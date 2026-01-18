// v8: QPSMAP を別途用意したファイルに置いたまま，全体を読み込まずに，一つずつ読み込んで kNN search (reranking) を行うオプションを追加．
//     QPSMAP_ON_SECONDARY_MEMORY を定義すると，全体は読み込まないで検索する．
// v7: QPSMAP を別途用意したファイルから読み込む
// v6:  v5_2 で使用した部分のみをできるだけ残し，不要部分をある程度削除して，マクロスイッチを減らしたバージョン
// v5_2: SISAP2025 Indexing Challenge で使用した最終バージョン
// v5_1: 立ち上がり（最初の質問に対する検索）が遅くなるという問題を解決
// v5: 検索時に対話型でなく，バッチ型で動かせるように変更
// v4: 解候補数をデータセットに対するppmdではなく，実際の候補数に変更．
// v3_1_1: latency (質問に対する応答時間) の分布を求める．
// v3_1: v2_1から派生：FTR_ON_MAIN_MEMORYに対応（最終段階検索で用いる特徴データをあらかじめ読み込んでおく）
// v2_1: 1st filtering で D~1 順の列挙に対応させる．
// v2: Filtering で interval_list （データ番号の区間（はじめ，長さ（or 終わり））の配列をスレッド数分）に対応させる．
// SMAP は圧縮形式のみ，射影距離計算も圧縮形式対応のみとし，それ以外は削除（未対応）
// UNPACKED_SMAP, USE_TABLE_FOR_UMPACKED に対応する部分を削除
// USE_INTERVAL が定義されているときは, 
// ともに INTERVAL_WITH_RUN, INTERVAL_WITH_PRIORITY, LOOP_CONTROL_BY_NUM_SKETCHESが定義されていることを前提とする．??
// そうでないときは，これまでのデータ番号のリストを用いる．
//
// v1: Filteringを並列化（ただし，2, 4, 8, 16, 32 スレッドのみ）
// sketch_and_smap_v0: sketch_filtering_v2 からの派生．SketchとSmapを用いた double-filtering
// ただし，v2にあった，D~inf 順と D~1 順の列挙は削除し，列挙パターン表を用いたハミング順列挙（Enhanced 版，2段階法対応）のみ
// 1st filtering: Enumeration in Hamming distance order（距離下限を考慮した2段階ハミング距離順列挙）
// 2nd filtering: Full-scan of candidates obtained by 1st filtering using qpsmap
// Enhanced 版は sketch_filtering_v2 で導入済．距離下限が小さいビットがONのパターンを先に列挙する．
// 例：下位6ビットのハミング距離順列挙で始めて，不足したら追加の上位ビットパターンの列挙を行う．
// sketch_filtering_v2: Filtering = parallel enumeration in D~inf order + selection by D~1
// フィルタリングを D~inf 順の並列列挙により上位のスケッチを選択し，それらの D~1 が上位のスケッチに絞り込む．
// FILTERING_BY_SKETCH_ENUMERATION_INF を定義しておく．
// 当面 NUM_K = 0 (フィルタリングのみ) or -2 (Filtering ので recall を求める) だけに対応し，実際の検索は行わない．
#include <stdio.h>
#include <string.h>
#include "parm.h"
#include "config.h"
#include "ftr.h"
#include "smap.h"
#include "sketch.h"
#include "quick.h"
#include "bit_op.h"
#include "e_time.h"
#include <stdlib.h>
#include <errno.h>
#include <unistd.h>
#include <sys/resource.h>

// ここは削除する
#ifdef UMPACKED_SMAP
// 非圧縮形式 qpsmap 射影像を用いる. Single-threadのみ
#error "This verson is not applicable for UMPACKED_SMAP."
#elif defined(USE_TABLE_FOR_UMPACKED)
// 表関数は非圧縮用を用いる．Single-threadのみ．qpsmap射影像は圧縮しているが，検索時には umpack して用いる
#error "This verson is not applicable for USE_TABLE_FOR_UMPACKED."
#endif

double filtering_cost_1st, filtering_cost_2nd; 
int ivl_count; // 1st filtering で求めた候補の区間数（候補数ではない）
void reset_filtering_cost(void) {
    filtering_cost_1st = filtering_cost_2nd = 0;
}

// #define QPSMAP_ON_SECONDARY_MEMORY

// /*
static interval_list *ivl = NULL;
static int intial_ivl_size_1st = 0;
static int intial_ivl_size_2nd = 0;
static struct_que_c2_n *que = NULL;
static answer_type *ans_buff = NULL;
static kNN_buffer *buff_p[NUM_THREADS] = {NULL};
//#if defined(QPSMAP_ON_SECONDARY_MEMORY) && defined(NUM_THREADS) && (NUM_THREADS > 1)
#if defined(QPSMAP_ON_SECONDARY_MEMORY) // && defined(NUM_THREADS) && (NUM_THREADS > 1)
// QPSMAP を全部読み込まずに，2次記憶に置いたままにするときで，マルチスレッド処理を行うときは，
// ある程度シングルスレッドでまとめ読みを行ってバッファに入れてから処理する
// シングルスレッドでは，この必要はないと思うが，使ってみよう．
#ifndef QPSMAP_BUFFER_SIZE
#define QPSMAP_BUFFER_SIZE 100
#endif

typedef struct {
    int num;            // バッファに読込まれているデータ数
//  int done;           // 2次フィルタリングでは，毎回すべてを用いるので不要．kNN検索の処理済みデータ数，バッファに読込んだ直後は 0  
    int interval_num;   // つぎにデータを読込む区間番号
    int next;           // つぎに読込む処理内の位置（QPSMAP_BUFFER_SIZEまで読んだら，残りは次回のまとめ読みで．
    tiny_int qpsmap_data[QPSMAP_BUFFER_SIZE][PACKED_QPSMAP_SIZE];
    int data_num[QPSMAP_BUFFER_SIZE];
} qpsmap_buffer_t;

static qpsmap_buffer_t *qpsmap_buffer = NULL;

void init_qpsmap_buffer(void) 
{
    if(qpsmap_buffer == NULL) {
        qpsmap_buffer = MALLOC(sizeof(qpsmap_buffer_t) * NUM_THREADS);
    }
}

void make_empty_qpsmap_buffer(void)
{
    for(int t = 0; t < NUM_THREADS; t++) {
        qpsmap_buffer[t].num = 0;
//      qpsmap_buffer[t].done = 0;
        qpsmap_buffer[t].interval_num = 0;
        qpsmap_buffer[t].next = 0;
    }
}

// ivl に未処理の候補が残っているかどうかを判定する．
int remain_ivl(void)
{
    for(int t = 0; t < NUM_THREADS; t++) {
        // qpsmap_buffer[t].interval_num は，つぎに処理する区間番号（interalリスト内での番号）
        // ivl->lg[t] は，スレッドに割り当てられている区間リストの長さ（区間数）
        if(qpsmap_buffer[t].interval_num < ivl->lg[t]) { // 残りがある．  
            return 1;
        }
    }
    return 0; // 全スレッドで処理が終わっている．
}

void read_qpsmap_to_buffer(void) 
{
    // 各スレッドのバッファに最大 QPSMAP_BUFFER_SIZE 個のデータを読み込む．並列処理にすると遅くなるので，ここでは，直列処理．
    // 残っているデータが QPSMAP_BUFFER_SIZE 個未満のときは，残っている分だけ．残りが無ければ何もしない．
    for(int t = 0; t < NUM_THREADS; t++) {
        interval *list = ivl->list + t * ivl->size; // 区間リスト（一本のlistをsize個ごとに区切ってそれぞれをスレッドで用いる）
        int lg = ivl->lg[t];                        // 区間数
        // スレッド t に割り当てられている区間は，list[0], ... , list[lg - 1]
        qpsmap_buffer_t *b = qpsmap_buffer + t;     // qpsmapバッファ
        int m = b->interval_num;                    // 次に読込むデータ番号が入っている区間番号，質問に対する2次フィルタリングが始まったときは，0（先頭）にしておく．
        // listのm番目の区間で表されるデータ番号は，list[m].start から list[m].run 個, 
        // つまり，list[m].start, list[m].start + 1, ... , list[m].start + list[m].run - 1
        // 前回の読込みが，区間の途中までで，残っているときは，つぎの読込みは，b->next番目から
        // バッファーにQPSMAP_BUFFER_SIZE個のqpsmapを読み込む．
        b->num = 0;
//      b->done = 0;
        while(1) {
            if(b->num >= QPSMAP_BUFFER_SIZE) { // バッファが満杯
                break;
            }
            if(m >= lg) { // 区間が残っていない
                break;
            }
            int n = list[m].run - b->next; // つぎに読込む区間に残っているデータ数
            if(n <= 0) { // 現在の区間は読み切っているので，次の区間に
                m++; 
                b->interval_num = m;
                b->next = 0;
                continue;
            }
            if(b->num + n > QPSMAP_BUFFER_SIZE) {   // 区間に残っているデータを全部読み込むとバッファがあふれる
                n = QPSMAP_BUFFER_SIZE - b->num;    // 丁度いっぱいになるだけ読み込む．
            }
            // n 個のデータを読み込む．その先頭は，スケッチ順で list[m].start + b->next
            read_qpsmap_record(list[m].start + b->next, n, (tiny_int *)(b->qpsmap_data + b->num));
            // 読み込んだデータ番号を記録する．
            for(int i = 0; i < n; i++) {
                b->data_num[b->num + i] = list[m].start + b->next + i;
            }
            b->num += n;
            b->next += n;
        }
    }
}
#endif

// size_1st = 1次フィルタリングの候補リストの大きさ（区間数）の最大値（interval_list の大きさ）
// size_2nd = 2次フィルタリンツが求める候補数の最大値（kNN_bufferの大きさ）
void init_space_for_double_filtering(int size_1st, int size_2nd)
{
//  int nnt = (1 << PARA_ENUM_INF); // 分割数（THREAD_PLUS > 0 のとき，nt * (1 << THREAD_PLUS)
    if(ivl != NULL) {
        fprintf(stderr, "init_space_for_double_filtering error: ivl != NULL\n");
        exit(1);
    }
    if(ans_buff != NULL) {
        fprintf(stderr, "init_space_for_double_filtering error: ans_buff[0] != NULL\n");
        exit(1);
    }
    if(buff_p[0] != NULL) {
        fprintf(stderr, "init_space_for_double_filtering error: buff_p[0] != NULL\n");
        exit(1);
    }
    if(que != NULL) {
        fprintf(stderr, "init_space_for_double_filtering error: que != NULL\n");
        exit(1);
    }
    ivl = new_interval_list(NUM_THREADS, size_1st);
    ivl->lg[0] = 0;
    ivl->list[0].start = 0;
    intial_ivl_size_1st = size_1st;
    intial_ivl_size_2nd = size_2nd;
    for(int t = 0; t < NUM_THREADS; t++) {
        buff_p[t] = new_kNN_buffer(size_2nd);
        buff_p[t]->buff[0].data_num = 0;
    }
    ans_buff = MALLOC(sizeof(answer_type) * size_2nd * NUM_THREADS);
    ans_buff[0].data_num = 0;
    fprintf(stderr, "MALLOC QUE: size = %lu\n", sizeof(struct_que_c2_n));
    use_system("VmSize");
    que = MALLOC(sizeof(struct_que_c2_n));
    fprintf(stderr, "MALLOC QUE: OK\n");
    use_system("VmSize");
    for(int i = 0; i < QSIZE; i += 1024) {
        que->element[i].key = 0;
        que->details[i].sk = 0;
    }
//    #if defined(QPSMAP_ON_SECONDARY_MEMORY) && defined(NUM_THREADS) && (NUM_THREADS > 1)
//    packed_qpsmap_data = MALLOC(sizeof(*qpsmap_buffer) * NUM_THREADS);
//
//    #endif
}
//#endif
// */

static void count_ivl(void) 
{
    ivl_count = 0;
    for(int t = 0; t < ivl->nt; t++) {
        ivl_count += ivl->lg[t];
    }
}

// qpsmap 射影像は圧縮したままで，表関数は圧縮用を用いる．
#if PARALLEL_ENUM == 0
// Double-Filtering（single-thread）
// 1st: スケッチ列挙(Hamming)によるフィルタリング（バケット（配列 idx と bkt）利用）（部分集合列挙の表を利用する）
// 2nd: qpsmap による射影距離を用いる（データの qpsmap 射影像は，ビット列を詰め合わせた圧縮表現を用いる）
int double_filtering_by_sketch_enumeration_hamming_and_qpsmap(
	struct_query_sketch *qs, 
    #if defined(USE_PACKED_3BIT) || defined(USE_PACKED_6BIT)
	unsigned int table_for_packed[][1 << QUANTIZE_BIT * TABLE_UNIT], 
    #else
	unsigned int table_for_packed[][256], 
    #endif
	struct_bucket *bucket, 
	tiny_int packed_qpsmap_data[][PACKED_QPSMAP_SIZE],          // QPSMAP_ON_SECONDARY_MEMORY のときは，NULL
	int data_num[], int num_candidates_1st, int num_candidates_2nd)
{
    #ifdef QPSMAP_ON_SECONDARY_MEMORY
    tiny_int packed_qpsmap_record[PACKED_QPSMAP_SIZE];
    #endif
    struct timespec tp1, tp2, tp3;
    clock_gettime(CLOCK_METHOD, &tp1);
    // Interval list を用いるときは，
    // ともに INTERVAL_WITH_RUN, INTERVAL_WITH_PRIORITY, LOOP_CONTROL_BY_NUM_SKETCHESが定義されていることを前提とする． 
    // まず，filtering_by_sketch_enumeration_hamming_interval を用いて，フィルタリングを行い，interval_list の形式で候補を求める．
    // 1st filtering の結果を用いて，2nd filtering を qpsmap を用いて行う．
    // int nt = 1; // single-thread
/*
    static interval_list *ivl = NULL;
    if(ivl == NULL) {
        int nnt = nt;
        #ifdef THREAD_PLUS
        nnt *= (1 << THREAD_PLUS);
        #endif
        ivl = new_interval_list(nnt, num_candidates_1st);            
    } else if(ivl->size < num_candidates_1st) {
        fprintf(stderr, "realloc_interval_list: ivl->size = %d, nc1 = %d\n", ivl->size, num_candidates_1st);
        realloc_interval_list(ivl, num_candidates_1st);
        fprintf(stderr, "realloc_interval_list: OK\n");
    }
*/
    #ifdef FILTERING_BY_SKETCH_ENUMERATION_HAMMING
//        int nc = filtering_by_sketch_enumeration_hamming_interval(qs, bucket, ivl, num_candidates_1st);
        filtering_by_sketch_enumeration_hamming_interval(qs, bucket, ivl, num_candidates_1st);
    #elif defined(FILTERING_BY_SKETCH_ENUMERATION_C2N)
/*
    static struct_que_c2_n *que = NULL;
    if(que == NULL) {
        que = MALLOC(sizeof(struct_que_c2_n));
        for(int i = 0; i < QSIZE; i += 1024) {
            que->element[i].key = 0;
            que->details[i].sk = 0;
        }
    }
*/
//        int nc = filtering_by_sketch_enumeration_c2_n_interval(qs, bucket, que, ivl, num_candidates_1st);
        filtering_by_sketch_enumeration_c2_n_interval(qs, bucket, que, ivl, num_candidates_1st);
    #else
    #error "FILTERING_BY_SKETCH_ENUMERATION_(HAMMING | C2N) should be defined"
    #endif
    clock_gettime(CLOCK_METHOD, &tp2);

    kNN_buffer *buff = new_kNN_buffer(num_candidates_2nd);
    int k = 0;

    #ifndef QPSMAP_BUFFER_SIZE
    for(int t = 0; t < ivl->nt; t++) {
        interval *list = ivl->list + t * ivl->size;
        int lg = ivl->lg[t];
        for(int n = 0; n < lg /* && k < num_candidates_1st */; n++) {
            int j = list[n].start;
            for(int r = 0; r < list[n].run /* && k < num_candidates_1st */; j++, r++, k++) {
                // ここで，スケッチ順に qpsmap が並んでいるので，質問と j 番目の qpsmap との部分復元射影距離を求める．
                #ifdef QPSMAP_ON_SECONDARY_MEMORY
                if(!read_qpsmap_record(j, 1, packed_qpsmap_record)) {
                    exit(1);
                }
                dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_record);
                #else
                dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_data[j]);
                #endif
                answer_type ans = (answer_type){ bucket->idx[j], p_dist };
                push_kNN_buffer(&ans, buff); // 元の順番でのデータ番号と射影距離の対を kNN_buffer に push
            }
        }
    }
    #else

    init_qpsmap_buffer();
    make_empty_qpsmap_buffer();

    while(remain_ivl()) {
        // qpsmap_buffer に候補データを読み込む
        read_qpsmap_to_buffer();
        // バッファに読み込んだ候補について，ｋNN選択をする．
        qpsmap_buffer_t *b = qpsmap_buffer;     // qpsmapバッファ
        for(int m = 0; m < b->num; m++, k++) {
            dist_type p_dist = projected_dist_packed_table(table_for_packed, b->qpsmap_data[m]);
            answer_type ans = (answer_type){ bucket->idx[b->data_num[m]], p_dist };
            push_kNN_buffer(&ans, buff); // 元の順番でのデータ番号と射影距離の対を kNN_buffer に push
        }
    }

    flush_kNN_buffer(buff);

    #endif

    flush_kNN_buffer(buff);
    for(int i = 0; i < num_candidates_2nd; i++) {
        data_num[i] = buff->buff[i].data_num;
    }
    free_kNN_buffer(buff);

    clock_gettime(CLOCK_METHOD, &tp3);
    filtering_cost_1st += e_time(&tp1, &tp2);
    filtering_cost_2nd += e_time(&tp2, &tp3);

    //#ifdef LOG_IVL_COUNT
    count_ivl();    
    //#endif

    return k;
}

#else // PARALLEL_ENUM == 0
/*
static interval_list *ivl = NULL;
static int intial_ivl_size_1st = 0;
static int intial_ivl_size_2nd = 0;
static struct_que_c2_n *que = NULL;
static answer_type *ans_buff = NULL;
static kNN_buffer *buff_p[NUM_THREADS] = {NULL};

// size_1st = 1次フィルタリングの候補リストの大きさ（区間数）の最大値（interval_list の大きさ）
// size_2nd = 2次フィルタリンツが求める候補数の最大値（kNN_bufferの大きさ）
void init_space_for_double_filtering(int size_1st, int size_2nd)
{
//  int nnt = (1 << PARA_ENUM_INF); // 分割数（THREAD_PLUS > 0 のとき，nt * (1 << THREAD_PLUS)
    if(ivl != NULL) {
        fprintf(stderr, "init_space_for_double_filtering error: ivl != NULL\n");
        exit(1);
    }
    if(ans_buff != NULL) {
        fprintf(stderr, "init_space_for_double_filtering error: ans_buff[0] != NULL\n");
        exit(1);
    }
    if(buff_p[0] != NULL) {
        fprintf(stderr, "init_space_for_double_filtering error: buff_p[0] != NULL\n");
        exit(1);
    }
    if(que != NULL) {
        fprintf(stderr, "init_space_for_double_filtering error: que != NULL\n");
        exit(1);
    }
    ivl = new_interval_list(NUM_THREADS, size_1st);
    ivl->lg[0] = 0;
    ivl->list[0].start = 0;
    intial_ivl_size_1st = size_1st;
    intial_ivl_size_2nd = size_2nd;
    for(int t = 0; t < NUM_THREADS; t++) {
        buff_p[t] = new_kNN_buffer(size_2nd);
        buff_p[t]->buff[0].data_num = 0;
    }
    ans_buff = MALLOC(sizeof(answer_type) * size_2nd * NUM_THREADS);
    ans_buff[0].data_num = 0;
    fprintf(stderr, "MALLOC QUE: size = %lu\n", sizeof(struct_que_c2_n));
    use_system("VmSize");
    que = MALLOC(sizeof(struct_que_c2_n));
    fprintf(stderr, "MALLOC QUE: OK\n");
    use_system("VmSize");
    for(int i = 0; i < QSIZE; i += 1024) {
        que->element[i].key = 0;
        que->details[i].sk = 0;
    }
}
//#endif
*/

// Double-Filtering（multi-thread）
// 1st: スケッチ列挙(Hamming)によるフィルタリング（バケット（配列 idx と bkt）利用）（部分集合列挙の表を利用する）
// 2nd: qpsmap による射影距離を用いる（データの qpsmap 射影像は，ビット列を詰め合わせた圧縮表現を用いる）
int double_filtering_by_sketch_enumeration_hamming_and_qpsmap(
	struct_query_sketch *qs, 
    #if defined(USE_PACKED_3BIT) || defined(USE_PACKED_6BIT)
	unsigned int table_for_packed[][1 << QUANTIZE_BIT * TABLE_UNIT], 
    #else
	unsigned int table_for_packed[][256], 
    #endif
	struct_bucket *bucket, 
	tiny_int packed_qpsmap_data[][PACKED_QPSMAP_SIZE],          // QPSMAP_ON_SECONDARY_MEMORY のときは，NULL 
	int data_num[], int num_candidates_1st, int num_candidates_2nd)
{
	int n = PARALLEL_ENUM;
//  int nnt = (1 << PARA_ENUM_INF); // 分割数（THREAD_PLUS > 0 のとき，nt * (1 << THREAD_PLUS)
    #ifdef THREAD_PLUS
        nnt *= (1 << THREAD_PLUS);
        #if PARA_ENUM_INF == 0
            if(n > THREAD_PLUS) {
                n = THREAD_PLUS;
            }
        #endif
    #endif

    #if defined(FILTERING_BY_SKETCH_ENUMERATION_C2N) && PARA_ENUM_INF == 0
    n = PARALLEL_ENUM;
    nnt = (1 << PARALLEL_ENUM);
    #endif

	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif

    #ifdef QPSMAP_ON_SECONDARY_MEMORY
    tiny_int packed_qpsmap_record[PACKED_QPSMAP_SIZE];
    #endif

    struct timespec tp1, tp2, tp3;

    clock_gettime(CLOCK_METHOD, &tp1);
    // Interval list を用いるときは，
    // ともに INTERVAL_WITH_RUN, INTERVAL_WITH_PRIORITY, LOOP_CONTROL_BY_NUM_SKETCHESが定義されていることを前提とする． 
    // まず，filtering_by_sketch_enumeration_hamming_interval を用いて，フィルタリングを行い，interval_list の形式で候補を求める．
    // 1st filtering の結果を用いて，2nd filtering を qpsmap を用いて行う．
    if(num_candidates_1st > intial_ivl_size_1st) {
        fprintf(stderr, "too large nc1 = %d > %d\n", num_candidates_1st, intial_ivl_size_1st);
        exit(1);
    }
    ivl->size = num_candidates_1st;
    static int first = 1;
    #ifdef FILTERING_BY_SKETCH_ENUMERATION_HAMMING
        if(first) {fprintf(stderr, "FILTERING_BY_SKETCH_ENUMERATION_HAMMING, USE_INTERVAL\n"); first = 0;}
        int nc = filtering_by_sketch_enumeration_hamming_interval(qs, bucket, ivl, num_candidates_1st);
        if(nc == 0) {fprintf(stderr, "nc = %d, num_candidates_1st = %d\n", nc, num_candidates_1st); getchar(); }
    #elif defined(FILTERING_BY_SKETCH_ENUMERATION_C2N)
        if(first) {fprintf(stderr, "FILTERING_BY_SKETCH_ENUMERATION_C2N, USE_INTERVAL\n"); first = 0;}
        int nc = filtering_by_sketch_enumeration_c2_n_interval(qs, bucket, que, ivl, num_candidates_1st);
    #else
        #error "FILTERING_BY_SKETCH_ENUMERATION_(HAMMING | C2N) should be defined"
    #endif
    clock_gettime(CLOCK_METHOD, &tp2);

    int num_data_1st = nc / nt;	// スレッドが 1st filtering で求めるデータ数
    num_candidates_1st = num_data_1st * nt;     // 1st filtering で求めるデータ数の合計（元の num_candidates_1st がスレッド数で割り切れないときに端数を切り捨てる）

    // 2nd filtering を個々のスレッドで kNN_buffer 法で行う
    int num_data_2nd = num_candidates_2nd; // スレッドに分けても同じ候補数を選ばせる
    if(num_data_2nd > intial_ivl_size_2nd) {
        fprintf(stderr, "too large nc2 = %d > %d\n", num_candidates_2nd, intial_ivl_size_2nd);
        exit(1);
    }
    for(int t = 0; t < nt; t++) {
        make_empty_kNN_buffer(buff_p[t]);
        buff_p[t]->k = num_data_2nd;
    }

    // 2nd filtering ...
    #ifdef QPSMAP_ON_SECONDARY_MEMORY

    init_qpsmap_buffer();
    make_empty_qpsmap_buffer();

    while(remain_ivl()) {
        // qpsmap_buffer に候補データを読み込む
        read_qpsmap_to_buffer();
        // バッファに読み込んだ候補について，ｋNN選択をする．
        #pragma omp parallel
        {
            int t = omp_get_thread_num(); // スレッド番号
            kNN_buffer *buff = buff_p[t];
            qpsmap_buffer_t *b = qpsmap_buffer + t;     // qpsmapバッファ

            // （注意）THREAD_PLUS には未対応
            for(int m = 0; m < b->num; m++) {
                dist_type p_dist = projected_dist_packed_table(table_for_packed, b->qpsmap_data[m]);
                answer_type ans = (answer_type){ bucket->idx[b->data_num[m]], p_dist };
                push_kNN_buffer(&ans, buff); // 元の順番でのデータ番号と射影距離の対を kNN_buffer に push
            }
        }
    }

    #pragma omp parallel
    {
        int t = omp_get_thread_num(); // スレッド番号
        kNN_buffer *buff = buff_p[t];
        answer_type *a_buff = ans_buff + t * num_data_2nd;
        flush_kNN_buffer(buff);
        for(int i = 0; i < num_data_2nd; i++) {
            if(i < buff->num) {
                a_buff[i] = buff->buff[i];
            } else {
                a_buff[i] = (answer_type) {0, INT_MAX};
            }
        }
    }
    #else
    #pragma omp parallel
    {
        int t = omp_get_thread_num(); // スレッド番号
        int m = 0; // 列挙したスケッチ数（パターン番号）
        int k = 0; // 1st filtering で求めたデータ数
        kNN_buffer *buff = buff_p[t];
        answer_type *a_buff = ans_buff + t * num_data_2nd;

        #ifndef THREAD_PLUS
        int tt = t;
        #elif PARA_ENUM_INF > 0
        for(int tt = t * (1 << THREAD_PLUS); tt < (t + 1) * (1 << THREAD_PLUS); tt++) 
        #else
        int tt = t;
        #endif
        {
            interval *list = ivl->list + tt * ivl->size;
            int lg = ivl->lg[tt];
            for(m = 0; m < lg /* && k < num_data_1st * FACTOR_INF3 */; m++) {
                int j = list[m].start;
                for(int r = 0; r < list[m].run /* && k < num_data_1st * FACTOR_INF3 */; j++, r++, k++) {
                    // ここで，スケッチ順に qpsmap が並んでいるので，質問と j 番目の qpsmap との部分復元射影距離を求める．
                    #ifdef QPSMAP_ON_SECONDARY_MEMORY
                    if(!read_qpsmap_record(j, packed_qpsmap_record)) {
                        exit(1);
                    }
                    dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_record);
                    #else
                    dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_data[j]);
                    #endif
                    //dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_data[j]);
                    answer_type ans = (answer_type){ bucket->idx[j], p_dist };
                    push_kNN_buffer(&ans, buff); // 元の順番でのデータ番号と射影距離の対を kNN_buffer に push
                }
            }
        }

        flush_kNN_buffer(buff);
        for(int i = 0; i < num_data_2nd; i++) {
            if(i < buff->num) {
                a_buff[i] = buff->buff[i];
            } else {
                a_buff[i] = (answer_type) {0, INT_MAX};
            }
        }
    }
    #endif

    quick_select_k_answer(ans_buff, 0, num_data_2nd * nt - 1, num_candidates_2nd);
    for(int i = 0; i < num_candidates_2nd; i++) {
        data_num[i] = ans_buff[i].data_num;
    }
    clock_gettime(CLOCK_METHOD, &tp3);
    filtering_cost_1st += e_time(&tp1, &tp2);
    filtering_cost_2nd += e_time(&tp2, &tp3);
//    #endif
    //#ifdef LOG_IVL_COUNT
    count_ivl();    
    //#endif

    return num_candidates_2nd;
}

#endif // PARALLEL_ENUM

#define NUM_NN 30
int main(int argc, char *argv[])
{
	int num_ftr_files = argc - 1;
	char **dataset_ftr_filename = argv + 1;
	char *pivot_file = PIVOT_FILE;
	char *smap_pivot_file = SMAP_PIVOT_FILE;
	char *bucket_filename = BUCKET_FILE;
	char *qpsmap_file = QPSMAP_FILE;
	char *qr_file[] = {QUERY_FILE, QUERY_2ND_FILE, QUERY_3RD_FILE};
    #ifdef SELF_EVAL
	char *an_file[] = {ANSWER_FILE, ANSWER_2ND_FILE, ANSWER_3RD_FILE};
    #endif
	char *result_filename = RESULT_FILE;
    #ifdef INPUT_HYPER_PARAMETER
    FILE *fp = fopen(INPUT_HYPER_PARAMETER, "r");
    #else
    FILE *fp = stdin;
    #endif
    #ifdef SUMMARY_FILE
    FILE *fp_summary = fopen(SUMMARY_FILE, "w");
    fprintf(stderr, "SUMMARY_FILE = %s\n", SUMMARY_FILE);
    if(fp_summary == NULL) {
        fprintf(stderr, "open error\n");
    } else {
        fprintf(stderr, "open OK\n");
    }
    #else
    FILE *fp_summary = NULL;
    #endif
    #ifdef PRINT_SEARCH_COST
    FILE *fp_search_cost = fopen(PRINT_SEARCH_COST, "w");
    #else
    FILE *fp_search_cost = NULL;
    #endif

	make_bitcnt_tbl(8);

	fprintf(stderr, "PJT_DIM = %d, SMAP_DIM = %d\n", PJT_DIM, SMAP_DIM);

	int num_query_files = sizeof(qr_file) / sizeof(qr_file[0]);
	query_type *qr[num_query_files];
    #ifdef SELF_EVAL
	answer_type_NN *correct_answer[num_query_files];
    #endif
	int num_qr[num_query_files];
	struct_dataset *ds_query[num_query_files];
	for(int m = 0; m < num_query_files; m++) {
		if(strcmp(qr_file[m], "NONE") == 0) {
			num_query_files = m;
			break;
		}
		use_system("VmSize");
		ds_query[m] = read_dataset_n(1, &qr_file[m]);
		num_qr[m] = ds_query[m]->num_data;
		printf("read query file (%s) OK. the number of queries = %d: ", qr_file[m], num_qr[m]);
		use_system("VmSize");
		qr[m] = (query_type *)malloc(sizeof(query_type) * num_qr[m]);
		for(int i = 0; i < num_qr[m]; i++) {
			qr[m][i] = (query_type) { i, ds_query[m]->ftr_id[i].ftr };
		}
        #ifdef SELF_EVAL
		correct_answer[m] = read_correct_answer_NN(an_file[m], num_qr[m]);
		printf("read correct answer (%s) OK. ", an_file[m]);
		use_system("VmSize");
        #endif
	}

	#if defined(PARTITION_TYPE_QBP)
	pivot_type *pivot = new_pivot(QBP);
	#elif defined(PARTITION_TYPE_PQBP)
	pivot_type *pivot = new_pivot(PQBP);
	#endif
	read_pivot(pivot_file, pivot);
	fprintf(stderr, "read pivot OK\n");

	smap_pivot_type *smap_pivot = new_smap_pivot(PQBP);
    read_smap_pivot(smap_pivot_file, smap_pivot);
	printf("read smap pivot OK. ");

	#if !defined(FILTERING_BY_SKETCH_ENUMERATION_HAMMING) && !defined(FILTERING_BY_SKETCH_ENUMERATION_C2N)
    #error "This program is only for 1st filtering by sketch enumeration in Hamming distance or c2_n."
	#endif
	struct_bucket *bucket_ds = read_bucket(bucket_filename);
	fprintf(stderr, "read bucket OK, ");
	int num_data = bucket_ds->num_data;
	fprintf(stderr, "number of data = %d\n", num_data);

    // NUM_K = 0 -> filtering only, -1 -> scoreing only, -2 recall by filtering only without search
    #if NUM_K <= 0
        #if NUM_K == -2
            int num_top_k = 1;
        #else
    	    int num_top_k = NUM_K;
        #endif
    #else
        int num_top_k = NUM_K;
    #endif // NUM_K

    // 複数に分かれた特徴データファイル対応の読込の準備（qpsmap像を作るときと最終段階での実距離計算による検索のときに使用）
    fprintf(stderr, "open ftr files (filename = %s, num_files = %d)\n", dataset_ftr_filename[0], num_ftr_files);
    #if !defined(_OPENMP) || NUM_THREADS < 1
    struct_multi_ftr *mf = open_multi_ftr(num_ftr_files, dataset_ftr_filename, BLOCK_SIZE);
    if(num_data != mf->num_data) {
        fprintf(stderr, "bucket (filename = %s, num_data = %d) is not compatible with datasets (filename = %s, ... , num_data = %d)\n", bucket_filename, num_data, dataset_ftr_filename[0], mf->num_data);
        return -1;
    }
    #else

    #endif

	dataset_handle dh;

    #ifndef NUM_THREADS_PREPERATION
    #define NUM_THREADS_PREPERATION 16
    #endif

    #ifdef _OPENMP
        #if NUM_THREADS > NUM_THREADS_PREPERATION
		dh.num_threads = NUM_THREADS;
        #else
		dh.num_threads = NUM_THREADS_PREPERATION;
        #endif
	#else
		dh.num_threads = 1;
	#endif

    dh.sorted = 0; // FTR is NOT sorted. Arrangement is as is for double filtering 
	dh.ftr_on = SECONDARY_MEMORY;
	dh.mf = (struct_multi_ftr **)malloc(sizeof(struct_multi_ftr *) * dh.num_threads);
	for(int t = 0; t < dh.num_threads; t++) {
		dh.mf[t] = open_multi_ftr(num_ftr_files, dataset_ftr_filename, BLOCK_SIZE);
	}
	dh.ds = NULL;
	if(num_data != dh.mf[0]->num_data) {
		fprintf(stderr, "bucket (filename = %s, num_data = %d) is not compatible with datasets (filename = %s, ... , num_data = %d)\n", bucket_filename, num_data, dataset_ftr_filename[0], dh.mf[0]->num_data);
		return -1;
	}

    //  検索をする
	#ifdef _OPENMP
		#ifdef NUM_THREADS
            #if NUM_THREADS > NUM_THREADS_PREPERATION
        		omp_set_num_threads(NUM_THREADS);
            #else
        		omp_set_num_threads(NUM_THREADS_PREPERATION);
            #endif
		#endif
	#endif

    #ifdef QPSMAP_ON_SECONDARY_MEMORY
    // データの圧縮表現のqpsmapは，2次記憶に置いたままで，2nd filtering を行う．
    tiny_int (*packed_qpsmap_data)[PACKED_QPSMAP_SIZE] = NULL; 
    #else
    // データを圧縮表現のqpsmapとして求めておく．ただし，1st filtering で用いるスケッチ順にソートしておく．
	fprintf(stderr, "before malloc packed quantized images of data.\n");
	use_system("VmSize");
    tiny_int (*packed_qpsmap_data)[PACKED_QPSMAP_SIZE]; // 圧縮形式の量子化射影像はデータのみ．質問は量子化していない射影像だけを使用する．
    packed_qpsmap_data = MALLOC(sizeof(tiny_int) * PACKED_QPSMAP_SIZE * num_data);
	fprintf(stderr, "malloc packed quantized images of data OK.\n");
	use_system("VmSize");
    #endif

    double offset[SMAP_DIM], slice[SMAP_DIM]; // 2nd filtering で用いる qpsmap のためのパラメタ
    #ifdef QPSMAP_ON_SECONDARY_MEMORY
    // このときは，qpsmap_fileのためのFILEやqpsmap_headerは，smap.c内の外部変数にしておく（将来は，書き換えた方がよいかも）．
    if(!open_qpsmap_file(qpsmap_file, offset, slice)) {
        fprintf(stderr, "cannot open qpsmap file = %s\n", qpsmap_file);
        return -1;
    }
    #else
    qpsmap_header qpsmap_hd;
    FILE *fp_qpsmap = fopen(qpsmap_file, "r");
    if(!fp_qpsmap) {
        fprintf(stderr, "cannot open qpsmap file = %s\n", qpsmap_file);
        return -1;
    }
    if(!read_qpsmap(&qpsmap_hd, offset, slice, (tiny_int *)packed_qpsmap_data, num_data, fp_qpsmap)) {
        return -1;
    }
    #endif

#ifndef MAX_NUM_CANDIDATES_1ST
#define MAX_NUM_CANDIDATES_1ST 4000000
#endif

#ifndef MAX_NUM_CANDIDATES_2ND
#define MAX_NUM_CANDIDATES_2ND 3000
#endif

use_system("VmSize");
	fprintf(stderr, "malloc data_num_candidates and data_num_org OK: num_data = %d.\n", num_data);
	int *data_num_candidates = MALLOC(sizeof(int) * MAX_NUM_CANDIDATES_1ST); // double-filtering で求めた候補データのデータ番号を格納する配列
    #ifdef FTR_ON_MAIN_MEMORY
	int *data_num_org = MALLOC(sizeof(int) * num_data); // データ番号をそのまま順番に格納する配列（データをファイルから読み込むときに使用する．
//	fprintf(stderr, "malloc data_num_candidates and data_num_org OK: num_data = %d.\n", num_data);
    #endif
	use_system("VmSize");

    #ifdef FTR_ON_MAIN_MEMORY
    struct_ftr_id *ftr_id = MALLOC(sizeof(struct_ftr_id) * num_data);
	fprintf(stderr, "malloc ftr_id OK: num_data = %d.\n", num_data);
	use_system("VmSize");
    #endif

	// いったん，元の順序のままで読み込んで，圧縮形式のSMAPに変換して，その後で並べ替える．
    #ifdef FTR_ON_MAIN_MEMORY
    for(int i = 0; i < num_data; i++) { data_num_org[i] = i; } // 特徴データを元の順序で読み込むためのデータ番号の配列としても使用する．
    #endif
    for(int i = 0; i < MAX_NUM_CANDIDATES_1ST; i++) { data_num_candidates[i] = i; } // 配列をRAMに置くためのダミーアクセス．

    struct timespec tread1, tread2;
    clock_gettime(CLOCK_METHOD, &tread1);

    #ifdef FTR_ON_MAIN_MEMORY
    // 特徴データをRAMに読み込む
    #ifdef _OPENMP
    omp_set_num_threads(NUM_THREADS_PREPERATION);
    #pragma omp parallel for
    #endif
    for(int i = 0; i < num_data; i++) {
        // ここは，マルチスレッドで並列化した方がよいかも．
        // ある程度まとめて（シングルスレッドで）連続に読み込んでから，並列処理する．（未実装11/19時点）
//        unsigned char quantized_projected_data[SMAP_DIM]; // 非圧縮形式のデータのqpsmap（作業用：つぎつぎに圧縮形式に変換するので，1個分のみで，ループ内の局所変数にする）
        #ifndef _OPENMP
 		struct_ftr_id *ftr_id_p = get_next_ftr_id_from_multi_ftr(mf, data_num_org, i, num_data);
        #else
        int t = omp_get_thread_num();
 		struct_ftr_id *ftr_id_p = get_next_ftr_id_from_multi_ftr(dh.mf[t], data_num_org, i, num_data);
        #endif

        #ifdef FTR_ON_MAIN_MEMORY
        memcpy(&ftr_id[i], ftr_id_p, sizeof(struct_ftr_id));
        #endif

        #ifdef _OPENMP
        int nt = omp_get_num_threads();
        int nd = num_data / nt;
        #else
        int t = 0; 
        int nd = num_data;
        #endif
        if(t == 0 && (i + 1) % (nd / 100) == 0) {
            fprintf(stderr, "read ftr data %4d%%\r", (i + 1) / (nd / 100));
        }
    }
    clock_gettime(CLOCK_METHOD, &tread2);
    fprintf(stderr, "\ndone: %.4lf (sec)\n", e_time(&tread1, &tread2));
    #endif 

// 2段階検索で実際に特徴データを2次記憶から読み込んで処理を行う場合は，ランダムアクセスになるので，ブロック読込みは非効率．
// ブロック読込みをやめるために，block_size = 1 に変更する．（一旦すべてcloseしてから，再度openする方法では，なぜか速度が落ちる）
	for(int t = 0; t < dh.num_threads; t++) {
        FREE(dh.mf[t]->ftr_id, sizeof(struct_ftr_id) * dh.mf[t]->block_size);
        FREE(dh.mf[t]->data_num, sizeof(int) * dh.mf[t]->block_size);
		dh.mf[t]->block_size = 1;
        dh.mf[t]->ftr_id = MALLOC(sizeof(struct_ftr_id));
        dh.mf[t]->data_num = NULL;
	}

    clock_gettime(CLOCK_METHOD, &tread2);
    fprintf(stderr, "\ndone: %.4lf (sec)\n", e_time(&tread1, &tread2));

	char line[4000], *cmdline, *file_name;
	FILE *fp2;

    #if NUM_K > 0 || NUM_K == -2
        int max_num_qr = 0;
        for(int f = 0; f < num_query_files; f++) {
            if(num_qr[f] > max_num_qr) {
                max_num_qr = num_qr[f];
            }
        }
        kNN_buffer **top_k = MALLOC(sizeof(kNN_buffer *) * max_num_qr);
        for(int q = 0; q < max_num_qr; q++) { top_k[q] = new_kNN_buffer(num_top_k); }
        for(int q = 0; q < max_num_qr; q++) { top_k[q]->buff[0].dist = 0; }
    #endif

    struct_query_sketch query_sketch;
    smap_element_type q_smap[SMAP_DIM]; // 質問のSMAP射影像 q_smap （質問ごとに射影像を作成する．∵検索コストに含める必要がある）
    unsigned int table[SMAP_DIM][1 << QUANTIZE_BIT];  // 非圧縮表現のデータの射影距離計算のための表．
    #ifdef QUANTIZE_MIXED_MOD3
        #ifdef USE_PACKED_6BIT
        unsigned int table_for_packed[44][64]; // 圧縮形式のqpsmapと質問の射影距離のための表関数
        #else
        unsigned int table_for_packed[(SMAP_DIM + 2) / 3][256]; // 圧縮形式のqpsmapと質問の射影距離のための表関数
        #endif
    #elif defined(USE_PACKED_3BIT) || defined(USE_PACKED_6BIT)
        unsigned int table_for_packed[(SMAP_DIM * QUANTIZE_BIT + QUANTIZE_BIT * TABLE_UNIT - 1)/ (QUANTIZE_BIT * TABLE_UNIT)][1 << QUANTIZE_BIT * TABLE_UNIT]; // 圧縮形式のqpsmapと質問の射影距離のための表関数
    #else
        unsigned int table_for_packed[SMAP_DIM * QUANTIZE_BIT / 8][256]; // 圧縮形式のqpsmapと質問の射影距離のための表関数
    #endif

    for(int m = 0; m < num_query_files; m++) {
        // 1st filtering のための query_sketch を作る．
        set_query_sketch(&query_sketch, &qr[m][0], pivot);

        // 2nd filtering のためのデータの圧縮形式の qpsmap との射影距離を計算するための表関数 table_for_packed を作る．
        psmap(q_smap, qr[m][0].ftr, smap_pivot); // 質問（ftr）の smap 射影像 q_smap を求める．
        make_table_for_query_p(q_smap, table, offset, slice, 0.1 * SCORE_P_2ND); // 射影像を用いて，射影距離計算の表 table を作成する．
        make_table_for_packed_data(table, table_for_packed, offset, slice); // 圧縮表現データのための表 table_for_packed を作成する．
    }

    #ifdef STATIC_KNN_BUFFER_FOR_SEARCH
	init_search_kNN_on_ram(num_top_k);
    #endif

//#if PARALLEL_ENUM != 0
    init_space_for_double_filtering(MAX_NUM_CANDIDATES_1ST, MAX_NUM_CANDIDATES_2ND);
	use_system("VmSize");
//#endif

    if(fp_summary != NULL) {
        #ifdef SELF_EVAL
        fprintf(fp_summary, "trial, query, width, q_bit, ftr_on, nc1, nc2, recall@1, recall@30, filtering, 1st(sec), 2nd(sec), kNN(sec), ave(ms/q), stdev(ms/q), min(ms/q), max(ms/q)\n");
        #else
        fprintf(fp_summary, "trial, query, width, q_bit, ftr_on, nc1, nc2, filtering, 1st(sec), 2nd(sec), kNN(sec), ave(ms/q), stdev(ms/q), min(ms/q), max(ms/q)\n");
        #endif
    }
    int trial = 0;

    // NUM_Q が質問ファイルの質問数の10分の1以下のときは，trial が進むたびに，つぎの質問に切り替えていく．
    // たとえば，num_qr = 1000, NUM_Q = 100 のときは，
    // trial = 1 => q = 0, ... , 99
    // trial = 2 => q = 100, ... , 199,
    // ...
    // trial = 10 => q = 900, ... , 999
    // trial = 11 => q = 0, ... , 99
    int q_1st = 0; // 最初の質問
    while(1) {
        if(fp == stdin) fprintf(stderr, "nc1 ? ");
        cmdline = fgets(line, 900, fp);
		if(fp != stdin) fprintf(stderr, "line = %s\n", line);
        if(cmdline == NULL) {		// EOF 
			if(fp == stdin) break;	// stdin なら終了
			fclose(fp);
            #ifdef INPUT_HYPER_PARAMETER
            break;
            #else
			fp = stdin;				// ファイルから stdin に戻す
			continue;
            #endif
		}
		if(fp == stdin && line[0] == '<') {	// 入力を切り替える
			for(file_name = line + 1; *file_name == ' '; file_name++);	// '<' に続く空白をスキップ
			file_name[strlen(file_name) - 1] = 0;						// 行末の '\n' を除去
			if((fp2 = fopen(file_name, "r")) == NULL) {
				fprintf(stderr, "Cannot open, file = %s\n", file_name);
			} else {
				fp = fp2;
			}
			continue;
		}
        int nc1 = atoi(line);
        if(fp == stdin) fprintf(stderr, "nc1 = %d\n", nc1);
        if(nc1 < 0) {				// 負の値が来たら終了
			if(fp == stdin) break;	// stdin のときは全体を終了
			fclose(fp);
            #ifdef INPUT_HYPER_PARAMETER
            break;
            #else
			fp = stdin;				// ファイルから stdin に戻す
			continue;
            #endif
		} else if(nc1 == 0) {
            continue;
        } else if(nc1 > MAX_NUM_CANDIDATES_1ST) {
            fprintf(stderr, "too large nc1 (%d) > %d\n", nc1, MAX_NUM_CANDIDATES_1ST);
            continue;
        }
        if(fp == stdin) fprintf(stderr, "nc2 ? ");
        cmdline = fgets(line, 900, fp);
        if(cmdline == NULL) {		// EOF
			if(fp == stdin) break;
			fclose(fp);
			fp = stdin;
			continue;
		}
        int nc2 = atoi(line);
        if(fp == stdin) fprintf(stderr, "nc2 = %d\n", nc2);
        if(nc2 > MAX_NUM_CANDIDATES_2ND) {
            fprintf(stderr, "too large nc2 (%d) > %d\n", nc2, MAX_NUM_CANDIDATES_2ND);
            continue;
        }        
        if(fp == stdin) fprintf(stderr, "OK -> <enter>, reset -> -1");
        cmdline = fgets(line, 900, fp);
        if(cmdline == NULL) {		// EOF
			if(fp == stdin) break;
			fclose(fp);
            #ifdef INPUT_HYPER_PARAMETER
            break;
            #else
			fp = stdin;				// ファイルから stdin に戻す
			continue;
            #endif
		}
        if(line[0] != '\n') continue;	// 空行をスキップ

        trial++;

//      int num_candidates_1st = nc1 * 0.000001 * num_data;
//		int num_candidates_2nd = nc2 * 0.000001 * num_data;
        int num_candidates_1st = nc1;
		int num_candidates_2nd = nc2;
//      if(num_candidates_2nd <= 30) num_candidates_2nd = 30;
        fprintf(stderr, "num_data = %d, num_candidates_1st = %d, num_candidates_2nd = %d\n", num_data, num_candidates_1st, num_candidates_2nd);
        #ifdef SELF_EVAL
        double total_filtering = 0, total_kNN = 0, total_total = 0, total_recall = 0;
        #else
        double total_filtering = 0, total_kNN = 0, total_total = 0;
        #endif
 
        for(int m = 0; m < num_query_files; m++) {
			use_system("VmSize");
			int num_queries = num_qr[m];
			#if NUM_K > 0 || NUM_K == -2
				for(int q = 0; q < num_queries; q++) { make_empty_kNN_buffer(top_k[q]); }
			#endif

			#ifdef NUM_Q
			if(NUM_Q != 0 && NUM_Q < num_queries) { num_queries = NUM_Q; }
			#endif
        
            #ifdef QPSMAP_ON_SECONDARY_MEMORY
            // QPSMAPがSSDのときは，検索が遅いので，質問数を100個に減らして実行
            //if(num_queries > 100) {
            //    num_queries = 100;
            //}
            #endif
            struct timespec tp1, tp2, tp3;
            clock_gettime(CLOCK_METHOD, &tp1);
            double e_time_filtering = 0, e_time_kNN = 0, e_time_total = 0;
            reset_filtering_cost();

            double trial_search_cost[num_queries];
            double trial_filtering_cost[num_queries];
            int trial_ivl_count[num_queries];
            #ifdef SELF_EVAL
			int found = 0;
            #endif

//            for(int q = 0; q < num_queries; q++) {
            int q; // もともとの質問ファイル内での質問番号 （q_1stから始まる．Trialが進むと，q_1stが0でなくなる）
            int q_local; // Trial内での局所的な質問番号（毎回0から始まる）
//            fprintf(stderr, "search starts (num_queries = %d): from %d to %d\n", num_queries, q_1st, q_1st + num_queries);
            for(q_local = 0; q_local < num_queries; q_local++) {
                q = q_local + q_1st;
                // 1st filtering のための query_sketch を作る．
	            clock_gettime(CLOCK_METHOD, &tp1);
                set_query_sketch(&query_sketch, &qr[m][q], pivot);

                // 2nd filtering のためのデータの圧縮形式の qpsmap との射影距離を計算するための表関数 table_for_packed を作る．
                psmap(q_smap, qr[m][q].ftr, smap_pivot); // 質問（ftr）の smap 射影像 q_smap を求める．
                make_table_for_query_p(q_smap, table, offset, slice, 0.1 * SCORE_P_2ND); // 射影像を用いて，射影距離計算の表 table を作成する．
                make_table_for_packed_data(table, table_for_packed, offset, slice); // 圧縮表現データのための表 table_for_packed を作成する．

                // double-filtering でデータ番号を求める．
                double_filtering_by_sketch_enumeration_hamming_and_qpsmap(&query_sketch, table_for_packed, bucket_ds, packed_qpsmap_data, data_num_candidates, num_candidates_1st, num_candidates_2nd);
                trial_ivl_count[q_local] = ivl_count;
                clock_gettime(CLOCK_METHOD, &tp2);
				e_time_filtering += e_time(&tp1, &tp2);
                trial_filtering_cost[q_local] = e_time(&tp1, &tp2);

                #if NUM_K == 0
                    // Filtering のみで，そのコストを求める．
                #elif NUM_K == -2
                    // 実際には検索を行わず，正解情報との照合で recall を求める．
                    #ifdef SELF_EVAL
    				found += answer_check(&correct_answer[m][q], num_candidates_2nd, data_num_candidates, top_k[q]);
                    #endif
                //    #define CHECK_NOT_FOUND
                    // Filteringで見つけた解候補に正解が見つからないとき，候補の実距離，射影距離（D~1）を調べる．
                    // ざっと見た限りではあるが，正解の射影距離（D~1）がかなり大きい傾向にあるようだ．
                    // Hammingの列挙がD~1が小さいものからになるようにしているので，あえて，D~1が大きいものも列挙されるようにしない限り，正解が候補に含まれることはないと思われる．
                    // QBP（PQBPの直積分解は用いないQBP）の射影距離は，D~infであれば距離下限が保証されるが，D~1やD~2とすると，射影距離が大きくなることは避けられない．
                    // それにも関わらず，全体としては，D~pの方が精度（recall）が高くなる理由については，再考が必要である．
                    #ifdef CHECK_NOT_FOUND
                        if(top_k[q]->k_nearest == UINT_MAX) {
                            fprintf(stderr, "query[%d] is not found. bd's sorted are:\n", q);
                            // ピボットと分割境界の最小距離を昇順にソートして表示．
                            for(int j = 0; j < PJT_DIM; j++) {
                                fprintf(stderr, " %d", query_sketch.bd[query_sketch.idx[j]]);
                            }
                            fprintf(stderr, "\n");
                            int correct = correct_answer[m][q].data_num;
                            sketch_type sk = 0;
                            mf->read_in = 0; mf->next = 0; // mf のバッファをキャンセル・リセット
                    		struct_ftr_id *ftr_id_p = get_next_ftr_id_from_multi_ftr(mf, data_num_org, correct, num_data);
                            sk = data_to_sketch(ftr_id_p->ftr, pivot);
                            fprintf(stderr, "p_dist of correct answer = %u, data_ID = %d, r_dist = %u\n", priority(sk, &query_sketch), ftr_id_p->data_id, dist_L2(ftr_id_p->ftr, qr[m][q].ftr, FTR_DIM));
                            dist_type min_p_dist = UINT_MAX, p_dist, min_r_dist = UINT_MAX, r_dist;
                            mf->read_in = 0; mf->next = 0; // mf のバッファをキャンセル・リセット
                            for(int k = 0; k < num_candidates_2nd; k++) {
                                ftr_id_p = get_next_ftr_id_from_multi_ftr(mf, data_num_candidates, k, num_candidates_2nd);
                                sk = data_to_sketch(ftr_id_p->ftr, pivot);
                                if((p_dist = priority(sk, &query_sketch)) < min_p_dist) {
                                    min_p_dist = p_dist;
                                }
                                if((r_dist = dist_L2(ftr_id_p->ftr, qr[m][q].ftr, FTR_DIM)) < min_r_dist) {
                                    min_r_dist = r_dist;
                                }
                            }
                            fprintf(stderr, "min_p_dist of candidates = %u\n", min_p_dist);
                            fprintf(stderr, "min_r_dist of candidates = %u\n", min_r_dist);
                            fprintf(stderr, "r_dist of correct_answer = %u\n", correct_answer[m][q].dist);
                            getchar();
                        }
                    #endif
                #elif NUM_K > 0
                    // k = NUM_K として，double-filtering で求めた nc2 個の候補から k-NN 検索を行う．
                    if(top_k[q] == NULL) {
                        fprintf(stderr, "もしこれが表示されたら変なので，停止します\n"); exit(1);
                        top_k[q] = new_kNN_buffer(num_top_k);
                    }
                    if(num_top_k == 1) {
                        #ifndef FTR_ON_MAIN_MEMORY
                        search_NN(&dh, &qr[m][q], num_candidates_2nd, data_num_candidates, top_k[q]);
                        #else
                        search_NN_on_ram(ftr_id, &qr[m][q], num_candidates_2nd, data_num_candidates, top_k[q]);
                        #endif
                        #ifdef SELF_EVAL
                        found += correct_answer[m][q].dist[0] == top_k[q]->buff[0].dist;
//                      trial_found[q] = correct_answer[m][q].dist[0] == top_k[q]->buff[0].dist;
//                      trial_dist[q] = top_k[q]->buff[0].dist; 
                        #endif
                    } else {
                        #ifndef FTR_ON_MAIN_MEMORY
                        search_kNN(&dh, &qr[m][q], num_candidates_2nd, data_num_candidates, top_k[q]);
                        #else
                        search_kNN_on_ram(ftr_id, &qr[m][q], num_candidates_2nd, data_num_candidates, top_k[q]);
                        #endif
                        #ifdef SELF_EVAL
                        for(int d = 0; d < num_top_k; d++) {
                            if(correct_answer[m][q].dist[0] == top_k[q]->buff[d].dist) {
                                found++;
                                break;
                            }
                        }
//                      trial_found[q] = correct_answer[m][q].dist[0] == top_k[q]->buff[0].dist;
//                      trial_dist[q] = top_k[q]->buff[0].dist; 
                        #endif
                    }
                #endif
				clock_gettime(CLOCK_METHOD, &tp3);
				e_time_kNN += e_time(&tp2, &tp3);
				e_time_total += e_time(&tp1, &tp3);
//                trial_search_cost[q] = e_time(&tp1, &tp3);
                trial_search_cost[q_local] = e_time(&tp1, &tp3);
            }
//          printf("trial, f, query, found,  dist,  filtering, total_cost, average, stdev\n");
            double sum = 0, sum2 = 0, ave, stdev, cost_min = 1000, cost_max = 0;
            for(int q = 0; q < num_queries; q++) {
                sum += trial_search_cost[q];
                sum2 += trial_search_cost[q] * trial_search_cost[q];
                if(trial_search_cost[q] < cost_min) cost_min = trial_search_cost[q];
                if(trial_search_cost[q] > cost_max) cost_max = trial_search_cost[q];
            }
            ave = sum / num_queries;
            stdev = sqrt(sum2 / num_queries - ave * ave);
            #ifdef PRINT_SEARCH_COST
            fprintf(fp_search_cost, "trial, file, query_num, filtering_cost, search_cost, ivl_count\n");
            for(int q = 0; q < num_queries; q++) {
                fprintf(fp_search_cost, "%d, %d, %d, %.4lf, %.4lf, %d\n", trial, m, q, trial_filtering_cost[q], trial_search_cost[q], trial_ivl_count[q]);
            }
            fprintf(fp_search_cost, "%d, %d, summary, %.4lf, %.4lf, %.4lf, %.4lf\n", trial, m, ave * 1000, stdev * 1000, cost_min * 1000, cost_max * 1000);
            #endif

//            #ifdef PRINT_KNN_RESULT
//            for(int q = 0; q < num_queries; q++) {
//                printf("%5d, %1d, %5d, %5d, %5d, %10.3le, %10.3le\n", trial, m, q, trial_found[q], trial_dist[q], trial_filtering_cost[q], trial_search_cost[q]);
//            }
//            #endif

            #ifdef SELF_EVAL
//            for(int x = 0; x < 5; x++) {
//                printf("q = %d, dist = %u\n", x, correct_answer[m][x].dist); 
//            }
//            getchar();
            double recall_1 = recall_kNN_1(num_queries, correct_answer[m] + q_1st, top_k + q_1st);
            double recall_k = recall_kNN(num_queries, correct_answer[m] + q_1st, top_k + q_1st);
            printf("filtering, %.4lf, kNN, %.4lf, total, %.4lf, ave = %.4lf (ms/q), stdev = %.4lf (ms/q), recall_1, %.1lf, found = %d, recall_k, %.1lf\n", 
                e_time_filtering, e_time_kNN, e_time_total, ave * 1000, stdev * 1000, recall_1, found, recall_k);
            printf("filtering cost: 1st = %.4lf (ms/q), 2nd = %.4lf (ms/q), kNN (reranking) cost: %.4lf (ms/q)\n", 
                (double)filtering_cost_1st / num_queries * 1000, (double)filtering_cost_2nd / num_queries * 1000, e_time_kNN / num_queries * 1000);
            total_filtering += e_time_filtering; total_kNN += e_time_kNN, total_total += e_time_total, total_recall += recall_1;
            if(fp_summary != NULL) {
                #ifdef FTR_ON_MAIN_MEMORY
                fprintf(fp_summary, "%d, %d, %d, %d, RAM, %d, %d, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf \n", 
                    trial, m, PJT_DIM, QUANTIZE_BIT, nc1, nc2, recall_1, recall_k, e_time_filtering, filtering_cost_1st, filtering_cost_2nd, e_time_kNN, ave * 1000, stdev * 1000, cost_min * 1000, cost_max * 1000);
                #else
                fprintf(fp_summary, "%d, %d, %d, %d, SSD, %d, %d, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf \n", 
                    trial, m, PJT_DIM, QUANTIZE_BIT, nc1, nc2, recall_1, recall_k, e_time_filtering, filtering_cost_1st, filtering_cost_2nd, e_time_kNN, ave * 1000, stdev * 1000, cost_min * 1000, cost_max * 1000);
                #endif
            }
            #else
            printf("filtering, %.4lf, kNN, %.4lf, total, %.4lf, ave = %.4lf (ms/q), stdev = %.4lf (ms/q)\n", 
                e_time_filtering, e_time_kNN, e_time_total, ave * 1000, stdev * 1000);
            printf("filtering cost: 1st = %.4lf (ms/q), 2nd = %.4lf (ms/q)\n", (double)filtering_cost_1st / num_queries * 1000, (double)filtering_cost_2nd / num_queries * 1000);
            total_filtering += e_time_filtering; total_kNN += e_time_kNN, total_total += e_time_total;
            if(fp_summary != NULL) {
                #ifdef FTR_ON_MAIN_MEMORY
                fprintf(fp_summary, "%d, %d, %d, %d, RAM, %d, %d, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf \n", 
                    trial, m, PJT_DIM, QUANTIZE_BIT, nc1, nc2, e_time_filtering, filtering_cost_1st, filtering_cost_2nd, e_time_kNN, ave * 1000, stdev * 1000, cost_min * 1000, cost_max * 1000);
                #else
                fprintf(fp_summary, "%d, %d, %d, %d, SSD, %d, %d, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf, %.4lf \n", 
                    trial, m, PJT_DIM, QUANTIZE_BIT, nc1, nc2, e_time_filtering, filtering_cost_1st, filtering_cost_2nd, e_time_kNN, ave * 1000, stdev * 1000, cost_min * 1000, cost_max * 1000);
                #endif
            }
            #endif
            #ifdef SELF_EVAL
            out_result_NN(result_filename, num_queries, correct_answer[m] + q_1st, top_k + q_1st);
            #else
            out_result_NN(result_filename, num_queries, NULL, top_k);
            #endif

//            strcpy(result_filename2, result_filename);
//            result_filename2 = strtok(result_filename2, ".");
//            result_filename2 = strcat(result_filename2, "_result.csv");
//            double recall = out_result_double(result_filename2, PJT_DIM, SMAP_DIM, FTR_DIM, e_time_1st, e_time_score, e_time_2nd, e_time_kNN, e_time(&tp1, &tp2), nc1, nc2,
//                                                num_queries, correct_answer[m], top_k, qr_file[m]);

//			#if NUM_K > 0 || NUM_K == -2
//				for(int q = 0; q < num_queries; q++) { free_kNN_buffer(top_k[q]); }
//				FREE(top_k, sizeof(kNN_buffer *) * num_queries);
//			#endif
			use_system("VmSize");
            if(NUM_Q < num_qr[m]) {
                q_1st += NUM_Q;
                if(q_1st + NUM_Q > num_qr[m]) {
                    q_1st = 0;
                }
            }
		}
        #ifdef SELF_EVAL
        #ifdef FILTERING_BY_SKETCH_ENUMERATION_HAMMING
        printf("filtering, %.4lf, kNN, %.4lf, total, %.4lf, recall, %.1lf, conjunctive enumeration (enum = %d, supp = %d), nt = %d\n", 
                total_filtering / num_query_files, total_kNN / num_query_files, total_total / num_query_files, total_recall / num_query_files, ENUM_DIM, SPP_BIT, 1 << PARALLEL_ENUM);
        #else
        printf("filtering, %.4lf, kNN, %.4lf, total, %.4lf, recall, %.1lf, sketch enumeration c2n, nt = %d\n", 
                total_filtering / num_query_files, total_kNN / num_query_files, total_total / num_query_files, total_recall / num_query_files, 1 << PARALLEL_ENUM);
        #endif
        #else
        #ifdef FILTERING_BY_SKETCH_ENUMERATION_HAMMING
        printf("filtering, %.4lf, kNN, %.4lf, total, %.4lf, conjunctive enumeration (enum = %d, supp = %d), nt = %d\n", 
                total_filtering / num_query_files, total_kNN / num_query_files, total_total / num_query_files, ENUM_DIM, SPP_BIT, 1 << PARALLEL_ENUM);
        #else
        printf("filtering, %.4lf, kNN, %.4lf, total, %.4lf, sketch enumeration c2n, nt = %d\n", 
                total_filtering / num_query_files, total_kNN / num_query_files, total_total / num_query_files, 1 << PARALLEL_ENUM);
        #endif
        #endif
	}

// 後始末
    for(int m = 0; m < num_query_files; m++) {
        free_dataset(ds_query[m]);
        free(qr[m]);
        #ifdef SELF_EVAL
        free(correct_answer[m]);
        #endif
    }

    for(int t = 0; t < dh.num_threads; t++) {
        close_multi_ftr(dh.mf[t]);
    }
	free(dh.mf);

	#if NUM_K > 0 || NUM_K == -2
		for(int q = 0; q < max_num_qr; q++) { free_kNN_buffer(top_k[q]); }
		FREE(top_k, sizeof(kNN_buffer *) * max_num_qr);
	#endif

    free_pivot(pivot);
    free_smap_pivot(smap_pivot);
    free_bucket(bucket_ds);

    if(fp_summary != NULL) {
        fclose(fp_summary);
    }

    if(fp_search_cost != NULL) {
        fclose(fp_search_cost);
    }

	return 0;
}
