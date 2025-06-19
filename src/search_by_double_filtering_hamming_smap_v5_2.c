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
// ともに INTERVAL_WITH_RUN, INTERVAL_WITH_PRIORITY, LOOP_CONTROL_BY_NUM_SKETCHESが定義されていることを前提とする．
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

#ifdef UMPACKED_SMAP
// 非圧縮形式 qpsmap 射影像を用いる. Single-threadのみ
#error "This verson is not applicable for UMPACKED_SMAP."
#elif defined(USE_TABLE_FOR_UMPACKED)
// 表関数は非圧縮用を用いる．Single-threadのみ．qpsmap射影像は圧縮しているが，検索時には umpack して用いる
#error "This verson is not applicable for USE_TABLE_FOR_UMPACKED."
#endif

double filtering_cost_1st, filtering_cost_2nd; 

void reset_filtering_cost(void) {
    filtering_cost_1st = filtering_cost_2nd = 0;
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
	tiny_int packed_qpsmap_data[][PACKED_QPSMAP_SIZE], 
	int data_num[], int num_candidates_1st, int num_candidates_2nd)
{
    struct timespec tp1, tp2, tp3;
    #ifdef USE_INTERVAL
        clock_gettime(CLOCK_METHOD, &tp1);
        // Interval list を用いるときは，
        // ともに INTERVAL_WITH_RUN, INTERVAL_WITH_PRIORITY, LOOP_CONTROL_BY_NUM_SKETCHESが定義されていることを前提とする． 
        // まず，filtering_by_sketch_enumeration_hamming_interval を用いて，フィルタリングを行い，interval_list の形式で候補を求める．
        // 1st filtering の結果を用いて，2nd filtering を qpsmap を用いて行う．
        int nt = 1; // single-thread
    	static interval_list *ivl = NULL;
        if(ivl == NULL) {
            int nnt = nt;
            #ifdef THREAD_PLUS
            nnt *= (1 << THREAD_PLUS);
            #endif
    	    ivl = new_interval_list(nnt, num_candidates_1st);            
        } else if(ivl->size < num_candidates_1st) {
            realloc_interval_list(ivl, num_candidates_1st);
        }
        #ifdef FILTERING_BY_SKETCH_ENUMERATION_HAMMING
        int nc = filtering_by_sketch_enumeration_hamming_interval(qs, bucket, ivl, num_candidates_1st);
        #elif defined(FILTERING_BY_SKETCH_ENUMERATION_C2N)
		static struct_que_c2_n *que = NULL;
        if(que == NULL) {
            que = MALLOC(sizeof(struct_que_c2_n));
            for(int i = 0; i < QSIZE; i += 1024) {
                que->element[i].key = 0;
                que->details[i].sk = 0;
            }
        }
        int nc = filtering_by_sketch_enumeration_c2_n_interval(qs, bucket, que, ivl, num_candidates_1st);
        #else
        #error "FILTERING_BY_SKETCH_ENUMERATION_(HAMMING | C2N) should be defined"
        #endif
        clock_gettime(CLOCK_METHOD, &tp2);

        #ifdef SECOND_FILTERING_KNN_BUFFER
        // 2nd filtering で 1st filtering で求めた nc 個の候補から qpsmap 射影距離で上位（近い方が上位） num_candidates_2nd を選択する．
        kNN_buffer *buff = new_kNN_buffer(num_candidates_2nd);
        #elif defined(SECOND_FILTERING_SELECT)
        // 1st filtering で求めた nc 個の候補のデータ番号と qpsmap 射影距離の対を配列に格納し，最後に quick_select_k_answer で num_candidates_2nd 個を選択する．
        answer_type *ans_buff = MALLOC(sizeof(answer_type) * nc * 2);
        #endif

        int k = 0;
        for(int t = 0; t < ivl->nt; t++) {
            interval *list = ivl->list + t * ivl->size;
            int lg = ivl->lg[t];
            for(int n = 0; n < lg /* && k < num_candidates_1st */; n++) {
                int j = list[n].start;
                for(int r = 0; r < list[n].run /* && k < num_candidates_1st */; j++, r++, k++) {
                    // ここで，スケッチ順に qpsmap が並んでいるので，質問と j 番目の qpsmap との部分復元射影距離を求める．
                    dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_data[j]);
                    answer_type ans = (answer_type){ bucket->idx[j], p_dist };
                    #ifdef SECOND_FILTERING_KNN_BUFFER
                        push_kNN_buffer(&ans, buff); // 元の順番でのデータ番号と射影距離の対を kNN_buffer に push
                    #else
                        ans_buff[k] = ans;
                    #endif
                }
            }
        }

        #ifdef SECOND_FILTERING_KNN_BUFFER
            flush_kNN_buffer(buff);
            for(int i = 0; i < num_candidates_2nd; i++) {
                data_num[i] = buff->buff[i].data_num;
            }
            free_kNN_buffer(buff);
        #elif defined(SECOND_FILTERING_SELECT)
            quick_select_k_answer(ans_buff, 0, k - 1, num_candidates_2nd);
            for(int i = 0; i < num_candidates_2nd; i++) {
                data_num[i] = ans_buff[i].data_num;
            }
            FREE(ans_buff, sizeof(answer_type) * nc * 2);
        #endif

        clock_gettime(CLOCK_METHOD, &tp3);
        filtering_cost_1st += e_time(&tp1, &tp2);
        filtering_cost_2nd += e_time(&tp2, &tp3);
        
        #ifdef SECOND_FILTERING
        return num_candidates_2nd;
        #elif defined(SECOND_FILTERING_SELECT)
        return num_candidates_2nd;
        #else
        return k;
        #endif

    #elif defined(FILTERING_BY_SKETCH_ENUMERATION_HAMMING) // WITHOUT_USING_INTERVAL (single-thread)
        static sub_dimension *table = NULL;
        static int table_size = 0;
        #ifdef ENUM_DIM
        int enum_dim = ENUM_DIM;
        #else
        int enum_dim = PJT_DIM - 19; // ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
        #endif
        if(table == NULL) {
            table_size = (1 << enum_dim) + 1;
            table = make_table_for_enumeration_hamming(table_size, enum_dim);
        }

        #define TABLE_SPP
        // tableを用いた列挙では不足するときに，追加のビットパターンを求めるための表
        #ifdef TABLE_SPP
        static sub_dimension *table_spp = NULL;
        #ifdef SPP_BIT
        int spp_bit = SPP_BIT; // 追加のビット数
        #else
        int spp_bit = 16; // 追加のビット数
        #endif
        static int table_spp_size = 0;
        if(table_spp == NULL) {
            table_spp_size = (1 << spp_bit) + 1;
            table_spp = make_table_for_enumeration_hamming(table_spp_size, spp_bit);
            fprintf(stderr, "make table for spplementary enumeation: spp_bit = %d, table_spp_size = %d\n", spp_bit, table_spp_size);
        }
        #endif

        static int first = 1;
        if(first) {
            #ifdef SECOND_FILTERING_SELECT
            fprintf(stderr, "double-filtering by enum_hamm and qpsmap (single-thread). using table, quick_select, enum_dim = %d, spp_bit = %d\n", enum_dim, spp_bit);
            #elif defined(SECOND_FILTERING_KNN_BUFFER)
            fprintf(stderr, "double-filtering by enum_hamm and qpsmap (single-thread). using table, kNN_buffer, enum_dim = %d, spp_bit = %d\n", enum_dim, spp_bit);
            #elif defined(SKETCH_ENUMERATION_ONLY)
            fprintf(stderr, "double-filtering by enum_hamm and qpsmap (single-thread). using table, enumeration only, enum_dim = %d, spp_bit = %d\n", enum_dim, spp_bit);
            #else
            fprintf(stderr, "double-filtering by enum_hamm and qpsmap (single-thread). using table, 1st only, enum_dim = %d, spp_bit = %d\n", enum_dim, spp_bit);
            #endif
            first = 0;
        }

        int *bd_idx = qs->idx;
        int *bkt = bucket->bkt;

        int k = 0; // 列挙したスケッチから得られたデータ番号の個数
        int n; // 列挙したスケッチの個数
        sketch_type sk, base_mask = 0, base_mask2 = 0;
        #ifdef TABLE_SPP
        int n_spp = 1; // つぎに使用する追加パターン番号（0 のものは空集合なので，1から使用する）
        int n_spp2 = 1; // 追加パターンも使い切ったら，さらに追加する．
        #endif

        #ifdef SECOND_FILTERING_KNN_BUFFER
        // 2nd filtering で 1st filtering で求めた num_candidates_1st 個の候補から qpsmap 射影距離で上位（近い方が上位） num_candidates_2nd を選択する．
        kNN_buffer *buff = new_kNN_buffer(num_candidates_2nd);
        #elif defined(SECOND_FILTERING_SELECT)
        // 1st filtering で求めた num_candidates_1st 個の候補のデータ番号と qpsmap 射影距離の対を配列に格納し，最後に quick_select_k_answer で num_candidates_2nd 個を選択する．
        answer_type *ans_buff = MALLOC(sizeof(answer_type) * num_candidates_1st);
        #endif

        // table を使って mask を作って，つぎのスケッチを列挙する．n は，最初のパターンで列挙されるパターン番号
        for(n = 0; k < num_candidates_1st; n++) {
            sketch_type mask = base_mask;
            for(int m = 0; m < table[n].num; m++) {
                mask |= (1 << bd_idx[(int)table[n].dim[m]]);
            }
            sk = qs->sketch ^ mask;
            for(int j = bkt[sk]; j < bkt[sk + 1] && k < num_candidates_1st; j++, k++) {
                #if defined(SECOND_FILTERING_KNN_BUFFER) || defined(SECOND_FILTERING_SELECT)
                    // ここで，スケッチ順に qpsmap が並んでいるときは，質問と j 番目の qpsmap との部分復元射影距離を求めて
                    dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_data[j]);
                    answer_type ans = (answer_type){ bucket->idx[j], p_dist };
                    #ifdef SECOND_FILTERING_KNN_BUFFER
                    push_kNN_buffer(&ans, buff); // 元の順番でのデータ番号と射影距離の対を kNN_buffer に push
                    #else
                    ans_buff[k] = ans;
                    #endif
                #elif defined(SKETCH_ENUMERATION_ONLY)
                    // 何もしない
                #else
                    data_num[k] = bucket->idx[j];
                #endif
            }
            if(table[n + 1].num == 0) { // 用意したパターンがなくなった．
                #ifdef TABLE_SPP
                if(table_spp[n_spp].num == 0) { // 追加分のパターンも使い切った．
                    if(table_spp[n_spp2].num == 0) break; // 2回目の追加分のパターンも使い切った．
                    n_spp = 0;
                    n_spp2++;
                    base_mask2 = 0;
                    for(int m = 0; m < table_spp[n_spp2].num; m++) {
                        base_mask2 |= (1 << bd_idx[(int)table_spp[n_spp2].dim[m] + enum_dim + spp_bit]);
                    }
                }
                base_mask = base_mask2;
                for(int m = 0; m < table_spp[n_spp].num; m++) {
                    base_mask |= (1 << bd_idx[(int)table_spp[n_spp].dim[m] + enum_dim]);
                }
                n = -1; // for の再初期化で n++ となって，0 になる．
                n_spp++;
                #else
                break;
                #endif
            }
        }
        #ifdef SECOND_FILTERING_KNN_BUFFER
        flush_kNN_buffer(buff);
        for(int i = 0; i < num_candidates_2nd; i++) {
            data_num[i] = buff->buff[i].data_num;
        }
        free_kNN_buffer(buff);
        #elif defined(SECOND_FILTERING_SELECT)
        quick_select_k_answer(ans_buff, 0, num_candidates_1st - 1, num_candidates_2nd);
        for(int i = 0; i < num_candidates_2nd; i++) {
            data_num[i] = ans_buff[i].data_num;
        }
        FREE(ans_buff, sizeof(answer_type) * num_candidates_1st);
        #endif

        #ifdef ENUM_CHECK
        if(table[n].num > max_H) {
            max_H = table[n].num;
            fprintf(stderr, "max Hamming = %d\n", max_H);
        }
        if(n > max_n) {
            max_n = n;
            fprintf(stderr, "max n (number of used table entries) = %d\n", max_n);
        }
        #endif
        #ifdef SECOND_FILTERING_KNN_BUFFER
        return num_candidates_2nd;
        #elif defined(SECOND_FILTERING_SELECT)
        return num_candidates_2nd;
        #else
        return k;
        #endif
    #else // FILTERING_BY_SKETCH_ENUMERATION_C2N (WITHOUT_USING_INTERVAL) (single-thread)
        sketch_type s;
        QUE_c2 qu, qu2;
        int *bd = qs->bd, *bd_idx = qs->idx;
        int *bkt = bucket->bkt;
		static struct_que_c2_n *que = NULL;
        if(que == NULL) que = MALLOC(sizeof(struct_que_c2_n));

        #ifdef SECOND_FILTERING_KNN_BUFFER
        // 2nd filtering で 1st filtering で求めた num_candidates_1st 個の候補から qpsmap 射影距離で上位（近い方が上位） num_candidates_2nd を選択する．
        kNN_buffer *buff = new_kNN_buffer(num_candidates_2nd);
        #elif defined(SECOND_FILTERING_SELECT)
        // 1st filtering で求めた num_candidates_1st 個の候補のデータ番号と qpsmap 射影距離の対を配列に格納し，最後に quick_select_k_answer で num_candidates_2nd 個を選択する．
        answer_type *ans_buff = MALLOC(sizeof(answer_type) * num_candidates_1st);
        #endif

        s = qs->sketch; // 先頭は質問のスケッチ
        int k = 0;

        for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates_1st; j++, k++) {
            #if defined(SECOND_FILTERING_KNN_BUFFER) || defined(SECOND_FILTERING_SELECT)
                // ここで，スケッチ順に qpsmap が並んでいるときは，質問と j 番目の qpsmap との部分復元射影距離を求めて
                dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_data[j]);
                answer_type ans = (answer_type){ bucket->idx[j], p_dist };
                #ifdef SECOND_FILTERING_KNN_BUFFER
                push_kNN_buffer(&ans, buff); // 元の順番でのデータ番号と射影距離の対を kNN_buffer に push
                #else
                ans_buff[k] = ans;
                #endif
            #elif defined(SKETCH_ENUMERATION_ONLY)
                // 何もしない
            #else
                data_num[k] = bucket->idx[j];
            #endif
        }

        s = s ^ (1 <<  bd_idx[0]); // 先頭の次は、質問のスケッチと距離下限が最小のビットだけが異なるもの

        for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates_1st; j++, k++) {
            #if defined(SECOND_FILTERING_KNN_BUFFER) || defined(SECOND_FILTERING_SELECT)
                // ここで，スケッチ順に qpsmap が並んでいるときは，質問と j 番目の qpsmap との部分復元射影距離を求めて
                dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_data[j]);
                answer_type ans = (answer_type){ bucket->idx[j], p_dist };
                #ifdef SECOND_FILTERING_KNN_BUFFER
                push_kNN_buffer(&ans, buff); // 元の順番でのデータ番号と射影距離の対を kNN_buffer に push
                #else
                ans_buff[k] = ans;
                #endif
            #elif defined(SKETCH_ENUMERATION_ONLY)
                // 何もしない
            #else
                data_num[k] = bucket->idx[j];
            #endif
        }

        make_empty_que_c2_n(que);

        // enq pattern of 0...10
        qu.cursor = new_que_e2_n(que);
        qu.key = bd[bd_idx[1]];
        que->details[qu.cursor].sk = qs->sketch ^ (1 << bd_idx[1]);
        que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
        enq_c2_n(&qu, que);		

        while(deq_c2_n(&qu, que) && k < num_candidates_1st) {
            s = que->details[qu.cursor].sk; // 列挙のつぎのスケッチ
            for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates_1st; j++, k++) {
                #if defined(SECOND_FILTERING_KNN_BUFFER) || defined(SECOND_FILTERING_SELECT)
                    // ここで，スケッチ順に qpsmap が並んでいるときは，質問と j 番目の qpsmap との部分復元射影距離を求めて
                    dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_data[j]);
                    answer_type ans = (answer_type){ bucket->idx[j], p_dist };
                    #ifdef SECOND_FILTERING_KNN_BUFFER
                    push_kNN_buffer(&ans, buff); // 元の順番でのデータ番号と射影距離の対を kNN_buffer に push
                    #else
                    ans_buff[k] = ans;
                    #endif
                #elif defined(SKETCH_ENUMERATION_ONLY)
                    // 何もしない
                #else
                    data_num[k] = bucket->idx[j];
                #endif
            }

            switch(que->details[qu.cursor].pt & 15) {
            case 0: // X0000 -> enq(X0001) and enq(Y10^{m+1}) if X0000 = Y010^m
            case 8: // X1000 -> enq(X1001) and enq(Y10^{m+1}) if X0000 = Y010^m
                {
                    int m = lsb_pos(que->details[qu.cursor].pt);
                    if(m > 0 && m < PJT_DIM - 1 && !(que->details[qu.cursor].pt & (1 << (m + 1)))) {
                        // Y010^m -> Y10^{m+1}
                        qu2.cursor = new_que_e2_n(que);
                        qu2.key = qu.key + bd[bd_idx[m + 1]] - bd[bd_idx[m]];
                        que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[m + 1])) ^ (1 << bd_idx[m]);
                        que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
                        // Y010^m -> Y010^{m-1}1
                        qu.key = qu.key + bd[bd_idx[0]];
                        que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0]);
                        que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
                        enq_c2_n(&qu, que);
                        enq_c2_n(&qu2, que);
                    } else {
                        qu.key = qu.key + bd[bd_idx[0]];
                        que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0]);
                        que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
                        enq_c2_n(&qu, que);
                    }
                }
                break;
            case 4:  // X0100 -> enq(X0101) and enq(X1000)
                // X1000
                qu2.cursor = new_que_e2_n(que);
                qu2.key = qu.key + bd[bd_idx[3]] - bd[bd_idx[2]];
                que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[3])) ^ (1 << bd_idx[2]);
                que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
                // X0101
                qu.key = qu.key + bd[bd_idx[0]];
                que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[0]);
                que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
                enq_c2_n(&qu, que);
                enq_c2_n(&qu2, que);
                break;
            case 1:  // X0001 -> enq(X0010)
            case 5:  // X0101 -> enq(X0110)
            case 9:  // X1001 -> enq(X1010)
            case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
                qu.key = qu.key + bd[bd_idx[1]] - bd[bd_idx[0]];
                que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[1])) ^ (1 << bd_idx[0]);
                que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
                enq_c2_n(&qu, que);
                break;
            case 2:  // X0010 -> enq(X0011) and enq(X0100)
            case 10: // X1010 -> enq(X1011) and enq(X1100)
                // X0100 and X1100
                qu2.cursor = new_que_e2_n(que);
                qu2.key = qu.key +  bd[bd_idx[2]] -  bd[bd_idx[1]];
                que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[2])) ^ (1 << bd_idx[1]);
                que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
                // X0011 and X1011
                qu.key = qu.key + bd[bd_idx[0]];
                que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[0]);
                que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
                enq_c2_n(&qu, que);
                enq_c2_n(&qu2, que);
                break;
            case 6:  // X0110 -> enq(X0111)
            case 12: // X1100 -> enq(X1101)
            case 14: // X1110 -> enq(10111)
                qu.key = qu.key + bd[bd_idx[0]];
                que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[0]);
                que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
                enq_c2_n(&qu, que);
                break;
            case 3:  // X0011
            case 7:  // X0111
            case 11: // X1011
            case 15: // X1111 -> nothing to do
                break;
            }
        }
        #ifdef SECOND_FILTERING_KNN_BUFFER
        flush_kNN_buffer(buff);
        for(int i = 0; i < num_candidates_2nd; i++) {
            data_num[i] = buff->buff[i].data_num;
        }
        free_kNN_buffer(buff);
        #elif defined(SECOND_FILTERING_SELECT)
        quick_select_k_answer(ans_buff, 0, num_candidates_1st - 1, num_candidates_2nd);
        for(int i = 0; i < num_candidates_2nd; i++) {
            data_num[i] = ans_buff[i].data_num;
        }
        FREE(ans_buff, sizeof(answer_type) * num_candidates_1st);
        #endif
        #ifdef SECOND_FILTERING_KNN_BUFFER
        return num_candidates_2nd;
        #elif defined(SECOND_FILTERING_SELECT)
        return num_candidates_2nd;
        #else
        return k;
        #endif
    #endif
}
#else

#ifdef STATIC_DF_WORK
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
    que = MALLOC(sizeof(struct_que_c2_n));
    for(int i = 0; i < QSIZE; i += 1024) {
        que->element[i].key = 0;
        que->details[i].sk = 0;
    }
}
#endif
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
	tiny_int packed_qpsmap_data[][PACKED_QPSMAP_SIZE], 
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
    struct timespec tp1, tp2, tp3;

    #ifdef USE_INTERVAL
        clock_gettime(CLOCK_METHOD, &tp1);
        // Interval list を用いるときは，
        // ともに INTERVAL_WITH_RUN, INTERVAL_WITH_PRIORITY, LOOP_CONTROL_BY_NUM_SKETCHESが定義されていることを前提とする． 
        // まず，filtering_by_sketch_enumeration_hamming_interval を用いて，フィルタリングを行い，interval_list の形式で候補を求める．
        // 1st filtering の結果を用いて，2nd filtering を qpsmap を用いて行う．
        #ifndef STATIC_DF_WORK
        static interval_list *ivl = NULL;
        if(ivl == NULL) {
    	    ivl = new_interval_list(nnt, num_candidates_1st);
        } else if(ivl->size < num_candidates_1st) {
            realloc_interval_list(ivl, num_candidates_1st);
        }
        #else
        if(num_candidates_1st > intial_ivl_size_1st) {
            fprintf(stderr, "too large nc1 = %d > %d\n", num_candidates_1st, intial_ivl_size_1st);
            exit(1);
        }
        ivl->size = num_candidates_1st;
        #endif
        static int first = 1;
        #ifdef FILTERING_BY_SKETCH_ENUMERATION_HAMMING
            if(first) {fprintf(stderr, "FILTERING_BY_SKETCH_ENUMERATION_HAMMING, USE_INTERVAL\n"); first = 0;}
            int nc = filtering_by_sketch_enumeration_hamming_interval(qs, bucket, ivl, num_candidates_1st);
            if(nc == 0) {fprintf(stderr, "nc = %d, num_candidates_1st = %d\n", nc, num_candidates_1st); getchar(); }
        #elif defined(FILTERING_BY_SKETCH_ENUMERATION_C2N)
            if(first) {fprintf(stderr, "FILTERING_BY_SKETCH_ENUMERATION_C2N, USE_INTERVAL\n"); first = 0;}
            #ifndef STATIC_DF_WORK
		    static struct_que_c2_n *que = NULL;
            if(que == NULL) {
                que = MALLOC(sizeof(struct_que_c2_n));
                for(int i = 0; i < QSIZE; i += 1024) {
                    que->element[i].key = 0;
                    que->details[i].sk = 0;
                }
            }
            #endif
            int nc = filtering_by_sketch_enumeration_c2_n_interval(qs, bucket, que, ivl, num_candidates_1st);
        #else
            #error "FILTERING_BY_SKETCH_ENUMERATION_(HAMMING | C2N) should be defined"
        #endif
        clock_gettime(CLOCK_METHOD, &tp2);

        int num_data_1st = nc / nt;	// スレッドが 1st filtering で求めるデータ数
        num_candidates_1st = num_data_1st * nt;     // 1st filtering で求めるデータ数の合計（元の num_candidates_1st がスレッド数で割り切れないときに端数を切り捨てる）

        #ifdef SECOND_FILTERING_KNN_BUFFER
           // 2nd filtering を個々のスレッドで kNN_buffer 法で行う
            int num_data_2nd = num_candidates_2nd; // スレッドに分けても同じ候補数を選ばせる
            #ifndef STATIC_DF_WORK
            answer_type *ans_buff = MALLOC(sizeof(answer_type) * num_data_2nd * nt); // 最後に各スレッドが求めてものを一つにまとめて，quick_selectする．
            kNN_buffer *buff_p[nt]; // kNN_buffer のプール（それぞれのスレッドで用いる）
            for(int t = 0; t < nt; t++) { buff_p[t] = new_kNN_buffer(num_data_2nd); }
            #else
            if(num_data_2nd > intial_ivl_size_2nd) {
                fprintf(stderr, "too large nc1 = %d > %d\n", num_candidates_1st, intial_ivl_size_1st);
                exit(1);
            }
            for(int t = 0; t < nt; t++) {
                make_empty_kNN_buffer(buff_p[t]);
                buff_p[t]->k = num_data_2nd;
            }
            #endif
        #elif defined(SECOND_FILTERING_SELECT)
            // 1st filtering で求めた num_candidates_1st 個の候補のデータ番号と qpsmap 射影距離の対を配列に格納し，最後に quick_select_k_answer で num_candidates_2nd 個を選択する．
//fprintf(stderr, "MALLOC for ans_buff: size = %d\n", sizeof(answer_type) * num_candidates_1st * nt);
            static int allocated_nc1 = 0;
            static answer_type *ans_buff = NULL;
            if(allocated_nc1 < num_candidates_1st) {
                if(ans_buff != NULL) {
                    FREE(ans_buff, sizeof(answer_type) * allocated_nc1 * nt);
                }
                ans_buff = MALLOC(sizeof(answer_type) * num_candidates_1st * nt);
                allocated_nc1 = num_candidates_1st;
            }
//fprintf(stderr, "MALLOC OK\n"); getchar();
            int k_th[nt];
        #endif

        // 2nd filtering ...
        #pragma omp parallel
        {
            int t = omp_get_thread_num(); // スレッド番号
            int m = 0; // 列挙したスケッチ数（パターン番号）
            int k = 0; // 1st filtering で求めたデータ数
            #ifdef SECOND_FILTERING_KNN_BUFFER
                kNN_buffer *buff = buff_p[t];
                answer_type *a_buff = ans_buff + t * num_data_2nd;
            #elif defined(SECOND_FILTERING_SELECT)
                answer_type *buff = ans_buff + t * num_candidates_1st;
            #elif defined(SKETCH_ENUMERATION_ONLY)
                // 何もしない
            #else
                int *buff = data_num + t * num_data_1st;
            #endif

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
                        dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_data[j]);
                        answer_type ans = (answer_type){ bucket->idx[j], p_dist };
                        #ifdef SECOND_FILTERING_KNN_BUFFER
                        push_kNN_buffer(&ans, buff); // 元の順番でのデータ番号と射影距離の対を kNN_buffer に push
                        #elif defined(SKETCH_ENUMERATION_ONLY)
                            // 何もしない
                        #else
                        buff[k] = ans;
                        #endif
                    }
                }
            }

            #ifdef SECOND_FILTERING_KNN_BUFFER
            flush_kNN_buffer(buff);
            for(int i = 0; i < num_data_2nd; i++) {
                if(i < buff->num) {
                    a_buff[i] = buff->buff[i];
                } else {
                    a_buff[i] = (answer_type) {0, INT_MAX};
                }
            }
            #elif defined(SECOND_FILTERING_SELECT)
            k_th[t] = k;
            #endif
        }

    #else
        clock_gettime(CLOCK_METHOD, &tp1);
        #ifdef FILTERING_BY_SKETCH_ENUMERATION_HAMMING
            static sub_dimension *table = NULL;
            static int table_size = 0;
            #ifdef ENUM_DIM
            int enum_dim = ENUM_DIM;
            #else
            int enum_dim = PJT_DIM - 19; // ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
            #endif
            if(table == NULL) {
                table_size = (1 << enum_dim) + nt; // 各スレッド用に最後に空集合を与えるので，nt 個余分に確保
                table = make_table_for_enumeration_hamming(table_size, enum_dim);
                rearrange_table(nt, table, table_size);
            }

            // tableを用いた列挙では不足するときに，追加のビットパターンを求めるための表
            static sub_dimension *table_spp = NULL;
            #ifdef SPP_BIT
            int spp_bit = SPP_BIT; // 追加のビット数
            #else
            int spp_bit = 16; // 追加のビット数
            #endif
            static int table_spp_size = 0;
            if(table_spp == NULL) {
                table_spp_size = (1 << spp_bit) + 1;
                table_spp = make_table_for_enumeration_hamming(table_spp_size, spp_bit);
                fprintf(stderr, "make table for spplementary enumeation: spp_bit = %d, table_spp_size = %d\n", spp_bit, table_spp_size);
            }

            static int first = 1;
            if(first) {
                #ifdef SECOND_FILTERING_SELECT
                fprintf(stderr, "double-filtering by enum_hamm and qpsmap (%d-thread). using table, quick_select, enum_dim = %d, spp_bit = %d\n", nt, enum_dim, spp_bit);
                #elif defined(SECOND_FILTERING_KNN_BUFFER)
                fprintf(stderr, "double-filtering by enum_hamm and qpsmap (%d-thread). using table, kNN_buffer, enum_dim = %d, spp_bit = %d\n", nt, enum_dim, spp_bit);
                #elif defined(SKETCH_ENUMERATION_ONLY)
                fprintf(stderr, "double-filtering by enum_hamm and qpsmap (%d-thread). using table, enumeration only, enum_dim = %d, spp_bit = %d\n", nt, enum_dim, spp_bit);
                #else
                fprintf(stderr, "double-filtering by enum_hamm and qpsmap (%d-thread). using table, 1st only, enum_dim = %d, spp_bit = %d\n", nt, enum_dim, spp_bit);
                #endif
                fprintf(stderr, "QUANTIZE_BIT = %d, TABLE_UNIT = %d\n", QUANTIZE_BIT, TABLE_UNIT);
                first = 0;
            }

            int *bd_idx = qs->idx;
            int *bkt = bucket->bkt;

            int num_data_1st = num_candidates_1st / nt;	// スレッドが 1st filtering で求めるデータ数
            num_candidates_1st = num_data_1st * nt;     // 1st filtering で求めるデータ数の合計（元の num_candidates_1st がスレッド数で割り切れないときに端数を切り捨てる）

            #ifdef SECOND_FILTERING_KNN_BUFFER
                // 2nd filtering を個々のスレッドで kNN_buffer 法で行う
                int num_data_2nd = num_candidates_2nd; // スレッドに分けても同じ候補数を選ばせる
                answer_type *ans_buff = MALLOC(sizeof(answer_type) * num_data_2nd * nt); // 最後に各スレッドが求めてものを一つにまとめて，quick_selectする．
                kNN_buffer *buff_p[nt];
                for(int t = 0; t < nt; t++) { buff_p[t] = new_kNN_buffer(num_data_2nd); }
            #elif defined(SECOND_FILTERING_SELECT)
                // 1st filtering で求めた num_candidates_1st 個の候補のデータ番号と qpsmap 射影距離の対を配列に格納し，最後に quick_select_k_answer で num_candidates_2nd 個を選択する．
                answer_type *ans_buff = MALLOC(sizeof(answer_type) * num_candidates_1st);
            #endif

            #pragma omp parallel
            {
                int t = omp_get_thread_num(); // スレッド番号
                int m = 0; // 列挙したスケッチ数（パターン番号）
                int k = 0; // 1st filtering で求めたデータ数
                sub_dimension *tbl = table + t * table_size / nt; // スレッドが使用する部分集合の列挙パターン配列
                sketch_type base_mask = 0, base_mask2 = 0;
                int n_spp = 1, n_spp2 = 1; // 追加のパターン番号
                #ifdef SECOND_FILTERING_KNN_BUFFER
        //        	kNN_buffer *buff = new_kNN_buffer(num_data_2nd);
                    kNN_buffer *buff = buff_p[t];
                    answer_type *a_buff = ans_buff + t * num_data_2nd;
                #elif defined(SECOND_FILTERING_SELECT)
                    answer_type *buff = ans_buff + t * num_data_1st;
                #elif defined(SKETCH_ENUMERATION_ONLY)
                    // 何もしない
                #else
                    int *buff = data_num + t * num_data_1st;
                #endif
                for(m = 0; k < num_data_1st ; m++) {
                    sketch_type sk, mask = base_mask;
                    for(int j = 0; j < tbl[m].num; j++) {
                        mask |= (1 << bd_idx[(int)tbl[m].dim[j]]);
                    }
                    sk = qs->sketch ^ mask;
                    for(int j = bkt[sk]; j < bkt[sk + 1] && k < num_data_1st; j++, k++) {
                        #if defined(SECOND_FILTERING_KNN_BUFFER) || defined(SECOND_FILTERING_SELECT)
                            // ここで，スケッチ順に qpsmap が並んでいるときは，質問と j 番目の qpsmap との部分復元射影距離を求めて
                            dist_type p_dist = projected_dist_packed_table(table_for_packed, packed_qpsmap_data[j]);
                            answer_type ans = (answer_type){ bucket->idx[j], p_dist };
                            #ifdef SECOND_FILTERING_KNN_BUFFER
                            push_kNN_buffer(&ans, buff); // 元の順番でのデータ番号と射影距離の対を kNN_buffer に push
                            #elif defined(SKETCH_ENUMERATION_ONLY)
                                // 何もしない
                            #else
                            buff[k] = ans;
                            #endif
                        #elif defined(SKETCH_ENUMERATION_ONLY)
                            // 何もしない
                        #else
                            buff[k] = bucket->idx[j];
                        #endif
                    }
                    if(tbl[m + 1].num == 0) { // 用意したパターンがなくなった．
                        if(table_spp[n_spp].num == 0) { // 追加分のパターンも使い切った．
                            if(table_spp[n_spp2].num == 0) break; // 2回目の追加分のパターンも使い切った．
                            n_spp = 0;
                            n_spp2++;
                            base_mask2 = 0;
                            for(int m = 0; m < table_spp[n_spp2].num; m++) {
                                base_mask2 |= (1 << bd_idx[(int)table_spp[n_spp2].dim[m] + enum_dim + spp_bit]);
                            }
                        }
                        base_mask = base_mask2;
                        for(int m = 0; m < table_spp[n_spp].num; m++) {
                            base_mask |= (1 << bd_idx[(int)table_spp[n_spp].dim[m] + enum_dim]);
                        }
                        m = -1; // forの再初期化で m++ となって，0になる．
                        n_spp++;
                    }
                }
                #ifdef SECOND_FILTERING_KNN_BUFFER
                flush_kNN_buffer(buff);
                for(int i = 0; i < num_data_2nd; i++) {
                    a_buff[i] = buff->buff[i];
                }
                #endif
            }

        #else
            #error "FILTERING_BY_SKETCH_ENUMERATION_C2N WITHOUT_INTERVAL not supported for multi-thread"
            return num_candidates_1st;
        #endif
        clock_gettime(CLOCK_METHOD, &tp2);
    #endif

    #if defined(USE_INTERVAL) || defined(FILTERING_BY_SKETCH_ENUMERATION_HAMMING)
    #if defined(SECOND_FILTERING_KNN_BUFFER)
        #ifndef STATIC_DF_WORK
        for(int t = 0; t < nt; t++) { free_kNN_buffer(buff_p[t]); }
        #endif
        quick_select_k_answer(ans_buff, 0, num_data_2nd * nt - 1, num_candidates_2nd);
        for(int i = 0; i < num_candidates_2nd; i++) {
            data_num[i] = ans_buff[i].data_num;
        }
        #ifndef STATIC_DF_WORK
        FREE(ans_buff, sizeof(answer_type) * num_data_2nd * nt);
        #endif
        clock_gettime(CLOCK_METHOD, &tp3);
        filtering_cost_1st += e_time(&tp1, &tp2);
        filtering_cost_2nd += e_time(&tp2, &tp3);
    #elif defined(SECOND_FILTERING_SELECT) && !defined(USE_INTERVAL)
        quick_select_k_answer(ans_buff, 0, num_candidates_1st - 1, num_candidates_2nd);
        for(int i = 0; i < num_candidates_2nd; i++) {
            data_num[i] = ans_buff[i].data_num;
        }
        FREE(ans_buff, sizeof(answer_type) * num_candidates_1st);
    #else
/*
fprintf(stderr, "before quick_select_k_answer_mt: nt = %d, num_candidates_1st = %d, num_candidates_2nd = %d\n", nt, num_candidates_1st, num_candidates_2nd);
for(int t = 0; t < nt; t++) {
    fprintf(stderr, "k_th[%d] = %d", t, k_th[t]);
    answer_type *list = ans_buff + t * num_candidates_1st;
    for(int i = 0; i < 5; i++) {
        fprintf(stderr, ", id = %10d, dist = %5d", list[i].data_num, list[i].dist);
    }
    fprintf(stderr, "\n");
}
*/
        int selected = quick_select_k_answer_mt(nt, k_th, num_candidates_1st, ans_buff, num_candidates_2nd);
/*
fprintf(stderr, "after quick_select_k_answer_mt\n");
int sel = 0;
for(int t = 0; t < nt; t++) {
    fprintf(stderr, "k_th[%d] = %d", t, k_th[t]);
    sel += k_th[t];
    answer_type *list = ans_buff + t * num_candidates_1st;
    for(int i = 0; i < k_th[t]; i++) {
        fprintf(stderr, ", id = %10d, dist = %5d", list[i].data_num, list[i].dist);
    }
    fprintf(stderr, "\n");
}
fprintf(stderr, "selected = %d, sel = %d\n", selected, sel); getchar();
*/
        int i = 0;
        nc = 0;
//fprintf(stderr, "selected = ");
        for(int t = 0; t < nt; t++) {
            answer_type *list = ans_buff + t * num_candidates_1st;
            for(int j = 0; j < k_th[t] && nc < num_candidates_2nd; j++, nc++) {
//fprintf(stderr, "%12d", list[j].data_num);
                data_num[i++] = list[j].data_num;
            }
        }
//getchar();
    #endif
    #endif

    #ifdef SECOND_FILTERING_KNN_BUFFER
    return num_candidates_2nd;
    #elif defined(SECOND_FILTERING_SELECT)
    return num_candidates_2nd;
    #else
    return num_candidates_1st;
    #endif
}



#endif

#define NUM_NN 30
int main(int argc, char *argv[])
{
	int num_ftr_files = argc - 1;
	char **dataset_ftr_filename = argv + 1;
	char *pivot_file = PIVOT_FILE;
	char *smap_pivot_file = SMAP_PIVOT_FILE;
	char *bucket_filename = BUCKET_FILE;
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

    // 2nd filtering で用いる qpsmap のためのパラメタを準備する． 
    double ave[SMAP_DIM], stdev[SMAP_DIM], offset[SMAP_DIM], slice[SMAP_DIM];
    compute_parmeters_for_qpsmap(smap_pivot, ds_query[0], ave, stdev, offset, slice);

    // データを圧縮表現のqpsmapとして求めておく．ただし，1st filtering で用いるスケッチ順にソートしておく．
    tiny_int (*packed_qpsmap_data)[PACKED_QPSMAP_SIZE]; // 圧縮形式の量子化射影像はデータのみ．質問は量子化していない射影像だけを使用する．
    packed_qpsmap_data = MALLOC(sizeof(tiny_int) * PACKED_QPSMAP_SIZE * num_data);
	fprintf(stderr, "malloc packed quantized images of data OK.\n");
	use_system("VmSize");

	int *data_num_candidates = MALLOC(sizeof(int) * num_data); // double-filtering で求めた候補データのデータ番号を格納する配列
	int *data_num_org = MALLOC(sizeof(int) * num_data); // データ番号をそのまま順番に格納する配列（データをファイルから読み込むときに使用する．
	fprintf(stderr, "malloc data_num_candidates and data_num_org OK: num_data = %d.\n", num_data);
	use_system("VmSize");

    #ifdef FTR_ON_MAIN_MEMORY
    struct_ftr_id *ftr_id = MALLOC(sizeof(struct_ftr_id) * num_data);
	fprintf(stderr, "malloc ftr_id OK: num_data = %d.\n", num_data);
	use_system("VmSize");
    #endif

    fprintf(stderr, "make packed quantized images of data ... \n");
    // スケッチ順に特徴データを読み込むと，ランダムアクセス的になって遅くなるかもしれないので，
	// いったん，元の順序のままで読み込んで，圧縮形式のSMAPに変換して，その後で並べ替える．
    for(int i = 0; i < num_data; i++) { data_num_org[i] = i; } // 特徴データを元の順序で読み込むためのデータ番号の配列としても使用する．
    for(int i = 0; i < num_data; i++) { data_num_candidates[i] = i; } // 配列をRAMに置くためのダミーアクセス．

    struct timespec tread1, tread2;
    clock_gettime(CLOCK_METHOD, &tread1);

    #ifdef _OPENMP
    omp_set_num_threads(NUM_THREADS_PREPERATION);
    #pragma omp parallel for
    #endif
    for(int i = 0; i < num_data; i++) {
        // ここは，マルチスレッドで並列化した方がよいかも．
        // ある程度まとめて（シングルスレッドで）連続に読み込んでから，並列処理する．（未実装11/19時点）
        unsigned char quantized_projected_data[SMAP_DIM]; // 非圧縮形式のデータのqpsmap（作業用：つぎつぎに圧縮形式に変換するので，1個分のみで，ループ内の局所変数にする）
        #ifndef _OPENMP
 		struct_ftr_id *ftr_id_p = get_next_ftr_id_from_multi_ftr(mf, data_num_org, i, num_data);
        #else
        int t = omp_get_thread_num();
 		struct_ftr_id *ftr_id_p = get_next_ftr_id_from_multi_ftr(dh.mf[t], data_num_org, i, num_data);
        #endif

        #ifdef FTR_ON_MAIN_MEMORY
        memcpy(&ftr_id[i], ftr_id_p, sizeof(struct_ftr_id));
        #endif

        q_uchar_psmap(quantized_projected_data, ftr_id_p->ftr, smap_pivot, offset, slice);
        #ifdef QUANTIZE_MIXED_MOD3
        char2tiny_mod3(quantized_projected_data, packed_qpsmap_data[i], SMAP_DIM);
        #else
        char2tiny(quantized_projected_data, packed_qpsmap_data[i], SMAP_DIM, QUANTIZE_BIT);
        #endif
        #ifdef _OPENMP
        int nt = omp_get_num_threads();
        int nd = num_data / nt;
        #else
        int t = 0; 
        int nd = num_data;
        #endif
        if(t == 0 && (i + 1) % (nd / 100) == 0) {
            fprintf(stderr, "made packed quantized images of data %4d%%\r", (i + 1) / (nd / 100));
        }
    }
    clock_gettime(CLOCK_METHOD, &tread2);
    fprintf(stderr, "\ndone: %.4lf (sec)\n", e_time(&tread1, &tread2));

// 2段階検索で実際に特徴データを2次記憶から読み込んで処理を行う場合は，ランダムアクセスになるので，ブロック読込みは非効率．
// ブロック読込みをやめるために，block_size = 1 に変更する．（一旦すべてcloseしてから，再度openする方法では，なぜか速度が落ちる）
	for(int t = 0; t < dh.num_threads; t++) {
        FREE(dh.mf[t]->ftr_id, sizeof(struct_ftr_id) * dh.mf[t]->block_size);
        FREE(dh.mf[t]->data_num, sizeof(int) * dh.mf[t]->block_size);
		dh.mf[t]->block_size = 1;
        dh.mf[t]->ftr_id = MALLOC(sizeof(struct_ftr_id));
        dh.mf[t]->data_num = NULL;
	}

	fprintf(stderr, "sort packed quantized images of data in sketch order ... \n");
	tiny_int *temp = MALLOC(sizeof(tiny_int) * PACKED_QPSMAP_SIZE);
	char *done = (char *)calloc(num_data, sizeof(char));
	for(int i = 0; i < num_data; i++) {
		if(done[i] || i == bucket_ds->idx[i]) {
			done[i] = 1;
			if((i + 1) % (num_data / 100) == 0) {
				fprintf(stderr, "sorted packed quantized images of data %4d%%\r", (i + 1) / (num_data / 100));
			}
			continue;
		}
		memcpy(temp, packed_qpsmap_data[i], sizeof(tiny_int) * PACKED_QPSMAP_SIZE);
		int j;
		for(j = i; i != bucket_ds->idx[j]; j = bucket_ds->idx[j]) {
			memcpy(packed_qpsmap_data[j], packed_qpsmap_data[bucket_ds->idx[j]], sizeof(tiny_int) * PACKED_QPSMAP_SIZE); 
			done[j] = 1;
		}
		memcpy(packed_qpsmap_data[j], temp, sizeof(tiny_int) * PACKED_QPSMAP_SIZE);
		done[j] = 1;
        if((i + 1) % (num_data / 100) == 0) {
            fprintf(stderr, "sorted packed quantized images of data %4d%%\r", (i + 1) / (num_data / 100));
        }
	}

    clock_gettime(CLOCK_METHOD, &tread2);
    fprintf(stderr, "\ndone: %.4lf (sec)\n", e_time(&tread1, &tread2));

	FREE(temp, sizeof(tiny_int) * PACKED_QPSMAP_SIZE);
	free(done);

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

    #ifdef STATIC_DF_WORK
    init_space_for_double_filtering(1000000, 2000);
    #endif

    if(fp_summary != NULL) {
        #ifdef SELF_EVAL
        fprintf(fp_summary, "trial, query, width, q_bit, ftr_on, nc1, nc2, recall@1, recall@30, filtering, 1st(sec), 2nd(sec), kNN(sec), ave(ms/q), stdev(ms/q), min(ms/q), max(ms/q)\n");
        #else
        fprintf(fp_summary, "trial, query, width, q_bit, ftr_on, nc1, nc2, filtering, 1st(sec), 2nd(sec), kNN(sec), ave(ms/q), stdev(ms/q), min(ms/q), max(ms/q)\n");
        #endif
    }
    int trial = 0;
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
        if(num_candidates_2nd <= 30) num_candidates_2nd = 30;
        fprintf(stderr, "num_data = %d, num_candidates_1st = %d, num_candidates_2nd = %d\n", num_data, num_candidates_1st, num_candidates_2nd);
        #ifdef SELF_EVAL
        double total_filtering = 0, total_kNN = 0, total_total = 0, total_recall = 0;
        #else
        double total_filtering = 0, total_kNN = 0, total_total = 0;
        #endif

        for(int m = 0; m < num_query_files; m++) {
			use_system("VmSize");
			int num_queries = num_qr[m];
			#ifdef NUM_Q
			if(NUM_Q != 0 && NUM_Q < num_queries) { num_queries = NUM_Q; }
			#endif
        
			#if NUM_K > 0 || NUM_K == -2
				for(int q = 0; q < num_queries; q++) { make_empty_kNN_buffer(top_k[q]); }
			#endif

            struct timespec tp1, tp2, tp3;
            clock_gettime(CLOCK_METHOD, &tp1);
            double e_time_filtering = 0, e_time_kNN = 0, e_time_total = 0;
            reset_filtering_cost();

            double trial_search_cost[num_queries];
            #ifdef SELF_EVAL
			int found = 0;
            #endif

            for(int q = 0; q < num_queries; q++) {

                // 1st filtering のための query_sketch を作る．
	            clock_gettime(CLOCK_METHOD, &tp1);
                set_query_sketch(&query_sketch, &qr[m][q], pivot);

                // 2nd filtering のためのデータの圧縮形式の qpsmap との射影距離を計算するための表関数 table_for_packed を作る．
                psmap(q_smap, qr[m][q].ftr, smap_pivot); // 質問（ftr）の smap 射影像 q_smap を求める．
                make_table_for_query_p(q_smap, table, offset, slice, 0.1 * SCORE_P_2ND); // 射影像を用いて，射影距離計算の表 table を作成する．
                make_table_for_packed_data(table, table_for_packed, offset, slice); // 圧縮表現データのための表 table_for_packed を作成する．

                // double-filtering でデータ番号を求める．
                double_filtering_by_sketch_enumeration_hamming_and_qpsmap(&query_sketch, table_for_packed, bucket_ds, packed_qpsmap_data, data_num_candidates, num_candidates_1st, num_candidates_2nd);

                clock_gettime(CLOCK_METHOD, &tp2);
				e_time_filtering += e_time(&tp1, &tp2);

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
                    }
                #endif
				clock_gettime(CLOCK_METHOD, &tp3);
				e_time_kNN += e_time(&tp2, &tp3);
				e_time_total += e_time(&tp1, &tp3);
                trial_search_cost[q] = e_time(&tp1, &tp3);
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
            fprintf(fp_search_cost, "trial, file, query_num, search_cost\n");
            for(int q = 0; q < num_queries; q++) {
                fprintf(fp_search_cost, "%d, %d, %d, %.4lf\n", trial, m, q, trial_search_cost[q]);
            }
            fprintf(fp_search_cost, "%d, %d, summary, %.4lf, %.4lf, %.4lf, %.4lf\n", trial, m, ave * 1000, stdev * 1000, cost_min * 1000, cost_max * 1000);
            #endif

//            #ifdef PRINT_KNN_RESULT
//            for(int q = 0; q < num_queries; q++) {
//                printf("%5d, %1d, %5d, %5d, %5d, %10.3le, %10.3le\n", trial, m, q, trial_found[q], trial_dist[q], trial_filtering_cost[q], trial_search_cost[q]);
//            }
//            #endif

            #ifdef SELF_EVAL
            double recall_1 = recall_kNN_1(num_queries, correct_answer[m], top_k);
            double recall_k = recall_kNN(num_queries, correct_answer[m], top_k);
            printf("filtering, %.4lf, kNN, %.4lf, total, %.4lf, ave = %.4lf (ms/q), stdev = %.4lf (ms/q), recall_1, %.1lf, found = %d, recall_k, %.1lf\n", 
                e_time_filtering, e_time_kNN, e_time_total, ave * 1000, stdev * 1000, recall_1, found, recall_k);
            printf("filtering cost: 1st = %.4lf, 2nd = %.4lf\n", filtering_cost_1st, filtering_cost_2nd);
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
            printf("filtering cost: 1st = %.4lf, 2nd = %.4lf\n", filtering_cost_1st, filtering_cost_2nd);
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
            out_result_NN(result_filename, num_queries, correct_answer[m], top_k);
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
