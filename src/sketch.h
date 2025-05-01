#pragma once
// スケッチの型定義と基本操作およびスケッチを用いるkNN検索に関するもののみ
// 距離関数は特徴データに関するものをまとめた ftr に置く
// ピボット選択（QBPやAIR, LS）に関するものは pivot_selection にまとめる

// 前提
// 定数マクロ
// FTR_DIM = 特徴データの次元
// PJT_DIM = 射影次元．ここではピボットの幅（ビット数）
// PIVOT_PARTITION = GHP | BP | QBP | PQBP （基礎分割関数）（とりあえずQBPのみなので無視）

#include <stdio.h>
#include <stdlib.h>
#include "parm.h"
#include "config.h"
#include "ftr.h"
#include <limits.h>

#ifndef PJT_DIM
#define PJT_DIM 16
#define NARROW_SKETCH
#endif
// SKETCH_SIZE = スケッチの大きさ（NARROW_SKETCH -> 32ビット単位, WIDE_SKETCH -> 64ビット単位, EXPANDED_SKETCH -> 64ビット単位）
// TABLE_SIZE = score計算のための表関数の表の数（8ビット毎に256の大きさの表）
// sketch_type:	スケッチの型（32ビットまで）
#if defined(NARROW_SKETCH)
	#define SKETCH_SIZE 1
	#define TABLE_SIZE 4
	typedef unsigned int sketch_type;
#elif defined(WIDE_SKETCH)
	#define SKETCH_SIZE 1
	#define TABLE_SIZE 8
	typedef unsigned long sketch_type;
#elif defined(EXPANDED_SKETCH)
	#define SKETCH_SIZE ((PJT_DIM + 63) / 64)
	#define TABLE_SIZE (SKETCH_SIZE * 8)
	typedef unsigned long sketch_type[SKETCH_SIZE];
#endif

/*
#ifdef PARTITION_TYPE_PQBP // 座標分割QBP（PQBP）のときは，FTR_DIMをPJT_DIM個に分割する．
// FTR_DIMがPJT_DIMで割り切れないときは，最初の部分空間に1個ずつ座標を加える．
// （例）FTR_DIM = 64, PJT_DIM = 6 のときは，10次元ずつで4次元分の余りがでる．
//  最初の4空間（0から3）までは，11次元，残りは10次元の部分空間とする．
//  0: 0 - 10, 1: 11 - 21, 2: 22 - 32, 3: 33 - 43, 4: 44 - 53, 6: 54 - 63  
#define PART_START(j) ((FTR_DIM / PJT_DIM) * j + (j < FTR_DIM % PJT_DIM ? j : 0)) // 部分座標空間の開始次元
#define PART_DIM(j) ((FTR_DIM / PJT_DIM) + (j < FTR_DIM % PJT_DIM ? 1 : 0)) // 部分座標空間の次元数
#endif
*/

#ifndef NUM_PART
#define NUM_PART 1
#endif
#define PART_PJT_MIN (PJT_DIM / NUM_PART)
#define PART_PJT_MAX ((PJT_DIM + NUM_PART - 1) / NUM_PART)

// PART_NUM(p) = 射影次元（p = 0, ... , PJT_DIM - 1）に対応する部分空間番号
#define PART_NUM(p) (((p) / PART_PJT_MAX) < (PJT_DIM % NUM_PART) ? ((p) / PART_PJT_MAX) : ((p) - PJT_DIM % NUM_PART) / PART_PJT_MIN)
// PART_START(p), PART_DIM(p), PART_END(p), PART_PJT_DIM(p) 射影次元（p = 0, ... , PJT_DIM - 1）に対応する部分空間の
// 開始次元番号, 次元数, 最終次元番号, 射影次元数
#define PART_START(j) ((FTR_DIM / NUM_PART) * PART_NUM(j) + (PART_NUM(j) < FTR_DIM % NUM_PART ? PART_NUM(j) : FTR_DIM % NUM_PART))
#define PART_DIM(j) ((FTR_DIM / NUM_PART) + (PART_NUM(j) < FTR_DIM % NUM_PART ? 1 : 0))
#define PART_END(j) (PART_START(j) + PART_DIM(j) - 1)
#define PART_PJT_DIM(j) ((PJT_DIM / NUM_PART) + (PART_NUM(j) < PJT_DIM % NUM_PART ? 1 : 0))

/*
p_dim   part_num        part_start      part_end        part_dim        part_pjt_dim
0 - 3   0               0               10              11              4
4 - 7   1               11              21              11              4
8 - 10  2               22              32              11              3
11 - 13 3               33              43              11              3
14 - 16 4               44              53              10              3
17 - 19 5               54              63              10              3
*/

typedef struct {
	ftr_type  p[PJT_DIM];	// BP, QBP の中心点
	dist_type r[PJT_DIM];	// BP, QBP の半径
	#ifdef PARTITION_TYPE_SQBP
	int used[PJT_DIM][FTR_DIM], num_used[PJT_DIM];
	#endif
	#ifdef PARTITION_TYPE_CSQBP
	int j_start[PJT_DIM], num_j[PJT_DIM]; // 使用する最初の座標軸番号 (0, ... , PJT_DIM - 1)，使用する座標軸数
	#endif
	#ifdef USE_PD_SKETCH
	int axis[PJT_DIM][FTR_DIM];
	int num_axis[PJT_DIM];
	#endif
	partition_type type;				// 基礎分割関数（GHP = 0, BP = 1, QBP = 2, PQBP = 3, SQBP = 4）（定数は enum で定義されている）
} pivot_type;

// query_num: 質問番号, ftr: 特徴データ（ポインタ）
typedef struct {
	int query_num;			// 質問番号（質問のftrファイル内での通し番号）
	ftr_type ftr;			// 特徴データ
	// auxi_type query_id;		// 質問ID（ftrの補助情報）（処理が重たくなるので削除）
} query_type;

// data_num: データ番号, dist: 質問との距離（フィルタリングでは，score_p）
#ifndef NUM_NN
#define NUM_NN 30
#endif
//	int data_num[NUM_NN];			// データ番号（ftrファイル内での通し番号）
//	dist_type dist[NUM_NN];			// 質問との距離（フィルタリングでは，score_p）
typedef struct {
	int data_num[NUM_NN];			// データ番号（ftrファイル内での通し番号）
	#ifndef ANSWER_DIST_FLOAT
	dist_type dist[NUM_NN];			// 質問との距離（フィルタリングでは，score_p）
	#else
	float dist[NUM_NN];
	#endif
} answer_type_NN;

typedef struct {
	int data_num;			// データ番号（ftrファイル内での通し番号）
	dist_type dist;			// 質問との距離（フィルタリングでは，score_p）
	#ifdef ANSWER_WITH_DATA_ID
	auxi_type data_id;		// データID（ftrの補助情報）（処理が重たくなるので削除）
	#endif
} answer_type;

// スケッチと優先度，データ数の構造体
typedef struct {
	sketch_type sk;			// スケッチ
	dist_type priority;		// 優先度（D~1）
	int num;				// スケッチを持つデータ数
} sketch_with_priority_num;

#ifdef CHECK_ENUMERATION
typedef struct {
	int num_enum_sketches;			// 正解が見つかるまでに列挙したすべてのスケッチ数
	int num_enum_nonempty_sketches;	// 正解が見つかるまでに列挙した空でないスケッチ数
	int num_enum_data;				// 正解が見つかるまでに列挙したデータ数
	int num_data[10];				// 列挙した最初の10個のスケッチまでのデータ数（累積）
} struct_check_result;
#endif

#ifndef USE_COMPACT_INTERVAL
typedef struct {
	#ifdef INTERVAL_WITH_PRIORITY
	dist_type priority;
	#endif
	int start;
	#ifdef INTERVAL_WITH_RUN
	int run;
	#else
	int end;
	#endif
} interval;
#else // COMPACT_INTERVAL
typedef struct {
	int start;
	unsigned short run;
	#ifdef INTERVAL_WITH_PRIORITY
	unsigned short priority;
	#endif
} interval;
#endif

//	int nt;		// nt 個のスレッドで分割
//	size;		// 各スレッドが用いる interval 配列 list の大きさ（inteval 数）
//	int *lg;	// lg[t] でスレッド t が分担する interval 数
//	interval *list;		// スレッド t が用いる部分: list[t * size], ... , list[(t + 1) * size - 1], ...
typedef struct {
	int nt;	// nt 個のスレッドで分割
	int size;		// 各スレッドが用いる interval 配列 list の大きさ（inteval 数）
	int *lg;	// num[t] でスレッド t が分担する interval 数
	interval *list;		// スレッド t が用いる部分: list[t * size], ... , list[(t + 1) * size - 1], ...
} interval_list;

typedef struct {
	int size, step;			// 配列の大きさ，拡張するときの増分
	int num_list;					// interval数
	int num_data;					// データ数（interval長の合計）
	interval *elm;	// interval を格納する配列
} vlist;

//	answer_type *buff;		// 解 (id, dist) を入れるバッファー（大きさは 2k）
//	int k, num;				// バッファーの大きさ k と，格納している解の個数
//	dist_type k_nearest;	// k-NN の距離（ただし，最後にソートして
typedef struct {
	answer_type *buff;		// 解 (id, dist) を入れるバッファー（大きさは 2k）
	int k, num;				// バッファーの大きさ k と，格納している解の個数
	dist_type k_nearest;	// k-NN の距離（ただし，最後にソートして
} kNN_buffer;

//	answer_type *answer;
//	ftr_type ftr;
//	sketch_type sketch;
// typedef struct {
//	answer_type *answer;
//	ftr_type ftr;
//	sketch_type sketch;
// } struct_answer_sketch;

//	query_type query;
//	sketch_type sketch;		// 質問点のスケッチ
//	int bd[PJT_DIM];		// 質問点と分割境界との距離（距離下限やscoreの計算に使用）
//	int idx[PJT_DIM];		// bd の順位表
//	dist_type tbl[TABLE_SIZE][256]; 	// scoreの計算に用いる表関数
//	answer_type answer; 	// 最近傍のデータ番号と距離
//	sketch_type answer_sketch; // 最近傍のスケッチ
//	#ifdef USE_AIR
//	#ifdef DUAL_SAMPLE_DS
//	answer_type answer_2; 	// 最近傍のデータ番号と距離
//	sketch_type answer_sketch_2; // 最近傍のスケッチ
//	#endif
//	int computed;			// AIRにおけるピボット評価のとき再サンプリングで選ばれていてquery_sketchが計算済みかどうか
//	#endif
typedef struct {
	query_type query;
	sketch_type sketch;		// 質問点のスケッチ
	#ifdef PLUS_SD
	int bd[PJT_DIM];		// 質問点と分割境界との距離（距離下限やscoreの計算に使用）
	double bd_plus[PJT_DIM];// bd をsmap射影像の標準偏差によって補正したもの
	int idx[PJT_DIM];		// bd_plus の順位表
	#else
	int bd[PJT_DIM];		// 質問点と分割境界との距離（距離下限やscoreの計算に使用）
	int idx[PJT_DIM];		// bd の順位表
	#endif
	dist_type tbl[TABLE_SIZE][256]; 	// scoreの計算に用いる表関数
	#ifndef NUM_NN
		answer_type answer; 				// 最近傍のデータ番号と距離（AIRではpriority）
		sketch_type answer_sketch;			// 最近傍のスケッチ
	#else
		answer_type_NN answer;			// answer[0] -> 最近傍のデータ番号と距離（AIRではpriority），answer[k] -> k + 1 近傍
		sketch_type answer_sketch[NUM_NN];	// k 近傍のスケッチ
		dist_type max_dist;				// 質問と近傍の射影距離の最大値
		int p_idx[NUM_NN];				// 質問と近傍の射影距離の相対ソート用インデックス
	#endif
	#ifdef USE_AIR
		#ifdef DUAL_SAMPLE_DS
			#ifndef NUM_NN
				answer_type answer_2; 				// 最近傍のデータ番号と距離
				sketch_type answer_sketch_2; 		// 最近傍のスケッチ
			#else
				answer_type_NN answer_2; 			// k-近傍のデータ番号と距離
				sketch_type answer_sketch_2[NUM_NN]; 	// k-近傍のスケッチ
			#endif
		#endif
		#ifdef USING_CORRELATION
			dist_type *distance, *priority;   // サンプルデータとサンプル質問間の距離，優先順位
		#endif
		int computed;			// AIRにおけるピボット評価のとき再サンプリングで選ばれていてquery_sketchが計算済みかどうか
	#endif
} struct_query_sketch;

//	sketch_type sk;		// sketch
//	int num;			// number of data whose sketch is s
//	int pos;			// offset position in sftr
typedef struct {
	sketch_type sk;		// sketch
	int num;			// number of data whose sketch are sk
	int pos;			// offset position in sftr or bucket
} sk_num_pair;

//	int num_data; 				// the number of data points in dataset
//	ftr_type *ftr_data;			// array of ftr data (NULL, if not used)
//	sketch_type *sk;			// array of sketches (NULL, if not used)
//	int *idx; 					// idx[i] = the original position of the i-th data in sketch order (NULL, if not used)
//	#ifdef NARROW_SKETCH
//	int *bkt;					// bkt[s] = The first position of the data whose sketch is s in the data arranged in the sketch order (NULL, if not used)
//	#endif
//	int num_nonempty_buckets;	// number of nonempty buckets = number of sk_num_pairs
//	sk_num_pair *sk_num;		// representing nonempty buckets (NULL, if not used)
typedef struct {
	int num_data; 				// the number of data points in dataset
	ftr_type *ftr_data;			// array of ftr data (NULL, if not used)
	sketch_type *sk;			// array of sketches (NULL, if not used)
	int *idx; 					// idx[i] = the original position of the i-th data in sketch order (NULL, if not used)
	#ifdef NARROW_SKETCH
	int *bkt;					// bkt[s] = The first position of the data whose sketch is s in the data arranged in the sketch order (NULL, if not used)
	#endif
	#ifdef REVERSE_IDX
	int *r_idx;					// r-idx[x] = the position in sketch order of the x-th data in the original dataset.
	#endif
	int num_nonempty_buckets;	// number of nonempty buckets = number of sk_num_pairs
	sk_num_pair *sk_num;		// representing nonempty buckets (NULL, if not used)
	#ifdef WITH_STAT
	int stat[PJT_DIM];			// stat[j] = number of data with 1s (ON-bit)
	#endif
} struct_bucket;

typedef struct {
	char *filename;				// filename of bucket file
	FILE *fp;					// bucket file handler
	int num_data; 				// the number of data points in dataset
	int num_nonempty_buckets;	// number of nonempty buckets = number of sk_num_pairs
	int processed_buckets;		// number of processed buckets (sk_num_pairs)
	sk_num_pair sk_num;			// next pair of sketch and number of records
} struct_bucket_sk_num;

//#ifdef NARROW_SKETCH // Priority queue for sketch enumeration only for NARROW sketches

#if defined(NARROW_SKETCH)

// 射影次元の部分集合を表す．
typedef struct {
	int num;			// 要素数
	char dim[PJT_DIM]; 	// 次元番号
} sub_dimension;

//#define QSIZE  BIT 
#define QSIZE  (1L << PJT_DIM) // 最悪の場合
// #define QSIZE 100000	// 実際には，m 個のスケッチを列挙するためには，queue の最大要素数は m．1個のスケッチに対する平均データ数は 22ビットで　900以上，26ビットでも20個以上
						// 2^w の割り当ては無駄で，減らした方がよいと考えてみたが，速度などへの悪影響はなかった．
typedef struct {
    dist_type key;
	sketch_type sk;
	unsigned idx;
} QUE;

typedef struct {
	QUE element[QSIZE + 1];
	int qsize;
} struct_que;

// カーソル使用（imamura）
typedef struct {
    dist_type key;
    unsigned cursor;
} QUE_c2;

// カーソル使用（imamura）mnodified by Takeshi (based on Naoya's idea)
typedef struct {
	sketch_type sk; // sketch (or sketch ^ query sketch)
	// ここには，スケッチ（または，スケッチと質問スケッチのXOR）
	sketch_type pt; // sketch ^ query sketch in rearranged bit order
	// ここには，スケッチと質問スケッチのXORをビットを距離下限順に並べ替えたもの
	// 最初から見ると，0, 1, 10, 11, 100, ... のようになるはず
} QUE_Detail_n;

typedef struct {
    int qsize;
    int detail_size;
    QUE_c2 element[QSIZE + 1];
    QUE_Detail_n details[QSIZE + 1];
} struct_que_c2_n;
#endif

#ifdef EXPANDED_SKETCH
void data_to_sketch(ftr_type o, pivot_type *pivot, sketch_type sk);
void data_to_sketch_1bit(ftr_type o, pivot_type *pivot, int dim, sketch_type sk);
#else
sketch_type data_to_sketch(ftr_type o, pivot_type *pivot);
void data_to_sketch_1bit(ftr_type o, pivot_type *pivot, int dim, sketch_type *sk);
#endif
pivot_type *new_pivot(int type);
void free_pivot(pivot_type *pivot);
// pivot_type *new_pivot_pjt_dim(int type, int pjt_dim);
void read_pivot(const char *filename, pivot_type *pivot);
// void read_pivot_pjt_dim(char *filename, pivot_type *pivot, int pjt_dim);
void write_pivot(char *filename, pivot_type *pivot);
#ifdef PRE_ROTATION
void read_pivot_with_rotation(char *filename, pivot_type *pivot, int *rotation);
void read_rotation(char *filename, int *rotation);
void write_pivot_with_rotation(char *filename, pivot_type *pivot, int *rotation);
#endif
struct_bucket *new_bucket(int num_data, sketch_type sketch[]);
void free_bucket(struct_bucket *b);
void write_bucket(char *filename, struct_bucket *b);
struct_bucket *read_bucket(char *filename);
struct_bucket *read_compact_bucket(char *filename);

struct_bucket_sk_num *open_bucket_sk_num(char *filename);
int read_next_bucket_sk_num(struct_bucket_sk_num *bsk);

#if defined(NARROW_SKETCH)
void min_heapify_p(int i, struct_que *que);
int deq_p(QUE *q, struct_que *que);
void enq_p(QUE *q, struct_que *que);

void make_empty_que_c2_n(struct_que_c2_n *que);
int new_que_e2_n(struct_que_c2_n *que);
void min_heapify_c2_n(int i, struct_que_c2_n *que);
int deq_c2_n(QUE_c2 *qe, struct_que_c2_n *que);
void deq_c2_n_del(struct_que_c2_n *que);
void enq_c2_n(QUE_c2 *qe, struct_que_c2_n *que);
void enq_c2_n_after_deq(QUE_c2 *qe, struct_que_c2_n *que);
#endif

// struct_query_sketch *make_query_sketch(struct_dataset *ds_query, pivot_type *pivot);

// query_type *new_query(void);

void set_query_sketch(struct_query_sketch *qs, query_type *query, pivot_type *pivot);
void compute_sketch_and_boundary_plus(dist_type bd_plus[][PJT_DIM], int num_queries, struct_query_sketch *qs_all, query_type *query_all, pivot_type *pivot);
void set_query_sketch_p_boundary_plus(dist_type bd_plus[PJT_DIM], struct_query_sketch *qs, query_type *query, pivot_type *pivot, double p);
void set_query_sketch_p(struct_query_sketch *qs, query_type *query, pivot_type *pivot, double p);
#ifdef PLUS_SD
void compute_query_sketch_and_ave_stdev_bd_of_pjt_dim(double ave[], double stdev[], int num_queries, struct_query_sketch query_sketch[], query_type query[], pivot_type *pivot);
void compute_query_sketch_and_ave0_ave1_of_pjt_dim(double ave0[], double ave1[], int num_queries, struct_query_sketch query_sketch[], query_type query[], pivot_type *pivot);
void set_query_sketch_p_plus_sd(double ave[], double stdev[], struct_query_sketch *query_sketch, double p);
void set_query_sketch_p_ave0_ave1(double ave0[], double ave1[], struct_query_sketch *query_sketch, double p);
#endif
dist_type priority(sketch_type s, struct_query_sketch *qs);
dist_type priority_inf(sketch_type s, struct_query_sketch *qs);
dist_type priority_partitioned(sketch_type s, struct_query_sketch *qs);
dist_type hamming(sketch_type s, sketch_type t);
void filtering_by_sequential_search_using_kNN_buffer(struct_query_sketch *qs, int num_data, sketch_type sketch[], kNN_buffer *buff, int data_num[], int num_candidates);
void filtering_by_sequential_search_using_quick_select_k(struct_query_sketch *qs, int num_data, sketch_type sk[], dist_type sc[], int data_num[], int num_candidates);

#if defined(NARROW_SKETCH)
interval_list *new_interval_list(unsigned int nt, unsigned int size);
void realloc_interval_list(interval_list *ivl, unsigned int size);
vlist *new_vlist(int size, int step);
void add_vlist(vlist *vl, interval i);
void makenull_vlist(vlist *vl);
void enum_sub_dimension_0(sub_dimension sub, int begin, int remain, sub_dimension tab[], int *k, int K, int dim);
void rearrange_table(int nt, sub_dimension *tab, int tab_size);
sub_dimension *make_table_for_enumeration_hamming(int n, int dim);
sketch_type Mu(int b, int lg, int idx[], int i, sub_dimension *sd);
int nextComb(int n);
void reset_max_H(void);
#ifdef CHECK_ENUMERATION
void check_sketch_enumeration_c2_n(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, struct_check_result result[], int num_candidates);
void check_sketch_enumeration_hamming(struct_query_sketch *qs, struct_bucket *bucket, struct_check_result result[], int num_candidates, int low, int add, int add_2);
#else
int filtering_by_sketch_enumeration_hamming(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates);
int filtering_by_sketch_enumeration_hamming_interval(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates);
int filtering_by_sketch_enumeration_hamming_interval_2(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates);
int filtering_by_sketch_enumeration_hamming_interval_3(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates);
int filtering_by_sketch_enumeration_hamming_interval_4(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates);
int filtering_by_sketch_enumeration_hamming_select(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates);
int filtering_by_sketch_enumeration(struct_query_sketch *qs, struct_bucket *bucket, struct_que *que, int data_num[], int num_candidates);
int filtering_by_sketch_enumeration_interval(struct_query_sketch *qs, struct_bucket *bucket, struct_que *que, interval_list *ivl, int num_candidates);
int filtering_by_sketch_enumeration_sketch(struct_query_sketch *qs, struct_bucket *bucket, struct_que *que, sketch_type sketch[], int num_candidates);
int filtering_by_sketch_enumeration_c2_n(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, int data_num[], int num_candidates);
#if !defined(USE_MU_COMMON) && PARA_ENUM_INF > 0
int filtering_by_sketch_enumeration_c2_n_interval(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que[], interval_list *ivl, int num_candidates);
#else
int filtering_by_sketch_enumeration_c2_n_interval(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, interval_list *ivl, int num_candidates);
#endif
int filtering_by_sketch_enumeration_c2_n_interval_st(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, interval_list *ivl, int num_candidates);
dist_type filtering_by_sketch_enumeration_c2_n_score(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, int num_candidates);
int filtering_by_sketch_enumeration_c2_n_sketch(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, sketch_type sketch[], int num_candidates);
#ifdef NEW_INF
int filtering_by_sketch_enumeration_inf(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates);
int filtering_by_sketch_enumeration_inf_interval(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates);
#else
int filtering_by_sketch_enumeration_inf(struct_query_sketch *qs, struct_bucket *bucket, sketch_type sketch[], int num_candidates);
//int filtering_by_sketch_enumeration_inf_serial(struct_query_sketch *qs, struct_bucket *bucket, sketch_type sketch[], int num_candidates);
int filtering_by_sketch_enumeration_inf_data(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates);
int filtering_by_sketch_enumeration_inf_data_select(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates);
int filtering_by_sketch_enumeration_inf_only(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates);
int filtering_by_sketch_enumeration_inf_only_2(struct_query_sketch *qs, struct_bucket *bucket, int *data_num_thread[], int num_data_thread[], int num_candidates);
int filtering_by_sketch_enumeration_inf_only_3(struct_query_sketch *qs, struct_bucket *bucket, int *data_num_thread[], int num_data_thread[], int num_candidates);
int filtering_by_sketch_enumeration_inf_only_4(struct_query_sketch *qs, struct_bucket *bucket, int *data_num_thread[], int num_data_thread[], int num_candidates);
int filtering_by_sketch_enumeration_inf_data_select_once(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates);
int filtering_by_sketch_enumeration_inf_only_sketch(struct_query_sketch *qs, struct_bucket *bucket, sketch_type sketch[], int num_candidates);
#endif
#endif
#endif
//void filtering_by_sequential_search_using_compact_bucket(struct_query_sketch *qs, int num_data, compact_bucket *bucket_ds, dist_type score[], int idx[], int data_num[], int num_candidates);
//void filtering_by_sequential_search_n(int num_queries, struct_query_sketch qs[], int num_data, sketch_type sketch[], kNN_buffer *can[]);
//void filtering_by_sequential_search_using_bucket();
//void filtering_by_sketch_enumeration();

int comp_sketch(sketch_type a, sketch_type b);
int comp_uint(const void *a, const void *b);
int find_pivot_for_sketch(int idx[], sketch_type sk[], int i, int j);
int partition_by_pivot_for_sketch(int idx[], sketch_type sk[], int i, int j, sketch_type piv);
void quick_sort_for_sketch(int idx[], sketch_type sk[], int i, int j);

// Metric space
typedef enum {EUCLID, L_1} metric;
// when L_P is used, ativate the following, add necessary program code, and compile with option like "-DL_P=1.3".
//#ifndef L_P
//#define L_P 1.3
//#endif

#ifndef MERTIC
#define METRIC EUCLID
// DeCAF descriptor assumes Euclidean
#endif

// typedef unsigned int dist_type;

#if METRIC == EUCLID
#define DISTANCE    dist_L2
#define DISTANCE_2  dist_L2_2
#define DISTANCE_22 dist_L2_22
#define SET_DIST    set_dist_L2_22
#define DISTANCE_PIVOT    dist_pivot_L2
#define DISTANCE_PIVOT_2  dist_pivot_L2_2
#define DISTANCE_PIVOT_22 dist_pivot_L2_22
#ifdef PARTITION_TYPE_PQBP
#define PART_DISTANCE    part_dist_L2
#define PART_DISTANCE_2  part_dist_L2_2
#define PART_DISTANCE_22 part_dist_L2_22
#define SET_PART_DIST    set_part_dist_L2_22
#define PART_DISTANCE_PIVOT    part_dist_pivot_L2
#define PART_DISTANCE_PIVOT_2  part_dist_pivot_L2_2
#define PART_DISTANCE_PIVOT_22 part_dist_pivot_L2_22
#endif
#else
#define DISTANCE    dist_L1
#define DISTANCE_2  dist_L1_2
#define DISTANCE_22 dist_L1_22
#define SET_DIST    set_dist_L2_22// set function is common for both L1 and L2
#ifdef PARTITION_TYPE_PQBP
#define PART_DISTANCE    part_dist_L1
#define PART_DISTANCE_2  part_dist_L1_2
#define PART_DISTANCE_22 part_dist_L1_22
#define SET_PART_DIST    set_part_dist_L2_22// set function is common to L2
#endif
#endif

// data_num: データ番号, dist: 質問との距離（フィルタリングでは，score_p）
//#ifdef NUM_NN
//	int data_num[NUM_NN];			// データ番号（ftrファイル内での通し番号）
//	dist_type dist[NUM_NN];			// 質問との距離（フィルタリングでは，score_p）
//typedef struct {
//	int data_num[NUM_NN];			// データ番号（ftrファイル内での通し番号）
//	dist_type dist[NUM_NN];			// 質問との距離（フィルタリングでは，score_p）
//} answer_type_NN;
//#endif

//typedef struct {
//	int data_num;			// データ番号（ftrファイル内での通し番号）
//	dist_type dist;			// 質問との距離（フィルタリングでは，score_p）
//	#ifdef ANSWER_WITH_DATA_ID
//	auxi_type data_id;		// データID（ftrの補助情報）（処理が重たくなるので削除）
//	#endif
//} answer_type;

// query_num: 質問番号, ftr: 特徴データ（ポインタ）
//typedef struct {
//	int query_num;			// 質問番号（質問のftrファイル内での通し番号）
//	ftr_type ftr;			// 特徴データ
//	// auxi_type query_id;		// 質問ID（ftrの補助情報）（処理が重たくなるので削除）
//} query_type;

//typedef struct {
//	answer_type *buff;		// 解 (id, dist) を入れるバッファー（大きさは 2k）
//	int k, num;				// バッファーの大きさ k と，格納している解の個数
//	dist_type k_nearest;	// k-NN の距離（ただし，最後にソートして
//} kNN_buffer;

dist_type dist_L1(ftr_type a, ftr_type b, int dim);
dist_type dist_L1_2(ftr_type a, ftr_type b, int dim, dist_type dist);
// void set_dist_L1_22(ftr_type a); // set function is common to L2
dist_type dist_L1_22(ftr_type b); 

dist_type dist_L2(ftr_type a, ftr_type b, int dim);
dist_type dist_L2_2(ftr_type a, ftr_type b, int dim, dist_type dist);
void set_dist_L2_22(ftr_type a);
dist_type dist_L2_22(ftr_type b);

#ifdef PARTITION_TYPE_PQBP
dist_type part_dist_L2(ftr_type a, ftr_type b, int dim_start, int dim);
dist_type part_dist_L2_2(ftr_type a, ftr_type b, int dim_start, int dim, dist_type dist);
void set_part_dist_L2_22(ftr_type a, int dim_start, int dim);
dist_type part_dist_L2_22(ftr_type b, int dim_start, int dim);
#endif

#if defined(PARTITION_TYPE_SQBP) || defined(PARTITION_TYPE_CSQBP)
dist_type sub_dist_L2(ftr_type a, pivot_type *pivot, int dim);
dist_type sub_dist_L2_2(ftr_type a, pivot_type *pivot, int dim, dist_type dist);
void set_sub_dist_L2_22(pivot_type *pivot, int dim);
dist_type sub_dist_L2_22(ftr_type b);
#endif

#ifndef USE_PD_SKETCH
dist_type dist_pivot_L2(ftr_type a, ftr_type piv_center, int dim);
dist_type dist_pivot_L2_2(ftr_type a, ftr_type piv_center, int dim, dist_type dist);
dist_type dist_pivot_L2_22(ftr_type b);
#else
dist_type dist_pivot_L2(ftr_type a, ftr_type piv_center, int num_axis, int axis[]);
dist_type dist_pivot_L2_2(ftr_type a, ftr_type piv_center, int num_axis, int axis[], dist_type dist);
dist_type dist_pivot_L2_22(ftr_type b, int num_axis, int axis[]);
#endif
#ifdef PARTITION_TYPE_PQBP
dist_type part_dist_pivot_L2(ftr_type a, ftr_type piv_center, int dim_start, int dim);
dist_type part_dist_pivot_L2_2(ftr_type a, ftr_type piv_center, int dim_start, int dim, dist_type dist);
dist_type part_dist_pivot_L2_22(ftr_type b, int dim_start, int dim);
#endif

int find_pivot_answer(answer_type ans[], int i, int j);
#ifdef PUBMED23
int partition_by_pivot_answer(answer_type ans[], int i, int j, answer_type piv);
#else
int partition_by_pivot_answer(answer_type ans[], int i, int j, dist_type piv);
#endif
void quick_sort_answer(answer_type ans[], int i, int j);
void quick_select_k_r_answer(answer_type ans[], int i, int j, int k);
void quick_select_k_answer(answer_type ans[], int i, int j, int k);

int find_pivot_sketch_with_priority_num(sketch_with_priority_num buff[], int i, int j);
int partition_by_pivot_sketch_with_priority_num(sketch_with_priority_num buff[], int i, int j, dist_type piv, int *sum);
int quick_select_sum_k_sketch_with_priority_num(sketch_with_priority_num buff[], int i, int j, int sum);

#ifdef INTERVAL_WITH_PRIORITY
void quick_sort_interval(interval buff[], int i, int j);
int find_pivot_interval(interval buff[], int i, int j);
int partition_by_pivot_interval(interval buff[], int i, int j, dist_type piv, int *sum);
int partition_by_pivot_interval_no_sum(interval buff[], int i, int j, dist_type piv);
int quick_select_sum_k_interval(interval buff[], int i, int j, int sum);
int quick_select_sum_k_interval_by_multi_thread(interval_list *ivl, int sum);
// for multi-thread
// int find_pivot_answer_mt(answer_type ans[], int i, int j);
#endif
int partition_by_pivot_answer_mt(answer_type ans[], int i, int j, dist_type piv);
int quick_select_k_answer_mt(int nt, int k_th[], int size, answer_type ans[], int k);

int balance_interval_list(interval_list *ivl);

kNN_buffer *new_kNN_buffer(int k);
void free_kNN_buffer(kNN_buffer *b);
dist_type make_empty_kNN_buffer(kNN_buffer *b);
dist_type push_kNN_buffer(answer_type *a, kNN_buffer *b);
int comp_answer(const void *a, const void *b);
dist_type flush_kNN_buffer(kNN_buffer *b);
dist_type final_flush_kNN_buffer(kNN_buffer *b);
dist_type merge_kNN_buffer(kNN_buffer *b, kNN_buffer *pool[], int n);
answer_type *get_kNN_buffer(kNN_buffer *b, int i);
answer_type *read_correct_answer(char *answer_csv_file, int num_queries);
//#ifdef NUM_NN
answer_type_NN *read_correct_answer_NN(char *answer_csv_file, int num_queries);
//#endif
void sort_answer(int num, answer_type ans[]);
int answer_check_interval(struct_bucket *bucket, answer_type *ans, interval_list *ivl, kNN_buffer *top_k);
int answer_check(answer_type *ans, int num_c, int candidate[], kNN_buffer *top_k);
void search_kNN(dataset_handle *dh, query_type *qr, int num_c, int candidate[], kNN_buffer *top_k);
void search_NN(dataset_handle *dh, query_type *qr, int num_c, int candidate[], kNN_buffer *top_k);
void search_kNN_on_ram(struct_ftr_id ftr_id[], query_type *qr, int num_c, int candidate[], kNN_buffer *top_k);
void search_NN_on_ram(struct_ftr_id ftr_id[], query_type *qr, int num_c, int candidate[], kNN_buffer *top_k);

void out_result(char *filename, int num_queries, answer_type ans[], kNN_buffer *top_k[]);
void out_result_NN(char *filename, int num_queries, answer_type_NN ans[], kNN_buffer *top_k[]);
void out_result2(char *filename, int w, int d, double etime, double num_c, int num_queries, answer_type ans[], kNN_buffer *top_k[], int sorted);
// #ifdef DOUBLE_FILTERING
// double out_result_double(char *filename, int n_w, int e_w, int d, double, double, double, double, double, double nc_1st, double nc_2nd, int num_queries, answer_type ans[], kNN_buffer *top_k[]);
double out_result_double(char *filename, int n_w, int e_w, int d, double e_time_1st, double e_time_score, double e_time_2nd, double e_time_kNN,  double etime, double nc_1st, double nc_2nd, 
                         int num_queries, answer_type ans[], kNN_buffer *top_k[], char *query_file);
double recall_without_search(int num_queries, answer_type ans[], kNN_buffer *top_k[]);
double recall_kNN_1(int num_queries, answer_type_NN ans[], kNN_buffer *top_k[]);
double recall_kNN(int num_queries, answer_type_NN ans[], kNN_buffer *top_k[]);
//#ifdef NUM_NN
double recall_without_search_NN(int num_queries, answer_type_NN ans[], kNN_buffer *top_k[]);
//#endif
// #endif

