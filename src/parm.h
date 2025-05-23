#ifndef COMPILE_TIME // 実際のコンパイル時にはCOMPILE_TIMEをdefineして，以下の定義を無効にすること．
                     // VSCODEでプログラムを作成している時に，未定義定数でエラーが出るのを防ぐために，
                     // 実際のコンパイル時に事前に定義して用いる定数を定義しておく． 

//#define DEEP1B
#define DECAF
#define NUM_K 10  // k-NN で求める近傍数
#define QUERY "query_xx.ftr"    // 質問データのftrファイル
#define SAMPLE_DATASET "base10M_00.ftr"     // ピボット選択時のサンプルデータセット
#define SAMPLE_DS_QBP "base_00.ftr"         // ピボットの中心点のQBPの元のデータセット
#define PIVOT_FILE "pivot.csv"
// #define INPUT_PIVOT_FILE "input_pivot.csv"  // 作成済みのピボットを改良するときに指定
#define PSMAP_PIVOT_FILE "psmap_pivot.csv"
#define SMAP_PIVOT_FILE "smap_pivot.csv"
#define PQBP_PIVOT_FILE "pqbp_pivot.csv"
#define BUCKET_FILE "pivot_range.bkt"
#define QUERY_FILE "query_xx.ftr"
#define ANSWER_FILE "query_xx.csv"
#define QUERY_2ND_FILE "query_yy.ftr"
#define ANSWER_2ND_FILE "query_yy.csv"
#define QUERY_3RD_FILE "query_zz.ftr"
#define ANSWER_3RD_FILE "query_zz.csv"
#define SFTR_FILE "base_yy_xx.ftr"
#define SEED 1 // random seed

#define WITHOUT_BKT
//#define PARTITION_TYPE_QBP
//#define SKETCH_PARTITION QBP // QBP
//#define SKETCH_PARTITION PQBP // PQBP
//#define SKETCH_PARTITION SQBP // SQBP
#define PARTITION_TYPE_PQBP
#define USE_PD
#define USE_PD_2
#define USE_PD_INC
#define USE_PD_SKETCH
#define NUM_AXIS 4
//#define PARTITION_TYPE_SQBP
//#define PARTITION_TYPE_CSQBP
//#define CS_PART_DIM 4 // CSQBP風のSMAP（連続する座標軸を用いる射影）で用いる部分空間の射影次元
//#define CS_NUM_PART (PJT_DIM / CS_PART_DIM) // CSQBP の部分空間の個数
//#define CS_PART_SKIP (FTR_DIM / CS_NUM_PART) // CSQBP の部分空間の開始次元の間隔
#define FIXED_NUM_J // CSQBPで部分空間の次元数を固定するオプション
#define FIXED_J_START
//#define SELECT_DIM 64 // number of selected dimensions by SQBP or CSQBP
//#define SELECT_SHIFT 5 // CSQBP の隣り合うピボットが使用する軸番号の差の10倍
//#define PRE_ROTATION
//#define TERNARY_QUANTIZATION 67
//#define IGNORE_MED
//#define USE_MEDIAN_AS_PIVOT
#define MAKE_PIVOT_FROM_OUTLIER

#define _OPENMP
#define NUM_THREADS 8
#define PARA_ENUM_INF 3
#define PARALLEL_ENUM 2
#define NARROW_SKETCH
//#define WIDE_SKETCH
//#define EXPANDED_SKETCH
#define REVERSE_IDX

#define PJT_DIM 26
#define SAMPLE_SIZE 10000
//#define FTR_ON_MAIN_MEMORY
#define FTR_ON_SECONDARY_MEMORY

#define FTR_ARRANGEMENT_ASIS

#define RESULT_FILE "result/${pivot}_${range}.csv"
#define NUM_CANDIDATES  100 // データセットに対する割合（パーミリアド（万分率））: 100 -> 1%
#define NUM_CANDIDATES1 200 // データセットに対する割合（パーミリアド（万分率））: 200 -> 2%
#define NUM_CANDIDATES2 300
#define NUM_CANDIDATES3 400
#define NUM_CANDIDATES4 500
#define NUM_CANDIDATES5 600

//#define FILTERING_BY_SKETCH_ENUMERATION_HAMMING
#define FILTERING_BY_SKETCH_ENUMERATION_C2N
//#define ALTERNATIVE_ENUMERATION
//#define ALTERNATIVE_5
#define WITH_ENUM_TABLE
#define WITH_MASK_TABLE
#define USE_INTERVAL
#define USE_COMPACT_INTERVAL
#define INTERVAL_WITH_PRIORITY
#define INTERVAL_WITH_RUN
//#define SELECT_BY_PRIORITY_AFTER_ENUMERATION
#define SELECT_SUM
#define SELECTION_BY_MULTI_THREAD
#define DATA_NUM_IN_SKETCH_ORDER
#define LOOP_CONTROL_BY_NUM_SKETCHES
#define SECOND_FILTERING_KNN_BUFFER
//#define SECOND_FILTERING_SELECT
#define STATIC_KNN_BUFFER_FOR_SEARCH
#define STATIC_DF_WORK

#define ENUM_DIM 5
#define ENUM_DIM_END 8
#define SPP_BIT 5
#define SPP_BIT_END 15
#define SPP_BIT_2 5
#define SPP_BIT_2_END 15
//#define CHECK_ENUMERATION
#define NEW_INF
#define INF_AND_HAMM
#define THREAD_PLUS 2

#define WITH_STAT

#define USE_MU
#define USE_MU_COMMON
//#define SELECT_BY_SINGLE
//#define SELECT_BY_PARA_MERGE
//#define FILTERING_BY_SKETCH_ENUMERATION
//#define FILTERING_BY_SKETCH_ENUMERATION_C2N
//#define FILTERING_BY_SKETCH_ENUMERATION_INF
//#define SEQUENTIAL_FILTERING
//#define SEQUENTIAL_FILTERING_USING_BUCKET
//#define SEQUENTIAL_FILTERING_USING_HAMMING
//#define DOUBLE_FILTERING

#define SELECT_K_BY_QUICK_SELECT_K_PARA_WORK

#define NUM_Q 100
#define NUM_D 10000
#define REPEAT 100
#define ACCESS_MODE "SEQUENTIAL"

//#define TINY_IN_CHAR
//#define QUANTIZE_MIXED_MOD3
//#define QUANTIZE_BIT_0 3
//#define QUANTIZE_BIT_1 3
//#define QUANTIZE_BIT_2 2
#define QUANTIZE_BIT 3
//#define USE_PACKED_3BIT
#define USE_PACKED_6BIT
//#define USE_PACKED_4BIT
#define QUANTIZE_RANGE 2.0

#define EXPANDED_PIVOT_FILE     "pivot_PQBP_t10_w192_sd00_ss10000_np6_sr1000_seed1.csv"
#define EXPANDED_BUCKET_FILE    "pivot_PQBP_t10_w192_sd00_ss10000_np6_sr1000_seed1_00_99.bkt"
#define NARROW_PIVOT_FILE       "pivot_QBP_t1000_w28_sd00_ss10000.csv"
#define NARROW_BUCKET_FILE      "pivot_QBP_t1000_w28_sd00_ss10000_00_99.bkt"

#define NUM_CANDIDATES_1ST      3000
#define NUM_CANDIDATES_1ST_END  8000
#define NUM_CANDIDATES_1ST_STEP 1000
#define NUM_CANDIDATES_2ND      1
#define NUM_CANDIDATES_2ND_END  5
#define NUM_CANDIDATES_2ND_STEP 1

#define SCORE_P_1ST 1.0
#define SCORE_P_2ND 1.0

#define PRIORITY 0 // Hamming

#define SCORE_1
#define PRIORITY 1

#define FACTOR_INF 100
#define FACTOR_INF2 10

#define PLUS_SD 0.67
// #define QUANTIZE_BIT 2
#define USE_LOWER_BOUND

#define USE_LOG_LOG_SCORE
//#define USE_LOG_LOG_SCORE_QUICK_SELECT
#define USE_LOG_LOG_SCORE_QUICK_SELECT_ANSWER
//#define USE_LOG_LOG_SCORE_KNN_BUFFER
#define SCORE_TOP
//#define SCORE_2
//#define PRIORITY 2

//#define SCORE_INF
//#define PRIORITY 3

#define NUM_TRIAL_QBP 100 // Number of Trials for each bit in incremental QBP
#define SAMPLE_SIZE_QBP 5000 // Number of samples to calculate collisions
#define SAMPLE_SIZE_RADIUS 1000 // Number of samples to calculate radius of pivot

#define SAMPLE_NUM_UNBALANCED_RAD 10 // Define this when making unbalanced partitioninig.
#define USE_MED // Define this when median is used to compute radius.
#define NUM_TRIAL_LS  0	  // Number of Trials in Local Search
// Macro for AIR
#define USE_AIR
#ifdef USE_AIR
#define NUM_TRIAL_AIR 500 // Number of Trials in AIR
//#define DUAL_SAMPLE_DS // Use another sample dataset to compare the evaluations of two sample datasets
#define TRUNCATE_AT 10    // Number of candidates (% to num_data) to use evaluate pivots in optimization by AIR 
#define SAMPLE_SIZE_AIR 10000 // NUmber of sample queries for AIR
#define RETRY 200
#define REUSE 10	//取り換え頻度
#define T_POWER 2
#define N0 100				//初期サンプルサイズ（データ）
#define Q0 50				//初期サンプルサイズ（質問）
#define P0 5				//初期サンプルサイズ（ピボット群）

#define INTERMEDIATE
//#define REPLACE_WHOLE
//#define REPLACE_TRIAL 50

//#define HYBRID_REPLACE_FLIP

#define EVAL_BY_SEQUENTIAL_FILTERING // AIRでピボットの評価に sequential filtering を用いるときに #define
//#define EVAL_BY_RANK_SCORE
#define USE_QUICK_SELECT
#define USE_QUICK_SELECT_PARA_WORK

//#define PARTITION_TYPE_PQBP
//#define USE_QUICK_SELECT_IN_EVALUATION

#define TRY_FLIP_IN_ROTATION
//#define USING_CORRELATION
//#define EVALUATE_BY_DISTANCE_PRESERVATION
#define EVALUATE_BY_RANKING
#define USE_TOP_PIVOT_SETS_FROM 100

// #define EVAL_PRECISION precision_resample_cursor_2_n
#define EVAL_PRECISION precision_resample_by_sequential_filtering
#define NUM_NN 10

#define K1 0.6645		    //95%到達…500-0.2802 100-0.5488
//#define K1 0.2802
#define K2 1.0
// #define K 131028

// #define SKIP_ENQ		// score による打ち切りを行うときに #define 
//#ifdef SKIP_ENQ
//double expand = 1.0 ; // 射影距離を引き延ばしたり，縮めたりする係数
//#define SKIP_CONDITION(pdist, dist) ((pdist) * expand > (dist))
// スケッチの列挙で用いる優先度付き待ち行列にスケッチをENQ（挿入）するときに，
// 射影距離 pdist × expand がその時点での最近傍暫定解との距離 dist を超えていたら，ENQ をしない．
// DISA の DeCAF記述子 においては，実距離に比べて非常に大きいので，expand は 0.1 のように小さくすること．
//#endif

#define SCORE_INF
#define SCORE_P 10
#define SCORE_FACTOR 0.4

//#define REPLACE_WHOLE
#define REPLACE_TRIAL 10

#define VAR_FLIP
#define NUM_FLIPS 1

#ifndef EVAL_MODE		// 評価にlogを使った重み付きスコア総和を用いるときは 1 にする．
#define EVAL_MODE 0
#endif

#define NUM_SAMPLE_DATA 100000
#define NUM_SAMPLE_QUERIES 1000
#define INCLUDE_KNN_FOR_EACH_QUERY

#endif // USE_AIR
#endif // COMPILE_TIME
// これより下には追加しないこと