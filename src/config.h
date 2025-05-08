#pragma once

#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <math.h>
#include <string.h>
#include <omp.h> 
#include <float.h>
#include <time.h>
#include "parm.h"

typedef enum {MAIN_MEMORY, SECONDARY_MEMORY} ftr_on_type;

#ifdef DOUBLE_FILTERING
//#define FTR_ON_SECONDARY_MEMORY
//#define FTR_ARRANGEMENT_ASIS
#define PARTITION_TYPE_PQBP
#define PARTITION_TYPE_QBP
#endif

#if !(defined(FTR_ON_MAIN_MEMORY) ^ defined(FTR_ON_SECONDARY_MEMORY))
#error "Just one of { FTR_ON_MAIN_MEMORY | FTR_ON_SECONDARY_MEMORY } should be defined."
#endif
 
#if !(defined(FTR_ARRANGEMENT_ASIS) ^ defined(FTR_SORT_BY_SKETCH) ^ defined(FTR_SORT_BY_SKETCH_IN_ADVANCE))
#error "Just one of { FTR_ARRANGEMENT_ASIS | FTR_SORT_BY_SKETCH | FTR_SORT_BY_SKETCH_IN_ADVANCE } should be defined."
#endif

#if !(defined(SEQUENTIAL_FILTERING) ^ defined(SEQUENTIAL_FILTERING_USING_BUCKET) ^ defined(SEQUENTIAL_FILTERING_USING_HAMMING) ^ defined(FILTERING_BY_SKETCH_ENUMERATION_HAMMING) ^ defined(FILTERING_BY_SKETCH_ENUMERATION) ^ defined(FILTERING_BY_SKETCH_ENUMERATION_C2N) ^ defined(FILTERING_BY_SKETCH_ENUMERATION_INF) ^ defined(WITHOUT_FILTERING) ^ defined(DOUBLE_FILTERING))
#error "Just one of { SEQUENTIAL_FILTERING | SEQUENTIAL_FILTERING_USING_BUCKET | SEQUENTIAL_FILTERING_USING_HAMMINGT | FILTERING_BY_SKETCH_ENUMERATION_HAMMING | FILTERING_BY_SKETCH_ENUMERATION | FILTERING_BY_SKETCH_ENUMERATION_C2N | WITHOUT_FILTERING | DOUBLE_FILTERING } should be defined."
#endif

#ifndef WITHOUT_FILTERING
#if !(defined(PARTITION_TYPE_GHP) ^ defined(PARTITION_TYPE_BP) ^ defined(PARTITION_TYPE_QBP) ^ defined(PARTITION_TYPE_PQBP) ^ defined(PARTITION_TYPE_SQBP) ^ defined(PARTITION_TYPE_CSQBP) ^ defined(PARTITION_TYPE_QBP_AND_PQBP))
#error "Just one of { PARTITION_TYPE_[GHP | BP | QBP | PQBP | SQBP | CSQBP| QBP_AND_PQBP] } should be defined when \"filtering by skech\"."
#endif
#endif

#ifndef BLOCK_SIZE
#define BLOCK_SIZE 200
#endif

typedef enum {GHP, BP, QBP, PQBP, SQBP, CSQBP} partition_type;

#define MAX_LEN 100000	// ファイル名の編集に用いる文字列ややcsvファイルの読込のバッファの大きさ
//#ifndef EXPANDED_SKETCH
//#define BIT (1 << PJT_DIM) 
//#endif

#define CLOCK_METHOD CLOCK_MONOTONIC // 経過時間と求めるときのパラメータ
//#define CLOCK_METHOD CLOCK_REALTIME // 経過時間と求めるときのパラメータ
