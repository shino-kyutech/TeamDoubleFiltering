#pragma once
//#include <stdio.h>
//#include <stdlib.h>
//#include <string.h>
//#include <math.h>
//#ifdef _OPENMP
//#include <omp.h>
//#endif
//#include "parm.h"
//#include "config.h"
#include "ftr.h"
#include "bit_op.h"
//#include "e_time.h"
#include "sketch.h"
//#include "quick.h"
//#include "pivot_selection.h"

#ifndef SMAP_DIM
#define SMAP_DIM 96
#endif

#ifndef SMAP_NUM_PART
#define SMAP_NUM_PART 1
#endif
#define PART_SMAP_MIN (SMAP_DIM / SMAP_NUM_PART)
#define PART_SMAP_MAX ((SMAP_DIM + SMAP_NUM_PART - 1) / SMAP_NUM_PART)

// 以下は，USE_PDが定義されている時は，初期ピボット割り当て時のみ使用する．
// SMAP_PART_NUM(p) = 射影次元（p = 0, ... , SMAP_DIM - 1）に対応する部分空間番号
#define SMAP_PART_NUM(p) (((p) / PART_SMAP_MAX) < (SMAP_DIM % SMAP_NUM_PART) ? ((p) / PART_SMAP_MAX) : ((p) - SMAP_DIM % SMAP_NUM_PART) / PART_SMAP_MIN)
// SMAP_PART_START(p), SMAP_PART_DIM(p), SMAP_PART_END(p), PART_SMAP_DIM(p) 射影次元（p = 0, ... , SMAP_DIM - 1）に対応する部分空間の
// 開始次元番号, 次元数, 最終次元番号, 射影次元数
#define SMAP_PART_START(j) ((FTR_DIM / SMAP_NUM_PART) * SMAP_PART_NUM(j) + (SMAP_PART_NUM(j) < FTR_DIM % SMAP_NUM_PART ? SMAP_PART_NUM(j) : FTR_DIM % SMAP_NUM_PART))
#define SMAP_PART_DIM(j) ((FTR_DIM / SMAP_NUM_PART) + (SMAP_PART_NUM(j) < FTR_DIM % SMAP_NUM_PART ? 1 : 0))
#define SMAP_PART_END(j) (SMAP_PART_START(j) + SMAP_PART_DIM(j) - 1)
#define PART_SMAP_DIM(j) ((SMAP_DIM / SMAP_NUM_PART) + (SMAP_PART_NUM(j) < SMAP_DIM % SMAP_NUM_PART ? 1 : 0))

typedef struct {
	ftr_element_type p[SMAP_DIM][FTR_DIM];	// SMAP の中心点
	#ifdef USE_PD
	int axis[SMAP_DIM][FTR_DIM];
	int num_axis[SMAP_DIM];
	#endif
	int axis_pool[FTR_DIM];
	int num_used_axis;
	partition_type type;				// 基礎分割関数（GHP = 0, BP = 1, QBP = 2, PQBP = 3, SQBP = 4）（定数は enum で定義されている）
} smap_pivot_type;

typedef unsigned short smap_element_type, *smap_type; // SMAPの射影像の各次元の型（16-bit）

#ifndef P_DISTANCE
#define P_DISTANCE projected_dist_L2
#endif

#ifndef Q_P_DISTANCE
#define Q_P_DISTANCE q_projected_dist_L2
#endif

smap_pivot_type *new_smap_pivot(int type);
void free_smap_pivot(smap_pivot_type *pivot);
void read_smap_pivot(char *filename, smap_pivot_type *pivot);
void write_smap_pivot(char *filename, smap_pivot_type *pivot);

#ifdef TERNARY_QUANTIZATION
void select_random_pivot_for_psmap(smap_pivot_type *pivot, double ave[], double stdev[], struct_dataset *ds);
#else
void select_random_pivot_for_psmap(smap_pivot_type *pivot, ftr_type median, struct_dataset *ds);
#endif

#ifdef USE_PD
dist_type partition_dist_L2(ftr_type a, ftr_type b, int num_axis, int axis[]);
#else
dist_type partition_dist_L2(ftr_type a, ftr_type b, int dim_start, int dim);
#endif
dist_type real_dist_L2(ftr_type a, ftr_type b);
#ifdef USE_PD
dist_type partition_dist_pivot_L2(ftr_type a, ftr_type piv_center, int num_axis, int axis[]);
#else
dist_type partition_dist_pivot_L2(ftr_type a, ftr_type piv_center, int dim_start, int dim);
#endif
dist_type real_dist_pivot_L2(ftr_type a, ftr_type piv_center);
dist_type projected_dist_L1(smap_type a, smap_type b);
dist_type projected_dist_L2(smap_type a, smap_type b);
dist_type projected_dist_L2_max(smap_type a, smap_type b);
dist_type projected_dist_Lp(smap_type a, smap_type b, double p);

#ifdef QUANTIZE_BIT

#ifndef TABLE_UNIT
#ifdef USE_PACKED_3BIT
#define TABLE_UNIT 3
#elif defined(USE_PACKED_6BIT)
#define TABLE_UNIT 2
#else
#define TABLE_UNIT (8 / QUANTIZE_BIT)
#endif
#endif

#ifdef QUANTIZE_MIXED_MOD3
	#define PACKED_QPSMAP_SIZE ((SMAP_DIM + 2) / 3)
#else
	#define PACKED_QPSMAP_SIZE ((SMAP_DIM * QUANTIZE_BIT + (sizeof(tiny_int) * 8 - 1)) / (sizeof(tiny_int) * 8))
#endif

void compute_parmeters_for_qpsmap(smap_pivot_type *pivot, struct_dataset *ds_query, double ave[], double stdev[], double offset[], double slice[]);
// void make_table_for_query(smap_type smap, unsigned int table[][1 << QUANTIZE_BIT], double offset[], double slice[], int range);
// void make_table_for_query_p(smap_type smap, unsigned int table[][1 << QUANTIZE_BIT], double offset[], double slice[], int range, double p);
void make_table_for_query(smap_type smap, unsigned int table[][1 << QUANTIZE_BIT], double offset[], double slice[]);
void make_table_for_query_p(smap_type smap, unsigned int table[][1 << QUANTIZE_BIT], double offset[], double slice[], double p);
#if defined(USE_PACKED_3BIT) || (defined(USE_PACKED_6BIT) && !defined(QUANTIZE_MIXED_MOD3))
void make_table_for_packed_data(unsigned int table[][1 << QUANTIZE_BIT], unsigned int p_table[][1 << QUANTIZE_BIT * TABLE_UNIT], double offset[], double slice[]);
dist_type projected_dist_packed_table(unsigned int p_table[][1 << QUANTIZE_BIT * TABLE_UNIT], tiny_int *qpsmap);
#elif defined(QUANTIZE_MIXED_MOD3) && defined(USE_PACKED_6BIT)
void make_table_for_packed_data(unsigned int table[][1 << QUANTIZE_BIT], unsigned int p_table[][64], double offset[], double slice[]);
dist_type projected_dist_packed_table(unsigned int p_table[][64], tiny_int *qpsmap);
#else
// void make_table_for_packed_data(unsigned int table[][1 << QUANTIZE_BIT], unsigned int p_table[][256], double offset[], double slice[]);
// dist_type projected_dist_packed_table(unsigned int p_table[][256], tiny_int *qpsmap);
void make_table_for_packed_data(unsigned int table[][1 << QUANTIZE_BIT], unsigned int p_table[][1 << QUANTIZE_BIT * TABLE_UNIT], double offset[], double slice[]);
dist_type projected_dist_packed_table(unsigned int p_table[][1 << QUANTIZE_BIT * TABLE_UNIT], tiny_int *qpsmap);
#endif
void make_table_for_packed_data_inf(unsigned int table[][1 << QUANTIZE_BIT], unsigned int p_table[][256], double offset[], double slice[]);
dist_type projected_dist_table(unsigned int table[][1 << QUANTIZE_BIT], unsigned char qpsmap[]);
dist_type projected_dist_table_inf(unsigned int table[][1 << QUANTIZE_BIT], unsigned char qpsmap[]);
dist_type projected_dist_packed_table_inf(unsigned int p_table[][256], tiny_int *qpsmap);
#endif

smap_element_type psmap_dim(int dim, ftr_type ftr, smap_pivot_type *pivot);
void psmap(smap_type projected, ftr_type ftr, smap_pivot_type *pivot);
smap_element_type q_psmap_dim(int dim, ftr_type ftr, smap_pivot_type *pivot, double offset, double slice);
void q_psmap(smap_type q_projected, ftr_type ftr, smap_pivot_type *pivot, double offset[], double slice[]);
void q_uchar_psmap(unsigned char q_projected[], ftr_type ftr, smap_pivot_type *pivot, double offset[], double slice[]);
void psmap2uchar_qpsmap(smap_type sm, unsigned char *uchar_smap, double offset[], double slice[], int range);

#ifndef PLUS_HALF
#define PLUS_HALF 0.5
#endif
