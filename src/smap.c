//#pragma once
#include "parm.h"
#include "bit_op.h"
#include "smap.h"
#include "e_time.h"
#ifdef USE_PD_INC
#include "pivot_selection.h"
#endif

// ピボット（中心点の配列とタイプ = smap_pivot_type）の構造体を動的に確保する．
// type は，基礎分割関数（GHP = 0, BP = 1, QBP = 2, PQBP = 3）を指定する（現状ではQBPとPQBPのみ）．
// スケッチ幅は，SMAP_DIMで指定したものを用いる．
smap_pivot_type *new_smap_pivot(int type)
{
	smap_pivot_type *pivot = MALLOC(sizeof(smap_pivot_type));
    #ifdef USE_PD
	for(int j = 0; j < SMAP_DIM; j++) {
        pivot->num_axis[j] = 0;
        for(int i = 0; i < FTR_DIM; i++) {
            pivot->axis[j][i] = i;
        }
	}
    #endif
    for(int i = 0; i < FTR_DIM; i++) {
        pivot->axis_pool[i] = i;
    }
    pivot->num_used_axis = 0;
	pivot->type = type;
	return pivot;
}

void free_smap_pivot(smap_pivot_type *pivot)
{
    FREE(pivot, sizeof(smap_pivot_type));
}

// ピボットをファイルから読み込む．
// 作成済みのピボットを用いて検索などを行うときに使用する．
// filename = ピボットを格納したファイル名（csv形式）．ただし，ファイル内の形式はスケッチのピボットと共通なので，半径は無視．
// pivot = ピボットの構造体（割り当ては，new_smap_pivotを用いて別途行っておくこと）
void read_smap_pivot(char *filename, smap_pivot_type *pivot)
{
	FILE *pfp ;
    char buf[100000] = {0};
	int i, j;

	pfp=fopen(filename, "r");
	if(pfp == NULL){
		fprintf(stderr, "cannot open pivot file = %s\n", filename);
		exit(0);
	}

	// 基礎分割関数
	if(fgets(buf, MAX_LEN, pfp) != NULL)
		pivot->type = atoi(strtok(buf, ","));

	printf("pivot-type = %d\n", pivot->type);

	// 半径
	if(fgets(buf, MAX_LEN, pfp) != NULL) {
        // smap のピボットは半径がないので，読み捨てる．
	}
	
	// 中心点
    for(i = 0; fgets(buf, MAX_LEN, pfp) != NULL; i++) {
        pivot->p[i][0] = atoi(strtok(buf, ","));
        for(j = 1; j < FTR_DIM; j++) {
            pivot->p[i][j] = atoi(strtok(NULL, ","));
        }
        #ifdef USE_PD
        int na = 0;
        for(j = 0; j < FTR_DIM; j++) {
            if(pivot->p[i][j] == FTR_MIN || pivot->p[i][j] == FTR_MAX) {
                ftr_element_type t = pivot->axis[i][na];
                pivot->axis[i][na++] = j;
                pivot->axis[i][j] = t;
            }
        }
        pivot->num_axis[i] = na;
        fprintf(stderr, "num_axis[%d] = %d\n", i, pivot->num_axis[i]);
        #endif
    }
/*
    #ifdef USE_PD
    for(int dim = 0; dim < SMAP_DIM; dim++) {
        printf("dim = %3d: num_axis = %d, ", dim, pivot->num_axis[dim]);
        for(int i = 0; i < pivot->num_axis[dim]; i++) {
            printf("axis[%d] = %d, ", i, pivot->axis[dim][i]);
        }
        printf("\n");
        printf("unused = ");
        for(int i = pivot->num_axis[dim]; i < FTR_DIM; i++) {
            printf("%4d", pivot->axis[dim][i]);
        }
        printf("\n");
        getchar();
    }
    getchar();
    #endif
*/

	fclose(pfp); 
}

void write_smap_pivot(char *filename, smap_pivot_type *pivot)
{
	int i, j;
	FILE *fp = fopen(filename, "w");
	if(fp == NULL) {
		printf("cannot write open file = %s\n", filename);
	} else {
		if(pivot->type == PQBP) {
			fprintf(fp, "%d, %d\n", pivot->type, SMAP_NUM_PART);		// 基礎分割関数, 分割数
		} else {
			fprintf(fp, "%d, 1\n", pivot->type);				// 基礎分割関数，分割数 = 1
		}
		for(i = 0; i < SMAP_DIM - 1; i++)	 	// 半径
			fprintf(fp, "%d,", 0); // smap_pivot には半径がないので 0 を書き込む
		fprintf(fp, "%d\n", 0);
//		fprintf(fp, "\n"); // smap_pivot には半径がないので空行を書き込む
		for(i = 0; i < SMAP_DIM; i++) {			// 中心点
			#ifdef PARTITION_TYPE_SQBP
				for(j = pivot->num_used[i]; j < FTR_DIM; j++) {
					pivot->p[i][pivot->used[i][j]] = FTR_IGN;		// 未使用部分をFTR_IGNに変更
				}
			#elif defined(PARTITION_TYPE_CSQBP)
				int used[FTR_DIM] = {0};
				for(int c = 0, j = pivot->j_start[i]; c < pivot->num_j[i]; c++) {
					used[j++] = 1;
					if(j == FTR_DIM) j = 0;			
				}
				for(j = 0; j < FTR_DIM; j++) {
					if(!used[j]) pivot->p[i][j] = FTR_IGN;
				}
			#endif
			for(j = 0; j < FTR_DIM - 1; j++)
				fprintf(fp, "%d,", pivot->p[i][j]);
			fprintf(fp, "%d\n", pivot->p[i][j]);
		}
	}
	fclose(fp);
}

#ifdef TERNARY_QUANTIZATION
void select_random_pivot_for_psmap(smap_pivot_type *pivot, double ave[], double stdev[], struct_dataset *ds)
{
	double tq = 0.01 * TERNARY_QUANTIZATION;
	for(int dim = 0; dim < SMAP_DIM; dim++) {
		int c = random() % ds->num_data;
		for(int i = 0; i < FTR_DIM; i++) {
            // psmapではSMAP_PART_START(dim)からSMAP_PART_DIM(dim)次元分しか用いないので，その部分は量子化して，未使用部分は平均値にする．
            if(i >= SMAP_PART_START(dim) && i < SMAP_PART_START(dim) + SMAP_PART_DIM(dim)) {
                if(ave[i] - tq * stdev[i] > 0 && ds->ftr_id[c].ftr[i] < ave[i] - tq * stdev[i]) {
                    pivot->p[dim][i] = FTR_MIN;
                } else if(ds->ftr_id[c].ftr[i] <= (ftr_element_type)(ave[i] + tq * stdev[i] + 0.5)) {
                    pivot->p[dim][i] = (ftr_element_type)(ave[i] + 0.5);
                } else {
                    pivot->p[dim][i] = FTR_MAX;
                }
            } else {
                pivot->p[dim][i] = (ftr_element_type)(ave[i] + 0.5);
            }
		}
    }
}
#else
// データセットdsから乱択した特徴データを2値量子化（中央値medianを閾値とする）した点を中心点にしたピボットを求める
#ifndef USE_PD_INC
#if defined(USE_PD_2) && NUM_AXIS > 1
// SMAP_NUM_PARTに基づいて各射影軸に割当てられた特徴次元の半分を使用する．
void select_random_pivot_for_psmap(smap_pivot_type *pivot, ftr_type median, struct_dataset *ds)
{
	for(int dim = 0; dim < SMAP_DIM; dim++) {
        #ifdef USE_PD
        pivot->num_axis[dim] = 0;
        for(int i = 0; i < FTR_DIM; i++) {
            pivot->axis[dim][i] = i;
        }
        #endif
		int c = random() % ds->num_data;
        // psmapではSMAP_PART_START(dim)からSMAP_PART_DIM(dim)次元分しか用いないので，その部分は量子化して，未使用部分は中央値にする．
        for(int i = 0; i < FTR_DIM; i++) {
            if(i >= SMAP_PART_START(dim) && i < SMAP_PART_START(dim) + SMAP_PART_DIM(dim)) {
                // 割当てられた特徴次元の本数 = SMAP_PART_DIM
                // そのうち，平均して NUM_AXIS 本を使用するように乱択
                if(random() % SMAP_PART_DIM(dim) < NUM_AXIS) {
                    pivot->p[dim][i] = ds->ftr_id[c].ftr[i] < median[i] ? FTR_MIN : FTR_MAX;
//      [0], [num_axis[dim]], ... , i 
                    #ifdef USE_PD
                    int u = pivot->axis[dim][pivot->num_axis[dim]];
                    pivot->axis[dim][pivot->num_axis[dim]] = i;
                    pivot->axis[dim][i] = u;
                    pivot->num_axis[dim]++;
                    #endif
                } else {
                    pivot->p[dim][i] = median[i];
                }
            } else {
                pivot->p[dim][i] = median[i];
            }
		}
        printf("num_axis[%d] = %d\n", dim, pivot->num_axis[dim]);
/*      #ifdef USE_PD
        for(int i = 0; i < FTR_DIM; i++) {
            pivot->axis[dim][i] = i;
        }
        for(int i = 0; i < pivot->num_axis[dim]; i++) {
            int u = SMAP_PART_START(dim) + i;
            int t = pivot->axis[dim][i];
            pivot->axis[dim][i] = u;
            pivot->axis[dim][u] = t;
        }
        #endif
*/
    }
}
#else
void select_random_pivot_for_psmap(smap_pivot_type *pivot, ftr_type median, struct_dataset *ds)
{
	for(int dim = 0; dim < SMAP_DIM; dim++) {
		int c = random() % ds->num_data;
        // psmapではSMAP_PART_START(dim)からSMAP_PART_DIM(dim)次元分しか用いないので，その部分は量子化して，未使用部分は中央値にする．
        for(int i = 0; i < FTR_DIM; i++) {
            if(i >= SMAP_PART_START(dim) && i < SMAP_PART_START(dim) + SMAP_PART_DIM(dim)) {
                pivot->p[dim][i] = ds->ftr_id[c].ftr[i] < median[i] ? FTR_MIN : FTR_MAX;
            } else {
                pivot->p[dim][i] = median[i];
            }
		}
        #ifdef USE_PD
        for(int i = 0; i < FTR_DIM; i++) {
            pivot->axis[dim][i] = i;
        }
        for(int i = 0; i < SMAP_PART_DIM(dim); i++) {
            int u = SMAP_PART_START(dim) + i;
            int t = pivot->axis[dim][i];
            pivot->axis[dim][i] = u;
            pivot->axis[dim][u] = t;
        }
        pivot->num_axis[dim] = SMAP_PART_DIM(dim);
        #endif
    }
/*
    #ifdef USE_PD
    for(int dim = 0; dim < SMAP_DIM; dim++) {
        printf("dim = %3d: num_axis = %d, ", dim, pivot->num_axis[dim]);
        for(int i = 0; i < pivot->num_axis[dim]; i++) {
            printf("axis[%d] = %d, ", i, pivot->axis[dim][i]);
        }
        printf("\n");
        printf("unused = ");
        for(int i = pivot->num_axis[dim]; i < FTR_DIM; i++) {
            printf("%4d", pivot->axis[dim][i]);
        }
        printf("\n");
        getchar();
    }
    getchar();
    #endif
*/
}
#endif
#else // USE_PD_INC
// FTR_DIM の中から乱択した SMAP_DIM 個の特徴軸1本を射影軸とするピボット群を求める．
// FTR_DIM > SMAP_DIM のときは，未使用の特徴軸がある．
// FTR_DIM < SMAP_DIM のときは，当面未対応
void select_random_pivot_for_psmap(smap_pivot_type *pivot, ftr_type median, struct_dataset *ds)
{
    shuffle(pivot->axis_pool, FTR_DIM, SMAP_DIM); // 0, 1, ... , FTR_DIM - 1 から SMAP_DIM 個を乱択する．
    // dim 次元のピボットの中心点
    // num_axis[dim] = 1; // 特徴軸を 1 本だけ使用する．その特徴軸番号 = axis_pool[dim]
    // 元々，new_smap_pivot で作成したときには，axis[dim][i] = i (i = 0, ... , FTR_DIM - 1) に初期化されている．
    // axis[dim][0] と axis[dim][axis_pool[dim]] を交換する．
    // p[dim][ax[j]]  (j = 0, )= FTR_MIN or FTR_MAX とし，ax[dim] 以外のところの p[dim][i] は median[i] にしておく（未使用を表す）．
	for(int dim = 0; dim < SMAP_DIM; dim++) {
        for(int i = 0; i < FTR_DIM; i++) {
            pivot->p[dim][i] = median[i]; // まず，すべて未使用にしておく．
		}
        int u = pivot->axis[dim][0];
        int t = pivot->axis[dim][pivot->axis_pool[dim]];
        pivot->axis[dim][pivot->axis_pool[dim]] = u;
        pivot->axis[dim][0] = t;
        pivot->p[dim][pivot->axis[dim][0]] = random() % 2 ? FTR_MIN : FTR_MAX;
        pivot->num_axis[dim] = 1;
    }
    pivot->num_used_axis = SMAP_DIM;
/*
    for(int dim = 0; dim < SMAP_DIM; dim++) {
        printf("num_axis[%d] = %d\n", dim, pivot->num_axis[dim]);
        for(int j = 0; j < FTR_DIM; j++) {
            if(pivot->axis[dim][j] != j) {
                printf("axis[%d][%d] = %d\n", dim, j, pivot->axis[dim][j]);
            }
        }
    }
    getchar();
    printf("num_used = %d\n", pivot->num_used_axis);
    for(int j = 0; j < FTR_DIM; j++) {
        printf("axis_pool[%d] = %d\n", j, pivot->axis_pool[j]);
    }
    getchar();
*/
}
#endif
#endif

// 座標分割距離関数（L2, dim_startからdim次元分のみ用いる）（SMAP_NUM_PART = 1 のときは，全次元を用いる通常の距離関数になる）
// 注意：平方根を取る．
#ifdef USE_PD
dist_type partition_dist_L2(ftr_type a, ftr_type b, int num_axis, int axis[])
{
	dist_type s = 0;
	for(int i = 0; i < num_axis; i++)  {
        int j = axis[i];
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
	}
	return sqrt(s);
}
#else
dist_type partition_dist_L2(ftr_type a, ftr_type b, int dim_start, int dim)
{
	dist_type s = 0;
	for(int j = dim_start; j < dim_start + dim; j++)  {
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
	}
	return sqrt(s);
}
#endif

//特徴データ間の実距離関数（全体を用いるので rotation は無関係）
// 注意：平方根を取る．（いまのところ rotation は考えない）（そのうち ftr.h に統合）
dist_type real_dist_L2(ftr_type a, ftr_type b)
{
	int j;
	dist_type s = 0;
	
	for(j = 0; j < FTR_DIM; j++)  {
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
	}
	return sqrt(s);
}

#ifdef USE_PD
dist_type partition_dist_pivot_L2(ftr_type a, ftr_type piv_center, int num_axis, int axis[])
{
    static int first = 1;
    if(first) {
		#ifdef IGNORE_MED
        fprintf(stderr, "partition_dist_pivot_L2 (IGNORE_MED): num_axis = %d, axis[0] = %d\n", num_axis, axis[0]);
        #else
        fprintf(stderr, "partition_dist_pivot_L2 (!IGNORE_MED): num_axis = %d, axis[0] = %d\n", num_axis, axis[0]);
        #endif
        first = 0;
    }
	dist_type s = 0;
	for(int i = 0; i < num_axis; i++) {
        int j = axis[i];
		#ifdef IGNORE_MED
		if(piv_center[j] == FTR_MIN) {
			s += ((int)a[j] - FTR_MIN) * ((int)a[j] - FTR_MIN);
		} else if(piv_center[j] == FTR_MAX) {
			s += ((int)a[j] - FTR_MAX) * ((int)a[j] - FTR_MAX);
		}
		#else
		s += ((int)a[j] - (int)piv_center[j]) * ((int)a[j] - (int)piv_center[j]);
		#endif
	}
	return sqrt(s);
}
#else
dist_type partition_dist_pivot_L2(ftr_type a, ftr_type piv_center, int dim_start, int dim)
{
    static int first = 1;
    if(first) {
		#ifdef IGNORE_MED
        fprintf(stderr, "partition_dist_pivot_L2 (IGNORE_MED): dim_start = %d, dim = %d\n", dim_start, dim);
        #else
        fprintf(stderr, "partition_dist_pivot_L2 (!IGNORE_MED): dim_start = %d, dim = %d\n", dim_start, dim);
        #endif
        first = 0;
    }
	dist_type s = 0;
	for(int j = dim_start; j < dim_start + dim; j++) {
		#ifdef IGNORE_MED
		if(piv_center[j] == FTR_MIN) {
			s += ((int)a[j] - FTR_MIN) * ((int)a[j] - FTR_MIN);
		} else if(piv_center[j] == FTR_MAX) {
			s += ((int)a[j] - FTR_MAX) * ((int)a[j] - FTR_MAX);
		}
		#else
		s += ((int)a[j] - (int)piv_center[j]) * ((int)a[j] - (int)piv_center[j]);
		#endif
	}
	return sqrt(s);
}
#endif

dist_type real_dist_pivot_L2(ftr_type a, ftr_type piv_center)
{
    static int first = 1;
    if(first) {
		#ifdef IGNORE_MED
        fprintf(stderr, "real_dist_pivot_L2 (IGNORE_MED)\n");
        #else
        fprintf(stderr, "real_dist_pivot_L2 (!IGNORE_MED)\n");
        #endif
        first = 0;
    }
	int j;
	dist_type s = 0;
	
	for(j = 0; j < FTR_DIM; j++) {
		#ifdef IGNORE_MED
		if(piv_center[j] == FTR_MIN) {
			s += ((int)a[j] - FTR_MIN) * ((int)a[j] - FTR_MIN);
		} else if(piv_center[j] == FTR_MAX) {
			s += ((int)a[j] - FTR_MAX) * ((int)a[j] - FTR_MAX);
		}
		#else
		s += ((int)a[j] - (int)piv_center[j]) * ((int)a[j] - (int)piv_center[j]);
		#endif
	}
	return sqrt(s);
}

// psmapで射影された射影像間の射影距離(L1)
dist_type projected_dist_L1(smap_type a, smap_type b)
{
	dist_type sum = 0;
	for(int j = 0; j < SMAP_DIM; j++)  {
		sum += abs((int)a[j] - (int)b[j]);
	}
	return sum;
}

// psmapで射影された射影像間の射影距離(L2)
dist_type projected_dist_L2(smap_type a, smap_type b)
{
	double sum = 0;
	for(int j = 0; j < SMAP_DIM; j++)  {
		sum += (double)((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]) ;
	}
	return sqrt(sum);
}

// psmapで射影された射影像間の射影距離(L2-MAX)
dist_type projected_dist_L2_max(smap_type a, smap_type b)
{
	double sum = 0, d_max = 0;
    int n = -1;
	for(int p = 0; p < SMAP_DIM; p++)  {
        if(SMAP_PART_NUM(p) != n) { // 部分空間が切り替わった
            sum += d_max * d_max;
            d_max = 0;
            n = SMAP_PART_NUM(p);
        }
        double pd = abs((int)a[p] - (int)b[p]);
        if(pd > d_max) d_max = pd;
	}
	sum += d_max * d_max;
	return sqrt(sum);
}

// psmapで射影された射影像間の射影距離(Lp)
dist_type projected_dist_Lp(smap_type a, smap_type b, double p)
{
	double sum = 0;
	for(int j = 0; j < SMAP_DIM; j++)  {
		sum += pow(abs((int)a[j] - (int)b[j]), p) ;
	}
	return pow(sum, 1 / p);
}

#ifdef QUANTIZE_BIT
/* 質問の表関数を固定してみたが，高速化にはならなかった．
unsigned int setted_table[NUM_THREADS][SMAP_DIM][1 << QUANTIZE_BIT];

void set_projected_dist_table(unsigned int table[][1 << QUANTIZE_BIT], int t)
{
    for(int j = 0; j < SMAP_DIM; j++) {
        for(int i = 0; i < (1 << QUANTIZE_BIT); i++) {
            setted_table[t][j][i] = table[j][i];
        }
    }
}

dist_type projected_dist_table_using_setted(unsigned char qpsmap[], int t)
{
	dist_type sum = 0;
	for(int j = 0; j < SMAP_DIM; j++)  {
		sum += setted_table[t][j][qpsmap[j]];
	}
	return sum;
}
*/

// 量子化のための部分距離の事前調査（各射影次元の平均，標準偏差，量子化のオフセット，量子化の刻み幅を求める）．
// ここでは，質問のデータセットのsmap射影像を求めているが，これは残さない．検索時には，別途，個々の質問の射影像を求める（検索コストに含めるため）．
void compute_parmeters_for_qpsmap(smap_pivot_type *pivot, struct_dataset *ds_query, double ave[], double stdev[], double offset[], double slice[])
{
    #ifdef QUANTIZE_MIXED_MOD3
    int range[] = {1 << QUANTIZE_BIT_0, 1 << QUANTIZE_BIT_1, 1 << QUANTIZE_BIT_2};
    #else
    int range = 1 << QUANTIZE_BIT; 
    #endif
    double sum[SMAP_DIM] = {0}, sum2[SMAP_DIM] = {0};
    int num_queries = ds_query->num_data;
    smap_element_type q_smap[SMAP_DIM];

    for(int q = 0; q < num_queries; q++) {
        psmap(q_smap, ds_query->ftr_id[q].ftr, pivot);
        for(int j = 0; j < SMAP_DIM; j++) {
            sum[j] += q_smap[j];
            sum2[j] += q_smap[j] * q_smap[j];
        }
    }

    for(int j = 0; j < SMAP_DIM; j++) {
        ave[j] = sum[j] / num_queries;
        stdev[j] = sqrt(sum2[j] / num_queries - ave[j] * ave[j]);
        #ifdef QUANTIZE_MIXED_MOD3
        slice[j] = stdev[j] * 2 * QUANTIZE_RANGE / range[j % 3];
        #else
        slice[j] = stdev[j] * 2 * QUANTIZE_RANGE / range;
        #endif
        offset[j] = ave[j] - QUANTIZE_RANGE * stdev[j];
        if(j <= 1) {
            fprintf(stderr, "dim = %d, ave = %.2lf, stdev = %.2lf, slice = %.2lf, offset = %.2lf\n", j, ave[j], stdev[j], slice[j], offset[j]);
        }
    }

}

// 質問とデータの量子化射影像間の射影距離計算のための表関数を用意する．
// D~2 用
/*
void make_table_for_query(smap_type smap, unsigned int table[][1 << QUANTIZE_BIT], double offset[], double slice[], int range)
{
    for(int j = 0; j < SMAP_DIM; j++) {
        for(int k = 0; k < range; k++) {
            #ifdef PLUS_HALF
            table[j][k] = pow(abs(smap[j] - ((k + PLUS_HALF) * slice[j] + offset[j])), 2); 
            #else
            table[j][k] = pow(abs(smap[j] - (k * slice[j] + offset[j])), 2); 
            #endif
        }
    }
}
*/

// 質問とデータの量子化射影像間の射影距離計算のための表関数を用意する．
// D~p 用 (p = SCORE_P / 10.0)
#ifndef SCORE_P
#define SCORE_P 10
#endif
void make_table_for_query(smap_type smap, unsigned int table[][1 << QUANTIZE_BIT], double offset[], double slice[])
{
    #ifdef QUANTIZE_MIXED_MOD3
    int range[] = {1 << QUANTIZE_BIT_0, 1 << QUANTIZE_BIT_1, 1 << QUANTIZE_BIT_2};
    #else
    int range = 1 << QUANTIZE_BIT; 
    #endif
    for(int j = 0; j < SMAP_DIM; j++) {
        #ifdef QUANTIZE_MIXED_MOD3
        for(int k = 0; k < range[j % 3]; k++) {
            #ifdef PLUS_HALF
            table[j][k] = pow(abs(smap[j] - ((k + PLUS_HALF) * slice[j] + offset[j])), SCORE_P / 10.0); 
            #else
            table[j][k] = pow(abs(smap[j] - (k * slice[j] + offset[j])), SCORE_P / 10.0); 
            #endif
        }
        #else
        for(int k = 0; k < range; k++) {
            #ifdef PLUS_HALF
            table[j][k] = pow(abs(smap[j] - ((k + PLUS_HALF) * slice[j] + offset[j])), SCORE_P / 10.0); 
            #else
            table[j][k] = pow(abs(smap[j] - (k * slice[j] + offset[j])), SCORE_P / 10.0); 
            #endif
        }
        #endif
    }
}

// 質問とデータの量子化射影像間の射影距離計算のための表関数を用意する．（D~p のpを指定．D~infはp=DBL_MAXを用いる）
/*
void make_table_for_query_p(smap_type smap, unsigned int table[][1 << QUANTIZE_BIT], double offset[], double slice[], int range, double p)
{
    if(p == DBL_MAX) p = 1.0; // D~infはD~1の表を用いる
    for(int j = 0; j < SMAP_DIM; j++) {
        for(int k = 0; k < range; k++) {
            #ifdef PLUS_HALF
            table[j][k] = pow(abs(smap[j] - ((k + PLUS_HALF) * slice[j] + offset[j])), p); 
            #else
            table[j][k] = pow(abs(smap[j] - (k * slice[j] + offset[j])), p); 
            #endif
        }
    }
}
*/

// 質問とデータの量子化射影像間の射影距離計算のための表関数を用意する．（D~p のpを指定．D~infはp=DBL_MAXを用いる）
void make_table_for_query_p(smap_type smap, unsigned int table[][1 << QUANTIZE_BIT], double offset[], double slice[], double p)
{
    #ifdef QUANTIZE_MIXED_MOD3
    int range[] = {1 << QUANTIZE_BIT_0, 1 << QUANTIZE_BIT_1, 1 << QUANTIZE_BIT_2};
    #else
    int range = 1 << QUANTIZE_BIT; 
    #endif
    #ifdef USE_LOWER_BOUND
    unsigned char uchar_q_smap[SMAP_DIM];
    psmap2uchar_qpsmap(smap, uchar_q_smap, offset, slice, range);
//    fprintf(stderr, "set table starts using lower bound\n"); exit(0);
    #endif
    if(p == DBL_MAX) p = 1.0; // D~infはD~1の表を用いる
    for(int j = 0; j < SMAP_DIM; j++) {
        #ifdef QUANTIZE_MIXED_MOD3
        for(int k = 0; k < range[j % 3]; k++) {
            #if defined(USE_LOWER_BOUND)
            if(k < uchar_q_smap[j]) {
                table[j][k] = fabs(smap[j] - ((k + 1) * slice[j] + offset[j])); 
            } else if(k == uchar_q_smap[j]) {
                table[j][k] = 0;
            } else {
                table[j][k] = fabs(smap[j] - (k * slice[j] + offset[j])); 
            }
            #elif defined(PLUS_HALF)
            table[j][k] = pow(abs(smap[j] - ((k + PLUS_HALF) * slice[j] + offset[j])), p); 
            #else
            table[j][k] = pow(abs(smap[j] - (k * slice[j] + offset[j])), p); 
            #endif
        }
        #else
        for(int k = 0; k < range; k++) {
            #if defined(USE_LOWER_BOUND)
            if(k < uchar_q_smap[j]) {
                table[j][k] = fabs(smap[j] - ((k + 1) * slice[j] + offset[j])); 
            } else if(k == uchar_q_smap[j]) {
                table[j][k] = 0;
            } else {
                table[j][k] = fabs(smap[j] - (k * slice[j] + offset[j])); 
            }
            #elif defined(PLUS_HALF)
            table[j][k] = pow(abs(smap[j] - ((k + PLUS_HALF) * slice[j] + offset[j])), p); 
            #else
            table[j][k] = pow(abs(smap[j] - (k * slice[j] + offset[j])), p); 
            #endif
        }
        #endif
    }
}

// 質問とデータの量子化射影像間の射影距離計算のための表関数（非圧縮射影像用の表を元に圧縮表現用の表を作成する）．（D~p を想定．D~inf は除く）
// D~p (p = 1, 1.5, or 2.0 ... ) 用に作成した table を用いれば，packed用もD~p用になる．
// ただし，p = DBL_MAX のときは，この関数は使用してはいけない．*_inf を用いること．
// table: 非圧縮表現のデータのqpsmap用の表関数（make_table_for_queryまたはmake_table_for_query_pを用いて，事前に作成しておく）
// p_table: 圧縮表現のデータのqpsmap用の表関数（出力）
#ifdef QUANTIZE_MIXED_MOD3
#ifdef USE_PACKED_6BIT
void make_table_for_packed_data(unsigned int t[][1 << QUANTIZE_BIT], unsigned int p[][64], double offset[], double slice[])
{
    static int first = 1;
    if(first) {
    fprintf(stderr, "make_table_for_packed_data (1): QUANTIZE_MIXED_MOD3 && USE_PACKED_6BIT\n"); first = 0;
    }

    for(int j = 0; j < 43; j++) {
        for(int k = 0; k < 64; k++) {
            p[j][k] = 0;
        }
    }
    int i, j;
	for(i = j = 0; j + 9 < SMAP_DIM; i += 4, j += 9) {
        //  0 1 2  3 4 5  6 7 8
        // <3-3-2><3-3-2><3-3-2>
        for(int a = 0; a < 8; a++) {
            for(int b = 0; b < 8; b++) {
                // p[0][ab] : t[0][a] + t[1][b]
                p[i][(a << 3) | b] += t[j][a] + t[j + 1][b];
                // p[1][ab] : t[3][a] + t[4][b]
                p[i + 1][(a << 3) | b] += t[j + 3][a] + t[j + 4][b];
                // p[2][ab] : t[6][a] + t[7][b]
                p[i + 2][(a << 3) | b] += t[j + 6][a] + t[j + 7][b];
            }
        }
        for(int a = 0; a < 4; a++) {
            for(int b = 0; b < 4; b++) {
                for(int c = 0; c < 4; c++) {
                    // p[3][ab] : t[2][a] + t[5][b] + t[8][c]
                    p[i + 3][(a << 4) | (b << 2) | c] += t[j + 2][a] + t[j + 5][b] + t[j + 8][c];
                }
            }
        }
	}
    //  0 1 2  3 4 5
    // <3-3-2><3-3-2>
    for(int a = 0; a < 8; a++) {
        for(int b = 0; b < 8; b++) {
            // p[0][ab] : t[0][a] + t[1][b]
            p[i][(a << 3) | b] += t[j][a] + t[j + 1][b];
            // p[1][ab] : t[3][a] + t[4][b]
            p[i + 1][(a << 3) | b] += t[j + 3][a] + t[j + 4][b];
        }
    }
    for(int a = 0; a < 4; a++) {
        for(int b = 0; b < 4; b++) {
            // p[2][ ab] : t[2][a] + t[5][b]
            p[i + 2][(a << 2) | b] += t[j + 2][a] + t[j + 5][b];
        }
    }
}
#else
void make_table_for_packed_data(unsigned int table[][1 << QUANTIZE_BIT], unsigned int p_table[][256], double offset[], double slice[])
{
    static int first = 1;
    if(first) {
    fprintf(stderr, "make_table_for_packed_data (2): QUANTIZE_MIXED_MOD3 && !USE_PACKED_6BIT\n"); first = 0;
    }

    for(int j = 0; j < (SMAP_DIM + 2) / 3; j++) {
        for(int k = 0; k < 256; k++) {
            p_table[j][k] = 0;
        }
    }
    for(int j = 0; j < (SMAP_DIM + 2) / 3; j++) {
        for(int k = 0; k < 256; k++) {
            int idx = (k >> (QUANTIZE_BIT_1 + QUANTIZE_BIT_2)); 
            p_table[j][k] += table[j * 3][idx];
            idx = (k >> QUANTIZE_BIT_2) & ((1 << QUANTIZE_BIT_1) - 1); 
            p_table[j][k] += table[j * 3 + 1][idx];
            idx = k & ((1 << QUANTIZE_BIT_2) - 1); 
            p_table[j][k] += table[j * 3 + 2][idx];
        }
    }
}
#endif
#elif defined(USE_PACKED_3BIT) || defined(USE_PACKED_6BIT)
void make_table_for_packed_data(unsigned int table[][1 << QUANTIZE_BIT], unsigned int p_table[][1 << QUANTIZE_BIT * TABLE_UNIT], double offset[], double slice[])
{
    static int first = 1;
    if(first) {
    fprintf(stderr, "make_table_for_packed_data (3): !QUANTIZE_MIXED_MOD3 && (USE_PACKED_3BIT || USE_PACKED_6BIT)\n"); first = 0;
    }

    int range = 1 << QUANTIZE_BIT; 

    for(int j = 0; j < (SMAP_DIM + TABLE_UNIT - 1) / TABLE_UNIT; j++) {
        for(int k = 0; k < 1 << QUANTIZE_BIT * TABLE_UNIT; k++) {
            p_table[j][k] = 0;
        }
    }
    if(TABLE_UNIT == 1) {
        for(int j = 0; j < SMAP_DIM; j++) {
            for(int k = 0; k < range; k++) {
                p_table[j][k] += table[j][k];
            }
        }
    } else {
        #ifdef TINY_IN_CHAR
        for(int j = 0; j < (SMAP_DIM + TABLE_UNIT - 1) / TABLE_UNIT; j++) {
            for(int k = 0; k < 256; k++) {
                for(int kk = 0; kk < TABLE_UNIT; kk++) {
                    int idx = (k >> ((TABLE_UNIT - kk - 1) * QUANTIZE_BIT)) & ((1 << QUANTIZE_BIT) - 1); 
                    p_table[j][k] += table[j * TABLE_UNIT + kk][idx];
                }
            }
        }
        #else
        for(int j = 0; j < SMAP_DIM; j++) {
            for(int k = 0; k < range; k++) {
                int jj = j / TABLE_UNIT; // table[j][k] を加える p_table[jj]
                int jr = TABLE_UNIT - 1 - j % TABLE_UNIT;
                int mask = ((1 << QUANTIZE_BIT) - 1) << (jr * QUANTIZE_BIT);
                for(int kk = 0; kk < 1 << QUANTIZE_BIT * TABLE_UNIT; kk++) {
                    if(((kk & mask) >> (jr * QUANTIZE_BIT)) == k) {
                        p_table[jj][kk] += table[j][k];
                    }
                }
            }
        }
        #endif
    }
}
#else
void make_table_for_packed_data(unsigned int table[][1 << QUANTIZE_BIT], unsigned int p_table[][1 << QUANTIZE_BIT * TABLE_UNIT], double offset[], double slice[])
{
    static int first = 1;
    if(first) {
    fprintf(stderr, "make_table_for_packed_data (4): !QUANTIZE_MIXED_MOD3 && !USE_PACKED_3BIT && !USE_PACKED_6BIT\n"); first = 0;
    fprintf(stderr, "QUANTIZE_BIT = %d, TABLE_UNIT = %d\n", QUANTIZE_BIT, TABLE_UNIT);
    }

    int range = 1 << QUANTIZE_BIT; 

//    fprintf(stderr, "clear p_table ... ");
    for(int j = 0; j < (SMAP_DIM + TABLE_UNIT - 1) / TABLE_UNIT; j++) {
//        for(int k = 0; k < 256; k++) {
        for(int k = 0; k < (1 << QUANTIZE_BIT * TABLE_UNIT); k++) {
            p_table[j][k] = 0;
        }
    }
//    fprintf(stderr, "OK\n");
    if(TABLE_UNIT == 1) {
        for(int j = 0; j < SMAP_DIM; j++) {
            for(int k = 0; k < range; k++) {
                p_table[j][k] += table[j][k];
            }
        }
    } else {
        #ifdef TINY_IN_CHAR
        int ngrp = (SMAP_DIM + TABLE_UNIT - 1) / TABLE_UNIT;
        // SMAP_DIM 個の次元を TABLE_UNIT 個ずつ束ねたグループを作るときの，グループの個数
        // (例1) SMAP_DIM = 24, TABLE_UNIT = 8 --> ngrp = 3 （ちょうど割り切れるので，どのグループもTABLE_UNIT個の次元が属する）
        // (例2) SMAP_DIM = 26, TABLE_UNIT = 8 --> ngrp = 4 （8個の次元からなるグループは3個，最後のグループは2次元しか属していない）
        int max_k = (1 << QUANTIZE_BIT * TABLE_UNIT);
        // QUANTIXE_BIT ビット整数を TABLE_UNIT 個まとめた整数の個数 2 ^ (QUANTIZE_BIT * TABLE_UNIT)
        // (例1) QUANTIZE_BIT = 2, TABLE_UNIT = 4 --> max_k = 1 << 8 = 256 (8-bit にちょうどいっぱいにできる)
        // (例2) QUANTIZE_BIT = 3, TABLE_UNIT = 2 --> mak_k = 1 << 6 = 64 
        for(int j = 0; j < ngrp; j++) {
            for(int k = 0; k < max_k; k++) {
                for(int kk = 0; kk < TABLE_UNIT; kk++) {
                    if(j * TABLE_UNIT + kk >= SMAP_DIM) {
                        // つぎに p_table に加える射影次元が SMAP_DIM に達したら，打切り．
                        break;
                    }
                    // kk = つぎに加える射影次元がグループ内で左から数えて何番目（0から数える）を表している．
                    int nk = (j + 1) * TABLE_UNIT > SMAP_DIM ? SMAP_DIM - j * TABLE_UNIT : TABLE_UNIT; 
                    // nk = グループに属する射影次元数．端数が出る最後の部分以外では，TABLE_UNIT
                    // (例1) SMAP_DIM = 26, TABLE_UNIT = 8, QUANTIZE_BIT = 1 
                    //       j = 0, 1, 2 --> nk = 8
                    //       j = 3       --> nk = 2 （最後のグループでは端数が生じる．この例では2次元しか残っていない）
                    // (例2) SMAP_DIM = 26, TABLE_UNIT = 3, QUANTIZE_BIT = 2
                    //       j = 0, ... , 7 --> nk = 3
                    //       j = 8          --> nk = 2
                    int kr = nk - kk - 1;
                    // kr = つぎに加える射影次元がグループ内で右から数えて何番目（0から数える）を表している．
                    // (例1) SMAP_DIM = 26, TABLE_UNIT = 8, QUANTIZE_BIT = 1 
                    //       j = 0, 1, 2 --> kr = 8 - kk - 1 (7, 6, 5, ... )
                    //       j = 3       --> kr = 2 - kk - 1 (1, 0)
                    // int idx = (k >> ((TABLE_UNIT - kk - 1) * QUANTIZE_BIT)) & ((1 << QUANTIZE_BIT) - 1); 
                    int idx = (k >> (kr * QUANTIZE_BIT)) & ((1 << QUANTIZE_BIT) - 1); 
                    // idx = QUANTIZE_BIT 整数を nk 個まとめた整数 k から，右から kr 番目を取り出したもの
                    p_table[j][k] += table[j * TABLE_UNIT + kk][idx];
                }
            }
        }
        #else
        for(int j = 0; j < SMAP_DIM; j++) {
            for(int k = 0; k < range; k++) {
                int jj = j / TABLE_UNIT; // table[j][k] を加える p_table[jj]
                int jr = TABLE_UNIT - 1 - j % TABLE_UNIT;
                int mask = ((1 << QUANTIZE_BIT) - 1) << (jr * QUANTIZE_BIT);
                for(int kk = 0; kk < 256; kk++) {
                    if(((kk & mask) >> (jr * QUANTIZE_BIT)) == k) {
                        p_table[jj][kk] += table[j][k];
                    }
                }
            }
        }
        #endif
    }
}
#endif
// D~inf 用（table は D~1 用に作成したものを用いる）
void make_table_for_packed_data_inf(unsigned int table[][1 << QUANTIZE_BIT], unsigned int p_table[][256], double offset[], double slice[])
{
    int range = 1 << QUANTIZE_BIT; 
    for(int j = 0; j < SMAP_DIM / TABLE_UNIT; j++) {
        for(int k = 0; k < 256; k++) {
            p_table[j][k] = 0;
        }
    }
    for(int j = 0; j < SMAP_DIM; j++) {
        for(int k = 0; k < range; k++) {
        	int jj = j / TABLE_UNIT; // table[j][k] が関係するのは p_table[jj]
            int jr = TABLE_UNIT - 1 - j % TABLE_UNIT; // QUANTIZE_BIT;
            int mask = ((1 << QUANTIZE_BIT) - 1) << (jr * QUANTIZE_BIT);
            for(int kk = 0; kk < 256; kk++) {
                if(((kk & mask) >> (jr * QUANTIZE_BIT)) == k) {
                    if(p_table[jj][kk] < table[j][k]) {
                        p_table[jj][kk] = table[j][k];
                    }
                }
            }
        }
    }
}

// 質問の射影距離計算tableを用いたucharの量子化射影像間の射影距離(D~p) (p乗根は取らない)
// tableに設定されたものをそのまま用いた総和を求める（D~p共通．ただし，D~infは除く）
dist_type projected_dist_table(unsigned int table[][1 << QUANTIZE_BIT], unsigned char qpsmap[])
{
	dist_type sum = 0;
	for(int j = 0; j < SMAP_DIM; j++)  {
		sum += table[j][qpsmap[j]];
//        printf("table[%d][%u], %d, sum, %d\n", j, qpsmap[j], table[j][qpsmap[j]], sum);
	}
//    printf("sum = %u\n", sum); getchar();
	return sum;
}

// 質問の射影距離計算tableを用いたucharの量子化射影像間の射影距離(Lp) (p乗根は取らない)
// D~inf 用（table は D~1 用に作成したものを用いる）
dist_type projected_dist_table_inf(unsigned int table[][1 << QUANTIZE_BIT], unsigned char qpsmap[])
{
	unsigned int max = 0;
	for(int j = 0; j < SMAP_DIM; j++)  {
        if(max < table[j][qpsmap[j]]) max = table[j][qpsmap[j]];
	}
	return max;
}

// 質問の射影距離計算のための表 packed_table を用いた圧縮表現されたデータのqpsmapの射影距離(D~p) (p乗根は取らない)
// tableに設定されたものをそのまま用いた総和を求める（D~p共通．ただし，D~infは除く）
/*
dist_type projected_dist_packed_table(unsigned int p_table[][256], tiny_int *qpsmap)
{
	dist_type sum = 0;
	for(int j = 0; j < SMAP_DIM * QUANTIZE_BIT / 32 ; j++) {
    	sum += (p_table[j * 4][qpsmap[j] >> 24] + p_table[j * 4 + 1][(qpsmap[j] >> 16) & 0xff] +p_table[j * 4 + 2][(qpsmap[j] >> 8) & 0xff] +p_table[j * 4 + 3][qpsmap[j] & 0xff]);
	}
	return sum;
}
*/
#ifdef TINY_IN_CHAR
#ifdef QUANTIZE_MIXED_MOD3
#ifdef USE_PACKED_6BIT
dist_type projected_dist_packed_table(unsigned int t[][64], tiny_int *q)
{
//    static int first = 1;
//    if(first) {
//    fprintf(stderr, "projected_dist_packed_table (1): TINY_IN_CHAR, QUANTIZE_MIXED_MOD3, USE_PACKED_6BIT\n"); first = 0;
//    }
	dist_type sum = 0;
    int i, j;
	for(i = 0, j = 0; j + 3 < (SMAP_DIM + 2) / 3; i += 4, j += 3) {
    // <3-bit(a)><3-bit(b)><2-bit(g)> <3-bit(c)><3-bit(d)><2-bit(h)> <3-bit(e)><3-bit(f)><2-bit(i)> 
        sum += t[i][q[j] >> 2] + t[i + 1][q[j + 1] >> 2] + t[i + 2][q[j + 2] >> 2] + t[i + 3][((q[j] & 3) << 4) | ((q[j + 1] & 3) << 2) | (q[j + 2] & 3)];
	}
    // <3-bit(a)><3-bit(b)><2-bit(g)> <3-bit(c)><3-bit(d)><2-bit(h)>  
        sum += t[i][q[j] >> 2] + t[i + 1][q[j + 1] >> 2] + t[i + 2][((q[j] & 3) << 2) | (q[j + 1] & 3)];
    return sum;
}
#else
dist_type projected_dist_packed_table(unsigned int p_table[][256], tiny_int *qpsmap)
{
//    static int first = 1;
//    if(first) {
//    fprintf(stderr, "projected_dist_packed_table (2): TINY_IN_CHAR, QUANTIZE_MIXED_MOD3, !USE_PACKED_6BIT\n"); first = 0;
//    }
	dist_type sum = 0;
	for(int j = 0; j < (SMAP_DIM + 2) / 3; j++) {
		sum += p_table[j][qpsmap[j]];
	}
    return sum;
}
#endif
#else
dist_type projected_dist_packed_table(unsigned int p_table[][1 << QUANTIZE_BIT * TABLE_UNIT], tiny_int *qpsmap)
// dist_type projected_dist_packed_table(unsigned int p_table[][256], tiny_int *qpsmap)
{
//    static int first = 1;
//    if(first) {
//    fprintf(stderr, "projected_dist_packed_table (3): TINY_IN_CHAR, !QUANTIZE_MIXED_MOD3"); first = 0;
//    fprintf(stderr, ", QUANTIZE_BIT = %d, TABLE_UNIT = %d, sizeof(tiny_int) = %ld\n", QUANTIZE_BIT, TABLE_UNIT, sizeof(tiny_int));
//    }
	dist_type sum = 0;
	for(int j = 0; j < (SMAP_DIM + TABLE_UNIT - 1)/ (TABLE_UNIT); j++) {
		sum += p_table[j][qpsmap[j]];
	}
    return sum;
}
#endif
#elif defined(USE_PACKED_3BIT) || defined(USE_PACKED_6BIT) || defined(USE_PACKED_4BIT)

/*
dist_type projected_dist_packed_table(unsigned int p_table[][1 << QUANTIZE_BIT * TABLE_UNIT], tiny_int *qpsmap)
{
// mask = (1 << length) - 1 = 00000011111 （下位のlength-bitだけが1，それ以外が0のビットパターン）
// （例）残りがあるとき
// length = 9 = (QUANTIZE_BIT(3) * TABLE_UNIT(3)), leftover = 5  
//     <- 5-bit ->     <- 4-bit -><- 9-bit -><- 9-bit -><- 9-bit -><- 1-bit ->    
//     xxxxxxxxxxx     yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww???????????XXXXXXXXXXX
// 最初の見出し = xxxxxxxxxxyyyyyyyyy = (xxxxxxxxxxx << 4) | (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 4))
// つぎの見出し = zzzzzzzzzzz = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> ((32 - 4) - 9)) & mask
// つぎの見出し = wwwwwwwwwww = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> ((32 - 4) - 18)) & mask
// つぎの見出し = ??????????? = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> ((32 - 4) - 27)) & mask
// つぎの残りビット数 = 1 = (32 - 4 - 27) 
// つぎの残りビットは = XXXXXXXXXXX = yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... & 0000001 (下位の残りビット）

// （例）残りがないとき
// length = 9 = (QUANTIZE_BIT(3) * TABLE_UNIT(3)), leftover = 0  
//     <- 0-bit ->     <- 9-bit -><- 9-bit -><- 9-bit -><- 5-bit ->    
//     xxxxxxxxxxx     yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwwwXXXXXXXXXXX
// 最初の見出し = yyyyyyyyyyy = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 9)
// つぎの見出し = zzzzzzzzzzz = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 18) & mask
// つぎの見出し = wwwwwwwwwww = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 27) & mask
// 残りビット数 = 5 = (32 - 27) 
// 残りビットは = XXXXXXXXXXX = yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... & 0011111 (下位の残りビット）


	dist_type sum = 0;
	int leftover = 0; // 前の32-bitを処理した残り（値．右寄せ）
	int leftover_bit = 0; // 前の32-bitを処理した残り（ビット数）
	int length = QUANTIZE_BIT * TABLE_UNIT; // 表見出しのビット数（表を引くときに用いる整数のビット長）
	int mask = (1 << length) - 1; // 下位のQUANTIZE_BITを取り出すときに用いるマスク（これとビット毎の＆を取る）
	
	for(int j = 0, i = 0, m = 0; (j < (SMAP_DIM * QUANTIZE_BIT + 31) / 32) && (m < SMAP_DIM); j++) {
		int shift = length - leftover_bit;
		int idx = (leftover << shift) | (qpsmap[j] >> (32 - shift));
		sum += p_table[i++][idx];
        m += TABLE_UNIT;
        shift = 32 - shift;
		while(m < SMAP_DIM && shift >= length) {
            m += TABLE_UNIT;
			shift -= length;
			idx = (qpsmap[j] >> shift) & mask;
            if(m > SMAP_DIM) {
                idx >>= (m - SMAP_DIM) * QUANTIZE_BIT; // 右シフトで余分な部分を消して右寄せにする．
            }
			sum += p_table[i++][idx];
		}
		leftover_bit = shift;
		leftover = qpsmap[j] & ((1 << leftover_bit) - 1);
	}
    return sum;
}
*/
dist_type projected_dist_packed_table(unsigned int p_table[][1 << QUANTIZE_BIT * TABLE_UNIT], tiny_int *qpsmap)
{
// mask = (1 << length) - 1 = 00000011111 （下位のlength-bitだけが1，それ以外が0のビットパターン）
// （例）残りがあるとき
// length = 9 = (QUANTIZE_BIT(3) * TABLE_UNIT(3)), leftover = 5  
//     <- 5-bit ->     <- 4-bit -><- 9-bit -><- 9-bit -><- 9-bit -><- 1-bit ->    
//     xxxxxxxxxxx     yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww???????????XXXXXXXXXXX
// 最初の見出し = xxxxxxxxxxyyyyyyyyy = (xxxxxxxxxxx << 4) | (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 4))
// つぎの見出し = zzzzzzzzzzz = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> ((32 - 4) - 9)) & mask
// つぎの見出し = wwwwwwwwwww = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> ((32 - 4) - 18)) & mask
// つぎの見出し = ??????????? = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> ((32 - 4) - 27)) & mask
// つぎの残りビット数 = 1 = (32 - 4 - 27) 
// つぎの残りビットは = XXXXXXXXXXX = yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... & 0000001 (下位の残りビット）

// 最初の 1-word = 32-bit
// length = 9 = (QUANTIZE_BIT(3) * TABLE_UNIT(3)), leftover = 0  
//     <- 0-bit ->     <- 9-bit -><- 9-bit -><- 9-bit ->            <- 5-bit -> 9次元分    
//     xxxxxxxxxxx     yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwwwXXXXXXXXXXX
// 最初の見出し = yyyyyyyyyyy = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 9)
// つぎの見出し = zzzzzzzzzzz = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 18) & mask
// つぎの見出し = wwwwwwwwwww = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 27) & mask
// 残りビット数 = 5 = (32 - 27) 
// 残りビットは = XXXXXXXXXXX = yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... & 0011111 (下位の残りビット）
//     <- 5-bit ->     <- 4-bit -><- 9-bit -><- 9-bit -><- 9-bit -><- 1-bit -> 3 + 9 = 12 次元分 (21次元)
//     <- 1-bit ->     <- 8-bit -><- 9-bit -><- 9-bit ->           <- 6-bit -> 3 + 6 =  9 次元分 (30次元)
//     <- 6-bit ->     <- 3-bit -><- 9-bit -><- 9-bit -><- 9-bit -><- 2-bit -> 3 + 9 = 12 次元分 (42次元)
//     <- 2-bit ->     <- 7-bit -><- 9-bit -><- 9-bit ->           <- 7-bit -> 3 + 6 =  9 次元分 (51次元)
//     <- 7-bit ->     <- 2-bit -><- 9-bit -><- 9-bit -><- 9-bit -><- 3-BIT -> 3 + 9 = 12 次元分 (63次元)
//     <- 3-bit ->     <- 7-bit -><- 9-bit -><- 9-bit ->           <- 7-bit -> 3 + 6 =  9 次元分 (72次元)
//     <- 8-bit ->     <- 2-bit -><- 9-bit -><- 9-bit -><- 9-bit -><- 4-BIT -> 3 + 9 = 12 次元分 (84次元)
//     <- 4-bit ->     <- 7-bit -><- 9-bit -><- 9-bit -><- 9-bit ->            3 + 6 = 12 次元分 (96次元)

#ifdef USE_PACKED_3BIT
//    static int first = 1;
//    if(first) {
//    fprintf(stderr, "projected_dist_packed_table (5): !TINY_IN_CHAR, USE_PACKED_3BIT\n"); first = 0;
//    }
// 96次元3-bitで9-bitの表専用
	dist_type sum = 0;
    tiny_int *q = qpsmap;
    unsigned int (*t)[1 << QUANTIZE_BIT * TABLE_UNIT] = p_table;
    unsigned int m = (1 << 9) - 1; // 111 111 111b
//     <- 0-bit ->     <- 9-bit -><- 9-bit -><- 9-bit ->            <- 5-bit -> 9次元分    
    sum +=                                          t[0][q[0] >> 23] + t[1][(q[0] >> 14) & m] + t[2][(q[0] >> 5) & m];
//     <- 5-bit ->     <- 4-bit -><- 9-bit -><- 9-bit -><- 9-bit -><- 1-bit -> 3 + 9 = 12 次元分 (21次元)
    sum += t[3][((q[0] << 4) & m) | (q[1] >> 28)] + t[4][(q[1] >> 19) & m] + t[5][(q[1] >> 10) & m] + t[6][(q[1] >> 1) & m];
//     <- 1-bit ->     <- 8-bit -><- 9-bit -><- 9-bit ->           <- 6-bit -> 3 + 6 =  9 次元分 (30次元)
    sum += t[7][((q[1] << 8) & m) | (q[2] >> 24)] + t[8][(q[2] >> 15) & m] + t[9][(q[2] >>  6) & m];
//     <- 6-bit ->     <- 3-bit -><- 9-bit -><- 9-bit -><- 9-bit -><- 2-bit -> 3 + 9 = 12 次元分 (42次元)
    sum += t[10][((q[2] << 3) & m) | (q[3] >> 29)] + t[11][(q[3] >> 20) & m] + t[12][(q[3] >> 11) & m] + t[13][(q[3] >> 2) & m];
//     <- 2-bit ->     <- 7-bit -><- 9-bit -><- 9-bit ->           <- 7-bit -> 3 + 6 =  9 次元分 (51次元)
    sum += t[14][((q[3] << 7) & m) | (q[4] >> 25)] + t[15][(q[4] >> 16) & m] + t[16][(q[4] >>  7) & m];
//     <- 7-bit ->     <- 2-bit -><- 9-bit -><- 9-bit -><- 9-bit -><- 3-BIT -> 3 + 9 = 12 次元分 (63次元)
    sum += t[17][((q[4] << 2) & m) | (q[5] >> 30)] + t[18][(q[5] >> 21) & m] + t[19][(q[5] >> 12) & m] + t[20][(q[5] >> 3) & m];
//     <- 3-bit ->     <- 6-bit -><- 9-bit -><- 9-bit ->           <- 8-bit -> 3 + 6 =  9 次元分 (72次元)
    sum += t[21][((q[5] << 6) & m) | (q[6] >> 26)] + t[22][(q[6] >> 17) & m] + t[23][(q[6] >>  8) & m];
//     <- 8-bit ->     <- 1-bit -><- 9-bit -><- 9-bit -><- 9-bit -><- 4-BIT -> 3 + 9 = 12 次元分 (84次元)
    sum += t[24][((q[6] << 1) & m) | (q[7] >> 31)] + t[25][(q[7] >> 22) & m] + t[26][(q[7] >> 13) & m] + t[27][(q[7] >> 4) & m];
//     <- 4-bit ->     <- 5-bit -><- 9-bit -><- 9-bit -><- 9-bit ->            3 + 9 = 12 次元分 (96次元)
    sum += t[28][((q[7] << 5) & m) | (q[8] >> 27)] + t[29][(q[8] >> 18) & m] + t[30][(q[8] >>  9) & m] + t[31][q[8] & m];
#elif defined(USE_PACKED_6BIT)
    static int first = 1;
    if(first) {
    fprintf(stderr, "projected_dist_packed_table (6-2): !TINY_IN_CHAR, USE_PACKED_6BIT\n"); first = 0;
    }
// 96次元3-bitで6-bitの表専用
	dist_type sum = 0;
    tiny_int *q = qpsmap;
    unsigned int (*t)[1 << QUANTIZE_BIT * TABLE_UNIT] = p_table;
    unsigned int m = (1 << QUANTIZE_BIT * TABLE_UNIT) - 1; // 111 111b
    int rest_bit = 0;
    unsigned int rest = 0;
    int dim;
    for(dim = 0; dim < SMAP_DIM; q++) {
//        printf("dim = %d: ", dim); print_bin(q[0]);
        if((dim <= SMAP_DIM - 10 && rest_bit != 4) || (dim <= SMAP_DIM - 12 && rest_bit == 4)) {
//            printf(": rest_bit = %d: \n", rest_bit);
            switch(rest_bit) {
                case 0:
                    //  <- 0-bit ->     <- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->            <- 2-bit -> 10次元分
//                    print_bin(q[0] >> 26); printf(" -> add %d\n", t[0][q[0] >> 26]);
//                    print_bin((q[0] >> 20) & m); printf(" -> add %d\n", t[1][(q[0] >> 20) & m]);
//                    print_bin((q[0] >> 14) & m); printf(" -> add %d\n", t[2][(q[0] >> 14) & m]);
//                    print_bin((q[0] >> 8) & m); printf(" -> add %d\n", t[3][(q[0] >> 8) & m]);
//                    print_bin((q[0] >> 2) & m); printf(" -> add %d\n", t[4][(q[0] >> 2) & m]);
//                    printf("rest = "); print_bin(q[0] & 3); printf(" -> "); print_bin((q[0] << 4) & m); printf("\n");
                    sum += t[0][q[0] >> 26] + t[1][(q[0] >> 20) & m] + t[2][(q[0] >> 14) & m] + t[3][(q[0] >>  8) & m] + t[4][(q[0] >>  2) & m];
                    t += 5;
                    rest_bit = 2;
                    rest = (q[0] << 4) & m; // 4-bit left shifted of rest 2-bit
                    dim += 10;
                    break;
                case 2:
                    //  <- 2-bit ->     <- 4-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->            <- 4-bit -> 10次元分    
//                    print_bin(rest | (q[0] >> 28)); printf(" -> add %d\n", t[0][rest | (q[0] >> 28)]);
//                    print_bin((q[0] >> 22) & m); printf(" -> add %d\n", t[1][(q[0] >> 22) & m]);
//                    print_bin((q[0] >> 16) & m); printf(" -> add %d\n", t[2][(q[0] >> 16) & m]);
//                    print_bin((q[0] >> 10) & m); printf(" -> add %d\n", t[3][(q[0] >> 10) & m]);
//                    print_bin((q[0] >> 4) & m); printf(" -> add %d\n", t[4][(q[0] >> 4) & m]);
//                    printf("rest = "); print_bin(q[0] & 15); printf(" -> "); print_bin((q[0] << 2) & m); printf("\n");
                    sum += t[0][rest | (q[0] >> 28)] 
                        + t[1][(q[0] >> 22) & m] + t[2][(q[0] >> 16) & m] + t[3][(q[0] >> 10) & m] + t[4][(q[0] >>  4) & m];
                    t += 5;
                    rest_bit = 4;
                    rest = (q[0] << 2) & m; // 2-bit left shifted of rest 4-bit
                    dim += 10;
                    break;
                case 4:
                    //  <- 4-bit ->     <- 2-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->             12次元分    
//                    print_bin(rest | (q[0] >> 30)); printf(" -> add %d\n", t[0][rest | (q[0] >> 30)]);
//                    print_bin((q[0] >> 24) & m); printf(" -> add %d\n", t[1][(q[0] >> 24) & m]);
//                    print_bin((q[0] >> 19) & m); printf(" -> add %d\n", t[2][(q[0] >> 18) & m]);
//                    print_bin((q[0] >> 12) & m); printf(" -> add %d\n", t[3][(q[0] >> 12) & m]);
//                    print_bin((q[0] >> 6) & m); printf(" -> add %d\n", t[4][(q[0] >> 6) & m]);
//                    print_bin(q[0] & m); printf(" -> add %d\n", t[5][q[0] & m]);
//                    printf("rest = "); print_bin(q[0] & 15); printf(" -> "); print_bin((q[0] << 2) & m); printf("\n");
                    sum += t[0][rest | (q[0] >> 30)] 
                        + t[1][(q[0] >> 24) & m] + t[2][(q[0] >> 18) & m] + t[3][(q[0] >> 12) & m] + t[4][(q[0] >>  6) & m] + t[5][q[0] & m];
                    t += 6;
                    rest_bit = 0;
                    rest = 0;
                    dim += 12;
                    break;
            }
        } else {
            // 最後に 32-bit を使い切らないとき（デバッグ不十分）
            int sp;
            switch(rest_bit) {
                case 0:
                    //  <- 0-bit ->     <- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->            <- 2-bit -> 10次元分
                    sp = 26;
                    while(dim < SMAP_DIM) {
                        sum += t[0][(q[0] >> sp) & m];
                        dim++;
                        sp -= 6;
                    }   
                    break;
                case 2:
                    //  <- 2-bit ->     <- 4-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->            <- 4-bit -> 10次元分    
                    sum += t[0][rest | (q[0] >> 28)];
                    dim += 2;
                    sp = 22;
                    while(dim < SMAP_DIM) {
                        sum += t[0][(q[0] >> sp) & m];
                        dim += 2;
                        sp -= 6;
                    }
                    break;
                case 4:
                    //  <- 4-bit ->     <- 2-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->             12次元分    
                    sum += t[0][rest | (q[0] >> 30)];
                    dim += 2;
                    sp = 24;
                    while(dim < SMAP_DIM) {
                        sum += t[0][(q[0] >> sp) & m];
                        dim += 2;
                        sp -= 6;
                    }
                    break;
            }

        }
    }
/*
    static int first = 1;
    if(first) {
    fprintf(stderr, "projected_dist_packed_table (6): !TINY_IN_CHAR, USE_PACKED_6BIT\n"); first = 0;
    }
// 96次元3-bitで6-bitの表専用
	dist_type sum = 0;
    tiny_int *q = qpsmap;
    unsigned int (*t)[1 << QUANTIZE_BIT * TABLE_UNIT] = p_table;
    unsigned int m = (1 << QUANTIZE_BIT * TABLE_UNIT) - 1; // 111 111b
//     <- 0-bit ->     <- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->            <- 2-bit -> 10次元分    
    sum +=             t[0][q[0] >> 26] + t[1][(q[0] >> 20) & m] + t[2][(q[0] >> 14) & m] + t[3][(q[0] >>  8) & m] + t[4][(q[0] >>  2) & m];
//     <- 2-bit ->     <- 4-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->            <- 4-bit -> 10次元分    
    sum += t[5][((q[0] << 4) & m) | (q[1] >> 28)] 
                     + t[6][(q[1] >> 22) & m] + t[7][(q[1] >> 16) & m] + t[8][(q[1] >> 10) & m] + t[9][(q[1] >>  4) & m];
//     <- 4-bit ->     <- 2-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->             12次元分    
    sum += t[10][((q[1] << 2) & m) | (q[2] >> 30)] 
                     + t[11][(q[2] >> 24) & m] + t[12][(q[2] >> 18) & m] + t[13][(q[2] >> 12) & m] + t[14][(q[2] >>  6) & m] + t[15][q[2] & m];

//     <- 0-bit ->     <- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->            <- 2-bit -> 10次元分    
    sum +=             t[16][q[3] >> 26] + t[17][(q[3] >> 20) & m] + t[18][(q[3] >> 14) & m] + t[19][(q[3] >>  8) & m] + t[20][(q[3] >>  2) & m];
//     <- 2-bit ->     <- 4-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->            <- 4-bit -> 10次元分    
    sum += t[21][((q[3] << 4) & m) | (q[4] >> 28)] 
                     + t[22][(q[4] >> 22) & m] + t[23][(q[4] >> 16) & m] + t[24][(q[4] >> 10) & m] + t[25][(q[4] >>  4) & m];
//     <- 4-bit ->     <- 2-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->             12次元分    
    sum += t[26][((q[4] << 2) & m) | (q[5] >> 30)] 
                     + t[27][(q[5] >> 24) & m] + t[28][(q[5] >> 18) & m] + t[29][(q[5] >> 12) & m] + t[30][(q[5] >>  6) & m] + t[31][q[5] & m];
//     <- 0-bit ->     <- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->            <- 2-bit -> 10次元分    
    sum +=             t[32][q[6] >> 26] + t[33][(q[6] >> 20) & m] + t[34][(q[6] >> 14) & m] + t[35][(q[6] >>  8) & m] + t[36][(q[6] >>  2) & m];
//     <- 2-bit ->     <- 4-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->            <- 4-bit -> 10次元分    
    sum += t[37][((q[6] << 4) & m) | (q[7] >> 28)] 
                     + t[38][(q[7] >> 22) & m] + t[39][(q[7] >> 16) & m] + t[40][(q[7] >> 10) & m] + t[41][(q[7] >>  4) & m];
//     <- 4-bit ->     <- 2-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit -><- 6-bit ->             12次元分    
    sum += t[42][((q[7] << 2) & m) | (q[8] >> 30)] 
                     + t[43][(q[8] >> 24) & m] + t[44][(q[8] >> 18) & m] + t[45][(q[8] >> 12) & m] + t[46][(q[8] >>  6) & m] + t[47][q[8] & m];
*/
#else 
//    static int first = 1;
//    if(first) {
//    fprintf(stderr, "projected_dist_packed_table (8): !TINY_IN_CHAR, USE_PACKED_4BIT"); first = 0;
//    fprintf(stderr, ", QUANTIZE_BIT = %d, TABLE_UNIT = %d, sizeof(tiny_int) = %ld\n", QUANTIZE_BIT, TABLE_UNIT, sizeof(tiny_int));
//    }
// 96次元4-bitで8-bitの表専用
	dist_type sum = 0;
//    tiny_int *q = qpsmap;
    unsigned char *q = (unsigned char *)qpsmap;
    unsigned int (*t)[1 << QUANTIZE_BIT * TABLE_UNIT] = p_table;
    unsigned int m = (1 << 8) - 1; // 1111 1111b
//         <- 8-bit ->        <- 8-bit ->              <- 8-bit ->             <- 8-bit -> 8次元分    
    sum = t[ 0][q [ 3]] + t[ 1][q[ 2]] + t[ 2][q[ 1]] + t[ 3][q[ 0]] +
          t[ 4][q [ 7]] + t[ 5][q[ 6]] + t[ 6][q[ 5]] + t[ 7][q[ 4]] +
          t[ 8][q [11]] + t[ 9][q[10]] + t[10][q[ 9]] + t[11][q[ 8]] +
          t[12][q [15]] + t[13][q[14]] + t[14][q[13]] + t[15][q[12]] +
          t[16][q [19]] + t[17][q[18]] + t[18][q[17]] + t[19][q[16]] +
          t[20][q [23]] + t[21][q[22]] + t[22][q[21]] + t[23][q[20]] +
          t[24][q [27]] + t[25][q[26]] + t[26][q[25]] + t[27][q[24]] +
          t[28][q [31]] + t[29][q[30]] + t[30][q[29]] + t[31][q[28]] +
          t[32][q [35]] + t[33][q[34]] + t[34][q[33]] + t[35][q[32]] +
          t[36][q [39]] + t[37][q[38]] + t[38][q[37]] + t[39][q[36]] +
          t[40][q [43]] + t[41][q[42]] + t[42][q[41]] + t[43][q[40]] +
          t[44][q [47]] + t[45][q[46]] + t[46][q[45]] + t[47][q[44]];
/*
    sum = t[0][q[0] >> 24] + t[1][(q[0] >> 16) & m] + t[2][(q[0] >> 8) & m] + t[3][q[0] & m] +
          t[4][q[1] >> 24] + t[5][(q[1] >> 16) & m] + t[6][(q[1] >> 8) & m] + t[7][q[1] & m] +
          t[8][q[2] >> 24] + t[9][(q[2] >> 16) & m] + t[10][(q[2] >> 8) & m] + t[11][q[2] & m] +
          t[12][q[3] >> 24] + t[13][(q[3] >> 16) & m] + t[14][(q[3] >> 8) & m] + t[15][q[3] & m] +
          t[16][q[4] >> 24] + t[17][(q[4] >> 16) & m] + t[18][(q[4] >> 8) & m] + t[19][q[4] & m] +
          t[20][q[5] >> 24] + t[21][(q[5] >> 16) & m] + t[22][(q[5] >> 8) & m] + t[23][q[5] & m] +
          t[24][q[6] >> 24] + t[25][(q[6] >> 16) & m] + t[26][(q[6] >> 8) & m] + t[27][q[6] & m] +
          t[28][q[7] >> 24] + t[29][(q[7] >> 16) & m] + t[30][(q[7] >> 8) & m] + t[31][q[7] & m] +
          t[32][q[8] >> 24] + t[33][(q[8] >> 16) & m] + t[34][(q[8] >> 8) & m] + t[35][q[8] & m] +
          t[36][q[9] >> 24] + t[37][(q[9] >> 16) & m] + t[38][(q[9] >> 8) & m] + t[39][q[9] & m] +
          t[40][q[10] >> 24] + t[41][(q[10] >> 16) & m] + t[42][(q[10] >> 8) & m] + t[43][q[10] & m] +
          t[44][q[11] >> 24] + t[45][(q[11] >> 16) & m] + t[46][(q[11] >> 8) & m] + t[47][q[11] & m];
*/
#endif
    return sum;
}

#else
dist_type projected_dist_packed_table(unsigned int p_table[][256], tiny_int *qpsmap)
{
// mask = (1 << length) - 1 = 00000011111 （下位のlength-bitだけが1，それ以外が0のビットパターン）
// （例）残りがあるとき
// length = 9 = (QUANTIZE_BIT(3) * TABLE_UNIT(3)), leftover = 5  
//     <- 5-bit ->     <- 4-bit -><- 9-bit -><- 9-bit -><- 9-bit -><- 1-bit ->    
//     xxxxxxxxxxx     yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww???????????XXXXXXXXXXX
// 最初の見出し = xxxxxxxxxxyyyyyyyyy = (xxxxxxxxxxx << 4) | (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 4))
// つぎの見出し = zzzzzzzzzzz = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> ((32 - 4) - 9)) & mask
// つぎの見出し = wwwwwwwwwww = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> ((32 - 4) - 18)) & mask
// つぎの見出し = ??????????? = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> ((32 - 4) - 27)) & mask
// つぎの残りビット数 = 1 = (32 - 4 - 27) 
// つぎの残りビットは = XXXXXXXXXXX = yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... & 0000001 (下位の残りビット）

// （例）残りがないとき
// length = 9 = (QUANTIZE_BIT(3) * TABLE_UNIT(3)), leftover = 0  
//     <- 0-bit ->     <- 9-bit -><- 9-bit -><- 9-bit -><- 5-bit ->    
//     xxxxxxxxxxx     yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwwwXXXXXXXXXXX
// 最初の見出し = yyyyyyyyyyy = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 9)
// つぎの見出し = zzzzzzzzzzz = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 18) & mask
// つぎの見出し = wwwwwwwwwww = (yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... >> (32 - 27) & mask
// 残りビット数 = 5 = (32 - 27) 
// 残りビットは = XXXXXXXXXXX = yyyyyyyyyyyzzzzzzzzzzzwwwwwwwwwww... & 0011111 (下位の残りビット）

//    static int first = 1;
//    if(first) {
//    fprintf(stderr, "projected_dist_packed_table (7): !TINY_IN_CHAR, !USE_PACKED_3BIT, !USE_PACKED_6BIT"); first = 0;
//    fprintf(stderr, ", QUANTIZE_BIT = %d, TABLE_UNIT = %d, sizeof(tiny_int) = %ld\n", QUANTIZE_BIT, TABLE_UNIT, sizeof(tiny_int));
//    }
	dist_type sum = 0;
	int leftover = 0; // 前の32-bitを処理した残り（値．右寄せ）
	int leftover_bit = 0; // 前の32-bitを処理した残り（ビット数）
	int length = QUANTIZE_BIT * TABLE_UNIT; // 表見出しのビット数（表を引くときに用いる整数のビット長）
	int mask = (1 << length) - 1; // 下位のQUANTIZE_BITを取り出すときに用いるマスク（これとビット毎の＆を取る）
	
	for(int j = 0, i = 0, m = 0; (j < (SMAP_DIM * QUANTIZE_BIT + 31) / 32) && (m < SMAP_DIM); j++) {
		int shift = length - leftover_bit;
		int idx = (leftover << shift) | (qpsmap[j] >> (32 - shift));
		sum += p_table[i++][idx];
        m += TABLE_UNIT;
        shift = 32 - shift;
		while(m < SMAP_DIM && shift >= length) {
            m += TABLE_UNIT;
			shift -= length;
			idx = (qpsmap[j] >> shift) & mask;
            if(m > SMAP_DIM) {
                idx >>= (m - SMAP_DIM) * QUANTIZE_BIT; // 右シフトで余分な部分を消して右寄せにする．
            }
			sum += p_table[i++][idx];
		}
		leftover_bit = shift;
		leftover = qpsmap[j] & ((1 << leftover_bit) - 1);
	}
    return sum;
}
#endif

// D~inf用．p_tableはD~inf用のものを用いる．
dist_type projected_dist_packed_table_inf(unsigned int p_table[][256], tiny_int *qpsmap)
{
	unsigned int max = 0;
	for(int j = 0; j < SMAP_DIM * QUANTIZE_BIT / 32 ; j++) {
    	if(max < p_table[j * 4][qpsmap[j] >> 24])               max = p_table[j * 4][qpsmap[j] >> 24];
        if(max < p_table[j * 4 + 1][(qpsmap[j] >> 16) & 0xff])  max = p_table[j * 4 + 1][(qpsmap[j] >> 16) & 0xff];
        if(max < p_table[j * 4 + 2][(qpsmap[j] >> 8) & 0xff])   max = p_table[j * 4 + 2][(qpsmap[j] >> 8) & 0xff];
        if(max < p_table[j * 4 + 3][qpsmap[j] & 0xff])          max = p_table[j * 4 + 3][qpsmap[j] & 0xff];
	}
	return max;
}

#endif

// 特徴データftrの射影像のdim次元を求める．
smap_element_type psmap_dim(int dim, ftr_type ftr, smap_pivot_type *pivot) {
//    #ifdef SMAP_PARTITION_TYPE_PQBP
    #if defined(SMAP_PARTITION_TYPE_PQBP) || defined(PARTITION_TYPE_PQBP) || defined(PARTITION_TYPE_QBP)
//    fprintf(stderr, "dim = %d, start = %d, part_dim = %d\n",dim, SMAP_PART_START(dim), SMAP_PART_DIM(dim));
//    return partition_dist_L2(pivot->p[dim], ftr, SMAP_PART_START(dim), SMAP_PART_DIM(dim));
    #ifdef USE_PD
    return partition_dist_pivot_L2(ftr, pivot->p[dim], pivot->num_axis[dim], pivot->axis[dim]);
    #else
    return partition_dist_pivot_L2(ftr, pivot->p[dim], SMAP_PART_START(dim), SMAP_PART_DIM(dim));
    #endif
    #elif defined(SMAP_PARTITION_TYPE_CSQBP)
    return sub_dist_L2(ftr, pivot, dim);
    #endif
}

// 特徴データftrの射影像projectedを求める．
void psmap(smap_type projected, ftr_type ftr, smap_pivot_type *pivot) {
    for(int dim = 0; dim < SMAP_DIM; dim++) {
        projected[dim] = psmap_dim(dim, ftr, pivot);
    }
}

// 特徴データftrの量子化射影像のdim次元を求める．
smap_element_type q_psmap_dim(int dim, ftr_type ftr, smap_pivot_type *pivot, double offset, double slice) {
    #ifdef QUANTIZE_MIXED_MOD3
    int rtbl[] = {1 << QUANTIZE_BIT_0, 1 << QUANTIZE_BIT_1, 1 << QUANTIZE_BIT_2};
    int range = rtbl[dim % 3];
    #else
    int range = 1 << QUANTIZE_BIT; 
    #endif

    #if defined(SMAP_PARTITION_TYPE_PQBP) || defined(PARTITION_TYPE_PQBP) || defined(PARTITION_TYPE_QBP)
        #ifdef USE_PD
            dist_type d = partition_dist_pivot_L2(ftr, pivot->p[dim], pivot->num_axis[dim], pivot->axis[dim]);
        #else
            dist_type d = partition_dist_pivot_L2(ftr, pivot->p[dim], SMAP_PART_START(dim), SMAP_PART_DIM(dim));
        #endif
    #elif defined(SMAP_PARTITION_TYPE_CSQBP)
        dist_type d = sub_dist_L2(ftr, pivot, dim);
    #endif

    d = (d < offset ? 0 : d - offset) / slice;
    if(d >= range) d = range - 1;
    return (smap_element_type)d;
} 

// 特徴データftrの量子化射影像q_projectedを求める．
void q_psmap(smap_type q_projected, ftr_type ftr, smap_pivot_type *pivot, double offset[], double slice[]) {
    for(int dim = 0; dim < SMAP_DIM; dim++) {
        q_projected[dim] = q_psmap_dim(dim, ftr, pivot, offset[dim], slice[dim]);
    }
}

// 特徴データftrの量子化射影像q_projected（各次元は，unsigned char）を求める．
void q_uchar_psmap(unsigned char q_projected[], ftr_type ftr, smap_pivot_type *pivot, double offset[], double slice[]) {
    for(int dim = 0; dim < SMAP_DIM; dim++) {
        q_projected[dim] = q_psmap_dim(dim, ftr, pivot, offset[dim], slice[dim]);
    }
}

// psmap を量子化したunsigned char現に変換する．
void psmap2uchar_qpsmap(smap_type sm, unsigned char *uchar_smap, double offset[], double slice[], int range)
{
    // sm を量子化して各次元をucharで表現する．
    for(int j = 0; j < SMAP_DIM; j++) {
        smap_element_type d = sm[j];
        d = (d < offset[j] ? 0 : d - offset[j]) / slice[j];
        #ifndef QUANTIZE_MIXED
        if(d >= range) d = range - 1;
        #else
        if(d >= (j % 2 == 0 ? 4 : range)) d = (j % 2 == 0 ? 4 : range) - 1;
        #endif
        uchar_smap[j] = d;
    }
}
