//#pragma once
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
#include <string.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include "parm.h"
#include "bit_op.h"
#include "sketch.h"
#include "quick.h"
#include "e_time.h"

#define WRITE_BIT write_bit

#ifndef SMAP_DIM
#define SMAP_DIM 96
#endif

#define count_bits(n) (0x0000ffff&((0x00ff00ff&((0x0f0f0f0f&((0x33333333&((0x55555555&(n))\
	+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))\
	+((0xf0f0f0f0&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+\
	+((0xaaaaaaaa&(n))>>1)))>>2)))>>4)))+((0xff00ff00&((0x0f0f0f0f&((0x33333333&((0x55555555&(n))\
	+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))\
	+((0xf0f0f0f0&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))>>4)))>>8)))\
	+((0xffff0000&((0x00ff00ff&((0x0f0f0f0f&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))\
	+((0xf0f0f0f0&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))>>4)))\
	+((0xff00ff00&((0x0f0f0f0f&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))\
	+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))+((0xf0f0f0f0&((0x33333333&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))\
	+((0xcccccccc&((0x55555555&(n))+((0xaaaaaaaa&(n))>>1)))>>2)))>>4)))>>8)))>>16)

sketch_type data_to_sketch(ftr_type o, pivot_type *pivot)
{
	sketch_type sk = 0;
	for(int j = 0; j < PJT_DIM; j++) {
		#ifdef USE_PD_SKETCH
		write_bit(j, dist_pivot_L2_2(o, pivot->p[j], pivot->num_axis[j], pivot->axis[j], pivot->r[j]) < pivot->r[j], &sk);
		#elif defined(PARTITION_TYPE_PQBP)
		WRITE_BIT(j, PART_DISTANCE_PIVOT_2(o, pivot->p[j], PART_START(j), PART_DIM(j), pivot->r[j]) < pivot->r[j], &sk);
		#elif defined(PARTITION_TYPE_SQBP) || defined(PARTITION_TYPE_CSQBP)
		WRITE_BIT(j, sub_dist_L2_2(o, pivot, j, pivot->r[j]) < pivot->r[j], &sk);
		#else
		WRITE_BIT(j, DISTANCE_PIVOT_2(o, pivot->p[j], FTR_DIM, pivot->r[j]) < pivot->r[j], &sk);
		#endif
	}
	return sk;
}

void data_to_sketch_1bit(ftr_type o, pivot_type *pivot, int dim, sketch_type *sp)
{
	#ifdef USE_PD_SKETCH
	write_bit(dim, dist_pivot_L2_2(o, pivot->p[dim], pivot->num_axis[dim], pivot->axis[dim], pivot->r[dim]) < pivot->r[dim], sp);
	#elif defined(PARTITION_TYPE_PQBP)
	WRITE_BIT(dim, PART_DISTANCE_PIVOT_2(o, pivot->p[dim], PART_START(dim), PART_DIM(dim), pivot->r[dim])  < pivot->r[dim], sp);
	#elif defined(PARTITION_TYPE_SQBP) || defined(PARTITION_TYPE_CSQBP)
	WRITE_BIT(dim, sub_dist_L2_2(o, pivot, dim, pivot->r[dim])  < pivot->r[dim], sp);
	#else
	WRITE_BIT(dim, DISTANCE_PIVOT_2(o, pivot->p[dim], FTR_DIM, pivot->r[dim])  < pivot->r[dim], sp);
	#endif
}

// ピボット（中心点と半径の対 = pivot_type）の配列を動的に確保する．
// type は，基礎分割関数（GHP = 0, BP = 1, QBP = 2, PQBP = 3, SQBP = 4, CSQBP = 5）を指定する（現状ではQBPとPQBPのみ）．
// スケッチ幅は，PJT_DIMで指定したものを用いる．
pivot_type *new_pivot(int type)
{
	pivot_type *pivot = MALLOC(sizeof(pivot_type));
	for(int j = 0; j < PJT_DIM; j++) {
		pivot->p[j] = MALLOC(sizeof(ftr_element_type) * FTR_DIM);
		pivot->r[j] = 0;
	}
	pivot->type = type;
	return pivot;
}

void free_pivot(pivot_type *pivot)
{
	for(int j = 0; j < PJT_DIM; j++) {
		FREE(pivot->p[j], sizeof(ftr_element_type) * FTR_DIM);
	}
	FREE(pivot, sizeof(pivot_type));
}

// ピボットをファイルから読み込む．
// 作成済みのピボットを用いて検索などを行うときに使用する．
// filename = ピボットを格納したファイル名（csv形式）
// pivot = ピボットの配列（割り当ては，new_pivotなどを用いて別途行っておくこと）
void read_pivot(const char *filename, pivot_type *pivot)
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

	fprintf(stderr, "pivot-type = %d\n", pivot->type);
	if(pivot->type == CSQBP) {
		pivot->type = SQBP;
		fprintf(stderr, "CSQBP pivot is used as SQBP pivot.\n");
	}

	// 半径
	if(fgets(buf, MAX_LEN, pfp) != NULL) {
		pivot->r[0] = atoi(strtok(buf, ","));
		for(i = 1; i < PJT_DIM; i++)
			pivot->r[i] = atoi(strtok(NULL, ","));
	}
	
	// 中心点
	#ifndef PARTITION_TYPE_SQBP
		for(i = 0; fgets(buf, MAX_LEN, pfp) != NULL; i++) {
			pivot->p[i][0] = atoi(strtok(buf, ","));
			for(j = 1; j < FTR_DIM; j++) {
				pivot->p[i][j] = atoi(strtok(NULL, ","));
			}
			#ifdef USE_PD_SKETCH
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
	#else
		// SQBPではピボットの使用次元（used, num_used）を求めておく．
		// CSQBP で作成したピボットを用いるときは，PARTITION_TYPE_SQBP を定義しておくこと．
		for(i = 0; fgets(buf, MAX_LEN, pfp) != NULL; i++) {
			pivot->p[i][0] = atoi(strtok(buf, ","));
			pivot->num_used[i] = 0;
			if(pivot->p[i][0] != FTR_IGN) { pivot->used[i][pivot->num_used[i]++] = 0; }
			for(j = 1; j < FTR_DIM; j++) {
				pivot->p[i][j] = atoi(strtok(NULL, ","));
				if(pivot->p[i][j] != FTR_IGN) { pivot->used[i][pivot->num_used[i]++] = j; }
			}
		}
    #endif

	fclose(pfp); 
}

#ifdef PRE_ROTATION
void read_pivot_with_rotation(char *filename, pivot_type *pivot, int *rotation)
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

	// 分割数
	int num_p = atoi(strtok(NULL, ","));

	if(num_p <= 0) {
		fprintf(stderr, "invalid number of partition = %d\n", num_p);
	}
	// rotation table
	for(i = 0; i < FTR_DIM; i++) {
		rotation[i] = atoi(strtok(NULL, ","));
	}

	// 半径
	if(fgets(buf, MAX_LEN, pfp) != NULL) {
		pivot->r[0] = atoi(strtok(buf, ","));
		for(i = 1; i < PJT_DIM; i++)
			pivot->r[i] = atoi(strtok(NULL, ","));
	}
	
	// 中心点
    for(i = 0; fgets(buf, MAX_LEN, pfp) != NULL; i++){
		pivot->p[i][0] = atoi(strtok(buf, ","));
    	for(j = 1; j < FTR_DIM; j++) {
    		pivot->p[i][j] = atoi(strtok(NULL, ","));
    	}
    }
	fclose(pfp); 
}

void read_rotation(char *filename, int *rotation)
{
	FILE *pfp ;
    char buf[100000] = {0};

	pfp=fopen(filename, "r");
	if(pfp == NULL){
		fprintf(stderr, "cannot open pivot file = %s\n", filename);
		exit(0);
	}

	// 基礎分割関数
	int p_type;
	if(fgets(buf, MAX_LEN, pfp) != NULL)
		p_type = atoi(strtok(buf, ","));

	if(p_type != PQBP) {
		fprintf(stderr, "invalid pivot type = %d\n", p_type);
		exit(1);
	}
	// 分割数
	int num_p = atoi(strtok(NULL, ","));

	if(num_p <= 0) {
		fprintf(stderr, "invalid number of partition = %d\n", num_p);
	}
	// rotation table
	for(int i = 0; i < FTR_DIM; i++) {
		rotation[i] = atoi(strtok(NULL, ","));
	}

	fclose(pfp); 
}

void write_pivot_with_rotation(char *filename, pivot_type *pivot, int *rotation)
{
	int i, j;
	if(pivot->type != PQBP) {
		fprintf(stderr, "Pivot with rotation is valid only for PQBP!\n");
		exit(1);
	}
	FILE *fp = fopen(filename, "w");
	if(fp == NULL) {
		printf("cannot write open file = %s\n", filename);
	} else {
		fprintf(fp, "%d, %d", pivot->type, NUM_PART);		// 基礎分割関数, 分割数
		for(i = 0; i < FTR_DIM; i++) {
			fprintf(fp, ", %d", rotation[i]);				// rotation 表
		}
		fprintf(fp, "\n");
		for(i = 0; i < PJT_DIM - 1; i++)	 	// 半径
			fprintf(fp, "%d,", pivot->r[i]);
		fprintf(fp, "%d\n", pivot->r[i]);
		for(i = 0; i < PJT_DIM; i++) {			// 中心点
			for(j = 0; j < FTR_DIM - 1; j++)
				fprintf(fp, "%d,", pivot->p[i][j]);
			fprintf(fp, "%d\n", pivot->p[i][j]);
		}
	}
	fclose(fp);
}
#endif

// ピボットをファイルに書き出す．
// QBPで作成したり，AIRで最適化して作成したピボットをcsv形式で書き出す．
// 1行目は，QBPやPQBPのtype
// 2行目は，PJT_DIM個の半径をコンマ区切りで
// 3行目以降は，ピボットの中心点（FTR_DIM次元）
// pivot_property.sh は，この形式から分割type(QBP or PQBP)やPJT_DIM，FTR_DIMを調べるスクリプト．
void write_pivot(char *filename, pivot_type *pivot)
{
	int i, j;
	FILE *fp = fopen(filename, "w");
	if(fp == NULL) {
		printf("cannot write open file = %s\n", filename);
	} else {
		if(pivot->type == PQBP) {
			fprintf(fp, "%d, %d\n", pivot->type, NUM_PART);		// 基礎分割関数, 分割数
		} else {
			fprintf(fp, "%d, 1\n", pivot->type);				// 基礎分割関数，分割数 = 1
		}
		for(i = 0; i < PJT_DIM - 1; i++)	 	// 半径
			fprintf(fp, "%d,", pivot->r[i]);
		fprintf(fp, "%d\n", pivot->r[i]);
		for(i = 0; i < PJT_DIM; i++) {			// 中心点
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

// 全データのスケッチ配列からバケット表（idx と bkt）を含む構造体（struct_bucket）を作成する．
// ただし，PJT_DIM（スケッチ幅）< 32 でのみ使用すること．
// PJT_DIMが大きいと，バケット表のための配列bktが非常に大きくなるので，使用できない．
struct_bucket *new_bucket(int num_data, sketch_type sk[])
// version 2.0 (using arrays only idx and bucket)
{
	struct_bucket *b = (struct_bucket *)malloc(sizeof(struct_bucket));
	b->num_data = num_data;
	b->ftr_data = NULL;
	b->sk = NULL;

	b->idx = (int *)malloc(sizeof(int) * num_data);
	b->bkt = (int *)calloc((1L << PJT_DIM) + 2, sizeof(int));

	int *bucket = b->bkt + 1;
	for(int i = 0; i < num_data; i++) bucket[sk[i] + 1]++;
	for(long i = 0; i < (1L << PJT_DIM); i++) bucket[i + 1] += bucket[i];
	for(int i = 0; i < num_data; i++) b->idx[(long)bucket[sk[i]]++] = i;

	b->num_nonempty_buckets = 0;
	for(long s = 0; s < (1L << PJT_DIM); s++) {
		b->num_nonempty_buckets += (b->bkt[s + 1] - b->bkt[s] > 0);
	}
	b->sk_num = (sk_num_pair *)malloc(sizeof(sk_num_pair) * b->num_nonempty_buckets);
	int num_elements, j = 0;
	for(long s = 0; s < (1L << PJT_DIM); s++) {
		num_elements = b->bkt[s + 1] - b->bkt[s];
		if(num_elements > 0) {
			b->sk_num[j++] = (sk_num_pair) {(sketch_type)s, num_elements, b->bkt[s]};
		}
	}

	return b;
}

int comp_sketch(sketch_type a, sketch_type b)
{
	if(a < b)
		return -1;
	else if(a == b)
		return 0;
	else
		return 1;
}

int comp_uint(const void *a, const void *b) {
	if(*((unsigned int *) a) < *((unsigned int *) b))
		return -1;
	else if(*((unsigned int *) a) == *((unsigned int *) b))
		return 0;
	else
		return 1;
}

// バケット表（idx と bkt）を含む構造体（struct_bucket）をコンパクトな形式でファイルに書き出す
// bkt は PJT_DIM (width) に対して指数オーダーなので，書き出さない．空でないバケット情報のみ書き出す．
void write_bucket(char *filename, struct_bucket *b)
{
	int num_data = b->num_data, *idx = b->idx, num_nonempty_buckets = b->num_nonempty_buckets;
	sk_num_pair *sk_num = b->sk_num;

	FILE *fp;

	if((fp = fopen(filename, "wb"))  == NULL) {
		fprintf(stderr, "Write open bucket file error, file name = %s\n", filename);
		exit(0);
	}
	if(fwrite(&num_data, sizeof(int), 1, fp) != 1) {  // num_data を書き出す
		fprintf(stderr, "fwrite error (num_data) file = %s\n", filename);
		exit(0);
	}
	if(fwrite(idx, sizeof(int) * num_data, 1, fp) != 1) {  // idx[data_num] を書き出す
		fprintf(stderr, "fwrite error (idx) file = %s\n", filename);
		exit(0);
	}

	if(fwrite(&num_nonempty_buckets, sizeof(int), 1, fp) != 1) {  // num_nonempty_buckets を書き出す
		fprintf(stderr, "fwrite error (num_nonempty_buckets) file = %s\n", filename);
		exit(0);
	}
	fprintf(stderr, "num_nonempty_buckets = %d, average number of elements in nonempty buckets = %.2lf\n", num_nonempty_buckets, (double)num_data / num_nonempty_buckets);
	for(int i = 0; i < num_nonempty_buckets; i++) {
		if(fwrite(&sk_num[i].sk, sizeof(sketch_type), 1, fp) != 1) {  // sketch s を書き出す
			fprintf(stderr, "fwrite error (sketch) file = %s\n", filename);
			exit(0);
		}
		if(fwrite(&sk_num[i].num, sizeof(int), 1, fp) != 1) {  // num_elements を書き出す
			fprintf(stderr, "fwrite error (num_elements) file = %s\n", filename);
			exit(0);
		}
	}
	fclose(fp);
	return;
}

void free_bucket(struct_bucket *b)
{
	int num_data = b->num_data;
	FREE(b->sk, sizeof(sketch_type) * num_data);
	#ifdef REVERSE_IDX
	FREE(b->r_idx, sizeof(int) * num_data);
	#endif
	free(b->sk_num);
	free(b->bkt);
	FREE(b->idx, sizeof(int) * num_data);
	free(b);
}

// コンパクトな形式でファイルに保存されたバケット表を読み込んで（idx と bkt）を含む構造体（struct_bucket）に展開する
struct_bucket *read_bucket(char *filename)
{
	FILE *fp;
	if((fp = fopen(filename, "rb"))  == NULL) {
		fprintf(stderr, "Read open bucket file error, file name = %s\n", filename);
		exit(0);
	}

	int num_data;
	if(fread(&num_data, sizeof(int), 1, fp) != 1) {  // ファイルに書かれている num_data を読み込む
		fprintf(stderr, "fread error (num_data) file = %s\n", filename);
		exit(0);
	}

	fprintf(stderr, "num_data in file %s = %d\n", filename, num_data);
	use_system("VmSize");

	struct_bucket *b = (struct_bucket *)malloc(sizeof(struct_bucket));
	b->num_data = num_data;
	b->ftr_data = NULL;  // ここではデータセットは読み込まない．オンメモリ検索をするときには，必要に応じて ftr または sftr を読み込む
	b->idx = MALLOC(sizeof(int) * num_data);
	if(fread(b->idx, sizeof(int) * num_data, 1, fp) != 1) {  // idx[num_data] を読み込む
		fprintf(stderr, "fread error (idx, size = %ld) file = %s\n", sizeof(int) * num_data, filename);
		exit(0);
	}

	if(fread(&b->num_nonempty_buckets, sizeof(int), 1, fp) != 1) {  // ファイルに書かれている num_nonempty_buckets を読み込む
		fprintf(stderr, "fread error (data_num) file = %s\n", filename);
		exit(0);
	}
	use_system("VmSize");

	b->bkt = (int *)calloc((1L << PJT_DIM) + 2, sizeof(int));
	if(b->bkt != NULL) {
	} else {
		fprintf(stderr, "calloc bkt failed. exit!\n");
	}
	b->sk_num = (sk_num_pair *)malloc(sizeof(sk_num_pair) * b->num_nonempty_buckets);
	int num, j;
	sk_num_pair snp;
	long int s_next = 0;
	num = j = 0;
	for(int i = 0; i < b->num_nonempty_buckets; i++) {
		if(fread(&snp.sk, sizeof(sketch_type), 1, fp) != 1) {  // sketch を読み込む
			fprintf(stderr, "fread error (sketch) file = %s\n", filename);
			exit(0);
		}
		if(fread(&snp.num, sizeof(int), 1, fp) != 1) {  // num_elements を読み込む
			fprintf(stderr, "fread error (num_elements) file = %s\n", filename);
			exit(0);
		}
		snp.pos = num;
		b->sk_num[j++] = snp; // (sk_num_pair) {s, num_elements, num};
		while(s_next <= snp.sk) {
			b->bkt[s_next++] = num;
		}
		num += snp.num;
	}
	while(s_next < (1L << PJT_DIM) + 2) {
		b->bkt[s_next++] = num;
	}

	#ifdef REVERSE_IDX
	b->r_idx = MALLOC(sizeof(int) * num_data);
	for(int i = 0; i < num_data; i++) {
		b->r_idx[b->idx[i]] = i;
	}
	fprintf(stderr, "make reverse idx OK\n");
	#endif

	fprintf(stderr, "num_nonempty_buckets = %d, average number of elements in nonempty buckets = %.2lf\n", b->num_nonempty_buckets, (double)num_data / b->num_nonempty_buckets);
	fclose(fp);

	b->sk = MALLOC(sizeof(sketch_type) * num_data);
	for(int i = 0; i < b->num_nonempty_buckets; i++) {
		for(int j = b->sk_num[i].pos; j < b->sk_num[i].pos + b->sk_num[i].num; j++) {
			memcpy(&b->sk[b->idx[j]], &b->sk_num[i].sk, sizeof(sketch_type));
		}
	}
	use_system("VmSize");

	#ifdef WITH_STAT
	for(int j = 0; j < PJT_DIM; j++) b->stat[j] = 0;
	for(sketch_type s = 0; s < (1 << PJT_DIM); s++) {
		if(b->bkt[s] == b->bkt[s + 1]) continue;
		sketch_type t = s;
		for(int j = 0; j < PJT_DIM; j++) {
			b->stat[j] += (t & 1) * (b->bkt[s + 1] - b->bkt[s]);
			t >>= 1;
		}
	}

	printf("stat: data_num = %d\n", num_data);
	for(int j = 0; j < PJT_DIM; j++) {
		printf("%12d", b->stat[j]);
	}
	getchar();
	#endif

	return b;
}

// bkt ファイルに保存されているバケットを読み込む．ただし，bkt は，配列に展開しない．
// コンパクト表現 compact_bucket，つまり sk_num_pair（空でないスケッチとデータ数の対＋オフセット）の配列を用いる．
struct_bucket *read_compact_bucket(char *filename)
{
	FILE *fp;
	if((fp = fopen(filename, "rb"))  == NULL) {
		fprintf(stderr, "Read open bucket file error, file name = %s\n", filename);
		exit(0);
	}

	int num_data;
	if(fread(&num_data, sizeof(int), 1, fp) != 1) {  // ファイルに書かれている num_data を読み込む
		fprintf(stderr, "fread error (num_data) file = %s\n", filename);
		exit(0);
	}

	struct_bucket *b = (struct_bucket *)malloc(sizeof(struct_bucket));
	b->num_data = num_data;
	b->ftr_data = NULL;  // ここではデータセットは読み込まない．オンメモリ検索をするときには，必要に応じて ftr または sftr を読み込む
	b->bkt = NULL; // bkt は配列に展開しない
	b->idx = (int *)malloc(sizeof(int) * num_data);
	if(fread(b->idx, sizeof(int) * num_data, 1, fp) != 1) {  // idx[num_data] を読み込む
		fprintf(stderr, "fread error (idx, size = %ld) file = %s\n", sizeof(int) * num_data, filename);
		exit(0);
	}

	if(fread(&b->num_nonempty_buckets, sizeof(int), 1, fp) != 1) {  // ファイルに書かれている num_nonempty_buckets を読み込む
		fprintf(stderr, "fread error (num_nonempty_buckets) file = %s\n", filename);
		exit(0);
	}

	b->sk_num = NULL; // おそらく不要なので，動作確認したら，このメンバを削除できるかも
	b->sk = (sketch_type *)malloc(sizeof(sketch_type) * num_data);
	int offset = 0;
	sk_num_pair snp;
	for(int i = 0; i < b->num_nonempty_buckets; i++) {
		if(fread(&snp.sk, sizeof(sketch_type), 1, fp) != 1) {  // sketch を読み込む
			fprintf(stderr, "fread error (sketch) file = %s\n", filename);
			exit(0);
		}
		if(fread(&snp.num, sizeof(int), 1, fp) != 1) {  // num_elements を読み込む
			fprintf(stderr, "fread error (num_elements) file = %s\n", filename);
			exit(0);
		}
		snp.pos = offset;
		for(int j = snp.pos; j < snp.pos + snp.num; j++) {
			b->sk[b->idx[j]] = snp.sk;
		}
		offset += snp.num;
	}
	fprintf(stderr, "num_nonempty_buckets = %d, average number of elements in nonempty buckets = %.2lf\n", b->num_nonempty_buckets, (double)num_data / b->num_nonempty_buckets);
	fclose(fp);
	return b;
}

// bkt ファイルに保存されているsk_numを読み込むためにopenする．
struct_bucket_sk_num *open_bucket_sk_num(char *filename)
{
	struct_bucket_sk_num *bsk = (struct_bucket_sk_num *)malloc(sizeof(struct_bucket_sk_num));
	bsk->filename = filename;
	if((bsk->fp = fopen(filename, "rb"))  == NULL) {
		fprintf(stderr, "Read open bucket file error, file name = %s\n", filename);
		exit(0);
	}
	int num_data;
	if(fread(&num_data, sizeof(int), 1, bsk->fp) != 1) {  // ファイルに書かれている num_data を読み込む
		fprintf(stderr, "fread error (num_data) file = %s\n", filename);
		exit(0);
	}
	bsk->num_data = num_data;
	if(fseek(bsk->fp, (long)(sizeof(int) * num_data), SEEK_CUR) != 0) {  // idx[num_data] を読み飛ばす
		fprintf(stderr, "fseek error (to skip ibk, size = %ld, file = %s)\n", sizeof(int) * num_data, filename);
		exit(0);
	}
	if(fread(&bsk->num_nonempty_buckets, sizeof(int), 1, bsk->fp) != 1) {  // ファイルに書かれている num_nonempty_buckets を読み込む
		fprintf(stderr, "fread error (num_nonempty_buckets) file = %s\n", filename);
		exit(0);
	}
	if(fread(&bsk->sk_num.sk, sizeof(sketch_type), 1, bsk->fp) != 1) {  // sketch を読み込む
		fprintf(stderr, "fread error (sketch) file = %s\n", filename);
		exit(0);
	}
	if(fread(&bsk->sk_num.num, sizeof(int), 1, bsk->fp) != 1) {  // num_elements を読み込む
		fprintf(stderr, "fread error (num_elements) file = %s\n", filename);
		exit(0);
	}
	fprintf(stderr, "num_nonempty_buckets = %d, average number of elements in nonempty buckets = %.2lf\n", bsk->num_nonempty_buckets, (double)num_data / bsk->num_nonempty_buckets);
	bsk->processed_buckets = 0;

	return bsk;
}

int read_next_bucket_sk_num(struct_bucket_sk_num *bsk)
{
	if(++bsk->processed_buckets >= bsk->num_nonempty_buckets) {
		return 0;	// EOF (all sk_num_pairs are processed)
	}
	if(fread(&bsk->sk_num.sk, sizeof(sketch_type), 1, bsk->fp) != 1) {  // sketch を読み込む
		fprintf(stderr, "fread error (sketch) file = %s, num_nonempty_buckets = %d, processed_buckets = %d\n", bsk->filename, bsk->num_nonempty_buckets, bsk->processed_buckets);
		exit(0);
	}
	if(fread(&bsk->sk_num.num, sizeof(int), 1, bsk->fp) != 1) {  // num_elements を読み込む
		fprintf(stderr, "fread error (num_elements) file = %s, num_nonempty_buckets = %d, processed_buckets = %d\n", bsk->filename, bsk->num_nonempty_buckets, bsk->processed_buckets);
		exit(0);
	}

	return 1;
}

#define PARENT(i) ((i)>>1)
#define LEFT(i)   ((i)<<1)
#define RIGHT(i)  (((i)<<1)+1)

void min_heapify_p(int i, struct_que *que)
{
    int l, r;
    int smallest;

    l = LEFT(i);
	r = RIGHT(i);
    if (l < que->qsize && que->element[l].key < que->element[i].key) smallest = l; else smallest = i;
    if (r < que->qsize && que->element[r].key < que->element[smallest].key) smallest = r;
    if (smallest != i) {
        QUE t = que->element[i]; que->element[i] = que->element[smallest]; que->element[smallest] = t;
        min_heapify_p(smallest, que);
    }
}

int deq_p(QUE *q, struct_que *que)
{
    if (que->qsize == 0) return 0;
    memcpy(q, &(que->element[0]), sizeof(QUE));
	que->element[0] = que->element[--(que->qsize)];
    min_heapify_p(0, que);
    return 1;
}

void enq_p(QUE *q, struct_que *que)
{
    int i, ii;

	i = (que->qsize)++;
    memcpy(&(que->element[i]), q, sizeof(QUE));
	while (i > 0 && que->element[ii = PARENT(i)].key > que->element[i].key) {
        QUE t = que->element[i];
		que->element[i] = que->element[ii]; 
		que->element[ii] = t;
        i = ii;
    }
}

void make_empty_que_c2_n(struct_que_c2_n *que)
{
	que->qsize = 0;
	que->detail_size = 0;
}

int new_que_e2_n(struct_que_c2_n *que)
{
	int i = que->detail_size;
	que->detail_size++;
	return i;
}

void min_heapify_c2_n(int i, struct_que_c2_n *que)
{
	QUE_c2 *element = que->element;
    QUE_c2 t = element[i];

    while(1) {
        int l = LEFT(i);
        if (que->qsize <= l) {
            break;
        }
        if (element[l].key < t.key) {
            int r = RIGHT(i);
            if (r < que->qsize && element[r].key < element[l].key) {
                element[i] = element[r];
                i = r;
            }
            else {
                element[i] = element[l];
                i = l;
            }
        }
        else {
            int r = RIGHT(i);
            if (!(r < que->qsize && element[r].key < t.key)) {
                break;
            }
            element[i] = element[r];
            i = r;
        }
    }
    element[i] = t;
}

int deq_c2_n(QUE_c2 *qe, struct_que_c2_n *que)
{
    if (que->qsize == 0) return 0;
    *qe = que->element[0];
    que->element[0] = que->element[que->qsize - 1];
	--(que->qsize);
    min_heapify_c2_n(0, que);
    return 1;
}

void deq_c2_n_del(struct_que_c2_n *que) // Deq で取り除かなかったものを取り除く 
{
    if (que->qsize == 0) return;
    que->element[0] = que->element[que->qsize - 1];
	--(que->qsize);
    min_heapify_c2_n(0, que);
}

void enq_c2_n(QUE_c2 *qe, struct_que_c2_n *que)
{
	QUE_c2 *element = que->element;
    int i, ii;
	QUE_c2 t;

    i = (que->qsize)++;
    element[i] = *qe;
    while (i > 0 && element[ii = PARENT(i)].key > element[i].key) {
        t = element[i];
        element[i] = element[ii]; 
        element[ii] = t;
        i = ii;
    }
}

void enq_c2_n_after_deq(QUE_c2 *qe, struct_que_c2_n *que)
{
    que->element[0] = *qe;
    min_heapify_c2_n(0, que);
}

// 質問 query (質問番号，ftr) のための構造体（スケッチなどを準備する
// 正解（最近傍）に関するメンバ(answerとanswer_sketch)の設定は行わない
void set_query_sketch(struct_query_sketch *qs, query_type *query, pivot_type *pivot)
{
	qs->query = *query;
	qs->sketch = 0;
	int tbl_size = 4;
	for(int p = 0; p < tbl_size; p++) {
		for(int n = 0; n < 256; n++) {
			qs->tbl[p][n] = 0;
		}
	}
	#if PRIORITY == 0 // hamming
		for(int j = 0; j < PJT_DIM; j++) {
			#ifndef PARTITION_TYPE_PQBP
			dist_type dist = DISTANCE(query->ftr, pivot->p[j], FTR_DIM);
			#else
			dist_type dist = PART_DISTANCE(query->ftr, pivot->p[j], PART_START(j), PART_DIM(j));
			#endif
			WRITE_BIT(j, dist < pivot->r[j], &qs->sketch);
		}
	#else // score_1, score_2, score_inf
		for(int j = 0; j < PJT_DIM; j++) qs->idx[j] = j;
		for(int j = 0; j < PJT_DIM; j++) {
			#ifndef PARTITION_TYPE_PQBP
			dist_type dist = DISTANCE(query->ftr, pivot->p[j], FTR_DIM);
			#else
			dist_type dist = PART_DISTANCE(query->ftr, pivot->p[j], PART_START(j), PART_DIM(j));
			#endif
			WRITE_BIT(j, dist < pivot->r[j], &qs->sketch);
			#ifndef SQRT_FTR
				#ifdef SCORE_2
				qs->bd[j] = fabs(sqrt(dist) - sqrt(pivot->r[j])) * fabs(sqrt(dist) - sqrt(pivot->r[j])); // score_2
				#else
				qs->bd[j] = fabs(sqrt(dist) - sqrt(pivot->r[j])); // score_1, score_inf
				#endif
			#else
				#ifdef SCORE_2
				fprintf(stderr, "(1) set_query_sketch\n"); exit(0);
				qs->bd[j] = fabs((int)dist - (int)pivot->r[j]) * fabs((int)dist - (int)pivot->r[j]); // score_2
				#else
				qs->bd[j] = fabs((int)dist - (int)pivot->r[j]); // score_1, score_inf
				#endif
			#endif

			// bd をソート（idxを用いた相対ソート）（挿入法：bd[j] を挿入）
			int l;
			for(l = j - 1; l >= 0 && qs->bd[qs->idx[l]] > qs->bd[j]; l--) {
				qs->idx[l + 1] = qs->idx[l];
			}
			qs->idx[l + 1] = j;

			// 表関数の更新（bd[j] を対応するところに加える）
			int tn = j / 8, bp = j % 8; // 表関数の変更すべき表番号が tn，そのビット位置が bp
			for(int n = 0; n < 256; n++) {
				if(n & (1 << bp)) {
					qs->tbl[tn][n] += qs->bd[j];
				}
			}
		}
	#endif
}

// 各質問のquery_sketchのsketchを求めて，質問点と分割境界との最小距離に反対側にある点と分割境界の最小距離の平均距離＋2σを求める
void compute_sketch_and_boundary_plus(dist_type bd_plus[][PJT_DIM], int num_queries, struct_query_sketch *qs_all, query_type *query_all, pivot_type *pivot)
{
	for(int q = 0; q < num_queries; q++) {
		struct_query_sketch *qs = &qs_all[q];
		query_type *query = &query_all[q];
		qs->query = *query;
		// 質問のsketchをゼロクリアする
		qs->sketch = 0;
		for(int j = 0; j < PJT_DIM; j++) {
			// dist = j-bit のピボット中心と質問[q]の距離（PQBPでは部分距離）
			#if defined(PARTITION_TYPE_PQBP)
			dist_type dist = PART_DISTANCE(query->ftr, pivot->p[j], PART_START(j), PART_DIM(j));
			#elif defined(PARTITION_TYPE_SQBP) || defined(PARTITION_TYPE_CSQBP)
			dist_type dist = sub_dist_L2(query->ftr, pivot, j); // 未対応
			#else
			dist_type dist = DISTANCE(query->ftr, pivot->p[j], FTR_DIM); // 未対応
			#endif
			// 質問のsketchの j-bit 目をセットする．
			WRITE_BIT(j, dist < pivot->r[j], &qs->sketch);
			#ifndef SQRT_FTR
			double bd = fabs(sqrt(dist) - sqrt(pivot->r[j])); // query と j-bit の分割境界との最小距離
			#else
			double bd = fabs(dist - pivot->r[j]); // query と j-bit の分割境界との最小距離
			#endif
			// query と j-bit の分割境界と反対側にある質問の分割境界との最小距離の平均を求める．
			double sum = 0, sum2 = 0; int count = 0; dist_type min_dist = UINT_MAX;
			for(int x = 0; x < num_queries; x++) {
				// xdist = j-bit のピボット中心と質問[x]の距離（PQBPでは部分距離）
				#if defined(PARTITION_TYPE_PQBP)
				dist_type xdist = PART_DISTANCE(query_all[x].ftr, pivot->p[j], PART_START(j), PART_DIM(j));
				#elif defined(PARTITION_TYPE_SQBP) || defined(PARTITION_TYPE_CSQBP)
				dist_type xdist = sub_dist_L2(query_all[x].ftr, pivot, j); // 未対応
				#else
				dist_type xdist = DISTANCE(query_all[q].ftr, pivot->p[j], FTR_DIM); // 未対応
				#endif
				if((dist < pivot->r[j]) != (xdist < pivot->r[j])) { // query と query_all[x] は，分割境界の反対側
					#ifndef SQRT_FRT
						sum += fabs(sqrt(xdist) - sqrt(pivot->r[j])); sum2 += (sqrt(xdist) - sqrt(pivot->r[j])) * (sqrt(xdist) - sqrt(pivot->r[j]));
						count++;
						if(fabs(sqrt(xdist) - sqrt(pivot->r[j])) < min_dist) min_dist = fabs(sqrt(xdist) - sqrt(pivot->r[j]));
					#else
						sum += fabs(xdist - pivot->r[j]); sum2 += (xdist - pivot->r[j]) * (xdist - pivot->r[j]);
						count++;
						if(fabs(xdist - pivot->r[j]) < min_dist) min_dist = fabs(xdist - pivot->r[j]);
					#endif
				}
			}
			double ave = sum / count;
			double stdev = sqrt(sum2 / count - ave * ave);
			bd_plus[q][j] = sqrt(bd * bd + (ave - 1.2 * stdev) * (ave - 1.2 * stdev));
		}
	}
}

void set_query_sketch_p_boundary_plus(dist_type bd_plus[PJT_DIM], struct_query_sketch *qs, query_type *query, pivot_type *pivot, double p)
{
	qs->query = *query;
	int tbl_size = 4;
	for(int p = 0; p < tbl_size; p++) {
		for(int n = 0; n < 256; n++) {
			qs->tbl[p][n] = 0;
		}
	}

	for(int j = 0; j < PJT_DIM; j++) qs->idx[j] = j;
	for(int j = 0; j < PJT_DIM; j++) {
		qs->bd[j] = pow(bd_plus[j], p) ; // score_p

		// bd をソート（idxを用いた相対ソート）（挿入法：bd[j] を挿入）
		int l;
		for(l = j - 1; l >= 0 && qs->bd[qs->idx[l]] > qs->bd[j]; l--) {
			qs->idx[l + 1] = qs->idx[l];
		}
		qs->idx[l + 1] = j;

		// 表関数の更新（bd[j] を対応するところに加える）
		int tn = j / 8, bp = j % 8; // 表関数の変更すべき表番号が tn，そのビット位置が bp
		for(int n = 0; n < 256; n++) {
			if(n & (1 << bp)) {
				qs->tbl[tn][n] += qs->bd[j];
			}
		}
	}
}

#ifdef PLUS_SD
// 質問query[q]のスケッチquery_sketch[q]に求めて，質問をサンプルとして，各射影次元の射影像（質問とピボットの中心点との距離）の平均と標準偏差を求める．
void compute_query_sketch_and_ave_stdev_bd_of_pjt_dim(double ave[], double stdev[], int num_queries, struct_query_sketch query_sketch[], query_type query[], pivot_type *pivot)
{
	double sum[PJT_DIM] = {0}, sum2[PJT_DIM] = {0};
	for(int q = 0; q < num_queries; q++) {
		query_sketch[q].query = query[q];
		// 質問のsketchをゼロクリアする
		query_sketch[q].sketch = 0;
		for(int j = 0; j < PJT_DIM; j++) {
			// dist = j-bit のピボット中心と質問[q]の距離（PQBPでは部分距離）
			#if defined(PARTITION_TYPE_PQBP)
			dist_type dist = PART_DISTANCE(query[q].ftr, pivot->p[j], PART_START(j), PART_DIM(j));
			#elif defined(PARTITION_TYPE_SQBP) || defined(PARTITION_TYPE_CSQBP)
			dist_type dist = sub_dist_L2(query[q].ftr, pivot, j); // 未対応
			#else
			dist_type dist = DISTANCE(query[q].ftr, pivot->p[j], FTR_DIM); // 未対応
			#endif
			// 質問のsketchの j-bit 目をセットする．
			WRITE_BIT(j, dist < pivot->r[j], &query_sketch[q].sketch);
			#ifndef SQRT_FTR
				sum[j] += sqrt(dist);
				sum2[j] += dist;
				query_sketch[q].bd[j] = fabs(sqrt(dist) - sqrt(pivot->r[j])); // query と j-bit の分割境界との最小距離
			#else
				sum[j] += dist;
				sum2[j] += dist * dist;
				query_sketch[q].bd[j] = fabs((int)dist - (int)pivot->r[j]); // query と j-bit の分割境界との最小距離
			#endif
		}
	}
	for(int j = 0; j < PJT_DIM; j++) {
		ave[j] = sum[j] / num_queries;
		stdev[j] = sqrt(sum2[j] / num_queries - ave[j] * ave[j]);
	}
}

// 質問query[q]のスケッチquery_sketch[q]に求めて，質問をサンプルとして，各射影次元のスケッチ0および1となる射影像（質問とピボットの中心点との距離）の平均を求める．
void compute_query_sketch_and_ave0_ave1_of_pjt_dim(double ave0[], double ave1[], int num_queries, struct_query_sketch query_sketch[], query_type query[], pivot_type *pivot)
{
	double sum0[PJT_DIM] = {0}, sum1[PJT_DIM] = {0};
	int num0[PJT_DIM] = {0}, num1[PJT_DIM] = {0};
	for(int q = 0; q < num_queries; q++) {
		query_sketch[q].query = query[q];
		// 質問のsketchをゼロクリアする
		query_sketch[q].sketch = 0;
		for(int j = 0; j < PJT_DIM; j++) {
			// dist = j-bit のピボット中心と質問[q]の距離（PQBPでは部分距離）
			#if defined(PARTITION_TYPE_PQBP)
			dist_type dist = PART_DISTANCE(query[q].ftr, pivot->p[j], PART_START(j), PART_DIM(j));
			#elif defined(PARTITION_TYPE_SQBP) || defined(PARTITION_TYPE_CSQBP)
			dist_type dist = sub_dist_L2(query[q].ftr, pivot, j); // 未対応
			#else
			dist_type dist = DISTANCE(query[q].ftr, pivot->p[j], FTR_DIM); // 未対応
			#endif
			// 質問のsketchの j-bit 目をセットする．
			WRITE_BIT(j, dist < pivot->r[j], &query_sketch[q].sketch);
			#ifndef SQRT_FTR
				if(dist < pivot->r[j]) {
					num1[j]++;
					sum1[j] += sqrt(dist);
				} else {
					num0[j]++;
					sum0[j] += sqrt(dist);
				}
				query_sketch[q].bd[j] = sqrt(dist); // query のピボット中心との距離，つまりsmap射影像
			#else
				if(dist < pivot->r[j]) {
					num1[j]++;
					sum1[j] += dist;
				} else {
					num0[j]++;
					sum0[j] += dist;
				}
				query_sketch[q].bd[j] = dist; // query のピボット中心との距離，つまりsmap射影像
			#endif
		}
	}
	for(int j = 0; j < PJT_DIM; j++) {
		double ave = (sum0[j] + sum1[j]) / num_queries;
		ave0[j] = sum0[j] / num0[j];
		ave1[j] = sum1[j] / num1[j];
		printf("ave1[%d] = %8.2lf, ave[%d] = %8.2lf, ave0[%d] = %8.2lf\n", j, ave1[j], j, ave, j, ave0[j]);
	}
}

// query_sketch の 表関数を作成する（plus stdev）
void set_query_sketch_p_plus_sd(double ave[], double stdev[], struct_query_sketch *query_sketch, double p)
{
	int tbl_size = 4;

	for(int p = 0; p < tbl_size; p++) {
		for(int n = 0; n < 256; n++) {
			query_sketch->tbl[p][n] = 0;
		}
	}

	for(int j = 0; j < PJT_DIM; j++) query_sketch->idx[j] = j;
	for(int j = 0; j < PJT_DIM; j++) {
		if(PLUS_SD < 0) {
			if(p != DBL_MAX) {
				query_sketch->bd_plus[j] = pow(query_sketch->bd[j], p); // score_p
			} else {
				query_sketch->bd_plus[j] = query_sketch->bd[j]; // score_inf
			}
		} else {
			if(p != DBL_MAX) {
				query_sketch->bd_plus[j] = pow(query_sketch->bd[j] + PLUS_SD * stdev[j], p) - pow(fabs(query_sketch->bd[j] - PLUS_SD * stdev[j]), p); // score_p
			} else {
				query_sketch->bd_plus[j] = query_sketch->bd[j] + PLUS_SD * stdev[j] - fabs(query_sketch->bd[j] - PLUS_SD * stdev[j]); // score_inf
			}
		}
		// bd_plus をソート（idxを用いた相対ソート）（挿入法：bd_plus[j] を挿入）
		int l;
		for(l = j - 1; l >= 0 && query_sketch->bd_plus[query_sketch->idx[l]] > query_sketch->bd_plus[j]; l--) {
			query_sketch->idx[l + 1] = query_sketch->idx[l];
		}
		query_sketch->idx[l + 1] = j;
		// 表関数の更新（bd[j] を対応するところに加える）
		int tn = j / 8, bp = j % 8; // 表関数の変更すべき表番号が tn，そのビット位置が bp
		for(int n = 0; n < 256; n++) {
			if(n & (1 << bp)) {
				if(p != DBL_MAX) {
					query_sketch->tbl[tn][n] += query_sketch->bd_plus[j];
				} else {
					if(query_sketch->tbl[tn][n] < query_sketch->bd_plus[j]) {
						query_sketch->tbl[tn][n] = query_sketch->bd_plus[j];
					}
				}
			}
		}
	}
}

// query_sketch の 表関数を作成する（using ave0 and ave1）
void set_query_sketch_p_ave0_ave1(double ave0[], double ave1[], struct_query_sketch *query_sketch, double p)
{
	int tbl_size = 4;

	for(int p = 0; p < tbl_size; p++) {
		for(int n = 0; n < 256; n++) {
			query_sketch->tbl[p][n] = 0;
		}
	}

	for(int j = 0; j < PJT_DIM; j++) query_sketch->idx[j] = j;
	for(int j = 0; j < PJT_DIM; j++) {
		double dave0 = PLUS_SD * ave0[j] + (1 - PLUS_SD) * ave1[j];
		double dave1 = PLUS_SD * ave1[j] + (1 - PLUS_SD) * ave0[j];
		if(dave0 < dave1) {
			double ta = dave0;
			dave0 = dave1;
			dave1 = ta;
		}
		double dave2 = 0.5 * ave1[j] + 0.5 * ave0[j]; // 半径の代わり
		if(query_sketch->bd[j] <  dave1) {
			query_sketch->bd_plus[j] = pow(dave2 - query_sketch->bd[j], p);
		} else if(query_sketch->bd[j] >  dave0) {
			query_sketch->bd_plus[j] = pow(query_sketch->bd[j] - dave2, p);
		} else {
			query_sketch->bd_plus[j] = 0;
		}
		int l;
		for(l = j - 1; l >= 0 && query_sketch->bd_plus[query_sketch->idx[l]] > query_sketch->bd_plus[j]; l--) {
			query_sketch->idx[l + 1] = query_sketch->idx[l];
		}
		query_sketch->idx[l + 1] = j;
		// 表関数の更新（bd[j] を対応するところに加える）
		int tn = j / 8, bp = j % 8; // 表関数の変更すべき表番号が tn，そのビット位置が bp
		for(int n = 0; n < 256; n++) {
			if(n & (1 << bp)) {
				query_sketch->tbl[tn][n] += query_sketch->bd_plus[j];
			}
		}
	}
}
#endif

// for score_p
void set_query_sketch_p(struct_query_sketch *qs, query_type *query, pivot_type *pivot, double p)
{
	qs->query = *query;
	qs->sketch = 0;
	int tbl_size = 4;
	for(int p = 0; p < tbl_size; p++) {
		for(int n = 0; n < 256; n++) {
			qs->tbl[p][n] = 0;
		}
	}
	for(int j = 0; j < PJT_DIM; j++) qs->idx[j] = j;
	for(int j = 0; j < PJT_DIM; j++) {
		#ifdef USE_PD_SKETCH
		dist_type dist = dist_pivot_L2(query->ftr, pivot->p[j], pivot->num_axis[j], pivot->axis[j]);
		#elif defined(PARTITION_TYPE_PQBP)
		dist_type dist = PART_DISTANCE_PIVOT(query->ftr, pivot->p[j], PART_START(j), PART_DIM(j));
		#elif defined(PARTITION_TYPE_SQBP) || defined(PARTITION_TYPE_CSQBP)
		dist_type dist = sub_dist_L2(query->ftr, pivot, j);
		#else
		dist_type dist = DISTANCE_PIVOT(query->ftr, pivot->p[j], FTR_DIM);
		#endif
		WRITE_BIT(j, dist < pivot->r[j], &qs->sketch);

		#ifndef SQRT_FTR
			if(p != DBL_MAX) {
				qs->bd[j] = pow(fabs(sqrt(dist) - sqrt(pivot->r[j])), p); // score_p
			} else {
				qs->bd[j] = fabs(sqrt(dist) - sqrt(pivot->r[j])); // score_inf
			}
		#else
			if(p != DBL_MAX) {
				qs->bd[j] = pow(fabs((int)dist - (int)pivot->r[j]), p); // score_p
			} else {
				qs->bd[j] = fabs((int)dist - (int)pivot->r[j]); // score_inf
			}
		#endif

		// bd をソート（idxを用いた相対ソート）（挿入法：bd[j] を挿入）
		int l;
		for(l = j - 1; l >= 0 && qs->bd[qs->idx[l]] > qs->bd[j]; l--) {
			qs->idx[l + 1] = qs->idx[l];
		}
		qs->idx[l + 1] = j;

		// 表関数の更新（bd[j] を対応するところに加える）
		int tn = j / 8, bp = j % 8; // 表関数の変更すべき表番号が tn，そのビット位置が bp
		for(int n = 0; n < 256; n++) {
			if(n & (1 << bp)) {
				if(p != DBL_MAX) {
					qs->tbl[tn][n] += qs->bd[j];
				} else {
					if(qs->tbl[tn][n] < qs->bd[j]) {
						qs->tbl[tn][n] = qs->bd[j];
					}
				}
			}
		}
	}

}

dist_type priority(sketch_type s, struct_query_sketch *q)
{
	#ifndef SCORE_INF// score_1 or score_2
		sketch_type d = s ^ q->sketch;
		return q->tbl[0][d & 0xff] + q->tbl[1][(d >> 8) & 0xff] + q->tbl[2][(d >> 16) & 0xff] + q->tbl[3][d >> 24];
	#else // score_inf
		sketch_type d = s ^ q->sketch;
		dist_type inf = 0;

		inf = q->tbl[0][d & 0xff];
		inf = inf >= q->tbl[1][(d >>  8) & 0xff] ? inf : q->tbl[1][(d >> 8 ) & 0xff];
		inf = inf >= q->tbl[2][(d >> 16) & 0xff] ? inf : q->tbl[2][(d >> 16) & 0xff];
		inf = inf >= q->tbl[3][(d >> 24) & 0xff] ? inf : q->tbl[3][(d >> 24) & 0xff];
		return inf;
	#endif	
}

dist_type priority_inf(sketch_type s, struct_query_sketch *q)
{
	sketch_type d = s ^ q->sketch;
	dist_type inf = 0;

	inf = q->tbl[0][d & 0xff];
	inf = inf >= q->tbl[1][(d >>  8) & 0xff] ? inf : q->tbl[1][(d >> 8 ) & 0xff];
	inf = inf >= q->tbl[2][(d >> 16) & 0xff] ? inf : q->tbl[2][(d >> 16) & 0xff];
	inf = inf >= q->tbl[3][(d >> 24) & 0xff] ? inf : q->tbl[3][(d >> 24) & 0xff];
	return inf;

}

dist_type hamming(sketch_type s, sketch_type t)
{
	sketch_type d = s ^ t;
	return bit_count(d);
}

#ifndef WITHOUT_FILTERING
void filtering_by_sequential_search_using_kNN_buffer(struct_query_sketch *qs, int num_data, sketch_type sketch[], kNN_buffer *buff, int data_num[], int num_candidates)
{
	static answer_type *ans = NULL;
	if(ans == NULL) ans = (answer_type *)malloc(sizeof(answer_type) * num_data);
	// 優先順位（priority）の計算だけ並列処理する（kNN_buffer による選択はフィルタリングでは効果がないので，直列処理する）
	#ifdef _OPENMP
	#if defined(NUM_THREADS) && NUM_THREADS > 1
	int nt = (num_data < NUM_THREADS ? num_data : NUM_THREADS);
	omp_set_num_threads(nt);
	#pragma omp parallel for
	#endif
	#endif
	for(int i = 0; i < num_data; i++) {
		ans[i].data_num = i;
		ans[i].dist = priority(sketch[i], qs);
	}
	if(num_candidates == 0) {
		return; // scoring (calculation of priorities) only
	}
	dist_type k_nearest = make_empty_kNN_buffer(buff);
	for(int i = 0; i < num_data; i++) {
		if(ans[i].dist < k_nearest) k_nearest = push_kNN_buffer(&ans[i], buff);
	}
	k_nearest = flush_kNN_buffer(buff);
	for(int i = 0; i < num_candidates; i++) {
		data_num[i] = buff->buff[i].data_num;
	}
}

void filtering_by_sequential_search_using_quick_select_k(struct_query_sketch *qs, int num_data, sketch_type sketch[], dist_type score[], int data_num[], int num_candidates)
{
	#ifdef _OPENMP
	#if defined(NUM_THREADS) && NUM_THREADS > 1
	int nt = (num_data < NUM_THREADS ? num_data : NUM_THREADS);
	omp_set_num_threads(nt);
	#pragma omp parallel for
	#endif
	#endif
	for(int i = 0; i < num_data; i++) {
		data_num[i] = i;
		score[i] = priority(sketch[i], qs);
	}
	if(num_candidates == 0) {
		return; // scoring (calculation of priorities) only
	}
	#ifdef QUICK_SELECT
	quick_select_k(data_num, score, 0, num_data - 1, num_candidates); 
	#else
	quick_sort(data_num, score, 0, num_data - 1); 
	#endif
}

interval_list *new_interval_list(unsigned int nt, unsigned int size)
{
	interval_list *il = MALLOC(sizeof(interval_list));
	il->nt 		= nt;
	il->size 	= size;
	il->lg 		= MALLOC(sizeof(unsigned int) * nt);
	il->list 	= MALLOC(sizeof(interval) * nt * size);

	il->lg[0] = 0;
	for(int i = 0; i < nt; i++) {
		il->list[i * size].start = 0;
	}

	return il;
}

void realloc_interval_list(interval_list *ivl, unsigned int size) {
	FREE(ivl->list, sizeof(interval) * ivl->nt * ivl->size);
	ivl->size = size;
	ivl->list = MALLOC(sizeof(interval) * ivl->nt * ivl->size);
}

vlist *new_vlist(int size, int step) {
	vlist *vl = MALLOC(sizeof(vlist));
	vl->size = size;
	vl->step = step;
	vl->num_list = 0;
	vl->num_data = 0;
	vl->elm = MALLOC(sizeof(interval) * size);
	return vl;
}

void add_vlist(vlist *vl, interval i) {
	if(vl->num_list >= vl->size) {
		int psize = vl->size;
		interval *pelm = vl->elm;
		vl->size += vl->step;
		vl->elm = MALLOC(sizeof(interval) * vl->size);
		memcpy(vl->elm, pelm, sizeof(interval) * psize);
		FREE(pelm, sizeof(interval) * psize);
	}
	vl->elm[vl->num_list++] = i;
	vl->num_data += i.run;
}

void makenull_vlist(vlist *vl) {
	vl->num_list = 0;
	vl->num_data = 0;
}

void enum_sub_dimension_0(sub_dimension sub, int begin, int remain, sub_dimension tab[], int *k, int K, int dim)
{
    if(remain!=0) {
    	int i;
        for(i = dim - remain; i >= begin; i--) {
        	sub_dimension sub2 = sub;
        	sub2.dim[sub2.num++] = dim - i - 1;
            enum_sub_dimension_0(sub2, i + 1, remain - 1, tab, k, K, dim);
            if(*k >= K) {
                return;
            }
        }
    } else {
		tab[(*k)++] = sub;
    	if(*k >= K) {
            return;
        }
    }
}

void rearrange_table(int nt, sub_dimension *tab, int tab_size)
{
	sub_dimension *temp = malloc(sizeof(sub_dimension) * tab_size);
	memcpy(temp, tab, sizeof(sub_dimension) * tab_size);
	int x = tab_size / nt;
	for(int i = 0; i < tab_size; i++) {
		tab[i / nt + i % nt * x] = temp[i];
	}
	free(temp);
}

int nextComb(int n)
{
	int x = n & -n;
	int y = x + n;
	return y | (((n & ~y) / x) >> 1);
}

// ハミング距離順の列挙のために要素数の昇順に同じ要素数内では辞書順に射影次元番号集合　{0, ... , dim -1} の部分集合を求める
sub_dimension *make_table_for_enumeration_hamming(int n, int dim)
{
	sub_dimension *table = MALLOC(sizeof(sub_dimension) * n);
	sub_dimension sub;
	int k = 0;

	sub.num = 0;
	table[k++] = sub; // 先頭に空集合を入れる．
    for(int n_elm = 1; n_elm <= dim  && k < n - 1; n_elm++) {
		for(int bit = (1 << n_elm) - 1; bit < (1 << dim) && k < n; bit = nextComb(bit)) {
			sub.num = 0;
			int temp = bit;
			while(temp != 0) {
				int p =  lsb_pos((unsigned)temp);
				sub.dim[sub.num++] = p;
				temp ^= (1 << p);
			}
			table[k++] = sub;
		}
	}
	sub.num = 0;
	while(k < n) {
		table[k++] = sub; // 最後に空集合を入れる．
	}
	return table;
}

void diff(sub_dimension *a, sub_dimension *b, sub_dimension *c) {
	int i, j;
	sub_dimension temp;
	temp.num = 0;
	for(i = 0; i < PJT_DIM; i++) {
		int x = 0, y = 0;
		for(j = 0; j < a->num; j++) {
			if(a->dim[j] == i) {
				x = 1;
			}
		}
		for(j = 0; j < b->num; j++) {
			if(b->dim[j] == i) {
				y = 1;
			}
		}
		if(x != y) {
			temp.dim[temp.num++] = i;
		}
	}
	*c = temp;
}

#define USE_MU
// #define USE_AVE_ALL

#ifdef USE_MU
inline sketch_type Mu(int b, int lg, int idx[], int i, sub_dimension *sd)
{
	sketch_type mask = 0;
	for(int m = 0; m < sd[i].num; m++) {
		int j = sd[i].dim[m];
		mask |= 1 << idx[b + j];
	}
	return mask;
}
#endif

#ifndef PARA_ENUM_INF
#define PARA_ENUM_INF 0
#endif

#ifdef WITH_ENUM_TABLE

	#ifdef ENUM_CHECK
	int max_H = 0, max_n = 0;

	void reset_max_H(void) {
		max_H = max_n = 0;
	}
	#endif

	#if PARA_ENUM_INF == 0  
	// スケッチ列挙(Hamming)によるフィルタリング（バケット（配列 idx と bkt）利用）（部分集合列挙の表を利用する）（single-threasd）
	int filtering_by_sketch_enumeration_hamming(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates)
	{
		double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
		int num_sketch = num_candidates / ave_num; // 求めるスケッチの個数（期待値）
		int num_sketch_hamm = num_sketch * 5 /* FACTOR_INF */; // Hamming距離順の列挙によって求めるスケッチの個数（十分多めにしておく）
		#ifdef SELECT_SUM
			static sketch_with_priority_num *buff = NULL;
			static int buff_size = 0;
			if(buff == NULL) {
				fprintf(stderr, "malloc answer buffer\n");
				buff = MALLOC(sizeof(sketch_with_priority_num) * num_sketch_hamm);
				fprintf(stderr, "malloc answer buffer OK: num_sketch = %d, num_sketch_hamm = %d\n", num_sketch, num_sketch_hamm);
				buff_size = num_sketch_hamm;
			} else if(buff_size < num_sketch_hamm) {
				fprintf(stderr, "free and malloc answer buffer\n");
				FREE(buff, sizeof(sketch_with_priority_num) * buff_size);
				buff = MALLOC(sizeof(sketch_with_priority_num) * num_sketch_hamm);
				fprintf(stderr, "malloc answer buffer OK: num_sketch = %d, num_sketch_hamm = %d\n", num_sketch, num_sketch_hamm);
				buff_size = num_sketch_hamm;
			}
		#else
			static answer_type *buff = NULL;
			static int buff_size = 0;
			if(buff == NULL) {
				fprintf(stderr, "malloc answer buffer\n");
				buff = MALLOC(sizeof(answer_type) * num_sketch_hamm);
				fprintf(stderr, "malloc answer buffer OK: num_sketch = %d, num_sketch_hamm = %d\n", num_sketch, num_sketch_hamm);
				buff_size = num_sketch_hamm;
			} else if(buff_size < num_sketch_hamm) {
				fprintf(stderr, "free and malloc answer buffer\n");
				FREE(buff, sizeof(answer_type) * buff_size);
				buff = MALLOC(sizeof(answer_type) * num_sketch_hamm);
				fprintf(stderr, "malloc answer buffer OK: num_sketch = %d, num_sketch_hamm = %d\n", num_sketch, num_sketch_hamm);
				buff_size = num_sketch_hamm;
			}
		#endif
		static sub_dimension *table = NULL;
		static int table_size = 0;
		#ifdef ENUM_DIM
		int enum_dim = ENUM_DIM;
		#else
		int enum_dim = PJT_DIM - 22; // ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
		#endif
		if(table == NULL) {
			table_size = (1 << enum_dim) + 1;
			table = make_table_for_enumeration_hamming(table_size, enum_dim);
			#ifdef USE_DIFF_TABLE
			for(int i = 0; i < table_size - 1; i++) {
				diff(&table[i], &table[i + 1], &table[i]);
				if(table[i + 1].num == 0) break; // 次が空集合になったら終了
			}
			#endif
		} else if(table_size < 8 * num_candidates) { // この部分は不要？
			FREE(table, sizeof(sub_dimension) * table_size);
			table_size = 8 * num_candidates;
			table = make_table_for_enumeration_hamming(table_size, enum_dim);
			#ifdef USE_DIFF_TABLE
			for(int i = 0; i < table_size - 1; i++) {
				diff(&table[i], &table[i + 1], &table[i]);
				if(table[i + 1].num == 0) break; // 次が空集合になったら終了
			}
			#endif
		}

		static int first = 1;
		if(first) {
			#ifdef USE_DIFF_TABLE
			fprintf(stderr, "enum_hamm_with_after_selection single-thread. using diff table, enum_dim = %d\n", enum_dim);
			#else
			fprintf(stderr, "enum_hamm_with_after_selection single-thread. using table, enum_dim = %d\n", enum_dim);
			#endif
			first = 0;
		}

		#define TABLE_SPP
		// tableを用いた列挙では不足するときに，追加のビットパターンを求めるための表
		#ifdef TABLE_SPP
		static sub_dimension *table_spp = NULL;
		#ifdef SPP_BIT
		int spp_bit = SPP_BIT; // 追加のビット数
		#else
		int spp_bit = 17; // 追加のビット数
		#endif
		static int table_spp_size = 0;
		if(table_spp == NULL) {
			table_spp_size = (1 << spp_bit) + 1;
			table_spp = make_table_for_enumeration_hamming(table_spp_size, spp_bit);
			fprintf(stderr, "make table for spplementary enumeation: spp_bit = %d, table_spp_size = %d\n", spp_bit, table_spp_size);
			#ifdef USE_DIFF_TABLE
			for(int i = 0; i < table_spp_size - 1; i++) {
				diff(&table_spp[i], &table_spp[i + 1], &table_spp[i]);
				if(table_spp[i + 1].num == 0) break; // 次が空集合になったら終了
			}
			#endif
		}
		#endif

		#ifndef WITHOUT_IDX
		int *bd_idx = qs->idx;
		#endif
		int *bkt = bucket->bkt;

		int n;
		sketch_type sk = qs->sketch;
		sketch_type base_mask = 0, base_mask2 = 0;
		#ifdef TABLE_SPP
		int n_spp = 1; // つぎに使用する追加パターン番号（0 のものは空集合なので，1から使用する）
		int n_spp2 = 1; // 追加パターンも使い切ったら，さらに追加する．
		#endif

		int num_enum_data = 0; // 列挙されたスケッチを持つデータ総数
		int num_nonempty = 0; // 列挙によって求めた空でないスケッチの個数
		for(n = 0; num_enum_data < num_candidates * FACTOR_INF && num_nonempty < buff_size; n++) { // table を使って mask を作って，つぎのスケッチを列挙する
			sketch_type mask = base_mask;
			for(int m = 0; m < table[n].num; m++) {
				#ifdef WITHOUT_IDX // bucket->idx を用いた相対位置によるビット操作を行わない
				mask |= (1 << table[n].dim[m]);
				#else
				mask |= (1 << bd_idx[(int)table[n].dim[m]]);
				#endif
			}
			#ifdef USE_DIFF_TABLE
			sk = sk ^ mask;
			#else
			sk = qs->sketch ^ mask;
			#endif
			if(bkt[sk + 1] > bkt[sk]) { // sk が空でなければ，buffに追加
				#ifdef SELECT_SUM
				buff[num_nonempty++] = (sketch_with_priority_num){sk, priority(sk, qs), bkt[sk + 1] - bkt[sk]};
				num_enum_data += bkt[sk + 1] - bkt[sk];
				#else
				buff[num_nonempty++] = (answer_type){sk, priority(sk, qs)};
				num_enum_data += bkt[sk + 1] - bkt[sk];
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
						#ifdef WITHOUT_IDX // bucket->idx を用いた相対位置によるビット操作を行わない
						base_mask2 |= (1 << (table_spp[n_spp2].dim[m] + enum_dim + spp_bit));
						#else
						base_mask2 |= (1 << bd_idx[(int)table_spp[n_spp2].dim[m] + enum_dim + spp_bit]);
						#endif
					}
				}
				base_mask = base_mask2;
				for(int m = 0; m < table_spp[n_spp].num; m++) {
					#ifdef WITHOUT_IDX // bucket->idx を用いた相対位置によるビット操作を行わない
					base_mask |= (1 << (table_spp[n_spp].dim[m] + enum_dim));
					#else
					base_mask |= (1 << bd_idx[(int)table_spp[n_spp].dim[m] + enum_dim]);
					#endif
				}
				n = -1; // forの再初期化で n++ となって，0になる．
				n_spp++;
				#else
				break;
				#endif
			}
		}

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

		#ifdef SELECT_SUM
			int num_selected_sketches, check_sum = 0;
			num_selected_sketches =quick_select_sum_k_sketch_with_priority_num(buff, 0, num_nonempty - 1, num_candidates);
			int k = 0, i; // 出力したデータ番号の個数
			for(i = 0; k < num_candidates; i++) {
				sketch_type sk = (sketch_type)(buff[i].sk);
				for(int j = bkt[sk]; j < bkt[sk + 1] && k < num_candidates; j++, k++) {
					#ifdef DATA_NUM_IN_SKETCH_ORDER
					data_num[k] = j;
					#else
					data_num[k] = bucket->idx[j];
					#endif
				}
			}
		#else
			int sel = (num_sketch < num_nonempty ? num_sketch : num_nonempty);
			quick_select_k_answer(buff, 0, num_nonempty - 1, sel * 0.4);
			int k = 0, i; // 出力したデータ番号の個数
			for(i = 0; k < num_candidates && i < sel * 0.4; i++) {
				sketch_type sk = (sketch_type)(buff[i].data_num);
				for(int j = bkt[sk]; j < bkt[sk + 1] && k < num_candidates; j++, k++) {
					#ifdef DATA_NUM_IN_SKETCH_ORDER
					data_num[k] = j;
					#else
					data_num[k] = bucket->idx[j];
					#endif
				}
			}
			if(k == num_candidates) return num_candidates;
			quick_select_k_answer(buff, sel * 0.4 - 1, num_nonempty - 1, sel * 0.3);
			for(; k < num_candidates && i < num_nonempty; i++) {
				sketch_type sk = (sketch_type)(buff[i].data_num);
				for(int j = bkt[sk]; j < bkt[sk + 1] && k < num_candidates; j++, k++) {
					#ifdef DATA_NUM_IN_SKETCH_ORDER
					data_num[k] = j;
					#else
					data_num[k] = bucket->idx[j];
					#endif
				}
			}
		#endif

		return k;
	}

	// single-thread
	#ifndef USE_MU
	int filtering_by_sketch_enumeration_hamming_interval(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates)
	{
		double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
		int num_sketch = num_candidates / ave_num; // 求めるスケッチの個数（期待値）
		int num_sketch_hamm = num_sketch * 5 /* FACTOR_INF */; // Hamming距離順の列挙によって求めるスケッチの個数（十分多めにしておく）
		static sketch_with_priority_num *buff = NULL;
		static int buff_size = 0;
		if(buff == NULL) {
			fprintf(stderr, "malloc answer buffer\n");
			buff = MALLOC(sizeof(sketch_with_priority_num) * num_sketch_hamm);
			fprintf(stderr, "malloc answer buffer OK: num_sketch = %d, num_sketch_hamm = %d\n", num_sketch, num_sketch_hamm);
			buff_size = num_sketch_hamm;
		} else if(buff_size < num_sketch_hamm) {
			fprintf(stderr, "free and malloc answer buffer\n");
			FREE(buff, sizeof(sketch_with_priority_num) * buff_size);
			buff = MALLOC(sizeof(sketch_with_priority_num) * num_sketch_hamm);
			fprintf(stderr, "malloc answer buffer OK: num_sketch = %d, num_sketch_hamm = %d\n", num_sketch, num_sketch_hamm);
			buff_size = num_sketch_hamm;
		}

		static sub_dimension *table = NULL;
		static int table_size = 0;
		#ifdef ENUM_DIM
		int enum_dim = ENUM_DIM;
		#else
		int enum_dim = PJT_DIM - 22; // ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
		#endif
		if(table == NULL) {
			table_size = (1 << enum_dim) + 1;
			table = make_table_for_enumeration_hamming(table_size, enum_dim);
			#ifdef USE_DIFF_TABLE
			for(int i = 0; i < table_size - 1; i++) {
				diff(&table[i], &table[i + 1], &table[i]);
				if(table[i + 1].num == 0) break; // 次が空集合になったら終了
			}
			#endif
		} else if(table_size < 8 * num_candidates) { // この部分は不要？
			FREE(table, sizeof(sub_dimension) * table_size);
			table_size = 8 * num_candidates;
			table = make_table_for_enumeration_hamming(table_size, enum_dim);
			#ifdef USE_DIFF_TABLE
			for(int i = 0; i < table_size - 1; i++) {
				diff(&table[i], &table[i + 1], &table[i]);
				if(table[i + 1].num == 0) break; // 次が空集合になったら終了
			}
			#endif
		}

		static int first = 1;
		if(first) {
			#ifdef USE_DIFF_TABLE
			fprintf(stderr, "enum_hamm_with_after_selection_interval single-thread. using diff table, enum_dim = %d\n", enum_dim);
			#else
			fprintf(stderr, "enum_hamm_with_after_selection_interval single-thread. using table, enum_dim = %d\n", enum_dim);
			#endif
			first = 0;
		}

		#define TABLE_SPP
		// tableを用いた列挙では不足するときに，追加のビットパターンを求めるための表
		#ifdef TABLE_SPP
		static sub_dimension *table_spp = NULL;
		#ifdef SPP_BIT
		int spp_bit = SPP_BIT; // 追加のビット数
		#else
		int spp_bit = 17; // 追加のビット数
		#endif
		static int table_spp_size = 0;
		if(table_spp == NULL) {
			table_spp_size = (1 << spp_bit) + 1;
			table_spp = make_table_for_enumeration_hamming(table_spp_size, spp_bit);
			fprintf(stderr, "make table for spplementary enumeation: spp_bit = %d, table_spp_size = %d\n", spp_bit, table_spp_size);
			#ifdef USE_DIFF_TABLE
			for(int i = 0; i < table_spp_size - 1; i++) {
				diff(&table_spp[i], &table_spp[i + 1], &table_spp[i]);
				if(table_spp[i + 1].num == 0) break; // 次が空集合になったら終了
			}
			#endif
		}
		#endif

		#ifndef WITHOUT_IDX
		int *bd_idx = qs->idx;
		#endif
		int *bkt = bucket->bkt;

		int n;
		sketch_type sk = qs->sketch;
		sketch_type base_mask = 0, base_mask2 = 0;
		#ifdef TABLE_SPP
		int n_spp = 1; // つぎに使用する追加パターン番号（0 のものは空集合なので，1から使用する）
		int n_spp2 = 1; // 追加パターンも使い切ったら，さらに追加する．
		#endif

		int num_enum_data = 0; // 列挙されたスケッチを持つデータ総数
		int num_nonempty = 0; // 列挙によって求めた空でないスケッチの個数
//			int num_all_low = 0;
//			#ifdef LOOP_CONTROL_BY_NUM_SKETCHES
//			#if defined(INCLUDE_EMPTY_SKETCHES)
//			int enum_sketches = 0;
//			for(n = 0; enum_sketches < num_sketch /* && k < num_data_thread * FACTOR_INF */; n++)
//			#else
//			for(n = 0; num_enum_data < num_candidates * FACTOR_INF2 && num_nonempty < num_sketch; n++) 
//			#endif
//			#else
//			#define BUCKET_SIZE
		#ifdef BUCKET_SIZE
		int visited_sketches = 0;
		double visited_sum = 0, visited_sum2 = 0; 
		#endif

		for(n = 0; num_enum_data < num_candidates * FACTOR_INF && num_nonempty < buff_size; n++) 
//			#endif
		{
			// table を使って mask を作って，つぎのスケッチを列挙する．
			sketch_type mask = base_mask;
			for(int m = 0; m < table[n].num; m++) {
				#ifdef WITHOUT_IDX // bucket->idx を用いた相対位置によるビット操作を行わない
				int bp = table[n].dim[m];
				#else
				int bp = bd_idx[(int)table[n].dim[m]];
				#endif
				mask |= (1 << bp);
			}
			#ifdef USE_DIFF_TABLE
			sk = sk ^ mask;
			#else
			sk = qs->sketch ^ mask;
			#endif
			if(bkt[sk + 1] > bkt[sk]) { // sk が空でなければ，buffに追加
				buff[num_nonempty++] = (sketch_with_priority_num){sk, priority(sk, qs), bkt[sk + 1] - bkt[sk]};
				num_enum_data += bkt[sk + 1] - bkt[sk];
			}
			#ifdef BUCKET_SIZE
			visited_sketches++;
			visited_sum += (bkt[sk + 1] - bkt[sk]);
			visited_sum2 += (bkt[sk + 1] - bkt[sk]) * (bkt[sk + 1] - bkt[sk]);
			#endif
			if(table[n + 1].num == 0) { // 用意したパターンがなくなった．
				// fprintf(stderr, "LOWのパターンを使い切った．n = %d, count = %d\n", n, ++num_all_low);
				#ifdef TABLE_SPP
				if(table_spp[n_spp].num == 0) { // 追加分のパターンも使い切った．
					// fprintf(stderr, "追加分のパターンを使い切った．Hit enter !\n"); getchar();
					if(table_spp[n_spp2].num == 0) break; // 2回目の追加分のパターンも使い切った．
					n_spp = 0;
					n_spp2++;
					base_mask2 = 0;
					for(int m = 0; m < table_spp[n_spp2].num; m++) {
						#ifdef WITHOUT_IDX // bucket->idx を用いた相対位置によるビット操作を行わない
						int bp = table_spp[n_spp2].dim[m] + enum_dim + spp_bit;
						#else
						int bp = bd_idx[(int)table_spp[n_spp2].dim[m] + enum_dim + spp_bit];
						#endif
						base_mask2 |= (1 << bp);
					}
				}
				base_mask = base_mask2;
				for(int m = 0; m < table_spp[n_spp].num; m++) {
					#ifdef WITHOUT_IDX // bucket->idx を用いた相対位置によるビット操作を行わない
					int bp = table_spp[n_spp].dim[m] + enum_dim;
					#else
					int bp = bd_idx[(int)table_spp[n_spp].dim[m] + enum_dim];
					#endif
					base_mask |= (1 << bp);
				}
				n = -1; // forの再初期化で n++ となって，0になる．
				n_spp++;
				#else
				break;
				#endif
			}
		}
		#ifdef BUCKET_SIZE
		double ave = visited_sum / visited_sketches;
		double ave2 = visited_sum2 / visited_sketches;
		double stdev = sqrt(ave2 - ave * ave);
		printf("visited_sketches = %d, average = %.2lf, stdev = %.2lf\n", visited_sketches, ave, stdev);
		// getchar();
		#endif

		quick_select_sum_k_sketch_with_priority_num(buff, 0, num_nonempty - 1, num_candidates);
		int k = 0, num_sk = 0; // k = 出力したデータ番号の個数, num_sk = スケッチの個数
		for(int i = 0; k < num_candidates; i++) {
			sketch_type sk = (sketch_type)(buff[num_sk].sk);
			ivl->list[num_sk++] = (interval){bkt[sk], bkt[sk + 1] - bkt[sk]};
			k += bkt[sk + 1] - bkt[sk];
		}
		ivl->lg[0] = num_sk;
		return k;

	}
	#else // USE_MU
	int filtering_by_sketch_enumeration_hamming_interval(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates)
	{
		double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
		int num_sketch = num_candidates / ave_num; // 求めるスケッチの個数（期待値）
		int num_sketch_hamm = num_sketch * 5 /* FACTOR_INF */; // Hamming距離順の列挙によって求めるスケッチの個数（十分多めにしておく）
		static sketch_with_priority_num *buff = NULL;
		static int buff_size = 0;
		if(buff == NULL) {
			fprintf(stderr, "malloc answer buffer\n");
			buff = MALLOC(sizeof(sketch_with_priority_num) * num_sketch_hamm);
			fprintf(stderr, "malloc answer buffer OK: num_sketch = %d, num_sketch_hamm = %d\n", num_sketch, num_sketch_hamm);
			buff_size = num_sketch_hamm;
		} else if(buff_size < num_sketch_hamm) {
			fprintf(stderr, "free and malloc answer buffer\n");
			FREE(buff, sizeof(sketch_with_priority_num) * buff_size);
			buff = MALLOC(sizeof(sketch_with_priority_num) * num_sketch_hamm);
			fprintf(stderr, "malloc answer buffer OK: num_sketch = %d, num_sketch_hamm = %d\n", num_sketch, num_sketch_hamm);
			buff_size = num_sketch_hamm;
		}

		static sub_dimension *table = NULL;
		static int table_size = 0;
		#ifdef ENUM_DIM
		int enum_dim = ENUM_DIM;
		#else
		int enum_dim = PJT_DIM - 22; // ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
		#endif
		if(table == NULL) {
			table_size = (1 << enum_dim) + 1;
			table = make_table_for_enumeration_hamming(table_size, enum_dim);
			#ifdef USE_DIFF_TABLE
			for(int i = 0; i < table_size - 1; i++) {
				diff(&table[i], &table[i + 1], &table[i]);
				if(table[i + 1].num == 0) break; // 次が空集合になったら終了
			}
			#endif
		} else if(table_size < 8 * num_candidates) { // この部分は不要？
			FREE(table, sizeof(sub_dimension) * table_size);
			table_size = 8 * num_candidates;
			table = make_table_for_enumeration_hamming(table_size, enum_dim);
			#ifdef USE_DIFF_TABLE
			for(int i = 0; i < table_size - 1; i++) {
				diff(&table[i], &table[i + 1], &table[i]);
				if(table[i + 1].num == 0) break; // 次が空集合になったら終了
			}
			#endif
		}

		static int first = 1;
		if(first) {
			#ifdef USE_DIFF_TABLE
			fprintf(stderr, "enum_hamm_with_after_selection_interval single-thread. using diff table, enum_dim = %d\n", enum_dim);
			#else
			fprintf(stderr, "enum_hamm_with_after_selection_interval single-thread. using table, enum_dim = %d\n", enum_dim);
			#endif
			first = 0;
		}

		#define TABLE_SPP
		// tableを用いた列挙では不足するときに，追加のビットパターンを求めるための表
		#ifdef TABLE_SPP
		static sub_dimension *table_spp = NULL;
		#ifdef SPP_BIT
		int spp_bit = SPP_BIT; // 追加のビット数
		#else
		int spp_bit = 17; // 追加のビット数
		#endif
		static int table_spp_size = 0;
		if(table_spp == NULL) {
			table_spp_size = (1 << spp_bit) + 1;
			table_spp = make_table_for_enumeration_hamming(table_spp_size, spp_bit);
			fprintf(stderr, "make table for spplementary enumeation: spp_bit = %d, table_spp_size = %d\n", spp_bit, table_spp_size);
			#ifdef USE_DIFF_TABLE
			for(int i = 0; i < table_spp_size - 1; i++) {
				diff(&table_spp[i], &table_spp[i + 1], &table_spp[i]);
				if(table_spp[i + 1].num == 0) break; // 次が空集合になったら終了
			}
			#endif
		}
		#endif

		#ifndef WITHOUT_IDX
		int *bd_idx = qs->idx;
		#endif
		int *bkt = bucket->bkt;

		#ifdef INF_AND_HAMM
		sketch_type mu_low = qs->sketch;
		#else
		sketch_type mu_low = 0;
		#endif
		sketch_type mu_add = 0;		// add-bit 部分のマスクパターン
		int num_enum_data = 0;		// 列挙されたスケッチを持つデータ総数
		int num_nonempty = 0;		// 列挙によって求めた空でないスケッチの個数
		int low_c = 0, add_c = 0;	// low-bit と add-bit のカウンタ

		while(num_enum_data < num_candidates * FACTOR_INF && num_nonempty < buff_size) {
			#ifndef INF_AND_HAMM
			mu_low = Mu(0, enum_dim, bd_idx, low_c, table);
			sketch_type sk = qs->sketch ^ mu_low ^ mu_add;
			#else
			sketch_type sk = mu_low ^ mu_add;
			#endif
			if(bkt[sk + 1] > bkt[sk]) { // sk が空でなければ，buffに追加
				buff[num_nonempty++] = (sketch_with_priority_num){sk, priority(sk, qs), bkt[sk + 1] - bkt[sk]};
				num_enum_data += bkt[sk + 1] - bkt[sk];
			}
			#ifdef INF_AND_HAMM
			mu_low = mu_low ^ (1 << bd_idx[bit_count(low_c ^ (low_c + 1)) - 1]);
			#endif
			low_c++;
			if(low_c == (1 << enum_dim)) {
				low_c = 0;
				#ifdef INF_AND_HAMM
				mu_low = qs->sketch;
				#endif
				add_c++;
				mu_add = Mu(enum_dim, spp_bit, bd_idx, add_c, table_spp);
			}
		}
		quick_select_sum_k_sketch_with_priority_num(buff, 0, num_nonempty - 1, num_candidates);
		int k = 0, num_sk = 0; // k = 出力したデータ番号の個数, num_sk = スケッチの個数
		for(int i = 0; k < num_candidates; i++) {
			sketch_type sk = (sketch_type)(buff[num_sk].sk);
			ivl->list[num_sk++] = (interval){bkt[sk], bkt[sk + 1] - bkt[sk]};
			k += bkt[sk + 1] - bkt[sk];
		}
		ivl->lg[0] = num_sk;
		return k;
	}
	#endif // USE_MU
	#else // PARA_ENUM_INF > 0 (multi-thread)
		#ifndef SELECT_BY_SINGLE
			#ifndef SELECT_BY_PARA_MERGE
			// スケッチ列挙(Hamming)によるフィルタリング（バケット（配列 idx と bkt）利用）（部分集合列挙の表を利用する）（multi-threasd）
			int filtering_by_sketch_enumeration_hamming(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates)
			{
				int n = PARA_ENUM_INF;
				int nt = (1 << n); // スレッド数
				#ifdef _OPENMP
				omp_set_num_threads(nt);
				#endif

				// ハミング距離順の列挙のための部分集合の表を用意する．
				static sub_dimension *table = NULL;
				static int table_size = 0;
				#ifdef ENUM_DIM
				int enum_dim = ENUM_DIM;
				#else
				int enum_dim = PJT_DIM - 22; // ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
				#endif
				if(table == NULL) {
					table_size = (1 << enum_dim) + nt;
					table = make_table_for_enumeration_hamming(table_size, enum_dim);
					rearrange_table(nt, table, table_size);
				fprintf(stderr, "made table for hamming enumeration: enum_dim = %d\n", enum_dim);
				use_system("VmSize");
				}

				static int first = 1;
				if(first) {
					#ifdef SELECT_SUM
					fprintf(stderr, "(2) ");
					#else
					fprintf(stderr, "(1) ");
					#endif
					#ifdef USE_DIFF_TABLE
					fprintf(stderr, "enum_hamm multi-thread (%d-thread). using diff table, enum_dim = %d\n", nt, enum_dim);
					#else
					fprintf(stderr, "enum_hamm multi-thread (%d-thread). using table, enum_dim = %d\n", nt, enum_dim);
					#endif
					first = 0;
				use_system("VmSize");
				}

				#define TABLE_SPP
				// tableを用いた列挙では不足するときに，追加のビットパターンを求めるための表
				#ifdef TABLE_SPP
				static sub_dimension *table_spp = NULL;
				#ifdef SPP_BIT
				int spp_bit = SPP_BIT; // 追加のビット数
				#else
				int spp_bit = 17; // 追加のビット数
				#endif
				static int table_spp_size = 0;
				if(table_spp == NULL) {
					table_spp_size = (1 << spp_bit) + 1;
					table_spp = make_table_for_enumeration_hamming(table_spp_size, spp_bit);
					fprintf(stderr, "make table for spplementary enumeation: spp_bit = %d, table_spp_size = %d\n", spp_bit, table_spp_size);
					#ifdef USE_DIFF_TABLE
					for(int i = 0; i < table_spp_size - 1; i++) {
						diff(&table_spp[i], &table_spp[i + 1], &table_spp[i]);
						if(table_spp[i + 1].num == 0) break; // 次が空集合になったら終了
					}
					#endif
				fprintf(stderr, "made table for hamming enumeration: spp_bit = %d\n", spp_bit);
				use_system("VmSize");
				}
				#endif

				#ifndef WITHOUT_IDX
				int *bd_idx = qs->idx;
				#endif
				int *bkt = bucket->bkt;

				// num_candidates  = 求めるデータ数（全体）
				// data_num[] = 求めたデータ番号を格納する配列
				// table[] = 部分集合を要素数の昇順に列挙するパターンを格納した配列
				// bd_idx[] = 相対位置でビット操作をするための距離下限の順位表
				// bkt[] = バケット表（bkt[s] = データをスケッチ順に並べたときに，スケッチ S のデータの先頭位置）
				// nt = スレッド数
				int num_data_thread = num_candidates / nt;	// スレッド求めるデータ番号数
				int t;										// スレッド番号
				int m;  									// スレッドが列挙したスケッチ数（パターン番号）
				int num_nonempty;							// スレッドが実際に列挙した空でないスケッチ数
				int k;  									// スレッドが求めたデータ番号数
				int table_size_thread = table_size / nt;	// スレッドが用いる列挙する部分集合のパターン配列の大きさ
				sub_dimension *q;							// スレッドが用いる列挙する部分集合のパターン配列
				int *d;										// スレッドが求めたデータ番号を格納する配列
				sketch_type base_mask = 0, base_mask2 = 0;	// 追加パターンで作成する mask
				int n_spp, n_spp2;							// 追加で使用するパターン番号
				sketch_type sk, mask;

				double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
				int num_sketch = num_candidates / ave_num; // 求めるスケッチの個数
				int num_sketch_thread = num_sketch / nt * 5; // 各スレッドが求めるスケッチの個数（十分多めにしておく）
				static int buff_pool_size = 0;
				#ifdef SELECT_SUM
					static sketch_with_priority_num *buff_pool = NULL;
					if(buff_pool == NULL || buff_pool_size < num_sketch_thread * nt) {
						if(buff_pool != NULL) {
							FREE(buff_pool, sizeof(sketch_with_priority_num) * buff_pool_size);
						}
						buff_pool_size = num_sketch_thread * nt;
						fprintf(stderr, "malloc answer buffer, ");
						buff_pool = MALLOC(sizeof(sketch_with_priority_num) * buff_pool_size);
						fprintf(stderr, "OK: buff_pool_size = %d, num_thread = %d\n", buff_pool_size, nt);
					}
					sketch_with_priority_num *buff;
				#else
					static answer_type *buff_pool = NULL;
					if(buff_pool == NULL || buff_pool_size < num_sketch_thread * nt) {
						if(buff_pool != NULL) {
							FREE(buff_pool, sizeof(answer_type) * buff_pool_size);
						}
						buff_pool_size = num_sketch_thread * nt;
						fprintf(stderr, "malloc answer buffer, ");
						buff_pool = MALLOC(sizeof(answer_type) * buff_pool_size);
						fprintf(stderr, "OK: buff_pool_size = %d, num_thread = %d\n", buff_pool_size, nt);
					}
					answer_type *buff;
				#endif

				#pragma omp parallel private(t, m, num_nonempty, k, q, d, base_mask, base_mask2, n_spp, n_spp2, buff, sk, mask)
				{
					t = omp_get_thread_num(); // スレッド番号の取得
					m = 0;
					num_nonempty = 0;
					k = 0;
					q = table + t * table_size_thread;
					d = data_num + t * num_data_thread;
					base_mask = base_mask2 = 0;
					n_spp = n_spp2 = 1;
					buff = buff_pool + t * num_sketch_thread;
					int enum_sketch = 0;

					for(m = 0; num_nonempty < num_sketch_thread && k < num_data_thread * FACTOR_INF; m++) {
						mask = base_mask;
						for(int j = 0; j < q[m].num; j++) {
							#ifndef WITHOUT_IDX
							mask |= (1 << bd_idx[(int)q[m].dim[j]]);
							#else
							mask |= (1 << (int)q[m].dim[j]);
							#endif
						}
						sk = qs->sketch ^ mask;
						// sk が空でなければ，buffに追加．ここでは，データ番号には展開しない
						if(bkt[sk + 1] > bkt[sk]) {
							#ifdef SELECT_SUM
							buff[num_nonempty++] = (sketch_with_priority_num){sk, priority(sk, qs), bkt[sk + 1] - bkt[sk]};
							#else
							buff[num_nonempty++] = (answer_type){sk, priority(sk, qs)};
							#endif
							k += bkt[sk + 1] - bkt[sk];
							if(enum_sketch == 0 && k >= num_candidates / nt) {
								enum_sketch = num_nonempty;
							}
						}
						if(q[m + 1].num == 0) { // 用意したパターンがなくなった．
							#ifdef TABLE_SPP
								if(table_spp[n_spp].num == 0) { // 追加分のパターンも使い切った．
									if(table_spp[n_spp2].num == 0) break; // 2回目の追加分のパターンも使い切った．
									n_spp = 0;
									n_spp2++;
									base_mask2 = 0;
									for(int m = 0; m < table_spp[n_spp2].num; m++) {
										#ifndef WITHOUT_IDX
										base_mask2 |= (1 << bd_idx[(int)table_spp[n_spp2].dim[m] + enum_dim + spp_bit]);
										#else
										base_mask2 |= (1 << ((int)table_spp[n_spp2].dim[m] + enum_dim + spp_bit));
										#endif
									}
								}
								base_mask = base_mask2;
								for(int m = 0; m < table_spp[n_spp].num; m++) {
									#ifndef WITHOUT_IDX
									base_mask |= (1 << bd_idx[(int)table_spp[n_spp].dim[m] + enum_dim]);
									#else
									base_mask |= (1 << ((int)table_spp[n_spp].dim[m] + enum_dim));
									#endif
								}
								m = -1; // forの再初期化で m++ となって，0になる．
								n_spp++;
							#else
								while(k < num_data_thread) {
									d[k++] = 0;
								}
							#endif
						}
					}
					#ifdef SELECT_SUM
					quick_select_sum_k_sketch_with_priority_num(buff, 0, num_nonempty - 1, num_data_thread);
					#else
					int sel = (num_sketch / nt < num_nonempty ? num_sketch / nt : num_nonempty);
					quick_select_k_answer(buff, 0, num_nonempty - 1, sel * 0.4);
					#endif
					int k = 0, i; // 出力したデータ番号の個数
					for(i = 0; k < num_data_thread /* && i < sel * 0.4 */; i++) {
						#ifdef SELECT_SUM
						sketch_type sk = (sketch_type)(buff[i].sk);
						#else
						sketch_type sk = (sketch_type)(buff[i].data_num);
						#endif
						for(int j = bkt[sk]; j < bkt[sk + 1] && k < num_data_thread; j++, k++) {
							#ifdef DATA_NUM_IN_SKETCH_ORDER
							d[k] = j;
							#else
							d[k] = bucket->idx[j];
							#endif
						}
					}
				}
				return num_candidates;
			}

			// データ数で制御．10分の1ずつを目標にして，少しずつ列挙して目標数を超えたら終了
			// スケッチ列挙(Hamming)によるフィルタリング（バケット（配列 idx と bkt）利用）（部分集合列挙の表を利用する）（multi-threasd）
			#ifndef USE_MU
			int filtering_by_sketch_enumeration_hamming_interval(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates)
			{
				// fprintf(stderr, "(5) filtering_by_sketch_enumeration_hamming （multi-threasd）\n"); exit(0);
				int n = PARA_ENUM_INF;
				int nt = (1 << n); // スレッド数
				#ifdef _OPENMP
				omp_set_num_threads(nt);
				#endif

				// ハミング距離順の列挙のための部分集合の表を用意する．
				static sub_dimension *table = NULL;
				static int table_size = 0;
				#ifdef ENUM_DIM
				int enum_dim = ENUM_DIM;
				#else
				int enum_dim = PJT_DIM - 22; // ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
				#endif
				if(table == NULL) {
					table_size = (1 << enum_dim) + nt;
					table = make_table_for_enumeration_hamming(table_size, enum_dim);
					rearrange_table(nt, table, table_size);
					fprintf(stderr, "made table for hamming enumeration: enum_dim = %d\n", enum_dim);
					use_system("VmSize");
				}

				static int first = 1;
				if(first) {
					#ifdef SELECT_SUM
					fprintf(stderr, "(2*) ");
					#else
					fprintf(stderr, "(1*+) ");
					#endif
					#ifdef USE_DIFF_TABLE
					fprintf(stderr, "enum_hamm_interval multi-thread (%d-thread). using diff table, enum_dim = %d\n", nt, enum_dim);
					#else
					fprintf(stderr, "enum_hamm_interval multi-thread (%d-thread). using table, enum_dim = %d\n", nt, enum_dim);
					#endif
					first = 0;
					use_system("VmSize");
				}

				static sub_dimension *table_spp = NULL;
				int spp_bit = SPP_BIT; // 追加のビット数
				static int table_spp_size = 0;
				if(table_spp == NULL) {
					table_spp_size = (1 << spp_bit) + 1;
					table_spp = make_table_for_enumeration_hamming(table_spp_size, spp_bit);
					fprintf(stderr, "make table for spplementary enumeation: spp_bit = %d, table_spp_size = %d\n", spp_bit, table_spp_size);
					fprintf(stderr, "made table for hamming enumeration: spp_bit = %d\n", spp_bit);
					use_system("VmSize");
				}

				int *bd_idx = qs->idx;
				int *bkt = bucket->bkt;

				// INTERVAL_WITH_PRIORITYが定義されているときは，interval_list をバッファとして用いることができる．
				interval *buff_pool = ivl->list, *buff;
				int buff_size = ivl->size;

				// 各スレッドで列挙するスケッチ数（空も含む）を同じにする．平均の10分の1くらいにする．
				// 実際にスレッドで求めたデータ数の合計を積算して，目標のデータ数を超えるまで繰り返す．
				// double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
				double ave_num = (double)bucket->num_data / (1 << PJT_DIM); // バケット（空も含む）の平均要素数
				int num_sketch = num_candidates * FACTOR_INF / ave_num; // 求めるスケッチの個数の期待値（FACTOR_INFで少し多めにする）
				int num_sketch_thread = num_sketch / nt / FACTOR_INF2; // 各スレッドが一回の操作で求めるスケッチの個数（目標の10分の1程度にしておく）
				int num_enum_data[nt], total_enum_data = 0;
				int total_enum_sketch = 0;
				// num_candidates  = 求めるデータ数（全体）
				// table[] = 部分集合を要素数の昇順に列挙するパターンを格納した配列
				// bd_idx[] = 相対位置でビット操作をするための距離下限の順位表
				// bkt[] = バケット表（bkt[s] = データをスケッチ順に並べたときに，スケッチ S のデータの先頭位置）
				// nt = スレッド数
				int t;										// スレッド番号
				int m, pn[nt];  							// スレッドが列挙したスケッチ数（パターン番号）
				int num_enum_sketches;						// スレッドが実際に列挙したスケッチ数（空も含む）
				int ne, num_nonempty[nt];							// スレッドが実際に列挙した空でないスケッチ数
				// int k[nt];  								// スレッドが求めたデータ番号数
				int table_size_thread = table_size / nt;	// スレッドが用いる列挙する部分集合のパターン配列の大きさ
				sub_dimension *q;							// スレッドが用いる列挙する部分集合のパターン配列
				sketch_type base_mask, base_mask2;
				sketch_type bm[nt], bm2[nt];					// 追加パターンで作成する mask
				int n_spp, n_spp2, ps[nt], ps2[nt];			// 追加で使用するパターン番号
				sketch_type sk, mask;
				for(t = 0; t < nt; t++) {
					pn[t] = 0;
					num_enum_data[t] = 0;
					num_nonempty[t] = 0;
					// k[t] = 0;
					bm[t] = 0;
					bm2[t] = 0;
					ps[t] = ps2[t] = 1;
				}
				do {
					#pragma omp parallel private(t, m, q, buff, base_mask, base_mask2, n_spp, n_spp2, sk, mask, num_enum_sketches)
					{
						t = omp_get_thread_num(); // スレッド番号の取得
						int en = 0;
						ne = num_nonempty[t];
						q = table + t * table_size_thread;
						buff = buff_pool + t * buff_size;
						m = pn[t];
						base_mask = bm[t]; base_mask2 = bm2[t];
						n_spp = ps[t]; n_spp2 = ps2[t];
						for(num_enum_sketches = 0; num_enum_sketches < num_sketch_thread; m++) // ここでは，データ数の制限はしない
						{
							mask = base_mask;
							for(int j = 0; j < q[m].num; j++) {
								mask |= (1 << bd_idx[(int)q[m].dim[j]]);
							}
							sk = qs->sketch ^ mask;
							if(bkt[sk + 1] > bkt[sk]) {// sk が空でなければ，buffに追加．ここでは，データ番号には展開しない
								buff[ne++] = (interval){priority(sk, qs), bkt[sk], bkt[sk + 1] - bkt[sk]};
								en += bkt[sk + 1] - bkt[sk];
							}
							num_enum_sketches++;
							if(q[m + 1].num == 0) { // 用意したパターンがなくなった．
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
						// quick_selection をまとめて並列処理する場合は，いったん，すべてのスレッドでの列挙が終わってから，
						// つまり，この並列ループが終わってから，別の並列ループで行う（並列処理の関数を呼び出す）
						num_nonempty[t] = ne;
						ivl->lg[t] = num_nonempty[t];
						pn[t] = m;
						bm[t] = base_mask; bm2[t] = base_mask2;
						ps[t] = n_spp; ps2[t] = n_spp2;
						num_enum_data[t] = en;
					}
					for(t = 0; t < nt; t++) {
						total_enum_data += num_enum_data[t];
					}
					total_enum_sketch += num_sketch;
					int num_rest = num_candidates * FACTOR_INF - total_enum_data; // 残りの候補数
					if(num_rest > 0) {
						#ifdef USE_AVE_ALL
						#else
						ave_num = (double)total_enum_data / total_enum_sketch;
						#endif
						num_sketch = num_rest / ave_num;
						if(num_sketch < nt) num_sketch = nt;
						num_sketch_thread = num_sketch / nt; 
					}
				} while (total_enum_data < num_candidates * FACTOR_INF);
				int nc2 = quick_select_sum_k_interval_by_multi_thread(ivl, num_candidates);
				if(nc2 == 0) {
					fprintf(stderr, "q = %d, num_candidates = %d, nc2 = %d\n", qs->query.query_num, num_candidates, nc2);
					getchar();
				}
				return nc2;
			}
			#else // USE_MU
#ifdef USE_MU_COMMON
			int filtering_by_sketch_enumeration_hamming_interval(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates)
			{
				int n = PARA_ENUM_INF;
				int nt = (1 << n); // スレッド数
				#ifdef _OPENMP
				omp_set_num_threads(nt);
				#endif
				int *bd_idx = qs->idx;
				int *bkt = bucket->bkt;

				#ifndef THREAD_PLUS
				int mu_thread_size = nt;
				#else
				int mu_thread_size = nt * (1 << THREAD_PLUS);
				#endif

				sketch_type mu_thread[mu_thread_size]; // スレッドに割り当てられた下位 n-bit のパターン（or 質問スケッチとのXOR）
				mu_thread[0] = 0;			// このパターンの順序は，おそらく，ほとんど recall に影響しないので，INF 用の grey code 生成ルールを用いて準備する
				for(int t = 0; t < mu_thread_size; t++) {
					mu_thread[t + 1] = mu_thread[t] ^ (1 << bd_idx[bit_count(t ^ (t + 1)) - 1]);
				}
				// #define MU_SORT
				#ifdef MU_SORT
					int idx_mu[mu_thread_size];
					dist_type pr_thread[mu_thread_size];
					for(int t = 0; t < mu_thread_size; t++) {
						pr_thread[t] = priority(qs->sketch ^ mu_thread[t], qs);
						idx_mu[t] = t;
					}
					quick_sort(idx_mu, pr_thread, 0, mu_thread_size - 1);
					sketch_type mu_temp[mu_thread_size];
					dist_type pr_temp[mu_thread_size];
					for(int t = 0; t < mu_thread_size; t++) {
						mu_temp[t] = mu_thread[idx_mu[t]];
						pr_temp[t] = pr_thread[idx_mu[t]];
					}
					for(int t = 0; t < mu_thread_size; t++) {
						mu_thread[t] = mu_temp[t];
						pr_thread[t] = pr_temp[t];
				//		printf("pr[%d] = %d\n", t, pr_thread[t]);
					}
				//	getchar();
				#endif

				// ハミング距離順の列挙のための部分集合の表を用意する．
				static sub_dimension *table_low = NULL;
				static int table_low_size = 0;
				#ifdef ENUM_DIM
				int enum_dim = ENUM_DIM;
				#else
				int enum_dim = PJT_DIM - 22;	// ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
				#endif
				#ifndef THREAD_PLUS
				int low_dim = enum_dim - n;		// スレッドで固定する割り当て分の n-bit を減らす
				#else
				int low_dim = enum_dim - n - THREAD_PLUS;		// スレッドで固定する割り当て分の n-bit を減らす
				#endif
				if(table_low == NULL) {
					table_low_size = (1 << low_dim) + nt;
					table_low = make_table_for_enumeration_hamming(table_low_size, low_dim);
					fprintf(stderr, "made table for hamming enumeration: low_dim = %d\n", low_dim);
					use_system("VmSize");
				}

				int add_dim = SPP_BIT;
				static int first = 1;
				if(first) {
					fprintf(stderr, "(*mu_common + post-selection) enum_hamm_interval multi-thread (%d-thread). low-add = %d-%d\n", nt, low_dim, add_dim);
					first = 0;
					use_system("VmSize");
				}

				static sub_dimension *table_add = NULL;
				static int table_add_size = 0;
				if(table_add == NULL) {
					table_add_size = (1 << add_dim) + 1;
					table_add = make_table_for_enumeration_hamming(table_add_size, add_dim);
					fprintf(stderr, "make table for spplementary enumeation: add_dim = %d\n", add_dim);
					use_system("VmSize");
				}
				
				#ifdef DECAF
				#define MAX_PRIORITY 1200
				#elif defined(DEEP1B)
				#define MAX_PRIORITY 500
				#else
				#define MAX_PRIORITY 200
				#endif

				#ifdef THREAD_PLUS
				static vlist *vl[(1 << PARA_ENUM_INF) * (1 << THREAD_PLUS)][MAX_PRIORITY + 1];
				#else
				static vlist *vl[(1 << PARA_ENUM_INF)][MAX_PRIORITY + 1];
				#endif
//					#define PR_FACTOR 20
//					dist_type max_priority = MAX_PRIORITY / PR_FACTOR;
				dist_type max_priority = MAX_PRIORITY;
//					for(int i = 0; i < low_dim + add_dim; i++) {
//						max_priority += qs->bd[qs->idx[i]];
//					}
//					dist_type max_priority = priority(-1 ^ qs->sketch, qs);
/*					static int max_priority_all = 0;
				if(qs->query.query_num == 0) {
					max_priority_all = max_priority;
					fprintf(stderr, "q = %d, max_priority = %d\n", qs->query.query_num, max_priority);
				} else if(max_priority > max_priority_all) {
					max_priority_all = max_priority;
					fprintf(stderr, "q = %d, max_priority = %d\n", qs->query.query_num, max_priority);
				}
				if(max_priority > MAX_PRIORITY) {
					fprintf(stderr, "Too large max_priority (= %d), q = %d\n", max_priority, qs->query.query_num);
					exit(1);
				}
*/					
				int num_data_of_priority[mu_thread_size][max_priority + 1];
				int lg_of_priority[mu_thread_size][max_priority + 1];
				int total_enum_data[max_priority + 1];

// fprintf(stderr, "max_priority = %d, max_low_priority = %d\n", max_priority, max_low_priority); getchar();
//					if(qs->query.query_num == 0) {
//						fprintf(stderr, "q = %d, max_priority = %d\n", qs->query.query_num, max_priority);
//					}
				#ifndef VLIST_SIZE
				#define VLIST_SIZE 10
				#define VLIST_STEP 10
				#endif

				static int initialized_vlist = 0;
				if(!initialized_vlist) {
					for(int t = 0; t < mu_thread_size; t++) {
						for(dist_type pr = 0; pr <= max_priority; pr++) {
							vl[t][pr] = new_vlist(VLIST_SIZE, VLIST_STEP);
						}
					}
					initialized_vlist = 1;
				} else {
					for(int t = 0; t < mu_thread_size; t++) {
						for(dist_type pr = 0; pr <= max_priority; pr++) {
							if(vl[t][pr] == NULL) {
								vl[t][pr] = new_vlist(VLIST_SIZE, VLIST_STEP);
							} else {
								makenull_vlist(vl[t][pr]);
							}
						}
					}
				}
				
				for(int t = 0; t < mu_thread_size; t++) {
					for(dist_type pr = 0; pr <= max_priority; pr++) {
						if(vl[t][pr]->num_list != 0 || vl[t][pr]->num_data != 0) {
							fprintf(stderr, "q = %d, invalid initialized vlist: num_list = %d, num_data = %d\n", qs->query.query_num, vl[t][pr]->num_list, vl[t][pr]->num_data);
							getchar();
						}
					}
				}
				int num_enum_data[nt];		// スレッドが実際に列挙した空でないスケッチ数とデータ数
				for(int t = 0; t < nt; t++) {
					num_enum_data[t] = 0;
				}
				int mask_id = 0;								// 列挙する共通マスクのid（先頭を0とした連番）
				int loop = FACTOR_INF2;									// まとめて共通マスク mu_mask を作成する個数
				sketch_type mu_common[loop];
				int total_data;
				#ifdef THREAD_PLUS
				int nn = n + THREAD_PLUS;
				#else
				int nn = n;
				#endif
				do {
//						#pragma omp parallel for
					for(int i = mask_id; i < mask_id + loop; i++) {
						sketch_type mu_low = Mu(nn, low_dim, bd_idx, i % (1 << low_dim), table_low);
						sketch_type mu_add = Mu(nn + low_dim, add_dim, bd_idx, i / (1 << low_dim), table_add);
						mu_common[i - mask_id] = mu_low ^ mu_add;
					}
					mask_id += loop;
					#pragma omp parallel 
//						for(int t = 0; t < nt; t++)
					{
						int t = omp_get_thread_num();	// スレッド番号
						#ifdef THREAD_PLUS
						for(int tn = t * (1 << THREAD_PLUS); tn < (t + 1) * (1 << THREAD_PLUS); tn++) 
						#else
						int tn = t;
						#endif
						{
							sketch_type sk0 = qs->sketch ^ mu_thread[tn];
							int nd = 0;						// 一回の処理で各スレッドで求めたデータ数
							for(int i = 0; i < loop; i++) {
								sketch_type sk = sk0 ^ mu_common[i];
								if(bkt[sk + 1] > bkt[sk]) {
//										dist_type pr = priority(sk, qs) / PR_FACTOR;
									dist_type pr = priority(sk, qs);
									add_vlist(vl[tn][pr], (interval){pr, bkt[sk], bkt[sk + 1] - bkt[sk]});
									nd += bkt[sk + 1] - bkt[sk];
								}
							}
							num_enum_data[t] += nd;
						}
					}
					total_data = 0;
					for(int t = 0; t < nt; t++) {
						total_data += num_enum_data[t];
					}
				} while (total_data < num_candidates * FACTOR_INF);
// printf("q, %d, max_priority, %d, => ,", qs->query.query_num, max_priority);
				int sd;
				while(1) {
					sd = 0;
					for(int t = 0; t < mu_thread_size; t++) {
						sd += vl[t][max_priority]->num_data;
					}
					if(sd != 0) break;
					max_priority--;
				}
//					static int max_priority_all = 0;
//					if(qs->query.query_num == 0) {
//						max_priority_all = max_priority;
//						fprintf(stderr, "q = %d, max_priority = %d\n", qs->query.query_num, max_priority);
//					} else if(max_priority > max_priority_all) {
//						max_priority_all = max_priority;
//						fprintf(stderr, "q = %d, max_priority = %d\n", qs->query.query_num, max_priority);
//					}
// printf("%d\n", max_priority); getchar();
				for(int pr = 0; pr <= max_priority; pr++) {
					total_enum_data[pr] = 0;
					for(int t = 0; t < mu_thread_size; t++) {
						num_data_of_priority[t][pr] = lg_of_priority[t][pr] = 0;
					}
				}

				// num_data_of_priority[t][pr] を pr = 0, ... , pr の累積に変更する．
				// lg_of_priority[t][pr] も累積
				for(int t = 0; t < mu_thread_size; t++) {
					num_data_of_priority[t][0] += vl[t][0]->num_data;
					lg_of_priority[t][0] += vl[t][0]->num_list;
					for(int pr = 1; pr <= max_priority; pr++) {
						num_data_of_priority[t][pr] += num_data_of_priority[t][pr - 1] + vl[t][pr]->num_data;
						lg_of_priority[t][pr] += lg_of_priority[t][pr - 1] + vl[t][pr]->num_list;
					}
				}
				// 全合計が目標のnum_candidatesを超える最小の優先度 final_priority を求める．
				int final_priority = 0;
				for(int pr = 0; pr <= max_priority; pr++) {
					total_enum_data[pr] = 0;
					for(int t = 0; t < mu_thread_size; t++) {
						total_enum_data[pr] += num_data_of_priority[t][pr];
					}
					if(total_enum_data[pr] < num_candidates) final_priority = pr;
				}
				final_priority++;

				// vl から ivl にまとめる
				#pragma omp parallel for
				for(int t = 0; t < mu_thread_size; t++) {
					interval *buff = ivl->list + t * ivl->size;
					int lg = 0;
					for(int pr = 0; pr <= final_priority; pr++) {
						if(vl[t][pr]->num_list == 0) continue;
						memcpy(buff + lg, vl[t][pr]->elm, sizeof(interval) * vl[t][pr]->num_list);
						lg += vl[t][pr]->num_list;
//							for(int i = 0; i < vl[t][pr]->num_list; i++) {
//								buff[lg++] = vl[t][pr]->elm[i];
//							}
					}
					ivl->lg[t] = lg;
				}
				// printf("total_enum_data = %d ===> ", total_enum_data[final_priority]);
				// final_priority の区間で余分なもの（データ数がnum_candidatesを超える部分）を削除
				for(int t = mu_thread_size - 1; t >= 0; t--) {
					interval *buff = ivl->list + t * ivl->size;
					if(total_enum_data[final_priority] - (num_data_of_priority[t][final_priority] - num_data_of_priority[t][final_priority - 1]) >= num_candidates) {
						ivl->lg[t] = lg_of_priority[t][final_priority - 1];
						total_enum_data[final_priority] -= (num_data_of_priority[t][final_priority] - num_data_of_priority[t][final_priority - 1]);
					} else {
						int j = lg_of_priority[t][final_priority] - 1;
						while(buff[j].priority >= final_priority && total_enum_data[final_priority] - buff[j].run >= num_candidates) {
							total_enum_data[final_priority] -= buff[j--].run;
						}
						ivl->lg[t] = j;
					}
				}
				return total_enum_data[final_priority];
			}
#else
			// スケッチ列挙(Hamming)によるフィルタリング（バケット（配列 idx と bkt）利用）（部分集合列挙の表を利用する）（マスク生成関数muを使用）（multi-threasd）
			int filtering_by_sketch_enumeration_hamming_interval(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates)
			{
				// fprintf(stderr, "(5) filtering_by_sketch_enumeration_hamming （multi-threasd）(Using Mu)\n"); exit(0);
				int n = PARA_ENUM_INF;
				int nt = (1 << n); // スレッド数
				#ifdef _OPENMP
				omp_set_num_threads(nt);
				#endif
				// ハミング距離順の列挙のための部分集合の表を用意する．
				static sub_dimension *table = NULL;
				static int table_size = 0;
				#ifdef ENUM_DIM
				int enum_dim = ENUM_DIM;
				#else
				int enum_dim = PJT_DIM - 22; // ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
				#endif
				if(table == NULL) {
					table_size = (1 << enum_dim) + nt;
					table = make_table_for_enumeration_hamming(table_size, enum_dim);
					//　rearrange_table(nt, table, table_size);　// mu では表の並べ替えをしない（遅くなるようなら後で修正）
					fprintf(stderr, "made table for hamming enumeration: enum_dim = %d\n", enum_dim);
					use_system("VmSize");
				}

				static int first = 1;
				if(first) {
					#ifdef SELECT_SUM
					fprintf(stderr, "(2*) ");
					#else
					fprintf(stderr, "(1*+mu) ");
					#endif
					#ifdef USE_DIFF_TABLE
					fprintf(stderr, "enum_hamm_interval multi-thread (%d-thread). using diff table, enum_dim = %d\n", nt, enum_dim);
					#else
					fprintf(stderr, "enum_hamm_interval multi-thread (%d-thread). using table, enum_dim = %d\n", nt, enum_dim);
					#endif
					first = 0;
					use_system("VmSize");
				}

				static sub_dimension *table_spp = NULL;
				int spp_bit = SPP_BIT; // 追加のビット数
				static int table_spp_size = 0;
				if(table_spp == NULL) {
					table_spp_size = (1 << spp_bit) + 1;
					table_spp = make_table_for_enumeration_hamming(table_spp_size, spp_bit);
					fprintf(stderr, "make table for spplementary enumeation: spp_bit = %d, table_spp_size = %d\n", spp_bit, table_spp_size);
					fprintf(stderr, "made table for hamming enumeration: spp_bit = %d\n", spp_bit);
					use_system("VmSize");
				}

				int *bd_idx = qs->idx;
				int *bkt = bucket->bkt;

				// INTERVAL_WITH_PRIORITYが定義されているときは，interval_list をバッファとして用いることができる．
				interval *buff_pool = ivl->list, *buff;
				// 各スレッドで列挙するスケッチ数（空も含む）を同じにする．平均の2分の1くらい（FACTOR_INF2）にする．
				// 実際にスレッドで求めたデータ数の合計を積算して，目標のデータ数を超えるまで繰り返す．
				double ave_num = (double)bucket->num_data / (1 << PJT_DIM); // バケット（空も含む）の平均要素数
//					typedef struct {
//						int num_enum_data; 						// スレッドで求めたデータ数
//						int num_nonempty;						// スレッドげ求めた空でないスケッチ数
//						sketch_type mu_add;						// マスクパターン
//						int add_p;								// add-bit カウンタの直前の値
//					} struct_work_enum;
//					struct_work_enum work[nt];
				int num_enum_data[nt];
				int total_enum_data = 0;
//					int total_enum_sketch = 0;
				// num_candidates  = 求めるデータ数（全体）
				// table[] = 部分集合を要素数の昇順に列挙するパターンを格納した配列
				// bd_idx[] = 相対位置でビット操作をするための距離下限の順位表
				// bkt[] = バケット表（bkt[s] = データをスケッチ順に並べたときに，スケッチ S のデータの先頭位置）
				// nt = スレッド数
				int t;										// スレッド番号
				int sj = 0;  									// スレッドが列挙したスレッド内でのスケッチ番号
				int num_sketch = num_candidates * FACTOR_INF / ave_num; // 求めるスケッチの個数の期待値（FACTOR_INFで少し多めにする）
				int ns = num_sketch / nt / FACTOR_INF2; 	// 各スレッドが一回の操作で求めるスケッチの個数（目標の2（FACTOR_INF2）分の1程度にしておく）
				int num_nonempty[nt];						// スレッドが実際に列挙した空でないスケッチ数
				sketch_type mu_low, mu_add[nt];
				int add_c, add_p[nt];						// add-bit のカウンタ，直前のカウンタ
				sketch_type sk;
				for(t = 0; t < nt; t++) {
					num_enum_data[t] = 0;
					num_nonempty[t] = 0;
					add_p[t] = -1;
					mu_add[t] = 0;
				}
				int i, j;
				int ne, nd, ap;
				sketch_type md;
				do {
					#pragma omp parallel private(t, ne, nd, ap, md, add_c, mu_low, buff, i, j, sk)
					{
						t = omp_get_thread_num(); // スレッド番号の取得
						nd = num_enum_data[t];
						ne = num_nonempty[t];
						ap = add_p[t];
						md = mu_add[t];
						buff = buff_pool + t * ivl->size;
						for(j = sj;  j < sj + ns; j++) {
							i = j * nt + t;
							mu_low = Mu(0, enum_dim, bd_idx, i % (1 << enum_dim), table);
							add_c = i / (1 << enum_dim);
							if(add_c != ap) {
								md = Mu(enum_dim, spp_bit, bd_idx, add_c, table_spp);
								ap = add_c;
							}
							sk = qs->sketch ^ mu_low ^ md;
							// sk が空でなければ，buffに追加．ここでは，データ番号には展開しない
							if(bkt[sk + 1] > bkt[sk]) {
								buff[ne++] = (interval){priority(sk, qs), bkt[sk], bkt[sk + 1] - bkt[sk]};
								nd += bkt[sk + 1] - bkt[sk];
							}
						}
						// quick_selection をまとめて並列処理する場合は，いったん，すべてのスレッドでの列挙が終わってから，
						// つまり，この並列ループが終わってから，別の並列ループで行う（並列処理の関数を呼び出す）
						ivl->lg[t] = num_nonempty[t] = ne;
						num_enum_data[t] = nd;
						add_p[t] = ap;
						mu_add[t] = md;
					}
					total_enum_data = 0;
					for(t = 0; t < nt; t++) {
						total_enum_data += num_enum_data[t];
					}
					sj += ns;
					int num_rest = num_candidates * FACTOR_INF - total_enum_data; // 残りの候補数
					if(num_rest > 0) {
						#ifdef USE_AVE_ALL
						#else
						ave_num = (double) total_enum_data / (sj * nt);
						#endif
						num_sketch = num_rest / ave_num;
						if(num_sketch < nt) num_sketch = nt;
						ns = num_sketch / nt; 
					}
				} while (total_enum_data < num_candidates * FACTOR_INF);
				int nc2 = quick_select_sum_k_interval_by_multi_thread(ivl, num_candidates);
				if(nc2 == 0) {
					fprintf(stderr, "q = %d, num_candidates = %d, nc2 = %d\n", qs->query.query_num, num_candidates, nc2);
					getchar();
				}
				return nc2;
			}
#endif
			int filtering_by_sketch_enumeration_hamming_interval_3(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates)
			{
				int n = PARA_ENUM_INF;
				int nt = (1 << n); // スレッド数
				#ifdef _OPENMP
				omp_set_num_threads(nt);
				#endif
				int *bd_idx = qs->idx;
				int *bkt = bucket->bkt;

				sketch_type mu_thread[nt]; // スレッドに割り当てられた下位 n-bit のパターン（or 質問スケッチとのXOR）
				// このパターンの順序は，おそらく，ほとんど recall に影響しないので，INF 用の grey code 生成ルールを用いて準備する
				mu_thread[0] = 0;
				for(int t = 0; t < nt; t++) {
					mu_thread[t + 1] = mu_thread[t] ^ (1 << bd_idx[bit_count(t ^ (t + 1)) - 1]);
				}
//					double ave_sample = (double)sample_nd / nt;

				// ハミング距離順の列挙のための部分集合の表を用意する．
				static sub_dimension *table_low = NULL;
				static int table_low_size = 0;
				#ifdef ENUM_DIM
				int enum_dim = ENUM_DIM;
				#else
				int enum_dim = PJT_DIM - 22; // ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
				#endif
				int low_dim = enum_dim - n; // スレッドで固定する割り当て分の n-bit を減らす
				if(table_low == NULL) {
					table_low_size = (1 << low_dim) + nt;
					table_low = make_table_for_enumeration_hamming(table_low_size, low_dim);
					fprintf(stderr, "made table for hamming enumeration: low_dim = %d\n", low_dim);
					use_system("VmSize");
				}

				int add_dim = SPP_BIT;
				static int first = 1;
				if(first) {
					fprintf(stderr, "(mu**3+ps) enum_hamm_interval multi-thread (%d-thread). low-add = %d-%d\n", nt, low_dim, add_dim);
					first = 0;
					use_system("VmSize");
				}

				static sub_dimension *table_add = NULL;
				static int table_add_size = 0;
				if(table_add == NULL) {
					table_add_size = (1 << add_dim) + 1;
					table_add = make_table_for_enumeration_hamming(table_add_size, add_dim);
					fprintf(stderr, "make table for spplementary enumeation: add_dim = %d\n", add_dim);
					use_system("VmSize");
				}

				int num_enum_data[nt], num_nonempty[nt];		// スレッドが実際に列挙した空でないスケッチ数とデータ数
				for(int t = 0; t < nt; t++) {
					num_enum_data[t] = num_nonempty[t] = 0;
				}
				int mask_id = 0;								// 列挙する共通マスクのid（先頭を0とした連番）
				int total_enum_data = 0;						// 求めたデータ数の合計
				#define MAX_LOOP 4000
				#define MIN_LOOP 200
				sketch_type mu_common[MAX_LOOP];
				int loop = 300;
				do {
//						#pragma omp parallel for
					for(int i = mask_id; i < mask_id + loop; i++) {
						sketch_type mu_low = Mu(n, low_dim, bd_idx, i % (1 << low_dim), table_low);
						sketch_type mu_add = Mu(n + low_dim, add_dim, bd_idx, i / (1 << low_dim), table_add);
						mu_common[i - mask_id] = mu_low ^ mu_add;
					}
					mask_id += loop;
					#pragma omp parallel 
					{
						int t = omp_get_thread_num();	// スレッド番号
						int ne = num_nonempty[t];		// 各スレッドで求めた空でないスケッチ数
						int nd = 0;						// 一回の処理で各スレッドで求めたデータ数
						interval *buff = ivl->list + t * ivl->size;
						for(int i = 0; i < loop; i++) {
							sketch_type sk = qs->sketch ^ mu_common[i] ^ mu_thread[t];
							if(bkt[sk + 1] > bkt[sk]) {
								buff[ne++] = (interval){priority(sk, qs), bkt[sk], bkt[sk + 1] - bkt[sk]};
								nd += bkt[sk + 1] - bkt[sk];
							}
						}
						num_nonempty[t] = ne;
						num_enum_data[t] += nd;
					}
					total_enum_data = 0;
					for(int t = 0; t < nt; t++) {
						ivl->lg[t] = num_nonempty[t];
						total_enum_data += num_enum_data[t];
					}
				} while (total_enum_data < num_candidates * FACTOR_INF);

				int nc2 = quick_select_sum_k_interval_by_multi_thread(ivl, num_candidates);
				if(nc2 == 0) {
					fprintf(stderr, "q = %d, num_candidates = %d, nc2 = %d\n", qs->query.query_num, num_candidates, nc2);
					getchar();
				}
				return nc2;
			}
			#endif // USE_MU
			#else // SELECT_BY_PARA_MERGE　（上位のスケッチ選択を並列処理で行った後でマージ処理する）
			// スケッチ列挙(Hamming)によるフィルタリング（バケット（配列 idx と bkt）利用）（部分集合列挙の表を利用する）（multi-threasd）
			int filtering_by_sketch_enumeration_hamming(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates)
			{
				int n = PARA_ENUM_INF;
				int nt = (1 << n); // スレッド数
				#ifdef _OPENMP
				omp_set_num_threads(nt);
				#endif

				// ハミング距離順の列挙のための部分集合の表を用意する．
				static sub_dimension *table = NULL;
				static int table_size = 0;
				#ifdef ENUM_DIM
				int enum_dim = ENUM_DIM;
				#else
				int enum_dim = PJT_DIM - 22; // ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
				#endif
				if(table == NULL) {
					table_size = (1 << enum_dim) + nt;
					table = make_table_for_enumeration_hamming(table_size, enum_dim);
					rearrange_table(nt, table, table_size);
				}

				static int first = 1;
				if(first) {
					#ifdef USE_DIFF_TABLE
					fprintf(stderr, "(4) enum_hamm multi-thread (%d-thread). using diff table, enum_dim = %d\n", nt, enum_dim);
					#else
					fprintf(stderr, "(4) enum_hamm multi-thread (%d-thread). using table, enum_dim = %d\n", nt, enum_dim);
					#endif
					first = 0;
				}

				#define TABLE_SPP
				// tableを用いた列挙では不足するときに，追加のビットパターンを求めるための表
				#ifdef TABLE_SPP
				static sub_dimension *table_spp = NULL;
				#ifdef SPP_BIT
				int spp_bit = SPP_BIT; // 追加のビット数
				#else
				int spp_bit = 17; // 追加のビット数
				#endif
				static int table_spp_size = 0;
				if(table_spp == NULL) {
					table_spp_size = (1 << spp_bit) + 1;
					table_spp = make_table_for_enumeration_hamming(table_spp_size, spp_bit);
					fprintf(stderr, "make table for spplementary enumeation: spp_bit = %d, table_spp_size = %d\n", spp_bit, table_spp_size);
					#ifdef USE_DIFF_TABLE
					for(int i = 0; i < table_spp_size - 1; i++) {
						diff(&table_spp[i], &table_spp[i + 1], &table_spp[i]);
						if(table_spp[i + 1].num == 0) break; // 次が空集合になったら終了
					}
					#endif
				}
				#endif

				#ifndef WITHOUT_IDX
				int *bd_idx = qs->idx;
				#endif
				int *bkt = bucket->bkt;

				// num_candidates  = 求めるデータ数（全体）
				// data_num[] = 求めたデータ番号を格納する配列
				// table[] = 部分集合を要素数の昇順に列挙するパターンを格納した配列
				// bd_idx[] = 相対位置でビット操作をするための距離下限の順位表
				// bkt[] = バケット表（bkt[s] = データをスケッチ順に並べたときに，スケッチ S のデータの先頭位置）
				// nt = スレッド数
				int num_data_thread = num_candidates / nt;	// スレッド求めるデータ番号数
				int t;										// スレッド番号
				int m;  									// スレッドが列挙したスケッチ数（パターン番号）
				int num_nonempty;							// スレッドが実際に列挙した空でないスケッチ数
				int k;  									// スレッドが求めたデータ番号数
				int table_size_thread = table_size / nt;	// スレッドが用いる列挙する部分集合のパターン配列の大きさ
				sub_dimension *q;							// スレッドが用いる列挙する部分集合のパターン配列
				int *d;										// スレッドが求めたデータ番号を格納する配列
				sketch_type base_mask = 0, base_mask2 = 0;	// 追加パターンで作成する mask
				int n_spp, n_spp2;							// 追加で使用するパターン番号
				sketch_type sk, mask;

				double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
				int num_sketch = num_candidates * FACTOR_INF / ave_num; // 求めるスケッチの個数（FACTOR_INFで少し多めにする）
				int num_sketch_thread = num_sketch / nt; // 各スレッドが求めるスケッチの個数（注意：このスケッチ数では，データ数が num_data_thread に達しない可能性がある）
				static int buff_pool_size = 0;
				int num_enumerated_data_thread[nt]; // スレッドで列挙したスケッチを持つデータ数
				int num_selected_sketches[nt]; // 列挙後に quick_select で選択した priority 上位のスケッチ数
				int num_selected_data[nt]; // 上のスケッチをもつデータ数
				#ifdef SELECT_SUM
					static sketch_with_priority_num *buff_pool = NULL;
					if(buff_pool == NULL || buff_pool_size < num_sketch_thread * nt) {
						if(buff_pool != NULL) {
							FREE(buff_pool, sizeof(sketch_with_priority_num) * buff_pool_size);
						}
						buff_pool_size = num_sketch_thread * nt;
						fprintf(stderr, "malloc answer buffer, ");
						buff_pool = MALLOC(sizeof(sketch_with_priority_num) * buff_pool_size);
						fprintf(stderr, "OK: buff_pool_size = %d, num_thread = %d\n", buff_pool_size, nt);
					}
					sketch_with_priority_num *buff;
				#else
					static answer_type *buff_pool = NULL;
					if(buff_pool == NULL || buff_pool_size < num_sketch_thread * nt) {
						if(buff_pool != NULL) {
							FREE(buff_pool, sizeof(answer_type) * buff_pool_size);
						}
						buff_pool_size = num_sketch_thread * nt;
						fprintf(stderr, "malloc answer buffer, ");
						buff_pool = MALLOC(sizeof(answer_type) * buff_pool_size);
						fprintf(stderr, "OK: buff_pool_size = %d, num_thread = %d\n", buff_pool_size, nt);
					}
					answer_type *buff;
				#endif

				#pragma omp parallel private(t, m, num_nonempty, k, q, d, base_mask, base_mask2, n_spp, n_spp2, buff, sk, mask)
				{
					t = omp_get_thread_num(); // スレッド番号の取得
					m = 0;
					num_nonempty = 0;
					k = 0;
					q = table + t * table_size_thread;
					d = data_num + t * num_data_thread;
					base_mask = base_mask2 = 0;
					n_spp = n_spp2 = 1;
					buff = buff_pool + t * num_sketch_thread;

					for(m = 0; num_nonempty < num_sketch_thread /* && k < num_data_thread * 3 // ここでは，データ数の制限はしない */; m++) {
						mask = base_mask;
						for(int j = 0; j < q[m].num; j++) {
							#ifndef WITHOUT_IDX
							mask |= (1 << bd_idx[(int)q[m].dim[j]]);
							#else
							mask |= (1 << (int)q[m].dim[j]);
							#endif
						}
						sk = qs->sketch ^ mask;
						if(bkt[sk + 1] > bkt[sk]) { // sk が空でなければ，buffに追加．ここでは，データ番号には展開しない
							#ifdef SELECT_SUM
							buff[num_nonempty++] = (sketch_with_priority_num){sk, priority(sk, qs), bkt[sk + 1] - bkt[sk]};
							#else
							buff[num_nonempty++] = (answer_type){sk, priority(sk, qs)};
							#endif
							k += bkt[sk + 1] - bkt[sk];
						}
						if(q[m + 1].num == 0) { // 用意したパターンがなくなった．
							#ifdef TABLE_SPP
								if(table_spp[n_spp].num == 0) { // 追加分のパターンも使い切った．
									if(table_spp[n_spp2].num == 0) break; // 2回目の追加分のパターンも使い切った．
									n_spp = 0;
									n_spp2++;
									base_mask2 = 0;
									for(int m = 0; m < table_spp[n_spp2].num; m++) {
										#ifndef WITHOUT_IDX
										base_mask2 |= (1 << bd_idx[(int)table_spp[n_spp2].dim[m] + enum_dim + spp_bit]);
										#else
										base_mask2 |= (1 << ((int)table_spp[n_spp2].dim[m] + enum_dim + spp_bit));
										#endif
									}
								}
								base_mask = base_mask2;
								for(int m = 0; m < table_spp[n_spp].num; m++) {
									#ifndef WITHOUT_IDX
									base_mask |= (1 << bd_idx[(int)table_spp[n_spp].dim[m] + enum_dim]);
									#else
									base_mask |= (1 << ((int)table_spp[n_spp].dim[m] + enum_dim));
									#endif
								}
								m = -1; // forの再初期化で m++ となって，0になる．
								n_spp++;
							#else
								while(k < num_data_thread) {
									d[k++] = 0;
								}
							#endif
						}
					}
					num_enumerated_data_thread[t] = k;
					#ifdef SELECT_SUM
						if(k > num_data_thread * FACTOR_INF) { // 列挙したスケッチを持つデータ数が十分多いときは，priority 上位のものを選択する．
							num_selected_sketches[t] = quick_select_sum_k_sketch_with_priority_num(buff, 0, num_nonempty - 1, num_data_thread * FACTOR_INF);
							num_selected_data[t] = num_data_thread * FACTOR_INF;
						} else {
							num_selected_sketches[t] = num_nonempty;
							num_selected_data[t] = num_enumerated_data_thread[t];
						}
					#else
						int sel = (num_sketch / nt < num_nonempty ? num_sketch / nt : num_nonempty);
						quick_select_k_answer(buff, 0, num_nonempty - 1, sel * 0.4);
					#endif
				}

				// スレッドでSELECTしたスケッチを詰め合わせ，それからSELECTする．
				int total_num_selected_sketches = num_selected_sketches[0];
				int total_num_data_enumerated = num_selected_data[0];
				for(int t = 1; t < nt; t++) {
					memcpy(buff_pool + total_num_selected_sketches, buff_pool + t * num_sketch_thread, sizeof(sketch_with_priority_num) * num_selected_sketches[t]);
					total_num_selected_sketches += num_selected_sketches[t];
					total_num_data_enumerated += num_selected_data[t];
				}
				if(total_num_data_enumerated < num_candidates * FACTOR_INF) { // 列挙したスケッチのデータ総数が少ないときは，目標の候補数を減らす．
					num_candidates = total_num_data_enumerated / FACTOR_INF;
				}

				int nk = quick_select_sum_k_sketch_with_priority_num(buff_pool, 0, total_num_selected_sketches - 1, num_candidates);

				#if defined(_OPENMP) && NUM_THREADS > 1
					omp_set_num_threads(nk < NUM_THREADS ? nk : NUM_THREADS);
					nt = omp_get_max_threads(); 	// スレッド数を求める
					int sk_th[nt + 1]; // スレッドが列挙するスケッチの先頭位置．スレッド t が列挙するスケッチは，sk_th[t] から sk_th[t + 1] - 1 まで． 
					int num[nt]; // スレッド 0 から スレッド t - 1 が列挙するスケッチのデータ数の合計．スレッド t が列挙したスケッチのデータは，data_num[num[t]] から data_num[num[t + 1] - 1] まで．
					int p = 0, sum = 0;
					for(int t = 0; t < nt; t++) {
						while(sum < (num_candidates / nt) * t) {
							sum += buff_pool[p++].num;
						}
						sk_th[t] = p; 
						num[t] = sum;
					}
					sk_th[nt] = nk;
					#pragma omp parallel
					{
						int t = omp_get_thread_num(); // スレッド番号の取得
						int k = 0;
						int *d = data_num + num[t];
						for(int p = sk_th[t]; p < sk_th[t + 1]; p++) {
							sketch_type sk = (sketch_type)(buff_pool[p].sk);
							for(int j = bkt[sk]; j < bkt[sk + 1]; j++, k++) {
								#ifdef DATA_NUM_IN_SKETCH_ORDER
								d[k] = j;
								#else
								d[k] = bucket->idx[j];
								#endif
							}
						}
					}
				#else
					k = 0;
					for(int i = 0; i < num_enumerated_sketches; i++) {
						s = buff_pool[i].sk;
						for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates; j++, k++) {
							#ifdef DATA_NUM_IN_SKETCH_ORDER
							data_num[k] = j;
							#else
							data_num[k] = bucket->idx[j];
							#endif
						}
					}
				#endif

				return num_candidates;

			}
			#endif // SELECT_BY_PARA_MERGE
		#else // !SELECT_BY_SINGLE 
		// スケッチ列挙(Hamming)によるフィルタリング（バケット（配列 idx と bkt）利用）（部分集合列挙の表を利用する）（multi-threasd）
		// スケッチのD1上位選択をsingle-threadで行う．
		int filtering_by_sketch_enumeration_hamming(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates)
		{
			int n = PARA_ENUM_INF;
			int nt = (1 << n); // スレッド数
			#ifdef _OPENMP
			omp_set_num_threads(nt);
			#endif

			// ハミング距離順の列挙のための部分集合の表を用意する．
			static sub_dimension *table = NULL;
			static int table_size = 0;
			#ifdef ENUM_DIM
			int enum_dim = ENUM_DIM;
			#else
			int enum_dim = PJT_DIM - 22; // ハミング距離順の列挙を求めるための次元数（射影次元より小さくする）
			#endif
			if(table == NULL) {
				table_size = (1 << enum_dim) + nt;
				table = make_table_for_enumeration_hamming(table_size, enum_dim);
				rearrange_table(nt, table, table_size);
			fprintf(stderr, "made table for hamming enumeration: enum_dim = %d\n", enum_dim);
			use_system("VmSize");
			}

			static int first = 1;
			if(first) {
				#ifdef USE_DIFF_TABLE
				fprintf(stderr, "(3) enum_hamm multi-thread (%d-thread). using diff table, enum_dim = %d\n", nt, enum_dim);
				#else
				fprintf(stderr, "(3) enum_hamm multi-thread (%d-thread). using table, enum_dim = %d\n", nt, enum_dim);
				#endif
				first = 0;
			use_system("VmSize");
			}

			#define TABLE_SPP
			// tableを用いた列挙では不足するときに，追加のビットパターンを求めるための表
			#ifdef TABLE_SPP
			static sub_dimension *table_spp = NULL;
			#ifdef SPP_BIT
			int spp_bit = SPP_BIT; // 追加のビット数
			#else
			int spp_bit = 17; // 追加のビット数
			#endif
			static int table_spp_size = 0;
			if(table_spp == NULL) {
				table_spp_size = (1 << spp_bit) + 1;
				table_spp = make_table_for_enumeration_hamming(table_spp_size, spp_bit);
				fprintf(stderr, "make table for spplementary enumeation: spp_bit = %d, table_spp_size = %d\n", spp_bit, table_spp_size);
				#ifdef USE_DIFF_TABLE
				for(int i = 0; i < table_spp_size - 1; i++) {
					diff(&table_spp[i], &table_spp[i + 1], &table_spp[i]);
					if(table_spp[i + 1].num == 0) break; // 次が空集合になったら終了
				}
				#endif
			fprintf(stderr, "made table for hamming enumeration: spp_bit = %d\n", spp_bit);
			use_system("VmSize");
			}
			#endif

			#ifndef WITHOUT_IDX
			int *bd_idx = qs->idx;
			#endif
			int *bkt = bucket->bkt;

			// num_candidates  = 求めるデータ数（全体）
			// data_num[] = 求めたデータ番号を格納する配列
			// table[] = 部分集合を要素数の昇順に列挙するパターンを格納した配列
			// bd_idx[] = 相対位置でビット操作をするための距離下限の順位表
			// bkt[] = バケット表（bkt[s] = データをスケッチ順に並べたときに，スケッチ S のデータの先頭位置）
			// nt = スレッド数
			int num_data_thread = num_candidates / nt;	// スレッドが求めるデータ番号数の目標（これに達しないかもしれない）
			int t;										// スレッド番号
			int m;  									// スレッドが列挙したスケッチ数（パターン番号）
			int num_nonempty;							// スレッドが実際に列挙した空でないスケッチ数
			int k;  									// スレッドが求めたデータ番号数
			int table_size_thread = table_size / nt;	// スレッドが用いる列挙する部分集合のパターン配列の大きさ
			sub_dimension *q;							// スレッドが用いる列挙する部分集合のパターン配列
			int *d;										// スレッドが求めたデータ番号を格納する配列
			sketch_type base_mask = 0, base_mask2 = 0;	// 追加パターンで作成する mask
			int n_spp, n_spp2;							// 追加で使用するパターン番号
			sketch_type sk, mask;

			double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
			int num_sketches_thread = num_candidates * FACTOR_INF / ave_num / nt; // スレッドが求めるスケッチ数（目標のデータ数より多めになるようにしておく）
			int num_sketches = num_sketches_thread * nt; // 列挙で求めるスケッチの総数
			static int buff_pool_size = 0;
			#ifdef SELECT_SUM
				static sketch_with_priority_num *buff_pool = NULL;
				if(buff_pool == NULL || buff_pool_size < num_sketches) {
					if(buff_pool != NULL) {
						FREE(buff_pool, sizeof(sketch_with_priority_num) * buff_pool_size);
					}
					buff_pool_size = num_sketches;
					fprintf(stderr, "malloc answer buffer, ");
					buff_pool = MALLOC(sizeof(sketch_with_priority_num) * buff_pool_size);
					fprintf(stderr, "OK: buff_pool_size = %d, num_thread = %d\n", buff_pool_size, nt);
				}
				int num_enumerated_data_thread[nt];
				sketch_with_priority_num *buff;
			#else
				static answer_type *buff_pool = NULL;
				if(buff_pool == NULL || buff_pool_size < num_sketch_thread * nt) {
					if(buff_pool != NULL) {
						FREE(buff_pool, sizeof(answer_type) * buff_pool_size);
					}
					buff_pool_size = num_sketch_thread * nt;
					fprintf(stderr, "malloc answer buffer, ");
					buff_pool = MALLOC(sizeof(answer_type) * buff_pool_size);
					fprintf(stderr, "OK: buff_pool_size = %d, num_thread = %d\n", buff_pool_size, nt);
				}
				answer_type *buff;
			#endif

			#pragma omp parallel private(t, m, num_nonempty, q, d, base_mask, base_mask2, n_spp, n_spp2, buff, sk, mask)
			{
				t = omp_get_thread_num(); // スレッド番号の取得
				m = 0;
				num_nonempty = 0;
				q = table + t * table_size_thread;
				d = data_num + t * num_data_thread;
				base_mask = base_mask2 = 0;
				n_spp = n_spp2 = 1;
				buff = buff_pool + t * num_sketches_thread;
				num_enumerated_data_thread[t] = 0;

				for(m = 0; num_nonempty < num_sketches_thread; m++) {
					mask = base_mask;
					for(int j = 0; j < q[m].num; j++) {
						#ifndef WITHOUT_IDX
						mask |= (1 << bd_idx[(int)q[m].dim[j]]);
						#else
						mask |= (1 << (int)q[m].dim[j]);
						#endif
					}
					sk = qs->sketch ^ mask;
					// sk が空でなければ，buffに追加．ここでは，データ番号には展開しない
					if(bkt[sk + 1] > bkt[sk]) {
						#ifdef SELECT_SUM
						buff[num_nonempty++] = (sketch_with_priority_num){sk, priority(sk, qs), bkt[sk + 1] - bkt[sk]};
						#else
						buff[num_nonempty++] = (answer_type){sk, priority(sk, qs)};
						#endif
						num_enumerated_data_thread[t] += bkt[sk + 1] - bkt[sk];
					}
					if(q[m + 1].num == 0) { // 用意したパターンがなくなった．
						#ifdef TABLE_SPP
							if(table_spp[n_spp].num == 0) { // 追加分のパターンも使い切った．
								if(table_spp[n_spp2].num == 0) break; // 2回目の追加分のパターンも使い切った．
								n_spp = 0;
								n_spp2++;
								base_mask2 = 0;
								for(int m = 0; m < table_spp[n_spp2].num; m++) {
									#ifndef WITHOUT_IDX
									base_mask2 |= (1 << bd_idx[(int)table_spp[n_spp2].dim[m] + enum_dim + spp_bit]);
									#else
									base_mask2 |= (1 << ((int)table_spp[n_spp2].dim[m] + enum_dim + spp_bit));
									#endif
								}
							}
							base_mask = base_mask2;
							for(int m = 0; m < table_spp[n_spp].num; m++) {
								#ifndef WITHOUT_IDX
								base_mask |= (1 << bd_idx[(int)table_spp[n_spp].dim[m] + enum_dim]);
								#else
								base_mask |= (1 << ((int)table_spp[n_spp].dim[m] + enum_dim));
								#endif
							}
							m = -1; // forの再初期化で m++ となって，0になる．
							n_spp++;
						#else

						#endif
					}
				}

			}

			#ifdef SELECT_SUM
				int total_num_data_enumerated = 0, num_selected_sketches;
				for(int t = 0; t < nt; t++) {
					total_num_data_enumerated += num_enumerated_data_thread[t];
				}
				if(total_num_data_enumerated < num_candidates * FACTOR_INF) { // 列挙したスケッチのデータ総数が少ないときは，目標の候補数を減らす．
					num_candidates = total_num_data_enumerated / FACTOR_INF;
				}
				num_selected_sketches = quick_select_sum_k_sketch_with_priority_num(buff_pool, 0, buff_pool_size - 1, num_candidates);
			#else
				int sel = (num_sketch / nt < num_nonempty ? num_sketch / nt : num_nonempty);
				quick_select_k_answer(buff, 0, num_nonempty - 1, sel * 0.4);
			#endif

			#ifdef SELECT_SUM
				#if defined(_OPENMP) && NUM_THREADS > 1
					omp_set_num_threads(num_selected_sketches < NUM_THREADS ? num_selected_sketches : NUM_THREADS);
					nt = omp_get_max_threads(); 	// スレッド数を求める
					int sk_th[nt + 1]; // スレッドが列挙するスケッチの先頭位置．スレッド t が列挙するスケッチは，sk_th[t] から sk_th[t + 1] - 1 まで． 
					int num[nt]; // スレッド 0 から スレッド t - 1 が列挙するスケッチのデータ数の合計．スレッド t が列挙したスケッチのデータは，data_num[num[t]] から data_num[num[t + 1] - 1] まで．
					int p = 0, sum = 0;
					for(int t = 0; t < nt; t++) {
						while(sum < (num_candidates / nt) * t) {
							sum += buff_pool[p++].num;
						}
						sk_th[t] = p; 
						num[t] = sum;
					}
					sk_th[nt] = num_selected_sketches;
					#pragma omp parallel
					{
						int t = omp_get_thread_num(); // スレッド番号の取得
						int k = 0;
						int *d = data_num + num[t];
						for(int p = sk_th[t]; p < sk_th[t + 1]; p++) {
							sketch_type sk = (sketch_type)(buff_pool[p].sk);
							for(int j = bkt[sk]; j < bkt[sk + 1]; j++, k++) {
								#ifdef DATA_NUM_IN_SKETCH_ORDER
								d[k] = j;
								#else
								d[k] = bucket->idx[j];
								#endif
							}
						}
					}
				#else
					k = 0;
					for(int i = 0; i < num_enumerated_sketches; i++) {
						s = buff_pool[i].sk;
						for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates; j++, k++) {
							#ifdef DATA_NUM_IN_SKETCH_ORDER
							data_num[k] = j;
							#else
							data_num[k] = bucket->idx[j];
							#endif
						}
					}
				#endif
			#endif

			return num_candidates;
		}
		#endif // SELECT_BY_SINGLE
	#endif // PARA_ENUM_INF

#elif defined(WITH_DIFF)
// スケッチ列挙(Hamming)によるフィルタリング（バケット（配列 idx と bkt）利用）(1)
int filtering_by_sketch_enumeration_hamming(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates)
{
	#ifndef WITHOUT_IDX
	sketch_type sk; // 求めたスケッチ
	int *bd_idx = qs->idx;
	#endif
	int *bkt = bucket->bkt;

	static int first = 1;
	if(first) {
		fprintf(stderr, "enum_hamm single-thread. (1)\n");
		first = 0;
	}

	sk = qs->sketch; // 列挙の先頭は，質問のスケッチ．
	int k = 0; // 列挙したスケッチから得られたデータ番号の個数
	for(int j = bkt[sk]; j < bkt[sk + 1] && k < num_candidates; j++, k++) {
		#ifdef DATA_NUM_IN_SKETCH_ORDER
		data_num[k] = j;
		#else
		data_num[k] = bucket->idx[j];
		#endif
	}

	int bit; // ハミング距離順に列挙するときの直前と現在のビットパターン（この違いを元に列挙するスケッチを求める）

	for(int n = 1; k < num_candidates; n++) { // ビットパターンでbit_countがnであるものを求める
		for(bit = (1L << n) - 1; bit < (1L << PJT_DIM) &&  k < num_candidates; bit = nextComb(bit)) {
			#ifdef WITHOUT_IDX // bucket->idx を用いた相対位置によるビット操作を行わない
			sk = qs->sketch ^ bit;
			#else
			sketch_type diff = (unsigned)bit;
			sketch_type mask = 0;
			while(diff != 0) {
				int p = lsb_pos(diff);
				mask |= (1 << bd_idx[p]);
				diff ^= (1 << p);
			}
			sk = qs->sketch ^ mask;
			#endif
			for(int j = bkt[sk]; j < bkt[sk + 1] && k < num_candidates; j++, k++) {
				#ifdef DATA_NUM_IN_SKETCH_ORDER
				data_num[k] = j;
				#else
				data_num[k] = bucket->idx[j];
				#endif
			}
		}
	}
}
#else
// スケッチ列挙(Hamming)によるフィルタリング（バケット（配列 idx と bkt）利用）(2)
int filtering_by_sketch_enumeration_hamming(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates)
{
	sketch_type sk_n; 
	#ifndef WITHOUT_IDX
	sketch_type sk_p; // ひとつ前のスケッチ
	int *bd_idx = qs->idx;
	#endif
	int *bkt = bucket->bkt;

	static int first = 1;
	if(first) {
		#ifndef WITHOUT_IDX
		fprintf(stderr, "enum_hamm single-thread WITH_IDX. (2)\n");
		#else
		fprintf(stderr, "enum_hamm single-thread WITHOUT_IDX. (2)\n");
		#endif
		first = 0;
	}

	sk_n = qs->sketch; // 列挙の先頭は，質問のスケッチ．
	int k = 0; // 列挙したスケッチから得られたデータ番号の個数
	for(int j = bkt[sk_n]; j < bkt[sk_n + 1] && k < num_candidates; j++, k++) {
		#ifdef DATA_NUM_IN_SKETCH_ORDER
		data_num[k] = j;
		#else
		data_num[k] = bucket->idx[j];
		#endif
	}

	int bit_p, bit_n; // ハミング距離順に列挙するときの直前と現在のビットパターン（この違いを元に列挙するスケッチを求める）

	#ifndef WITHOUT_IDX
	sk_p = sk_n;
	#endif
	bit_p = 0;

	int n;
	for(n = 1; k < num_candidates; n++) { // ビットパターンでbit_countがnであるものを求める
		for(bit_n = (1L << n) - 1; bit_n < (1L << PJT_DIM) &&  k < num_candidates; bit_n = nextComb(bit_p)) {
			#ifdef WITHOUT_IDX // bucket->idx を用いた相対位置によるビット操作を行わない
			sk_n = qs->sketch ^ bit_n;
			#else
			sketch_type diff = (unsigned)(bit_p ^ bit_n);
			sketch_type mask = 0;
			while(diff != 0) {
				int p = lsb_pos(diff);
				mask |= (1 << bd_idx[p]);
				diff ^= (1 << p);
			}
			sketch_type sk_n = sk_p ^ mask;
			#endif
			for(int j = bkt[sk_n]; j < bkt[sk_n + 1] && k < num_candidates; j++, k++) {
				#ifdef DATA_NUM_IN_SKETCH_ORDER
				data_num[k] = j;
				#else
				data_num[k] = bucket->idx[j];
				#endif
			}
			#ifndef WITHOUT_IDX
			sk_p = sk_n;
			#endif
			bit_p = bit_n;
		}
	}
	return k;
}

int filtering_by_sketch_enumeration_hamming_interval(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates)
{
	sketch_type sk_n; 
	#ifndef WITHOUT_IDX
	sketch_type sk_p; // ひとつ前のスケッチ
	int *bd_idx = qs->idx;
	#endif
	int *bkt = bucket->bkt;
	int nt = ivl->nt;

	static int first = 1;
	if(first) {
		#ifndef WITHOUT_IDX
		fprintf(stderr, "enum_hamm_interval single-thread WITH_IDX. (2*)a\n");
		#else
		fprintf(stderr, "enum_hamm_interval single-thread WITHOUT_IDX. (2*)b\n");
		#endif
		first = 0;
	}

	sk_n = qs->sketch; // 列挙の先頭は，質問のスケッチ．
	int k = 0; // 列挙したスケッチから得られたデータ番号の個数
	int num_sk = 0; // 列挙したスケッチ数
	if(bkt[sk_n] < bkt[sk_n + 1]) { // 空でないバケツの始めと終わりを追加
		num_sk++;
		int t = num_sk % nt;
		int lg = num_sk / nt;
		ivl->lg[t] = lg;
		ivl->list[t * ivl->size + lg] = (interval){bkt[sk_n], bkt[sk_n + 1] - 1};
		k += bkt[sk_n + 1] - bkt[sk_n];
	}
	if(k >= num_candidates) return k;

	int bit_p, bit_n; // ハミング距離順に列挙するときの直前と現在のビットパターン（この違いを元に列挙するスケッチを求める）

	#ifndef WITHOUT_IDX
	sk_p = sk_n;
	#endif
	bit_p = 0;

	int n;
	for(n = 1; k < num_candidates; n++) { // ビットパターンでbit_countがnであるものを求める
		for(bit_n = (1L << n) - 1; bit_n < (1L << PJT_DIM) &&  k < num_candidates; bit_n = nextComb(bit_p)) {
			#ifdef WITHOUT_IDX // bucket->idx を用いた相対位置によるビット操作を行わない
			sk_n = qs->sketch ^ bit_n;
			#else
			sketch_type diff = (unsigned)(bit_p ^ bit_n);
			sketch_type mask = 0;
			while(diff != 0) {
				int p = lsb_pos(diff);
				mask |= (1 << bd_idx[p]);
				diff ^= (1 << p);
			}
			sketch_type sk_n = sk_p ^ mask;
			#endif
			if(bkt[sk_n] < bkt[sk_n + 1]) { // 空でないバケツの始めと終わりを追加
				num_sk++;
				int t = num_sk % nt;
				int lg = num_sk / nt;
				ivl->lg[t] = lg;
				ivl->list[t * ivl->size + lg] = (interval){bkt[sk_n], bkt[sk_n + 1] - 1};
				k += bkt[sk_n + 1] - bkt[sk_n];
			}
			if(k >= num_candidates) return k;
			#ifndef WITHOUT_IDX
			sk_p = sk_n;
			#endif
			bit_p = bit_n;
		}
	}
	return k;
}
#endif

#ifndef FACTOR_INF
#define FACTOR_INF 2.0
#endif

// スケッチ列挙(Hamming -> D~1 (quick_select_k))によるフィルタリング（バケット（配列 idx と bkt）利用）(3)
int filtering_by_sketch_enumeration_hamming_select(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates)
{
//	fprintf(stderr, "enum_hamm single-thread-select. FACTOR_INF = %d\n", FACTOR_INF); exit(0);
	double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
	int num_sketch = num_candidates / ave_num; // 求めるスケッチの個数
	int num_sketch_hamm = num_sketch * FACTOR_INF; // Hamming距離順の列挙によって求めるスケッチの個数（D~1で求める個数の1.5倍）
	static answer_type *buff = NULL;
	if(buff == NULL) {
		fprintf(stderr, "malloc answer buffer\n");
		buff = MALLOC(sizeof(answer_type) * num_sketch_hamm);
		fprintf(stderr, "malloc answer buffer OK\n");
	}

	static int first = 1;
	if(first) {
		fprintf(stderr, "enum_hamm single-thread. (3)\n");
		first = 0;
	}

	sketch_type sk_n; 
	#ifndef WITHOUT_IDX
	sketch_type sk_p; // ひとつ前のスケッチ
	int *bd_idx = qs->idx;
	#endif
	int *bkt = bucket->bkt;

	sk_n = qs->sketch; // 列挙の先頭は，質問のスケッチ．

	int bit_p, bit_n; // ハミング距離順に列挙するときの直前と現在のビットパターン（この違いを元に列挙するスケッチを求める）

	#ifndef WITHOUT_IDX
	sk_p = sk_n;
	#endif
	bit_p = 0;

	int num_nonempty = 0; // 列挙によって求めた空でないスケッチの個数
	// sk_n が空でなければ，buffに追加
	if(bkt[sk_n + 1] > bkt[sk_n]) {
		buff[num_nonempty++] = (answer_type){sk_n, priority(sk_n, qs)};
	}

	for(int n = 1; num_nonempty < num_sketch_hamm; n++) { // ビットパターンでbit_countがnであるものを求める
		for(bit_n = (1L << n) - 1; bit_n < (1L << PJT_DIM) &&  num_nonempty < num_sketch_hamm; bit_n = nextComb(bit_p)) {
			#ifdef WITHOUT_IDX // bucket->idx を用いた相対位置によるビット操作を行わない
				sk_n = qs->sketch ^ bit_n;
			#else
				sketch_type diff = (unsigned)(bit_p ^ bit_n);
				sketch_type mask = 0;
				while(diff != 0) {
					int p = lsb_pos(diff);
					mask |= (1 << bd_idx[p]);
					diff ^= (1 << p);
				}
				sketch_type sk_n = sk_p ^ mask;
			#endif
			if(bkt[sk_n + 1] > bkt[sk_n]) {
				buff[num_nonempty++] = (answer_type){sk_n, priority(sk_n, qs)};
			}
			#ifndef WITHOUT_IDX
			sk_p = sk_n;
			#endif
			bit_p = bit_n;
		}
	}

	quick_select_k_answer(buff, 0, num_sketch_hamm - 1, num_sketch);
	int k = 0; // 出力したデータ番号の個数
	for(int i = 0; k < num_candidates; i++) {
		sketch_type sk = (sketch_type)(buff[i].data_num);
		for(int j = bkt[sk]; j < bkt[sk + 1] && k < num_candidates; j++, k++) {
			#ifdef DATA_NUM_IN_SKETCH_ORDER
			data_num[k] = j;
			#else
			data_num[k] = bucket->idx[j];
			#endif
		}
	}
	return k;
}

// スケッチ列挙によるフィルタリング（バケット（配列 idx と bkt）利用）
int filtering_by_sketch_enumeration(struct_query_sketch *qs, struct_bucket *bucket, struct_que *que, int data_num[], int num_candidates)
{
	sketch_type s;
	QUE qu, qu2;
	int *bd = qs->bd, *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	s = qs->sketch;
	int k = 0;
	for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates; j++, k++) {
		#ifdef DATA_NUM_IN_SKETCH_ORDER
		data_num[k] = j;
		#else
		data_num[k] = bucket->idx[j];
		#endif
	}
	que->qsize = 0;
	qu.key = bd[bd_idx[0]];
	qu.sk = s ^ (1 << bd_idx[0]);
	qu.idx = 1;
	enq_p(&qu, que);

	while(deq_p(&qu, que) && k < num_candidates) {
		for(int j = bkt[qu.sk]; j < bkt[qu.sk + 1] && k < num_candidates; j++, k++) {
			#ifdef DATA_NUM_IN_SKETCH_ORDER
			data_num[k] = j;
			#else
			data_num[k] = bucket->idx[j];
			#endif
		}
		if(qu.idx < PJT_DIM) {
			s = qu.sk ^ (1 << bd_idx[qu.idx]);
			qu2.key = qu.key + bd[bd_idx[qu.idx]];
			qu2.sk = s;
			qu2.idx = qu.idx + 1;
			enq_p(&qu2, que);

			qu2.key = qu.key + bd[bd_idx[qu.idx]] - bd[bd_idx[qu.idx - 1]];
			qu2.sk = s ^ (1 << bd_idx[qu.idx - 1]);
			qu2.idx = qu.idx + 1;
			enq_p(&qu2, que);
		}
	}
	return k;
}

// スケッチ列挙によるフィルタリング（バケット（配列 idx と bkt）利用）
int filtering_by_sketch_enumeration_interval(struct_query_sketch *qs, struct_bucket *bucket, struct_que *que, interval_list *ivl, int num_candidates)
{
	sketch_type s;
	QUE qu, qu2;
	int *bd = qs->bd, *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	s = qs->sketch;
	int k = 0;
	int num_sk = 0; // 列挙した空でないスケッチ数
	if(bkt[s] < bkt[s + 1]) { // 空でないバケツの始めと長さ（run）を追加
		ivl->list[num_sk++] = (interval){bkt[s], bkt[s + 1] - bkt[s]};
		k += bkt[s + 1] - bkt[s];
	}

	que->qsize = 0;
	qu.key = bd[bd_idx[0]];
	qu.sk = s ^ (1 << bd_idx[0]);
	qu.idx = 1;
	enq_p(&qu, que);

	while(deq_p(&qu, que) && k < num_candidates) {
		if(bkt[qu.sk] < bkt[qu.sk + 1]) { // 空でないバケツの始めと長さ（run）を追加
			ivl->list[num_sk++] = (interval){bkt[qu.sk], bkt[qu.sk + 1] - bkt[qu.sk]};
			k += bkt[qu.sk + 1] - bkt[qu.sk];
		}
		if(qu.idx < PJT_DIM) {
			s = qu.sk ^ (1 << bd_idx[qu.idx]);
			qu2.key = qu.key + bd[bd_idx[qu.idx]];
			qu2.sk = s;
			qu2.idx = qu.idx + 1;
			enq_p(&qu2, que);

			qu2.key = qu.key + bd[bd_idx[qu.idx]] - bd[bd_idx[qu.idx - 1]];
			qu2.sk = s ^ (1 << bd_idx[qu.idx - 1]);
			qu2.idx = qu.idx + 1;
			enq_p(&qu2, que);
		}
	}
	ivl->lg[0] = num_sk;
	return k;
}

int filtering_by_sketch_enumeration_sketch(struct_query_sketch *qs, struct_bucket *bucket, struct_que *que, sketch_type sketch[], int num_candidates)
{
	double ave_num = (double)bucket->num_data / (1L << PJT_DIM);
	int num_sketch = num_candidates / ave_num;
	sketch_type s;
	QUE qu, qu2;
	int *bd = qs->bd, *bd_idx = qs->idx;

	s = qs->sketch;
	int num = 0;
	sketch[num++] = s;

	que->qsize = 0;
	qu.key = bd[bd_idx[0]];
	qu.sk = s ^ (1 << bd_idx[0]);
	qu.idx = 1;
	enq_p(&qu, que);

	while(deq_p(&qu, que) && num < num_sketch) {
		sketch[num++] = qu.sk;
		if(qu.idx < PJT_DIM) {
			s = qu.sk ^ (1 << bd_idx[qu.idx]);
			qu2.key = qu.key + bd[bd_idx[qu.idx]];
			qu2.sk = s;
			qu2.idx = qu.idx + 1;
			enq_p(&qu2, que);

			qu2.key = qu.key + bd[bd_idx[qu.idx]] - bd[bd_idx[qu.idx - 1]];
			qu2.sk = s ^ (1 << bd_idx[qu.idx - 1]);
			qu2.idx = qu.idx + 1;
			enq_p(&qu2, que);
		}
	}
	return num;
}

int filtering_by_sketch_enumeration_c2_n(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, int data_num[], int num_candidates)
{
	static int first = 1;
	if(first) {
		#ifndef WITHOUT_IDX
		fprintf(stderr, "filtering_by_sketch_enumeration_c2_n single-thread WITH_IDX. \n");
		#else
		fprintf(stderr, "filtering_by_sketch_enumeration_c2_n single-thread WITHOUT_IDX. \n");
		#endif
		first = 0;
	}
	sketch_type s;
	QUE_c2 qu, qu2;
	int *bd = qs->bd, *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	#ifdef SELECT_SUM
	static sketch_with_priority_num *sk_buff = NULL;
	static int sk_buff_size = 0;
	if(sk_buff_size < num_candidates) {
		if(sk_buff != NULL) {
			FREE(sk_buff, sizeof(sketch_with_priority_num) * sk_buff_size);
		}
		sk_buff_size = num_candidates;
		sk_buff = MALLOC(sizeof(sketch_with_priority_num) * sk_buff_size);
	}
	dist_type priority;
	int num_enumerated_sketches = 0;
	#endif

	s = qs->sketch; // 先頭は質問のスケッチ
	int k = 0;

	#ifdef SELECT_SUM
	priority = 0;
	if(bkt[s + 1] > bkt[s]) { // sが空でなければバッファーに追加
		sk_buff[num_enumerated_sketches++] = (sketch_with_priority_num){s, priority, bkt[s + 1] - bkt[s]};
		k += bkt[s + 1] - bkt[s]; // k が目標の num_candidates を超えるかも？
	}
	#else
	for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates; j++, k++) {
		#ifdef DATA_NUM_IN_SKETCH_ORDER
		data_num[k] = j;
		#else
		data_num[k] = bucket->idx[j];
		#endif
	}
	#endif

	s = s ^ (1 <<  bd_idx[0]); // 先頭の次は、質問のスケッチと距離下限が最小のビットだけが異なるもの

	#ifdef SELECT_SUM
	if(k < num_candidates) {
		priority = qs->bd[bd_idx[0]];
		if(bkt[s + 1] > bkt[s]) { // sが空でなければバッファーに追加
			sk_buff[num_enumerated_sketches++] = (sketch_with_priority_num){s, priority, bkt[s + 1] - bkt[s]};
			k += bkt[s + 1] - bkt[s]; // k が目標の num_candidates を超えるかも？
		}
	}
	#else
	for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates; j++, k++) {
		#ifdef DATA_NUM_IN_SKETCH_ORDER
		data_num[k] = j;
		#else
		data_num[k] = bucket->idx[j];
		#endif
	}
	#endif

	make_empty_que_c2_n(que);

	// enq pattern of 0...10
	qu.cursor = new_que_e2_n(que);
	qu.key = bd[bd_idx[1]];
	que->details[qu.cursor].sk = qs->sketch ^ (1 << bd_idx[1]);
	que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
	enq_c2_n(&qu, que);		

	while(deq_c2_n(&qu, que) && k < num_candidates) {

		s = que->details[qu.cursor].sk; // 列挙のつぎのスケッチ

		#ifdef SELECT_SUM
		priority = qu.key;
		if(bkt[s + 1] > bkt[s]) { // sが空でなければバッファーに追加
			sk_buff[num_enumerated_sketches++] = (sketch_with_priority_num){s, priority, bkt[s + 1] - bkt[s]};
			k += bkt[s + 1] - bkt[s]; // k が目標の num_candidates を超えるかも？
		}
		#else
		for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates; j++, k++) {
			#ifdef DATA_NUM_IN_SKETCH_ORDER
			data_num[k] = j;
			#else
			data_num[k] = bucket->idx[j];
			#endif
		}
		#endif

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

	#ifdef SELECT_SUM
		#if defined(_OPENMP) && NUM_THREADS > 1
			omp_set_num_threads(num_enumerated_sketches < NUM_THREADS ? num_enumerated_sketches : NUM_THREADS);
			int nt = omp_get_max_threads(); 	// スレッド数を求める
			int sk_th[nt + 1]; // スレッドが列挙するスケッチの先頭位置．スレッド t が列挙するスケッチは，sk_th[t] から sk_th[t + 1] - 1 まで． 
			int num[nt]; // スレッド 0 から スレッド t - 1 が列挙するスケッチのデータ数の合計．スレッド t が列挙したスケッチのデータは，data_num[num[t]] から data_num[num[t + 1] - 1] まで．
			int p = 0, sum = 0;
			for(int t = 0; t < nt; t++) {
				while(sum < (num_candidates / nt) * t) {
					sum += sk_buff[p++].num;
				}
				sk_th[t] = p; 
				num[t] = sum;
//				fprintf(stderr, "sk_th[%d] = %d, num[%d] = %d\n", t, sk_th[t], t, num[t]);
			}
			sk_th[nt] = num_enumerated_sketches;
			#pragma omp parallel
			{
				int t = omp_get_thread_num(); // スレッド番号の取得
				int k = 0;
				int *d = data_num + num[t];
				for(int p = sk_th[t]; p < sk_th[t + 1]; p++) {
					sketch_type sk = (sketch_type)(sk_buff[p].sk);
					for(int j = bkt[sk]; j < bkt[sk + 1]; j++, k++) {
						#ifdef DATA_NUM_IN_SKETCH_ORDER
						d[k] = j;
						#else
						d[k] = bucket->idx[j];
						#endif
					}
				}
			}
		#else
			k = 0;
			for(int i = 0; i < num_enumerated_sketches; i++) {
				s = sk_buff[i].sk;
				for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates; j++, k++) {
					#ifdef DATA_NUM_IN_SKETCH_ORDER
					data_num[k] = j;
					#else
					data_num[k] = bucket->idx[j];
					#endif
				}
			}
		#endif
	#endif
	return k;
}

#if !defined(THREAD_PLUS) || THREAD_PLUS == 0
#if PARA_ENUM_INF == 0
int filtering_by_sketch_enumeration_c2_n_interval(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, interval_list *ivl, int num_candidates)
#else
int filtering_by_sketch_enumeration_c2_n_interval_st(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, interval_list *ivl, int num_candidates)
#endif
{
	static int first = 1;
	if(first) {
		#ifndef WITHOUT_IDX
		fprintf(stderr, "(*)filtering_by_sketch_enumeration_c2_n_interval single-thread WITH_IDX. \n");
		#else
		fprintf(stderr, "filtering_by_sketch_enumeration_c2_n_interval single-thread WITHOUT_IDX. \n");
		#endif
		first = 0;
	}

	sketch_type s;
//	dist_type priority;
	QUE_c2 qu, qu2;
	int *bd = qs->bd, *bd_idx = qs->idx;
	int *bkt = bucket->bkt;
	int nt = ivl->nt;

	s = qs->sketch; // 先頭は質問のスケッチ
	int k = 0;
	int num_sk = 0; // 列挙したスケッチ数
	if(bkt[s] < bkt[s + 1]) { // 空でないバケツの始めと終わりを追加
		num_sk++;
		int t = num_sk % nt;
		int lg = num_sk / nt;
		ivl->lg[t] = lg;
		ivl->list[t * ivl->size + lg] = (interval){bkt[s], bkt[s + 1] - bkt[s]};
		k += bkt[s + 1] - bkt[s];
	}

	s = s ^ (1 <<  bd_idx[0]); // 先頭の次は、質問のスケッチと距離下限が最小のビットだけが異なるもの
	if(bkt[s] < bkt[s + 1]) { // 空でないバケツの始めと終わりを追加
		num_sk++;
		int t = num_sk % nt;
		int lg = num_sk / nt;
		ivl->lg[t] = lg;
		ivl->list[t * ivl->size + lg] = (interval){bkt[s], bkt[s + 1] - bkt[s]};
		k += bkt[s + 1] - bkt[s];
	}

	make_empty_que_c2_n(que);

	// enq pattern of 0...10
	qu.cursor = new_que_e2_n(que);
	qu.key = bd[bd_idx[1]];
	que->details[qu.cursor].sk = qs->sketch ^ (1 << bd_idx[1]);
	que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
	enq_c2_n(&qu, que);		

	while(deq_c2_n(&qu, que) && k < num_candidates) {

		s = que->details[qu.cursor].sk; // 列挙のつぎのスケッチ
		if(bkt[s] < bkt[s + 1]) { // 空でないバケツの始めと終わりを追加
			num_sk++;
			int t = num_sk % nt;
			int lg = num_sk / nt;
			ivl->lg[t] = lg;
			ivl->list[t * ivl->size + lg] = (interval){bkt[s], bkt[s + 1] - bkt[s]};
			k += bkt[s + 1] - bkt[s];
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

	return k;
}
#else
#if PARA_ENUM_INF == 0
int filtering_by_sketch_enumeration_c2_n_interval(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, interval_list *ivl, int num_candidates)
// 列挙による 1st-filtering は，直列処理（PARA_ENUM_INF = 0）．
// ただし，PARALLEL_ENUM = 1, 2, ... のときは，double-filtering での 2nd-filtering では，並列処理を行うので，
// ivl は，その並列処理に対応できるようにしておく．
// ここでは，THREAD_PLUS = 1, 2, ... であるので，
// もし THREAD_PLUS < PARALLEL_ENUM のときは，ivl->nt = 2 ^ PARALLEL_ENUM として，ここでのフィルタリングの出力は，複数（2 ^ (PARALLEL_ENUM - THREAD_PLUS)) に分ける．
// もし THREAD_PLUS > PARALLEL_ENUM のときは，未対応．
{
	static int first = 1;
	if(first) {
		#ifndef WITHOUT_IDX
		fprintf(stderr, "filtering_by_sketch_enumeration_c2_n_interval single-thread WITH_IDX: thread_plus = %d", THREAD_PLUS);
		fprintf(stderr, ", !INTERVAL_WITH_PRIORITY\n");
		fprintf(stderr, "ivl->nt = %d, ivl->size = %d, num_candidates = %d\n", ivl->nt, ivl->size, num_candidates);
		#else
		fprintf(stderr, "filtering_by_sketch_enumeration_c2_n_interval %d-thread WITHOUT_IDX. \n", nt);
		#endif
		first = 0;
	}

	sketch_type s;
	QUE_c2 qu, qu2;
	int *bd = qs->bd, *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	int mu_thread_size = 1 << THREAD_PLUS;

	sketch_type mu_thread[mu_thread_size]; // スレッドに割り当てられた下位 n-bit のパターン（or 質問スケッチとのXOR）
	// このパターンの順序は，D_INF 用の grey code 生成ルールを用いて準備する
	mu_thread[0] = 0;
	for(int t = 0; t < mu_thread_size - 1; t++) {
		mu_thread[t + 1] = mu_thread[t] ^ (1 << bd_idx[bit_count(t ^ (t + 1)) - 1]);
	}
	#ifdef PARALLEL_ENUM
//	fprintf(stderr, "mu_thread OK\n");
//	fprintf(stderr, "lg[t]: t = 0, ... , %d ", 1 << PARALLEL_ENUM);
	for(int tt = 0; tt < (1 << PARALLEL_ENUM); tt++) {
		ivl->lg[tt] = 0;
	}
	#else
//	fprintf(stderr, "!defined(PARALLEL_ENUM), THREAD_PLUS = %d, mu_thread_size = %d, ivl->nt = %d\n", THREAD_PLUS, mu_thread_size, ivl->nt); exit(0);
	for(int t = 0; t < mu_thread_size; t++) {
		ivl->lg[t] = 0;
	}
	#endif
//	fprintf(stderr, "OK\n");

	// 下位 THREAD_PLUS ビットは，スレッド毎に割り当てた mu_thread[t] との XOR で求める．
	// それ以降の PJT_DIM - THREAD_PLUS ビットのみを変化させた，スケッチを列挙する（下位 THREAD_PLUS ビットは，質問のまま）

	int nn = THREAD_PLUS;
	int nc = 0;
	int nc0;

	s = qs->sketch; // 先頭は質問のスケッチ
	nc0 = 0;

// mu_thread_size = 1 << THREAD_PLUS 個に分割している．
// PARALELL_ENUM > THREAD_PLUS のときは，t 番目の分割で求められる区間の出力は，複数のリストに交互に格納していく．
// 複数リストの個数 num_list は，1 << (PARALLEL_ENUM - THREAD_PLUS) となる．

	#ifdef PARALLEL_ENUM
	int num_list = (1 << (PARALLEL_ENUM - THREAD_PLUS));
	char output_t[mu_thread_size];
	for(int t = 0; t < mu_thread_size; t++) {
		output_t[t] = 0;
	} 
	#endif
	for(int t = 0; t < mu_thread_size; t++) {
		sketch_type sk = s ^ mu_thread[t];
		if(bkt[sk + 1] > bkt[sk]) {
//			fprintf(stderr, "t = %d, output_t[%d] = %d, tt = %d, num_list = %d\n", t, t, output_t[t], t * num_list + output_t[t], num_list);
			#ifdef PARALLEL_ENUM
			int tt = t * num_list + output_t[t]++;
			if(output_t[t] == num_list) output_t[t] = 0;
			#else
			int tt = t;
			#endif
			interval *buff = ivl->list + tt * ivl->size;
			buff[ivl->lg[tt]++] = (interval){bkt[sk], bkt[sk + 1] - bkt[sk]};
			nc0 += bkt[sk + 1] - bkt[sk];
		}
	}
	nc += nc0;

	if(nc >= num_candidates) goto ret;
	s = s ^ (1 <<  bd_idx[nn]); // 先頭の次は、質問のスケッチと距離下限が最小のビットだけが異なるもの

	nc0 = 0;
	for(int t = 0; t < mu_thread_size; t++) {
		sketch_type sk = s ^ mu_thread[t];
		if(bkt[sk + 1] > bkt[sk]) {
			#ifdef PARALLEL_ENUM
			int tt = t * num_list + output_t[t]++;
			if(output_t[t] == num_list) output_t[t] = 0;
			#else
			int tt = t;
			#endif
			interval *buff = ivl->list + tt * ivl->size;
			buff[ivl->lg[tt]++] = (interval){bkt[sk], bkt[sk + 1] - bkt[sk]};
			nc0 += bkt[sk + 1] - bkt[sk];
		}
	}
	nc += nc0;
	if(nc >= num_candidates) goto ret;

	make_empty_que_c2_n(que);

	// enq pattern of 0...10
	qu.cursor = new_que_e2_n(que);
	qu.key = bd[bd_idx[nn + 1]];
	que->details[qu.cursor].sk = qs->sketch ^ (1 << bd_idx[nn + 1]);
	que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
	enq_c2_n(&qu, que);		

	while(deq_c2_n(&qu, que)) {
		s = que->details[qu.cursor].sk; // 列挙のつぎのスケッチ
		nc0 = 0;
		for(int t = 0; t < mu_thread_size; t++) {
			sketch_type sk = s ^ mu_thread[t];
			if(bkt[sk + 1] > bkt[sk]) {
				#ifdef PARALLEL_ENUM
				int tt = t * num_list + output_t[t]++;
				if(output_t[t] == num_list) output_t[t] = 0;
				#else
				int tt = t;
				#endif
				interval *buff = ivl->list + tt * ivl->size;
				buff[ivl->lg[tt]++] = (interval){bkt[sk], bkt[sk + 1] - bkt[sk]};
				nc0 += bkt[sk + 1] - bkt[sk];
			}
		}
		nc += nc0;
		if(nc >= num_candidates) goto ret;

		switch(que->details[qu.cursor].pt & 15) {
		case 0: // X0000 -> enq(X0001) and enq(Y10^{m+1}) if X0000 = Y010^m
		case 8: // X1000 -> enq(X1001) and enq(Y10^{m+1}) if X0000 = Y010^m
			{
				int m = lsb_pos(que->details[qu.cursor].pt);
				if(m > 0 && nn + m < PJT_DIM - 1 && !(que->details[qu.cursor].pt & (1 << (m + 1)))) {
					// Y010^m -> Y10^{m+1}
					qu2.cursor = new_que_e2_n(que);
					qu2.key = qu.key + bd[bd_idx[nn + m + 1]] - bd[bd_idx[nn + m]];
					que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + m + 1])) ^ (1 << bd_idx[nn + m]);
					que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
					// Y010^m -> Y010^{m-1}1
					qu.key = qu.key + bd[bd_idx[nn + 0]];
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 0]);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					enq_c2_n(&qu2, que);
				} else {
					qu.key = qu.key + bd[bd_idx[nn + 0]];
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 0]);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
				}
			}
			break;
		case 4:  // X0100 -> enq(X0101) and enq(X1000)
			// X1000
			qu2.cursor = new_que_e2_n(que);
			qu2.key = qu.key + bd[bd_idx[nn + 3]] - bd[bd_idx[nn + 2]];
			que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 3])) ^ (1 << bd_idx[nn + 2]);
			que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
			// X0101
			qu.key = qu.key + bd[bd_idx[nn + 0]];
			que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 0]);
			que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
			enq_c2_n(&qu, que);
			enq_c2_n(&qu2, que);
			break;
		case 1:  // X0001 -> enq(X0010)
		case 5:  // X0101 -> enq(X0110)
		case 9:  // X1001 -> enq(X1010)
		case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
			qu.key = qu.key + bd[bd_idx[nn + 1]] - bd[bd_idx[nn + 0]];
			que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 1])) ^ (1 << bd_idx[nn + 0]);
			que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
			enq_c2_n(&qu, que);
			break;
		case 2:  // X0010 -> enq(X0011) and enq(X0100)
		case 10: // X1010 -> enq(X1011) and enq(X1100)
			// X0100 and X1100
			qu2.cursor = new_que_e2_n(que);
			qu2.key = qu.key +  bd[bd_idx[nn + 2]] -  bd[bd_idx[nn + 1]];
			que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 2])) ^ (1 << bd_idx[nn + 1]);
			que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
			// X0011 and X1011
			qu.key = qu.key + bd[bd_idx[nn + 0]];
			que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[nn + 0]);
			que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
			enq_c2_n(&qu, que);
			enq_c2_n(&qu2, que);
			break;
		case 6:  // X0110 -> enq(X0111)
		case 12: // X1100 -> enq(X1101)
		case 14: // X1110 -> enq(10111)
			qu.key = qu.key + bd[bd_idx[nn + 0]];
			que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[nn + 0]);
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
ret:

	return nc;
}
#endif
#endif

#if PARA_ENUM_INF > 0
#if defined(USE_MU_COMMON)
// with post-selection (using mu_common)
int filtering_by_sketch_enumeration_c2_n_interval(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, interval_list *ivl, int num_candidates)
{
	int n = PARA_ENUM_INF;
	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif
	static int first = 1;
	if(first) {
		#ifndef WITHOUT_IDX
		fprintf(stderr, "(b with post-selection using mu_common) filtering_by_sketch_enumeration_c2_n_interval %d-thread WITH_IDX. \n", nt);
		#else
		fprintf(stderr, "filtering_by_sketch_enumeration_c2_n_interval %d-thread WITHOUT_IDX. \n", nt);
		#endif
		first = 0;
	}

	sketch_type s;
	QUE_c2 qu, qu2;
	int *bd = qs->bd, *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	#ifndef THREAD_PLUS
	int mu_thread_size = nt;
	#else
	int mu_thread_size = nt * (1 << THREAD_PLUS);
	#endif

	sketch_type mu_thread[mu_thread_size]; // スレッドに割り当てられた下位 n-bit のパターン（or 質問スケッチとのXOR）
	// このパターンの順序は，D_INF 用の grey code 生成ルールを用いて準備する
	mu_thread[0] = 0;
	for(int t = 0; t < mu_thread_size - 1; t++) {
		mu_thread[t + 1] = mu_thread[t] ^ (1 << bd_idx[bit_count(t ^ (t + 1)) - 1]);
	}
	dist_type pr_thread[mu_thread_size];
	for(int t = 0; t < mu_thread_size; t++) {
		pr_thread[t] = priority(mu_thread[t] ^ qs->sketch, qs);
	}
//#define MU_SORT
#ifdef MU_SORT
	int idx_mu[mu_thread_size];
	for(int t = 0; t < mu_thread_size; t++) {
		idx_mu[t] = t;
	}
	quick_sort(idx_mu, pr_thread, 0, mu_thread_size - 1);
	sketch_type mu_temp[mu_thread_size];
	dist_type pr_temp[mu_thread_size];
	for(int t = 0; t < mu_thread_size; t++) {
		mu_temp[t] = mu_thread[idx_mu[t]];
		pr_temp[t] = pr_thread[idx_mu[t]];
	}
	for(int t = 0; t < mu_thread_size; t++) {
		mu_thread[t] = mu_temp[t];
		pr_thread[t] = pr_temp[t];
	}
#endif

	dist_type max_priority = priority(-1 ^ qs->sketch, qs);

	int num_nonempty[mu_thread_size];								// スレッドが実際に列挙した空でないスケッチ数
	int num_data_of_priority[mu_thread_size][max_priority + 1];		// スレッドが求めた各priority毎のデータ数
	int lg_of_priority[mu_thread_size][max_priority + 1];			// スレッドが求めた各priority毎のリスト長（interval数）
	int total_enum_data[max_priority + 1];				// 求めた各priority毎のデータ数の合計
	int final_priority = 0;
	for(int t = 0; t < mu_thread_size; t++) {
		num_nonempty[t] = 0;
	}
	for(int p = 0; p <= max_priority; p++) {
		total_enum_data[p] = 0;
		for(int t = 0; t < mu_thread_size; t++) {
			num_data_of_priority[t][p] = lg_of_priority[t][p] = 0;
		}
	}

	static sketch_type *mu_common = NULL;			// 一回の処理で求めた共通マスクの配列
	static dist_type *pr_common = NULL;
	static int common_size = 0;
	if(common_size < num_candidates) {
		if(mu_common != NULL) {
			FREE(mu_common, sizeof(sketch_type) * common_size);
			FREE(pr_common, sizeof(dist_type) * common_size);
		}
		mu_common = MALLOC(sizeof(sketch_type) * num_candidates);
		pr_common = MALLOC(sizeof(dist_type) * num_candidates);
		common_size = num_candidates;
	}

	// 下位 n ビットは，スレッド毎に割り当てた mu_thread[t] との XOR で求める．
	// それ以降の PJT_DIM - n ビットのみを変化させた，スケッチを列挙する（下位 n ビットは，質問のまま）

	s = qs->sketch; // 先頭は質問のスケッチ
	int num_sk = 0; // 列挙したスケッチ数
	mu_common[num_sk] = s;
	pr_common[num_sk++] = 0;
	int chk_sum = 0;
	int num_checked = 0; // chk_sum を調べたスケッチ数

	#ifndef THREAD_PLUS
	int nn = n;
	#else
	int nn = n + THREAD_PLUS;
	#endif

	s = s ^ (1 <<  bd_idx[nn]); // 先頭の次は、質問のスケッチと距離下限が最小のビットだけが異なるもの
	mu_common[num_sk] = s;
	pr_common[num_sk++] = qs->bd[bd_idx[nn]];

	make_empty_que_c2_n(que);

	// enq pattern of 0...10
	qu.cursor = new_que_e2_n(que);
	qu.key = bd[bd_idx[nn + 1]];
	que->details[qu.cursor].sk = qs->sketch ^ (1 << bd_idx[nn + 1]);
	que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
	enq_c2_n(&qu, que);		

	while(deq_c2_n(&qu, que)) {

		s = que->details[qu.cursor].sk; // 列挙のつぎのスケッチ
		mu_common[num_sk] = s;
		pr_common[num_sk++] = qu.key;

		if(num_sk >= num_checked + FACTOR_INF2) {
			// チェック済みの後で求めた共通マスク数が FACTOR_INF2 を超えたので，区間リストを求める．
			int sum_part = 0; // 求めた区間リストの区間長の合計（＝候補データ数）
			#pragma omp parallel reduction (+:sum_part)
			{
				#ifndef THREAD_PLUS
				int t = omp_get_thread_num();	// スレッド番号
				int ne = num_nonempty[t];		// 各スレッドで求めた空でないスケッチ数
				int nd = 0;						// 一回の処理で各スレッドで求めたデータ数
				int *nd_p = num_data_of_priority[t];
				int *lg_p = lg_of_priority[t];
				interval *buff = ivl->list + t * ivl->size;		// スレッドに割当てられた区間リスト
				for(int j = num_checked; j < num_sk; j++) {
					sketch_type sk = mu_common[j] ^ mu_thread[t];
					if(bkt[sk + 1] > bkt[sk]) {
						dist_type pr = pr_common[j] + pr_thread[t];
						buff[ne++] = (interval){bkt[sk], bkt[sk + 1] - bkt[sk]};
						nd += bkt[sk + 1] - bkt[sk];
						sum_part += bkt[sk + 1] - bkt[sk];
						nd_p[pr] += bkt[sk + 1] - bkt[sk];
						lg_p[pr]++;
					}
				}
				num_nonempty[t] = ne;
				#else
				int t = omp_get_thread_num();	// スレッド番号
				for(int tn = t * (1 << THREAD_PLUS); tn < (t + 1) * (1 << THREAD_PLUS); tn++) {
					int ne = num_nonempty[tn];		// 各スレッドで求めた空でないスケッチ数
					int nd = 0;						// 一回の処理で各スレッドで求めたデータ数
					int *nd_p = num_data_of_priority[tn];
					int *lg_p = lg_of_priority[tn];
					interval *buff = ivl->list + tn * ivl->size;
					for(int j = num_checked; j < num_sk; j++) {
						sketch_type sk = mu_common[j] ^ mu_thread[tn];
						if(bkt[sk + 1] > bkt[sk]) {
							dist_type pr = pr_common[j] + pr_thread[tn];
							buff[ne++] = (interval){bkt[sk], bkt[sk + 1] - bkt[sk]};
							nd += bkt[sk + 1] - bkt[sk];
							sum_part += bkt[sk + 1] - bkt[sk];
							nd_p[pr] += bkt[sk + 1] - bkt[sk];
							lg_p[pr]++;
						}
					}
					num_nonempty[tn] = ne;
				}
				#endif
			}
			// num_checked までは，データ数を確認済み（chk_sum）．
			// num_sk までは，共通マスク（mu_common）作成済み．
			// num_checked 以降 num_sk までのデータ数は，sum_part
			if(chk_sum + sum_part >= num_candidates * FACTOR_INF) {
				// max_priority を求めなおす
				int sd;
				#ifndef THREAD_PLUS
				int NT = nt;
				#else
				int NT = nt * (1 << THREAD_PLUS);
				#endif
				while(1) {
					sd = 0;
					for(int t = 0; t < NT; t++) {
						sd += num_data_of_priority[t][max_priority];
					}
					if(sd != 0) break;
					max_priority--;
				}
				// num_data_of_priority[t][pr] を pr = 0, ... , pr の累積に変更する．
				// lg_of_priority[t][pr] も累積
				for(int t = 0; t < NT; t++) {
					for(int pr = 1; pr <= max_priority; pr++) {
						num_data_of_priority[t][pr] += num_data_of_priority[t][pr - 1];
						lg_of_priority[t][pr] += lg_of_priority[t][pr - 1];
					}
				}
				final_priority = 0;
				for(int pr = 0; pr <= max_priority; pr++) {
					total_enum_data[pr] = 0;
					for(int t = 0; t < NT; t++) {
						total_enum_data[pr] += num_data_of_priority[t][pr];
					}
					if(total_enum_data[pr] < num_candidates) final_priority = pr;
				}
				final_priority++;
				for(int t = 0; t < NT; t++) {
					ivl->lg[t] = lg_of_priority[t][final_priority];
				}
				for(int t = NT - 1; t >= 0; t--) {
					if(total_enum_data[final_priority] - (num_data_of_priority[t][final_priority] - num_data_of_priority[t][final_priority - 1]) >= num_candidates) {
						ivl->lg[t] = lg_of_priority[t][final_priority - 1];
						total_enum_data[final_priority] -= (num_data_of_priority[t][final_priority] - num_data_of_priority[t][final_priority - 1]);
					} else {

					}
				}
				break;
			}

			num_checked = num_sk;
			chk_sum += sum_part;
		}
		switch(que->details[qu.cursor].pt & 15) {
		case 0: // X0000 -> enq(X0001) and enq(Y10^{m+1}) if X0000 = Y010^m
		case 8: // X1000 -> enq(X1001) and enq(Y10^{m+1}) if X0000 = Y010^m
			{
				int m = lsb_pos(que->details[qu.cursor].pt);
				if(m > 0 && nn + m < PJT_DIM - 1 && !(que->details[qu.cursor].pt & (1 << (m + 1)))) {
					// Y010^m -> Y10^{m+1}
					qu2.cursor = new_que_e2_n(que);
					qu2.key = qu.key + bd[bd_idx[nn + m + 1]] - bd[bd_idx[nn + m]];
					que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + m + 1])) ^ (1 << bd_idx[nn + m]);
					que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
					// Y010^m -> Y010^{m-1}1
					qu.key = qu.key + bd[bd_idx[nn + 0]];
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 0]);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
					enq_c2_n(&qu2, que);
				} else {
					qu.key = qu.key + bd[bd_idx[nn + 0]];
					que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 0]);
					que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
					enq_c2_n(&qu, que);
				}
			}
			break;
		case 4:  // X0100 -> enq(X0101) and enq(X1000)
			// X1000
			qu2.cursor = new_que_e2_n(que);
			qu2.key = qu.key + bd[bd_idx[nn + 3]] - bd[bd_idx[nn + 2]];
			que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 3])) ^ (1 << bd_idx[nn + 2]);
			que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
			// X0101
			qu.key = qu.key + bd[bd_idx[nn + 0]];
			que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 0]);
			que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
			enq_c2_n(&qu, que);
			enq_c2_n(&qu2, que);
			break;
		case 1:  // X0001 -> enq(X0010)
		case 5:  // X0101 -> enq(X0110)
		case 9:  // X1001 -> enq(X1010)
		case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
			qu.key = qu.key + bd[bd_idx[nn + 1]] - bd[bd_idx[nn + 0]];
			que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 1])) ^ (1 << bd_idx[nn + 0]);
			que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
			enq_c2_n(&qu, que);
			break;
		case 2:  // X0010 -> enq(X0011) and enq(X0100)
		case 10: // X1010 -> enq(X1011) and enq(X1100)
			// X0100 and X1100
			qu2.cursor = new_que_e2_n(que);
			qu2.key = qu.key +  bd[bd_idx[nn + 2]] -  bd[bd_idx[nn + 1]];
			que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 2])) ^ (1 << bd_idx[nn + 1]);
			que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
			// X0011 and X1011
			qu.key = qu.key + bd[bd_idx[nn + 0]];
			que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[nn + 0]);
			que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
			enq_c2_n(&qu, que);
			enq_c2_n(&qu2, que);
			break;
		case 6:  // X0110 -> enq(X0111)
		case 12: // X1100 -> enq(X1101)
		case 14: // X1110 -> enq(10111)
			qu.key = qu.key + bd[bd_idx[nn + 0]];
			que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[nn + 0]);
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
	balance_interval_list(ivl);
	return total_enum_data[final_priority];
}
#else
// with post-selection (without using mu_common)
int filtering_by_sketch_enumeration_c2_n_interval(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que_pool[], interval_list *ivl, int num_candidates)
{
	int n = PARA_ENUM_INF;
	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif

	static int first = 1;
	if(first) {
		#ifndef WITHOUT_IDX
		fprintf(stderr, "(c: with post-selection by independent enumeration) filtering_by_sketch_enumeration_c2_n_interval %d-thread WITH_IDX. \n", nt);
		#else
		fprintf(stderr, "filtering_by_sketch_enumeration_c2_n_interval %d-thread WITHOUT_IDX. \n", nt);
		#endif
		first = 0;
	}

	#ifndef THREAD_PLUS
	int mu_thread_size = nt;
	#else
	int mu_thread_size = nt * (1 << THREAD_PLUS);
	#endif

	sketch_type mu_thread[mu_thread_size]; // スレッドに割り当てられた下位 n-bit のパターン（or 質問スケッチとのXOR）
	// このパターンの順序は，D_INF 用の grey code 生成ルールを用いて準備する
	mu_thread[0] = 0;
	for(int t = 0; t < mu_thread_size - 1; t++) {
		mu_thread[t + 1] = mu_thread[t] ^ (1 << qs->idx[bit_count(t ^ (t + 1)) - 1]);
	}
	dist_type pr_thread[mu_thread_size];
	for(int t = 0; t < mu_thread_size; t++) {
		pr_thread[t] = priority(mu_thread[t] ^ qs->sketch, qs);
	}
//#define MU_SORT
#ifdef MU_SORT
	int idx_mu[mu_thread_size];
	for(int t = 0; t < mu_thread_size; t++) {
		idx_mu[t] = t;
	}
	quick_sort(idx_mu, pr_thread, 0, mu_thread_size - 1);
	sketch_type mu_temp[mu_thread_size];
	dist_type pr_temp[mu_thread_size];
	for(int t = 0; t < mu_thread_size; t++) {
		mu_temp[t] = mu_thread[idx_mu[t]];
		pr_temp[t] = pr_thread[idx_mu[t]];
	}
	for(int t = 0; t < mu_thread_size; t++) {
		mu_thread[t] = mu_temp[t];
		pr_thread[t] = pr_temp[t];
//		printf("pr[%d] = %d\n", t, pr_thread[t]);
	}
//	getchar();
#endif

	dist_type max_priority = priority(-1 ^ qs->sketch, qs);

	#if !defined(THREAD_PLUS) || THREAD_PLUS == 0
//	int num_nonempty[nt];								// スレッドが実際に列挙した空でないスケッチ数
	int num_data_of_priority[nt][max_priority + 1];		// スレッドが求めた各priority毎のinterval数
	int lg_of_priority[nt][max_priority + 1];			// スレッドが求めた各priority毎のリスト長（interval数）
	int total_enum_data[max_priority + 1];				// 求めた各priority毎のデータ数の合計
	int final_priority = 0;
//	for(int t = 0; t < nt; t++) {
//		num_nonempty[t] = 0;
//	}
	for(int p = 0; p <= max_priority; p++) {
		total_enum_data[p] = 0;
		for(int t = 0; t < nt; t++) {
			num_data_of_priority[t][p] = lg_of_priority[t][p] = 0;
		}
	}
	#else
	int num_nonempty[nt * (1 << THREAD_PLUS)];								// スレッドが実際に列挙した空でないスケッチ数
	int num_data_of_priority[nt * (1 << THREAD_PLUS)][max_priority + 1];
	int lg_of_priority[nt * (1 << THREAD_PLUS)][max_priority + 1];
	int total_enum_data[max_priority + 1];						// 求めたデータ数の合計
	int final_priority = 0;
	for(int t = 0; t < nt * (1 << THREAD_PLUS); t++) {
		num_nonempty[t] = 0;
	}
	for(int p = 0; p <= max_priority; p++) {
		total_enum_data[p] = 0;
		for(int t = 0; t < nt * (1 << THREAD_PLUS); t++) {
			num_data_of_priority[t][p] = lg_of_priority[t][p] = 0;
		}
	}
	#endif

	int num_data_thread = num_candidates * FACTOR_INF / nt;
	#pragma omp parallel
	{
		// 下位 n ビットは，スレッド毎に割り当てた mu_thread[t] との XOR で求める．
		// それ以降の PJT_DIM - n ビットのみを変化させた，スケッチを列挙する（下位 n ビットは，質問のまま）
		int *bd = qs->bd, *bd_idx = qs->idx;
		int *bkt = bucket->bkt;
		sketch_type s, sk;
		int t = omp_get_thread_num();	// スレッド番号
		dist_type pr, pr_offset = pr_thread[t];
		QUE_c2 qu, qu2;
		struct_que_c2_n *que = que_pool[t];
		interval *buff = ivl->list + t * ivl->size;
		int ne = 0;		// 各スレッドで求めた空でないスケッチ数
		int nd = 0;		// 一回の処理で各スレッドで求めたデータ数
		int *nd_p = num_data_of_priority[t];
		int *lg_p = lg_of_priority[t];

		s = qs->sketch; // 先頭は質問のスケッチ
		sk = s ^ mu_thread[t];
		if(bkt[sk + 1] > bkt[sk]) {
			pr = pr_offset;
			buff[ne++] = (interval){pr, bkt[sk], bkt[sk + 1] - bkt[sk]};
			nd += bkt[sk + 1] - bkt[sk];
			nd_p[pr] += bkt[sk + 1] - bkt[sk];
			lg_p[pr]++;
		}

		#ifndef THREAD_PLUS
		int nn = n;
		#else
		int nn = n + THREAD_PLUS;
		#endif

		s = s ^ (1 <<  bd_idx[nn]); // 先頭の次は、質問のスケッチと距離下限が最小（ただし，nnビットは除く）のビットだけが異なるもの
		sk = s ^ mu_thread[t];
		if(bkt[sk + 1] > bkt[sk]) {
			pr = qs->bd[bd_idx[nn]] + pr_offset;
			buff[ne++] = (interval){pr, bkt[sk], bkt[sk + 1] - bkt[sk]};
			nd += bkt[sk + 1] - bkt[sk];
			nd_p[pr] += bkt[sk + 1] - bkt[sk];
			lg_p[pr]++;
		}

		make_empty_que_c2_n(que);

		// enq pattern of 0...10
		qu.cursor = new_que_e2_n(que);
		qu.key = bd[bd_idx[nn + 1]];
		que->details[qu.cursor].sk = qs->sketch ^ (1 << bd_idx[nn + 1]);
		que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
		enq_c2_n(&qu, que);		

		dist_type last_pr = 10000;
		while(deq_c2_n(&qu, que)) {
			s = que->details[qu.cursor].sk; // 列挙のつぎのスケッチ
			sk = s ^ mu_thread[t];
			pr = qu.key + pr_offset;
			if(pr > last_pr) break;
			if(bkt[sk + 1] > bkt[sk]) {
				buff[ne++] = (interval){pr, bkt[sk], bkt[sk + 1] - bkt[sk]};
				nd += bkt[sk + 1] - bkt[sk];
				nd_p[pr] += bkt[sk + 1] - bkt[sk];
				lg_p[pr]++;
				if(nd >= num_data_thread) {
					last_pr = pr;
				}
			}
			switch(que->details[qu.cursor].pt & 15) {
			case 0: // X0000 -> enq(X0001) and enq(Y10^{m+1}) if X0000 = Y010^m
			case 8: // X1000 -> enq(X1001) and enq(Y10^{m+1}) if X0000 = Y010^m
				{
					int m = lsb_pos(que->details[qu.cursor].pt);
					if(m > 0 && nn + m < PJT_DIM - 1 && !(que->details[qu.cursor].pt & (1 << (m + 1)))) {
						// Y010^m -> Y10^{m+1}
						qu2.cursor = new_que_e2_n(que);
						qu2.key = qu.key + bd[bd_idx[nn + m + 1]] - bd[bd_idx[nn + m]];
						que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + m + 1])) ^ (1 << bd_idx[nn + m]);
						que->details[qu2.cursor].pt = que->details[qu.cursor].pt + (1 << m);
						// Y010^m -> Y010^{m-1}1
						qu.key = qu.key + bd[bd_idx[nn + 0]];
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 0]);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
						enq_c2_n(&qu2, que);
					} else {
						qu.key = qu.key + bd[bd_idx[nn + 0]];
						que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 0]);
						que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
						enq_c2_n(&qu, que);
					}
				}
				break;
			case 4:  // X0100 -> enq(X0101) and enq(X1000)
				// X1000
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key + bd[bd_idx[nn + 3]] - bd[bd_idx[nn + 2]];
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 3])) ^ (1 << bd_idx[nn + 2]);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 4;
				// X0101
				qu.key = qu.key + bd[bd_idx[nn + 0]];
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 0]);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 1:  // X0001 -> enq(X0010)
			case 5:  // X0101 -> enq(X0110)
			case 9:  // X1001 -> enq(X1010)
			case 13: // X1101 -> enq(X1110) (note that X <> 0, because 0...00 and 0...01 is already processed before while loop)
				qu.key = qu.key + bd[bd_idx[nn + 1]] - bd[bd_idx[nn + 0]];
				que->details[qu.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 1])) ^ (1 << bd_idx[nn + 0]);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				break;
			case 2:  // X0010 -> enq(X0011) and enq(X0100)
			case 10: // X1010 -> enq(X1011) and enq(X1100)
				// X0100 and X1100
				qu2.cursor = new_que_e2_n(que);
				qu2.key = qu.key +  bd[bd_idx[nn + 2]] -  bd[bd_idx[nn + 1]];
				que->details[qu2.cursor].sk = (que->details[qu.cursor].sk ^ (1 << bd_idx[nn + 2])) ^ (1 << bd_idx[nn + 1]);
				que->details[qu2.cursor].pt = que->details[qu.cursor].pt + 2;
				// X0011 and X1011
				qu.key = qu.key + bd[bd_idx[nn + 0]];
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[nn + 0]);
				que->details[qu.cursor].pt = que->details[qu.cursor].pt + 1;
				enq_c2_n(&qu, que);
				enq_c2_n(&qu2, que);
				break;
			case 6:  // X0110 -> enq(X0111)
			case 12: // X1100 -> enq(X1101)
			case 14: // X1110 -> enq(10111)
				qu.key = qu.key + bd[bd_idx[nn + 0]];
				que->details[qu.cursor].sk = que->details[qu.cursor].sk ^ (1 <<  bd_idx[nn + 0]);
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
	}

	// max_priority を求め直す（実際にデータが存在する最大の priority，優先度としては最低）
	int sd;
	#ifndef THREAD_PLUS
	int NT = nt;
	#else
	int NT = nt * (1 << THREAD_PLUS);
	#endif
	while(1) {
		sd = 0;
		for(int t = 0; t < NT; t++) {
			sd += num_data_of_priority[t][max_priority];
		}
		if(sd != 0) break;
		max_priority--;
	}

	// num_data_of_priority[t][pr] を pr = 0, ... , pr の累積に変更する．
	// lg_of_priority[t][pr] も累積
	for(int t = 0; t < NT; t++) {
		for(int pr = 1; pr <= max_priority; pr++) {
			num_data_of_priority[t][pr] += num_data_of_priority[t][pr - 1];
			lg_of_priority[t][pr] += lg_of_priority[t][pr - 1];
		}
	}
	// 全合計が目標のnum_candidatesを超える最小の優先度 final_priority を求める．
	final_priority = 0;
	for(int pr = 0; pr <= max_priority; pr++) {
		total_enum_data[pr] = 0;
		for(int t = 0; t < NT; t++) {
			total_enum_data[pr] += num_data_of_priority[t][pr];
		}
		if(total_enum_data[pr] < num_candidates) final_priority = pr;
	}
	final_priority++;
	// 各スレッドで求めた final_priority までの区間数を ivl_lg[t] にセット．
	for(int t = 0; t < NT; t++) {
		ivl->lg[t] = lg_of_priority[t][final_priority];
	}
//			printf("total_enum_data = %d ===> ", total_enum_data[final_priority]);
	// final_priority の区間で余分なもの（データ数がnum_candidatesを超える部分）を削除
	for(int t = NT - 1; t >= 0; t--) {
		interval *buff = ivl->list + t * ivl->size;
		if(total_enum_data[final_priority] - (num_data_of_priority[t][final_priority] - num_data_of_priority[t][final_priority - 1]) >= num_candidates) {
			ivl->lg[t] = lg_of_priority[t][final_priority - 1];
			total_enum_data[final_priority] -= (num_data_of_priority[t][final_priority] - num_data_of_priority[t][final_priority - 1]);
		} else {
			int j = lg_of_priority[t][final_priority] - 1;
			while(buff[j].priority >= final_priority && total_enum_data[final_priority] - buff[j].run >= num_candidates) {
				total_enum_data[final_priority] -= buff[j--].run;
			}
			ivl->lg[t] = j;
		}
	}

	return total_enum_data[final_priority];
}
#endif

// score_1 (D~1) 順に列挙して，データ数が num_candidates を超える最初のスケッチの優先順位（score_1, D~1）を求める．
// 実際のデータ番号は求めない．
dist_type filtering_by_sketch_enumeration_c2_n_score(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, int num_candidates)
{
	sketch_type s;
	QUE_c2 qu, qu2;
	int *bd = qs->bd, *bd_idx = qs->idx;
	int *bkt = bucket->bkt;
	
	s = qs->sketch;
	int k = 0;
	for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates; j++, k++) {
	}
	if(k >= num_candidates) return 0;

	s = s ^ (1 <<  bd_idx[0]);
	for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates; j++, k++) {
	}
	if(k >= num_candidates) return bd[bd_idx[0]];

	make_empty_que_c2_n(que);

	// enq pattern of 0...10
	qu.cursor = new_que_e2_n(que);
	qu.key = bd[bd_idx[1]];
	que->details[qu.cursor].sk = qs->sketch ^ (1 << bd_idx[1]);
	que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
	enq_c2_n(&qu, que);		

	while(deq_c2_n(&qu, que) && k < num_candidates) {
		s = que->details[qu.cursor].sk;
		for(int j = bkt[s]; j < bkt[s + 1] && k < num_candidates; j++, k++) {
		}
		if(k >= num_candidates) return qu.key;

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
	return qu.key;
}

// D~p (p != inf) 順の列挙を用いたフィルタリング（空を含む）スケッチも含める
int filtering_by_sketch_enumeration_c2_n_sketch(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, sketch_type sketch[], int num_candidates)
{
	double ave_num = (double)bucket->num_data / (1L << PJT_DIM);
	int num_sketch = num_candidates / ave_num;
	sketch_type s;
	QUE_c2 qu, qu2;
	int *bd = qs->bd, *bd_idx = qs->idx;
	
	int num = 0;
	s = qs->sketch;
	sketch[num++] = s;

	s = s ^ (1 <<  bd_idx[0]);
	sketch[num++] = s;

	make_empty_que_c2_n(que);

	// enq pattern of 0...10
	qu.cursor = new_que_e2_n(que);
	qu.key = bd[bd_idx[1]];
	que->details[qu.cursor].sk = qs->sketch ^ (1 << bd_idx[1]);
	que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
	enq_c2_n(&qu, que);		

	while(deq_c2_n(&qu, que) && num < num_sketch) {
		s = que->details[qu.cursor].sk;
		sketch[num++] = s;

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
	return num;
}

#ifdef NEW_INF
#if PARA_ENUM_INF > 0
#else
	// D~inf順の列挙を用いる．single-thread after-selection
	int filtering_by_sketch_enumeration_inf(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates)
	{
		int *bd_idx = qs->idx;
		int *bkt = bucket->bkt;

		// 全体の D~inf 順のスケッチの列挙を a[0], a[1], ... とする．
		// a[0] = qs->sketch （質問のスケッチ，D~inf = 0）
		// a[i + 1] = a[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1])

		sketch_type s = qs->sketch; // 列挙するスケッチ
		int nc = 0;					// 求めた候補数

		for(int i = 0; nc < num_candidates; i++) { // このループでは，D~inf順にスケッチを列挙して，空でなければ，データ番号をdata_numに格納する．
			if(bkt[s + 1] > bkt[s]) { // s が空でない．
				for(int j = bkt[s]; j < bkt[s + 1] && nc < num_candidates; j++, nc++) {
					#ifdef DATA_NUM_IN_SKETCH_ORDER
					data_num[nc] = j;
					#else
					data_num[nc] = bucket->idx[j];
					#endif
				}
			}
			s = s ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1]); // 次のスケッチ
		}

		return nc;
	}

	// D~inf順の列挙を用いる．single-thread after-selection using interval
	int filtering_by_sketch_enumeration_inf_interval(struct_query_sketch *qs, struct_bucket *bucket, interval_list *ivl, int num_candidates)
	{
	//	fprintf(stderr, "enum_inf single-thread\n"); exit(0);
		int num_data_inf = num_candidates * FACTOR_INF; // D~inf順の列挙によって求める目標データ数（D~1で求める個数の1.5～3倍か？）
		int *bd_idx = qs->idx;
		int *bkt = bucket->bkt;

		#error "filtering_by_sketch_enumeration_inf_interval post-selection should be used with INTERVAL_WITH_PRIORITY"

		static int first = 1;
		if(first) {
			fprintf(stderr, "num_candidates = %d, num_data_inf = %d\n", num_candidates, num_data_inf);
			#ifndef WITHOUT_IDX
			fprintf(stderr, "enum_inf_interval single-thread WITH_IDX. INTERVAL_WITH_PRIORITY\n");
			#else
			fprintf(stderr, "enum_inf_interval single-thread WITHOUT_IDX. INTERVAL_WITH_PRIORITY\n");
			#endif
			first = 0;
		}


		// 全体の D~inf 順のスケッチの列挙を a[0], a[1], ... とする．
		// a[0] = qs->sketch （質問のスケッチ，D~inf = 0）
		// a[i + 1] = a[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1])

		sketch_type s = qs->sketch; // 列挙するスケッチ
		int num_nonempty = 0;		// 列挙した空でないスケッチ数（＝求めたinterval数）
		int nc = 0;					// 求めた候補数

		for(int i = 0; nc < num_data_inf; i++) { // このループでは，D~inf順にスケッチを列挙して，空でなければ，その区間（interval)を格納する．
			if(bkt[s + 1] > bkt[s]) { // s が空でない．
				buff[num_nonempty++] = (interval){priority(s, qs), bkt[s], bkt[s + 1] - bkt[s]};
				nc += (bkt[s + 1] - bkt[s]);
			}
			s = s ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1]); // 次のスケッチ
		}
		ivl->lg[0] = num_nonempty;
		ivl->lg[0] = quick_select_sum_k_interval(buff, 0, num_nonempty - 1, num_candidates);
		return num_candidates;
	}
#endif
#else // !defined(NEW_INF)
#if PARA_ENUM_INF > 0 
#define NO_EMPTY_CHECK
// D~inf 順に列挙した空でないスケッチで，データ数の合計が num_candidates 以上となるものを求める．
// 求めたスケッチ数を返す．求めたスケッチは sk[0], ... ．
// n は，マクロ PARA_ENUM_INF で定義しておくこと．
// スレッド数 2 ^ n (= 1, 2, 4, 8, 16, 32) の並列処理を行う．
int filtering_by_sketch_enumeration_inf(struct_query_sketch *qs, struct_bucket *bucket, sketch_type *sketch, int num_candidates)
{
	int n = PARA_ENUM_INF;
	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif
	double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
	int num_sketch_thread = num_candidates / ave_num / nt; // 各スレッドで求めるスケッチの個数（すべてのスレッドで同じ）
	int num_sketch = num_sketch_thread * nt; // 求めるスケッチの個数（合計）
	int num_sketch_inf = num_sketch * FACTOR_INF; // D~inf順の列挙によって求めるスケッチの個数（D~1で求める個数の3～4倍）
	static kNN_buffer *buff[1 << PARA_ENUM_INF] = {NULL};
	if(buff[0] == NULL) {
		fprintf(stderr, "new_kNN_buffer\n");
		for(int t = 0; t < nt; t++) {
			buff[t] = new_kNN_buffer(num_sketch_thread);
		}
	} else {
		for(int t = 0; t < nt; t++) {
			make_empty_kNN_buffer(buff[t]);
		}
	}

	unsigned bit = 1;
	bit <<= qs->idx[n - 1]; // 分割数に応じて固定で決定される操作ビット

	sketch_type s_0[nt]; // 各スレッドが列挙するスケッチの先頭
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	s_0[0] = qs->sketch;
	for(int i = 0; i < nt - 1; i++) {
		s_0[i + 1] = s_0[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1]);
	}

	int t; // スレッド番号
	sketch_type *p; // スレッド t が出力するスケッチの配列
	sketch_type s; // 列挙の先頭
	int r; // 列挙したスケッチの番号（先頭は0）
	kNN_buffer *b;

	#ifdef _OPENMP
	#pragma omp parallel private(t, p, s, r, b)
	#endif
	{
		#ifdef _OPENMP
		t = omp_get_thread_num(); // スレッド番号
		#else
		t = 0;
		#endif
		p = sketch + num_sketch_thread * t; // スレッド t が出力するスケッチの配列
		s = s_0[t]; // 列挙の先頭
		r = 0; // 列挙したスケッチの番号（先頭は0）
		b = buff[t];

		#ifdef _OPENMP
		#pragma omp for schedule(static, 1)
		#endif
		for(int i = 0; i < num_sketch_inf; i++) { // このループでは，D~inf順に空でないスケッチを1個求めて，kNN_bufferにpushする．
			int num_data_of_sk; // s は最後に列挙したスケッチ
			while((num_data_of_sk = bkt[s + 1] - bkt[s]) == 0) { // s が空ならば，空でないものが見つかるまでスケッチを列挙する．
				s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
				r++;
			}
			answer_type ans = (answer_type){(int)s, priority(s, qs)};
			push_kNN_buffer(&ans, b);
			// 次を列挙する．
			s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
			r++;
		}
		flush_kNN_buffer(b);
		for(int i = 0; i < num_sketch_thread; i++) {
			p[i] = (sketch_type)(b->buff[i].data_num);
		}
	}
	return num_sketch;
}

// D~inf 順に列挙した空でないスケッチで，データ数の合計が num_candidates 以上となるものを求める．
// 求めたスケッチ数を返す．求めたスケッチは sk[0], ... ．
// n は，マクロ PARA_ENUM_INF で定義しておくこと．
// スレッド数 2 ^ n (= 1, 2, 4, 8, 16, 32) の並列処理を行う．
int filtering_by_sketch_enumeration_inf_data(struct_query_sketch *qs, struct_bucket *bucket, int *data_num, int num_candidates)
{
	int n = PARA_ENUM_INF;
	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif
	double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
	int num_sketch_thread = num_candidates / ave_num / nt; // 各スレッドで求めるスケッチの個数（すべてのスレッドで同じ）
	int num_data_thread = num_candidates / nt; // 各スレッドで求めるデータの個数（すべてのスレッドで同じ）
	int num_sketch = num_sketch_thread * nt; // 求めるスケッチの個数（合計）
	int num_sketch_inf = num_sketch * FACTOR_INF; // D~inf順の列挙によって求めるスケッチの個数（D~1で求める個数の3～4倍）
	static kNN_buffer *buff[1 << PARA_ENUM_INF] = {NULL};
	if(buff[0] == NULL) {
		fprintf(stderr, "new_kNN_buffer\n");
		for(int t = 0; t < nt; t++) {
			buff[t] = new_kNN_buffer(num_sketch_thread);
		}
	} else {
		for(int t = 0; t < nt; t++) {
			make_empty_kNN_buffer(buff[t]);
		}
	}

	unsigned bit = 1;
	bit <<= qs->idx[n - 1]; // 分割数に応じて固定で決定される操作ビット

	sketch_type s_0[nt]; // 各スレッドが列挙するスケッチの先頭
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	s_0[0] = qs->sketch;
	for(int i = 0; i < nt - 1; i++) {
		s_0[i + 1] = s_0[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1]);
	}

	int t; // スレッド番号
	int *p; // スレッド t が出力するデータ番号の配列
	sketch_type s; // 列挙の先頭
	int r; // 列挙したスケッチの番号（先頭は0）
	int num; // 出力したデータ番号の個数
	kNN_buffer *b;

	#ifdef _OPENMP
	#pragma omp parallel private(t, p, s, r, num, b)
	#endif
	{
		#ifdef _OPENMP
		t = omp_get_thread_num(); // スレッド番号
		#else
		t = 0;
		#endif
		p = data_num + num_data_thread * t; // スレッド t が出力するデータ番号の配列
		s = s_0[t]; // 列挙の先頭
		r = 0; // 列挙したスケッチの番号（先頭は0）
		b = buff[t];

		#ifdef _OPENMP
		#pragma omp for schedule(static, 1)
		#endif
		for(int i = 0; i < num_sketch_inf; i++) { // このループでは，D~inf順に空でないスケッチを1個求めて，kNN_bufferにpushする．
			int num_data_of_sk; // s は最後に列挙したスケッチ
			while((num_data_of_sk = bkt[s + 1] - bkt[s]) == 0) { // s が空ならば，空でないものが見つかるまでスケッチを列挙する．
				s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
				r++;
			}
			answer_type ans = (answer_type){(int)s, priority(s, qs)};
			push_kNN_buffer(&ans, b);
			// 次を列挙する．
			s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
			r++;
		}
		flush_kNN_buffer(b);
		num = 0; // 出力したデータ番号の個数
		for(int i = 0; num < num_data_thread; i++) {
			sketch_type sk = (sketch_type)(b->buff[i].data_num);
			for(int j = bkt[sk]; j < bkt[sk + 1] && num < num_data_thread; j++) {
				p[num++] = j;
			}
		}
	}
	return num_data_thread * nt;
}

// D~inf 順に列挙した空でないスケッチで，データ数の合計が num_candidates 以上となるものを求める．(quick_select_k)
// 求めたスケッチ数を返す．求めたスケッチは sk[0], ... ．
// n は，マクロ PARA_ENUM_INF で定義しておくこと．
// スレッド数 2 ^ n (= 1, 2, 4, 8, 16, 32) の並列処理を行う．
int filtering_by_sketch_enumeration_inf_data_select(struct_query_sketch *qs, struct_bucket *bucket, int *data_num, int num_candidates)
{
	int n = PARA_ENUM_INF;
	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif
	double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
	int num_sketch_thread = num_candidates / ave_num / nt; // 各スレッドで求めるスケッチの個数（すべてのスレッドで同じ）
	int num_sketch_thread_inf = num_sketch_thread * FACTOR_INF; // 各スレッドでD~inf順の列挙によって求めるスケッチの個数（すべてのスレッドで同じ）
	int num_data_thread = num_candidates / nt; // 各スレッドで求めるデータの個数（すべてのスレッドで同じ）
	int num_sketch = num_sketch_thread * nt; // 求めるスケッチの個数（合計）
	int num_sketch_inf = num_sketch * FACTOR_INF; // D~inf順の列挙によって求めるスケッチの個数（D~1で求める個数の3～4倍）
	static answer_type *buff[1 << PARA_ENUM_INF] = {NULL};
	if(buff[0] == NULL) {
		fprintf(stderr, "malloc answer buffer\n");
		for(int t = 0; t < (1<< PARA_ENUM_INF); t++) {
			buff[t] = MALLOC(sizeof(answer_type) * num_sketch_thread_inf);
		}
		fprintf(stderr, "malloc answer buffer OK\n");
	}

	unsigned bit = 1;
	bit <<= qs->idx[n - 1]; // 分割数に応じて固定で決定される操作ビット

	sketch_type s_0[nt]; // 各スレッドが列挙するスケッチの先頭
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	s_0[0] = qs->sketch;
	for(int i = 0; i < nt - 1; i++) {
		s_0[i + 1] = s_0[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1]);
	}

	int t; // スレッド番号
	int *p; // スレッド t が出力するデータ番号の配列
	sketch_type s; // 列挙の先頭
	int r; // 列挙したスケッチの番号（先頭は0）
	int r_t; // スレッドで列挙したスケッチの個数
	int num; // 出力したデータ番号の個数
	answer_type *b;

	#ifdef _OPENMP
	#pragma omp parallel private(t, p, s, r, r_t, num, b)
	#endif
	{
		#ifdef _OPENMP
		t = omp_get_thread_num(); // スレッド番号
		#else
		t = 0;
		#endif
		p = data_num + num_data_thread * t; // スレッド t が出力するデータ番号の配列
		s = s_0[t]; // 列挙の先頭
		r = 0; // 列挙したスケッチの番号（先頭は0）
		r_t = 0;
		b = buff[t];

		#ifdef _OPENMP
		#pragma omp for schedule(static, 1)
		#endif
		for(int i = 0; i < num_sketch_inf; i++) { // このループでは，D~inf順に空でないスケッチを1個求めて，kNN_bufferにpushする．
			int num_data_of_sk; // s は最後に列挙したスケッチ
			while((num_data_of_sk = bkt[s + 1] - bkt[s]) == 0) { // s が空ならば，空でないものが見つかるまでスケッチを列挙する．
				s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
				r++;
			}
			answer_type ans = (answer_type){(int)s, priority(s, qs)};
			b[r_t++] = ans;
			// 次を列挙する．
			s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
			r++;
		}
		quick_select_k_answer(b, 0, num_sketch_thread_inf - 1, num_sketch_thread);
		num = 0; // 出力したデータ番号の個数
		for(int i = 0; num < num_data_thread; i++) {
			sketch_type sk = (sketch_type)(b[i].data_num);
			for(int j = bkt[sk]; j < bkt[sk + 1] && num < num_data_thread; j++) {
				#ifdef DATA_NUM_IN_SKETCH_ORDER
				p[num++] = j;
				#else
				p[num++] = bucket->idx[j];
				#endif
			}
		}
	}
	return num_data_thread * nt;
}

// D~inf 順に列挙した空でないスケッチで，データ数の合計が num_candidates 以上となるものを求める．
// D~1による選択は行わない
// 求めたスケッチ数を返す．求めたスケッチは sk[0], ... ．
// n は，マクロ PARA_ENUM_INF で定義しておくこと．
// スレッド数 2 ^ n (= 1, 2, 4, 8, 16, 32) の並列処理を行う．
int filtering_by_sketch_enumeration_inf_only(struct_query_sketch *qs, struct_bucket *bucket, int *data_num, int num_candidates)
{
	int n = PARA_ENUM_INF;
	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif
	int num_data_thread = num_candidates / nt; // 各スレッドで求めるデータの個数（すべてのスレッドで同じ．端数が出るときは合計はnum_candidatesに達しない）
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	unsigned bit = 1;
	bit <<= bd_idx[n - 1]; // 分割数に応じて固定で決定される操作ビット

	sketch_type s_0[nt]; // 各スレッドが列挙するスケッチの先頭
	s_0[0] = qs->sketch;
	for(int i = 0; i < nt - 1; i++) {
		s_0[i + 1] = s_0[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1]);
	}

	int t; // スレッド番号
	int *p; // スレッド t が出力するデータ番号の配列
	sketch_type s; // 列挙の先頭
	int r; // 列挙したスケッチの番号（先頭は0）
	int num; // 出力したデータ番号の個数

	#ifdef _OPENMP
	#pragma omp parallel private(t, p, s, r, num)
	#endif
	{
		#ifdef _OPENMP
		t = omp_get_thread_num(); // スレッド番号
		#else
		t = 0;
		#endif
		p = data_num + num_data_thread * t; // スレッド t が出力するデータ番号の配列
		s = s_0[t]; // 列挙の先頭
		r = 0; // 列挙したスケッチのスレッド内の番号（先頭は0）
		num = 0;

		while(num < num_data_thread) {
			while(bkt[s + 1] == bkt[s]) { // s が空ならば，空でないものが見つかるまでスケッチを列挙する．
				s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
				r++;
			}
			for(int j = bkt[s]; j < bkt[s + 1] && num < num_data_thread; j++) {
				#ifdef DATA_NUM_IN_SKETCH_ORDER
					p[num++] = j;
				#else
					p[num++] = bucket->idx[j];
				#endif
			}
			// 次を列挙する．
			s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
			r++;
		}
	}
	return num_data_thread * nt;
}

int filtering_by_sketch_enumeration_inf_only_2(struct_query_sketch *qs, struct_bucket *bucket, int *data_num_thread[], int num_data_thread[], int num_candidates)
{
	int n = PARA_ENUM_INF;
	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif
	double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
	int num_sketch_thread = num_candidates / ave_num / nt; // 各スレッドで求めるスケッチの個数（すべてのスレッドで同じ）
	int nc_thread = num_candidates / nt; // 各スレッドで求めるデータの個数（すべてのスレッドで同じ）
	int num_sketch = num_sketch_thread * nt; // 求めるスケッチの個数（合計）
	static sketch_type *buff[1 << PARA_ENUM_INF] = {NULL};
	if(buff[0] == NULL) {
		for(int t = 0; t < (1<< PARA_ENUM_INF); t++) {
			buff[t] = MALLOC(sizeof(sketch_type) * num_sketch_thread);
		}
	}

	unsigned bit = 1;
	bit <<= qs->idx[n - 1]; // 分割数に応じて固定で決定される操作ビット

	sketch_type s_0[nt]; // 各スレッドが列挙するスケッチの先頭
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	s_0[0] = qs->sketch;
	for(int i = 0; i < nt - 1; i++) {
		s_0[i + 1] = s_0[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1]);
	}

	int t; // スレッド番号
	int *p; // スレッド t が出力するデータ番号の配列
	sketch_type s; // 列挙の先頭
	int r; // 列挙したスケッチの番号（先頭は0）
	int r_t; // スレッドで列挙したスケッチの個数
	int num; // 出力したデータ番号の個数
	sketch_type *b;

	#ifdef _OPENMP
	#pragma omp parallel private(t, p, s, r, r_t, num, b)
	#endif
	{
		#ifdef _OPENMP
		t = omp_get_thread_num(); // スレッド番号
		#else
		t = 0;
		#endif
		p = data_num_thread[t]; // スレッド t が出力するデータ番号の配列
		s = s_0[t]; // 列挙の先頭
		r = 0; // 列挙したスケッチの番号（先頭は0）
		r_t = 0;
		b = buff[t];

		#ifdef _OPENMP
		#pragma omp for schedule(static, 1)
		#endif
		for(int i = 0; i < num_sketch; i++) { // このループでは，D~inf順に空でないスケッチを1個求めて，bufferに入れる．
			int num_data_of_sk; // s は最後に列挙したスケッチ
			while((num_data_of_sk = bkt[s + 1] - bkt[s]) == 0) { // s が空ならば，空でないものが見つかるまでスケッチを列挙する．
				s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
				r++;
			}
			b[r_t++] = s;
			// 次を列挙する．
			s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
			r++;
		}

		num = 0; // 出力したデータ番号の個数
		for(int i = 0; i < num_sketch_thread; i++) {
			sketch_type sk = b[i];
			for(int j = bkt[sk]; j < bkt[sk + 1] && num < nc_thread; j++) {
				p[num++] = j;
			}
		}
		num_data_thread[t] = num;
	}

	num = 0;
	for(t = 0; t < nt; t++) {
		num += num_data_thread[t];
	}
	
	return num;
}

int filtering_by_sketch_enumeration_inf_only_3(struct_query_sketch *qs, struct_bucket *bucket, int *data_num_thread[], int num_data_thread[], int num_candidates)
{
	int n = PARA_ENUM_INF;
	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif
	double ave_num = (double)bucket->num_data / (1L << PJT_DIM); // 空も含むすべてのバケットの平均要素数
	int num_sketch_thread = num_candidates / ave_num / nt; // 各スレッドで求めるスケッチの個数（すべてのスレッドで同じ）
	int nc_thread = num_candidates / nt; // 各スレッドで求めるデータの個数（すべてのスレッドで同じ）
	int num_sketch = num_sketch_thread * nt; // 求めるスケッチの個数（合計）
	static sketch_type *buff[1 << PARA_ENUM_INF] = {NULL};
	if(buff[0] == NULL) {
		for(int t = 0; t < (1<< PARA_ENUM_INF); t++) {
			buff[t] = MALLOC(sizeof(sketch_type) * num_sketch_thread);
		}
	}

	unsigned bit = 1;
	bit <<= qs->idx[n - 1]; // 分割数に応じて固定で決定される操作ビット

	sketch_type s_0[nt]; // 各スレッドが列挙するスケッチの先頭
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	s_0[0] = qs->sketch;
	for(int i = 0; i < nt - 1; i++) {
		s_0[i + 1] = s_0[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1]);
	}

	int t; // スレッド番号
	int *p; // スレッド t が出力するデータ番号の配列
	sketch_type s; // 列挙の先頭
	int r; // 列挙したスケッチの番号（先頭は0）
	int r_t; // スレッドで列挙したスケッチの個数
	int num; // 出力したデータ番号の個数
	sketch_type *b;

	#ifdef _OPENMP
	#pragma omp parallel private(t, p, s, r, r_t, num, b)
	#endif
	{
		#ifdef _OPENMP
		t = omp_get_thread_num(); // スレッド番号
		#else
		t = 0;
		#endif
		p = data_num_thread[t]; // スレッド t が出力するデータ番号の配列
		s = s_0[t]; // 列挙の先頭
		r = 0; // 列挙したスケッチの番号（先頭は0）
		r_t = 0;
		b = buff[t];

		#ifdef _OPENMP
		#pragma omp for schedule(static, 1)
		#endif
		for(int i = 0; i < num_sketch; i++) { // このループでは，D~inf順にスケッチを1個求めて，bufferに入れる．
			b[r_t++] = s;
			// 次を列挙する．
			s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
			r++;
		}

		num = 0; // 出力したデータ番号の個数
		int ii;
		for(ii = 0; ii < num_sketch_thread && num < nc_thread; ii++) {
			sketch_type sk = b[ii];
			for(int j = bkt[sk]; j < bkt[sk + 1] && num < nc_thread; j++) {
				p[num++] = j;
			}
		}
		num_data_thread[t] = num;
	}

	num = 0;
	for(t = 0; t < nt; t++) {
		num += num_data_thread[t];
	}
	return num;
}

int filtering_by_sketch_enumeration_inf_only_4(struct_query_sketch *qs, struct_bucket *bucket, int *data_num_thread[], int num_data_thread[], int num_candidates)
{
	int n = PARA_ENUM_INF;
	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif
	double ave_num = (double)bucket->num_data / (1L << PJT_DIM); // 空も含むすべてのバケットの平均要素数
	int num_sketch_thread = num_candidates / ave_num / nt; // 各スレッドで求めるスケッチの個数（すべてのスレッドで同じ）
	int nc_thread = num_candidates / nt; // 各スレッドで求めるデータの個数（すべてのスレッドで同じ）
	int num_sketch = num_sketch_thread * nt; // 求めるスケッチの個数（合計）
	static sketch_type *buff = NULL;
	if(buff == NULL) {
		buff = MALLOC(sizeof(sketch_type) * num_sketch);
	}

	unsigned bit = 1;
	bit <<= qs->idx[n - 1]; // 分割数に応じて固定で決定される操作ビット

	sketch_type s_0[nt]; // 各スレッドが列挙するスケッチの先頭
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	s_0[0] = qs->sketch;
	for(int i = 0; i < nt - 1; i++) {
		s_0[i + 1] = s_0[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1]);
	}

	int t; // スレッド番号
	int *p; // スレッド t が出力するデータ番号の配列
	sketch_type s; // 列挙の先頭
	int r; // 列挙したスケッチの番号（先頭は0）
	int r_t; // スレッドで列挙したスケッチの個数
	int num; // 出力したデータ番号の個数

	#ifdef _OPENMP
	#pragma omp parallel private(t, s, r_t, r)
	#endif
	{
		#ifdef _OPENMP
		t = omp_get_thread_num(); // スレッド番号
		#else
		t = 0;
		#endif
		s = s_0[t]; // 列挙の先頭
		r = 0; // 列挙したスケッチの番号（先頭は0）
		r_t = 0;

		#ifdef _OPENMP
		#pragma omp for schedule(static, 1)
		#endif
		for(int i = 0; i < num_sketch; i++) { // このループでは，D~inf順にスケッチを1個求めて，bufferに入れる．
			buff[r_t++ + t * num_sketch_thread] = s;
			// 次を列挙する．
			r++;
			s = s ^ (1 << bd_idx[bit_count((r - 1) ^ (r)) - 1 + n]);
		}
	}
	// データ番号を展開
	#ifdef _OPENMP
	#pragma omp parallel private(t, p, num)
	#endif
	{
		#ifdef _OPENMP
		t = omp_get_thread_num(); // スレッド番号
		#else
		t = 0;
		#endif
		p = data_num_thread[t]; // スレッド t が出力するデータ番号の配列
		num = 0; // 出力したデータ番号の個数
		int ii;
		for(ii = 0; ii < num_sketch_thread && num < nc_thread; ii++) {
			sketch_type sk = buff[ii + t * num_sketch_thread];
			for(int j = bkt[sk]; j < bkt[sk + 1] && num < nc_thread; j++) {
				p[num++] = j;
			}
		}
		num_data_thread[t] = num;
	}

	num = 0;
	for(t = 0; t < nt; t++) {
		num += num_data_thread[t];
	}
	
	return num;
}

int filtering_by_sketch_enumeration_inf_only_sketch(struct_query_sketch *qs, struct_bucket *bucket, sketch_type sketch[], int num_candidates)
{
	int n = PARA_ENUM_INF;
	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif
	double ave_num = (double)bucket->num_data / (1L << PJT_DIM); // 空も含むすべてのバケットの平均要素数
	int num_sketch_thread = num_candidates / ave_num / nt; // 各スレッドで求めるスケッチの個数（すべてのスレッドで同じ）
	int num_sketch = num_sketch_thread * nt; // 求めるスケッチの個数（合計）

	static int first = 1;
	if(first) {
		fprintf(stderr, "num_candidates = %d, num_sketch = %d, num_sketch_thread = %d\n", num_candidates, num_sketch, num_sketch_thread);
		first = 0;
	}

	unsigned bit = 1;
	bit <<= qs->idx[n - 1]; // 分割数に応じて固定で決定される操作ビット

	sketch_type s_0[nt]; // 各スレッドが列挙するスケッチの先頭
	int *bd_idx = qs->idx;

	s_0[0] = qs->sketch;
	for(int i = 0; i < nt - 1; i++) {
		s_0[i + 1] = s_0[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1]);
	}

	int t; // スレッド番号
	sketch_type s; // 列挙の先頭
	int r_t; // スレッドで列挙したスケッチの個数
	int r; // 列挙したスケッチの番号（先頭は0）
	sketch_type *b;

	#ifdef _OPENMP
	#pragma omp parallel private(t, s, r, r_t, b)
	#endif
	{
		#ifdef _OPENMP
		t = omp_get_thread_num(); // スレッド番号
		#else
		t = 0;
		#endif
		s = s_0[t]; // 列挙の先頭
		r = 0; // 列挙したスケッチの番号（先頭は0）
		r_t = 0;
		b = sketch + (t * num_sketch_thread);

		#ifdef _OPENMP
		#pragma omp for schedule(static, 1)
		#endif
		for(int i = 0; i < num_sketch; i++) { // このループでは，D~inf順にスケッチを1個求めて，bufferに入れる．
			b[r_t++] = s;
			// 次を列挙する．
			s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
			r++;
		}
	}
	return num_sketch;
}

// D~inf 順に列挙した空でないスケッチで，データ数の合計が num_candidates 以上となるものを求める．(quick_select_k は最後に一回だけ)
// 求めたスケッチ数を返す．求めたスケッチは sk[0], ... ．
// n は，マクロ PARA_ENUM_INF で定義しておくこと．
// スレッド数 2 ^ n (= 1, 2, 4, 8, 16, 32) の並列処理を行う．
int filtering_by_sketch_enumeration_inf_data_select_once(struct_query_sketch *qs, struct_bucket *bucket, int *data_num, int num_candidates)
{
	int n = PARA_ENUM_INF;
	int nt = (1 << n); // スレッド数
	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif
	double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
	int num_sketch = num_candidates / ave_num / nt * nt; // 求めるD~1が上位のスケッチの個数（合計）（スレッド数の倍数になるようにしている）
	int num_sketch_inf = num_sketch * FACTOR_INF; // D~inf順の列挙によって求めるスケッチの合計個数（D~1で求める個数の3～4倍）
	int num_sketch_thread_inf = num_sketch_inf / nt; // 各スレッドでD~inf順の列挙によって求めるスケッチの個数（すべてのスレッドで同じ）
	int num_data_thread = num_candidates / nt; // 各スレッドで求めるデータの個数（すべてのスレッドで同じ）
	static answer_type *buff = NULL;
	if(buff == NULL) {
		fprintf(stderr, "malloc answer buffer\n");
		buff = MALLOC(sizeof(answer_type) * num_sketch_inf);
		fprintf(stderr, "malloc answer buffer OK\n");
	}

	unsigned bit = 1;
	bit <<= qs->idx[n - 1]; // 分割数に応じて固定で決定される操作ビット

	sketch_type s_0[nt]; // 各スレッドが列挙するスケッチの先頭
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	s_0[0] = qs->sketch;
	for(int i = 0; i < nt - 1; i++) {
		s_0[i + 1] = s_0[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1]);
	}

	int t; // スレッド番号
	sketch_type s; // 列挙の先頭
	int r; // 列挙したスケッチの番号（先頭は0）
	int r_t; // スレッドで列挙したスケッチの個数
	int num; // 出力したデータ番号の個数
	answer_type *b;

	#ifdef _OPENMP
	#pragma omp parallel private(t, s, r, r_t, b)
	#endif
	{
		#ifdef _OPENMP
		t = omp_get_thread_num(); // スレッド番号
		#else
		t = 0;
		#endif
		s = s_0[t]; // 列挙の先頭
		r = 0; // スレッドが列挙したスケッチの番号（先頭は0）
		r_t = 0; // スレッドがバッファに書き込んだスケッチ数
		b = buff + t * num_sketch_thread_inf; // スレッドがスケッチを書き込むバッファ（列挙は一つおきに行うが，書き込みは連続する）

		#ifdef _OPENMP
		#pragma omp for schedule(static, 1)
		#endif
		for(int i = 0; i < num_sketch_inf; i++) { // このループでは，D~inf順に空でないスケッチを1個求めて，buffに入れる．
			int num_data_of_sk; // s は最後に列挙したスケッチ
			while((num_data_of_sk = bkt[s + 1] - bkt[s]) == 0) { // s が空ならば，空でないものが見つかるまでスケッチを列挙する．
				s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
				r++;
			}
			answer_type ans = (answer_type){(int)s, priority(s, qs)};
			b[r_t++] = ans;
			// 次を列挙する．
			s = s ^ bit ^ (1 << bd_idx[bit_count((r * nt + nt - 1) ^ (r * nt + nt)) - 1]);
			r++;
		}
	}

	// 一旦，並列ループを終了して，一度だけで quick_select_k で buff 全体から D~1 の上位のものを選ぶ．
	quick_select_k_answer(buff, 0, num_sketch_inf - 1, num_sketch);

	// スケッチをデータ番号に展開する
	int *p; // スレッドがデータを書き込むときの配列
	int remain;

	#ifdef _OPENMP
	#pragma omp parallel private(t, p, num, b, remain)
	#endif
	{
		#ifdef _OPENMP
		t = omp_get_thread_num(); // スレッド番号
		#else
		t = 0;
		#endif
		p = data_num + t * num_data_thread; // スレッド t が出力するデータ番号の配列
		num = 0;
		#ifdef _OPENMP
		#pragma omp for
		#endif
		for(int i = 0; i < num_sketch; i++) {
			sketch_type sk = (sketch_type)(buff[i].data_num);
			for(int j = bkt[sk]; j < bkt[sk + 1] && num < num_data_thread; j++) {
				p[num++] = j;
			}
		}
		// データ数が目標に到達していないときは，
		if(num < num_data_thread) {
			b = buff + num_sketch + t * num_sketch / nt;
			remain = num_data_thread - num;
			for(int i = 0; i < remain; i++) {
				sketch_type sk = (sketch_type)(b[i].data_num);
				for(int j = bkt[sk]; j < bkt[sk + 1] && num < num_data_thread; j++) {
					p[num++] = j;
				}
			}
		}
	}
	return num_data_thread * nt;
}

#else // PARA_ENUM_INF == 0
// D~inf順の列挙を用いる．シングルスレッド版
int filtering_by_sketch_enumeration_inf(struct_query_sketch *qs, struct_bucket *bucket, sketch_type sketch[], int num_candidates)
{
//	fprintf(stderr, "enum_inf single-thread\n"); exit(0);
	double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
	int num_sketch = num_candidates / ave_num; // 求めるスケッチの個数
	kNN_buffer *buff = new_kNN_buffer(num_sketch);
	int num_sketch_inf = num_sketch * FACTOR_INF; // D~inf順の列挙によって求めるスケッチの個数（D~1で求める個数の3～4倍）
	sketch_type s; // 列挙するスケッチ
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	// 全体の D~inf 順のスケッチの列挙を a[0], a[1], ... とする．
	// a[0] = qs->sketch （質問のスケッチ，D~inf = 0）
	// a[i + 1] = a[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1])

	s = qs->sketch;
	int r = 0; // 列挙したスケッチの番号（先頭は0）

	for(int i = 0; i < num_sketch_inf; i++) { // このループでは，D~inf順に空でないスケッチを1個求めて，kNN_bufferにpushする．
		int num_data_of_sk; // s は最後に列挙したスケッチ
		while((num_data_of_sk = bkt[s + 1] - bkt[s]) == 0) { // s が空ならば，空でないものが見つかるまでスケッチを列挙する．
			s = s ^ (1 << bd_idx[bit_count(r ^ (r + 1)) - 1]);
			r++;
		}
		answer_type ans = (answer_type){(int)s, priority(s, qs)};
		push_kNN_buffer(&ans, buff);
		s = s ^ (1 << bd_idx[bit_count(r ^ (r + 1)) - 1]); // つぎを列挙
		r++;
	}
	flush_kNN_buffer(buff);
	for(int i = 0; i < num_sketch; i++) {
		sketch[i] = (sketch_type)(buff->buff[i].data_num);
	}
	free_kNN_buffer(buff);

	return num_sketch;
}

// D~inf順の列挙を用いる．シングルスレッド版(data_numのリストを返す)
int filtering_by_sketch_enumeration_inf_data(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates)
{
//	fprintf(stderr, "enum_inf single-thread\n"); exit(0);
	double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
	int num_sketch = num_candidates / ave_num; // 求めるスケッチの個数
	kNN_buffer *buff = new_kNN_buffer(num_sketch);
	int num_sketch_inf = num_sketch * FACTOR_INF; // D~inf順の列挙によって求めるスケッチの個数（D~1で求める個数の3～4倍）
	sketch_type s; // 列挙するスケッチ
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	// 全体の D~inf 順のスケッチの列挙を a[0], a[1], ... とする．
	// a[0] = qs->sketch （質問のスケッチ，D~inf = 0）
	// a[i + 1] = a[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1])

	s = qs->sketch;
	int r = 0; // 列挙したスケッチの番号（先頭は0）

	for(int i = 0; i < num_sketch_inf; i++) { // このループでは，D~inf順に空でないスケッチを1個求めて，kNN_bufferにpushする．
		int num_data_of_sk; // s は最後に列挙したスケッチ
		while((num_data_of_sk = bkt[s + 1] - bkt[s]) == 0) { // s が空ならば，空でないものが見つかるまでスケッチを列挙する．
			s = s ^ (1 << bd_idx[bit_count(r ^ (r + 1)) - 1]);
			r++;
		}
		answer_type ans = (answer_type){(int)s, priority(s, qs)};
		push_kNN_buffer(&ans, buff);
		s = s ^ (1 << bd_idx[bit_count(r ^ (r + 1)) - 1]); // つぎを列挙
		r++;
	}
	flush_kNN_buffer(buff);

	int num_data = 0;
	for(int i = 0; i < num_sketch && num_data < num_candidates; i++) {
		sketch_type sk = (sketch_type)(buff->buff[i].data_num);
		for(int j = bkt[sk]; j < bkt[sk + 1] && num_data < num_candidates; j++) {
			data_num[num_data++] = j;
		}
	}
	free_kNN_buffer(buff);

	return num_data;
}

// D~inf順の列挙を用いる．シングルスレッド版(data_numのリストを返す)(quick_select_k)
int filtering_by_sketch_enumeration_inf_data_select(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates)
{
	double ave_num = (double)bucket->num_data / bucket->num_nonempty_buckets; // 空でないバケットの平均要素数
	int num_sketch = num_candidates / ave_num; // 求めるスケッチの個数
	static answer_type *buff = NULL;
	int num_sketch_inf = num_sketch * FACTOR_INF; // D~inf順の列挙によって求めるスケッチの個数（D~1で求める個数の3～4倍）
	if(buff == NULL) { buff = MALLOC(sizeof(answer_type) * num_sketch_inf); }
	sketch_type s; // 列挙するスケッチ
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	// 全体の D~inf 順のスケッチの列挙を a[0], a[1], ... とする．
	// a[0] = qs->sketch （質問のスケッチ，D~inf = 0）
	// a[i + 1] = a[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1])

	s = qs->sketch;
	int r = 0; // 列挙したスケッチの番号（先頭は0）

	for(int i = 0; i < num_sketch_inf; i++) { // このループでは，D~inf順に空でないスケッチを1個求めて，kNN_bufferにpushする．
		int num_data_of_sk; // s は最後に列挙したスケッチ
		while((num_data_of_sk = bkt[s + 1] - bkt[s]) == 0) { // s が空ならば，空でないものが見つかるまでスケッチを列挙する．
			s = s ^ (1 << bd_idx[bit_count(r ^ (r + 1)) - 1]);
			r++;
		}
		answer_type ans = (answer_type){(int)s, priority(s, qs)};
		buff[i] = ans;
		s = s ^ (1 << bd_idx[bit_count(r ^ (r + 1)) - 1]); // つぎを列挙
		r++;
	}
	quick_select_k_answer(buff, 0, num_sketch_inf - 1, num_sketch);

	int num_data = 0;
	for(int i = 0; i < num_sketch && num_data < num_candidates; i++) {
		sketch_type sk = (sketch_type)(buff[i].data_num);
		for(int j = bkt[sk]; j < bkt[sk + 1] && num_data < num_candidates; j++) {
			#ifdef DATA_NUM_IN_SKETCH_ORDER
			data_num[num_data++] = j;
			#else
			data_num[num_data++] = bucket->idx[j];
			#endif
		}
	}

	return num_data;
}

// D~inf順の列挙を用いる．D~1での選択は行わないシングルスレッド版(data_numのリストを返す)
int filtering_by_sketch_enumeration_inf_only(struct_query_sketch *qs, struct_bucket *bucket, int data_num[], int num_candidates)
{
	sketch_type s; // 列挙するスケッチ
	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	// 全体の D~inf 順のスケッチの列挙を a[0], a[1], ... とする．
	// a[0] = qs->sketch （質問のスケッチ，D~inf = 0）
	// a[i + 1] = a[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1])

	s = qs->sketch;
	int r = 0; // 列挙したスケッチの番号（先頭は0）
	int num_data = 0;
	while(num_data < num_candidates) { // このループでは，D~inf順に空でないスケッチを1個求めて，data_numに展開する．
		// s が空ならば，空でないものが見つかるまでスケッチを列挙する．
		int num_data_sketch = bkt[s + 1] - bkt[s];
		while(num_data_sketch == 0) {
			s = s ^ (1 << bd_idx[bit_count(r ^ (r + 1)) - 1]);
			r++;
			num_data_sketch = bkt[s + 1] - bkt[s];
		}
		// 空でないスケッチが見つかったので，データ番号を格納する
		int m = num_data_sketch < (num_candidates - num_data) ? num_data_sketch : num_candidates - num_data;
		for(int j = 0; j < m; j++) {
			#ifdef DATA_NUM_IN_SKETCH_ORDER
				data_num[num_data++] = bkt[s] + j;
			#else
				data_num[num_data++] = bucket->idx[bkt[s] + j];
			#endif
		}
		s = s ^ (1 << bd_idx[bit_count(r ^ (r + 1)) - 1]); // つぎを列挙
		r++;
	}
	return num_data;
}

// D~inf順の列挙を用いる．空のスケッチの除去は行わないで，列挙したものを全部返す．D~1での選択は行わないシングルスレッド版(sketchのリストを返す)
int filtering_by_sketch_enumeration_inf_only_sketch(struct_query_sketch *qs, struct_bucket *bucket, sketch_type sketch[], int num_candidates)
{
//	fprintf(stderr, "enum_inf single-thread-inf_only_sketch. FACTOR_INF = %d, PARA_ENUM_INF = %d\n", FACTOR_INF, PARA_ENUM_INF); exit(0);
	double ave_num = (double)bucket->num_data / (1L << PJT_DIM); // バケットの平均要素数（スケッチがもつデータ数の平均）
	int num_sketch = num_candidates / ave_num; // 列挙するスケッチ数
	int *bd_idx = qs->idx;

	static int first = 1;
	if(first) {
		fprintf(stderr, "num_candidates = %d, num_sketch = %d\n", num_candidates, num_sketch);
		first = 0;
	}

	// 全体の D~inf 順のスケッチの列挙を a[0], a[1], ... とする．
	// a[0] = qs->sketch （質問のスケッチ，D~inf = 0）
	// a[i + 1] = a[i] ^ (1 << bd_idx[bit_count(i ^ (i + 1)) - 1])

	sketch[0] = qs->sketch;
	for(int r = 1; r < num_sketch; r++) { // このループでは，D~inf順にスケッチを1個求めて，sketchに格納する．
		sketch[r] = sketch[r - 1] ^ (1 << bd_idx[bit_count(r ^ (r + 1)) - 1]);
	}

	return num_sketch;
}
#endif // PARA_ENUM_INF
#endif // NEW_INF
#endif

#endif
// ! WITHOUT_FILTERING

#ifdef CHECK_ENUMERATION
void check_sketch_enumeration_c2_n(struct_query_sketch *qs, struct_bucket *bucket, struct_que_c2_n *que, struct_check_result *result, int num_candidates)
{
	sketch_type s;
	QUE_c2 qu, qu2;
	int *bd = qs->bd, *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	result->num_enum_sketches = result->num_enum_nonempty_sketches = result->num_enum_data = -1;	// 列挙する候補数が num_candidates 以下では正解が得られないとき．
	for(int n = 0; n < 10; n++) result->num_data[n] = 0;
	s = qs->sketch; // 先頭は質問のスケッチ
	int k = 0;		// 列挙したデータ数
	int num_sk = 0; // 列挙したスケッチ数 - 1
	int num_ne = 0; // 列挙した空でないスケッチ数
	if(bkt[s] < bkt[s + 1]) { // 空でないバケツの範囲に正解があるかをチェックして，見つかれば，列挙したスケッチ数とデータ数をresultに記録して復帰
		num_ne++;
		k += bkt[s + 1] - bkt[s];
		if(qs->answer.data_num >= bkt[s] && qs->answer.data_num < bkt[s + 1]) {
			result->num_enum_sketches = num_sk;
			result->num_enum_nonempty_sketches = num_ne;
			result->num_enum_data = k;
			result->num_data[num_sk++] = k;
			return;
		}
	}
	result->num_data[num_sk++] = k;
	s = s ^ (1 <<  bd_idx[0]); // 先頭の次は、質問のスケッチと距離下限が最小のビットだけが異なるもの
	if(bkt[s] < bkt[s + 1]) { // 空でないバケツの範囲に正解があるかをチェックして，見つかれば，列挙したスケッチ数とデータ数をresultに記録して復帰
		num_ne++;
		k += bkt[s + 1] - bkt[s];
		if(qs->answer.data_num >= bkt[s] && qs->answer.data_num < bkt[s + 1]) {
			result->num_enum_sketches = num_sk;
			result->num_enum_nonempty_sketches = num_ne;
			result->num_enum_data = k;
			result->num_data[num_sk++] = k;
			return;
		}
	}
	result->num_data[num_sk++] = k;

	make_empty_que_c2_n(que);

	// enq pattern of 0...10
	qu.cursor = new_que_e2_n(que);
	qu.key = bd[bd_idx[1]];
	que->details[qu.cursor].sk = qs->sketch ^ (1 << bd_idx[1]);
	que->details[qu.cursor].pt = 1 << 1; // pt = "0...00000010"
	enq_c2_n(&qu, que);		

	while(deq_c2_n(&qu, que) && k < num_candidates) {
		s = que->details[qu.cursor].sk; // 列挙のつぎのスケッチ
		if(bkt[s] < bkt[s + 1]) { // 空でないバケツの範囲に正解があるかをチェックして，見つかれば，列挙したスケッチ数とデータ数をresultに記録して復帰
			num_ne++;
			k += bkt[s + 1] - bkt[s];
			if(qs->answer.data_num >= bkt[s] && qs->answer.data_num < bkt[s + 1]) {
				result->num_enum_sketches = num_sk;
				result->num_enum_nonempty_sketches = num_ne;
				result->num_enum_data = k;
				if(num_sk < 10) result->num_data[num_sk++] = k;
				return;
			}
		}
		if(num_sk < 10) {
			result->num_data[num_sk++] = k;
		} else {
			num_sk++;
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

	return;
}

void check_sketch_enumeration_hamming(struct_query_sketch *qs, struct_bucket *bucket, struct_check_result *result, int num_candidates, int low, int add, int add_2)
{
	static sub_dimension *table[PJT_DIM + 1] = {NULL};
	static int table_size[PJT_DIM + 1] = {0};

	if(low == 0) {
//		fprintf(stderr, "\n");
//		for(int i = 0; i < PJT_DIM; i++) {
//			fprintf(stderr, "%3d", qs->bd[qs->idx[i]]);
//		}
//		fprintf(stderr, "\n");

//		fprintf(stderr, "low-add-add_2 = %d-%d-%d ", low, add, add_2);
		int b = qs->bd[qs->idx[0]];
		for(int i = 0; i < PJT_DIM; i++) {
			if(qs->bd[qs->idx[i]] > b) {
				low = i;
				break;
			}
		}
		b = qs->bd[qs->idx[low]];
		for(int i = 0; i + low < PJT_DIM; i++) {
			if(qs->bd[qs->idx[i + low]] > b) {
				add = i;
				break;
			}
		}
		b = qs->bd[qs->idx[low + add]];
		for(int i = 0; i + low + add < PJT_DIM; i++) {
			if(qs->bd[qs->idx[i + low + add]] > b) {
				add_2 = i;
				break;
			}
		}
//		fprintf(stderr, " -> %d-%d-%d\n", low, add, add_2);
//		getchar();
	}
//	fprintf(stderr, "%d-%d-%d -> ", low, add, add_2);
//	while(qs->bd[qs->idx[low - 1]] == qs->bd[qs->idx[low]]) {
//		low++; if(add > 0) add--;
//	}
//	fprintf(stderr, "%d-%d-%d -> ", low, add, add_2);
//	while(low + add < PJT_DIM && qs->bd[qs->idx[low + add - 1]] == qs->bd[qs->idx[low + add]]) {
//		add++; if(add_2 > 0) add_2--;
//	}
//	fprintf(stderr, "%d-%d-%d\n", low, add, add_2);
	if(low > PJT_DIM) low = PJT_DIM;
	if(low + add > PJT_DIM) add = PJT_DIM - low;
	if(low + add + add_2 > PJT_DIM) add_2 = PJT_DIM - low - add;
	int res = PJT_DIM - low - add - add_2;
	if(table[low] == NULL) {
		table_size[low] = (1 << low) + 1;
		fprintf(stderr, "make table for %d-bit, table_size = %d: ", low, table_size[low]);
		table[low] = make_table_for_enumeration_hamming(table_size[low], low);
		fprintf(stderr, "OK\n");
	}
	if(add > 0 && table[add] == NULL) {
		table_size[add] = (1 << add) + 1;
		fprintf(stderr, "make table for %d-bit, table_size = %d: ", add, table_size[add]);
		table[add] = make_table_for_enumeration_hamming(table_size[add], add);
		fprintf(stderr, "OK\n");
	}
	if(add_2 > 0 && table[add_2] == NULL) {
		table_size[add_2] = (1 << add_2) + 1;
		fprintf(stderr, "make table for %d-bit, table_size = %d: ", add_2, table_size[add_2]);
		table[add_2] = make_table_for_enumeration_hamming(table_size[add_2], add_2);
		fprintf(stderr, "OK\n");
	}
	if(res > 0 && table[res] == NULL) {
		table_size[res] = (1 << res) + 1;
		fprintf(stderr, "make table for %d-bit, table_size = %d: ", res, table_size[res]);
		table[res] = make_table_for_enumeration_hamming(table_size[res], res);
		fprintf(stderr, "OK\n");
	}

	int *bd_idx = qs->idx;
	int *bkt = bucket->bkt;

	result->num_enum_sketches = result->num_enum_nonempty_sketches = result->num_enum_data = -1;	// 列挙する候補数が num_candidates 以下では正解が得られないとき．

	sketch_type s, mu_low = 0, mu_add = 0, mu_add_2 = 0, mu_res = 0;
	int k = 0; 				// 列挙されたスケッチを持つデータ総数
	int num_sk = 0; 			// 列挙されたスケッチ総数（空も含む）
	int num_ne = 0; 		// 列挙によって求めた空でないスケッチの個数
	int add_c, add_p = -1; // add-bit のカウンタ，直前のカウンタ
	int add_2_c, add_2_p = -1; // add_2-bit のカウンタ，直前のカウンタ
	int res_c, res_p = -1; // res-bit のカウンタ，直前のカウンタ

	while(k < num_candidates * FACTOR_INF) {
		mu_low = Mu(0, low, bd_idx, num_sk % (1 << low), table[low]);
		if(add > 0) {
			add_c = (num_sk % (1 << (low + add))) / (1 << low);
			if(add_c != add_p) {
				mu_add = Mu(low, add, bd_idx, add_c, table[add]);
				add_p = add_c;
			}
			if(add_2 > 0) {
				add_2_c = (num_sk % (1 << (low + add + add_2))) / (1 << (low + add));
				if(add_2_c != add_2_p) {
					mu_add_2 = Mu(low + add, add_2, bd_idx, add_2_c, table[add_2]);
					add_2_p = add_2_c;
				}
				if(res > 0) {
					res_c = num_sk / (1 << (low + add + add_2));
					if(res_c != res_p) {
						mu_res = Mu(low + add + add_2, res, bd_idx, res_c, table[res]);
						res_p = res_c;
					}
				}
			}
		}
		s = qs->sketch ^ mu_low ^ mu_add ^ mu_add_2 ^ mu_res;
		if(bkt[s] < bkt[s + 1]) { // 空でないバケツの範囲に正解があるかをチェックして，見つかれば，列挙したスケッチ数とデータ数をresultに記録して復帰
			num_ne++;
			k += bkt[s + 1] - bkt[s];
			if(qs->answer.data_num >= bkt[s] && qs->answer.data_num < bkt[s + 1]) {
				result->num_enum_sketches = num_sk;
				result->num_enum_nonempty_sketches = num_ne;
				result->num_enum_data = k;
				return;
			}
		}
		num_sk++;
	}

	return;
}

#endif

/*
バケットの表現

0.	前提
	int data_num; // データ数
	ftr[0], ... , ftr[data_num - 1]; // 特徴データ（オリジナル）
	sk[0], ... , sk[data_num - 1]; // スケッチ
	sftr[0], ... , sftr[data_num - 1]; // 特徴データをスケッチ順に並べたもの
	idx[i] = "sftr[i] のオリジナルのインデックス" （sftr[i] = ftr[idx[i]]）
			「スケッチ順で i 番目のデータがオリジナルで何番目であるかと表す」
	bkt[s] = スケッチ順にデータを並べたときに，スケッチが s であるものの先頭の位置 (0 - data_num - 1)
	bkt[s + 1] - bkt[s] = スケッチが s であるもののデータ数

	スケッチが s であるデータ = sftr[j] (j = bkt[s], ... , bkt[s + 1] - 1) (sftr が使えるとき)
	スケッチが s であるデータ = ftr[idx[j]] (j = bkt[s], ... , bkt[s + 1] - 1) 
	スケッチが s であるデータのオリジナルインデックス = idx[j] (j = bkt[s], ... , bkt[s + 1] - 1)
   	
1. 	プログラム内部

	オプション1 (FTR_ON): 特徴データを主記憶に置くかどうか (FTR_ON = MAIN_MEMORY | SECONDARY_MEMORY)
	オプション2 (FTR_ARRANGEMANT): 特徴データをスケッチ順にソートしたものを用いるかどうか (FTR_ARRANGEMENT = SORT_BY_SKETCH | ASIS)
	オプション3 (FILTERING_BY): フィルタリング手法 (FILTERING_BY = SKETCH_ENUMERATION | SEQUENTIAL_SEARCH | SEQUENTIAL_SEARCH_USING_BUCKET)

1.1	オンメモリ（すべて主記憶上に読み込んでいる）
	FTR_ON = MAIN_MEMORY, FTR_ARRANGEMANT = SORT_BY_SKETCH | ASIS, FILTERING_BY = SKETCH_ENUMERATION
	database[data_num]; 			// 特徴データ（オリジナル）
	sketch[data_num];  				// データのスケッチ
	bkt[2^w + 2], idx[data_num];	// バケット表とデータのインデックス
	database2[data_num]; 			// スケッチ順の特徴データ

1.2	特徴データを2次記憶（HDD，SSD）に置いたまま
	FTR_ON = SECONDARY_MEMORY, FTR_ARRANGEMANT = SORT_BY_SKETCH | ASIS, FILTERING_BY = SKETCH_ENUMERATION
	bkt[2^w + 2], idx[data_num];	// バケット表とデータのインデックス
	特徴データは，sftrを2次記憶上に置いておく．
	検索は，スケッチ列挙法による．

1.3	Sequential Filtering （オンメモリ）
	FTR_ON = MAIN_MEMORY, FTR_ARRANGEMANT = ASIS, FILTERING_BY = SEQUENTIAL_SEARCH
	database[data_num]; 			// 特徴データ（オリジナル）
	sketch[data_num];  				// データのスケッチ

1.4	Sequential Filtering using bucket
	FTR_ON = SECONDARY_MEMORY, FTR_ARRANGEMANT = SORT_BY_SKETCH | ASIS, FILTERING_BY = SEQUENTIAL_SEARCH_USING_BUCKET
	database[data_num]; 			// 特徴データ（オリジナル）
	bkt[2^w + 2], idx[data_num];	// バケット表とデータのインデックス

		
2. 外部ファイルの保存時

*/

#ifndef TRIAL
#define TRIAL 5
#endif

// スケッチの配列を相対的にソートする．idxを入れ替える．
int find_pivot_for_sketch(int idx[], sketch_type sk[], int i, int j)
{
   int k;

	for(int t = 0; t < TRIAL; t++) {
		k = random() % (j - i + 1) + i;
		if(sk[idx[i]] != sk[idx[k]]) {
			return (sk[idx[i]] > sk[idx[k]] ? i : k);
		}
	}
	for(k = i + 1; k <= j; k++) {
		if(sk[idx[i]] != sk[idx[k]]) {
			return (sk[idx[i]] > sk[idx[k]] ? i : k);
		}
	}
	return -1;
}

int partition_by_pivot_for_sketch(int idx[], sketch_type sk[], int i, int j, sketch_type piv)
{
   int left, right;
   int temp;
   left = i;   right = j;
   do {
      while(sk[idx[left]] < piv) left++;
      while(sk[idx[right]] >= piv) right--;
      if (left < right) { 
         temp = idx[left];
         idx[left] = idx[right];
         idx[right] = temp;
      }
   } while(left <= right);
   return left;
}

void quick_sort_for_sketch(int idx[], sketch_type sk[], int i, int j)
{
	int pivotindex, k;
	pivotindex = find_pivot_for_sketch(idx, sk, i, j);
	if (pivotindex >= 0) {
		k = partition_by_pivot_for_sketch(idx, sk, i, j, sk[idx[pivotindex]]);
		quick_sort_for_sketch(idx, sk, i, k - 1);
		quick_sort_for_sketch(idx, sk, k, j);
	}
}

// 同じ点との距離を繰り返し計算するときに点データをコピーしておくための配列．
static int point_a[FTR_DIM];
#if defined(PARTITION_TYPE_SQBP)
static int used[FTR_DIM], num_used;
#elif defined(PARTITION_TYPE_CSQBP)
static int j_start, num_j;
#endif

void set_dist_L2_22(ftr_type a)
{
	int j;
	for(j = 0; j < FTR_DIM; j++) 
		point_a[j] = a[j];
}

void set_part_dist_L2_22(ftr_type a, int dim_start, int dim)
{
	int j;
	for(j = dim_start; j < dim_start + dim; j++) {
		point_a[j] = a[j];
	}
}

#ifdef PARTITION_TYPE_SQBP
void set_sub_dist_L2_22(pivot_type *pivot, int dim)
{
	ftr_type a = pivot->p[dim];
	num_used = pivot->num_used[dim];
	for(int j = 0; j < FTR_DIM; j++) {
		used[j] = pivot->used[dim][j];
		point_a[j] = a[j];
	}
}
#endif

#ifndef SQRT_FTR

dist_type dist_L1(ftr_type a, ftr_type b, int dim)
{
	dist_type s = 0;
	for(int j = 0; j < dim; j++)  {
		s += abs((int)a[j] - (int)b[j]);
	}
	return s;
}

dist_type dist_L1_2(ftr_type a, ftr_type b, int dim, dist_type dist)
{
	dist_type s = 0;
	for(int j = 0; j < dim; j++) {
		s += abs((int)a[j] - (int)b[j]);
		if(s >= dist) return s;
	}
	return s;
}

dist_type dist_L1_22(ftr_type b)
{
	dist_type s = 0;
	for(int j = 0; j < FTR_DIM; j++)  {
		s += abs(point_a[j] - (int)b[j]);
	}
	return s;
}

dist_type dist_L2(ftr_type a, ftr_type b, int dim)
// 注意：平方根を取っていない
{
	dist_type s = 0;
	for(int j = 0; j < dim; j++)  {
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
	}
	return s;
}

dist_type dist_L2_2(ftr_type a, ftr_type b, int dim, dist_type dist)
{
	dist_type s = 0;
	for(int j = 0; j < dim; j++)  {
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
		if(s >= dist) return s;
	}
	return s;
}

dist_type dist_L2_22(ftr_type b)
{
	dist_type s = 0;
	for(int j = 0; j < FTR_DIM; j++)  {
		s += (point_a[j] - (int)b[j]) * (point_a[j] - (int)b[j]);
	}
	return s;
}

#ifdef PARTITION_TYPE_PQBP
dist_type part_dist_L2(ftr_type a, ftr_type b, int dim_start, int dim)
// 注意：平方根を取っていない
{
	dist_type s = 0;
	for(int j = dim_start; j < dim_start + dim; j++)  {
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
	}
	return s;
}

dist_type part_dist_L2_2(ftr_type a, ftr_type b, int dim_start, int dim, dist_type dist)
{
	dist_type s = 0;
	for(int j = dim_start; j < dim_start + dim; j++)  {
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
		if(s >= dist) return s;
	}
	return s;
}

dist_type part_dist_L1_22(ftr_type b, int dim_start, int dim)
{
	dist_type s = 0;
	for(int j = dim_start; j < dim_start + dim; j++) {
		s += abs(point_a[j] - (int)b[j]);
	}
	return s;
}

dist_type part_dist_L2_22(ftr_type b, int dim_start, int dim)
{
	dist_type s = 0;
	for(int j = dim_start; j < dim_start + dim; j++) {
		s += (point_a[j] - (int)b[j]) * (point_a[j] - (int)b[j]);
	}
	return s;
}

#endif

#else // SQRT_FTR

dist_type dist_L1(ftr_type a, ftr_type b, int dim)
{
	dist_type s = 0;
	for(int j = 0; j < dim; j++)  {
		s += abs((int)a[j] * a[j] - (int)b[j] * b[j]);
	}
	return s;
}

dist_type dist_L1_2(ftr_type a, ftr_type b, int dim, dist_type dist)
{
	dist_type s = 0;
	for(int j = 0; j < dim; j++) {
		s += abs((int)a[j] * a[j] - (int)b[j] * b[j]);
		if(s >= dist) return s;
	}
	return s;
}

dist_type dist_L1_22(ftr_type b)
{
	dist_type s = 0;
	for(int j = 0; j < FTR_DIM; j++)  {
		s += abs((int)point_a[j] * point_a[j] - (int)b[j] * b[j]);
	}
	return s;
}

dist_type dist_L2(ftr_type a, ftr_type b, int dim)
// 注意：SQRT_FTRではintをオーバーフローするので平方根を取っている
{
	double s = 0;
	for(int j = 0; j < dim; j++) {
		s += ((double)a[j] - (double)b[j]) * ((double)a[j] + (double)b[j]) * ((double)a[j] - (double)b[j]) * ((double)a[j] + (double)b[j]);
	}
	return sqrt(s);
}

dist_type dist_L2_2(ftr_type a, ftr_type b, int dim, dist_type dist)
{
	double s = 0;
	for(int j = 0; j < dim; j++) {
		s += ((double)a[j] - (double)b[j]) * ((double)a[j] + (double)b[j]) * ((double)a[j] - (double)b[j]) * ((double)a[j] + (double)b[j]);
		if(s >= (double)dist * dist) return sqrt(s);
	}
	return sqrt(s);
}

dist_type dist_L2_22(ftr_type b)
{
	double s = 0;
	for(int j = 0; j < FTR_DIM; j++)  {
		s += ((double)point_a[j] * point_a[j] - (double)b[j] * b[j]) * ((double)point_a[j] * point_a[j] - (double)b[j] * b[j]);
	}
	return sqrt(s);
}

#ifdef PARTITION_TYPE_PQBP
dist_type part_dist_L2(ftr_type a, ftr_type b, int dim_start, int dim)
// 注意：SQRT_FTRではintをオーバーフローするので平方根を取っている
{
	double s = 0;
	for(int j = dim_start; j < dim_start + dim; j++)  {
		s += ((double)a[j] * a[j] - (double)b[j] * b[j]) * ((double)a[j] * a[j] - (double)b[j] * b[j]);
	}
	return sqrt(s);
}

dist_type part_dist_L2_2(ftr_type a, ftr_type b, int dim_start, int dim, dist_type dist)
{
	double s = 0;
	for(int j = dim_start; j < dim_start + dim; j++)  {
		s += ((double)a[j] * a[j] - (double)b[j] * b[j]) * ((double)a[j] * a[j] - (double)b[j] * b[j]);
		if(s >= (double)dist * dist) return sqrt(s);
	}
	return sqrt(s);
}

dist_type part_dist_L1_22(ftr_type b, int dim_start, int dim)
{
	dist_type s = 0;
	for(int j = dim_start; j < dim_start + dim; j++) {
		s += abs((int)point_a[j] * point_a[j] - (int)b[j] * b[j]);
	}
	return s;
}

dist_type part_dist_L2_22(ftr_type b, int dim_start, int dim)
{
	double s = 0;
	for(int j = dim_start; j < dim_start + dim; j++) {
		s += ((double)point_a[j] * point_a[j] - (double)b[j] * b[j]) * ((double)point_a[j] * point_a[j] - (double)b[j] * b[j]);
	}
	return sqrt(s);
}

#endif

#endif

#ifdef PARTITION_TYPE_SQBP

dist_type sub_dist_L2(ftr_type a, pivot_type *pivot, int dim)
// 注意：平方根を取っていない
{
	ftr_type b = pivot->p[dim];
	int *used = pivot->used[dim], num_used = pivot->num_used[dim];
	dist_type s = 0;
	
	for(int i = 0; i < num_used; i++)  {
		int j = used[i];
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
	}
	return s;
}

dist_type sub_dist_L2_2(ftr_type a, pivot_type *pivot, int dim, dist_type dist)
{
	ftr_type b = pivot->p[dim];
	int *used = pivot->used[dim], num_used = pivot->num_used[dim];
	dist_type s = 0;
	
	for(int i = 0; i < num_used; i++)  {
		int j = used[i];
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
		if(s >= dist) return s;
	}
	return s;
}

dist_type sub_dist_L2_22(ftr_type b)
{
	int j;
	dist_type s = 0;
	for(int i = 0; i < num_used; i++) {
		j = used[i];
		s += (point_a[j] - (int)b[j]) * (point_a[j] - (int)b[j]);
	}
	return s;
}

#endif

#ifdef PARTITION_TYPE_CSQBP

dist_type sub_dist_L2(ftr_type a, pivot_type *pivot, int dim)
// 注意：平方根を取っていない
{
	ftr_type b = pivot->p[dim];
	int j_start = pivot->j_start[dim], num_j = pivot->num_j[dim];
	dist_type s = 0;
	
	for(int i = j_start; i < j_start + num_j; i++)  {
		int j = i % PJT_DIM;
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
	}
	return s;
}

dist_type sub_dist_L2_2(ftr_type a, pivot_type *pivot, int dim, dist_type dist)
{
	ftr_type b = pivot->p[dim];
	int j_start = pivot->j_start[dim], num_j = pivot->num_j[dim];
	dist_type s = 0;
	
	for(int i = j_start; i < j_start + num_j; i++)  {
		int j = i % PJT_DIM;
		s += ((int)a[j] - (int)b[j]) * ((int)a[j] - (int)b[j]);
		if(s >= dist) return s;
	}
	return s;
}

dist_type sub_dist_L2_22(ftr_type b)
{
	dist_type s = 0;
	for(int i = j_start; i < j_start + num_j; i++)  {
		int j = i % PJT_DIM;
		s += (point_a[j] - (int)b[j]) * (point_a[j] - (int)b[j]);
	}
	return s;
}

#endif

#ifndef USE_PD_SKETCH
#ifndef SQRT_FTR
dist_type dist_pivot_L2(ftr_type a, ftr_type piv_center, int dim)
// 注意：平方根を取っていない
{
	#ifdef DEBUG_SQRT_FTR
	static int first = 1;
	if(first) {
		fprintf(stderr, "(1) dist_pivot_L2 (WITHOUT USE_PD_SKETCH)\n");
		first = 0;
	}
	#endif

	dist_type s = 0;
	for(int j = 0; j < dim; j++)  {
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
	return s;
}

dist_type dist_pivot_L2_2(ftr_type a, ftr_type piv_center, int dim, dist_type dist)
{
	#ifdef DEBUG_SQRT_FTR
	static int first = 1;
	if(first) {
		fprintf(stderr, "(2) dist_pivot_L2_2 (WITHOUT USE_PD_SKETCH)\n");
		first = 0;
	}
	#endif

	dist_type s = 0;
	for(int j = 0; j < dim; j++)  {
		#ifdef IGNORE_MED
			if(piv_center[j] == FTR_MIN) {
				s += ((int)a[j] - FTR_MIN) * ((int)a[j] - FTR_MIN);
			} else if(piv_center[j] == FTR_MAX) {
				s += ((int)a[j] - FTR_MAX) * ((int)a[j] - FTR_MAX);
			}
		#else
			s += ((int)a[j] - (int)piv_center[j]) * ((int)a[j] - (int)piv_center[j]);
		#endif
		if(s >= dist) return s;
	}
	return s;
}

// set_dist_L2 で pivot_center が point_a に set されているときに用いること
dist_type dist_pivot_L2_22(ftr_type b)
{
	#ifdef DEBUG_SQRT_FTR
	static int first = 1;
	if(first) {
		fprintf(stderr, "(3) dist_pivot_L2_22 (WITHOUT USE_PD_SKETCH)\n");
		first = 0;
	}
	#endif

	dist_type s = 0;
	for(int j = 0; j < FTR_DIM; j++)  {
		#ifdef IGNORE_MED
			if(point_a[j] == FTR_MIN) {
				s += ((int)b[j] - FTR_MIN) * ((int)b[j] - FTR_MIN);
			} else if(point_a[j] == FTR_MAX) {
				s += ((int)b[j] - FTR_MAX) * ((int)b[j] - FTR_MAX);
			}
		#else
			s += ((int)point_a[j] - (int)b[j]) * ((int)point_a[j] - (int)b[j]);
		#endif
	}
	return s;
}
#else
dist_type dist_pivot_L2(ftr_type a, ftr_type piv_center, int dim)
// 注意：SQRT_FTRのときはオーバーフローする可能性があるので平方根を取る
{
	#ifdef DEBUG_SQRT_FTR
	static int first = 1;
	if(first) {
		fprintf(stderr, "(4) dist_pivot_L2 SQRT_FTR (WITHOUT USE_PD_SKETCH)\n");
		first = 0;
	}
	#endif

	double s = 0;
	for(int j = 0; j < dim; j++)  {
		#ifdef IGNORE_MED
			if(piv_center[j] == FTR_MIN) {
				s += ((double)a[j] * a[j] - (double)FTR_MIN * FTR_MIN) * ((double)a[j] * a[j] - (double)FTR_MIN * FTR_MIN);
			} else if(piv_center[j] == FTR_MAX) {
				s += ((double)a[j] * a[j] - (double)FTR_MAX * FTR_MAX) * ((double)a[j] * a[j] - (double)FTR_MAX * FTR_MAX);
			}
		#else
			s += ((double)a[j] * a[j] - (double)piv_center[j] * piv_center[j]) * ((double)a[j] * a[j] - (double)piv_center[j] * piv_center[j]);
		#endif
	}
	return sqrt(s);
}

dist_type dist_pivot_L2_2(ftr_type a, ftr_type piv_center, int dim, dist_type dist)
{
	#ifdef DEBUG_SQRT_FTR
	static int first = 1;
	if(first) {
		fprintf(stderr, "(5) dist_pivot_L2_2 SQRT_FTR (WITHOUT USE_PD_SKETCH)\n");
		first = 0;
	}
	#endif

	double s = 0;
	for(int j = 0; j < dim; j++)  {
		#ifdef IGNORE_MED
			if(piv_center[j] == FTR_MIN) {
				s += ((double)a[j] * a[j] - (double)FTR_MIN * FTR_MIN) * ((double)a[j] * a[j] - (double)FTR_MIN * FTR_MIN);
			} else if(piv_center[j] == FTR_MAX) {
				s += ((double)a[j] * a[j] - (double)FTR_MAX * FTR_MAX) * ((double)a[j] * a[j] - (double)FTR_MAX * FTR_MAX);
			}
		#else
			s += ((double)a[j] * a[j] - (double)piv_center[j] * piv_center[j]) * ((double)a[j] * a[j] - (double)piv_center[j] * piv_center[j]);
		#endif
		if(s >= (double)dist * dist) return sqrt(s);
	}
	return sqrt(s);
}

// set_dist_L2 で pivot_center が point_a に set されているときに用いること
dist_type dist_pivot_L2_22(ftr_type b)
{
	#ifdef DEBUG_SQRT_FTR
	static int first = 1;
	if(first) {
		fprintf(stderr, "(6) dist_pivot_L2_22 SQRT_FTR (WITHOUT USE_PD_SKETCH)\n");
		first = 0;
	}
	#endif

	double s = 0;
	for(int j = 0; j < FTR_DIM; j++)  {
		#ifdef IGNORE_MED
			if(point_a[j] == FTR_MIN) {
				s += ((double)b[j] * b[j] - (double)FTR_MIN * FTR_MIN) * ((double)b[j] * b[j] - (double)FTR_MIN * FTR_MIN);
			} else if(point_a[j] == FTR_MAX) {
				s += ((double)b[j] * b[j] - (double)FTR_MAX * FTR_MAX) * ((double)b[j] * b[j] - (double)FTR_MAX * FTR_MAX);
			}
		#else
			s += ((double)point_a[j] * point_a[j] - (double)b[j] * b[j]) * ((double)point_a[j] * point_a[j] - (double)b[j] * b[j]);
		#endif
	}
	return sqrt(s);
}
#endif
#else
// #ifdef USE_PD_SKETCH
dist_type dist_pivot_L2(ftr_type a, ftr_type piv_center, int num_axis, int axis[])
// 注意：平方根を取っていない
{
	#ifdef DEBUG_SQRT_FTR
	static int first = 1;
	if(first) {
		fprintf(stderr, "(7) dist_pivot_L2 (USE_PD_SKETCH)\n");
		first = 0;
	}
	#endif

	dist_type s = 0;
	for(int i = 0; i < num_axis; i++)  {
		int j = axis[i];
		s += ((int)a[j] - (int)piv_center[j]) * ((int)a[j] - (int)piv_center[j]);
	}
	return s;
}

dist_type dist_pivot_L2_2(ftr_type a, ftr_type piv_center, int num_axis, int axis[], dist_type dist)
{
	#ifdef DEBUG_SQRT_FTR
	static int first = 1;
	if(first) {
		fprintf(stderr, "(8) dist_pivot_L2_2 (USE_PD_SKETCH)\n");
		first = 0;
	}
	#endif

	dist_type s = 0;
	for(int i = 0; i < num_axis; i++)  {
		int j = axis[i];
		s += ((int)a[j] - (int)piv_center[j]) * ((int)a[j] - (int)piv_center[j]);
		if(s >= dist) return s;
	}
	return s;
}

// set_dist_L2 で pivot_center が point_a に set されているときに用いること
dist_type dist_pivot_L2_22(ftr_type b, int num_axis, int axis[])
{
	#ifdef DEBUG_SQRT_FTR
	static int first = 1;
	if(first) {
		fprintf(stderr, "(9) dist_pivot_L2_22 (USE_PD_SKETCH)\n");
		first = 0;
	}
	#endif

	dist_type s = 0;
	for(int i = 0; i < num_axis; i++)  {
		int j = axis[i];
		s += ((int)point_a[j] - (int)b[j]) * ((int)point_a[j] - (int)b[j]);
	}
	return s;
}
#endif

#ifdef PARTITION_TYPE_PQBP
dist_type part_dist_pivot_L2(ftr_type a, ftr_type piv_center, int dim_start, int dim)
// 注意：平方根を取っていない
{
	#ifdef SQRT_FTR
	fprintf(stderr, "part_dist_pivot_L2 is invoked\n"); exit(0);
	#endif
	int j;
	dist_type s = 0;
	
	for(j = dim_start; j < dim_start + dim; j++)  {
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
	return s;
}

dist_type part_dist_pivot_L2_2(ftr_type a, ftr_type piv_center, int dim_start, int dim, dist_type dist)
{
	#ifdef SQRT_FTR
	fprintf(stderr, "part_dist_pivot_L2_2 is invoked\n"); exit(0);
	#endif
	int j;
	dist_type s = 0;
	
	for(j = dim_start; j < dim_start + dim; j++)  {
		#ifdef IGNORE_MED
		if(piv_center[j] == FTR_MIN) {
			s += ((int)a[j] - FTR_MIN) * ((int)a[j] - FTR_MIN);
		} else if(piv_center[j] == FTR_MAX) {
			s += ((int)a[j] - FTR_MAX) * ((int)a[j] - FTR_MAX);
		}
		#else
		s += ((int)a[j] - (int)piv_center[j]) * ((int)a[j] - (int)piv_center[j]);
		#endif
		if(s >= dist) return s;
	}
	return s;
}

dist_type part_dist_pivot_L2_22(ftr_type b, int dim_start, int dim)
{
	#ifdef SQRT_FTR
	fprintf(stderr, "part_dist_pivot_L2_22 is invoked\n"); exit(0);
	#endif
	int j;
	dist_type s = 0;
	for(j = dim_start; j < dim_start + dim; j++) {
		#ifdef IGNORE_MED
		if(point_a[j] == FTR_MIN) {
			s += ((int)b[j] - FTR_MIN) * ((int)b[j] - FTR_MIN);
		} else if(point_a[j] == FTR_MAX) {
			s += ((int)b[j] - FTR_MAX) * ((int)b[j] - FTR_MAX);
		}
		#else
		s += ((int)point_a[j] - (int)b[j]) * ((int)point_a[j] - (int)b[j]);
		#endif
	}
	return s;
}

#endif

// #endif

// answer_typeの配列をソートする．
// quick sort のための 関数群 find_pivot, partition_by _pivot, quick_sort

// PUBMED23 のときは，answer の比較に，data_numも用いる．
int find_pivot_answer(answer_type ans[], int i, int j)

{
// 初めに 5 回だけピボット候補をランダムに選ぶ．
// i 番目のデータ ans[i].dist と k（ただし，k = i + 1, ... , j）番目のデータ ans[k].dist を比較し，
// すべて等しいときには -1 を返し, 
// そうでないときには, ans[i].dist と異なる ans[k].dist  で最初に 
// 現れたもののうちで, 大きい方の位置(i または k) を返す. 
	int k;

	for(int t = 0; t < 5; t++) {
		k = random() % (j - i + 1) + i;
		if(ans[i].dist != ans[k].dist) {
			return (ans[i].dist > ans[k].dist ? i : k);
		}
	}
	for(k = i + 1; k <= j; k++) {
		if(ans[i].dist != ans[k].dist) {
			return (ans[i].dist > ans[k].dist ? i : k);
		}
		#ifdef PUBMED23
		else {
//			fprintf(stderr, "find_pivot_answer at tie\n");
//			fprintf(stderr, "ans[%d].data_num = %d, ans[%d].data_num = %d\n", i, ans[i].data_num, k, ans[k].data_num); 
			return (ans[i].data_num > ans[k].data_num ? i : k);
		}
		#endif
	}
	return -1;
}


#ifdef PUBMED23
int partition_by_pivot_answer(answer_type ans[], int i, int j, answer_type piv)
// ans[i], ... , ans[j] をそれらの dist と piv との大小によって分け， 
// piv より小さいものが ans[i], ... , ans[k-1] に，   
// そうでないものが ans[k], ... , ans[j] に来るようにする. 
// 右側のリストの先頭の位置(k)を返す. 
{
   int left, right;
   answer_type temp;
   left = i;   right = j;
   do {
//    while(ans[left].dist < piv) left++;
//fprintf(stderr, "comp_answer = ");
//int chk = comp_answer(&ans[left], &piv);
//fprintf(stderr, "%d\n", chk);
      while(comp_answer(&ans[left], &piv) < 0) left++;
//      while(ans[right].dist >= piv) right--;
      while(comp_answer(&ans[right], &piv) >= 0) right--;
	  if (left < right) { 
         temp = ans[left];
         ans[left] = ans[right];
         ans[right] = temp;
      }
   } while(left <= right);
   return left;
}
#else

int partition_by_pivot_answer(answer_type ans[], int i, int j, dist_type piv)
// ans[i], ... , ans[j] をそれらの dist と piv との大小によって分け， 
// piv より小さいものが ans[i], ... , ans[k-1] に，   
// そうでないものが ans[k], ... , ans[j] に来るようにする. 
// 右側のリストの先頭の位置(k)を返す. 
{
   int left, right;
   answer_type temp;
   left = i;   right = j;
   do {
      while(ans[left].dist < piv) left++;
      while(ans[right].dist >= piv) right--;
      if (left < right) { 
         temp = ans[left];
         ans[left] = ans[right];
         ans[right] = temp;
      }
   } while(left <= right);
   return left;
}
#endif


#ifdef PUBMED23
void quick_sort_answer(answer_type ans[], int i, int j)
{
	int pivotindex, k;
	pivotindex = find_pivot_answer(ans, i, j);
	if (pivotindex >= 0) {
		k = partition_by_pivot_answer(ans, i, j, ans[pivotindex]);
		quick_sort_answer(ans, i, k - 1);
		quick_sort_answer(ans, k, j);
	}
}
#else

void quick_sort_answer(answer_type ans[], int i, int j)
{
	int pivotindex, k;
	pivotindex = find_pivot_answer(ans, i, j);
	if (pivotindex >= 0) {
		k = partition_by_pivot_answer(ans, i, j, ans[pivotindex].dist);
		quick_sort_answer(ans, i, k - 1);
		quick_sort_answer(ans, k, j);
	}
}
#endif


#ifdef PUBMED23
void quick_select_k_r_answer(answer_type ans[], int i, int j, int k)
{
	int pivotindex, m;
	pivotindex = find_pivot_answer(ans, i, j);
	if (pivotindex >= 0) {
		m = partition_by_pivot_answer(ans, i, j, ans[pivotindex]);
		if(k < m) {
			quick_select_k_r_answer(ans, i, m - 1, k);
		} else if(k > m) {
			quick_select_k_r_answer(ans, m, j, k);
		}
	}
}
#else

void quick_select_k_r_answer(answer_type ans[], int i, int j, int k)
{
	int pivotindex, m;
	pivotindex = find_pivot_answer(ans, i, j);
	if (pivotindex >= 0) {
		m = partition_by_pivot_answer(ans, i, j, ans[pivotindex].dist);
		if(k < m) {
			quick_select_k_r_answer(ans, i, m - 1, k);
		} else if(k > m) {
			quick_select_k_r_answer(ans, m, j, k);
		}
	}
}
#endif


#ifdef PUBMED23
void quick_select_k_answer(answer_type ans[], int i, int j, int k)
{
	int pivotindex, m;
	while((pivotindex = find_pivot_answer(ans, i, j)) >= 0) {
		m = partition_by_pivot_answer(ans, i, j, ans[pivotindex]);
		if(k < m) {
			j = m - 1;
		} else if(k > m) {
			i = m;
		} else {
			break;
		}
	}
}
#else

void quick_select_k_answer(answer_type ans[], int i, int j, int k)
{
	int pivotindex, m;
	while((pivotindex = find_pivot_answer(ans, i, j)) >= 0) {
		m = partition_by_pivot_answer(ans, i, j, ans[pivotindex].dist);
		if(k < m) {
			j = m - 1;
		} else if(k > m) {
			i = m;
		} else {
			break;
		}
	}
}
#endif

// sketch_with_priority_num の配列を優先度でソートしたり，優先度が上位のものでデータ数の合計が k 以上になる最小個数のスケッチを選択する．
int find_pivot_sketch_with_priority_num(sketch_with_priority_num buff[], int i, int j)
{
// 初めに 5 回だけピボット候補をランダムに選ぶ．
// i 番目のデータ buff[i].priority と k（ただし，k = i + 1, ... , j）番目のデータ buff[k].priority を比較し，
// すべて等しいときには -1 を返し, 
// そうでないときには, buff[i].priority と異なる buff[k].priority で最初に 
// 現れたもののうちで, 大きい方の位置(i または k) を返す. 
	int k, t;

	for(t = 0; t < 5; t++) {
		k = random() % (j - i + 1) + i;
		if(buff[i].priority != buff[k].priority) {
			return (buff[i].priority > buff[k].priority ? i : k);
		}
	}
	for(k = i + 1; k <= j; k++) {
		if(buff[i].priority != buff[k].priority) {
			return (buff[i].priority > buff[k].priority ? i : k);
		}
	}
	return -1;
}

int partition_by_pivot_sketch_with_priority_num(sketch_with_priority_num buff[], int i, int j, dist_type piv, int *sum)
// buff[i], ... , buff[j] をそれらの priority と piv との大小によって分け， 
// piv より小さいものが buff[i], ... , buff[k-1] に，   
// そうでないものが buff[k], ... , buff[j] に来るようにする. 
// 右側のリストの先頭の位置(k)を返す. 
// *sum に左側のリストの num の合計を返す．
{
	int s = 0;
	int left, right;
	sketch_with_priority_num temp;
	left = i;   right = j;
	do {
		while(buff[left].priority < piv) {
			s += buff[left].num;
			left++;
		}
		while(buff[right].priority >= piv) right--;
		if (left < right) { 
			temp = buff[left];
			buff[left] = buff[right];
			buff[right] = temp;
		}
	} while(left <= right);
	*sum = s;
	return left;
}

int quick_select_sum_k_sketch_with_priority_num(sketch_with_priority_num buff[], int i, int j, int sum)
{
	int pivotindex, m = i;
	int s, done = 0; // 左側に priority が小さいもので num の合計が目標の sum 以下で，確定した部分の num の合計
	while((pivotindex = find_pivot_sketch_with_priority_num(buff, i, j)) >= 0) {
		m = partition_by_pivot_sketch_with_priority_num(buff, i, j, buff[pivotindex].priority, &s);
		if(s + done > sum) {
			j = m - 1;
		} else if(s + done < sum) {
			i = m;
			done += s;
		} else {
			break;
		}
	}
	return j;
}

int partition_by_pivot_answer_mt(answer_type ans[], int i, int j, dist_type piv)
// ans[i], ... , ans[j] をそれらの dist と piv との大小によって分け， 
// piv より小さいものが ans[i], ... , ans[k - 1] に，   
// そうでないものが ans[k], ... , ans[j] に来るようにする. 
// 右側のリストの先頭の位置(k)を返す. 
{
	if(i > j) {
		return i; // 区間 i, ... , j が空なので，piv より小さいものは一つもないことになるので，すべて右側になると考える．
	}
	int left, right;
	answer_type temp;
	left = i;   right = j;
	do {
		while(ans[left].dist < piv && left < right) left++;
		if(left >= right) { 
			if(ans[left].dist < piv) {
				return left + 1; // break; */
			} else {
				return left;
			}
		}
		while(ans[right].dist >= piv && right > left) right--;
		if (left < right) { 
			temp = ans[left];
			ans[left] = ans[right];
			ans[right] = temp;
		}
	} while(left <= right);
	return left;
}

#ifdef _OPENMP
int quick_select_k_answer_mt(int nt, int k_th[], int size, answer_type ans[], int k)
{
	int left[nt], right[nt];
	answer_type *list = NULL;
	int t;
	int num_c = 0; // list で調べる範囲のデータ数
	for(t = 0; t < nt; t++) {
		left[t] = 0; right[t] = k_th[t] - 1;
		num_c += k_th[t];
	}
	int pivotindex;
	dist_type pivot = UINT_MAX;
	int s = 0, done = 0; // 左側に dist が小さいもので num の合計が目標の sum 以下で，確定した部分の num の合計

	while(1) {
		// まず，ピボットを求める．各スレッドが求めた list から探す．
		for(t = 0; t < nt; t++) {
			list = ans + t * size;
			if(left[t] >= right[t]) continue; // スレッドtが求めた内，範囲（配列の添え字）がleft[t]からright[t]までのものが1個か0個
			if((pivotindex = find_pivot_answer(list, left[t], right[t])) >= 0) {
				pivot = list[pivotindex].dist;
				break;
			}
		}
		if(t == nt) { // find_pivot_answerでは，ピボットが見つからなかった．各list内に一つしか値が無い．
/*
fprintf(stderr, "find_pivot_answer could not find pivot: done = %d\n", done);
for(t = 0; t < nt; t++) {
	list = ans + t * size;
	fprintf(stderr, "t = %d, left[%d] = %d, right[%d] = %d\n", t, t, left[t], t, right[t]);
	if(s > 30) continue;
	for(int j = left[t]; j <= right[t]; j++) {
		fprintf(stderr, "dist[%d] = %6d ", j, list[j].dist);
	}
	fprintf(stderr, "\n");
}
*/
			do {
				dist_type min_pivot = UINT_MAX;
				int min_t = -1;
				// list内の最小値を求める．スレッドmin_tが持っている値min_pivot
				for(t = 0; t < nt; t++) {
					if(left[t] > right[t]) continue; // listが空
					list = ans + t * size;
					if(list[left[t]].dist < min_pivot) {
						min_pivot = list[left[t]].dist;
						min_t = t;
					}
				}
				if(min_t == -1) {
					fprintf(stderr, "all list empty!\n"); getchar();
				} else {
					// min_tのlistを求めるものに加える．
					if(left[min_t] <= right[min_t]) {
						done += right[min_t] + 1 - left[min_t];
						left[min_t] = right[min_t] + 1;
					}
				}
			} while(done < k);
//fprintf(stderr, "after adding data: done = %d\n", done);
//for(t = 0; t < nt; t++) {
//	list = ans + t * size;
//	fprintf(stderr, "t = %d, left[%d] = %d, right[%d] = %d\n", t, t, left[t], t, right[t]);
//	if(s > 30) continue;
//	for(int j = left[t]; j <= right[t]; j++) {
//		fprintf(stderr, "dist[%d] = %6d ", j, list[j].dist);
//	}
//	fprintf(stderr, "\n");
//}
			for(t = 0; t < nt; t++) {
				k_th[t] = left[t];
//				k_th[t] = right[t] < 0 ? 0 : right[t];
//fprintf(stderr, "k_th[%d] = %d, ", t, k_th[t]);
			}
//fprintf(stderr, "\n");
//			fprintf(stderr, "(0) pivot = %d, done = %d\n", pivot, done);
/*
			for(int s = 0; s < nt; s++) {
				fprintf(stderr, "k_th[%d] = %d\n", s, k_th[s]);
				answer_type *list = ans + s * size;
				for(int j = 0; j < size; j++) {
					if(j < k_th[s] && list[j].dist >= pivot) {
						fprintf(stderr, "(a) bad partition: list[%d].dist = %d\n", j, list[j].dist);
//						break;
					} else if(j >= k_th[s] && list[j].dist < pivot) {
						fprintf(stderr, "(b) bad partition: list[%d].dist = %d\n", j, list[j].dist);
						if(list[j].dist == 0) break;
					}
				}
			}
			getchar();
*/

// fprintf(stderr, "HERE(0): done = %d\n", done); getchar();
			return done;
		}
		int m[nt];
//		fprintf(stderr, "done = %d\n", done);
/*
		s = 0;
		for(t = 0; t < nt; t++) {
			fprintf(stderr, "left[%d] = %d, right[%d] = %d, k_th[%d] = %d\n", t, left[t], t, right[t], t, k_th[t]);
			s += left[t];
		}
		fprintf(stderr, "s = %d\n", s);
*/
		s = 0;
#define MIN_LIST 1000
		if(num_c > MIN_LIST) {
			#pragma omp parallel private(t, list) reduction(+:s)
			{
				t = omp_get_thread_num();
				list = ans + t * size;
				m[t] = partition_by_pivot_answer_mt(list, left[t], right[t], pivot);
				s = m[t] - left[t];
			}
		} else {
// fprintf(stderr, "num (= %d) < MIN_LIST\n", num_c); getchar();
			for(int t = 0; t < nt; t++)
			{
				list = ans + t * size;
				if(left[t] >= right[t]) {
					m[t] = list[left[t]].dist < pivot ? left[t] + 1 : left[t];
				} else {
					m[t] = partition_by_pivot_answer_mt(list, left[t], right[t], pivot);
				}
				s += m[t] - left[t];
			}
		}
//		fprintf(stderr, "after partition: pivot = %d, s = %d, done = %d\n", pivot, s, done);
//		for(t = 0; t < nt; t++) {
//			list = ans + t * size;
//			fprintf(stderr, "t = %d, m[%d] = %d, left[%d] = %d, right[%d] = %d\n", t, t, m[t], t, left[t], t, right[t]);
//			if(s > 30) continue;
//			for(int j = left[t]; j <= right[t]; j++) {
//				fprintf(stderr, "%6d", list[j].dist);
//			}
//			fprintf(stderr, "\n");
//		}
/*		for(t = 0; t < nt; t++) {
			int ea = 0, eb = 0;
			list = ans + t * size;
			for(int j = 0; j < m[t]; j++) {
				if(list[j].dist >= pivot) {
					ea++;
				}
			}
			for(int j = m[t]; j < k_th[t]; j++) {
				if(list[j].dist < pivot) {
					fprintf(stderr, "error: t = %d, m[%d] = %d, list[%d].dist = %d, pivot = %d", t, t, m[t], j, list[j].dist, pivot);
					eb++;
				}
			}
			fprintf(stderr, "t = %d, ea = %d, eb = %d\n", t, ea, eb);
		}
		getchar();
*/
		if(s + done > k) {
			num_c = 0;
			for(t = 0; t < nt; t++) {
				right[t] = m[t] - 1;
				num_c += left[t] <= right[t] ? right[t] + 1 - left[t] : 0;
			}
/*
			fprintf(stderr, "s is too large\n");
			for(t = 0; t < nt; t++) {
				list = ans + t * size;
				fprintf(stderr, "t = %d, m[%d] = %d, left[%d] = %d, right[%d] = %d\n", t, t, m[t], t, left[t], t, right[t]);
				if(s > 30) continue;
				for(int j = left[t]; j <= right[t]; j++) {
					fprintf(stderr, "%6d", list[j].dist);
				}
				fprintf(stderr, "\n");
			}
*/
		} else if(s + done < k) {
			num_c = 0;
			for(t = 0; t < nt; t++) {
				left[t] = m[t];
				num_c += left[t] <= right[t] ? right[t] + 1 - left[t] : 0;
			}
			done += s;
		} else {
//			for(t = 0; t < nt; t++) {
//				k_th[t] = right[t] < 0 ? 0 : right[t] + 1;
//			}
			done += s;
			for(t = 0; t < nt; t++) {
				k_th[t] = m[t];
			}
/*
			fprintf(stderr, "(1) pivot = %d, done = %d\n", pivot, done);
			for(int s = 0; s < nt; s++) {
				fprintf(stderr, "k_th[%d] = %d\n", s, k_th[s]);
				answer_type *list = ans + s * size;
				for(int j = 0; j < size; j++) {
					if(j < k_th[s] && list[j].dist >= pivot) {
						fprintf(stderr, "(a) bad partition: list[%d].dist = %d\n", j, list[j].dist);
						break;
					} else if(j >= k_th[s] && list[j].dist < pivot) {
						fprintf(stderr, "(b) bad partition: list[%d].dist = %d\n", j, list[j].dist);
						break;
					}
				}
			}
			getchar();
*/
// fprintf(stderr, "HERE(1): done = %d\n", done); getchar();
			return done;
		}
	}
	int selected;
	for(t = 0; t < nt; t++) {
		k_th[t] = right[t] < 0 ? 0 : right[t];
		selected += k_th[t];
	}
/*
	fprintf(stderr, "(2) pivot = %d, done = %d\n", pivot, done);
	for(int s = 0; s < nt; s++) {
		fprintf(stderr, "k_th[%d] = %d\n", s, k_th[s]);
		answer_type *list = ans + s * size;
		for(int j = 0; j < size; j++) {
			if(j < k_th[s] && list[j].dist >= pivot) {
				fprintf(stderr, "(a) bad partition: list[%d].dist = %d\n", j, list[j].dist);
				break;
			} else if(j >= k_th[s] && list[j].dist < pivot) {
				fprintf(stderr, "(b) bad partition: list[%d].dist = %d\n", j, list[j].dist);
				break;
			}
		}
	}
	getchar();
*/
//	return done;
//	fprintf(stderr, "HERE(2): selected = %d, done = %d\n", selected, done); getchar();
	return selected;
}
#endif

int balance_interval_list(interval_list *ivl)
{
	int nt = ivl->nt;
	interval *list = NULL;
	int sum_part[nt];
	int sum_all = 0;
	for(int t = 0; t < nt; t++) {
		list = ivl->list + t * ivl->size;
		sum_part[t] = 0;
		for(int j = 0; j < ivl->lg[t]; j++) {
			sum_part[t] += list[j].run;
		}
		sum_all += sum_part[t];
	}
	// 並列して選択すると，スレッドで選択してデータ数が不揃いになる可能性がある．
	// そのまま，第2段階検索を行うと，スレッドの処理時間にばらつきが出て，早く終わったスレッドが休み，遅くなるスレッドの終了を待つので，全体として遅くなる．
	// それを防ぐために，データ数ができるだけ均一になるように，区間を移動してバランスを取っておく．
	// 多すぎるスレッドから，少なすぎるスレッドに区間を移動する．
	// リストのデータ数で sort してからバランスを取る
	int ave_num = sum_all / nt;
	// sum_part[0], ... , sum_part[nt - 1] を idx を用いて降順に相対ソートする．
	// sum_part[idx[0]] >= sum_part[idx[1]] >= ... >= sum_part[ifx[nt - 1]] となるようにする．
	// 最初は，すべての i に対して，idx[i] = i としておく．
	int idx[nt], movable[nt];
	for(int t = 0; t < nt; t++) {
		idx[t] = t;
		movable[t] = 1;
	}
	for(int i = 1; i < nt; i++) {
		int temp = idx[i], j;
		for(j = i; j > 0 && sum_part[idx[j - 1]] < sum_part[temp]; j--) {
			idx[j] = idx[j - 1];
		}
		idx[j] = temp;
	}

	// 移動するものを集めておいて，まとめて移動する．∵大き過ぎて移動できない区間があると，移動するものが飛び飛びになるので．
	while(1) {
		// sum_part[t] が最大のものが ave_num を超えていて移動可能な interval を持っているなら，移動を試みる．
		int max_t = -1;
		for(int i = 0; i < nt; i++) {
			max_t = idx[i];
			if(ivl->lg[max_t] == 1) movable[max_t] = 0;
			if(movable[max_t]) break;
		}
		if(!movable[max_t]) {
			// fprintf(stderr, "no more movable list.\n");
			break;
		}
		interval *max_list = ivl->list + max_t * ivl->size;
		if(sum_part[max_t] < 1.01 * ave_num) break; // データ数が少ないので，終了．
		// リストの後部から順に移動できるもの後部（左側）に集める．
		int min_t = idx[nt - 1];
		interval *min_list = ivl->list + min_t * ivl->size;
		if(sum_part[min_t] > 0.99 * ave_num) break; // 不足するデータ数が少ないものしかないので，終了．
		int m = ivl->lg[max_t] - 1, sum_move = 0; 
		static int *not_movable = NULL, ivl_size = 0;
		int num_not_movable = 0;
		if(not_movable == NULL) {
			not_movable = MALLOC(sizeof(int) * ivl->size);
			ivl_size = ivl->size;
		} else if(ivl->size > ivl_size) {
			FREE(not_movable, sizeof(int) * ivl_size);
			not_movable = MALLOC(sizeof(int) * ivl->size);
			ivl_size = ivl->size;
		}
		while(m > 0 && sum_part[max_t] - sum_move >= 1.01 * ave_num) {
			while(m >= 0 && sum_part[max_t] - sum_move - max_list[m].run < ave_num) {
				// m 番目の interval の run が長過ぎる（移動すると max_t の run の合計が ave_num を切ってしまう．
				not_movable[num_not_movable++] = m--; // 移動できないものを覚えておいて，次（m--）を探す．
			}
			// 上のループを抜けたときは，m = -1 または m >= 0 && sum_part[max_t] - sum_move - max_list[m].run >= ave_num
			if(m < 0) break; // 右端まで到達して，移動できるものが見つからなかった．
			// m >= 0 であれば，m 番目のものは移動できる．
			sum_move += max_list[m--].run;
			if(sum_part[min_t] + sum_move >= ave_num) {
				// ここまでで見つけた移動可能なものを min_t に移動すると，初めて min_t が ave_num 以上になるので，ここまでを移してしまう．
				// そのために，ループを脱出する．
				break;
			}
		}
		if(sum_move == 0) { // 移動できるものが見つからなかった．∵どれを移動しても短くなり過ぎる．つまり，ここからは移動不可．
			movable[max_t] = 0;
			break;
		}
		m++;
		while(num_not_movable) {
			// 途中に移動できないものがあるので，それを右に移動して，左側に移動できるものを集める．
			// m のものは移動できる．not_movable[i] (i = 0, ... , num_not_movable - 1) のものは移動できない．
			interval temp = max_list[m];
			max_list[m] = max_list[not_movable[--num_not_movable]];
			max_list[not_movable[num_not_movable]] = temp;
			m++; // m に移動できないものを置いたので，m を左に．
			// まだ，移動できないものが右側にあるならば，m には移動できるものがあるはず．
		}
		// ここでは，m から ivl->lg[max_t] - 1 までに移動できるものが連続して集められている．
		// run の合計は sum_move であり，それを移動しても，残りの run の合計は ave_num 以上．
		// 集めたものを min_t に移す．
		memcpy(min_list + ivl->lg[min_t], max_list + m, sizeof(interval) * (ivl->lg[max_t] - m));
		sum_part[max_t] -= sum_move;
		sum_part[min_t] += sum_move;
		ivl->lg[min_t] += ivl->lg[max_t] - m;
		ivl->lg[max_t] = m;
		int j;
		for(j = 0; j < nt - 1 && sum_part[idx[j + 1]] > sum_part[max_t]; j++) idx[j] = idx[j + 1];
		idx[j] = max_t;
		for(j = nt - 1; j > 0 && sum_part[idx[j - 1]] < sum_part[min_t]; j--) idx[j] = idx[j - 1];
		idx[j] = min_t;
	}

	return sum_all;
}

#define KNN_BUFFER_FACTOR 3 // kNN_buffer の buff の大きさを k の何倍にするかを指定する．
#define KNN_BUFFER_FACTOR2 2 // タイ（同点）のものをここまで残す．

kNN_buffer *new_kNN_buffer(int k)
{
	kNN_buffer *b = MALLOC(sizeof(kNN_buffer));

	b->k = k;
	b->buff = MALLOC(sizeof(answer_type) * KNN_BUFFER_FACTOR * k); 
	if(b == NULL) {return b;}
	b->num = 0;
	b->k_nearest = INT_MAX;
	return b;
}

void free_kNN_buffer(kNN_buffer *b)
{
//	fprintf(stderr, "free_KNN_buffer: FREE = %lu + %lu\n", sizeof(kNN_buffer), sizeof(answer_type) * KNN_BUFFER_FACTOR * b->k);
	FREE(b->buff, sizeof(answer_type) * KNN_BUFFER_FACTOR * b->k);
	FREE(b, sizeof(kNN_buffer));
	b = NULL;
}

dist_type make_empty_kNN_buffer(kNN_buffer *b)
{
	b->num = 0;
	b->k_nearest = INT_MAX;
	return b->k_nearest;
}

dist_type push_kNN_buffer(answer_type *a, kNN_buffer *b)
{
	dist_type d;
	if(b->num >= KNN_BUFFER_FACTOR * b->k) {
		d = flush_kNN_buffer(b);
	} else {
		d = b->k_nearest;
	}
	if(a->dist > d) {return d;} // k_nearest より遠い解はバッファーには入れない
	b->buff[b->num] = *a;
	b->num++;
	return d;
}

/*
#ifndef PUBMED23
int comp_answer(const void *a, const void *b) {
	if(((answer_type *) a) -> dist < (((answer_type *) b) -> dist))
		return -1;
	else if(((answer_type *) a) -> dist == (((answer_type *) b) -> dist))
		return 0;
	else
		return 1;
}
#else
*/
int comp_answer(const void *a, const void *b) {
	if(((answer_type *) a) -> dist < (((answer_type *) b) -> dist)) {
		return -1;
	} else if(((answer_type *) a) -> dist == (((answer_type *) b) -> dist)) {
		if(((answer_type *) a) -> data_num < (((answer_type *) b) -> data_num)) {
			return -1;
		} else if(((answer_type *) a) -> data_num == (((answer_type *) b) -> data_num)) {
			return 0;
		} else {
			return 1;
		}
	} else {
		return 1;
	}
}
//#endif

dist_type flush_kNN_buffer(kNN_buffer *b)
{
	if(b->num == 0) return b->k_nearest;
//	if(b->num <= b->k) { fprintf(stderr, "kNN_buff: num <= k, num = %d, k = %d\n", b->num, b->k); }
	if(b->num > b->k) {
		quick_select_k_answer(b->buff, 0, b->num - 1, b->k);
//		qsort(b->buff, b->num, sizeof(answer_type), comp_answer);
		b->num = b->k;
	}
	b->k_nearest = b->buff[0].dist;
	for(int i = 1; i < b->k; i++) if(b->k_nearest < b->buff[i].dist) b->k_nearest = b->buff[i].dist;
	return b->k_nearest;
}

dist_type final_flush_kNN_buffer(kNN_buffer *b)
{
	qsort(b->buff, b->num, sizeof(answer_type), comp_answer);
	if(b->num > b->k) b->num = b->k;
	b->k_nearest = b->buff[b->num - 1].dist;
	return b->k_nearest;
}

dist_type merge_kNN_buffer(kNN_buffer *b, kNN_buffer *pool[], int n)
{
	int i, k = b->k, m;
	answer_type *ans = (answer_type *)malloc(sizeof(answer_type) * k * n);
	
	for(i = m = 0; m < n; m++) {
		for(int j = 0; j < (pool[m]->num < pool[m]->k ? pool[m]->num : pool[m]->k); j++) {
			ans[i++] = pool[m]->buff[j];
		}
	}
	quick_select_k_answer(ans, 0, i - 1, k);
	quick_sort_answer(ans, 0, k - 1);
	for(i = 0; i < k; i++) {
		b->buff[i] = ans[i];
	}
	b->k_nearest = b->buff[k - 1].dist;
	free(ans);
	return b->k_nearest;
}

answer_type *get_kNN_buffer(kNN_buffer *b, int i)
{
	return &(b->buff[i]);
}

// （現状）最近傍解のみ取り込んでいる．→100NNまでcsvに保存されている
answer_type *read_correct_answer(char *answer_csv_file, int num_queries)
{
	FILE *afp ;
    char buf[100000]={0};
	int i, q;
	answer_type *ans = (answer_type *)malloc(sizeof(answer_type) * num_queries);
	
	afp=fopen(answer_csv_file,"r");
	if(afp == NULL){
		fprintf(stderr, "cannot open answer csv file = %s\n", answer_csv_file);
		exit(0);
	}

	// 1行読み捨てる
	if(fgets(buf, MAX_LEN, afp) == NULL) {
		fprintf(stderr, "read error: answer csv file = %s\n", answer_csv_file);
		exit(0);
	}
	// もう1行読み捨てる
	if(fgets(buf, MAX_LEN, afp) == NULL) {
		fprintf(stderr, "read error: answer csv file = %s\n", answer_csv_file);
		exit(0);
	}
	// answer index, dist を読み込む
    for(i = 0; fgets(buf, MAX_LEN, afp) != NULL && i < num_queries; i++){
		q = atoi(strtok(buf, ","));
    	if(q != i) {
    		fprintf(stderr, "invalid answer file (line number != query number): answer csv file = %s\n", answer_csv_file);
			exit(0);
    	}
    	strtok(NULL, ","); // auxi_type query_id = atol(strtok(NULL, ","));
		ans[i].data_num = atoi(strtok(NULL, ","));
		ans[i].dist = atoi(strtok(NULL, ","));
		#ifdef ANSWER_BY_DATA_ID
		ans[i].data_num[k] = atoi(strtok(NULL, ","));
		#endif
}
	fclose(afp); 

	if(i != num_queries) { 
		fprintf(stderr, "query_num != answer_num\n");
		exit(0);
	}
    
	return ans;
}

#ifdef NUM_NN
// （改訂）NUM_NN が定義されていたら，NUM_NN近傍まで読み込む
answer_type_NN *read_correct_answer_NN(char *answer_csv_file, int num_queries)
{
	FILE *afp ;
    char buf[100000]={0};
	int i, q;
	answer_type_NN *ans = (answer_type_NN *)malloc(sizeof(answer_type_NN) * num_queries);
	
	afp=fopen(answer_csv_file,"r");
	if(afp == NULL){
		fprintf(stderr, "cannot open answer csv file = %s\n", answer_csv_file);
		exit(0);
	}

	// 1行読み捨てる
	if(fgets(buf, MAX_LEN, afp) == NULL) {
		fprintf(stderr, "read error: answer csv file = %s\n", answer_csv_file);
		exit(0);
	}
	// もう1行読み捨てる
	if(fgets(buf, MAX_LEN, afp) == NULL) {
		fprintf(stderr, "read error: answer csv file = %s\n", answer_csv_file);
		exit(0);
	}
	// answer index, dist を読み込む
    for(i = 0; fgets(buf, MAX_LEN, afp) != NULL && i < num_queries; i++){
		q = atoi(strtok(buf, ","));
    	if(q != i) {
    		fprintf(stderr, "invalid answer file (line number != query number): answer csv file = %s\n", answer_csv_file);
			exit(0);
    	}
		for(int k = 0; k < NUM_NN; k++) {
	    	strtok(NULL, ","); // auxi_type query_id = atol(strtok(NULL, ","));
			ans[i].data_num[k] = atoi(strtok(NULL, ","));
			#ifndef ANSWER_DIST_FLOAT
			ans[i].dist[k] = atoi(strtok(NULL, ","));
			#else
			ans[i].dist[k] = atof(strtok(NULL, ","));
			#endif
			#ifdef ANSWER_BY_DATA_ID
			ans[i].data_num[k] = atoi(strtok(NULL, ","));
			#endif
		}
    }
	fclose(afp); 

	if(i != num_queries) { 
		fprintf(stderr, "query_num != answer_num\n");
		exit(0);
	}
    
	return ans;
}
#endif

void sort_answer(int num, answer_type ans[])
{
	qsort(ans, num, sizeof(answer_type), comp_answer);
}

// 解候補から top-K を求める（第２段階検索）
void search_kNN(dataset_handle *dh, query_type *qr, int num_candidates, int data_num_of_candidate[], kNN_buffer *top_k)
{
//	fprintf(stderr, "search_kNN\n");
	#if defined(_OPENMP) && NUM_THREADS > 1
	omp_set_num_threads(num_candidates < NUM_THREADS ? num_candidates : NUM_THREADS);
	int num_k = top_k->k;
	int nt = omp_get_max_threads(); 	// スレッド数を求める
	static kNN_buffer *b_pool[NUM_THREADS] = {NULL};
	dist_type k_nearest[nt];
	if(b_pool[0] == NULL) {
		for(int t = 0; t < nt; t++) {
			if((b_pool[t] = new_kNN_buffer(num_k)) == NULL) {
				fprintf(stderr, "cannot allocate new kNN buffer\n");
				exit(0);
			}
		}
	}
	#else
	dist_type k_nearest;
	#endif

	SET_DIST(qr->ftr); // 距離はすべて質問からになるので，片方を質問に固定
	#if defined(_OPENMP) && NUM_THREADS > 1
	for(int t = 0; t < nt; t++) {
		k_nearest[t] = make_empty_kNN_buffer(b_pool[t]);
		if(dh->ftr_on == SECONDARY_MEMORY) dh->mf[t]->read_in = 0;
	}
	#pragma omp parallel for
	#else
	k_nearest = make_empty_kNN_buffer(top_k);
	if(dh->ftr_on == SECONDARY_MEMORY) dh->mf[0]->read_in = 0;
	#endif
	for(int i = 0; i < num_candidates; i++) {
		struct_multi_ftr *mf = NULL;
		#if defined(_OPENMP) && NUM_THREADS > 1
		int t = omp_get_thread_num();
		kNN_buffer *b = b_pool[t];
		if(dh->ftr_on == SECONDARY_MEMORY) mf = dh->mf[t];
		#else
		kNN_buffer *b = top_k;
		if(dh->ftr_on == SECONDARY_MEMORY) mf = dh->mf[0];
		#endif
		answer_type ans;
		ans.data_num = data_num_of_candidate[i];
		if(dh->ftr_on == MAIN_MEMORY) {
			ans.dist = DISTANCE_22(dh->ds->ftr_id[data_num_of_candidate[i]].ftr);
			#ifdef ANSWER_WITH_DATA_ID
			ans.data_id = dh->ds->ftr_id[data_num_of_candidate[i]].data_id;
			#endif
		} else {
			struct_ftr_id *ftr_id_p;
			ftr_id_p = get_next_ftr_id_from_multi_ftr(mf, data_num_of_candidate, i, num_candidates);
			ans.dist = DISTANCE_22(ftr_id_p->ftr);
			#ifdef ANSWER_WITH_DATA_ID
			ans.data_id = ftr_id_p->data_id;
			#endif
		}
		#if defined(_OPENMP) && NUM_THREADS > 1
		if(ans.dist < k_nearest[t]) {
			k_nearest[t] = push_kNN_buffer(&ans, b);
		}
		#else
		if(ans.dist < k_nearest) {
			k_nearest = push_kNN_buffer(&ans, b);
		}
		#endif
	}

	#if defined(_OPENMP) && NUM_THREADS > 1
		#pragma omp parallel for
		for(int t = 0; t < nt; t++) {
			k_nearest[t] = flush_kNN_buffer(b_pool[t]);
		}
		merge_kNN_buffer(top_k, b_pool, nt);
	#else
		k_nearest = final_flush_kNN_buffer(top_k);
	#endif
}

#ifdef STATIC_KNN_BUFFER_FOR_SEARCH

static kNN_buffer *b_pool[NUM_THREADS] = {NULL};

void init_search_kNN_on_ram(int num_k) 
{
	if(b_pool[0] == NULL) {
		for(int t = 0; t < NUM_THREADS; t++) {
			if((b_pool[t] = new_kNN_buffer(num_k)) == NULL) {
				fprintf(stderr, "cannot allocate new kNN buffer\n");
				exit(0);
			}
			for(int i = 0; i < num_k; i++) {
				b_pool[t]->buff[i].data_num = 0;
			}
		}
	}
}

// 解候補から top-K を求める（第２段階検索）特徴データは配列に読み込んだものを用いる
void search_kNN_on_ram(struct_ftr_id ftr_id[], query_type *qr, int num_candidates, int data_num_of_candidate[], kNN_buffer *top_k)
{
	#if defined(_OPENMP) && NUM_THREADS > 1
	omp_set_num_threads(num_candidates < NUM_THREADS ? num_candidates : NUM_THREADS);
	int num_k = top_k->k;
	int nt = omp_get_max_threads(); 	// スレッド数を求める
	for(int t = 0; t < nt; t++) {
		make_empty_kNN_buffer(b_pool[t]);
	}
	#else
	make_empty_kNN_buffer(top_k);
	#endif

	SET_DIST(qr->ftr); // 距離はすべて質問からになるので，片方を質問に固定
	#if defined(_OPENMP) && NUM_THREADS > 1
	#pragma omp parallel for
	#endif
	for(int i = 0; i < num_candidates; i++) {
		int k_nearest = INT_MAX;
		#if defined(_OPENMP) && NUM_THREADS > 1
		int t = omp_get_thread_num();
		kNN_buffer *b = b_pool[t];
		#else
		kNN_buffer *b = top_k;
		#endif
		answer_type ans;
		ans.data_num = data_num_of_candidate[i];
		ans.dist = DISTANCE_22(ftr_id[data_num_of_candidate[i]].ftr);
		#ifdef ANSWER_WITH_DATA_ID
		ans.data_id = 0;
		#endif
		if(ans.dist < k_nearest) {
			k_nearest = push_kNN_buffer(&ans, b);
		}
	}

	#if defined(_OPENMP) && NUM_THREADS > 1
		#pragma omp parallel for
		for(int t = 0; t < nt; t++) {
			flush_kNN_buffer(b_pool[t]);
		}
		merge_kNN_buffer(top_k, b_pool, nt);
	#else
		k_nearest = final_flush_kNN_buffer(top_k);
	#endif
}
#else
// 解候補から top-K を求める（第２段階検索）特徴データは配列に読み込んだものを用いる
void search_kNN_on_ram(struct_ftr_id ftr_id[], query_type *qr, int num_candidates, int data_num_of_candidate[], kNN_buffer *top_k)
{
	#if defined(_OPENMP) && NUM_THREADS > 1
	omp_set_num_threads(num_candidates < NUM_THREADS ? num_candidates : NUM_THREADS);
	int num_k = top_k->k;
	int nt = omp_get_max_threads(); 	// スレッド数を求める
	kNN_buffer *b_pool[nt];
	for(int t = 0; t < nt; t++) {
		if((b_pool[t] = new_kNN_buffer(num_k)) == NULL) {
			fprintf(stderr, "cannot allocate new kNN buffer\n");
			exit(0);
		}
	}
	#else
	make_empty_kNN_buffer(top_k);
	#endif

	SET_DIST(qr->ftr); // 距離はすべて質問からになるので，片方を質問に固定
	#if defined(_OPENMP) && NUM_THREADS > 1
	#pragma omp parallel for
	#endif
	for(int i = 0; i < num_candidates; i++) {
		int k_nearest = INT_MAX;
		#if defined(_OPENMP) && NUM_THREADS > 1
		int t = omp_get_thread_num();
		kNN_buffer *b = b_pool[t];
		#else
		kNN_buffer *b = top_k;
		#endif
		answer_type ans;
		ans.data_num = data_num_of_candidate[i];
		ans.dist = DISTANCE_22(ftr_id[data_num_of_candidate[i]].ftr);
		#ifdef ANSWER_WITH_DATA_ID
		ans.data_id = 0;
		#endif
		if(ans.dist < k_nearest) {
			k_nearest = push_kNN_buffer(&ans, b);
		}
	}

	#if defined(_OPENMP) && NUM_THREADS > 1
		#pragma omp parallel for
		for(int t = 0; t < nt; t++) {
			flush_kNN_buffer(b_pool[t]);
		}
		merge_kNN_buffer(top_k, b_pool, nt);
		for(int t = 0; t < nt; t++) {
			free_kNN_buffer(b_pool[t]);
		}
	#else
		final_flush_kNN_buffer(top_k);
	#endif
}
#endif

// answer_check と同様．ただし，解候補をデータ番号のリストではなく，区間リストで与える．
// 区間リストからデータ番号を求めるために，バケットを与える．
int answer_check_interval(struct_bucket *bucket, answer_type *ans, interval_list *ivl, kNN_buffer *top_k)
{
//	int nt = ivl->nt;
	int size = ivl->size;
	int *lg = ivl->lg;
	int found = 0;

	#if defined(_OPENMP) && NUM_THREADS > 1
	#ifndef THREAD_PLUS
		omp_set_num_threads(ivl->nt);
		#pragma omp parallel reduction (max: found)
		{
			int t = omp_get_thread_num(); // スレッド番号の取得
			found = 0;
			for(int i = 0; i < lg[t] && !found; i++) {
				interval *list = ivl->list + t * size;
				#ifdef DATA_NUM_IN_SKETCH_ORDER
					if(ans->data_num >= list[i].start && ans->data_num <= list[i].start + list[i].run - 1) {
						found = 1;
					}
				#else
					for(int j = list[i].start; j < list[i].start + list[i].run; j++) {
						int data_num_of_candidate = bucket->idx[j];
						if(ans->data_num == data_num_of_candidate) {
							found = 1;
							break;
						}
					}
				#endif
			}
		}
	#else
		omp_set_num_threads(NUM_THREADS);
		#pragma omp parallel reduction (max: found)
		{
			int tn = omp_get_thread_num(); // スレッド番号の取得
			found = 0;
			for(int t = tn * (1 << THREAD_PLUS); t < (tn + 1) * (1 << THREAD_PLUS); t++) {
				for(int i = 0; i < lg[t] && !found; i++) {
					interval *list = ivl->list + t * size;
					#ifdef DATA_NUM_IN_SKETCH_ORDER
						if(ans->data_num >= list[i].start && ans->data_num <= list[i].start + list[i].run - 1) {
							found = 1;
						}
					#else
						for(int j = list[i].start; j < list[i].start + list[i].run; j++) {
							int data_num_of_candidate = bucket->idx[j];
							if(ans->data_num == data_num_of_candidate) {
								found = 1;
								break;
							}
						}
					#endif
				}
			}
		}
	#endif
	#else
		for(int t = 0; t < ivl->nt && !found; t++) {
			for(int i = 0; i < lg[t] && !found; i++) {
				interval *list = ivl->list + t * size;
				#ifdef DATA_NUM_IN_SKETCH_ORDER
					if(ans->data_num >= list[i].start && ans->data_num <= list[i].start + list[i].run - 1) {
						found = 1;
						break;
					}
				#else
					for(int j = list[i].start; j < list[i].start + list[i].run; j++) {
						int data_num_of_candidate = bucket->idx[j];
						if(ans->data_num == data_num_of_candidate) {
							found = 1;
							break;
						}
					}
				#endif
			}
		}
	#endif
	
	top_k->num = 1;
	top_k->k_nearest = top_k->buff[0].dist;
	if(found) {
		top_k->buff[0].data_num = ans->data_num;
		top_k->buff[0].dist = ans->dist;
		return 1;
	} else {
		top_k->buff[0].data_num = 0;
		top_k->buff[0].dist = UINT_MAX;
		return 0;
	}
}

// 解放補のデータ番号 dana_bun_of_candidate[i] (0 <= i < num_candidates) の中に
// 正解のデータ番号 ans->data_num と一致するものがあるかどうかを調べて，
// 一致するものがあれば，top_k の先頭にデータ番号と正解の距離を与えて（実際に距離を計算していない），1 を返す．
// そうでなければ，top_k の先頭に距離の最大値 UINT_MAX を設定して，0 を返す．
int answer_check(answer_type *ans, int num_candidates, int data_num_of_candidate[], kNN_buffer *top_k)
{
	for(int i = 0; i < num_candidates; i++) {
		if(ans->data_num == data_num_of_candidate[i]) {
			top_k->buff[0].data_num = ans->data_num;
			top_k->buff[0].dist = ans->dist;
			top_k->num = 1;
			top_k->k_nearest = top_k->buff[0].dist;
			return 1;
		}
	}
	top_k->buff[0].data_num = data_num_of_candidate[0];
	top_k->buff[0].dist = UINT_MAX;
	top_k->num = 1;
	top_k->k_nearest = top_k->buff[0].dist;
	return 0;
}

void search_NN(dataset_handle *dh, query_type *qr, int num_candidates, int data_num_of_candidate[], kNN_buffer *top_k) // 解候補から top-K を求める（第２段階検索）
{
	#if defined(_OPENMP) && NUM_THREADS > 1
	omp_set_num_threads(NUM_THREADS);
	int nt = omp_get_max_threads(); 	// スレッド数を求める
	answer_type ans[nt];
	dist_type k_nearest[nt];
	#else
	answer_type ans;
	dist_type k_nearest;
	#endif

	SET_DIST(qr->ftr); // 距離はすべて質問からになるので，片方を質問に固定
	#if defined(_OPENMP) && NUM_THREADS > 1
	for(int t = 0; t < nt; t++) {
		k_nearest[t] = INT_MAX;
		if(dh->ftr_on != MAIN_MEMORY) {
			dh->mf[t]->read_in = 0;
		}
		ans[t].dist = INT_MAX;
	}
	#pragma omp parallel for
	#else
	k_nearest = INT_MAX;
	if(dh->ftr_on != MAIN_MEMORY) {
		dh->mf[0]->read_in = 0;
	}
	#endif
	for(int i = 0; i < num_candidates; i++) {
		#if defined(_OPENMP) && NUM_THREADS > 1
		int t = omp_get_thread_num();
		answer_type *a = &ans[t];
		struct_multi_ftr *mf = dh->mf[t];
		#else
		answer_type *a = &ans;
		struct_multi_ftr *mf = dh->mf[0];
		#endif
		int dist;
		if(dh->ftr_on == MAIN_MEMORY) {
			dist = DISTANCE_22(dh->ds->ftr_id[data_num_of_candidate[i]].ftr);
		} else {
			dist = DISTANCE_22(get_next_ftr_id_from_multi_ftr(mf, data_num_of_candidate, i, num_candidates)->ftr);
		}
		#if defined(_OPENMP) && NUM_THREADS > 1
		if(dist < k_nearest[t]) {
			a->data_num = data_num_of_candidate[i];
			a->dist = dist;
			k_nearest[t] = dist;
		}
		#else
		if(dist < k_nearest) {
			a->data_num = data_num_of_candidate[i];
			a->dist = dist;
			k_nearest = dist;
		}
		#endif
	}
	#if defined(_OPENMP) && NUM_THREADS > 1
	// 各スレッドが見つけた暫定解から距離最小のものを選ぶ
	int min_t = 0;
	dist_type min_dist = UINT_MAX;
	for(int t = 0; t < nt; t++) {
		if(ans[t].dist < min_dist) {
			min_dist = ans[t].dist;
			min_t = t;
		}
	}
	top_k->buff[0] = ans[min_t];
	#else
	top_k->buff[0] = ans;
	#endif
	top_k->num = 1;
	top_k->k_nearest = top_k->buff[0].dist;
}

void search_NN_on_ram(struct_ftr_id ftr_id[], query_type *qr, int num_candidates, int data_num_of_candidate[], kNN_buffer *top_k) // 解候補から top-K を求める（第２段階検索）
{
	#if defined(_OPENMP) && NUM_THREADS > 1
	omp_set_num_threads(NUM_THREADS);
	int nt = omp_get_max_threads(); 	// スレッド数を求める
	answer_type ans[nt];
	dist_type k_nearest[nt];
	#else
	answer_type ans;
	dist_type k_nearest;
	#endif

	SET_DIST(qr->ftr); // 距離はすべて質問からになるので，片方を質問に固定
	#if defined(_OPENMP) && NUM_THREADS > 1
	for(int t = 0; t < nt; t++) {
		k_nearest[t] = INT_MAX;
//		if(dh->ftr_on != MAIN_MEMORY) {
//			dh->mf[t]->read_in = 0;
//		}
		ans[t].dist = INT_MAX;
	}
	#pragma omp parallel for
	#else
	k_nearest = INT_MAX;
//	if(dh->ftr_on != MAIN_MEMORY) {
//		dh->mf[0]->read_in = 0;
//	}
	#endif
	for(int i = 0; i < num_candidates; i++) {
		#if defined(_OPENMP) && NUM_THREADS > 1
		int t = omp_get_thread_num();
		answer_type *a = &ans[t];
//		struct_multi_ftr *mf = dh->mf[t];
		#else
		answer_type *a = &ans;
//		struct_multi_ftr *mf = dh->mf[0];
		#endif
		int dist;
		dist = DISTANCE_22(ftr_id[data_num_of_candidate[i]].ftr);
//		if(dh->ftr_on == MAIN_MEMORY) {
//			dist = DISTANCE_22(dh->ds->ftr_id[data_num_of_candidate[i]].ftr);
//		} else {
//			dist = DISTANCE_22(get_next_ftr_id_from_multi_ftr(mf, data_num_of_candidate, i, num_candidates)->ftr);
//		}
		#if defined(_OPENMP) && NUM_THREADS > 1
		if(dist < k_nearest[t]) {
			a->data_num = data_num_of_candidate[i];
			a->dist = dist;
			k_nearest[t] = dist;
		}
		#else
		if(dist < k_nearest) {
			a->data_num = data_num_of_candidate[i];
			a->dist = dist;
			k_nearest = dist;
		}
		#endif
	}
	#if defined(_OPENMP) && NUM_THREADS > 1
	// 各スレッドが見つけた暫定解から距離最小のものを選ぶ
	int min_t = 0;
	dist_type min_dist = UINT_MAX;
	for(int t = 0; t < nt; t++) {
		if(ans[t].dist < min_dist) {
			min_dist = ans[t].dist;
			min_t = t;
		}
	}
	top_k->buff[0] = ans[min_t];
	#else
	top_k->buff[0] = ans;
	#endif
	top_k->num = 1;
	top_k->k_nearest = top_k->buff[0].dist;
}

void out_result(char *filename, int num_queries, answer_type ans[], kNN_buffer *top_k[])
{
	int i;
	FILE *fp;
	char name[MAX_LEN] = "";
	char temp[MAX_LEN];
	
	strcat(name, filename);
	strcpy(temp, name);
	char *temp2 = strtok(name, "." );
	i = 1;
	while((fp = fopen(temp, "r")) != NULL){   // resultファイル名を上書きしないようにする
		fclose(fp);
		char str[MAX_LEN] = {0};
		memset(temp, '\0', MAX_LEN);
		strcpy(temp, temp2);
		sprintf(str, "%d" , i);
		strcat(temp, "_");
		strcat(temp, str);
		strcat(temp, ".csv");
		i++;
	}

	fp = fopen(temp, "w");
	if(fp == NULL) {
		printf("ファイルが開けません（filename = %s, temp = %s）\n", filename, temp);
		return;
	}
	fprintf(fp, "nearest_idx, nearest_dist, ans_idx[0], ans_dist[0], ans_idx[1], ans_dist[1], ... \n");
	for(int i = 0; i < num_queries; i++) {
		fprintf(fp, "%d,%d,", ans[i].data_num, ans[i].dist);
		int k;
		for(k = 0; k < top_k[i]->k - 1; k++) {
			fprintf(fp, "%d, %d,", top_k[i]->buff[k].data_num, top_k[i]->buff[k].dist);
		}
		fprintf(fp, "%d, %d\n", top_k[i]->buff[k].data_num, top_k[i]->buff[k].dist);
	}
	fclose(fp);
}

// #ifdef NUM_NN
void out_result_NN(char *filename, int num_queries, answer_type_NN ans[], kNN_buffer *top_k[])
{
	int i;
	FILE *fp;
	char name[MAX_LEN] = "";
	char temp[MAX_LEN];
	
	strcat(name, filename);
	strcpy(temp, name);
	char *temp2 = strtok(name, "." );
	i = 1;
	while((fp = fopen(temp, "r")) != NULL){   // resultファイル名を上書きしないようにする
		fclose(fp);
		char str[MAX_LEN] = {0};
		memset(temp, '\0', MAX_LEN);
		strcpy(temp, temp2);
		sprintf(str, "%d" , i);
		strcat(temp, "_");
		strcat(temp, str);
		strcat(temp, ".csv");
		i++;
	}

	fp = fopen(temp, "w");
	if(fp == NULL) {
		printf("ファイルが開けません（filename = %s, temp = %s）\n", filename, temp);
		return;
	}
	#ifdef SELF_EVAL
	fprintf(fp, "nearest_idx, nearest_dist, ans_idx[0], ans_dist[0], ans_idx[1], ans_dist[1], ... \n");
	for(int i = 0; i < num_queries; i++) {
		#ifndef ANSWER_DIST_FLOAT
		fprintf(fp, "%d,%d,", ans[i].data_num[0], ans[i].dist[0]);
		#else
		fprintf(fp, "%d,%.8f,", ans[i].data_num[0], sqrt(ans[i].dist[0] * 2));
		#endif
		int k;
		for(k = 0; k < top_k[i]->k - 1; k++) {
			#ifndef ANSWER_DIST_FLOAT
			fprintf(fp, "%d, %d,", top_k[i]->buff[k].data_num, top_k[i]->buff[k].dist);
			#else
			fprintf(fp, "%d, %.8f,", top_k[i]->buff[k].data_num, sqrt(top_k[i]->buff[k].dist) / 500.0);
			#endif
		}
		#ifndef ANSWER_DIST_FLOAT
		fprintf(fp, "%d, %d\n", top_k[i]->buff[k].data_num, top_k[i]->buff[k].dist);
		#else
		fprintf(fp, "%d, %.8f\n", top_k[i]->buff[k].data_num, sqrt(top_k[i]->buff[k].dist) / 500.0);
		#endif
	}
	#else
	#ifdef KNN_WITH_DISTANCE
	fprintf(fp, "nearest_idx, nearest_dist, ans_idx[0], ans_dist[0], ans_idx[1], ans_dist[1], ... \n");
	for(int i = 0; i < num_queries; i++) {
		fprintf(fp, "secret, secret, ");
		int k;
		for(k = 0; k < top_k[i]->k - 1; k++) {
			#ifndef ANSWER_DIST_FLOAT
			fprintf(fp, "%d, %d,", top_k[i]->buff[k].data_num + 1, top_k[i]->buff[k].dist);
			#else
			fprintf(fp, "%d, %.8f,", top_k[i]->buff[k].data_num + 1, sqrt(top_k[i]->buff[k].dist) / 500.0);
			#endif
		}
		#ifndef ANSWER_DIST_FLOAT
		fprintf(fp, "%d, %d\n", top_k[i]->buff[k].data_num + 1, top_k[i]->buff[k].dist);
		#else
		fprintf(fp, "%d, %.8f\n", top_k[i]->buff[k].data_num + 1, sqrt(top_k[i]->buff[k].dist) / 500.0);
		#endif
	}
	#else
	// Output knns only
	for(int i = 0; i < num_queries; i++) {
		int k;
		for(k = 0; k < top_k[i]->k - 1; k++) {
			fprintf(fp, "%d, ", top_k[i]->buff[k].data_num + 1);
		}
		fprintf(fp, "%d\n", top_k[i]->buff[k].data_num + 1);
	}
	#endif
	#endif
	fclose(fp);
}
// #endif

void out_result2(char *filename, int w, int d, double etime, double num_c, int num_queries, answer_type ans[], kNN_buffer *top_k[], int sorted)
{
	int i;
	FILE *fp;
	char name[MAX_LEN] = "";
	char temp[MAX_LEN];
	double sum = 0;
	strcat(name, filename);
	strcpy(temp, name);
	char *temp2 = strtok(name, "." );
	i = 1;
	while((fp = fopen(temp, "r")) != NULL){   // resultファイル名を上書きしないようにする
		fclose(fp);
		char str[MAX_LEN] = {0};
		memset(temp, '\0', MAX_LEN);
		strcpy(temp, temp2);
		sprintf(str, "%d" , i);
		strcat(temp, "_");
		strcat(temp, str);
		strcat(temp, ".csv");
		i++;
	}

	fp = fopen(temp, "w");
	if(fp == NULL) {
		printf("ファイルが開けません（filename = %s, temp = %s）\n", filename, temp);
		return;
	}

	for(int i = 0; i < num_queries; i++) {
		if(ans[i].dist == top_k[i]->buff[0].dist) {
			sum++;
		}
	}

	#ifdef NUM_THREADS
	int nt = NUM_THREADS;
	#else
	int nt = 1;
	#endif
	
	char *m = "";
	#ifdef SEQUENTIAL_FILTERING
	m = "SF";
	#endif
	#ifdef SEQUENTIAL_FILTERING_USING_BUCKET
	m = "BKT";
	#endif
	#ifdef SEQUENTIAL_FILTERING_USING_HAMMING
	m = "HAMMING";
	#endif
	#ifdef FILTERING_BY_SKETCH_ENUMERATION_C2N
	m = "ENU";
	#endif
	#ifdef WITHOUT_FILTERING
	m = "NO";
	#endif
	

	fprintf(fp, "w, dim, nt, time, K, recall, method, sort\n");
	fprintf(fp, "%d, %d, %d, %.9lf, %f, %f, %s, %d, mark\n", w, d, nt, etime, num_c, sum / num_queries * 100, m, sorted);

	fclose(fp);
}

// #ifdef DOUBLE_FILTERING

double out_result_double(char *filename, int n_w, int e_w, int d, double e_time_1st, double e_time_score, double e_time_2nd, double e_time_kNN,  double etime, double nc_1st, double nc_2nd, 
                         int num_queries, answer_type ans[], kNN_buffer *top_k[], char *query_file)
{
	int i;
	FILE *fp;
	char name[MAX_LEN] = "";
	char temp[MAX_LEN];
	double sum = 0;
	strcat(name, filename);
	strcpy(temp, name);
	char *temp2 = strtok(name, "." );
	i = 1;
	while((fp = fopen(temp, "r")) != NULL){   // resultファイル名を上書きしないようにする
		fclose(fp);
		char str[MAX_LEN] = {0};
		memset(temp, '\0', MAX_LEN);
		strcpy(temp, temp2);
		sprintf(str, "%d" , i);
		strcat(temp, "_");
		strcat(temp, str);
		strcat(temp, ".csv");
		i++;
	}

	fp = fopen(temp, "w");
	if(fp == NULL) {
		printf("ファイルが開けません（filename = %s, temp = %s）\n", filename, temp);
		return 0;
	}

	for(int i = 0; i < num_queries; i++) {
		if(ans[i].dist == top_k[i]->buff[0].dist) {
			sum++;
		}
	}

	#ifdef NUM_THREADS
	int nt = NUM_THREADS;
	#else
	int nt = 1;
	#endif
	
	#ifdef SEQUENTIAL_FILTERING
	char *m = "SF";
	#endif
	#ifdef SEQUENTIAL_FILTERING_USING_BUCKET
	char *m = "BKT";
	#endif
	#ifdef SEQUENTIAL_FILTERING_USING_HAMMING
	char *m = "HAMMING";
	#endif
	#ifdef FILTERING_BY_SKETCH_ENUMERATION
	char *m = "ENU";
	#endif
	#ifdef FILTERING_BY_SKETCH_ENUMERATION_HAMMING
	char *m = "ENU_HAM";
	#endif
	#ifdef FILTERING_BY_SKETCH_ENUMERATION_C2N
	char *m = "ENU_C2N";
	#endif
	#ifdef FILTERING_BY_SKETCH_ENUMERATION_INF
	char *m = "ENU_INF";
	#endif
	#ifdef WITHOUT_FILTERING
	char *m = "NO";
	#endif
	#ifdef DOUBLE_FILTERING
	char *m = "DF";
	#endif
	
	#ifndef SCORE_P
	#if defined(SCORE_1)
	#define SCORE_P 10
	#elif defined(SCORE_2)
	#define SCORE_P 20
	#else
	#define SCORE_P 10
	#endif
	#endif

	#ifndef SCORE_P_1ST
	#define SCORE_P_1ST SCORE_P
	#endif
	#ifndef SCORE_P_2ND
	#define SCORE_P_2ND 0
	#endif

	if(SCORE_P_1ST == DBL_MAX) {
		fprintf(fp, "n_w, e_w, dim, nt, 1st, score, 2nd, kNN, time, p_1,  k_1'/n(ppm), p_2, k_2'/n(ppm), recall, method, query\n");
		fprintf(fp, "%d, %d, %d, %d, %.9lf, %.9lf, %.9lf, %.9lf, %.9lf, DBL_MAX, %.2lf, %lf, %.2lf, %lf, %s, %s, mark\n", 
		n_w, e_w, d, nt, e_time_1st, e_time_score, e_time_2nd, e_time_kNN, etime, nc_1st, SCORE_P_2ND / 10.0, nc_2nd, sum / num_queries * 100, m, query_file);
		printf(     "n_w, e_w, dim, nt, 1st, score, 2nd, kNN, time, p_1,  k_1'/n(ppm), p_2, k_2'/n(ppm), recall, method, query\n");
		printf(     "%d, %d, %d, %d, %.9lf, %.9lf, %.9lf, %.9lf, %.9lf, DBL_MAX, %.2lf, %lf, %.2lf, %lf, %s, %s, mark\n", 
		n_w, e_w, d, nt, e_time_1st, e_time_score, e_time_2nd, e_time_kNN, etime, nc_1st, SCORE_P_2ND / 10.0, nc_2nd, sum / num_queries * 100, m, query_file);
	} else {
		fprintf(fp, "n_w, e_w, dim, nt, 1st, score, 2nd, kNN, time, p_1,  k_1'/n(ppm), p_2, k_2'/n(ppm), recall, method, query\n");
		fprintf(fp, "%d, %d, %d, %d, %.9lf, %.9lf, %.9lf, %.9lf, %.9lf, %lf, %.2lf, %lf, %.2lf, %lf, %s, %s, mark\n", 
		n_w, e_w, d, nt, e_time_1st, e_time_score, e_time_2nd, e_time_kNN, etime, SCORE_P_1ST / 10.0, nc_1st, SCORE_P_2ND / 10.0, nc_2nd, sum / num_queries * 100, m, query_file);
		printf(     "n_w, e_w, dim, nt, 1st, score, 2nd, kNN, time, p_1,  k_1'/n(ppm), p_2, k_2'/n(ppm), recall, method, query\n");
		printf(     "%d, %d, %d, %d, %.9lf, %.9lf, %.9lf, %.9lf, %.9lf, %lf, %.2lf, %lf, %.2lf, %lf, %s, %s, mark\n", 
		n_w, e_w, d, nt, e_time_1st, e_time_score, e_time_2nd, e_time_kNN, etime, SCORE_P_1ST / 10.0, nc_1st, SCORE_P_2ND / 10.0, nc_2nd, sum / num_queries * 100, m, query_file);
	}

	fclose(fp);

	return sum / num_queries * 100;
}

// ここでは，NN 検索を行わないで，Recallを評価する．（double_filtering では，この関数を用いる前に，kNN_search で NN 検索を行ってる）
// 質問 q の正解情報: q の最近傍のデータ番号 ans[q].data_num と距離 ans[q].dist 
// フィルタリングで求めた，上位のものは top_k[q] に格納されている． k = top_k[q]->k
// top_k[q]->buff[i].data_num, top_k[q]->buff[i].dist (i = 0, ... , k - 1):
// 上位 i 番目のもののデータ番号と距離（ただし，フィルタリングの結果では，射影距離（優先順位））  
// Evaluate_filtering のときは，answer_check で，正解が見つかっているときには，正解情報から距離を top_k の先頭に与えている．
double recall_without_search(int num_queries, answer_type ans[], kNN_buffer *top_k[])
{
	int sum = 0;
	for(int q = 0; q < num_queries; q++) {
		if(ans[q].dist == top_k[q]->buff[0].dist) {
			sum++;
		}
	}
	return (double)sum / num_queries * 100;
}

// kNNの結果をデータ番号のみで評価する
double recall_kNN(int num_queries, answer_type_NN ans[], kNN_buffer *top_k[])
{
	int sum = 0;
	int k = top_k[0]->k;
	for(int q = 0; q < num_queries; q++) {
		for(int r = 0; r < k; r++) {
			for(int a = 0; a < k; a++) {
				if(top_k[q]->buff[r].data_num == ans[q].data_num[a]) {
//				if(top_k[q]->buff[r].dist <= ans[q].dist[k-1]) {
					sum++;
					break;
				}
			}
		}
	}
	return (double)sum / num_queries / k * 100;
}

// 1-NNの結果をdistで評価する
double recall_kNN_1(int num_queries, answer_type_NN ans[], kNN_buffer *top_k[])
{
	int sum = 0;
	for(int q = 0; q < num_queries; q++) {
		if(top_k[q]->buff[0].data_num == ans[q].data_num[0]) {
			sum++;
		}
	}
	return (double)sum / num_queries * 100;
}

#ifdef NUM_NN
double recall_without_search_NN(int num_queries, answer_type_NN ans[], kNN_buffer *top_k[])
{
//	fprintf(stderr, "recall_without_search_NN: NUM_NN = %d, k = %d\n", NUM_NN, top_k[0]->k); exit(0);
	int sum = 0;
	int k = top_k[0]->k;
	int sum_correct_NN = 0;
	int sum_correct_NN_k = 0;
	for(int q = 0; q < num_queries; q++) {
		int correct_NN = (ans[q].dist[0] == top_k[q]->buff[0].dist);
		sum_correct_NN += correct_NN;
		int correct_NN_k = 0;
		for(int j = 0; j < k; j++) {
			if(ans[q].dist[k - 1] >= top_k[q]->buff[j].dist) {
				correct_NN_k++;
			}
		}
		sum_correct_NN_k += correct_NN_k;
/*
		if(q < 30) {
			fprintf(stderr, "ans.dist = ");
			for(int j = 0; j < k; j++) {
				fprintf(stderr, "%6d", ans[q].dist[j]);
			}
			fprintf(stderr, "\n");
			fprintf(stderr, "top.dist = ");
			for(int j = 0; j < k; j++) {
				fprintf(stderr, "%6d", top_k[q]->buff[j].dist);
			}
			fprintf(stderr, "\n");
			fprintf(stderr, "%d, %2d\n", correct_NN, correct_NN_k);
//			getchar();
		}
*/

		for(int j = 0; j < k; j++) {
			if(ans[q].dist[k - 1] >= top_k[q]->buff[j].dist) {
				sum++;
			}
		}
	}
	fprintf(stderr, "sum_correct_NN = %d, sum_correct_NN_k = %d\n", sum_correct_NN, sum_correct_NN_k);
	return (double)sum / num_queries / k * 100;
}
#endif

// #endif


