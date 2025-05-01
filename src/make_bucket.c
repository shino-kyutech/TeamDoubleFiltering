// 特徴データのスケッチを作成してファイルに保存する．検索はしない．
// 引数1: ftrファイル1
// 引数2: ftrファイル2
// ...
// 以下のものはマクロで定義して与える
// ピボットファイル（csv形式），スケッチ幅（w），スケッチファイル（バイナリファイル）

#include "parm.h"
#include "config.h"
#include "bit_op.h"
#include "ftr.h"
//#include "kNN_search.h"
#include "sketch.h"

int main(int argc, char *argv[])
{
	int num_ftr_files = argc - 1;
	char **dataset_ftr_filename = argv + 1;
	dataset_handle dh;

	#if defined(_OPENMP) && defined(NUM_THREADS)
	dh.num_threads = NUM_THREADS;
	#else
	dh.num_threads = 1;
	#endif

	#ifdef FTR_SORT_BY_SKETCH_IN_ADVANCE
	fprintf(stderr, "FTR is sorted in advanve\n");
	dh.sorted = 1;
	#else
	fprintf(stderr, "FTR is NOT sorted. Arrangement is as is.\n");
	dh.sorted = 0;
	#endif

	#if defined(FTR_ON_SECONDARY_MEMORY)
	dh.ftr_on = SECONDARY_MEMORY;
	dh.mf = (struct_multi_ftr **)malloc(sizeof(struct_multi_ftr *) * dh.num_threads);
	for(int t = 0; t < dh.num_threads; t++) {
		dh.mf[t] = open_multi_ftr(num_ftr_files, dataset_ftr_filename, BLOCK_SIZE);
	}
	dh.ds = NULL;
	int num_data = dh.mf[0]->num_data;
	fprintf(stderr, "total number of data = %d\n", num_data);
	#else
	#error "FTR_ON_SECONDARY_MEMORY should be defined."
	#endif // FTR_ON

	#ifdef _OPENMP
	#if defined(NUM_THREADS) && NUM_THREADS > 1
	int nt = NUM_THREADS;
	omp_set_num_threads(NUM_THREADS);
	#endif
	#else
	int nt = 1;
	#endif

	char *pivot_file = PIVOT_FILE;
	char *bucket_file = BUCKET_FILE;

	fprintf(stderr, "PJT_DIM=%d\n", PJT_DIM);

	#if defined(PARTITION_TYPE_QBP)
	pivot_type *pivot = new_pivot(QBP);
	#elif defined(PARTITION_TYPE_PQBP)
	pivot_type *pivot = new_pivot(PQBP);
	#elif defined(PARTITION_TYPE_SQBP)
	pivot_type *pivot = new_pivot(SQBP);
	#elif defined(PARTITION_TYPE_CSQBP)
	pivot_type *pivot = new_pivot(SQBP); // 当面 CSQBP は SQBP で対応
	#endif

	#ifndef PRE_ROTATION
	read_pivot(pivot_file, pivot);
//	for(int j = 0; j < PJT_DIM; j++) {
//		printf("pivot->p[%d] = ", j);
//		for(int k = 0; k < pivot->num_used[j]; k++) {
//			printf(", [%d] = %d", pivot->used[j][k], pivot->p[j][pivot->used[j][k]]);
//		}
//		printf("\n");
//	}
	fprintf(stderr, "read pivot OK\n");
	#else
    int rotation[FTR_DIM];
	read_pivot_with_rotation(pivot_file, pivot, rotation);
	fprintf(stderr, "read pivot with rotation OK\n");
	#endif

	fprintf(stderr, "make sketch of data ... \n");
	sketch_type *sketch = (sketch_type *)malloc(sizeof(sketch_type) * num_data);

	int *data_num = (int *)malloc(sizeof(int) * num_data);
	#if defined(_OPENMP) && NUM_THREADS > 1
	#pragma omp parallel for
	for(int i = 0; i < num_data; i++) data_num[i] = i;
	#endif
	
	#if defined(_OPENMP) && NUM_THREADS > 1
	for(int t = 0; t < nt; t++) {
		if(dh.ftr_on == SECONDARY_MEMORY) dh.mf[t]->read_in = 0;
	}
	#pragma omp parallel for
	#else
	if(dh.ftr_on == SECONDARY_MEMORY) dh.mf[0]->read_in = 0;
	#endif
	for(int i = 0; i < num_data; i++) {
		struct_multi_ftr *mf = NULL;
		#if defined(_OPENMP) && NUM_THREADS > 1
		int t = omp_get_thread_num();
		if(dh.ftr_on == SECONDARY_MEMORY) mf = dh.mf[t];
		#else
		if(dh.ftr_on == SECONDARY_MEMORY) mf = dh.mf[0];
		#endif

		struct_ftr_id *ftr_id_p = get_next_ftr_id_from_multi_ftr(mf, data_num, i, num_data);
		#ifdef PRE_ROTATION
        pre_rotation_ftr(ftr_id_p->ftr, rotation);
		#endif

		#ifndef EXPANDED_SKETCH
		sketch[i] = data_to_sketch(ftr_id_p->ftr, pivot);
//		printf("i = %d, sk = ", i); print_bin(sketch[i]); printf("\n"); getchar();
		#else
		data_to_sketch(ftr_id_p->ftr, pivot, sketch[i]);
//		printf("sketch[%d] = ", i);
//		print_bin_expanded_2(sketch[i], SKETCH_SIZE);
//		printf("\n"); getchar();
		#endif
	}
	fprintf(stderr, "make sketch DONE.\n");

/*
	fprintf(stderr, "make bucket ... ");
	struct_bucket *bucket = new_bucket(num_data, sketch);
	fprintf(stderr, "OK\n");

	fprintf(stderr, "write bucket (file = %s) ... ", bucket_file);
	write_bucket(bucket_file, bucket);
	fprintf(stderr, "OK\n");

	free_bucket(bucket);
*/
	int *idx = (int *)malloc(sizeof(int) * num_data);
	for(int i = 0; i < num_data; i++) {
		idx[i] = i;
	}
	fprintf(stderr, "sort sketch starts: num_data = %d\n", num_data);
	quick_sort_for_sketch(idx, sketch, 0, num_data - 1);
	fprintf(stderr, "sort sketch ended: num_data = %d\n", num_data);

	int num_nonempty_buckets = 0;
	for(int i = 0; i < num_data; ) {
		int j;
		for(j = i + 1; j < num_data && comp_sketch(sketch[idx[i]], sketch[idx[j]]) == 0; j++);
		num_nonempty_buckets++;
		i = j;
	}

	fprintf(stderr, "write bucket (file = %s) ... ", bucket_file);
	FILE *fp;
	if((fp = fopen(bucket_file, "wb"))  == NULL) {
		fprintf(stderr, "Write open bucket file error, file name = %s\n", bucket_file);
		exit(0);
	}
	if(fwrite(&num_data, sizeof(int), 1, fp) != 1) {  // num_data を書き出す
		fprintf(stderr, "fwrite error (num_data) file = %s\n", bucket_file);
		exit(0);
	}
	if(fwrite(idx, sizeof(int) * num_data, 1, fp) != 1) {  // idx[data_num] を書き出す
		fprintf(stderr, "fwrite error (idx) file = %s\n", bucket_file);
		exit(0);
	}
	if(fwrite(&num_nonempty_buckets, sizeof(int), 1, fp) != 1) {  // num_nonempty_buckets を書き出す
		fprintf(stderr, "fwrite error (num_nonempty_buckets) file = %s\n", bucket_file);
		exit(0);
	}
	fprintf(stderr, "num_nonempty_buckets = %d, average number of elements in nonempty buckets = %.2lf\n", num_nonempty_buckets, (double)num_data / num_nonempty_buckets);
	for(int i = 0; i < num_data; ) {
		#ifndef EXPANDED_SKETCH
		if(fwrite(&sketch[idx[i]], sizeof(sketch_type), 1, fp) != 1) {  // sketch s を書き出す
			fprintf(stderr, "fwrite error (sketch) file = %s\n", bucket_file);
			exit(0);
		}
		#else
		if(fwrite(sketch[idx[i]], sizeof(sketch_type), 1, fp) != 1) {  // sketch s を書き出す
			fprintf(stderr, "fwrite error (sketch) file = %s\n", bucket_file);
			exit(0);
		}
		#endif
		int j;
		for(j = i + 1; j < num_data && comp_sketch(sketch[idx[i]], sketch[idx[j]]) == 0; j++);
		int num_elements = j - i;
		if(fwrite(&num_elements, sizeof(int), 1, fp) != 1) {  // num_elements を書き出す
			fprintf(stderr, "fwrite error (num_elements) file = %s\n", bucket_file);
			exit(0);
		}
		i = j;
	}
	fclose(fp);
	fprintf(stderr, "OK\n");

	return 0;
}
