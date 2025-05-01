/*
v2: データ番号を1からではなく0からに変更．付随情報には元の1開始の番号．
For PUBMED23 (23M * 384)
.h5ファイルから検索対象データ（train）に対する2種類の質問（itest, otest）の正解情報（1000近傍）
を読み込んで，make_answerで作成したものと同様の csv で出力する．

提供される h5 ファイルには，検索対象のデータセットと２種類の質問と正解情報が含まれている．
$ h5ls -r h5/benchmark-dev-pubmed23.h5
/                        Group
/itest                   Group
/itest/dists             Dataset {11000, 1000}
/itest/knns              Dataset {11000, 1000}
/itest/queries           Dataset {11000, 384}
/otest                   Group
/otest/dists             Dataset {11000, 1000}
/otest/knns              Dataset {11000, 1000}
/otest/queries           Dataset {11000, 384}
/train                   Dataset {23887701, 384}

コンパイル：h5cc -O3 ground_trueth_to_csv_PUBMED_v2.c -o ground_trueth_to_csv_PUBMED_v2
実行時パラメータ: file.h5 [itest | otest] num_queries k > answer_[i | o].csv
   1. file.h5       = input file of h5 format
   2. itest | otest = 
   3. num_queries
   4. k             = number of NN
*/

#include <stdio.h>
#include <hdf5.h>
#include <stdlib.h>
#include <string.h>
#include <limits.h>
#include <math.h>



void **read_dataset(char *h5file, char *dataset_name, int *num)
{
	// ファイルIDの取得
	hid_t file_id = H5Fopen(h5file, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(file_id == H5I_INVALID_HID) {
		fprintf(stderr, "cannot open file: %s\n", h5file);
		return NULL;
	}
	fprintf(stderr, "H5Fopen %s OK\n", h5file);

	// データセットIDの取得
	hid_t dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
	if(dataset_id == H5I_INVALID_HID) {
		fprintf(stderr, "cannot open dataset: %s\n", dataset_name);
		return NULL;
	}
	fprintf(stderr, "H5Dopen %s OK\n", dataset_name);

	// データ空間IDの取得
	hid_t filespace_id = H5Dget_space(dataset_id);
	if(filespace_id == H5I_INVALID_HID) {
		fprintf(stderr, "cannot get space\n");
		return NULL;
	}

	// データの次元数の取得（ここでは，2次元の配列を想定）
	int ndims = H5Sget_simple_extent_ndims(filespace_id);
	if (ndims != 2) {
		fprintf(stderr, "Error: Dimension of %s must be 2.\n", dataset_name);
		return NULL;
	}
	hsize_t *dims = malloc(sizeof(hsize_t) * ndims);
	H5Sget_simple_extent_dims(filespace_id, dims, NULL);
	fprintf(stderr, "dims[0] = %lld, dims[1] = %lld", dims[0], dims[1]);

	// 配列の行数（各行はベクトルとみなし，ベクトル数を取得，dims[1] は近傍数（おそらく1000）
	hssize_t npoints = H5Sget_simple_extent_npoints(filespace_id);
	int num_rows = npoints / dims[1];
	fprintf(stderr, "npoints = %lld, num_lines = %d\n", npoints, num_rows);
	if(*num == 0 || *num > num_rows) {
		*num = num_rows;
	}
	
	// データ型の取得
	hid_t type = H5Dget_type(dataset_id);
	if (!H5Tequal(type, H5T_NATIVE_FLOAT) && !H5Tequal(type, H5T_NATIVE_INT)) {
		fprintf(stderr, "type is not H5T_NATIVE_FLOAT or H5T_NATIVE_INT\n");
		return NULL;
	}
	// バッファの型の作成
	hsize_t dimsmr[1] = {dims[1]};
	hid_t memspace_id = H5Screate_simple(1, dimsmr, NULL);

	// 配列を用意
	void **data;
	data = malloc(*num * sizeof(void *));
	for(int i = 0; i < *num; i++) {
		data[i] = malloc(sizeof(int) * dims[1]);
	}
	fprintf(stderr, "malloc for data (total size = %lld) OK\n", sizeof(int) * dims[1] * *num);
	// バッファの読み取り位置
	int start = 0;
	hsize_t roffset[2] = {start, 0};

	// 読み取り数
	hsize_t rcount[2] = {1, dims[1]};

	herr_t status;
	
	for(int i = 0; i < *num; i++) {
		// 読み取り位置の設定
		status = H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, roffset, NULL, rcount, NULL);
		// 読み取り
		status = H5Dread(dataset_id, type, memspace_id, filespace_id, H5P_DEFAULT, data[i]);
		// 読み取り位置の更新
		roffset[0]++;
	}

	// 後始末
	status = H5Sclose(filespace_id);
	status = H5Sclose(memspace_id);
	status = H5Dclose(dataset_id);
	status = H5Fclose(file_id);
	free(dims);

	return data;
}

int main(int argc, char *argv[])
{
	if(argc != 5) {
		fprintf(stderr, "usage> %s h5_file [ itest | otest ] num_queries k > answer_file.csv\n", argv[0]);
		return 1;
	}

	char *h5file = argv[1];
	char *query = argv[2];
	int num_queries = atoi(argv[3]);
	int last_k = atoi(argv[4]);
	if(last_k > 1000) last_k = 1000;
	char *ds_dists, *ds_knns;
	int offset = 0;
	
	if(strcmp(query, "itest") == 0) {
		ds_dists = "itest/dists";
		ds_knns  = "itest/knns";
		offset   = 30000000;
	} else if(strcmp(query, "otest") == 0) {
		ds_dists = "otest/dists";
		ds_knns  = "otest/knns";
		offset   = 40000000;
	}

	// 距離のリストを読み込む

	// kNNの距離を読み込む．
	float **dist_list = (float **)read_dataset(h5file, ds_dists, &num_queries);
	if(dist_list == NULL) {
		return -1;
	} else {
		fprintf(stderr, "read_dataset %s OK. num_data = %d\n", ds_dists, num_queries);
	}

	// kNNの解（データ番号）を読み込む．
	int **knn_list = (int **)read_dataset(h5file, ds_knns, &num_queries);
	if(knn_list == NULL) {
		return -1;
	} else {
		fprintf(stderr, "read_dataset %s OK. num_data = %d\n", ds_knns, num_queries);
	}

	printf("answer for %s\n", query);
	printf("q_number, q_id");
	for(int k = 1; k <= last_k; k++) {
		printf(", NN[%d]_idx, NN[%d]_dist, NN[%d]_id", k, k, k);
	}
	printf("\n");
	for(int q = 0; q < num_queries; q++) {
		printf("%d,%d", q, q + offset);
		for(int k = 0; k < last_k; k++) {
			printf(", %d, %lf, %d", knn_list[q][k] - 1, dist_list[q][k], knn_list[q][k]);
		}
		printf("\n");
	}
	
	return 0;
}
