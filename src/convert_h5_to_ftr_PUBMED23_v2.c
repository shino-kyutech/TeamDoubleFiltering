//  v2 is used batch read to speed-up
//  float -> char only (float16 is not supported)
//  multi-threading is not used
//  .h5 -> .ftr
//  support on native float or 16-bit float data
//  modified for PUBMED23 used in SISAP Indexing Challenge 2025 
#include <stdlib.h>
#include <hdf5.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <limits.h>
#include <math.h>
#include <omp.h>
// PUBMED23 のデータを.h5（float）から.ftr（signed char）に変換するプログラム
// 元のデータにはデータ番号（ID）が付いていないので，読み込んだデータ順番に
// パラメタで指定するオフセットを加えたものをIDとする．
// 今回のチャレンジでは，オフセットはデータセットについては，1とする．
// 質問のIDは，itestにはオフセット300,000,001, otestには400,000,001とする．
// 提供される h5 ファイルには，検索対象のデータセットと２種類の質問と正解情報が含まれている．
// $ h5ls -r h5/benchmark-dev-pubmed23.h5
// /                        Group
// /itest                   Group
// /itest/dists             Dataset {11000, 1000}
// /itest/knns              Dataset {11000, 1000}
// /itest/queries           Dataset {11000, 384}
// /otest                   Group
// /otest/dists             Dataset {11000, 1000}
// /otest/knns              Dataset {11000, 1000}
// /otest/queries           Dataset {11000, 384}
// /train                   Dataset {23887701, 384}
//
// 後の取り扱いのために，1M（100万）個のデータ単位に分割する．若干端数が出るが捨てる（？）
// 使用法：一度に読み込むベクトル数を変更するときはマクロBATCH_SIZEを変更してコンパイルすること．
// コンパイル：h5cc -O3 -DBATCH_SIZE=500 convert_h5_to_ftr_PUBMED23.c -o convert_h5_to_ftr_PUBMED23
//   実行時パラメータ: file.h5 dataset_name start count offset convert_factor file.ftr
//   1. file.h5      = input file of h5 format
//   2. dataset_name = 今回は，データセットはtrain, 質問は，itest/queriesまたはotest/queries
//   start        = data number of the 1st input data
//   count        = number of data to convert (0 for all to end)
//   offset       = データ番号のオフセット．データセットでは1, 質問は300,000,001または400,000,001
//   file.ftr     = output file ftr format
// 

typedef signed char ftr_element_type; // 特徴データの要素の型
typedef ftr_element_type *ftr_type; // 特徴データベクトルの型
#define FTR_DIM 384       // 特徴データの次元数
#define DATA_NUM 1000000  // 出力するデータ総数（実際には，後で書換える）
#ifndef BATCH_SIZE
#define BATCH_SIZE 1000   // 一度にh5から読み込むベクトル数，まとめてftrに書き出すベクトル数
#endif

int CONVERT_FACTOR = 255; // 整数化する前にこの数を掛ける（実際には，実行時のパラメータで指定する） 

#define CONVERT(v) ((int)((v) * CONVERT_FACTOR)) // float -> int の変換

// auxi_type（特徴データの付随情報）の型
typedef unsigned int auxi_type; // PUBMED23では 32-bit 整数で十分（データID）

// 特徴データファイルのヘッダ情報
typedef struct {
	int element_size;	// 各次元のサイズ（おそらく1(char)，または2(short)）
	int data_dim;   	// 特徴次元数（image:64, music:96, colors:112, DeCAF:4096, Deep1B）
	int data_num;		// データ数（縮小するときは，後で変更する）
	int auxi_size;		// auxi_type（付随情報）のサイズ（PUBMED23では，データ数が約2,300万，32-bit intに収まるので 4）
} ftr_header_type;

// 特徴データと付随情報（データID）をまとめた構造体
typedef struct {
	ftr_element_type ftr[FTR_DIM];
	auxi_type data_id;
} struct_ftr_id;

ftr_header_type db_header = {sizeof(ftr_element_type), FTR_DIM, DATA_NUM, sizeof(auxi_type)};

// floatのvectorをsigned charのftrとdata_idに変換する 	
void convert_vectors(struct_ftr_id *ftr_id_buff, float *buff, auxi_type data_id, long *overflow)
{
  for(int i = 0; i < FTR_DIM; i++) {
	  int val = CONVERT(buff[i]);
		if(val > 127) {
			val = 127;
			(*overflow)++;
		} else if(val < -127) {
			val = -127;
			(*overflow)++;
		}
		ftr_id_buff->ftr[i] = val;
  }
	ftr_id_buff->data_id = data_id;
}

int main(int argc, const char **argv)
{
	if (argc < 8) {
		fprintf(stderr, "Usage: %s %s\n", argv[0], "file.h5 dataset_name start count offset convert_factor file.ftr");
		for(int i = 0; i < argc; i++) {
			fprintf(stderr, "argv[%d] = %s\n", i, argv[i]);
		}
		return -1;
	}

	const char *h5file = argv[1];
	const char *dataset_name = argv[2];
	int start = atoi(argv[3]);
	int count = atoi(argv[4]);
	int offset = atoi(argv[5]);
	CONVERT_FACTOR = atoi(argv[6]);
	const char *ftr_file = argv[7];

	// ファイルIDの取得
	hid_t file_id = H5Fopen(h5file, H5F_ACC_RDONLY, H5P_DEFAULT);
	if(file_id == H5I_INVALID_HID) {
		fprintf(stderr, "cannot open file: %s\n", h5file);
		return -1;
	}

	// データセットIDの取得
	hid_t dataset_id = H5Dopen(file_id, dataset_name, H5P_DEFAULT);
	if(dataset_id == H5I_INVALID_HID) {
		fprintf(stderr, "cannot open dataset: %s\n", dataset_name);
		return -1;
	}

	// データ空間IDの取得
	hid_t filespace_id = H5Dget_space(dataset_id);
	if(filespace_id == H5I_INVALID_HID) {
		fprintf(stderr, "cannot get space\n");
		return -1;
	}

	// データの次元数の取得（ここでは，2次元の配列を想定）
	int ndims = H5Sget_simple_extent_ndims(filespace_id);
	if (ndims != 2) {
		fprintf(stderr, "Error: Dimension of %s must be 2.\n", dataset_name);
		return -1;
	}
	hsize_t *dims = malloc(sizeof(hsize_t) * ndims);
	H5Sget_simple_extent_dims(filespace_id, dims, NULL);
	if(dims[1] != FTR_DIM) {
		fprintf(stderr, "data dimension error: expected = %d, stored = %lld", FTR_DIM, dims[1]);
		return -1;
	}

	// 配列の行数（各行はベクトルとみなし，ベクトル数を取得，dims[1] はベクトルの次元数）
	hssize_t npoints = H5Sget_simple_extent_npoints(filespace_id);
	int nvectors = npoints / dims[1];
	fprintf(stderr, "npoints = %lld, nvectors = %d\n", npoints, nvectors);

	// データ型の取得
	hid_t type = H5Dget_type(dataset_id);

	// バッファの型の作成
	hid_t memspace_id = H5Screate_simple(2, (hsize_t[]){BATCH_SIZE, FTR_DIM}, NULL);

	float (*rdata)[FTR_DIM];
	if(H5Tequal(type, H5T_NATIVE_FLOAT)){
		rdata = malloc(dims[1] * sizeof(float) * BATCH_SIZE);
		fprintf(stderr, "Data type = H5NATIVE_FLOAT (32-bit float)\n");
	} else {
		fprintf(stderr, "Data type != H5NATIVE_FLOAT\n");
		return -1;
	}

	// バッファの読み取り位置
	hsize_t roffset[2] = {start, 0};

	// 読み取り数
	hsize_t rcount[2] = {BATCH_SIZE, FTR_DIM};

	if(count == 0 || count > nvectors) {
		count = nvectors;
	}

	FILE *fp;
	struct_ftr_id *ftr_id_buff = malloc(sizeof(struct_ftr_id) * BATCH_SIZE);

	long overflow = 0;

	if((fp = fopen(ftr_file, "wb"))  == NULL) {
		fprintf(stderr, "cannot open ftr_file: file name = %s\n", ftr_file);
		exit(1);
	}
	fprintf(stderr, "open ftr_file OK: file name = %s\n", ftr_file);
	db_header.data_num = count;
	fwrite(&db_header, 1, sizeof(ftr_header_type), fp);  // ヘッダ情報を書き出す

	int written = 0;
	for (hsize_t i = start; i < start + count && roffset[0] < nvectors; i += BATCH_SIZE) {
		roffset[0] = i;
		if(i + BATCH_SIZE > nvectors) {
			rcount[0] = nvectors - i;
		}
		if(rcount[0] != BATCH_SIZE) {
            H5Sclose(memspace_id);
            memspace_id = H5Screate_simple(2, (hsize_t[]){nvectors - i, FTR_DIM}, NULL);
		}

		
		// 読み取り位置の設定
		H5Sselect_hyperslab(filespace_id, H5S_SELECT_SET, roffset, NULL, rcount, NULL);

		// 読み取り
		herr_t status = H5Dread(dataset_id, H5T_NATIVE_FLOAT, memspace_id, filespace_id, H5P_DEFAULT, rdata);
		if (status < 0) {
			fprintf(stderr, "Error: failed to read batch at index %llu\n", (unsigned long long)i);
			return -1;
		}
		
		for(int k = 0; k < (int)rcount[0]; k++) {
			// floatの特徴データをcharに変換して struct_ftr_id に格納する
			convert_vectors(&ftr_id_buff[k], rdata[k], (int)(i + k + offset), &overflow);
		}
		
		// struct_ftr_id のレコードを書き込む
		if(fwrite(ftr_id_buff, sizeof(struct_ftr_id), rcount[0], fp) != rcount[0]) {
			fprintf(stderr, "fwrite error (ftr_id): data id = %d, record number = %d\n",
				(int)(i + offset), (int)(i));
			return -1;
		}
		written += rcount[0];

		// 読み取り位置の更新
		roffset[0] += rcount[0];
	}

	fseek(fp, 0L, SEEK_SET);
	db_header.data_num = written;
	fwrite(&db_header, 1, sizeof(ftr_header_type), fp);  // ヘッダ情報を書き出す

	fprintf(stderr, "CONVERT_FACTOR = %d, num_vectors = %d, dim = %d, overflow = %ld\n", CONVERT_FACTOR, count, (int)dims[1], overflow);
	fclose(fp);

	// 後始末
	H5Sclose(filespace_id);
	H5Sclose(memspace_id);
	H5Dclose(dataset_id);
	H5Fclose(file_id);
	free(ftr_id_buff);
	free(rdata);
	free(dims);

	return 0;
}
