//#pragma once
// 特徴データに関する基本操作（距離関数，データセットファイルの入出力）
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "ftr.h"
#include "e_time.h"

void print_header(FILE *f, char *filename, ftr_header_type *header)
{
	fprintf(f, "header of %s:\n", filename);
	fprintf(f, "size of feature element = %10d\n", header->element_size);
	fprintf(f, "data dimension          = %10d\n", header->data_dim);
	fprintf(f, "number of data          = %10d\n", header->num_data);
	fprintf(f, "size of auxiarily data  = %10d\n", header->auxi_size);
}

file_handle open_ftr_file(char *filename, ftr_header_type *header)
// 特徴ファイル（データベース用，質問データ用）
// ファイルを open して，ヘッダ情報を得る
{
	file_handle fh;

	if((fh = READ_OPEN(filename)) == OPEN_ERROR) {
		fprintf(stderr, "open_ftr_file error, file name = %s\n", filename);
		return OPEN_ERROR;
	}
	if(!READ(fh, header, sizeof(ftr_header_type))) {
		fprintf(stderr, "open_ftr_file error: file name = %s\n", filename);
		return OPEN_ERROR;
	}

	return fh;
}

int get_ftr_id(file_handle fh, ftr_header_type *header, int data_num, struct_ftr_id *ftr_id)
// 特徴データとデータIDを読む．retuen value (0 -> error, 1 -> no error) 
{
	if(data_num >= header->num_data) return 0;
	SEEK(fh, sizeof(ftr_header_type) + (long)data_num * sizeof(struct_ftr_id), SEEK_SET);
	return READ(fh, ftr_id, sizeof(struct_ftr_id)); 
}

// 特徴データとデータIDをnum_data個まとめて読む．retuen value (error -> 0, no error -> number of data read-in) 
// あまり大量のデータをまとめて読むとエラーになるみたい（4GBが限界？）
int get_ftr_id_bulk(file_handle fh, ftr_header_type *header, int data_num, int num_data, struct_ftr_id *ftr_id)
{
	// data_num = 0, num_data = 100 ならば，先頭から100個まとめて読み込む
	// このとき，header->num_data = 90 ならば，90個しか読み込めないので，num_data = 90 とする
	// data_num = 50 のときは，51番目のデータから読み込むが，data_num >= header->data_num であれば，EOFを通り過ぎている
	// まだ読み込むデータが残っているのであれば，残りの個数は，header->num_data - data_num である
//	printf("get_ftr_id_bulk: data_num = %d, num_data = %d\n", data_num, num_data);
	if(data_num >= header->num_data) return 0;
	if(data_num + num_data >= header->num_data) {
		num_data = header->num_data - data_num;
	}
	SEEK(fh, sizeof(ftr_header_type) + (long)data_num * sizeof(struct_ftr_id), SEEK_SET);
	return READ(fh, ftr_id, sizeof(struct_ftr_id) * num_data) ? num_data : 0; 
}

struct_dataset *read_dataset_n(int num_ftr_files, char *filename[])
{
	// ヘッダ情報の整合性を確認して，struct_dataset内に読み込んでおくようにした方がよいかも
	file_handle fh[num_ftr_files];
	ftr_header_type hd[num_ftr_files];
	int num_data = 0;

	for(int i = 0; i < num_ftr_files; i++) {
		if((fh[i] = open_ftr_file(filename[i], &hd[i])) == OPEN_ERROR) {
			fprintf(stderr, "cannot open ftr file = %s\n", filename[i]);
			exit(0);
		}
		num_data += hd[i].num_data;
	}
	fprintf(stderr, "open OK: num_data = %d, file_name = %s\n", num_data, filename[0]);
	struct_dataset *ds = MALLOC(sizeof(struct_dataset));
	ds->num_data = num_data;
	ds->sorted = 0;
	ds->ftr_id = MALLOC(sizeof(struct_ftr_id) * num_data);
	if(ds->ftr_id == NULL) {
		fprintf(stderr, "cannot allocate memory for %d data (probably too many)\n", num_data);
		exit(0);
	}
	int n = 0, r;
	for(int i = 0; i < num_ftr_files; i++) {
//		for(int j = 0; j < hd[i].num_data; j++) {
//			get_ftr_id(fh[i], &hd[i], j, &ds->ftr_id[n++]);
//		}
		int j;
		for(j = 0; j < hd[i].num_data; j += r) {
			if(!(r = get_ftr_id_bulk(fh[i], &hd[i], j, 1000, &ds->ftr_id[n + j]))) {
				fprintf(stderr, "READ error: file_name = %s, file_num = %d, num_read = %d, FTR_DIM = %d\n", filename[i], i, r, FTR_DIM);
				exit(0);
			}
		}
		if(j != hd[i].num_data) {
			fprintf(stderr, "READ error: file_name = %s, file_num = %d, num_read = %d, FTR_DIM = %d\n", filename[i], i, j, FTR_DIM);
			exit(0);
		}
		n += j;
//		if(!(r = get_ftr_id_bulk(fh[i], &hd[i], 0, hd[i].num_data, &ds->ftr_id[n]))) {
//			fprintf(stderr, "READ error: file_name = %s, file_num = %d, data_num = %d, FTR_DIM = %d\n", filename[i], i, n, FTR_DIM);
//			exit(0);
//		}
//		n += r;
		CLOSE(fh[i]);
	}
	return ds;
}

void free_dataset(struct_dataset *ds)
{
	free(ds->ftr_id);
	free(ds);
}

void sort_dataset(struct_dataset *ds, int idx[])
// データセット内の特徴データ（付随情報付き）を idx で指定した順に並べ直す．
{
	struct_ftr_id *temp = (struct_ftr_id *)MALLOC(sizeof(struct_ftr_id));
	char *done = (char *)calloc(ds->num_data, sizeof(char));

	if(ds->sorted != 0) {
		fprintf(stderr, "Dataset is already sorted.\n");
		return;
	}
	for(int i = 0; i < ds->num_data; i++) {
		if(done[i]) continue;
		if(i == idx[i]) {
			done[i] = 1;
			continue;
		}
//		memcpy(temp, database[i], FTR_DIM); temp_sub_fb = sub_fb[i];
		memcpy(temp, &ds->ftr_id[i], sizeof(struct_ftr_id));
		int j;
		for(j = i; i != idx[j]; j = idx[j]) {
//			memcpy(database[j], database[idx[j]], FTR_DIM); sub_fb[j] = sub_fb[idx[j]];
			memcpy(&ds->ftr_id[j], &ds->ftr_id[idx[j]], sizeof(struct_ftr_id)); 
			done[j] = 1;
		}
//		memcpy(database[j], temp, FTR_DIM); sub_fb[j] = temp_sub_fb;
		memcpy(&ds->ftr_id[j], temp, sizeof(struct_ftr_id));
		done[j] = 1;
	}
	ds->sorted = 1;
	free(temp);
	free(done);
}

// 複数の ftr ファイルをひとまとめに扱えるように open する．
// block_size > 1 のときは，まとめ読みに対応する．
struct_multi_ftr *open_multi_ftr(int num_ftr_files, char *ftr_filename[], int block_size)
{
	struct_multi_ftr *mf = (struct_multi_ftr *)MALLOC(sizeof(struct_multi_ftr));
	mf->num_ftr_files = num_ftr_files;
	mf->fh = (file_handle *)MALLOC(sizeof(file_handle) * num_ftr_files);
	mf->hd = (ftr_header_type *)MALLOC(sizeof(ftr_header_type) * num_ftr_files);
	mf->num_data = 0;
	mf->offset = (int *)MALLOC(sizeof(int) * num_ftr_files);
	for(int i = 0; i < num_ftr_files; i++) {
		if((mf->fh[i] = open_ftr_file(ftr_filename[i], &mf->hd[i])) == OPEN_ERROR) {
			fprintf(stderr, "cannot open ftr file = %s\n", ftr_filename[i]);
			exit(0);
		}
		mf->offset[i] = mf->num_data;
		mf->num_data += mf->hd[i].num_data;
	}
	mf->block_size = block_size;
	if(block_size <= 1) { // まとめ読みしないので，ftr_id は配列にしないでデータ1個分だけメモリ確保する．data_num は不要．
		mf->ftr_id = (struct_ftr_id *)MALLOC(sizeof(struct_ftr_id));
		mf->data_num = NULL;
	} else {
		mf->ftr_id = (struct_ftr_id *)MALLOC(sizeof(struct_ftr_id) * block_size);
		mf->data_num = (int *)MALLOC(sizeof(int) * block_size);
	}
	mf->read_in = 0;
	mf->next = 0;
	return mf;
}

void close_multi_ftr(struct_multi_ftr *mf)
{
	for(int i = 0; i < mf->num_ftr_files; i++) {
		CLOSE(mf->fh[i]);
	}
	free(mf->fh);
	free(mf->hd);
	free(mf->offset);
	if(mf->ftr_id != NULL) free(mf->ftr_id);
	if(mf->data_num != NULL) free(mf->data_num);
	free(mf);
}
/*
static int run(int a[], int lg)
{
	int i, r;
	for(r = 0, i = 0; i < lg - 1; i++, r++) {
		if(a[i] + 1 != a[i + 1]) break;
	}
	return r + 1;
}
*/
/*
static int select_file(struct_multi_ftr *mf, int data_num)
{
	if(data_num >= mf->num_data) {
		fprintf(stderr, "data_num = %d is out of range: number of data in muli ftr files = %d\n", data_num, mf->num_data);
		exit(0);
	}
	int right = mf->num_ftr_files;
	int left = 0;
//	if(data_num >= 102000000) { printf("num_ftr_files = %d, data_num = %d\n", mf->num_ftr_files, data_num); }
	while(left < right - 1) {
		int m = (left + right) / 2;
//	if(data_num >= 102000000) { printf("left = %d, right  %d, m = %d, offset = %d\n", left, right, m, mf->offset[m]); } // getchar();
		if(data_num < mf->offset[m]) {
			right = m;
		} else if(data_num > mf->offset[m]) {
			left = m;
		} else {
			return m;
		}
	}
	return left;
}
*/
// マルチftrでの通し番号から，ファイル番号を求める
static int select_file(struct_multi_ftr *mf, int data_num)
{
	int f;
	for(f = 0; f < mf->num_ftr_files - 1; f++) {
//		printf("offset[%d] = %d\n", f + 1, mf->offset[f + 1]);
		if(data_num < mf->offset[f + 1]) {
//			printf("selected file num = %d for data_num = %d\n", f, data_num); getchar();
//			if(f != select_file_2(mf, data_num)) {
//				fprintf(stderr, "data_num = %d, select_file = %d, select_file_2 = %d\n", data_num, f, select_file_2(mf,data_num));
//				exit(0);
//			}
			return f;
		}
	}
	if(data_num < mf->num_data) {
//		printf("selected file num = %d for data_num = %d\n", f, data_num); getchar();
//		if(f != select_file_2(mf, data_num)) {
//			fprintf(stderr, "data_num = %d, select_file = %d, select_file_2 = %d\n", data_num, f, select_file_2(mf,data_num));
//			exit(0);
//		}
		return f;
	}
	fprintf(stderr, "data_num = %d is out of range: number of data in muli ftr files = %d\n", data_num, mf->num_data);
	exit(0);
}


static int run(struct_multi_ftr *mf, int a[], int lg, int f)
{
	int i, r;
	for(r = 0, i = 0; i < lg - 1; i++, r++) {
		if(a[i] + 1 != a[i + 1]) break;
		if(select_file(mf, a[i + 1]) != f) break;
	}
	return r + 1;
}

int get_ftr_block(struct_multi_ftr *mf, int data_num_list[], int start, int list_length)
{
// 特徴データをブロック読み込みする．
// 読み込むデータのデータ番号（ftrファイル中の通し番号）は，data_num_list[start + b] (b = 0, 1, ... , block_size - 1)
// ただし，start + b < list_length であること．
// return value = number of records read in
// 0: error
	int b, lg;
//	fprintf(stderr, "get_ftr_block: start = %d, length = %d\n", start, list_length);
	lg = mf->block_size;
	if(start + mf->block_size >= list_length) { // start から始めて block_size 分を読み込むと範囲をはみ出してしまう
		lg = list_length - start; 
	}
	b = 0; // 読み込み済みデータ数
	while(b < lg) {
//	fprintf(stderr, "lg = %d, b = %d\n", lg, b);
		int f = select_file(mf, data_num_list[start + b]); // 次に読み込むデータ番号（data_num_list[start + b]）のファイル番号
		int r = run(mf, data_num_list + start + b, lg - b, f); // data_num_list[start + b] から連続読み込みが可能なデータ番号が何個続いているか
		int data_num_in_f = data_num_list[start + b] - mf->offset[f]; // 次に読み込むデータのファイル内のデータ番号
//		fprintf(stderr, "start = %d, b = %d, data_num = %d, file = %d, data_num_in_f = %d, r = %d\n", start, b, data_num_list[start + b], f, data_num_in_f, r);
//		getchar();
//		if(mf->next_data_num != data_num_list[start + b]) 	// seek が必要か？
//	#pragma omp critical
//	{
		SEEK(mf->fh[f], (long)data_num_in_f * sizeof(struct_ftr_id) + sizeof(ftr_header_type), SEEK_SET);
		if(!READ(mf->fh[f], mf->ftr_id + b, sizeof(struct_ftr_id) * r)) { // r 個まとめて読み込む
			fprintf(stderr, "cannot read ftr data (data_num = %d, file_num = %d, data_num_in_f = %d)\n", data_num_list[start + b], f, data_num_in_f);
			exit(0);
		}
//	}
//		fprintf(stderr, "read OK\n");
		for(int j = 0; j < r; j++) {
			mf->data_num[b + j] = data_num_list[start + b + j];
		}
		b += r;
	}
	mf->read_in = b;
	mf->next = 0;

	return b;
}

struct_ftr_id *get_next_ftr_id_from_multi_ftr(struct_multi_ftr *mf, int data_num_list[], int start, int list_length)
{
	struct_ftr_id *ftr_id;
	if(mf->block_size <= 1) {
		int f = select_file(mf, data_num_list[start]);
		int data_num_in_f = data_num_list[start] - mf->offset[f];
//		fprintf(stderr, "start = %d, data_num = %d, file = %d, data_num_in_f = %d\n", start, data_num_list[start], f, data_num_in_f);
		SEEK(mf->fh[f], (long)data_num_in_f * sizeof(struct_ftr_id) + sizeof(ftr_header_type), SEEK_SET);
		if(!READ(mf->fh[f], mf->ftr_id, sizeof(struct_ftr_id))) { // 1 個だけ読み込む
			fprintf(stderr, "cannot read ftr data (data_num = %d, file_num = %d, data_num_in_f = %d)\n", data_num_list[start], f, data_num_in_f);
			exit(0);
		}
		ftr_id = mf->ftr_id;
//		fprintf(stderr, "data_id = %ld\n", ftr_id->data_id);  getchar();
	} else {
//		fprintf(stderr, "get_next: start = %d, next = %d, read_in = %d\n", start, mf->next, mf->read_in);
		// ブロック内に読み込んだデータがなくなったら，つぎのブロックをバッファーに読み込む
		if(mf->next >= mf->read_in) {
			get_ftr_block(mf, data_num_list, start, list_length);
		}
		// バッファーの次のデータとデータ番号を受け取る
		ftr_id = &mf->ftr_id[mf->next];
		mf->next++;
	}
	return ftr_id;
}

void pre_rotation_ftr(ftr_type ftr, int *rotation)
{
	ftr_element_type temp[FTR_DIM];
	for(int i = 0; i < FTR_DIM; i++) {
		temp[i] = ftr[rotation[i]];
	}
	memcpy(ftr, temp, sizeof(ftr_element_type) * FTR_DIM);
}

/*
int open_ftr_file_n(int num_files, char *filename[], ftr_handle *handle[], int num_data[], int block_size)
// n 個の ftr ファイルを new_ftr_handle で open し，
// num_data[i] = sum(handle[0].header.num_data + ... handle[i].header.num_data), つまり，i 番目のファイルまでのデータ数の合計
{
	int n, num = 0;
	for(n = 0; n < num_files; n++) {
		handle[n] = new_ftr_handle(filename[n], block_size);
		num += handle[n].header.num_data;
		num_data[n] = num;
	}
	return num;
}

static int select_file(int num_files, int num_data[], int data_num)
{
	for(int f = 0; f < num_files; f++) {
		if(data_num < num_data[f]) {
			return f;
		}
	}
	fprintf(stderr, "data_num = %d is out of range\n", data_num);
	exit(0);
}

static int offset(int num_files, int num_data[], int f)
{
	return f == 0 ? 0 : num_data[f - 1];
}

void close_ftr_file_n(int num_files, ftr_handle *handle[])
{
	int n;
	for(n = 0; n < num_files; n++) {
		free_ftr_handle(handle[n]);
	}
}

struct_ftr_id *get_ftr_id_n(int num_files, ftr_handle *handle[], int num_data[], int data_num)
// 特徴データを読む
// 1: no error
// 0: error
{
	int f = select_file(num_files, num_data, data_num);
	int data_num_in_f = data_num - offset(num_files, num_data, f);

	return get_next_ftr_id_from_block(handle[f], 

	SEEK(fh_n[f], sizeof(ftr_header_type) + (long)idx * (db_header->element_size * db_header->data_dim + db_header->sub_fb), SEEK_SET);
	return READ(fh_n[f], ftr, db_header->element_size * db_header->data_dim); 
}
*/
