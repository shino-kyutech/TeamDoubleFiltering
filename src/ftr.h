#pragma once
// 特徴データに関する型定義とデータセットファイルの入出力
#include "parm.h"
#include "config.h"

// dataset setting (for Image dataset, #define DATASE IMAGE and switch using DATASET macro)
// typedef enum {IMAGE, DECAF, DEEP1B} dataset;

typedef unsigned int dist_type;

#if defined(DECAF)
// 特徴データの型
typedef unsigned char ftr_element_type, *ftr_type; 
#ifndef FTR_DIM
#define FTR_DIM 4096
#endif
#define FTR_MIN 0
#define FTR_IGN 127 // SQBP のときに無視する
#define FTR_MAX 255
typedef long auxi_type; // 特徴データの付随情報の型: DISA_h の 100M データセットでは long (おそらく，なんらかのデータ識別番号)
#elif defined(DEEP1B)
// 特徴データの型
typedef signed char ftr_element_type, *ftr_type; 
#ifndef FTR_DIM
#define FTR_DIM 96
#endif
#ifndef FTR_MIN
#define FTR_MIN -127
#endif
#define FTR_IGN 0 // SQBP のとき無視する
#ifndef FTR_MAX
#define FTR_MAX 127
#endif
typedef unsigned int auxi_type; // 特徴データの付随情報の型: Deep1B では unsigned int (元のデータにはないので，データ変換時に通し番号を割振った)
#elif defined(LAION2B)
// 特徴データの型
typedef signed char ftr_element_type, *ftr_type; 
#ifndef FTR_DIM
#define FTR_DIM 768
#endif
#ifndef FTR_MIN
#define FTR_MIN -127
#endif
#define FTR_IGN 0 // SQBP のとき無視する
#ifndef FTR_MAX
#define FTR_MAX 127
#endif
typedef unsigned int auxi_type; // 特徴データの付随情報の型: Deep1B では unsigned int (元のデータにはないので，データ変換時に通し番号を割振った)
#elif defined(PUBMED23)
// 特徴データの型
typedef signed char ftr_element_type, *ftr_type; 
#ifndef FTR_DIM
#define FTR_DIM 384
#endif
#ifndef FTR_MIN
#define FTR_MIN -127
#endif
#define FTR_IGN 0 // SQBP のとき無視する
#ifndef FTR_MAX
#define FTR_MAX 127
#endif
typedef unsigned int auxi_type; // 特徴データの付随情報の型: Deep1B では unsigned int (元のデータにはないので，データ変換時に通し番号を割振った)
#endif

// 特徴データファイルのヘッダ情報
typedef struct {
	int element_size;	// 各次元のサイズ, つまり，正しければ，sizeof(ftr_element_type)
	int data_dim;		// 特徴次元数 (FTR_DIMのはず) 
	int num_data;		// データ数
	int auxi_size;		// 付随情報）のサイズ sizeof(auxi_type)
} ftr_header_type;

// 特徴データと付随情報（データID）をまとめた構造体（DeCAF記述子用）
typedef struct {
	ftr_element_type ftr[FTR_DIM];
	auxi_type data_id;
} struct_ftr_id;

// データセット（特徴データ）
//	int num_data;			// 特徴データ数
//	struct_ftr_id *ftr_id;	// 特徴データ（データID付き）の配列
//	int sorted;				// 特徴データがスケッチ順にソートされているかどうか
typedef struct {
	int num_data;			// 特徴データ数
	struct_ftr_id *ftr_id;	// 特徴データ（データID付き）の配列
	int sorted;				// 特徴データがスケッチ順にソートされているかどうか
} struct_dataset;

// 特徴データのファイル入出力

typedef enum {LOW_LEVEL, HIGH_LEVEL} file_io;

#ifndef FILE_IO
#define FILE_IO HIGH_LEVEL 
#endif

#if FILE_IO == LOW_LEVEL

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>

typedef int file_handle;

#define READ_OPEN(fn) open(fn, O_RDONLY)
#define OPEN_ERROR (-1)
#define CLOSE(fh) close(fh)
#define READ(fh, bf, bs) (read(fh, bf, bs) == (bs))
#define SEEK(fh, offset, origin) lseek(fh, offset, origin)

#else // HIGH_LEVEL

typedef FILE *file_handle;

#define READ_OPEN(fn) fopen(fn, "rb")
#define OPEN_ERROR NULL
#define CLOSE(fh) fclose(fh)
#define READ(fh, bf, bs) (fread(bf, bs, 1, fh) == 1)
#define SEEK(fh, offset, origin) fseek(fh, offset, origin)

#endif // FILE_IO

// データセット群（複数のデータセットをひとまとめに取り扱う）（まとめ読み）
typedef struct {
	int num_ftr_files;		// 特徴データファイル数
	file_handle *fh;		// ハンドルの配列
	ftr_header_type *hd;	// ヘッダの配列
	int num_data;			// 特徴データ数の総合計
	int *offset;			// offset[i] = i 番目のデータセットの先頭のデータの通し番号
	int block_size;			// まとめ読みの最大個数
	struct_ftr_id *ftr_id;	// 特徴データ（データID付き）のバッファ
	int *data_num;			// 読み込む前: ftr_id に読み込むデータ番号列，読み込み後: 読み込まれているデータ番号列
	int read_in;			// ブロックリードで読み込んだデータ数
	int next;				// つぎに返すデータの添え字（バッファの添え字）
} struct_multi_ftr;

// オンメモリのときは srtuct_dataset，2次記憶に置いたままのときは，struct_multi_ftr を使用する．
typedef struct {
	ftr_on_type ftr_on; // 0 -> main memory, 1 -> secondary
	int sorted;
	struct_dataset *ds;
	int num_threads;		// 並列アクセスのスレッド数
	struct_multi_ftr **mf;	// struct_mlti_ftr の配列（要素数はスレッド数） 
} dataset_handle;

void print_header(FILE *f, char *filename, ftr_header_type *header);
file_handle open_ftr_file(char *filename, ftr_header_type *header);
int get_ftr_id(file_handle fh, ftr_header_type *header, int data_num, struct_ftr_id *ftr_id);
int get_ftr_id_bulk(file_handle fh, ftr_header_type *header, int data_num, int num_data, struct_ftr_id *ftr_id);
struct_dataset *read_dataset_n(int num_ftr_files, char *filename[]);
void free_dataset(struct_dataset *ds);
void sort_dataset(struct_dataset *ds, int idx[]);

struct_multi_ftr *open_multi_ftr(int num_ftr_files, char *ftr_filename[], int block_size);
void close_multi_ftr(struct_multi_ftr *mf);
int get_ftr_block(struct_multi_ftr *mf, int data_num_list[], int start, int list_length);
struct_ftr_id *get_next_ftr_id_from_multi_ftr(struct_multi_ftr *mf, int data_num_list[], int start, int list_length);

void pre_rotation_ftr(ftr_type ftr, int *rotation);
