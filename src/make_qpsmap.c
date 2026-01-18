#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "config.h"
#include "ftr.h"
#include "smap.h"
#include "sketch.h"
#include "bit_op.h"
#include "e_time.h"

int main(int argc, char *argv[]) {
	int num_ftr_files = argc - 1;               // パラメータは，データセットのファイル群 (ds_00.ftr, ds_01.ftr, ... )
	char **dataset_ftr_file = argv + 1;
    char *smap_pivot_file = SMAP_PIVOT_FILE;    // SMAP用のピボット
    char *bucket_file = BUCKET_FILE;            // スケッチで作ったバケット表
    char *qpsmap_file = QPSMAP_FILE;            // 作成したQPSMAPを保存するファイル

	make_bitcnt_tbl(8); // bit_countを使うならば，最初に呼び出しておくこと．

   	fprintf(stderr, "Start make_qpsmap: PJT_DIM = %d, SMAP_DIM = %d\n", PJT_DIM, SMAP_DIM);    // PJT_DIM は不要かな？　SMAP_DIMは？
	use_system("VmSize");

    FILE *fp_out = fopen(qpsmap_file, "w");
    if(!fp_out) {
        perror("fopen");
        return 1;
    }
    fprintf(stderr, "Fopen QPSMAP file OK.\n");

    // バケット表の読み込み
    struct_bucket *bucket_ds = read_bucket(bucket_file);            // 変数名は bucket_ds ではなく，bucket　とした方が自然？　（後で確認）
    int num_data = bucket_ds->num_data;
	fprintf(stderr, "Read bucket OK: number of data = %d\n", num_data);
	use_system("VmSize");

    // QPSMAP ファイルのヘッダ
    qpsmap_header hd;
    hd.smap_dim = SMAP_DIM;
    hd.quantize_bit = QUANTIZE_BIT;
    hd.num_data = num_data;

    // 複数に分かれた特徴データファイル対応の読込の準備（qpsmap像を作るときと最終段階での実距離計算による検索のときに使用）
    fprintf(stderr, "open ftr files (filename = %s, num_files = %d)\n", dataset_ftr_file[0], num_ftr_files);
    #if !defined(_OPENMP) || NUM_THREADS < 1
    struct_multi_ftr *mf = open_multi_ftr(num_ftr_files, dataset_ftr_file, BLOCK_SIZE);
    if(num_data != mf->num_data) {
        fprintf(stderr, "bucket (filename = %s, num_data = %d) is not compatible with datasets (filename = %s, ... , num_data = %d)\n", 
            bucket_file, num_data, dataset_ftr_file[0], mf->num_data);
        return -1;
    }
    #else

    #endif

	dataset_handle dh;

    #ifdef _OPENMP
		dh.num_threads = NUM_THREADS;
	#else
		dh.num_threads = 1;
	#endif

    dh.sorted = 0; // FTR is NOT sorted. Arrangement is as is for double filtering 
	dh.ftr_on = SECONDARY_MEMORY;
	dh.mf = (struct_multi_ftr **)malloc(sizeof(struct_multi_ftr *) * dh.num_threads);
	for(int t = 0; t < dh.num_threads; t++) {
		dh.mf[t] = open_multi_ftr(num_ftr_files, dataset_ftr_file, BLOCK_SIZE);
	}
	dh.ds = NULL;
	if(num_data != dh.mf[0]->num_data) {
		fprintf(stderr, "bucket (filename = %s, num_data = %d) is not compatible with datasets (filename = %s, ... , num_data = %d)\n", bucket_file, num_data, dataset_ftr_file[0], dh.mf[0]->num_data);
		return -1;
	}

    // QPSMAPのピボット読込み
    smap_pivot_type *smap_pivot = new_smap_pivot(PQBP);
    read_smap_pivot(smap_pivot_file, smap_pivot);

    // QPSMAP 用パラメタ計算
    double ave[SMAP_DIM], stdev[SMAP_DIM], offset[SMAP_DIM], slice[SMAP_DIM];
    struct_dataset *ds_sample = read_dataset_n(1, &dataset_ftr_file[0]); // 先頭のデータセットフィルをサンプルとして用いる
    compute_parmeters_for_qpsmap(smap_pivot, ds_sample, ave, stdev, offset, slice);
    free_dataset(ds_sample);

    // QPSMAP 作成用メモリ確保
    tiny_int (*packed_qpsmap_data)[PACKED_QPSMAP_SIZE];
    packed_qpsmap_data = malloc(sizeof(tiny_int) * PACKED_QPSMAP_SIZE * num_data);
	fprintf(stderr, "malloc packed quantized images of data OK.\n");
	use_system("VmSize");

    // 元のデータ番号
    int *data_num_org = malloc(sizeof(int) * num_data);
    for(int i = 0; i < num_data; i++) data_num_org[i] = i;

    // QPSMAP 作成
    #ifdef _OPENMP
    omp_set_num_threads(NUM_THREADS);
    #pragma omp parallel for
    #endif
    for(int i = 0; i < num_data; i++) {
        unsigned char quantized_projected_data[SMAP_DIM];
        #ifndef _OPENMP
        int t = 100;
        int nd = num_data;
 		struct_ftr_id *ftr_id_p = get_next_ftr_id_from_multi_ftr(mf, data_num_org, i, num_data);
        #else
        int t = omp_get_thread_num();
        int nd = num_data / NUM_THREADS;
 		struct_ftr_id *ftr_id_p = get_next_ftr_id_from_multi_ftr(dh.mf[t], data_num_org, i, num_data);
        #endif
        q_uchar_psmap(quantized_projected_data, ftr_id_p->ftr, smap_pivot, offset, slice);
        char2tiny(quantized_projected_data, packed_qpsmap_data[i], SMAP_DIM, QUANTIZE_BIT);

        if(t == 0 && (i + 1) % (nd / 100) == 0) {
            fprintf(stderr, "packed qpsmap (t = %d) %4d%%\r", t, (i + 1) / (nd / 100));
        }
    }
    fprintf(stderr, "\nQPSMAP creation done.\n");

    // スケッチ順に並べ替え
    fprintf(stderr, "\nQPSMAP sorting ... ");
    tiny_int *temp = MALLOC(sizeof(tiny_int) * PACKED_QPSMAP_SIZE);
    char *done = calloc(num_data, sizeof(char));
    for(int i = 0; i < num_data; i++) {
        if(done[i] || i == bucket_ds->idx[i]) {
            done[i] = 1;
            continue;
        }
        memcpy(temp, packed_qpsmap_data[i], sizeof(tiny_int) * PACKED_QPSMAP_SIZE);
        int j;
        for(j = i; i != bucket_ds->idx[j]; j = bucket_ds->idx[j]) {
            memcpy(packed_qpsmap_data[j], packed_qpsmap_data[bucket_ds->idx[j]], sizeof(tiny_int) * PACKED_QPSMAP_SIZE);
            done[j] = 1;
        }
        memcpy(packed_qpsmap_data[j], temp, sizeof(tiny_int) * PACKED_QPSMAP_SIZE);
        done[j] = 1;
    }
	FREE(temp, sizeof(tiny_int) * PACKED_QPSMAP_SIZE);
    free(done);
    fprintf(stderr, "done.\n");

    // バイナリファイルに保存
    write_qpsmap(&hd, offset, slice, (tiny_int *)packed_qpsmap_data, num_data, fp_out);
    fclose(fp_out);

    fprintf(stderr, "QPSMAP saved to \"%s\".\n", qpsmap_file);

    free(packed_qpsmap_data);
    free(data_num_org);

    return 0;
}
