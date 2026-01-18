// QPSMAP ファイルのヘッダー情報を得る

#include <stdio.h>
#include <string.h>
// QPSMAP ファイルの先頭
typedef struct {
	int smap_dim;			// 射影次元数
	int quantize_bit;		// 量子化ビット数
	int num_data;			// データ数
} qpsmap_header;

int main(int argc, char *argv[]) {
    if(argc != 3 && strcmp(argv[2], "DIM") && strcmp(argv[2], "BIT") && strcmp(argv[2], "NUM")) {
        fprintf(stderr, "usage> %s <qpsmap_file> [ DIM | BIT | NUM ]\n", argv[0]);
        fprintf(stderr, "argv[1] = %s\n", argv[1]);
        fprintf(stderr, "argv[2] = %s\n", argv[2]);
        return 1;
    }

    char *qpsmap_file = argv[1];
    FILE *fp = fopen(qpsmap_file, "r");
    if(!fp) {
        fprintf(stderr, "Cannot fopen %s\n", qpsmap_file);
        return 1;
    }
    qpsmap_header hd;
    if(fread(&hd, sizeof(qpsmap_header), 1, fp) != 1) {
        fprintf(stderr, "fread error\n");
        return 1;
    }
	int smap_dim;			// 射影次元数
	int quantize_bit;		// 量子化ビット数
	int num_data;			// データ数

    char *name = argv[2];
    if(!strcmp(name, "DIM")) {
        printf("%d\n", hd.smap_dim);
    } else if(!strcmp(name, "BIT")) {
        printf("%d\n", hd.quantize_bit);
    } else if(!strcmp(name, "NUM")) {
        printf("%d\n", hd.num_data);
    }

    return 0;
}
