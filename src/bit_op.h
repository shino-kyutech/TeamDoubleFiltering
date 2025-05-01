// ビット操作
void print_bin_w(unsigned int s, int w);
void print_bin(unsigned int s);
void print_bin_err(unsigned int s);
void write_bit(int offset, int on_off, unsigned int *d);
void print_bin_long(unsigned long s);
void write_bit_long(int offset, int on_off, unsigned long *d);
void print_bin_expanded(unsigned long *s, int lg);
void print_bin_expanded_2(unsigned long *s, int lg);
void write_bit_expanded(int offset, int on_off, unsigned long *d);
int msb_pos(unsigned int x);
int lsb_pos(unsigned int x);
int msb_pos_long(unsigned long x);
int lsb_pos_long(unsigned long x);
void make_bitcnt_tbl(int width);
int bit_count(unsigned int a);
int bit_count_long(unsigned long a);
int bit_count_expanded(unsigned long *a, int lg);
int bit_check(unsigned a, int pos);
int bit_check_long(unsigned long a, int pos);
int bit_check_expanded(unsigned long a[], int pos);
int comp_bit(const void *a, const void *b);
void make_bitnum_pat(int width);

// 2-bit 整数などの操作
// typedef struct {
//     unsigned int *mem;  // w-bit 整数を詰め合わせて格納する領域
//     int num, w;         // 要素数，整数が使用するビット数
// } tiny_int;

#ifdef TINY_IN_CHAR  // tiny_int を 8-bit 単位で詰め合わせる．右寄せで，端数が出ても使わない．
typedef unsigned char tiny_int;
#else
typedef unsigned int tiny_int;
#endif
tiny_int *new_tiny_int(int num, int w);
void char2tiny(unsigned char a[], tiny_int p[], int num, int w);
void modify_tiny_int(tiny_int mem[], int i, int w, unsigned char v);
#ifdef QUANTIZE_MIXED_MOD3
void char2tiny_mod3(unsigned char a[], tiny_int mem[], int num);
#endif
void tiny2char(tiny_int p[], unsigned char a[], int num, int w);
