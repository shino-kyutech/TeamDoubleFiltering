#include <stdio.h>
#include <stdlib.h>
#include "e_time.h"
#include "bit_op.h"
#include "parm.h"

#ifdef DEBUG
int main(void)
{
	int num = 20, w = 3;
	tiny_int *packed = new_tiny_int(num, w);
	unsigned char a[num];
	srandom(3);
	for(int i = 0; i < num; i++) {
		a[i] = random() % (1<<w);
		printf("%3u", a[i]);
	}
	printf("\n");
	char2tiny(a, packed, num, w);
	for(int i = 0; i < 10; i++) {
		print_bin_w(packed[i], 8); printf(" ");
	}
	printf("\n");
	unsigned char b[num];
	tiny2char(packed, b, num, w);
	for(int i = 0; i < num; i++) {
		printf("%3d", b[i]);
	}
	printf("\n");
	modify_tiny_int(packed, 12, w, 5);
	tiny2char(packed, b, num, w);
	for(int i = 0; i < num; i++) {
		printf("%3d", b[i]);
	}
	printf("\n");
	modify_tiny_int(packed, 5, w, 3);
	tiny2char(packed, b, num, w);
	for(int i = 0; i < num; i++) {
		printf("%3d", b[i]);
	}
	printf("\n");
	modify_tiny_int(packed, 6, w, 6);
	tiny2char(packed, b, num, w);
	for(int i = 0; i < num; i++) {
		printf("%3d", b[i]);
	}
	printf("\n");
	modify_tiny_int(packed, 6, w, 0);
	tiny2char(packed, b, num, w);
	for(int i = 0; i < num; i++) {
		printf("%3d", b[i]);
	}
	printf("\n");
/*
	unsigned long i;
	unsigned long a, b, c, d;
	unsigned long e[3] = {0}, f[3] = {0xa0a0a0a0a0a0a0a0, 0xffffffffffffffff, 0x5050505050505050};

	for(int i = 0; i < 192; i++) {
		printf("i = %d, bit(e) = %d, bit(f) = %d\n", i, bit_check_expanded(e, i), bit_check_expanded(f, i));
	}

	make_bitcnt_tbl(8);
	for(int j = 0; j < 5; j++) {
		a = random() & 0xffff; b = random() & 0xffff; c = random() & 0xffff; d = random() & 0xffff;
		i = a | (b << 16) | (c << 32) | (d << 48);
		printf("i = %20lu (", i); print_bin_long(i); printf("), bit_count = %2d\n", bit_count_long(i));
		e[0] = i;
		a = random() & 0xffff; b = random() & 0xffff; c = random() & 0xffff; d = random() & 0xffff;
		i = a | (b << 16) | (c << 32) | (d << 48);
		printf("i = %20lu (", i); print_bin_long(i); printf("), bit_count = %2d\n", bit_count_long(i));
		e[1] = i;
		printf("####");
		print_bin_expanded(e, 2);
		printf("####\n");
		printf("bit_count_expand = %d\n", bit_count_expanded(e, 2));
//		int on = random() % 128;
		int on = 127;
		printf("on = %d\n", on);
		write_bit_expanded(0, 0, f);
		write_bit_expanded(on, 0, f);
		printf("####");
		print_bin_expanded(f, 2);
		printf("####\n");
		write_bit_expanded(on, 1, e);
		printf("####");
		print_bin_expanded(e, 2);
		printf("####\n");
		write_bit_expanded(on, 0, e);
		printf("####");
		print_bin_expanded(e, 2);
		printf("####\n");
	}
*/	
	return 0;
}
#endif

// ビット操作
void print_bin_w(unsigned int s, int w)
{
	int i;
	for(i = 0; i < w; i++)
		printf("%1d",(s >> (w - 1 - i)) & 1);
}

void print_bin(unsigned int s)
{
	int i;
	for(i = 0; i < 32; i++)
		printf("%1d",(s >> (31 - i)) & 1);
}

void print_bin_err(unsigned int s)
{
	int i;
	for(i = 0; i < 32; i++)
		fprintf(stderr, "%1d",(s >> (31 - i)) & 1);
}

void write_bit(int offset, int on_off, unsigned int *d)
{
	if(on_off == 1) {
		*d = *d | (1 << offset);
	} else {
		*d = *d & ~(1 << offset);
	}
}

// ビット操作（long用）
void print_bin_long(unsigned long s)
{
	int i;
	for(i = 0; i < 64; i++)
		printf("%1ld",(s >> (63 - i)) & 1);
}

void write_bit_long(int offset, int on_off, unsigned long *d)
{
	if(on_off == 1) {
		*d = *d | ((unsigned long)1 << offset);
	} else {
		*d = *d & ~((unsigned long)1 << offset);
	}
}

// ビット操作（expanded用）
void print_bin_expanded(unsigned long *s, int lg)
{
	for(int i = 0; i < lg; i++)
		print_bin_long(s[i]);
}

void print_bin_expanded_2(unsigned long *s, int lg)
{
	for(int i = 0; i < lg; i++) {
		if(i != 0) printf(":");
		print_bin_long(s[i]);
	}
}

void write_bit_expanded(int offset, int on_off, unsigned long *d)
{
	int pos = offset / 64;
	offset = offset - pos * 64;
	if(on_off == 1) {
		d[pos] = d[pos] | ((unsigned long)1 << offset);
	} else {
		d[pos] = d[pos] & ~((unsigned long)1 << offset);
	}
}
int msb_pos(unsigned int x) {
  int pos = -1;
 
  if (x != 0) {
    __asm__("bsrl %1, %0": "=r" (pos): "m" (x));
  }
 
  return pos;
}

int lsb_pos(unsigned int x) {
  int pos = -1;
 
  if (x != 0) {
    __asm__("bsfl %1, %0": "=r" (pos): "m" (x));
  }
 
  return pos;
}

int msb_pos_long(unsigned long x) {
  long pos = -1;
 
  if (x != 0) {
    __asm__("bsrq %1, %0": "=r" (pos): "m" (x));
  }
 
  return pos;
}

int lsb_pos_long(unsigned long x) {
  long pos = -1;
 
  if (x != 0) {
    __asm__("bsfq %1, %0": "=r" (pos): "m" (x));
  }
 
  return pos;
}

// ビットカウントのため（ハミング距離で使用する）
int *bitcnt_tbl = NULL;

// width = 8 で十分
void make_bitcnt_tbl(int width)
{
	int BIT = (1 << width);
	int i;
	unsigned int v;

	if(bitcnt_tbl == NULL) {
		bitcnt_tbl = (int *)malloc(sizeof(int) * BIT);
	}
	for(i = 0; i < BIT; i++) {
		v = i;
		for (bitcnt_tbl[i] = 0; v; v >>= 1) {
			bitcnt_tbl[i] += v & 1;
		}
	}
}

int bit_count(unsigned int a)
{
	int c = 0;

	// if(sizeof(sketch_type) == sizeof(short))
	//	c = bitcnt_tbl[a]; // 16bit
	// else {
		//for(; a; a >>= 8)  // 64bit等にも対応　ただし遅い
		//	c += bitcnt_tbl[a & 0xff]; 
		c = bitcnt_tbl[a & 0xff] + bitcnt_tbl[(a >> 8) & 0xff] + bitcnt_tbl[(a >> 16) & 0xff] + bitcnt_tbl[(a >> 24) & 0xff]; // 32bit
	//}
	
	return c;
}

int bit_count_long(unsigned long a)
{
	int c = 0;

	c = bitcnt_tbl[a & 0xff] + bitcnt_tbl[(a >> 8) & 0xff] + bitcnt_tbl[(a >> 16) & 0xff] + bitcnt_tbl[(a >> 24) & 0xff]
	  + bitcnt_tbl[(a >> 32) & 0xff] + bitcnt_tbl[(a >> 40) & 0xff] + bitcnt_tbl[(a >> 48) & 0xff] + bitcnt_tbl[(a >> 56) & 0xff]; // 64bit

	return c;
}

int bit_count_expanded(unsigned long *a, int lg)
{
	int c = 0;
	for(int i = 0; i < lg; i++) {
		c += bit_count_long(a[i]);
	}

	return c;
}

int comp_bit(const void *a, const void *b)
{
	if(bitcnt_tbl[*((int *) a)] < bitcnt_tbl[*((int *) b)])
		return -1;
	else if(bitcnt_tbl[*((int *) a)] == bitcnt_tbl[*((int *) b)])
		return 0;
	else
		return 1;
}

int bit_check(unsigned a, int pos)
{
	unsigned m = 1;
	return (a & (m << pos) ? 1 : 0);
}

int bit_check_long(unsigned long a, int pos)
{
	unsigned long m = 1;
	return (a & (m << pos) ? 1 : 0);
}

int bit_check_expanded(unsigned long a[], int pos)
{
	int n = pos / 64;
	pos = pos - n * 64;
	unsigned long m = 1;
	return (a[n] & (m << pos) ? 1 : 0);
}

// ハミング距離順の列挙で使用するときに使用する
// ON_BITの昇順に並べておく
int *bitnum_pat = NULL;

void make_bitnum_pat(int width)
{
	int BIT = (1 << width);
	int i;

	if(bitnum_pat == NULL) {
		bitnum_pat = (int *)malloc(sizeof(int) * BIT);
	}
	for(i = 0; i < BIT; i++) {
		bitnum_pat[i] = i;
	}
	qsort(bitnum_pat, BIT, sizeof(int), comp_bit);
}

#ifdef TINY_IN_CHAR
// w-bit 整数 を num 個分の領域を確保する．ただし，8-bit 毎に端数の詰め合わせをしない．
// TABLE_UNIT が定義されているときは，TABLE_UNIT 個まで詰め合わせる．
tiny_int *new_tiny_int(int num, int w)
{
	#ifdef TABLE_UNIT
		int capa = 8 / w; // 8-bit に格納できる w-bit 整数の個数 
		if(capa < TABLE_UNIT) {
			fprintf(stderr, "too large tiny int or table unit: tiny int * table unit should be smaller than 8-bit\n");
			exit(0);
		}
		return MALLOC((num / TABLE_UNIT * sizeof(unsigned char)));
	#else
		int capa = 8 / w; // 8-bit に格納できる w-bit 整数の個数 
		if(capa == 0) {
			fprintf(stderr, "too large tiny int: tiny int should be smaller than 8-bit\n");
			exit(0);
		}
		return MALLOC((num / capa * sizeof(unsigned char)));
	#endif
}

// 8-bit 整数 a[0], ... , a[num - 1] を，w-bit 整数として詰め合わせて，unsigned char の mem に格納する．
// 8 % w != 0 のとき，つまり，端数が残っても詰め合わせしないで，つぎに回す．
// TABLE_UNIT が定義されているときは，TABLE_UNIT 個まで詰め合わせる．
void char2tiny(unsigned char a[], tiny_int mem[], int num, int w)
{
//	static int first = 1;
//	if(first) {
//		fprintf(stderr, "char2tiny: TINY_IN_CHAR, num = %d, w = %d\n", num, w); first = 0;
//	}
	int capa = (8 / w);
	#ifdef TABLE_UNIT
	if(capa < TABLE_UNIT) { // このチェックは別の所で行うようにした方がよい．
		fprintf(stderr, "too large TABLE_UNIT %d (QUANTIZE_BIT = %d)\n", TABLE_UNIT, QUANTIZE_BIT); 
	}
	capa = TABLE_UNIT;
	#endif

	for(int i = 0, j = 0; i < num; j++) {
		mem[j] = a[i++]; // 最初のものを格納する．
		for(int m = 1; m < capa && i < num; m++) { // 8-bit に入る個数（8 / w）だけ，ただし，全体で num 個まで（i < num）
			mem[j] = (mem[j] << w) | a[i++];
		}
	}
}

// 8-bit 整数 a[0], ... , a[num - 1] を，w-bit 整数として詰め合わせて，unsigned char の mem に格納されている．
// そのとき，a[i] を v に修正する．
void modify_tiny_int(tiny_int mem[], int i, int w, unsigned char v)
{
	int capa = (8 / w);	//											(例) w = 2 のときは，capa = 4
	#ifdef TABLE_UNIT
//		if(capa < TABLE_UNIT) { // このチェックは別の所で行うようにした方がよい．
//		fprintf(stderr, "too large TABLE_UNIT %d (QUANTIZE_BIT = %d)\n", TABLE_UNIT, QUANTIZE_BIT); 
//	}
	capa = TABLE_UNIT;
	#endif
	int pos = i / capa; // a[i] が格納されているのが mem[pos]		(例) i = 6 のときは，pos = 1
	int off = i % capa; // a[i] が mem[pos] の中で何番目のものか	(例) i = 6 のときは, off = 2
//	printf("i = %d, w = %d, capa = %d, pos = %d, off = %d, v = %d\n", i, w, capa, pos, off, v);
	unsigned char mask = ((1 << w) - 1) << ((capa - off - 1) * w); // 修正するビット位置 00001100
//	printf("mem  = "); print_bin_w(mem[pos], 8); printf("\n");
//	printf("mask = "); print_bin_w(mask, 8); printf("\n");
//	printf("mem & ~mask = "); print_bin_w(mem[pos] & ~mask, 8); printf("\n");
//	printf("(v << (capa - off - 1) * w = "); print_bin_w((v << (capa - off - 1) * w), 8); printf("\n");
	mem[pos] = (mem[pos] & ~mask) | (v << (capa - off - 1) * w);
//	printf("mem  = "); print_bin_w(mem[pos], 8); printf("\n");
}

#ifdef QUANTIZE_MIXED_MOD3
// 8-bit 整数 a[0], ... , a[num - 1] を，QUANTIZE_BIT_0, QUANTIZE_BIT_1, QUANTIZE_BIT_2 で指定するビット幅の整数として詰め合わせて，
// unsigned char の mem に格納する．
void char2tiny_mod3(unsigned char a[], tiny_int mem[], int num)
{
	for(int i = 0, j = 0; i < num; j++) {
		mem[j] = a[i++]; // 最初のものを格納する．
		mem[j] = (mem[j] << QUANTIZE_BIT_1) | a[i++];
		mem[j] = (mem[j] << QUANTIZE_BIT_2) | a[i++];
	}
}
#endif
// char2tiny を用いて　mem に格納された num 個の w-bit 整数を uchar a[0], ... , a[num - 1] に展開する．
void tiny2char(tiny_int mem[], unsigned char a[], int num, int w)
{
	unsigned char mask = (1 << w) - 1; // mask = 0001111 （下位のw-bitだけが1）
	for(int i = 0, j = 0; i < num; j++) {
		int mm = 8 / w < num - i ? 8 / w : num - i;
		a[i++] = mem[j] >> (mm - 1) * w;
		for(int m = 1; m < mm; m++) {
			a[i++] = mem[j] >> (mm - m - 1) * w & mask;
		}
	}
}
#else // TINY_IN_CHAR
// w-bit 整数 を num 個分の領域を確保する．
tiny_int *new_tiny_int(int num, int w)
{
	return MALLOC((num * w + sizeof(int) - 1) / sizeof(int));
}

// 8-bit 整数 a[0], ... , a[num - 1] を，w-bit 整数として詰め合わせて，mem に格納する．
void char2tiny(unsigned char a[], tiny_int mem[], int num, int w)
{
//	static int first = 1;
//	if(first) {
//		fprintf(stderr, "char2tiny: !TINY_IN_CHAR, num = %d, w = %d\n", num, w); first = 0;
//	}
//	int num = p->num, w = p->w;
//	unsigned int *mem = p->mem;
	unsigned int buff = 0, mask = (1 << w) - 1; // mask = 0001111 （下位のw-bitだけが1）
	int used = 0, j = 0;

	for(int i = 0; i < num; i++) {
		if(used + w > sizeof(int) * 8) { // つぎに w-bit を入れると，はみ出す
			int r = sizeof(int) * 8 - used; // 残っているビット数
			buff = (buff << r) | ((a[i] & mask) >> (w - r));
			mem[j++] = buff;
			buff = 0;
			buff = buff | (a[i] & ((1 << (w - r)) - 1));
			used = w - r;
		} else {
			buff = (buff << w) | (a[i] & mask);
			used += w;
			if(used == sizeof(int) * 8) { // ちょうど int 1個分（32-bit）になった．
				mem[j++] = buff;
				buff = 0;
				used = 0;
			}
		}
	}
	if(used) {
		int r = sizeof(int) * 8 - used; // 残っているビット数
		mem[j] = (buff << r); // 左づめにする．
	}
}

// mem に格納された num 個の w-bit 整数を uchar a[0], ... , a[num - 1] に展開する．
void tiny2char(tiny_int mem[], unsigned char a[], int num, int w)
{
//	int num = p->num, w = p->w;
//	unsigned int *mem = p->mem;
	unsigned int buff = 0, mask = (1 << w) - 1;
	int j = (num * w + sizeof(int) * 8 - 1) / (sizeof(int) * 8) - 1;
	buff = mem[j];
	int r = num * w % (sizeof(int) * 8);
	int x = sizeof(int) * 8;

	if(r) { // 端数がある
		buff = buff >> (x - r);
		x = r;
	}
	for(int i = num - 1; i >= 0; i--) {
		if(x < w) { // つぎの範囲にまたがっている．
			a[i] = buff;
			buff = mem[--j];
			a[i] = a[i] | ((buff & ((1 << (w - x)) - 1)) << x);
			buff = buff >> (w - x);
			x = sizeof(int) * 8 - (w - x);
		} else {
			a[i] = buff & mask;
			x -= w;
			if(x == 0) {
				buff = mem[--j];
				x = sizeof(int) * 8;
			} else {
				buff = buff >> w;
			}
		}
	}
}
#endif
