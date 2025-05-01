#include <stdio.h>
#include "sketch.h"
#include "quick.h"
#include <stdlib.h>
#include <limits.h>
#include <omp.h>

#ifdef DEBUG
#include <time.h>
#define NUM_DATA 99990000
//#define NUM_DATA 10000

int topidx[NUM_DATA], num_top, max_top, num_rest, min_rest;

int main(void)
{
	int num_data = NUM_DATA;
	unsigned *sc = (unsigned *)malloc(sizeof(unsigned) * (num_data + 1));
	int *idx = (int *)malloc(sizeof(int) * (num_data + 1));
	int *idx2 = (int *)malloc(sizeof(int) * (num_data + 1));
	int *idx3 = (int *)malloc(sizeof(int) * (num_data + 1));
	int *selected = (int *)malloc(sizeof(int) * (num_data + 1));
	int *selected2 = (int *)malloc(sizeof(int) * (num_data + 1));
	struct timespec tp1, tp2;
	long sec, nsec;
	#ifndef SEED
	#define SEED 1
	#endif
	srandom(SEED);
	printf("SEED = %d\n", SEED);
	for(int i = 0; i < num_data; i++) {
		sc[i] = random() % (num_data / 5);
		idx[i] = idx2[i] = idx3[i] = i;
	}
	sc[num_data] = 10000;
	idx[num_data] = idx2[num_data] = idx3[num_data] = num_data;
	int nc = num_data / 10;

	printf("sort starts\n");
	clock_gettime(CLOCK_REALTIME, &tp1);
	quick_sort(idx2, sc, 0, num_data - 1);
	clock_gettime(CLOCK_REALTIME, &tp2);
	sec = tp2.tv_sec - tp1.tv_sec;
	nsec = tp2.tv_nsec - tp1.tv_nsec;
	if(nsec < 0){
		sec--;
		nsec += 1000000000L;
	}
	printf("sort ends: %ld.%09ld\n", sec, nsec);

	int k;
//	for(k = nc - 5; k < nc + 5; k++) {
//		printf("sc[idx2[%10d]] = %10d\n", k, sc[idx2[k]]);
//	}
	for(k = nc - 1; k < num_data && sc[idx2[k]] == sc[idx2[k + 1]]; k++);
	printf("k = %d, sc[idx2[%d] = %d\n", k, k + 1, sc[idx2[k + 1]]);

	printf("select starts\n");
	clock_gettime(CLOCK_REALTIME, &tp1);
//	quick_select_k(idx, sc, 0, num_data - 1, nc);
	#ifdef _OPENMP
	int nt = NUM_THREADS;
	printf("NUM_THREADS = %d\n", nt);
	omp_set_num_threads(nt);
	#endif
	int nc_by_para = quick_select_k_para(idx, sc, num_data, nc, nt);
	clock_gettime(CLOCK_REALTIME, &tp2);
	sec = tp2.tv_sec - tp1.tv_sec;
	nsec = tp2.tv_nsec - tp1.tv_nsec;
	if(nsec < 0){
		sec--;
		nsec += 1000000000L;
	}
	printf("select ends: %ld.%09ld: #selected = %d\n", sec, nsec, nc_by_para);
	for(int i = 0; i < nc_by_para; i++) {
		selected[i] = idx[i];
	}
	quick_sort(selected, sc, 0, nc_by_para - 1);

	printf("select work starts\n");
	work_select_para *wp = new_work_select_para(nt);
	clock_gettime(CLOCK_REALTIME, &tp1);
//	quick_select_k(idx, sc, 0, num_data - 1, nc);
	#ifdef _OPENMP
	printf("NUM_THREADS = %d\n", nt);
	omp_set_num_threads(nt);
	#endif
	int nc_by_para2 = quick_select_k_para_work(idx3, sc, num_data, nc, wp);
	clock_gettime(CLOCK_REALTIME, &tp2);
	sec = tp2.tv_sec - tp1.tv_sec;
	nsec = tp2.tv_nsec - tp1.tv_nsec;
	if(nsec < 0){
		sec--;
		nsec += 1000000000L;
	}
	printf("select work ends: %ld.%09ld: #selected = %d\n", sec, nsec, nc_by_para2);
	int k2 = get_top_k_from_work_para(selected2, idx3, wp); 
//	for(int i = 0; i < nc_by_para2; i++) {
//		selected2[i] = get_next_from_top_k(idx3, wp);
//	}
	quick_sort(selected2, sc, 0, nc_by_para - 1);

	for(int i = 0; i < nc_by_para; i++) {
		if(sc[selected[i]] != sc[idx2[i]]) {
			fprintf(stderr, "selected != idx2 at %d\n", i);
			break;
		}
	}

	for(int i = 0; i < nc_by_para; i++) {
		if(sc[selected[i]] != sc[selected2[i]]) {
			fprintf(stderr, "selected != selected2 at %d\n", i);
			break;
		}
	}
	return 0;
}
#endif

#ifndef TRIAL
#define TRIAL 5
#endif

// スケッチの配列を相対的にソートする．idxをソートする．
// quick sort のための 関数群 find_pivot, partition_by_pivot, quick_sort
// マルチスレッドによる並列処理のための関数群
int find_pivot(int idx[], unsigned sc[], int i, int j)

{
// idx[i] と idx[j] の比較 = sc[idx[i]] と sc[idx[j]] の比較
// 初めに TRIAL 回だけピボット候補 k をランダムに選ぶ．
// i 番目のデータ sc[idx[i]] と k（ただし，k = i + 1, ... , j）番目のデータ sc[idx[k]] を比較し，
// すべて等しいときには -1 を返し, 
// そうでないときには, sc[idx[i]] と異なる sc[idx[k]]  で最初に 
// 現れたもののうちで, 大きい方の位置(i または k) を返す. 
   int k;

	for(int t = 0; t < TRIAL; t++) {
		k = random() % (j - i + 1) + i;
		if(sc[idx[i]] != sc[idx[k]]) {
			return (sc[idx[i]] > sc[idx[k]] ? i : k);
		}
	}
	for(k = i + 1; k <= j; k++) {
		if(sc[idx[i]] != sc[idx[k]]) {
			return (sc[idx[i]] > sc[idx[k]] ? i : k);
		}
	}
	return -1;
}
/*
{
// idx[i] と idx[j] の比較 = sc[idx[i]] と sc[idx[j]] の比較
// i 番目のデータ sc[idx[i]] と k（ただし，k = i + 1, ... , j）番目のデータ sc[idx[k]] を比較し，
// すべて等しいときには -1 を返し, 
// そうでないときには, sc[idx[i]] と異なる sc[idx[k]]  で最初に 
// 現れたもののうちで, 大きい方の位置(i または k) を返す. 
   int k;

   for(k = i + 1; k <= j; k++)
      if(sc[idx[i]] != sc[idx[k]])
         return (sc[idx[i]] > sc[idx[k]] ? i : k);

   return -1;
}
*/

int partition_by_pivot(int idx[], unsigned sc[], int i, int j, unsigned piv)
// sc[idx[i]], ... , sc[idx[j]] を piv との大小によって分け， 
// piv より小さいものが sc[idx[i]], ... , sc[idx[k-1]] に，   
// そうでないものが sc[idx[k]], ... , sc[idx[j]] に来るようにする. 
// 右側のリストの先頭の位置(k)を返す. 
{
   int left, right;
   int temp;
   left = i;   right = j;
   do {
      while(sc[idx[left]] < piv) left++;
      while(sc[idx[right]] >= piv) right--;
      if (left < right) { 
         temp = idx[left];
         idx[left] = idx[right];
         idx[right] = temp;
      }
   } while(left <= right);
   return left;
}

void insertion_sort(unsigned data[], int n)
{
	int i, j;
	unsigned temp;
	for (i = 1; i < n; i++) {
		temp = data[i];
		for (j = i; j > 0 && data[j-1] > temp; j--) {
			data[j] = data[j - 1];
		}
		data[j] = temp;
	}
}

void quick_sort(int idx[], unsigned sc[], int i, int j)
{
	int pivotindex, k;
	pivotindex = find_pivot(idx, sc, i, j);
	if (pivotindex >= 0) {
		k = partition_by_pivot(idx, sc, i, j, sc[idx[pivotindex]]);
		quick_sort(idx, sc, i, k - 1);
		quick_sort(idx, sc, k, j);
	}
}

void quick_sort_para(int idx[], unsigned sc[], int i, int j)
{
	int pivotindex, k;
	pivotindex = find_pivot(idx, sc, i, j);
	if (pivotindex >= 0) {
		k = partition_by_pivot(idx, sc, i, j, sc[idx[pivotindex]]);
		#ifdef _OPENMP
		#pragma omp parallel sections
		#endif
		{
			#ifdef _OPENMP
			#pragma omp section
			#endif
			quick_sort_para(idx, sc, i, k - 1);
			#ifdef _OPENMP
			#pragma omp section
			#endif
			quick_sort_para(idx, sc, k, j);
		}
	}
}

void quick_select_k_r(int idx[], unsigned sc[], int i, int j, int k)
{
	int pivotindex, m;
	pivotindex = find_pivot(idx, sc, i, j);
	if (pivotindex >= 0) {
		m = partition_by_pivot(idx, sc, i, j, sc[idx[pivotindex]]);
		if(k < m) {
			quick_select_k_r(idx, sc, i, m - 1, k);
		} else if(k > m) {
			quick_select_k_r(idx, sc, m, j, k);
		}
	}
}

void quick_select_k(int idx[], unsigned sc[], int i, int j, int k)
{
	int pivotindex, m;
	while((pivotindex = find_pivot(idx, sc, i, j)) >= 0) {
		m = partition_by_pivot_para(idx, sc, i, j, sc[idx[pivotindex]]);
		if(k < m) {
			j = m - 1;
		} else if(k > m) {
			i = m;
		} else {
			break;
		}
	}
}

/*
int partition_by_pivot_bkt(int idx[], unsigned sc[], sk_num_pair sk_num[], int i, int j, unsigned piv, int *num)
// sc[idx[i]], ... , sc[idx[j]] を piv との大小によって分け， 
// piv より小さいものが sc[idx[i]], ... , sc[idx[k-1]] に，   
// そうでないものが sc[idx[k]], ... , sc[idx[j]] に来るようにする. 
// 右側のリストの先頭の位置(k)を返す. 
// *num には左側のリストのバケット要素数の合計を代入する．
{
	int left, right;
	int temp;
	int nn = 0;
	
	left = i;   right = j;
	do {
		while(sc[idx[left]] < piv) { nn += sk_num[idx[left]].num; left++; }
		while(sc[idx[right]] >= piv) right--;
		if (left < right) { 
			temp = idx[left];
			idx[left] = idx[right];
			idx[right] = temp;
		}
	} while(left <= right);
	*num = nn;
	return left;
}

int quick_select_bkt(int idx[], unsigned sc[], sk_num_pair sk_num[], int i, int j, int k)
// sc[i], ... , sc[j] のうちの上位が sc[idx[i]], ... , sc[idx[n]] となるように idx を入れ替えて，n + 1 を返す．
// ただし，sk_num[idx[i]].num + ... + sk_num[idx[n]] >= k となる最小の n にする．
{
	int pivotindex, m;
	int num = 0, num_temp; // 左に集めたスケッチの要素数の合計（これが k を超えるぎりぎりに左に集める）
	while((pivotindex = find_pivot(idx, sc, i, j)) >= 0) {
		m = partition_by_pivot_bkt(idx, sc, sk_num, i, j, sc[idx[pivotindex]], &num_temp);
		if(k < num + num_temp) {
			j = m - 1;
		} else if(k > num + num_temp) {
			num += num_temp;
			i = m;
		} else {
//			fprintf(stderr, "just k = %d found\n", num + num_temp);
			return m;
		}
	}
//	fprintf(stderr, "%d found in loop\n", num);
	while(num + sk_num[idx[i]].num < k) {
//		fprintf(stderr, "added %d, sc = %d\n", sk_num[idx[i]].num, sc[idx[i]]);
		num += sk_num[idx[i]].num;
		i++;
	}
//	fprintf(stderr, "found = %d\n", num + sk_num[idx[i]].num);
	return i + 1;
}
*/

int partition_by_pivot_para(int idx[], unsigned sc[], int i, int j, unsigned piv)
// 並列処理による quick_select_para のときは，空のパーティションができる可能性があるので，
// while の判定に添字範囲のチェック（left <= j および right >= i）を入れている
// sc[idx[i]], ... , sc[idx[j]] を piv との大小によって分け， 
// piv より小さいものが sc[idx[i]], ... , sc[idx[k-1]] に，   
// そうでないものが sc[idx[k]], ... , sc[idx[j]] に来るようにする. 
// 右側のリストの先頭の位置(k)を返す. 
{
   int left, right;
   int temp;
   left = i;   right = j;
   do {
      while(left <= j && sc[idx[left]] < piv) left++;
      while(right >= i && sc[idx[right]] >= piv) right--;
      if (left < right) { 
         temp = idx[left];
         idx[left] = idx[right];
         idx[right] = temp;
      }
   } while(left <= right);
   return left;
}

/*
0 = 先頭データの位置
n - 1 = 最終データの位置

n = データ数 (j - i + 1)
nt = スレッド数
s_t = スレッド t の分担範囲の要素数
s_t = (n + nt - t - 1) / nt = (n - t - 1) / nt + 1


ii[0] = i

L_i = L_{i-1} + s_{i-1} (i = 1, ... )
jj[t] = ii[t] + s_t - 1

*/

// n 個のデータ sc[0], ... , sc[n - 1] の上位 k 個を求める．
// ただし，sc そのものは並べ替えたりしないで，添字の配列 idx を並べ替えて相対的に求める．
// つまり，sc[idx[0]], sc[idx[1]], ... , sc[idx[k - 1]] が求める上位 k 個となるように idx を並べ替える
// ただし，上位 k 個は，値の順に並んでいる必要はない（保証はない）．
// 同点の k 位があるものも含めた個数（k 以上）を返す．
int quick_select_k_para(int idx[], unsigned sc[], int n, int k, int nt)
{
	int t, pivotindex[nt], ii[nt], jj[nt], m[nt], L[nt], R[nt];
	int target = k; // 求めている上位の個数
	int found = 0;		// 現状で見つかっている上位の個数

	// 分担範囲の設定
	// ii[t] = t の分担範囲の先頭位置
	// jj[t] = t の分担範囲の末尾位置
	
	ii[0] = 0;
	for(t = 0; t < nt; t++) {
		int s_t = (n - t - 1) / nt + 1; // スレッド t の分担範囲の要素数
		jj[t] = ii[t] + s_t - 1;
		L[t] = ii[t];
		R[t] = jj[t];
		m[t] = L[t];
		if(t < nt - 1) ii[t + 1] = ii[t] + s_t;
	}

	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif

	int pivot, pivot_L = INT_MAX;
	while(1) {
		int pi = -1;
		pivot = INT_MAX;
		#ifdef _OPENMP
		#pragma omp parallel for
		// 真の2分割が可能なピボット（その区間に2種類以上の値が含まれているときその2種類の内の大きい方）が存在すれば，その最小のものを求める．
		#endif
		for(t = 0; t < nt; t++) {
			if(L[t] >= R[t]) {
				pivotindex[t] = -1;
			} else {
				pivotindex[t] = find_pivot(idx, sc, L[t], R[t]);
				if(pivotindex[t] >= 0 && sc[idx[pivotindex[t]]] < pivot) {
					pi = pivotindex[t];
					pivot = sc[idx[pi]];
				}
			}
		}
		if(pi < 0) { // 真の2分割が確実なピボットが見つからなかった（すべてのスレッドの受け持ち区間が大きさ1以下かすべて同じ値だけになっている）
			// 空でない区間がもつ値の最小値を求めてそれをピボットにして分割をした状態にする．同じ値を持つ区間が複数ある可能性もあるので注意．
			// ピボット未満のものが targer 以上になったら終了．
			// すべての区間が空になってしまうことはないはず．（念のためすべて空になったら「エラー」終了するようにしておく）
			// 空でない区間がもつ異なる値が2個以上あれば，その大きい方を pivot にする．
//			fprintf(stderr, "no proper pivot found\n");
			do {
				pivot = INT_MAX;
				int pivot_min = INT_MAX;
				found = 0;
				for(t = 0; t < nt; t++) {
					if(L[t] > R[t]) {
						continue; // 区間が空
					}
					found += L[t] - ii[t];
					if(sc[idx[L[t]]] < pivot_min) {
						pivot = pivot_min;
						pivot_min = sc[idx[L[t]]];
					} else if(sc[idx[L[t]]] > pivot_min && sc[idx[L[t]]] < pivot) {
						pivot = sc[idx[L[t]]];
					}
				}
				if(pivot == INT_MAX) { // 空でない区間がもつ値がすべて同じだった．
					pivot = pivot_L;
					found = 0;
					for(t = 0; t < nt; t++) {
						if(L[t] <= R[t]) {
							L[t] = R[t] + 1;
						}
						found += L[t] - ii[t];
					}
					break;
				} else {
					found = 0;
					for(t = 0; t < nt; t++) {
						if(L[t] <= R[t] && sc[idx[L[t]]] < pivot) {
							L[t] = R[t] + 1;
						}
						found += L[t] - ii[t];
					}
				}
			} while(found < target);
			break; // while(1)
		}

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(t = 0; t < nt; t++) {
			m[t] = partition_by_pivot_para(idx, sc, L[t], R[t], pivot);
		}
		found = 0;
		for(t = 0; t < nt; t++) {
			found += m[t] - ii[t];
		}
		if(found > target) {
			pivot_L = pivot;
			for(t = 0; t < nt; t++) {
				R[t] = m[t] - 1;
			}
		} else if(found < target) {
			for(t = 0; t < nt; t++) {
				L[t] = m[t];
			}
		} else {
			pivot_L = pivot;
			for(t = 0; t < nt; t++) {
				L[t] = m[t];
			}
			break; // while(1)
		}
	}

// それぞれで見つけたもの（pivot 未満のもの）
//      0 : ii[0], ... , L[0] - 1
//      1 : ii[1], ... , L[1] - 1
// ...
// nt - 1 : ii[nt - 1], ... , L[nt - 1] - 1
//
// pivot 以上のものが入っている場所（詰め合わせで使える空いている場所）
//      0 : L[0], ... , jj[0]
//      t : L[t], ... , jj[t]
// ...
// nt - 1 :  

//	return 0;
	
	int lt = 0, lp = L[0], rt = nt - 1, rp = L[nt - 1] - 1, temp;
	while(lp < rp) {
		temp = idx[lp];
		idx[lp] = idx[rp];
		idx[rp] = temp;
		lp = (lp == jj[lt] ? L[++lt] : lp + 1);
		rp = (rp == ii[rt] ? L[--rt] - 1 : rp - 1);
	}

//	lp--;
#ifdef DEBUG
	num_top = 0; max_top = 0;
	max_top = 0;
	int num_less = 0, num_ge = 0;
	for(int x = 0; x < lp; x++) {
		if(sc[idx[x]] > max_top) max_top = sc[idx[x]];
		if(sc[idx[x]] < pivot) {
			num_less++;
		} else {
			num_ge++;
		}
	}
	min_rest = INT_MAX;
	int num_less2 = 0, num_ge2 = 0;
	for(int x = lp; x < n; x++) {
		if(sc[idx[x]] < min_rest) min_rest = sc[idx[x]];
		if(sc[idx[x]] < pivot) {
			num_less2++;
		} else {
			num_ge2++;
		}
	}
	printf("k = %d, lp = %d, max_top = %d, min_rest = %d, pivot = %d, pivot_L = %d\n", k, lp, max_top, min_rest, pivot, pivot_L);
	printf("num_less = %d, num_ge = %d, num_less2 = %d, num_ge2 = %d\n", num_less, num_ge, num_less2, num_ge2);
	for(int x = lp - 2; x < lp + 2; x++) {
		printf("sc[idx[%d]] = %d(%s), ", x, sc[idx[x]], sc[idx[x]] < pivot_L ? "<" : ">=");
	}
	printf("\n");
#endif
	return lp;
}

int partition_by_pivot_bkt_para(int idx[], unsigned sc[], sk_num_pair sk_num[], int i, int j, unsigned piv, int *num)
// 並列処理による quick_select_para のときは，空のパーティションができる可能性があるので，
// while の判定に添字範囲のチェック（left <= j および right >= i）を入れている
// sc[idx[i]], ... , sc[idx[j]] を piv との大小によって分け， 
// piv より小さいものが sc[idx[i]], ... , sc[idx[k-1]] に，   
// そうでないものが sc[idx[k]], ... , sc[idx[j]] に来るようにする. 
// 右側のリストの先頭の位置(k)を返す. 
{
	int left, right;
	int temp;
	int nn = 0;

	left = i;   right = j;
	do {
		while(left <= j && sc[idx[left]] < piv) { nn += sk_num[idx[left]].num; left++; }
		while(right >= i && sc[idx[right]] >= piv) right--;
		if (left < right) { 
			temp = idx[left];
			idx[left] = idx[right];
			idx[right] = temp;
		}
	} while(left <= right);
	*num = nn;
//	if(i == j) fprintf(stderr, "smallest partition: pivot = %d, i = %d, j = %d, num = %d\n", piv, i, j, nn); 
	return left;
}

//typedef struct {
//	int pivotindex, i, j, m, left, right, num;
//} work_select_bkt;

int quick_select_bkt_para(int idx[], unsigned sc[], sk_num_pair sk_num[], int n, int k, int nt)
// n 個のデータ sc[0], ... , sc[n - 1] の上位 k 個を求める．
// ただし，sc そのものは並べ替えたりしないで，添字の配列 idx を並べ替えて相対的に求める．
// つまり，sc[idx[0]], sc[idx[1]], ... , sc[idx[k - 1]] が求める上位 k 個となるように idx を並べ替える
// ただし，上位 k 個は，値の順に並んでいる必要はない（保証はない）．
// 同点の k 位があるものも含めた個数（k 以上）を返す．
{
	int t;
	work_select_bkt work[nt];
//	, pivotindex[nt], ii[nt], jj[nt], m[nt], L[nt], R[nt];
//	, num_temp[nt];
	int target = k; // 求めている上位の num の合計
	int found = 0;		// 現状で見つかっている上位の num の合計

	// 分担範囲の設定
	// work[t].i = t の分担範囲の先頭位置
	// work[t].j = t の分担範囲の末尾位置

	work[0].i = 0;
	//ii[0] = 0;
	for(t = 0; t < nt; t++) {
		int s_t = (n - t - 1) / nt + 1; // スレッド t の分担範囲の要素数
		work[t].j = work[t].i + s_t - 1;
		// jj[t] = ii[t] + s_t - 1;
		work[t].left = work[t].i;
		work[t].right = work[t].j;
		work[t].m = work[t].left;
		work[t].pivotindex = work[t].left;
		//L[t] = ii[t];
		//R[t] = jj[t];
		//m[t] = L[t];
		if(t < nt - 1) work[t + 1].i = work[t].i + s_t; // ii[t + 1] = ii[t] + s_t;
	}

	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif

	int pi, pivot, too_large = INT_MAX;
	int num = 0;
	while(1) {
//		#pragma omp parallel for
		pi = -1;
		// 2分割（未満も以上も空にならない）ができるピボットがあれば，それを採用する
		for(t = 0; t < nt; t++) {
			work_select_bkt *w = &work[t];
			if(w->pivotindex < 0) continue;
			if(w->left >= w->right) {
				w->pivotindex = -1;
			} else {
				w->pivotindex = find_pivot(idx, sc, w->left, w->right);
			}
			if(w->pivotindex >= 0) {
				pivot = sc[idx[w->pivotindex]];
				pi = w->pivotindex;
				break;
			}
		}
		if(pi < 0) { // 2分割が確実なピボットが見つからなかったので，
			pivot = INT_MAX;
			for(t = 0; t < nt; t++) {
				work_select_bkt *w = &work[t];
				int size = w->right - w->left + 1;
				if(size > 0) { // 分担区間が空でない
					if(pivot == INT_MAX) { // 最初はそれをピボットの候補にする
						pivot = sc[idx[w->left]];
						pi = w->left;
					} else if(pivot < sc[idx[w->left]]) { // 最初の候補より大きいものが見つかったので，それを採用（未満には最初の候補があり，以上も空にならない）
						pivot = sc[idx[w->left]];
						pi = w->left;
						break;
					} else if(pivot > sc[idx[w->left]]) { // 最初のものより小さい値が残っているので，最初のものを採用（未満が空にならない）
						break;
					}
				}
			}
			if(t == nt) pi = -1; // 2分割ができない（残っている値がすべて等しい）
		}
				
		if(pi < 0) { // 2分割できないので，最後の too_large で分割をしておく
			#ifdef _OPENMP
			#pragma omp parallel for
			#endif
			for(t = 0; t < nt; t++) {
				work[t].m = partition_by_pivot_bkt_para(idx, sc, sk_num, work[t].left, work[t].right, too_large, &work[t].num);
			}
			break;
		}
		// 2分割ができるピボットで分割を進める
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(t = 0; t < nt; t++) {
			work[t].m = partition_by_pivot_bkt_para(idx, sc, sk_num, work[t].left, work[t].right, pivot, &work[t].num);
		}
		found = 0; // 最後の分割で見つかったピボット未満のものの個数
		for(t = 0; t < nt; t++) {
			found += work[t].num;
		}
		if(num + found > target) { // ピボット未満が求めているものより多過ぎる => ピボット以上側は確定するので，未満側を調べる
			too_large = pivot;
			for(t = 0; t < nt; t++) {
				work[t].right = work[t].m - 1;
			}
		} else if(num + found < target) { // ピボット未満のものが求めているものより少ない．ピボット未満側は確定して，以上側を調べる
			num += found; // 確定した小さいもの個数
			for(t = 0; t < nt; t++) {
				work[t].left = work[t].m;
			}
		} else { // ピボット未満のものが「ちょうど」求めている個数になった．これで終了．
			break;
		}
	}

// それぞれで見つけたもの（pivot 未満のもの）
//      0 : ii[0], ... , m[0] - 1
//      1 : ii[1], ... , m[1] - 1
// ...
// nt - 1 : ii[nt - 1], ... , m[nt - 1] - 1
//
// pivot 以上のものが入っている場所（詰め合わせで使える空いている場所）
//      0 : m[0], ... , jj[0]
//      t : m[t], ... , jj[t]
// ...
// nt - 1 :  

	int num2 = 0, num_sk = 0;
	for(t = 0; t < nt; t++) {
		num_sk += work[t].m - work[t].i;
		for(int j = work[t].i; j < work[t].m; j++) {
			num2 += sk_num[idx[j]].num;
		}
	}
//	fprintf(stderr, "Before collect: sum = %d, num_sk = %d\n", num2, num_sk);
//	for(t = 0; t < nt; t++) {
//		work_select *w = &work[t];
//		fprintf(stderr, "t = %2d: i = %10d, j = %10d, L = %10d, R = %10d, m = %10d, found = %10d\n", t, w->i, w->j, w->left, w->right, w->m, w->m - w->i);
//	}
	
	int lt = 0, lp = work[0].m, rt = nt - 1, rp = work[nt - 1].m - 1, temp;

	int done = 0;
	for(lt = 0; lt < nt; lt++) {
		done += work[lt].m;
		if(work[lt].m <= work[lt].j) {
			lp = work[lt].m;
			break;
		}
	}
	if(lt >= nt) {
		fprintf(stderr, "(1) ERROR: no empty space to select\n"); // 空きがない
		exit(0);
	}
//	fprintf(stderr, "start collection: lt = %d, lp = %d\n", lt, lp);
	for(rt = nt - 1; rt >= 0; rt--) {
		if(work[rt].m > work[rt].i) {
			rp = work[rt].m - 1;
			break;
		}
	}
	if(rt < 0) {
		fprintf(stderr, "(1) ERROR: no small data\n"); // 空きがない
		exit(0);
	}
//	fprintf(stderr, "start collection: rt = %d, rp = %d\n", rt, rp);
//	int lt = 0, lp = m[0], rt = nt - 1, rp = m[nt - 1] - 1, temp;
	while(lp < rp && done < num_sk) {
		temp = idx[lp];
		idx[lp] = idx[rp];
		idx[rp] = temp;
		done++;
		if(done >= num_sk) break;
//		lp = (lp == work[lt].j ? work[++lt].m : lp + 1);
		if(lp < work[lt].j) {
			lp++;
		} else {
			lt++;
			while(lt < nt) {
				if(work[lt].m <= work[lt].j) {
					lp = work[lt].m;
					break;
				}
				lt++;
			}
			if(lt >= nt) {
				fprintf(stderr, "(2) ERROR: no empty space to select\n"); // 空きがない
				exit(0);
			}
//			fprintf(stderr, "lt changed to %d, lp = %d\n", lt, lp);
		}		
//		rp = (rp == work[rt].i ? work[--rt].m - 1 : rp - 1);
		if(rp > work[rt].i) {
			rp--;
		} else {
			rt--;
			while(rt >= 0) {
				if(work[rt].m > work[rt].i) {
					rp = work[rt].m - 1;
					break;
				}
				rt--;
			}
			if(rt < 0) {
				fprintf(stderr, "(2) ERROR: no small data, done = %d\n", done); // 空きがない
				exit(0);
			}
		}
	}

	return lp + 1;
}

work_select_para *new_work_select_para(int nt)
{
	work_select_para *ws = (work_select_para *)malloc(sizeof(work_select_para));
	ws->nt = nt;
	ws->work = (work_select *)malloc(sizeof(work_select) * nt);
	return ws;
}

int quick_select_k_para_work(int idx[], unsigned sc[], int n, int k, work_select_para *wp)
// n 個のデータ sc[0], ... , sc[n - 1] の上位 k 個を求める．
// ただし，sc そのものは並べ替えたりしないで，添字の配列 idx を並べ替えて相対的に求める．
// つまり，sc[idx[0]], sc[idx[1]], ... , sc[idx[k - 1]] が求める上位 k 個となるように idx を並べ替える
// ただし，上位 k 個は，値の順に並んでいる必要はない（保証はない）．
// 同点の k 位があるものも含めた個数（k 以上）を返す．
{
	int t;
	work_select *work = wp->work;
	int nt = wp->nt;
	int target = k; // 求めている上位の個数
	int found = 0;		// 現状で見つかっている上位の個数

	// 分担範囲の設定
	// ii[t] = t の分担範囲の先頭位置
	// jj[t] = t の分担範囲の末尾位置
	
	work[0].i = 0;
	for(t = 0; t < nt; t++) {
		work_select *w = &work[t];
		int s_t = (n - t - 1) / nt + 1; // スレッド t の分担範囲の要素数
		w->j = w->i + s_t - 1;
		w->left = w->i;
		w->right = w->j;
		w->m = w->left;
		if(t < nt - 1) work[t + 1].i = w->i + s_t;
	}

	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif

	int pivot, pivot_L = INT_MAX;
	while(1) {
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		// 分担範囲を真に2分割できるピボット（その区間に2種類以上の値が含まれているときその2種類の内の大きい方）を探す
		for(t = 0; t < nt; t++) {
			work_select *w = &work[t];
			if(w->left >= w->right) {
				w->pivotindex = -1;
			} else {
				w->pivotindex = find_pivot(idx, sc, w->left, w->right);
			}
		}
		// そうしたピボットが存在すれば，その最小のものを求める．
		int pi = -1;
		pivot = INT_MAX;
		for(t = 0; t < nt; t++) {
			work_select *w = &work[t];
			if(w->pivotindex >= 0 && sc[idx[w->pivotindex]] < pivot) {
				pi = w->pivotindex;
				pivot = sc[idx[pi]];
			}
		}
		if(pi < 0) { // 真の2分割が確実なピボットが見つからなかった（すべてのスレッドの受け持ち区間が大きさ1以下かすべて同じ値だけになっている）
			// 空でない区間がもつ値の最小値を求めてそれをピボットにして分割をした状態にする．同じ値を持つ区間が複数ある可能性もあるので注意．
			// ピボット未満のものが targer 以上になったら終了．
			// すべての区間が空になってしまうことはないはず．（念のためすべて空になったら「エラー」終了するようにしておく）
			// 空でない区間がもつ異なる値が2個以上あれば，その大きい方を pivot にする．
//			fprintf(stderr, "no proper pivot found\n");
			do {
				pivot = INT_MAX;
				int pivot_min = INT_MAX;
				found = 0;
				for(t = 0; t < nt; t++) {
					work_select *w = &work[t];
					if(w->left > w->right) {
						continue; // 区間が空
					}
					found += w->left - w->i;
					if(sc[idx[w->left]] < pivot_min) {
						pivot = pivot_min;
						pivot_min = sc[idx[w->left]];
					} else if(sc[idx[w->left]] > pivot_min && sc[idx[w->left]] < pivot) {
						pivot = sc[idx[w->left]];
					}
				}
				if(pivot == INT_MAX) { // 空でない区間がもつ値がすべて同じだった．
					pivot = pivot_L;
					found = 0;
					for(t = 0; t < nt; t++) {
						work_select *w = &work[t];
						if(w->left <= w->right) {
							w->left = w->right + 1;
						}
						found += w->left - w->i;
					}
					break;
				} else {
					found = 0;
					for(t = 0; t < nt; t++) {
						work_select *w = &work[t];
						if(w->left <= w->right && sc[idx[w->left]] < pivot) {
							w->left = w->right + 1;
						}
						found += w->left - w->i;
					}
				}
			} while(found < target);
			break; // while(1)
		}

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(t = 0; t < nt; t++) {
			work[t].m = partition_by_pivot_para(idx, sc, work[t].left, work[t].right, pivot);
		}
		found = 0;
		for(t = 0; t < nt; t++) {
			found += work[t].m - work[t].i;
		}
		if(found > target) {
			pivot_L = pivot;
			for(t = 0; t < nt; t++) {
				work[t].right = work[t].m - 1;
			}
		} else if(found < target) {
			for(t = 0; t < nt; t++) {
				work[t].left = work[t].m;
			}
		} else {
			pivot_L = pivot;
			for(t = 0; t < nt; t++) {
				work[t].left = work[t].m;
			}
			break; // while(1)
		}
	}
	wp->selected = found;
	wp->next_t = 0;
	wp->next = wp->work[0].i;
	
//	for(t = 0; t < nt; t++) {
//		work_select *w = &work[t];
//		printf("t, %d, i, %d, j, %d, m, %d, left, %d, right, %d\n", t, wp->work[t].i, wp->work[t].j, wp->work[t].m, wp->work[t].left, wp->work[t].right);
//	}
//	printf("found, %d\n", found);

	return found;
}

int get_top_k_from_work_para(int data_num[], int idx[], work_select_para *wp)
{
	int k = 0;
	while(1) {
		while(wp->next >= wp->work[wp->next_t].left) {
			if(wp->next_t == wp->nt - 1) return k;
			wp->next_t++;
			wp->next = wp->work[wp->next_t].i;
		}
//		printf("wp->next = %d, wp->next_t = %d, wp->work[wp->next_t].left = %d\n", wp->next, wp->next_t, wp->work[wp->next_t].left);
		data_num[k++] = idx[wp->next];
		wp->next++;
	}
	return k;
}
/*
int get_top_k_from_work_para(int data_num[], int idx[], work_select_para *wp)
{
	int k = 0;
	while(1) {
		while(wp->next >= wp->work[wp->next_t].left) {
			if(wp->next_t == wp->nt - 1) return k;
			wp->next_t++;
			wp->next = wp->work[wp->next_t].i;
		}
		data_num[k++] = idx[wp->next];
		wp->next++;
	}
	return k;
}
*/

// idx には，0, 1, ... , n - 1 を入れておく．
// n 個のデータ sc[0], ... , sc[n - 1] の上位 k 番目の値を求める．
// つまり，top k - 1 を求めるときの最後の sc を求める．
int get_threshold_k(int idx[], unsigned sc[], int n, int k, int nt)
{
	int t, pivotindex[nt], ii[nt], jj[nt], m[nt], L[nt], R[nt];
	int target = k; // 求めている上位の個数
	int found = 0;		// 現状で見つかっている上位の個数
/*	
	char chk[n];
	for(int i = 0; i < n; i++) chk[i] = 0;
	for(int i = 0; i < n; i++) chk[idx[i]]++;
	int ng = 0;
	for(int i = 0; i < n; i++) {
		if(chk[i] != 1) {
			printf("chk[%d] = %d\n", i, chk[i]);
			ng = 1;
		}
	}
	if(ng) exit(0);
*/
	// 分担範囲の設定
	// ii[t] = t の分担範囲の先頭位置
	// jj[t] = t の分担範囲の末尾位置
	
	ii[0] = 0;
	for(t = 0; t < nt; t++) {
		int s_t = (n - t - 1) / nt + 1; // スレッド t の分担範囲の要素数
		jj[t] = ii[t] + s_t - 1;
		L[t] = ii[t];
		R[t] = jj[t];
		m[t] = L[t];
		if(t < nt - 1) ii[t + 1] = ii[t] + s_t;
	}

	#ifdef _OPENMP
	omp_set_num_threads(nt);
	#endif

	int pivot, pivot_L = INT_MAX;
	while(1) {
		int pi = -1;
		pivot = INT_MAX;
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		// 真の2分割が可能なピボット（その区間に2種類以上の値が含まれているときその2種類の内の大きい方）が存在すれば，その最小のものを求める．
		for(t = 0; t < nt; t++) {
			if(L[t] >= R[t]) {
				pivotindex[t] = -1;
			} else {
				pivotindex[t] = find_pivot(idx, sc, L[t], R[t]);
//				if(pivotindex[t] >= 0 && sc[idx[pivotindex[t]]] < pivot) {
//					pi = pivotindex[t];
//					pivot = sc[idx[pi]];
//				}
			}
		}
		for(t = 0; t < nt; t++) {
			if(pivotindex[t] >= 0 && sc[idx[pivotindex[t]]] < pivot) {
				pi = pivotindex[t];
				pivot = sc[idx[pi]];
			}
		}
		if(pi < 0) { // 真の2分割が確実なピボットが見つからなかった（すべてのスレッドの受け持ち区間が大きさ1以下かすべて同じ値だけになっている）
			// 空でない区間がもつ値の最小値を求めてそれをピボットにして分割をした状態にする．同じ値を持つ区間が複数ある可能性もあるので注意．
			// ピボット未満のものが targer 以上になったら終了．
			// すべての区間が空になってしまうことはないはず．（念のためすべて空になったら「エラー」終了するようにしておく）
			// 空でない区間がもつ異なる値が2個以上あれば，その大きい方を pivot にする．
//			fprintf(stderr, "no proper pivot found\n");
			do {
				pivot = INT_MAX;
				int pivot_min = INT_MAX;
				found = 0;
				for(t = 0; t < nt; t++) {
					if(L[t] > R[t]) {
						continue; // 区間が空
					}
					found += L[t] - ii[t];
					if(sc[idx[L[t]]] < pivot_min) {
						pivot = pivot_min;
						pivot_min = sc[idx[L[t]]];
					} else if(sc[idx[L[t]]] > pivot_min && sc[idx[L[t]]] < pivot) {
						pivot = sc[idx[L[t]]];
					}
				}
				if(pivot == INT_MAX) { // 空でない区間がもつ値がすべて同じだった．
					pivot = pivot_L;
					found = 0;
					for(t = 0; t < nt; t++) {
						if(L[t] <= R[t]) {
							L[t] = R[t] + 1;
						}
						found += L[t] - ii[t];
					}
					break;
				} else {
					found = 0;
					for(t = 0; t < nt; t++) {
						if(L[t] <= R[t] && sc[idx[L[t]]] < pivot) {
							L[t] = R[t] + 1;
						}
						found += L[t] - ii[t];
					}
				}
			} while(found < target);
			break; // while(1)
		}

		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(t = 0; t < nt; t++) {
			m[t] = partition_by_pivot_para(idx, sc, L[t], R[t], pivot);
		}
		found = 0;
		for(t = 0; t < nt; t++) {
			found += m[t] - ii[t];
		}
		if(found > target) {
			pivot_L = pivot;
			for(t = 0; t < nt; t++) {
				R[t] = m[t] - 1;
			}
		} else if(found < target) {
			for(t = 0; t < nt; t++) {
				L[t] = m[t];
			}
		} else {
			pivot_L = pivot;
			for(t = 0; t < nt; t++) {
				L[t] = m[t];
			}
			break; // while(1)
		}
	}

	return pivot;
}
/*
#if defined(WIDE_SKETCH)
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

#elif defined(EXPANDED_SKETCH)

static int comp_sketch(sketch_type a, sketch_type b)
{
	for(int j = 0; j < SKETCH_SIZE; j++) {
		if(a[j] < b[j])
			return -1;
		else if(a[j] == b[j])
			continue;
		else
			return 1;
	}
	return 0;
}

// スケッチの配列を相対的にソートする．idxを入れ替える．
int find_pivot_for_sketch(int idx[], sketch_type sk[], int i, int j)
{
	int k, cmp;

	for(int t = 0; t < TRIAL; t++) {
		k = random() % (j - i + 1) + i;
		if((cmp = comp_sketch(sk[idx[i]], sk[idx[k]])) != 0) {
			return (cmp > 0 ? i : k);
		}
	}
	for(k = i + 1; k <= j; k++) {
		if((cmp = comp_sketch(sk[idx[i]], sk[idx[k]])) != 0) {
			return (cmp > 0 ? i : k);
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
		while(comp_sketch(sk[idx[left]], piv) < 0) left++;
		while(comp_sketch(sk[idx[right]], piv) >= 0) right--;
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
#endif
*/
