typedef struct {
	int pivotindex, i, j, m, left, right;
} work_select;

typedef struct {
	int nt; 
	work_select *work;
	int selected;
	int next, next_t;
} work_select_para;

typedef struct {
	int pivotindex, i, j, m, left, right, num;
} work_select_bkt;

int find_pivot(int idx[], unsigned sc[], int i, int j);
int partition_by_pivot(int idx[], unsigned sc[], int i, int j, unsigned piv);
void insertion_sort(unsigned data[], int n);
void quick_sort(int idx[], unsigned sc[], int i, int j);
void quick_sort_para(int idx[], unsigned sc[], int i, int j);
void quick_select_k_r(int idx[], unsigned sc[], int i, int j, int k);
void quick_select_k(int idx[], unsigned sc[], int i, int j, int k);

// int partition_by_pivot_bkt(int idx[], unsigned sc[], sk_num_pair sk_num[], int i, int j, unsigned piv, int *num);
// int quick_select_bkt(int idx[], unsigned sc[], sk_num_pair sk_num[], int i, int j, int k);
int partition_by_pivot_para(int idx[], unsigned sc[], int i, int j, unsigned piv);
int quick_select_k_para(int idx[], unsigned sc[], int n, int k, int nt);

int partition_by_pivot_bkt_para(int idx[], unsigned sc[], sk_num_pair sk_num[], int i, int j, unsigned piv, int *num);
int quick_select_bkt_para(int idx[], unsigned sc[], sk_num_pair sk_num[], int n, int k, int nt);

work_select_para *new_work_select_para(int nt);
int quick_select_k_para_work(int idx[], unsigned sc[], int n, int k, work_select_para *wp);
int get_top_k_from_work_para(int data_num[], int idx[], work_select_para *wp);
int get_threshold_k(int idx[], unsigned sc[], int n, int k, int nt);

//#if defined(WIDE_SKETCH) || defined(EXPANDED_SKETCH)
//int find_pivot_for_sketch(int idx[], sketch_type sk[], int i, int j);
//int partition_by_pivot_for_sketch(int idx[], sketch_type sk[], int i, int j, sketch_type piv);
//void quick_sort_for_sketch(int idx[], sketch_type sk[], int i, int j);
//#endif
