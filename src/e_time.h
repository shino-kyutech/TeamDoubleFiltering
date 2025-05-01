#pragma once
// 経過時間計測および使用メモリ表示
#include <sys/resource.h>
#include <time.h>
#include <stdlib.h>
#include <sys/types.h>
#include <unistd.h>

double e_time(struct timespec *start, struct timespec *end);

// VmSize or VmRSS
void use_system(char *rs);


// double allocated_memory = 0;

// #ifndef MEMORY_LIMIT
// #define MEMORY_LIMIT 48e9
// #endif

void *MALLOC(size_t size);
void FREE(void *m, size_t size);
