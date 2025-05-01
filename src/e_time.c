// 経過時間計測および使用メモリ表示
#include <stdio.h>
#include <stdlib.h>
#include "e_time.h"

double e_time(struct timespec *start, struct timespec *end)
{
	long sec = end->tv_sec - start->tv_sec;
	long nsec = end->tv_nsec - start->tv_nsec;
	if(nsec < 0){
		sec--;
		nsec += 1000000000L;
	}
	return (double)sec + (double)nsec/1000000000;
}

size_t allocated_memory = 0;
// double allocated_memory = 0;

// VmSize or VmRSS
void use_system(char *rs)
{
	char command[500];

	fflush(stdout);
	sprintf(command, "grep %s /proc/%d/status", rs, getpid());
	if(system(command) == -1) {
		fprintf(stderr, "system() error.\n");
		exit(1);
	}
	printf("allocated memory = %10.4e\n", (double)allocated_memory);
	fflush(stdout);
}


#ifndef MEMORY_LIMIT
#define MEMORY_LIMIT 48e9
#endif

void *MALLOC(size_t size) {
    if(allocated_memory + size > MEMORY_LIMIT) {
        fprintf(stderr, "Cannot allocate memory: current allocated = %10.4e, requested = %10.4e, LIMIT = %10.4e\n", (double)allocated_memory, (double)size, MEMORY_LIMIT);
        exit(0);
    }
    allocated_memory += size;
    return malloc(size);
}

void FREE(void *m, size_t size) {
    free(m);
    allocated_memory -= size;
}

