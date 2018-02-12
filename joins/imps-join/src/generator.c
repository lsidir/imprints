

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#endif
#include <stdio.h>              /* perror */
#include <stdlib.h>             /* posix_memalign */
#include <math.h>               /* fmod, pow */
#include <time.h>               /* time() */
#include <unistd.h>             /* getpagesize() */
#include <string.h>             /* memcpy() */

#include "generator.h"          /* create_relation_*() */

/* return a random number in range [0,N] */
#define RAND_RANGE(N) ((double)rand() / ((double)RAND_MAX + 1) * (N))
#define RAND_RANGE48(N,STATE) ((double)nrand48(STATE)/((double)RAND_MAX+1)*(N))
#define MALLOC(SZ) alloc_aligned(SZ+RELATION_PADDING) /*malloc(SZ+RELATION_PADDING)*/ 
#define FREE(X,SZ) free(X)

/**
 * Generate tuple IDs -> random distribution
 * column must have been allocated
 */
#define RANDOM_GEN(column, maxid, T)								\
	do {											\
		uint32_t _i;										\
		T  *restrict C = (T *) column->col;				\
		for (_i = 0; _i < column->colcount; _i++)		\
			C[_i] = RAND_RANGE(maxid);		\
	} while (0)

#define RANDOM_UNIQUE_GEN(column, T)	\
	do {								\
		uint32_t i;						\
		T  *restrict C = (T *) column->col;	\
	    for (i = 0; i < column->colcount; i++) 	\
	    	C[i] = (i+1);						\
	    /* randomly shuffle elements */			\
	    KNUTH_SHUFFLE(C, column->colcount, T);						\
    } while (0)

/**
 * Shuffle column values using Knuth shuffle.
 */
#define KNUTH_SHUFFLE(col, count, T)		\
	do {							\
		int32_t i, j;					\
		T tmp;						\
		for (i = count - 1; i > 0; i--) {						\
			j = RAND_RANGE(i);						\
			tmp = col[i];							\
			col[i] = col[j];						\
			col[j] = tmp;							\
		}											\
	} while (0)

#define GEN_FK_FROM_PK(fkcol, pkcol, T)				\
	do {											\
		int i, iters;								\
		int64_t remainder;							\
		T  *restrict C = (T *) fkcol->col;					\
													\
		iters = fkcol->colcount / pkcol->colcount;	\
		for(i = 0; i < iters; i++){					\
			memcpy(C + i * pkcol->colcount, 	\
					pkcol->col,								\
					pkcol->colcount * sizeof(T));			\
		}													\
		/* if fkcol->colcount is not an exact multiple of pkcol->colcount*/	\
		remainder = fkcol->colcount % pkcol->colcount;						\
		if(remainder > 0) {													\
			memcpy(C + i * pkcol->colcount, 					\
					pkcol->col,												\
					remainder * sizeof(T));									\
		}																	\
		KNUTH_SHUFFLE(C, fkcol->colcount, T);												\
    } while(0)

/* Uncomment the following to persist input relations to disk. */
/* #define PERSIST_RELATIONS 1 */

/** An experimental feature to allocate input relations numa-local */
int numalocalize;
int nthreads;

static int seeded = 0;
static unsigned int seedValue;

void 
seed_generator(unsigned int seed) 
{
    srand(seed);
    seedValue = seed;
    seeded = 1;
}

/** Check wheter seeded, if not seed the generator with current time */
static void
check_seed()
{
    if(!seeded) {
        seedValue = time(NULL);
        srand(seedValue);
        seeded = 1;
    }
}

int 
create_relation_pk(Column * column)
{
    check_seed();

	switch (column->coltype) {
		case TYPE_sht: RANDOM_UNIQUE_GEN(column, short); break;
		case TYPE_int: RANDOM_UNIQUE_GEN(column, int); break;
		case TYPE_lng: RANDOM_UNIQUE_GEN(column, long); break;
		default: break;
	}

#ifdef PERSIST_RELATIONS
    write_relation(relation, "R.tbl");
#endif

    return 0;
}

/** 
 * Create a foreign-key column using the given primary-key column and
 * foreign-key column size. Keys in pkcol is randomly distributed in the full
 * integer range.
 * 
 * @param fkcol [output] foreign-key column
 * @param pkcol [input] primary-key column
 * 
 * @return 
 */
int 
create_column_fk_from_pk(Column *fkcol, Column *pkcol)
{
	switch (fkcol->coltype) {
		case TYPE_sht: GEN_FK_FROM_PK(fkcol, pkcol, short); break;
		case TYPE_int: GEN_FK_FROM_PK(fkcol, pkcol, int); break;
		case TYPE_lng: GEN_FK_FROM_PK(fkcol, pkcol, long); break;
		default: break;
	}

    return 0;
}


int create_column_nonunique(Column *column, const int32_t maxid)
{
    check_seed();

	switch (column->coltype) {
		case TYPE_sht: RANDOM_GEN(column, maxid, short); break;
		case TYPE_int: RANDOM_GEN(column, maxid, int); break;
		case TYPE_lng: RANDOM_GEN(column, maxid, long); break;
		default: break;
	}

    return 0;
}

#if 0
int
ReadColumnFromFile(Column *column, char * filename)
{
	int  *restrict column_values = (int *) column->col;	// all regard as integers for now
	FILE* file = fopen(filename, "r");
	if (file == NULL) {
		printf("[ERROR] Failed to open column file %s\n", filename);
		exit(EXIT_SUCCESS);
	}
	char line[256];
	unsigned long count = 0;

	while (fgets(line, sizeof(line), file)) {
		column_values[count] = atoi(line);
		count++;
	}
	assert(count == column->colcount);

	fclose(file);

	return 0;
}
#endif

int
ReadColumnFromFile(Column *column, char * filename)
{
	long  *restrict column_values = (long *) column->col;	// all regard as integers for now
	FILE* file = fopen(filename, "r");
	if (file == NULL) {
		printf("[ERROR] Failed to open column file %s\n", filename);
		exit(EXIT_SUCCESS);
	}
	char line[256];
	unsigned long count = 0;

	while (fgets(line, sizeof(line), file)) {
		column_values[count] = atol(line);
		count++;
	}
	assert(count == column->colcount);

	fclose(file);

	return 0;
}

