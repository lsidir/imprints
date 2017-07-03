/*
 * join.h
 *
 *  Created on: 4 Jun 2017
 *      Author: zeroxwj
 */

#ifndef JOIN_H_
#define JOIN_H_

#include "types.h"
#include "common.h"
#include "barrier.h"
#include "hashtable.h"

typedef struct {
	unsigned int threadId; /* id of each thread */
	unsigned int threadCount;
	//ParallelJoin* join;
	Column * colL;
	Column * colR;
	HashTable * ht;
	pthread_barrier_t * barrier;	/* for thread sync */

	unsigned int imps_bits;			/* # of bits used in hash value to encode imprints */
} JoinThreadArgs;


#define PREFETCH_OFFSET 16
#define MAX_IMPS_BITS 6


/**
 *
 */
query_result_t hashjoin(Column * colL, Column * colR, unsigned int threadCount, joinconf_t joinCfg);

//query_result_t imps_hashjoin(Column * colL, Column * colR, unsigned int threadCount);

#endif /* JOIN_H_ */
