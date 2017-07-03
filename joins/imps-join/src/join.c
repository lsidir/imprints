/*
 * join.c
 *
 *  Created on: 5 Jun 2017
 *      Author: zeroxwj
 */

#include "join.h"
#include "utils.h"

const uint64_t oid_nil = 1ULL << (sizeof(uint64_t) * 8 - 1);

#define INSERT_POS(H, P, O) 														\
do {																				\
	uint64_t pos = P;																\
    while (ht->table[pos].hash														\
            || (!__sync_bool_compare_and_swap(&ht->table[pos].hash, 0, H))) {		\
        pos = (pos + 1) & ht->mask;													\
    }																				\
    ht->table[pos].oid = O;															\
} while (0)

#define HASH(KEY, CACHE)						\
do {											\
        CACHE.hash = hashKey(KEY);				\
        CACHE.pos = CACHE.hash & ht->mask;		\
} while (0)

#define HASH_IMPS(KEY, CACHE, LMASK, HMASK)						\
do {											\
        CACHE.hash = hashKey_imps(KEY, LMASK, HMASK);				\
        CACHE.pos = CACHE.hash & ht->mask;		\
} while (0)

#define LOOKUP_POS(P, K) 														\
do {																			\
	uint64_t pos = P;															\
    while ((ht->table[pos].hash) && (colRValues[ht->table[pos].oid] != K)) {	\
        pos = (pos + 1) & ht->mask;												\
    }																			\
																				\
    found_tuple_oid = ht->table[pos].hash ? ht->table[pos].oid : oid_nil;		\
} while (0)

#define GETBIN(Z,V,B,S,X)				\
do {								\
	int _i;							\
	Z = 0;							\
	for (_i = S; _i < B; _i += S)	\
		Z += ((V) >= bounds[_i].X);	\
} while (0)

void* run_multipass(void *arg) {

	JoinThreadArgs* threadArgs = (JoinThreadArgs*)(arg);

	//affinitizer.affinitize(threadArgs->threadNum);
	//threadArgs->join->runInParallel(threadArgs->threadNum);
	unsigned int threadId = threadArgs->threadId;
	unsigned int threadCount = threadArgs->threadCount;
	Column* colL = threadArgs->colL;
	Column* colR = threadArgs->colR;
	HashTable * ht = threadArgs->ht;
	pthread_barrier_t * barrier = threadArgs->barrier;

	unsigned int imps_bits = threadArgs->imps_bits;

	unsigned long colLCount = colL->colcount;
	unsigned long colRCount = colR->colcount;

	int rv;

	/* initialize hash table */
	uint64_t htSize = ht->nBuckets * sizeof(HTEntry);
	if (threadId == 0) {
		VERBOSE fprintf(stdout, "HT: initializing memory in thread 0\n");
		ht->table = (HTEntry *) malloc_aligned(htSize);
		memset(ht->table, 0, htSize);
	}
	BARRIER_ARRIVE(barrier, rv);

#define BUILD_MULTIPASS(T, X) {																				\
		const unsigned int bin_merge_bits = MAX_IMPS_BITS - imps_bits;										\
		const unsigned int step = (unsigned int)1 << bin_merge_bits;										\
        const uint64_t per_thread = colRCount / threadCount;												\
        const uint64_t startIndex = threadId * per_thread;													\
        const uint64_t endIndex = ((threadId == threadCount - 1) ? colRCount : (threadId + 1) * per_thread);\
        																									\
        T  *restrict colRValues = (T *) colR->col;															\
        																									\
		Imprints_index * r_imps_idx = colR->imps_idx;														\
		const ValRecord *restrict bounds = r_imps_idx->bounds;												\
		const int B = r_imps_idx->bins;																		\
		uint64_t lowermask = (ht->mask) >> imps_bits;														\
		unsigned int left_shift = __builtin_popcount(ht->mask) - imps_bits;									\
		uint64_t highermask = 0; 																			\
		unsigned long bin = 0; 																				\
																											\
        PrefetchCache *cache = (PrefetchCache *) malloc_aligned(PREFETCH_OFFSET * sizeof(PrefetchCache));	\
        const uint64_t cacheMask = PREFETCH_OFFSET - 1;														\
        																									\
        /* fill cache*/																						\
        for (uint32_t i = startIndex; i < startIndex + PREFETCH_OFFSET; i++) {								\
            const uint64_t cacheIndex = i & cacheMask;														\
            																								\
            GETBIN(bin, colRValues[i], B, step, X);															\
			highermask = bin << left_shift;																	\
            HASH_IMPS(colRValues[i], cache[cacheIndex], lowermask, highermask);								\
            __builtin_prefetch(&ht->table[cache[cacheIndex].pos], 1, 1);									\
        }																									\
																											\
        /* build - [startIndex ... endIndex - prefetchOffset)	*/											\
        uint64_t prefetchIndex = startIndex + PREFETCH_OFFSET;    											\
        const uint64_t endIndexPrefetch = endIndex - PREFETCH_OFFSET;										\
        for (uint64_t i = startIndex; i < endIndexPrefetch; i++) {											\
            /* re-use hash values and positions, which have been already computed*/							\
            const uint64_t cacheIndex = i & cacheMask;														\
            INSERT_POS(cache[cacheIndex].hash, cache[cacheIndex].pos, i);									\
            																								\
            /* prefetching*/																				\
			GETBIN(bin, colRValues[prefetchIndex], B, step, X);												\
			highermask = bin << left_shift;																	\
            HASH_IMPS(colRValues[prefetchIndex], cache[cacheIndex], lowermask, highermask);					\
            __builtin_prefetch(&ht->table[cache[cacheIndex].pos], 1, 1);									\
            prefetchIndex++;																				\
        }																									\
        for (uint64_t i = endIndexPrefetch; i < endIndex; i++) {											\
            const uint64_t cacheIndex = i & cacheMask;														\
            INSERT_POS(cache[cacheIndex].hash, cache[cacheIndex].pos, i);									\
        }																									\
        free(cache);																						\
}

    	struct timeval start_timeval, build_timeval;

        if (threadId == 0) {
        	gettimeofday(&start_timeval, NULL);
        }

    	/* build phase ------------------------------------------------------------------*/
    	switch (colR->coltype) {
    		case TYPE_sht: BUILD_MULTIPASS(short, sval); break;
    		case TYPE_int: BUILD_MULTIPASS(int, ival); break;
    		case TYPE_lng: BUILD_MULTIPASS(long, lval); break;
    		default: break;
    	}

        //syncThreads();
        BARRIER_ARRIVE(barrier, rv);

        if (threadId == 0) {
        	gettimeofday(&build_timeval, NULL);
        }

#define PROBE_MULTIPASS(T) {																						\
    	T  *restrict colLValues = (T *) colL->col;																	\
    	T  *restrict colRValues = (T *) colR->col;																	\
        const uint64_t per_thread = colLCount / threadCount;														\
        const uint64_t startIndex = threadId * per_thread;															\
        const uint64_t endIndex = ((threadId == threadCount - 1) ? colLCount : (threadId + 1) * per_thread);		\
        /*uint64_t count = 0;*/																						\
        /*uint64_t preJoinCount = 0;			*/																	\
        /*double checksum = 0;					*/																	\
																											\
		/*Imprints_index * l_imps_idx = colL->imps_idx;	*/													\
		/*const ValRecord *restrict bounds = l_imps_idx->bounds;*/												\
		/*const int B = l_imps_idx->bins;					*/													\
		/*uint64_t lowermask = (ht->mask) >> imps_bits;		*/												\
		/*unsigned int left_shift = __builtin_popcount(ht->mask) - imps_bits;		*/							\
		/*uint64_t highermask = 0; 												*/							\
		/*unsigned long bin = 0; 												*/								\
        																											\
        PrefetchCache *cache = (PrefetchCache *) malloc_aligned(PREFETCH_OFFSET * sizeof(PrefetchCache));			\
        const uint64_t cacheMask = (1 << unsignedLog2(PREFETCH_OFFSET)) - 1;										\
        																											\
		/* fill cache*/																								\
		for (uint32_t i = startIndex; i < startIndex + PREFETCH_OFFSET; i++) {										\
			const uint64_t cacheIndex = i & cacheMask;																\
			HASH(colLValues[i], cache[cacheIndex]);																	\
			__builtin_prefetch(&ht->table[cache[cacheIndex].pos], 0, 1);											\
		}																											\
																													\
		/* probe [0 .. n - prefetchOffset)*/																		\
		uint64_t prefetchIndex = startIndex + PREFETCH_OFFSET; 														\
		const uint64_t endIndexPrefetch = endIndex - PREFETCH_OFFSET;												\
		for (uint64_t i = startIndex; i < endIndexPrefetch; i++) {													\
			/* re-use hash values and position, which have been already computed during prefetching */				\
			const uint64_t cacheIndex = i & cacheMask;																\
			uint64_t found_tuple_oid = oid_nil;																		\
			LOOKUP_POS(cache[cacheIndex].pos, colLValues[i]);														\
			count += (found_tuple_oid != oid_nil);																	\
																													\
			/* prefetching	*/																						\
			HASH(colLValues[prefetchIndex], cache[cacheIndex]);														\
			__builtin_prefetch(&ht->table[cache[cacheIndex].pos], 0, 1);											\
			prefetchIndex++;																						\
		}																											\
		for (uint64_t i = endIndexPrefetch; i < endIndex; i++) {													\
			const uint64_t cacheIndex = i & cacheMask;																\
			uint64_t found_tuple_oid = oid_nil;																		\
			LOOKUP_POS(cache[cacheIndex].pos, colLValues[i]);														\
			count += (found_tuple_oid != oid_nil);																	\
		}																											\
		free(cache);																								\
    }


    BARRIER_ARRIVE(barrier, rv);

    uint64_t count = 0;
    // probe phase ------------------------------------------------------------------
	switch (colL->coltype) {
		case TYPE_sht: PROBE_MULTIPASS(short); break;
		case TYPE_int: PROBE_MULTIPASS(int); break;
		case TYPE_lng: PROBE_MULTIPASS(long); break;
		default: break;
	}


	BARRIER_ARRIVE(barrier, rv);

    // print report ------------------------------------------------------------------
    {
		fprintf(stdout, "count:%ld\n", count);
    }

    free(threadArgs);

	if (threadId == 0) {
		free(ht->table);
	}
	return NULL;
}

void* run_default(void *arg) {

	JoinThreadArgs* threadArgs = (JoinThreadArgs*)(arg);

	//affinitizer.affinitize(threadArgs->threadNum);
	//threadArgs->join->runInParallel(threadArgs->threadNum);
	unsigned int threadId = threadArgs->threadId;
	unsigned int threadCount = threadArgs->threadCount;
	Column* colL = threadArgs->colL;
	Column* colR = threadArgs->colR;
	HashTable * ht = threadArgs->ht;
	pthread_barrier_t * barrier = threadArgs->barrier;

	unsigned long colLCount = colL->colcount;
	unsigned long colRCount = colR->colcount;

	int rv;

	/* initialize hash table */
	uint64_t htSize = ht->nBuckets * sizeof(HTEntry);
	if (threadId == 0) {
		VERBOSE fprintf(stdout, "HT: initializing memory in thread 0\n");
		ht->table = (HTEntry *) malloc_aligned(htSize);
		memset(ht->table, 0, htSize);
	}
	BARRIER_ARRIVE(barrier, rv);

#define BUILD(T) {																							\
        const uint64_t per_thread = colRCount / threadCount;												\
        const uint64_t startIndex = threadId * per_thread;													\
        const uint64_t endIndex = ((threadId == threadCount - 1) ? colRCount : (threadId + 1) * per_thread);\
        																									\
        T  *restrict colRValues = (T *) colR->col;															\
																											\
        PrefetchCache *cache = (PrefetchCache *) malloc_aligned(PREFETCH_OFFSET * sizeof(PrefetchCache));	\
        const uint64_t cacheMask = PREFETCH_OFFSET - 1;														\
        																									\
        /* fill cache*/																						\
        for (uint32_t i = startIndex; i < startIndex + PREFETCH_OFFSET; i++) {								\
            const uint64_t cacheIndex = i & cacheMask;														\
            HASH(colRValues[i], cache[cacheIndex]);															\
            __builtin_prefetch(&ht->table[cache[cacheIndex].pos], 1, 1);									\
        }																									\
																											\
        /* build - [startIndex ... endIndex - prefetchOffset)	*/											\
        uint64_t prefetchIndex = startIndex + PREFETCH_OFFSET;    											\
        const uint64_t endIndexPrefetch = endIndex - PREFETCH_OFFSET;										\
        for (uint64_t i = startIndex; i < endIndexPrefetch; i++) {											\
            /* re-use hash values and positions, which have been already computed*/							\
            const uint64_t cacheIndex = i & cacheMask;														\
            INSERT_POS(cache[cacheIndex].hash, cache[cacheIndex].pos, i);									\
            																								\
            /* prefetching*/																				\
            HASH(colRValues[prefetchIndex], cache[cacheIndex]);												\
            __builtin_prefetch(&ht->table[cache[cacheIndex].pos], 1, 1);									\
            prefetchIndex++;																				\
        }																									\
        for (uint64_t i = endIndexPrefetch; i < endIndex; i++) {											\
            const uint64_t cacheIndex = i & cacheMask;														\
            INSERT_POS(cache[cacheIndex].hash, cache[cacheIndex].pos, i);									\
        }																									\
        free(cache);																						\
}

    	struct timeval start_timeval, build_timeval;

        if (threadId == 0) {
        	gettimeofday(&start_timeval, NULL);
        }

    	/* build phase ------------------------------------------------------------------*/
    	switch (colR->coltype) {
    		case TYPE_sht: BUILD(short); break;
    		case TYPE_int: BUILD(int); break;
    		case TYPE_lng: BUILD(long); break;
    		default: break;
    	}

        //syncThreads();
        BARRIER_ARRIVE(barrier, rv);

        if (threadId == 0) {
        	gettimeofday(&build_timeval, NULL);
        }

#define PROBE(T) {																									\
    	T  *restrict colLValues = (T *) colL->col;																	\
    	T  *restrict colRValues = (T *) colR->col;																	\
        const uint64_t per_thread = colLCount / threadCount;														\
        const uint64_t startIndex = threadId * per_thread;															\
        const uint64_t endIndex = ((threadId == threadCount - 1) ? colLCount : (threadId + 1) * per_thread);		\
        /*uint64_t count = 0;*/																						\
        /*uint64_t preJoinCount = 0;			*/																	\
        /*double checksum = 0;					*/																	\
        																											\
        PrefetchCache *cache = (PrefetchCache *) malloc_aligned(PREFETCH_OFFSET * sizeof(PrefetchCache));			\
        const uint64_t cacheMask = (1 << unsignedLog2(PREFETCH_OFFSET)) - 1;										\
        																											\
		/* fill cache*/																								\
		for (uint32_t i = startIndex; i < startIndex + PREFETCH_OFFSET; i++) {										\
			const uint64_t cacheIndex = i & cacheMask;																\
			HASH(colLValues[i], cache[cacheIndex]);																	\
			__builtin_prefetch(&ht->table[cache[cacheIndex].pos], 0, 1);											\
		}																											\
																													\
		/* probe [0 .. n - prefetchOffset)*/																		\
		uint64_t prefetchIndex = startIndex + PREFETCH_OFFSET; 														\
		const uint64_t endIndexPrefetch = endIndex - PREFETCH_OFFSET;												\
		for (uint64_t i = startIndex; i < endIndexPrefetch; i++) {													\
			/* re-use hash values and position, which have been already computed during prefetching */				\
			const uint64_t cacheIndex = i & cacheMask;																\
			uint64_t found_tuple_oid = oid_nil;																		\
			LOOKUP_POS(cache[cacheIndex].pos, colLValues[i]);														\
			count += (found_tuple_oid != oid_nil);																	\
																													\
			/* prefetching	*/																						\
			HASH(colLValues[prefetchIndex], cache[cacheIndex]);														\
			__builtin_prefetch(&ht->table[cache[cacheIndex].pos], 0, 1);											\
			prefetchIndex++;																						\
		}																											\
		for (uint64_t i = endIndexPrefetch; i < endIndex; i++) {													\
			const uint64_t cacheIndex = i & cacheMask;																\
			uint64_t found_tuple_oid = oid_nil;																		\
			LOOKUP_POS(cache[cacheIndex].pos, colLValues[i]);														\
			count += (found_tuple_oid != oid_nil);																	\
		}																											\
		free(cache);																								\
    }


    BARRIER_ARRIVE(barrier, rv);

    uint64_t count = 0;
    // probe phase ------------------------------------------------------------------
	switch (colL->coltype) {
		case TYPE_sht: PROBE(short); break;
		case TYPE_int: PROBE(int); break;
		case TYPE_lng: PROBE(long); break;
		default: break;
	}


	BARRIER_ARRIVE(barrier, rv);

    // print report ------------------------------------------------------------------
    {
		fprintf(stdout, "count:%ld\n", count);
    }

    free(threadArgs);

	if (threadId == 0) {
		free(ht->table);
	}
	return NULL;
}

query_result_t hashjoin(Column * colL, Column * colR, unsigned int threadCount, joinconf_t joinCfg) {
	query_result_t result;

	/* step 1: initialize the hash table */
	unsigned long colRCount = colR->colcount;
    uint64_t htSize = nextPowerOfTwo(colRCount) << 1;
    uint64_t htBits = unsignedLog2(htSize);

    HashTable * ht = (HashTable *) malloc_aligned(sizeof(HashTable));
    //initHashTable(ht, htBits, threadCount);
	ht->nBuckets = 1ull << htBits;
	ht->mask = ht->nBuckets - 1;

    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, threadCount);

    /* step 2: multi-threaded build and probe */
	pthread_t threads[1000];

	for (unsigned i = 0; i < threadCount; i++) {
		JoinThreadArgs* threadArgs = (JoinThreadArgs *)malloc_aligned(sizeof(JoinThreadArgs));
		//threadArgs->join = parallelJoin;
		threadArgs->threadId = i;
		threadArgs->threadCount = threadCount;
		threadArgs->colL = colL;
		threadArgs->colR = colR;
		threadArgs->ht = ht;
		threadArgs->barrier = &barrier;
		threadArgs->imps_bits = joinCfg.imps_bits;
		switch(joinCfg.joinalgo) {
		case 0:	//defaults
			pthread_create(&threads[i], NULL, run_default, (void*)threadArgs);
			break;
		case 1:	//multi-pass
			pthread_create(&threads[i], NULL, run_multipass, (void*)threadArgs);
			break;
		default:
			break;
		}
	}
	for (unsigned i = 0; i < threadCount; i++) {
		pthread_join(threads[i], NULL);
	}
	//query_result_t result = parallelJoin->getResult();

	free(ht);

	return result;
}

#if 0
query_result_t hashjoin_visitlist(Column * colL, Column * colR, unsigned int threadCount) {

	query_result_t result;

	/* step 1: build the hash table */
	unsigned long colRCount = colR->colcount;
    uint64_t htSize = nextPowerOfTwo(colRCount) << 1;
    uint64_t htBits = unsignedLog2(htSize);

    HashTable * ht = (HashTable *) malloc_aligned(sizeof(HashTable));
    createHashTable(ht, htBits, threadCount);


	return result;
}
#endif
