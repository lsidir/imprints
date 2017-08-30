/*
 * join.c
 *
 *  Created on: 5 Jun 2017
 *      Author: zeroxwj
 */

#include "join.h"
#include "utils.h"
#include "hybrid_timer.h"
#include "vector.h"

const uint64_t oid_nil = 1ULL << (sizeof(uint64_t) * 8 - 1);
const uint32_t hash_nil = 0xFFFFFFFF;

#define INSERT_POS(H, P, O) 														\
do {																				\
	uint64_t pos = P;																\
    while ((ht->table[pos].hash != hash_nil)														\
            || (!__sync_bool_compare_and_swap(&ht->table[pos].hash, hash_nil, H))) {		\
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
    while ((ht->table[pos].hash != hash_nil) && (colRValues[ht->table[pos].oid] != K)) {	\
        pos = (pos + 1) & ht->mask;												\
    }																			\
																				\
    found_tuple_oid = (ht->table[pos].hash != hash_nil) ? ht->table[pos].oid : oid_nil;		\
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

    ThreadDctPointer* dctPtrs = threadArgs->dctPointers;
    uint64_t* matches = threadArgs->matches;
    unsigned int num_cacheline_per_thread = threadArgs->num_cacheline_per_thread;

	unsigned long colLCount = colL->colcount;
	unsigned long colRCount = colR->colcount;

	int rv;

	/* initialize hash table */
	uint64_t htSize = ht->nBuckets * sizeof(HTEntry);
	if (threadId == 0) {
		VERBOSE fprintf(stdout, "HT: initializing memory in thread 0\n");
		ht->table = (HTEntry *) malloc_aligned(htSize);
		memset(ht->table, 0xFF, htSize);
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

    	//struct timeval start_timeval, build_timeval;
		HybridTimer build_timer;
		HybridTimer probe_timer;

        if (threadId == 0) {
        	//gettimeofday(&start_timeval, NULL);
        	timer_start(&build_timer);
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
        	//gettimeofday(&build_timeval, NULL);
        	timer_stop(&build_timer);
        	printf("[INFO] BUILD phase takes %lf ms\n", GetMilliSeconds(&build_timer));
        	timer_start(&probe_timer);
        }
#if 0
		int i;
		for (i = 0; i < 64; ++i) {
        	printf("%d\n", bounds[i].ival);
        }
		int * value_checked = (int *)malloc_aligned((endIndex - startIndex) * sizeof(int));
		for (i = 0; i < (endIndex-startIndex); ++i)
			value_checked[i] = 0;
		fprintf(stdout, "probenum:%ld, checknum:%ld\n", probenum, checknum);\
		for (i = 0; i < (endIndex-startIndex); ++i) {\
			if (value_checked[i] == 0) printf("unchecked: %d, value: %d\n", i, colLValues[i]);\
		}\
		printf("bounds[31]: %d, bounds[32]: %d\n", bounds[31].X, bounds[32].X);\
		printf("imps[95074]: %.4lx\n", im[95074]);\
		free(value_checked);
#endif

#define PROBE_MULTIPASS(T, X) {																						\
    	T  *restrict colLValues = (T *) colL->col;																	\
    	T  *restrict colRValues = (T *) colR->col;																	\
        /*const uint64_t per_thread = colLCount / threadCount;*/	\
		unsigned int vpc = IMPS_PAGE / colL->typesize;\
        const uint64_t startIndex = threadId * num_cacheline_per_thread * vpc;															\
        const uint64_t endIndex = ((threadId == threadCount - 1) ? colLCount : (threadId + 1) * num_cacheline_per_thread * vpc);		\
        /*printf("startindex: %ld, endindex: %ld\n", startIndex, endIndex);	*/\
        /*uint64_t probenum = 0, checknum = 0;*/																	\
        /*uint64_t preJoinCount = 0;			*/																	\
        /*double checksum = 0;					*/																	\
		size_t dentry = 0;		\
		unsigned int maxpart = IMPS_WIDTH;	\
		unsigned int partid_lb = 0;		\
		unsigned int partid_ub = 0;		\
		unsigned int super_partid = 0;	\
		/*size_t lcur;*/\
		uint64_t icnt = 0;	\
		uint64_t top_icnt = 0;\
		uint64_t cache_cnt = 0;\
		uint64_t valueid;\
		uint64_t valueid_lim;\
		Imprints_index * l_imps_idx = colL->imps_idx;														\
		Dct *restrict dct = l_imps_idx->dct;					\
		uint64_t *restrict im = (uint64_t *)l_imps_idx->imprints;	\
        const ValRecord *restrict bounds = l_imps_idx->bounds;												\
		/*const int B = l_imps_idx->bins;*/																		\
		uint64_t lowermask = (ht->mask) >> imps_bits;														\
		unsigned int left_shift = __builtin_popcount(ht->mask) - imps_bits;									\
		unsigned int pattern_bit =  (unsigned int)1 << (MAX_IMPS_BITS - imps_bits) ;		\
		uint64_t pattern = (imps_bits != 0) ? (((uint64_t) 1 << pattern_bit) - 1) : 0xFFFFFFFFFFFFFFFF;\
        uint64_t highermask = 0; 																			\
		T lbound, hbound;	\
        unsigned long bin = 0; 																				\
        uint64_t position = 0; 																				\
        uint64_t found_tuple_oid = 0; 																				\
        unsigned int valid_blks = 0;																											\
		/*printf("inital pattern:%ld\n", pattern);*/\
        for (; super_partid < maxpart/pattern_bit; ++super_partid) {			\
			partid_lb = super_partid * pattern_bit;\
			partid_ub = (super_partid + 1) * pattern_bit - 1;\
			/*printf("partid_lb: %d, partid_ub: %d\n", partid_lb, partid_ub);*/\
			icnt = dctPtrs[threadId].imps_start;	/* the start index in the imprints bit vectors */\
			assert(icnt < l_imps_idx->imps_cnt);\
			cache_cnt = startIndex / vpc;	/* the start index of the cachelines */\
			top_icnt = 0;\
			for (dentry = dctPtrs[threadId].dct_start_entry; dentry <= dctPtrs[threadId].dct_end_entry; dentry++) {\
				\
				/*determine the valid # of imprints bv for current dct entry*/\
				if (dentry == dctPtrs[threadId].dct_end_entry) /* the very last dct entry */ {\
					valid_blks = dctPtrs[threadId].end_valid_bv_num;\
				} else if (dentry == dctPtrs[threadId].dct_start_entry) /* the very first dct entry */ {\
					valid_blks = dct[dentry].blks - dctPtrs[threadId].start_bv_index;\
				} else {\
					valid_blks = dct[dentry].blks;\
				}\
				/*printf("valid_blks:%d, dct[dentry].blks:%d\n", valid_blks, dct[dentry].blks);*/\
				if (!dct[dentry].repeated) { /* not repeated */	\
					top_icnt = icnt + valid_blks;\
					for (; icnt < top_icnt; cache_cnt++, icnt++) {\
						/*if (icnt == 95074) printf("partid: %d, icnt#95074: %ld\n", super_partid, (im[icnt] & pattern));*/\
						if (im[icnt] & pattern) {\
							valueid = cache_cnt * vpc;\
							valueid_lim = valueid + vpc;\
							valueid_lim = valueid_lim > endIndex ? endIndex : valueid_lim;\
							for (; valueid < valueid_lim; valueid++) {\
								PROBE_HT(X);\
								/*checknum++;*/\
								/*if (icnt == 95074) printf("target value:%d, ", colLValues[valueid]);*/\
							}\
						}\
					}\
				} 	\
				else {	/* repeated imprints */ \
					if (im[icnt] & pattern) {\
						valueid = cache_cnt * vpc;\
						valueid_lim = valueid + vpc * valid_blks;\
						valueid_lim = valueid_lim > endIndex ? endIndex : valueid_lim;\
						for (; valueid < valueid_lim; valueid++) {\
							PROBE_HT(X);\
						}\
					}\
					icnt++;\
					cache_cnt += valid_blks;\
				}\
			}\
			/*printf("icnt:%ld\n", icnt);*/\
			pattern = pattern << pattern_bit;\
		}\
    }

#define PROBE_HT(X)\
		do {/*need to replace INT_MIN and INT_MAX somehow */\
			lbound = (partid_lb == 0) ? LONG_MIN : bounds[partid_lb].X;\
			hbound = (partid_ub == 63) ? LONG_MAX : bounds[partid_ub + 1].X;\
			if ((colLValues[valueid] >= lbound) && (colLValues[valueid] < hbound)) {		\
				/*GETBIN_RANGE(bin, *(const TYPE*)v, partid_lb, partid_ub);*/	\
				/*value_checked[valueid] = 1;*/\
				/*probenum++;*/\
				bin = super_partid;					\
				highermask = bin << left_shift;		\
				position = hashKey_imps(colLValues[valueid], lowermask, highermask) & ht->mask;	\
				found_tuple_oid = oid_nil;\
				LOOKUP_POS(position, colLValues[valueid]);\
				/*assert(found_tuple_oid != oid_nil);*/\
				count += (found_tuple_oid != oid_nil);\
			}\
		} while (0)


    BARRIER_ARRIVE(barrier, rv);

    uint64_t count = 0;
    // probe phase ------------------------------------------------------------------
	switch (colL->coltype) {
		//case TYPE_sht: PROBE_MULTIPASS(short, sval); break;
		//case TYPE_int: PROBE_MULTIPASS(int, ival); break;
		case TYPE_lng: PROBE_MULTIPASS(long, lval); break;
		default: break;
	}
	matches[threadId] = count;

	BARRIER_ARRIVE(barrier, rv);


    free(threadArgs);

	if (threadId == 0) {
		free(ht->table);

		timer_stop(&probe_timer);
		printf("[INFO] PROBE phase takes %lf ms\n", GetMilliSeconds(&probe_timer));

		//fprintf(stdout, "count:%ld\n", count);
	}
	return NULL;
}

void* run_orderlist(void *arg) {

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

    ThreadDctPointer* dctPtrs = threadArgs->dctPointers;
    uint64_t* matches = threadArgs->matches;
    unsigned int num_cacheline_per_thread = threadArgs->num_cacheline_per_thread;

	unsigned long colLCount = colL->colcount;
	unsigned long colRCount = colR->colcount;

	int rv;

	/* initialize hash table */
	uint64_t htSize = ht->nBuckets * sizeof(HTEntry);
	if (threadId == 0) {
		VERBOSE fprintf(stdout, "HT: initializing memory in thread 0\n");
		ht->table = (HTEntry *) malloc_aligned(htSize);
		memset(ht->table, 0xFF, htSize);
	}
	BARRIER_ARRIVE(barrier, rv);

#define BUILD_ORDERLIST(T, X) {																				\
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

    	//struct timeval start_timeval, build_timeval;
		HybridTimer build_timer;
		HybridTimer probe_timer;
		HybridTimer lvisit_timer;

        if (threadId == 0) {
        	//gettimeofday(&start_timeval, NULL);
        	timer_start(&build_timer);
        }

    	/* build phase ------------------------------------------------------------------*/
    	switch (colR->coltype) {
    		case TYPE_sht: BUILD_ORDERLIST(short, sval); break;
    		case TYPE_int: BUILD_ORDERLIST(int, ival); break;
    		case TYPE_lng: BUILD_ORDERLIST(long, lval); break;
    		default: break;
    	}

        //syncThreads();
        BARRIER_ARRIVE(barrier, rv);

        if (threadId == 0) {
        	//gettimeofday(&build_timeval, NULL);
        	timer_stop(&build_timer);
        	printf("[INFO] BUILD phase takes %lf ms\n", GetMilliSeconds(&build_timer));
        	timer_start(&probe_timer);
        }

#define PROBE_ORDERLIST(T, X) {																						\
    	T  *restrict colLValues = (T *) colL->col;																	\
    	T  *restrict colRValues = (T *) colR->col;																	\
        /*const uint64_t per_thread = colLCount / threadCount;*/	\
		unsigned int vpc = IMPS_PAGE / colL->typesize;\
        const uint64_t startIndex = threadId * num_cacheline_per_thread * vpc;															\
        const uint64_t endIndex = ((threadId == threadCount - 1) ? colLCount : (threadId + 1) * num_cacheline_per_thread * vpc);		\
        /*printf("startindex: %ld, endindex: %ld\n", startIndex, endIndex);	*/\
        /*uint64_t probenum = 0, checknum = 0;*/																	\
        /*uint64_t preJoinCount = 0;			*/																	\
        /*double checksum = 0;					*/																	\
		size_t dentry = 0;		\
		/*unsigned int maxpart = IMPS_WIDTH;*/	\
		unsigned int partid_lb = 0;		\
		unsigned int partid_ub = 0;		\
		unsigned int super_partid = 0;	\
		/*size_t lcur;*/\
		uint64_t icnt = 0;	\
		uint64_t top_icnt = 0;\
		uint32_t cache_cnt = 0;\
		uint32_t start_cache_cnt = startIndex / vpc;\
		uint64_t valueid;\
		uint64_t valueid_lim;\
		Imprints_index * l_imps_idx = colL->imps_idx;														\
		Dct *restrict dct = l_imps_idx->dct;					\
		uint64_t *restrict im = (uint64_t *)l_imps_idx->imprints;	\
        const ValRecord *restrict bounds = l_imps_idx->bounds;												\
		/*const int B = l_imps_idx->bins;*/																		\
		uint64_t lowermask = (ht->mask) >> imps_bits;														\
		unsigned int left_shift = __builtin_popcount(ht->mask) - imps_bits;									\
		unsigned int pattern_bit =  (unsigned int)1 << (MAX_IMPS_BITS - imps_bits) ;		\
		uint64_t pattern = (imps_bits != 0) ? (((uint64_t) 1 << pattern_bit) - 1) : 0xFFFFFFFFFFFFFFFF;\
		uint64_t ori_pattern = pattern;	\
        uint64_t highermask = 0; 																			\
        uint32_t partn = (unsigned int)1 << imps_bits;\
		T lbound, hbound;	\
        unsigned long bin = 0; 																				\
        uint64_t position = 0; 																				\
        uint64_t found_tuple_oid = 0; 																				\
        unsigned int valid_blks = 0;																											\
        int i, cache_cnt_iter, lpos, exclusive_cnt = 0;\
        timer_start(&lvisit_timer);\
        uint8_t* is_exclusive = (uint8_t *) malloc_aligned(num_cacheline_per_thread * sizeof(uint8_t));\
        Array * lvisit = (Array *) malloc_aligned(partn * sizeof(Array));\
        for (i = 0; i < partn; ++i) {\
        	initArray(lvisit + i, num_cacheline_per_thread / partn);\
        }\
		/*printf("inital pattern:%ld\n", pattern);*/\
		icnt = dctPtrs[threadId].imps_start;	/* the start index in the imprints bit vectors */\
		cache_cnt = startIndex / vpc;	/* the start index of the cachelines */\
		top_icnt = 0;\
		for (dentry = dctPtrs[threadId].dct_start_entry; dentry <= dctPtrs[threadId].dct_end_entry; dentry++) {\
			\
			/*determine the valid # of imprints bv for current dct entry*/\
			if (dentry == dctPtrs[threadId].dct_end_entry) /* the very last dct entry */ {\
				valid_blks = dctPtrs[threadId].end_valid_bv_num;\
			} else if (dentry == dctPtrs[threadId].dct_start_entry) /* the very first dct entry */ {\
				valid_blks = dct[dentry].blks - dctPtrs[threadId].start_bv_index;\
			} else {\
				valid_blks = dct[dentry].blks;\
			}\
			/*printf("valid_blks:%d, dct[dentry].blks:%d\n", valid_blks, dct[dentry].blks);*/\
			if (!dct[dentry].repeated) { /* not repeated */	\
				top_icnt = icnt + valid_blks;\
				for (; icnt < top_icnt; cache_cnt++, icnt++) {\
					pattern = ori_pattern;\
					for (i = 0; i < partn; ++i) {\
						if (pattern & im[icnt]) {\
							insertArray(lvisit + i, cache_cnt);\
							/* judge is exclusive for this pattern */\
							assert((cache_cnt >= start_cache_cnt) && ((cache_cnt - start_cache_cnt) < num_cacheline_per_thread));\
							/*assert(!((~pattern) & im[icnt]));*/\
							/*if ((~pattern) & im[icnt]) {printf("not exclusive:%d: %04lx\n", cache_cnt, im[icnt]);}*/\
							exclusive_cnt += (!((~pattern) & im[icnt]));\
        					is_exclusive[cache_cnt - start_cache_cnt] = !((~pattern) & im[icnt]);\
						}\
						pattern = pattern << pattern_bit;\
					}\
				}\
			} 	\
			else {	/* repeated imprints */ \
				pattern = ori_pattern;\
				for (i = 0; i < partn; ++i) {\
					if (im[icnt] & pattern) {\
						for (cache_cnt_iter = cache_cnt; \
								cache_cnt_iter < cache_cnt + valid_blks;\
								++cache_cnt_iter) {\
							insertArray(lvisit + i, cache_cnt_iter);\
							is_exclusive[cache_cnt_iter - start_cache_cnt] = !((~pattern) & im[icnt]);\
						}\
					}\
					pattern = pattern << pattern_bit;\
				}\
				icnt++;\
				cache_cnt += valid_blks;\
			}\
		}\
		timer_stop(&lvisit_timer);\
		printf("[INFO] Create visit list takes %lf ms\n", GetMilliSeconds(&lvisit_timer));\
		DEBUG printf("exclusive percentage: %lf\n", 1.0*exclusive_cnt / num_cacheline_per_thread);\
		\
        /* probing according to the visit list */\
		icnt = 0;\
        for (super_partid = 0; super_partid < partn; ++super_partid) {\
        	partid_lb = super_partid * pattern_bit;\
        	partid_ub = (super_partid + 1) * pattern_bit - 1;\
        	for (lpos = 0; lpos < lvisit[super_partid].used; ++lpos) {\
        		icnt = lvisit[super_partid].array[lpos];\
        		valueid = icnt * vpc;\
        		valueid_lim = valueid + vpc;\
        		valueid_lim = valueid_lim > endIndex ? endIndex : valueid_lim;\
        		if (is_exclusive[icnt - start_cache_cnt]) {\
					for (; valueid < valueid_lim; valueid++) {\
						PROBE_HT_EXCLUSIVE(X);\
					}\
        		} else {\
					for (; valueid < valueid_lim; valueid++) {\
						PROBE_HT(X);\
					}\
        		}\
        	}\
        }\
		\
		for (i = 0; i < partn; ++i) {\
			freeArray(lvisit + i);\
		}\
		free(lvisit);\
		free(is_exclusive);\
    }

#define PROBE_HT_EXCLUSIVE(X)\
		do {/*need to replace INT_MIN and INT_MAX somehow */\
			/*GETBIN_RANGE(bin, *(const TYPE*)v, partid_lb, partid_ub);*/	\
			/*value_checked[valueid] = 1;*/\
			/*probenum++;*/\
			bin = super_partid;					\
			highermask = bin << left_shift;		\
			position = hashKey_imps(colLValues[valueid], lowermask, highermask) & ht->mask;	\
			found_tuple_oid = oid_nil;\
			LOOKUP_POS(position, colLValues[valueid]);\
			/*assert(found_tuple_oid != oid_nil);*/\
			count += (found_tuple_oid != oid_nil);\
		} while (0)

    BARRIER_ARRIVE(barrier, rv);

    uint64_t count = 0;
    // probe phase ------------------------------------------------------------------
	switch (colL->coltype) {
		//case TYPE_sht: PROBE_ORDERLIST(short, sval); break;
		//case TYPE_int: PROBE_ORDERLIST(int, ival); break;
		case TYPE_lng: PROBE_ORDERLIST(long, lval); break;
		default: break;
	}
	matches[threadId] = count;

	BARRIER_ARRIVE(barrier, rv);


    free(threadArgs);

	if (threadId == 0) {
		free(ht->table);

		timer_stop(&probe_timer);
		printf("[INFO] PROBE phase takes %lf ms\n", GetMilliSeconds(&probe_timer));

		//fprintf(stdout, "count:%ld\n", count);
	}
	return NULL;
}

#if 0
for (; super_partid < maxpart/pattern_bit; ++super_partid) {			\
	partid_lb = super_partid * pattern_bit;\
	partid_ub = (super_partid + 1) * pattern_bit - 1;\
	/*lcur = lstart;*/\
	icnt = 0;\
	cache_cnt = 0;\
	top_icnt = 0;\
	for (dentry = 0; dentry < l_imps_idx->dct_cnt; dentry++) {\
		if (!dct[dentry].repeated) {	\
			top_icnt = icnt + dct[dentry].blks;\
			for (; icnt < top_icnt; cache_cnt++, icnt++) {\
				if (im[icnt] & pattern) {\
					valueid = cache_cnt * vpc;\
					valueid_lim = valueid + vpc;\
					valueid_lim = valueid_lim > total_cnt ? total_cnt : valueid_lim;\
					for (; valueid < valueid_lim; valueid++) {\
						PROBE_HT(X);\
					}\
				}\
			}\
		} 	\
		else {	\
			if (im[icnt] & pattern) {\
				valueid = cache_cnt * vpc;\
				valueid_lim = valueid + vpc * dct[dentry].blks;\
				valueid_lim = valueid_lim > total_cnt ? total_cnt : valueid_lim;\
				for (; valueid < valueid_lim; valueid++) {\
					PROBE_HT(X);\
				}\
			}\
			icnt++;\
			cache_cnt += dct[dentry].blks;\
		}\
	}\
	assert(icnt == l_imps_idx->imps_cnt);\
	pattern = pattern << pattern_bit;\
}\
}

#define PROBE_HT(X)\
do {/*need to replace INT_MIN and INT_MAX somehow */\
	lbound = (partid_lb == 0) ? INT_MIN : bounds[partid_lb].X;\
	hbound = (partid_ub == 63) ? INT_MAX : bounds[partid_ub+1].X;\
	if ((colLValues[startIndex+valueid] >= lbound) && (colLValues[startIndex+valueid] < hbound)) {		\
		/*GETBIN_RANGE(bin, *(const TYPE*)v, partid_lb, partid_ub);*/	\
		/*probenum++;*/\
		bin = super_partid;					\
		highermask = bin << left_shift;		\
		position = hashKey_imps(colLValues[startIndex+valueid], lowermask, highermask) & ht->mask;	\
		found_tuple_oid = oid_nil;\
		LOOKUP_POS(position, colLValues[startIndex+valueid]);\
		count += (found_tuple_oid != oid_nil);\
	}\
} while (0)
#endif


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

	uint64_t* matches = threadArgs->matches;

	unsigned long colLCount = colL->colcount;
	unsigned long colRCount = colR->colcount;

	int rv;

	/* initialize hash table */
	uint64_t htSize = ht->nBuckets * sizeof(HTEntry);
	if (threadId == 0) {
		VERBOSE fprintf(stdout, "HT: initializing memory in thread 0\n");
		ht->table = (HTEntry *) malloc_aligned(htSize);
		memset(ht->table, 0xFF, htSize);
	}
	DEBUG printf("first value: %04x, hash_nil: %04x\n", ht->table[0].hash, hash_nil);
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

    	//struct timeval start_timeval, build_timeval;
		HybridTimer build_timer;
		HybridTimer probe_timer;
        if (threadId == 0) {
        	//gettimeofday(&start_timeval, NULL);
        	timer_start(&build_timer);
        }

    	/* build phase ------------------------------------------------------------------*/
    	switch (colR->coltype) {
    		//case TYPE_sht: BUILD(short); break;
    		case TYPE_int: BUILD(int); break;
    		case TYPE_lng: BUILD(long); break;
    		default: break;
    	}

        //syncThreads();
        BARRIER_ARRIVE(barrier, rv);

        if (threadId == 0) {
        	//gettimeofday(&build_timeval, NULL);
        	timer_stop(&build_timer);
        	printf("[INFO] BUILD phase takes %lf ms\n", GetMilliSeconds(&build_timer));
        	timer_start(&probe_timer);
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
			/*DEBUG if (found_tuple_oid == oid_nil) fprintf(stdout, "probe failed:%d\n", colLValues[i]);*/\
        \
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
		//case TYPE_sht: PROBE(short); break;
		case TYPE_int: PROBE(int); break;
		case TYPE_lng: PROBE(long); break;
		default: break;
	}
	matches[threadId] = count;

	BARRIER_ARRIVE(barrier, rv);

    free(threadArgs);

	if (threadId == 0) {
		free(ht->table);

    	timer_stop(&probe_timer);
    	printf("[INFO] PROBE phase takes %lf ms\n", GetMilliSeconds(&probe_timer));

    	DEBUG fprintf(stdout, "count:%ld\n", count);
	}
	return NULL;
}

void* run_default_no_prefetching(void *arg) {

	JoinThreadArgs* threadArgs = (JoinThreadArgs*)(arg);

	unsigned int threadId = threadArgs->threadId;
	unsigned int threadCount = threadArgs->threadCount;
	Column* colL = threadArgs->colL;
	Column* colR = threadArgs->colR;
	HashTable * ht = threadArgs->ht;
	pthread_barrier_t * barrier = threadArgs->barrier;

	uint64_t* matches = threadArgs->matches;

	unsigned long colLCount = colL->colcount;
	unsigned long colRCount = colR->colcount;

	int rv;

	/* initialize hash table */
	uint64_t htSize = ht->nBuckets * sizeof(HTEntry);
	if (threadId == 0) {
		VERBOSE fprintf(stdout, "HT: initializing memory in thread 0\n");
		ht->table = (HTEntry *) malloc_aligned(htSize);
		memset(ht->table, 0xFF, htSize);
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

    	//struct timeval start_timeval, build_timeval;
		HybridTimer build_timer;
		HybridTimer probe_timer;
        if (threadId == 0) {
        	//gettimeofday(&start_timeval, NULL);
        	timer_start(&build_timer);
        }

    	/* build phase ------------------------------------------------------------------*/
    	switch (colR->coltype) {
    		//case TYPE_sht: BUILD(short); break;
    		case TYPE_int: BUILD(int); break;
    		//case TYPE_lng: BUILD(long); break;
    		default: break;
    	}

        //syncThreads();
        BARRIER_ARRIVE(barrier, rv);

        if (threadId == 0) {
        	//gettimeofday(&build_timeval, NULL);
        	timer_stop(&build_timer);
        	printf("[INFO] BUILD phase takes %lf ms\n", GetMilliSeconds(&build_timer));
        	timer_start(&probe_timer);
        }

#define PROBE_NO_PREFETCH(T) {																									\
    	T  *restrict colLValues = (T *) colL->col;																	\
    	T  *restrict colRValues = (T *) colR->col;																	\
        const uint64_t per_thread = colLCount / threadCount;														\
        const uint64_t startIndex = threadId * per_thread;															\
        const uint64_t endIndex = ((threadId == threadCount - 1) ? colLCount : (threadId + 1) * per_thread);		\
        /*uint64_t count = 0;*/																						\
        /*uint64_t preJoinCount = 0;			*/																	\
        /*double checksum = 0;					*/																	\
        PrefetchCache cacheEntry;																					\
																											\
																													\
		\
		for (uint64_t i = startIndex; i < endIndex; i++) {													\
			/* re-use hash values and position, which have been already computed during prefetching */				\
																			\
			uint64_t found_tuple_oid = oid_nil;																		\
			HASH(colLValues[i], cacheEntry);																		\
			\
			LOOKUP_POS(cacheEntry.pos, colLValues[i]);														\
			DEBUG if (found_tuple_oid == oid_nil) fprintf(stdout, "probe failed:%d\n", colLValues[i]);				\
																													\
        	count += (found_tuple_oid != oid_nil);																	\
																													\
		}																											\
    }


    BARRIER_ARRIVE(barrier, rv);

    uint64_t count = 0;
    // probe phase ------------------------------------------------------------------
	switch (colL->coltype) {
		//case TYPE_sht: PROBE(short); break;
		case TYPE_int: PROBE_NO_PREFETCH(int); break;
		//case TYPE_lng: PROBE(long); break;
		default: break;
	}
	matches[threadId] = count;

	BARRIER_ARRIVE(barrier, rv);

    free(threadArgs);

	if (threadId == 0) {
		free(ht->table);

    	timer_stop(&probe_timer);
    	printf("[INFO] PROBE phase takes %lf ms\n", GetMilliSeconds(&probe_timer));

    	DEBUG fprintf(stdout, "count:%ld\n", count);
	}
	return NULL;
}

query_result_t hashjoin(Column * colL, Column * colR, unsigned int threadCount, joinconf_t joinCfg) {
	query_result_t result = {0, 0, 0.0, 0, 0, 0};

	/* step 1: initialize the hash table */
	unsigned long colRCount = colR->colcount;
    uint64_t htSize = nextPowerOfTwo(colRCount) << 1;	/* why move one more bits here ?*/
	//uint64_t htSize = nextPowerOfTwo(colRCount);
    uint64_t htBits = unsignedLog2(htSize);

    DEBUG printf("htBits: %ld, sizeof(HTEntry):%ld\n", htBits, sizeof(HTEntry));

    HashTable * ht = (HashTable *) malloc_aligned(sizeof(HashTable));
    //initHashTable(ht, htBits, threadCount);
	ht->nBuckets = 1ull << htBits;
	ht->mask = ht->nBuckets - 1;

    pthread_barrier_t barrier;
    pthread_barrier_init(&barrier, NULL, threadCount);

    uint64_t * matches = (uint64_t *) malloc_aligned(threadCount * sizeof(uint64_t));/* to store result num*/

    /* step 2: initialize pointer to imprints for each partition of colL (the probe column) */
    unsigned long colLCount = colL->colcount;
    unsigned int vpc = IMPS_PAGE / colL->typesize;
    //unsigned int num_cacheline_per_thread = ceil(colLCount*1.0 / threadCount / vpc);
    unsigned int num_cacheline_per_thread = ceil(ceil(colLCount*1.0 / vpc) / threadCount);
    unsigned int num_cacheline_last_thread = ceil(colLCount*1.0 / vpc) - num_cacheline_per_thread * (threadCount - 1);
    unsigned int accu_cacheline_num = 0;
    unsigned int cur_thread = 0;
	uint64_t icnt = 0;		/* index of imprints bit vector */
	uint64_t dct_entry_remain_cnt = 0;
	uint64_t thread_remain_cnt = 0;
	uint64_t bv_consumed_dct_entry = 0;
	Imprints_index * l_imps_idx = colL->imps_idx;
	size_t dentry;
	Dct *restrict dct = l_imps_idx->dct;

    ThreadDctPointer * dctPointers = (ThreadDctPointer *) malloc_aligned(
    		threadCount * sizeof(ThreadDctPointer));

	dctPointers[0].dct_start_entry = 0;
	dctPointers[0].start_bv_index = 0;
	dctPointers[0].imps_start = 0;

	for (dentry = 0; dentry < l_imps_idx->dct_cnt; dentry++) {
		if (!dct[dentry].repeated) /* not repeated */{
			//top_icnt = icnt + dct[dentry].blks;

			dct_entry_remain_cnt = dct[dentry].blks;

			bv_consumed_dct_entry = 0;	/* # of bit vectors consumed within current entry */

			thread_remain_cnt = ((cur_thread == threadCount-1) ?
					num_cacheline_last_thread:
					num_cacheline_per_thread) - accu_cacheline_num;

			while (thread_remain_cnt <= dct_entry_remain_cnt) /* current dct entry could fulfill requirement...*/{
				dctPointers[cur_thread].dct_end_entry = dentry;
				dctPointers[cur_thread].end_valid_bv_num = thread_remain_cnt;
				dct_entry_remain_cnt -= thread_remain_cnt;
				accu_cacheline_num = 0;
				bv_consumed_dct_entry += thread_remain_cnt;
				cur_thread++;

				if (cur_thread < threadCount) {
					dctPointers[cur_thread].dct_start_entry =
							(dct_entry_remain_cnt == 0) ? (dentry + 1) : dentry;
					dctPointers[cur_thread].start_bv_index =
							(dct_entry_remain_cnt == 0) ? 0 : bv_consumed_dct_entry;
					dctPointers[cur_thread].imps_start = icnt + bv_consumed_dct_entry;
				}
			}

			accu_cacheline_num += dct_entry_remain_cnt;
			icnt += dct[dentry].blks;
		}
		else {	/* repeated */
			dct_entry_remain_cnt = dct[dentry].blks;

			bv_consumed_dct_entry = 0;	/* # of bit vectors consumed within current entry */

			thread_remain_cnt = ((cur_thread == threadCount-1) ?
					num_cacheline_last_thread:
					num_cacheline_per_thread) - accu_cacheline_num;

			while (thread_remain_cnt <= dct_entry_remain_cnt) {
				dctPointers[cur_thread].dct_end_entry = dentry;
				dctPointers[cur_thread].end_valid_bv_num = thread_remain_cnt;
				dct_entry_remain_cnt -= thread_remain_cnt;
				accu_cacheline_num = 0;
				bv_consumed_dct_entry += thread_remain_cnt;
				cur_thread++;

				if (cur_thread < threadCount) {
					dctPointers[cur_thread].dct_start_entry =
							(dct_entry_remain_cnt == 0) ? (dentry + 1) : dentry;
					dctPointers[cur_thread].start_bv_index =
							(dct_entry_remain_cnt == 0) ? 0 : bv_consumed_dct_entry;
					dctPointers[cur_thread].imps_start = icnt;	/* different with <not repeated> */
				}
			}

			accu_cacheline_num += dct_entry_remain_cnt;
			icnt += 1;	/* different with <not repeated> */
		}
	}

    /* step 3: multi-threaded build and probe */
	pthread_t threads[50];

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
		threadArgs->dctPointers = dctPointers;
		threadArgs->matches	= matches;
		threadArgs->num_cacheline_per_thread = num_cacheline_per_thread;
		switch(joinCfg.joinalgo) {
		case 0:	//default
			printf("[INFO] start default hash-join\n");
			pthread_create(&threads[i], NULL, run_default, (void*)threadArgs);
			break;
		case 1:	//multi-pass
			printf("[INFO] start imps-hashjoin with multipass\n");
			pthread_create(&threads[i], NULL, run_multipass, (void*)threadArgs);
			break;
		case 2:	//order-list
			printf("[INFO] start imps-hashjoin with orderlist\n");
			pthread_create(&threads[i], NULL, run_orderlist, (void*)threadArgs);
			break;
		case 3: //default with no prefetching
			printf("[INFO] start default hash-join without prefetching\n");
			pthread_create(&threads[i], NULL, run_default_no_prefetching, (void*)threadArgs);
		default:
			break;
		}
	}
	for (unsigned i = 0; i < threadCount; i++) {
		pthread_join(threads[i], NULL);
	}
	//query_result_t result = parallelJoin->getResult();
	/* to aggregate the join result */
	uint64_t total_num = 0;
	for (int i = 0; i < threadCount; ++i) {
		total_num += matches[i];
	}
	fprintf(stdout, "[INFO] Total Result Num: %ld\n", total_num);

	free(ht);
	free(dctPointers);
	free(matches);

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

#if 0
do {								\
	lcur = lstart + valueid;\
	lo = lcur + l->hseqbase;\
	v = FVALUE(l, lcur);				\
	nr = 0;						\
	lbound = (partid_lb == 0) ? INT_MIN : bins[partid_lb];\
	hbound = (partid_ub == 63) ? INT_MAX : bins[partid_ub+1];\
	if ((*(const TYPE*)v != TYPE##_nil) && (*(const TYPE*)v >= lbound) && (*(const TYPE*)v < hbound)) {		\
		probenum++;					\
		bin = PART;					\
		highermask = bin << left_shift;		\
		for (rb = HASHget##WIDTH(hsh, hash_imps_##TYPE(lowermask, v, highermask)); \
			 rb != hashnil;			\
			 /*rb = HASHgetlink##WIDTH(hsh, rb))*/	\
			 rb = HASHgetlink##WIDTH(hsh, rb), collisionlen++)	\
			if (rb >= rl && rb < rh &&	\
				* (const TYPE *) v == ((const TYPE *) base)[rb]) { \
				ro = (oid) (rb - rl + rseq); \
				HASHLOOPBODY();		\
				result_num++;		\
			}				\
	}
#endif

