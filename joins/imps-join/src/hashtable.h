/*
 * hashtable.h
 *
 *  Created on: 6 Jun 2017
 *      Author: zeroxwj
 */

#ifndef HASHTABLE_H_
#define HASHTABLE_H_

//#define _IDHASH_ 1
#define _FIBHASH_ 1

typedef struct  {
    uint32_t hash;
    uint64_t oid;
}HTEntry ;

typedef struct {
    uint64_t mask;
    HTEntry *table;
    uint64_t count;
    //bool lockFree;
    //Lock *locks;
    size_t locksLength;
    size_t lockMask;
    //size_t striping;
    uint64_t nBuckets;
	//PThreadLockCVBarrier *barrier;
} HashTable;

typedef struct {
    uint32_t hash;
    uint64_t pos;
    uint32_t rowID;	//??
} PrefetchCache;

#if defined(_IDHASH_)
	/** Identity Hashing */
	/*
	inline intkey_t hashKey(const intkey_t k) const {
		return (k & mask);
	}
	*/
#define hashKey(KEY) (KEY & ht->mask)
#define hashKey_imps(KEY, LMASK, HMASK) (((KEY) & LMASK) | HMASK)


#elif defined(_FIBHASH_)
	/** Fibonacci Hashing */
	/*
	inline intkey_t hashKey(const intkey_t k) const {
		return (k * 11400714819323198485ull) | (1ull<<((sizeof(intkey_t)*8-1)));
	}
	*/
#define hashKey(KEY) (((KEY * 11400714819323198485ull) | (1ull<<((sizeof(uint32_t)*8-1)))) & ht->mask)
#define hashKey_imps(KEY, LMASK, HMASK) ((((KEY * 11400714819323198485ull) | (1ull<<((sizeof(uint32_t)*8-1)))) & LMASK) | HMASK)

#elif defined(_CRCHASH_)
	/** CRC Hashing */
	inline intkey_t hashKey(const intkey_t k) const {
		return _mm_crc32_u64(0, k) | (1ull<<((sizeof(intkey_t)*8-1)));
	}
#else

    /** MurmurHash64A */
	/*
    inline intkey_t hashKey(intkey_t k) const {
        const intkey_t m = 0xc6a4a7935bd1e995;
        const int r = 47;
        intkey_t h = 0x8445d61a4e774912 ^(8 * m);
        k *= m;
        k ^= k >> r;
        k *= m;
        h ^= k;
        h *= m;
        h ^= h >> r;
        h *= m;
        h ^= h >> r;
        return h | (1ull << ((sizeof(intkey_t) * 8 - 1)));
    }
    */

#endif


#endif /* HASHTABLE_H_ */
