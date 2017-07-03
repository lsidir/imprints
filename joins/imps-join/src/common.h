/**
 */


#ifndef COMMON_H_
#define COMMON_H_

#include <stdint.h>
#include <stdlib.h>             /* posix_memalign, EXIT_FAILURE */
#include <stdio.h>              /* FILE */
#include <x86intrin.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>
//#include <ammintrin.h>
//#include <smmintrin.h>
//#include <immintrin.h>
#include <assert.h>
#include <time.h>
#include <sys/time.h>
#include <math.h>
#include <string.h>

#include "params.h"             /* macro parameters */


inline __attribute__((__always_inline__))
void *malloc_aligned(size_t size)
{
    void * ret;
    int rv;
    rv = posix_memalign((void**)&ret, CACHE_LINE_SIZE, size);

    if (rv) {
        printf("[ERROR] malloc_aligned() failed: out of memory\n");
        return 0;
    }

    return ret;
}

#if 0
/**
 * Makes a non-temporal write of 64 bytes from src to dst.
 * Uses vectorized non-temporal stores if available, falls
 * back to assignment copy.
 *
 * @param dst
 * @param src
 *
 * @return
 */
static inline void
store_nontemp_64B(void * dst, void * src)
{
//#ifdef __AVX__
    register __m256i * d1 = (__m256i*) dst;
    register __m256i s1 = *((__m256i*) src);
    register __m256i * d2 = d1+1;
    register __m256i s2 = *(((__m256i*) src)+1);

    _mm256_stream_si256(d1, s1);
    _mm256_stream_si256(d2, s2);
//#elif defined(__SSE2__)
#if 0
    register __m128i * d1 = (__m128i*) dst;
    register __m128i * d2 = d1+1;
    register __m128i * d3 = d1+2;
    register __m128i * d4 = d1+3;
    register __m128i s1 = *(__m128i*) src;
    register __m128i s2 = *((__m128i*)src + 1);
    register __m128i s3 = *((__m128i*)src + 2);
    register __m128i s4 = *((__m128i*)src + 3);

    _mm_stream_si128 (d1, s1);
    _mm_stream_si128 (d2, s2);
    _mm_stream_si128 (d3, s3);
    _mm_stream_si128 (d4, s4);

//#else
    /* just copy with assignment */
    *(cacheline_t *)dst = *(cacheline_t *)src;

#endif
}

template <class T>
inline __attribute__((__always_inline__)) void
swap(T ** A, T ** B)
{
    T * tmp = *A;
    *A = *B;
    *B = tmp;
}

inline __attribute__((__always_inline__)) uint32_t
mylog2(const uint32_t n) 
{
    register uint32_t res;
    __asm__ ( "\tbsr %1, %0\n" : "=r"(res) : "r"(n) );
    return res;
}

template <typename T>
inline T min(T a, T b) {
	return a < b ? a : b;
}

template <typename T>
inline T max(T a, T b) {
	return a > b ? a : b;
}
#endif

#endif /* COMMON_H_ */
