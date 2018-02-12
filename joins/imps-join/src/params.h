/**
 * @file    params.h
 * @author  Wenjian Xu
 *
 */
#ifndef PARAMS_H_
#define PARAMS_H_

#define IMPS_PAGE	64		/* the block size each bit vector represents */
#define IMPS_WIDTH	64
#define MAX_IMPS_BITS 6



#ifdef OID_WIDTH_8B /* 64-bit surrogate_id*/
typedef uint64_t surrogate_t;
#else /* 32-bit surrogate_id*/
typedef uint32_t surrogate_t;
#endif

/** The partitioning fan-out for the inital step of sort-merge joins */
#ifndef NRADIXBITS_DEFAULT
#define NRADIXBITS_DEFAULT 7
#endif

/** Default partitioning fan-out, can be adjusted from command line. */
#ifndef PARTFANOUT_DEFAULT
#define PARTFANOUT_DEFAULT (1<<NRADIXBITS_DEFAULT)
#endif

#define VALUESSPERCACHELINE(T) (CACHE_LINE_SIZE/sizeof(T))

/** System cache line size */
#ifndef CACHE_LINE_SIZE
#define CACHE_LINE_SIZE 64
#endif
//constexpr uint32_t CACHE_LINE_SIZE = 64;


#define LH_MIN_SIZE_LOG_2 6

#define LH_MIN_LOAD_FACTOR .15f
#define LH_MAX_LOAD_FACTOR .5f
//#define LH_MIN_LOAD_FACTOR .275f
//#define LH_MAX_LOAD_FACTOR .675f

//#define DYNAMIC_GROW

#endif /* PARAMS_H_ */
