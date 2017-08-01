/**
 * @file    types.h
 * 
 * @brief  Provides general type definitions used by all join algorithms.
 * 
 * 
 */
#ifndef TYPES_H
#define TYPES_H

#include <stdint.h>

/**
 * @defgroup Types Common Types
 * Common type definitions used by all join implementations.
 * @{
 */
#if 0
#ifdef KEY_8B /* 64-bit key/value, 16B tuples */
typedef int64_t intkey_t;
typedef int64_t value_t;
#else /* 32-bit key/value, 8B tuples */
typedef int32_t intkey_t;
typedef int32_t value_t;
#endif

typedef struct tuple_t    tuple_t;
typedef struct relation_t relation_t;

/** Type definition for a tuple, depending on KEY_8B a tuple can be 16B or 8B */
struct tuple_t {
    intkey_t key;
    value_t  payload;
};

/**
 * Type definition for a relation. 
 * It consists of an array of tuples and a size of the relation.
 */
struct relation_t {
  tuple_t * tuples;
  uint32_t  num_tuples;
};

typedef struct result_t result_t;
typedef struct joinconfig_t joinconfig_t;
typedef struct threadresult_t threadresult_t;

/** Holds the join results of a thread */
struct threadresult_t {
    int64_t  nresults;
    void *   results;
    uint32_t threadid;
};

/** Type definition for join results. */
struct result_t {
    int64_t          totalresults;
    threadresult_t * resultlist;
    int              nthreads;
};


/**
 * Various NUMA shuffling strategies for data shuffling phase of join
 * algorithms as also described by NUMA-aware data shuffling paper [CIDR'13].
 *
 * NUMA_SHUFFLE_RANDOM, NUMA_SHUFFLE_RING, NUMA_SHUFFLE_NEXT
 */
enum numa_strategy_t {RANDOM, RING, NEXT};

/** Join configuration parameters. */
struct joinconfig_t {
    int NTHREADS;
    int PARTFANOUT;
    int SCALARSORT;
    int SCALARMERGE;
    int MWAYMERGEBUFFERSIZE;
    enum numa_strategy_t NUMASTRATEGY;
};

/** @} */
#endif

/* printing options */
#define VERBOSE        if (1)
#define TIMING		   if (1)
#define DEBUG		   if (1)


/* Using MonetDB style types */
#define TYPE_void   0
#define TYPE_bte    3
#define TYPE_sht    4
#define TYPE_int    6
#define TYPE_oid    7
#define TYPE_flt    10
#define TYPE_dbl    11
#define TYPE_lng    12
#define TYPE_str    13

typedef struct {
    uint64_t matches;
    uint64_t preJoinTuples;
    double checksum;
    uint64_t time_cycles;
    uint64_t part_cycles;
    uint64_t join_cycles;
} query_result_t;


typedef struct {
	unsigned int joinalgo;
	unsigned int imps_bits;
} joinconf_t;

typedef union {
	char bval;
	short sval;
	int ival;
	long lval;
	unsigned long ulval;
	float fval;
	double dval;
} ValRecord;

typedef struct {
	unsigned long 	dct_start_entry;	/* dct start entry (start from 0) */
	unsigned long	dct_end_entry;		/* dct end entry */
	//unsigned int 	residual_blks;		/* remaining imps to be traversed for start entry */
	unsigned int	start_bv_index;		/* the start index of bv in dct_start_entry */
	unsigned int	end_valid_bv_num;	/* valid num of bv in dct_end_entry */
	unsigned long 	imps_start;			/* start index of imprints bit vectors */
	//unsigned int	flag_single_entry;	/* 1 indicate single entry, i.e., dct_start_entry == dct_end_entry, otherwise 0*/
} ThreadDctPointer;

/* The Column Imprints contains a binned imprint bitvector
 * to weed out blocks of no interest in a scan.
 * The imprint vectors are simply compressed, keeping track on
 * the number of blocks each imprint covers. From the blks we can
 * calculate the actual oid ranges, provided we scan only.
 */
#define MAXOFFSET 24
typedef struct {
	unsigned int   blks:MAXOFFSET;	/* blocks count						*/
	unsigned int   repeated:1;		/* same or unique imprints in blks	*/
	unsigned int   flgs:8 * sizeof(int) - MAXOFFSET-1; /* for future use	*/
} Dct;

typedef struct {
	Dct            *dct;		/* the dictionary structure				*/
	char           *imprints;	/* the imprint vectors					*/
	ValRecord      *bounds;		/* the bounds of the binning			*/
	unsigned long  dct_cnt;		/* count of dictionary entries			*/
	unsigned long  imps_cnt;	/* count of imprint vector entries		*/
	int            bins;		/* number of bins						*/
	int            imprintsize;	/* bytes per imprint					*/
	int            blocksize;	/* bytes per block						*/
} Imprints_index;

typedef struct {
	char          filename[1024];
	char          colname[1024];
	char          typename[64];
	int           coltype;			/* as given by TYPE_* defines	*/
	int           typesize;			/* size of type in bytes		*/
	unsigned long colcount;			/* count of values in column	*/
	char          *col;				/* heap of data					*/
	ValRecord     min;				/* min value in col				*/
	ValRecord     max;				/* max value in col				*/
	int           sorted;			/* 1 if sorted					*/
	Imprints_index *imps_idx;
} Column;


#endif /* TYPES_H */
