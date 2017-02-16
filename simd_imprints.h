#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#ifndef __APPLE__
#include <malloc.h>
#endif
#include <unistd.h>
#include <limits.h>
#include <float.h>
#include <ammintrin.h>
#include <smmintrin.h>
#include <immintrin.h>

#define BITS      64
#define SIMD_BITS 256
#define MAX_IMPS  50000000
#define SAMPLE_SZ 2048
#define REPETITION  10   /* how many times to repeat each experiment to get an avg */
#define BITRANGE    1    /* 1 = search by bin. 0 = by sampling the data*/

/* printing options */
#define STATS          if (0)
#define VERBOSE        if (1)
#define PRINT_HISTO    if (0)
#define PRINT_IMPRINTS if (0)
#define PRINT_QUERIES  if (0)

#undef COARSE
#ifdef COARSE
#define PAGESIZE sysconf(_SC_PAGESIZE) /* normal page size */
#else
#define PAGESIZE 64
#endif


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
/* size (in bytes) of supported types as defined by MonetDB */
int stride[14] = {0,0,0,1,2,0,4,8,0,0,4,8,8,0};


#define setBit(X,Y)      ((((long)1)<<Y) | ( ~(((long)1)<<Y) & X))
#define isSet(X,Y)       (((((long)1)<<Y) & X) ? 1 : 0)
#define COMPRESSION_MASK (~((~((unsigned long) 0))<<(bins)))
#define getMask(I)       (bins==64?imprints[I]:((((unsigned long) imprints[((I)*bins)/BITS])>>(((I)%(BITS/bins))*(bins))) & COMPRESSION_MASK))

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
	char          filename[1024];
	char          colname[1024];
	char          typename[64];
	int           coltype;			/* as given by TYPE_* defines	*/
	unsigned long colcount;			/* count of values in column	*/
	char          *col;				/* heap of data					*/
	ValRecord     min;				/* min value in col				*/
	ValRecord     max;				/* max value in col				*/
} Column;

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
	unsigned long  dct_cnt;		/* count of dictionary entries			*/
	unsigned long  imps_cnt;	/* count of imprint vector entries		*/
	int            bins;		/* number of bins						*/
	int            blocksize;	/* number of values covered per block	*/
} Imprints_index;

/* Zonemaps */
typedef struct {
	ValRecord min;
	ValRecord max;
} Zonemap;

typedef struct {
	Zonemap *zmaps;
	unsigned long zmaps_cnt;
} Zonemap_index;
