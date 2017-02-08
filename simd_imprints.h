#include <stdio.h>
#include <assert.h>
#include <string.h>
#include <stdlib.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>
#include <malloc.h>
#include <unistd.h>
#include <limits.h>
#include <float.h>

#define BITS      64
#define MAX_IMPS  50000000
#define SAMPLE_SZ 2048
#define REPETITION  10   /* how many times to repeat each experiment to get an avg */



#define BITRANGE    1    /* search by bin 1 or by sampling the data 0 */

/* printing options */
#define STATS if (1)
#define VERBOSE if (1)

#undef COARSE
#ifdef COARSE
#define PAGESIZE sysconf(_SC_PAGESIZE) /* normal page size */
#else
#define PAGESIZE 64
#endif


/* just to get a BATlike structure */
#define TYPE_void   0
#define TYPE_bte    3
#define TYPE_sht    4
#define TYPE_int    6
#define TYPE_oid    7
#define TYPE_flt    10
#define TYPE_dbl    11
#define TYPE_lng    12
#define TYPE_str    13

#define setBit(X,Y)      ((((long)1)<<Y) | ( ~(((long)1)<<Y) & X))
#define isSet(X,Y)       (((((long)1)<<Y) & X) ? 1 : 0)
#define COMPRESSION_MASK (~((~((unsigned long) 0))<<(BINS)))
#define getMask(I)       (BINS==64?bitmask[I]:((((unsigned long) bitmask[((I)*BINS)/BITS])>>(((I)%(BITS/BINS))*(BINS))) & COMPRESSION_MASK))

extern int stride[14];

typedef union {
	char bval;
	short sval;
	int ival;
	long lval;
	unsigned long ulval;
	float fval;
	double dval;
} ValRecord;

/* The Column Imprints contains a binned bit-mask structure
 * to weed out lines/blocks of no interest in a scan.
 * The mask vector is simply compressed, keeping track on
 * the number of blocks it covers. From the blks we can
 * calculate the actual oid ranges, provided we scan only.
*/
#define MAXOFFSET 24
typedef struct {
	unsigned int blks:MAXOFFSET;
	unsigned int repeated:1; /* the imprints in the range are all the same */
	unsigned int flgs:8 * sizeof(int) - MAXOFFSET-1; /* for future use, e.g. micro swaps */
} Imprint;


/* global bounds */
extern int equidistance;

extern Imprint *imprint;
extern long *bitmask;
extern int imptop;

/* each column imprint is characterized by a partition vector */

extern long histogram[BITS]; /* of bin fillers */
/* vector filling distribution */
extern long vectors[BITS+1];

extern long qhits[101];
