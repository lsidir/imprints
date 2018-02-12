/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

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

#define BITS 64

//#define COARSE
#ifdef COARSE
#define PAGESIZE sysconf(_SC_PAGESIZE) /*normal page size */
#else
#define PAGESIZE 64
#endif

#define COMPRESSLIST

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

/* the target column being analysed */
extern char filename[1024];
extern char column[1024];
extern char *col;
extern int coltype;
extern unsigned long colcount;
extern unsigned long pages;

/* the zonemaps are our frame of reference */
typedef struct{
	ValRecord min, max; /* zone marks */
} Zonemap;

/* The ColumnFingerprint contains a binned bit-mask structure
 * to weed out lines/blocks of no interest in a scan.
 * The mask vector is simply compressed, keeping track on
 * the number of blocks it covers. From the blks we can
 * calculate the actual oid ranges, provided we scan only.
*/
#define MAXOFFSET 24
typedef struct {
	unsigned int blks:MAXOFFSET;
	unsigned int repeated:1; /* the finger prints in the range are all the same */
	unsigned int flgs:8 * sizeof(int) - MAXOFFSET-1; /* for future use, e.g. micro swaps */
} Fingerprint;


/* global bounds */
extern ValRecord min, max, slice; /* 62 bins equidistant in value distribution */
extern ValRecord absmin, absmax, mxbins[64],mibins[64];
extern int equidistance;

#define MAXCFP 50000000
extern Zonemap *zmap;
extern long zonetop;

extern Fingerprint *cfp;
extern long *bitmask;
extern int cfptop;

/* each column fingerprint is characterized by a partition vector */

extern long histogram[BITS]; /* of bin fillers */
/* vector filling distribution */
extern long vectors[BITS+1];

extern int rpp; /* rows per page */

extern long qhits[101];

#define setBit(X,Y)  ((((long)1)<<Y) | ( ~(((long)1)<<Y) & X))
#define isSet(X,Y)   (((((long)1)<<Y) & X) ? 1 : 0)
#define COMPRESSION_MASK (~((~((unsigned long) 0))<<(BINS)))
#define getMask(I)   (BINS==64?bitmask[I]:((((unsigned long) bitmask[((I)*BINS)/BITS])>>(((I)%(BITS/BINS))*(BINS))) & COMPRESSION_MASK))

/* what type of histogram to build */
enum {EQWIDTH,EQHEIGHT};
