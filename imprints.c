#include "simd_imprints.h"

Imprints_index *
scalar_imprints(Column *column, Imprints_index *imps)
{
	long timer;
	unsigned long i;
	unsigned long colcnt = column->colcount;
	int bit;
	int k, values_per_block = imps->blocksize/column->typesize;

#define GETBIT(Z, X)								\
	do {											\
		int _i;										\
		Z = 0;										\
		for (_i = 0; _i < imps->bins-1; _i++)			\
			Z += (col[i] > imps->bounds[_i].X);	\
	} while (0)

#define SCALAR_IMPS(T, _T, X) {																		\
	T  *restrict col = (T *) column->col;															\
	_T *imprints = (_T *) imps->imprints;													\
	_T mask = 0, prevmask = 0;																		\
	for (i = 0; i < colcnt; i++) {																	\
		mask = 0;																					\
		for (k = 0; k < values_per_block; k++) {													\
			GETBIT(bit, X);																			\
			mask = setBit(mask, bit);																\
			i++;																					\
		}																							\
		if (mask == prevmask && imps->dct[imps->dct_cnt-1].blks < ((1<<MAXOFFSET)-1)) {				\
			if (imps->dct[imps->dct_cnt - 1].repeated == 0) {										\
				if (imps->dct[imps->dct_cnt - 1].blks > 1) {										\
					imps->dct[imps->dct_cnt - 1].blks--; /* reduce previous by 1 */					\
					imps->dct_cnt++;																\
					imps->dct[imps->dct_cnt-1].blks = 1;   /* the new is a repeat */				\
				}																					\
				imps->dct[imps->dct_cnt-1].repeated = 1;											\
			}																						\
			imps->dct[imps->dct_cnt - 1].blks++;													\
		} else {																					\
			/* new mask */																			\
			prevmask = mask;																		\
			imprints[imps->imps_cnt++] = mask;														\
			if (imps->dct_cnt > 0 && imps->dct[imps->dct_cnt - 1].repeated == 0						\
			     && imps->dct[imps->dct_cnt-1].blks < ((1<<MAXOFFSET)-1)) {							\
					imps->dct[imps->dct_cnt - 1].blks++;											\
			} else {																				\
				imps->dct[imps->dct_cnt].blks = 1;													\
				imps->dct[imps->dct_cnt].repeated = 0;												\
				imps->dct_cnt++;																	\
			}																						\
		}																							\
	}																								\
}

#define COLTYPE_SWITCH(_T)												\
	switch (column->coltype) {											\
		case TYPE_bte: SCALAR_IMPS(char, _T, bval); break;				\
		case TYPE_sht: SCALAR_IMPS(short, _T, sval); break;				\
		case TYPE_int: SCALAR_IMPS(int, _T, ival); break;				\
		case TYPE_lng: SCALAR_IMPS(long, _T, lval); break;				\
		case TYPE_oid: SCALAR_IMPS(unsigned long, _T, ulval); break;	\
		case TYPE_flt: SCALAR_IMPS(float, _T, fval); break;				\
		case TYPE_dbl: SCALAR_IMPS(double, _T, dval); break;			\
		default: break;													\
	}

	timer = usec();
	switch (imps->imprintsize) {
		case 1: COLTYPE_SWITCH(unsigned char); break;
		case 2: COLTYPE_SWITCH(unsigned short); break;
		case 4: COLTYPE_SWITCH(unsigned int); break;
		case 8: COLTYPE_SWITCH(unsigned long); break;
	}
	timer = usec() - timer;
	VERBOSE printf("%s imprints scalar creation time=%ld, %ld usec per thousand values\n",
			column->colname, timer, ((long)timer*1000)/column->colcount);

	return imps;
}

/*
Imprints_index *
simd_imprints(Column *column, Imprints_index *imps)
{
	return imps;
}*/

Imprints_index*
create_imprints(Column *column, int blocksize, int max_bins, int simd)
{
	Imprints_index *imps;
	unsigned long max_imprints;

	assert(max_bins <= 256);

	imps              = (Imprints_index *) malloc (sizeof(Imprints_index));
	imps->blocksize   = blocksize; /* in bytes */
	imps->bounds      = (ValRecord *) aligned_alloc(32, sizeof(ValRecord) * max_bins);
	max_imprints      = (column->colcount/(imps->blocksize/column->typesize))+1;
	imps->dct         = (Dct *)  aligned_alloc (32, sizeof(Dct)  * max_imprints);
	imps->dct_cnt     = 0;
	imps->imps_cnt    = 0;

	/* decide the bounds and the number of bins */
	binning(column, imps->bounds, &(imps->bins), max_bins);
	imps->imprintsize = imps->bins/8; /* in bytes */
	imps->imprints    = (char *) aligned_alloc (32, imps->imprintsize * max_imprints);

	if (simd) {
		//return simd_imprints(column, imps);
		return NULL;
	} else {
		if (max_bins > 64 || imps->bins > 64) return NULL;
		return scalar_imprints(column, imps);
	}

	return NULL;
}

#define cmpvalues(X)	\
static int	\
cmpvalues##X(const void *p1, const void *p2)\
{\
	ValRecord v1 = *(ValRecord *) p1;\
	ValRecord v2 = *(ValRecord *) p2;\
	if ( v1.X < v2.X ) return -1; \
	if ( v1.X > v2.X ) return 1; \
	return 0;\
}

cmpvalues(bval);
cmpvalues(sval);
cmpvalues(ival);
cmpvalues(lval);
cmpvalues(ulval);
cmpvalues(fval);
cmpvalues(dval);

void
binning(Column *column, ValRecord *bounds, int *bins, int max_bins) {
	ValRecord *sample;
	unsigned long pos;
	int i, smp;


	sample = (ValRecord *) malloc(sizeof(ValRecord)*SAMPLE_SZ);

#define sampleDistribution(X,T)												\
	/* draw the sample, */													\
	for (i = 0; i < SAMPLE_SZ; i++) {										\
		pos = (rand()*column->colcount)/RAND_MAX;							\
		sample[i].X = ((T*)column->col)[pos];								\
	}																		\
	/* sort the sample and make it unique*/									\
	qsort((char *) sample, SAMPLE_SZ, sizeof(ValRecord), cmpvalues##X);		\
	for (smp = 0, i = 1; i < SAMPLE_SZ; i++) {								\
		if (sample[i].X != sample[smp].X) {									\
			sample[++smp] = sample[i];										\
		}																	\
	}																		\
	smp++;

	switch(column->coltype){
	case TYPE_bte:
		sampleDistribution(bval, char);
		break;
	case TYPE_sht:
		sampleDistribution(sval, short);
		break;
	case TYPE_int:
		sampleDistribution(ival, int);
		break;
	case TYPE_lng:
		sampleDistribution(lval, long);
		break;
	case TYPE_oid:
		sampleDistribution(ulval, unsigned long);
		break;
	case TYPE_flt:
		sampleDistribution(fval, float);
		break;
	case TYPE_dbl:
		sampleDistribution(dval, double);
	}

	bounds[0] = sample[1];
	if (smp < max_bins-1) {
		for (i = 1; i < smp; i++) {
			bounds[i] = sample[i+1];
		}
		if (i<8)              *bins = 8;
		if (8<=i   && i<16)   *bins = 16;
		if (16<=i  && i<32)   *bins = 32;
		if (32<=i  && i<64)   *bins = 64;
		if (64<=i  && i<128)  *bins = 128;
		if (128<=i && i<256)  *bins = 256;
		for (; i < max_bins; i++) {
			bounds[i] = sample[smp-1];
		}
	} else {
		double y, ystep = (double)smp/(double)(max_bins-2);
		*bins = 64;
		for (i = 1, y = ystep; y < smp; y += ystep, i++) {
			bounds[i] = sample[(int)y];
		}
		if (i == max_bins - 2) { /* there is some leftovers */
			assert (y>=smp);
			assert ((y-ystep)<smp);
			bounds[i] = sample[smp-1];
		}
	}

	VERBOSE printf("%s binning gives %d unique values from %d and %d bins\n", column->colname, smp, SAMPLE_SZ, *bins);

	return;
}
