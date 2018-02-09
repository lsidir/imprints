/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "simd_imprints.h"

Imprints_index *
scalar_imprints(Column *column, Imprints_index *imps)
{
	long timer;
	unsigned long i, loops = 0, unique_imps = 0, compres_imps = 0;
	unsigned long colcnt = column->colcount;
	int bit;
	int k, values_per_block = imps->blocksize/column->typesize;

#define GETBIT(Z, X)								\
	do {											\
		int _i;										\
		Z = 0;										\
		for (_i = 0; _i < imps->bins-1; _i++)		\
			Z += (col[i] > imps->bounds[_i].X);		\
	} while (0)

#define SCALAR_IMPS(T, _T, X) {																		\
	T  *restrict col = (T *) column->col;															\
	_T *restrict imprints = (_T *) imps->imprints;													\
	_T mask = 0, prevmask = 0;																		\
	for (i = 0; i < colcnt;) {																		\
		mask = 0;																					\
		for (k = 0; k < values_per_block; k++) {													\
			GETBIT(bit, X);																			\
			mask = setBit(mask, bit);																\
			i++;																					\
		}																							\
		loops++;\
		if (mask == prevmask && imps->dct[imps->dct_cnt-1].blks < ((1<<MAXOFFSET)-1)) {				\
			compres_imps++;\
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
			unique_imps++;\
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
		case 1: COLTYPE_SWITCH(uint8_t); break;
		case 2: COLTYPE_SWITCH(uint16_t); break;
		case 4: COLTYPE_SWITCH(uint32_t); break;
		case 8: COLTYPE_SWITCH(uint64_t); break;
	}
	timer = usec() - timer;
	VERBOSE printf("%s imprints scalar creation time=%ld, %ld usec per thousand values, [loops=%lu unique=%lu compressed=%lu]\n",
			column->colname, timer, ((long)timer*1000)/column->colcount, loops, unique_imps, compres_imps);

	return imps;
}

Imprints_index *
simd_imprints(Column *column, Imprints_index *imps)
{
	long          timer;
	unsigned long i, loops = 0, unique_imps = 0, compres_imps = 0;
	unsigned long colcnt              = column->colcount;
	unsigned long e, k;
	int           values_per_block    = imps->blocksize/column->typesize;
	int           values_per_simd     = 32/column->typesize;
	int           simds_per_block     = values_per_block/values_per_simd;
	char *restrict imprints           = imps->imprints;
	/* simd stuff */
	//__m256i       *simd_imprints   = (__m256i *) imps->imprints;
	__m256i       zero             = _mm256_setzero_si256();
	__m256i       *restrict limits = aligned_alloc(32, (imps->bins) * sizeof(__m256i));
	__m256i       bitmasks[256];

	for (i = 0; i < 256; i++) {
		/* setbit_256(x, 0) will set the first bit, and so on */
		bitmasks[i] = setbit_256(zero, i);
	}

	#define MAKE_LIMITS(SIMDTYPE, X)											\
	for (int _i = 0; _i < imps->bins; _i++)										\
		limits[_i] = _mm256_set1_##SIMDTYPE(imps->bounds[_i].X);				\

	switch (column->coltype) {
		case TYPE_bte: MAKE_LIMITS(epi8, bval); break;
		case TYPE_sht: MAKE_LIMITS(epi16, sval); break;
		case TYPE_int: MAKE_LIMITS(epi32, ival); break;
		case TYPE_lng: MAKE_LIMITS(epi64x, lval); break;
		case TYPE_oid: MAKE_LIMITS(epi64x, ulval); break;
		case TYPE_flt: MAKE_LIMITS(ps, fval); break;
		case TYPE_dbl: MAKE_LIMITS(pd, dval); break;
		default: return NULL;
	}

	#define GETBIT_SIMD(SIMDTYPE)															\
		/* perform 2 bin comparisons per simd instruction until all bins are checked */		\
		for (int bin1 = 0, bin2 = 1; bin1 < imps->bins-2; bin1+=2, bin2+=2) {				\
			result = _mm256_add_##SIMDTYPE(result,											\
						_mm256_add_##SIMDTYPE(												\
							_mm256_cmpgt_##SIMDTYPE(values_v, limits[bin1]),				\
							_mm256_cmpgt_##SIMDTYPE(values_v, limits[bin2])));				\
		}																					\
		result = _mm256_sub_##SIMDTYPE(zero,result);										\
		for (int i1 = 0, i2=1; i1 < values_per_simd-1; i1+=2, i2+=2) {						\
			simd_mask = _mm256_or_si256(													\
				_mm256_or_si256(simd_mask, bitmasks[_mm256_extract_##SIMDTYPE(result, i1)]),\
				_mm256_or_si256(simd_mask, bitmasks[_mm256_extract_##SIMDTYPE(result, i2)])	\
			);																				\
		}

	#define GETBIT_SIMDD(SIMDTYPE) return NULL;
	#define GETBIT_SIMDF(SIMDTYPE) return NULL;

	#define SIMD_IMPS(T, X) {																			\
		T  *restrict col      = (T *) column->col;														\
		__m256i simd_mask     = zero;																	\
		__m256i simd_prevmask = zero;																	\
		__m256i check;																					\
		for (i = 0; i < colcnt;) {																		\
			simd_mask = zero;																			\
			for (int spb = 0; spb < simds_per_block && i < colcnt; spb++) {								\
				__m256i values_v = _mm256_load_si256((__m256i*) (col+i));								\
				__m256i result = zero;																	\
				GETBIT_SIMD(X);																			\
				i += values_per_simd;																	\
			}																							\
			loops++;\
			check = _mm256_xor_si256(simd_mask,simd_prevmask);											\
			if (_mm256_testz_si256(check,check) && imps->dct[imps->dct_cnt-1].blks < ((1<<MAXOFFSET)-1)) {\
				compres_imps++;\
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
				unique_imps++;\
				simd_prevmask = simd_mask;																\
				for (e = 0, k = imps->imps_cnt*imps->imprintsize; e < imps->imprintsize; e++) {			\
					imprints[k+e] = _mm256_extract_epi8(simd_mask, e);							\
				}																						\
				imps->imps_cnt++;																		\
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

	timer = usec();
	switch (column->coltype) {
		case TYPE_bte: SIMD_IMPS(char, epi8); break;
		case TYPE_sht: SIMD_IMPS(short, epi16); break;
		case TYPE_int: SIMD_IMPS(int, epi32); break;
		case TYPE_lng: SIMD_IMPS(long, epi64); break;
		case TYPE_oid: SIMD_IMPS(unsigned long, epi64); break;
		case TYPE_flt: break;
		case TYPE_dbl: break;
		default: break;
	}
	timer = usec() - timer;
	VERBOSE printf("%s imprints simd creation time=%ld, %ld usec per thousand values [loops=%lu unique=%lu compressed=%lu]\n",
			column->colname, timer, ((long)timer*1000)/column->colcount, loops, unique_imps, compres_imps);
	return imps;
}

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
	imps->imprints    = (char *) aligned_alloc (32, imps->imprintsize * (max_imprints+32));

	if (simd) {
		imps = simd_imprints(column, imps);
	} else {
		if (max_bins > 64 || imps->bins > 64) {
			printf("%s bin size larger than 64 for scalar version is not permitted\n", column->colname);
			return NULL;
		}
		imps = scalar_imprints(column, imps);
	}

	VERBOSE printf("%s %s imprints data_size=%ld(bytes) typesize=%d imprints_size=%ld(bytes)[%ld%%] imps_cnt=%lu dct_cnt=%lu bins=%d blocksize=%d(bytes) imprintsize=%d values_per_bloc=%d\n",
				   column->colname,
				   simd?"SIMD  ":"SCALAR",
				   column->colcount*column->typesize,
				   column->typesize,
				   imps->imps_cnt*imps->imprintsize + imps->dct_cnt*sizeof(Dct),
				   ((imps->imps_cnt*imps->imprintsize + imps->dct_cnt*sizeof(Dct))*100)/(column->colcount*column->typesize),
				   imps->imps_cnt,
				   imps->dct_cnt,
				   imps->bins,
				   imps->blocksize,
				   imps->imprintsize,
				   imps->blocksize/column->typesize
				   );

	return imps;
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
	ValRecord *sample, max;
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
		max.bval = 127;
		break;
	case TYPE_sht:
		sampleDistribution(sval, short);
		max.sval = 32767;
		break;
	case TYPE_int:
		sampleDistribution(ival, int);
		max.ival = INT_MAX;
		break;
	case TYPE_lng:
		sampleDistribution(lval, long);
		max.lval = LONG_MAX;
		break;
	case TYPE_oid:
		sampleDistribution(ulval, unsigned long);
		max.ulval = ULONG_MAX;
		break;
	case TYPE_flt:
		sampleDistribution(fval, float);
		max.fval = FLT_MAX;
		break;
	case TYPE_dbl:
		sampleDistribution(dval, double);
		max.dval = DBL_MAX;
	}

	printf("%s new binning with max_bins = %d\n", column->colname, max_bins);

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
			bounds[i] = max;
		}
	} else {
		double y, ystep = (double)smp/(double)(max_bins-2);
		*bins = max_bins;
		for (i = 1, y = ystep; y < smp; y += ystep, i++) {
			bounds[i] = sample[(int)y];
		}
		if (i == max_bins - 2) { /* there is some leftovers */
			assert (y>=smp);
			assert ((y-ystep)<smp);
			bounds[i++] = sample[smp-1];
		}
		for (; i < max_bins; i++) {
			bounds[i] = max;
		}
	}

	VERBOSE printf("%s binning gives %d unique values from %d and %d bins\n", column->colname, smp, SAMPLE_SZ, *bins);

	return;
}
