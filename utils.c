/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "simd_imprints.h"

/* check if a column is sorted */
void isSorted(Column *column)
{
	int i;

#define checksorted(T)													\
	for (column->sorted = 1, i = 1; i < column->colcount; i++)			\
		if (((T*)column->col)[i] < ((T*)column->col)[i-1]) {				\
			column->sorted = 0;											\
			break;														\
		}																\
	if (column->sorted == 0)												\
		for (column->sorted=1, i = 1; i<column->colcount; i++) {			\
			if (((T*)column->col)[i] > ((T*)column->col)[i-1]) {			\
				column->sorted = 0;										\
				break;													\
			}															\
		}

	switch (column->coltype) {
	case TYPE_bte:
		checksorted(char);
		break;
	case TYPE_sht:
		checksorted(short);
		break;
	case TYPE_int:
		checksorted(int);
		break;
	case TYPE_lng:
		checksorted(long);
		break;
	case TYPE_oid:
		checksorted(unsigned long);
		break;
	case TYPE_flt:
		checksorted(float);
		break;
	case TYPE_dbl:
		checksorted(double);
	}

	VERBOSE printf("%s sorted property %d\n", column->colname, column->sorted);
}

long usec()
{
	static struct timeval tpbase; /* automatically initialized to 0 */
	struct timeval tp;

	if (tpbase.tv_sec == 0)
		gettimeofday(&tpbase, 0);
	gettimeofday(&tp, 0);
	tp.tv_sec -= tpbase.tv_sec;
	return (long) tp.tv_sec * 1000000 + (long) tp.tv_usec;
}

__m256i setbit_256(__m256i x, int k) {
	// constants that will (hopefully) be hoisted out of a loop after inlining
	__m256i indices = _mm256_set_epi32(224,192,160,128,96,64,32,0);
	__m256i one = _mm256_set1_epi32(-1);
	one = _mm256_srli_epi32(one, 31);    // set1(0x1)
	__m256i kvec = _mm256_set1_epi32(k);
	// if 0<=k<=255 then kvec-indices has exactly one element with a value between 0 and 31
	__m256i shiftcounts = _mm256_sub_epi32(kvec, indices);
	__m256i kbit        = _mm256_sllv_epi32(one, shiftcounts);   // shift counts outside 0..31 shift the bit out of the element
	return _mm256_or_si256(kbit, x);                             // use _mm256_andnot_si256 to unset the k-th bit
}

