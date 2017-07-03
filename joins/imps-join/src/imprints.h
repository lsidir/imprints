/*
 * imprints.h
 *
 * alternative methods to build imprints for both columns
 *  Created on: 29 May 2017
 *      Author: zeroxwj
 */

#ifndef IMPRINTS_H_
#define IMPRINTS_H_

#include "types.h"
#include "common.h"

#define BITS      64
#define SIMD_BITS 256
#define MAX_IMPS  50000000
#define SAMPLE_SZ 2048
#define REPETITION  10   /* how many times to repeat each experiment to get an avg */
#define BITRANGE    1    /* 1 = search by bin. 0 = by sampling the data*/

#define setBit(X,Y)      ((X) | ((unsigned long long) 1 << (Y)))
#define isSet(X,Y)       (((((unsigned long long)1)<<Y) & X) ? 1 : 0)
#define COMPRESSION_MASK (~((~((unsigned long) 0))<<(imps->bins)))
#define getMask(I)       (imps->bins==64?imps->imprints[I]:((((unsigned long) imps->imprints[((I)*imps->bins)/64])>>(((I)%(64/imps->bins))*(imps->bins))) & COMPRESSION_MASK))


/**
 * binning (build) column R, build imprints for (probe) column L according to R's bin borders
 * @return 0  succeed to build
 * @return -1 fail to build
 */
int buildImpsR2L(Column *coLL, Column *colR, int max_bins, int blocksize);

int buildImpsL2R(Column *colL, Column *colR, int max_bins, int blocksize);

#endif /* IMPRINTS_H_ */
