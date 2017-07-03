
#ifndef UTILS_H_
#define UTILS_H_

#include "types.h"
#include "common.h"




void isSorted(Column *column);
long usec();
__m256i setbit_256(__m256i x, int k);
uint64_t unsignedLog2 (uint64_t val);
uint64_t nextPowerOfTwo(uint64_t n);

#endif /* UTILS_H_ */
