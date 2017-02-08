#include <ammintrin.h>
#include <smmintrin.h>
#include <immintrin.h>


#include <stdio.h>

char p[5] = {1, 1 << 1, 1 << 2, 1 << 3, 1 << 4};

void dump(__m256i val) {
		for (int i = 0; i < 8; i++) {
		printf("%i ", _mm256_extract_epi32(val, i));
	}
	printf("\n");
}

void printBits(size_t const size, void const * const ptr)
{
    unsigned char *b = (unsigned char*) ptr;
    unsigned char byte;
    int i, j;

    for (i=size-1;i>=0;i--)
    {
        for (j=7;j>=0;j--)
        {
            byte = b[i] & (1<<j);
            byte >>= j;
            printf("%u", byte);
        }
    }
    puts("");
}

inline __m256i setbit_256(__m256i x,int k){
// constants that will (hopefully) be hoisted out of a loop after inlining  
  __m256i indices = _mm256_set_epi32(224,192,160,128,96,64,32,0);
  __m256i one = _mm256_set1_epi32(-1);
  one = _mm256_srli_epi32(one, 31);    // set1(0x1)


  __m256i kvec = _mm256_set1_epi32(k);  
// if 0<=k<=255 then kvec-indices has exactly one element with a value between 0 and 31
  __m256i shiftcounts = _mm256_sub_epi32(kvec, indices);
  __m256i kbit        = _mm256_sllv_epi32(one, shiftcounts);   // shift counts outside 0..31 shift the bit out of the element
                                                               // kth bit set, all 255 other bits zero.
  return _mm256_or_si256(kbit, x);                             // use _mm256_andnot_si256 to unset the k-th bit
}


int main(int argc,char** argv) {
	int values[8] = {-2, 10, 20, 50, 100, 42, 140, 200};
	int limits[4] = {0, 5, 100, 150};

	__m256i zero        = _mm256_setzero_si256();
	__m256i one         = _mm256_set1_epi32(1);

	__m256i values_v    = _mm256_load_si256((__m256i*) values);
	__m256i result      = _mm256_setzero_si256();

	__m256i bitmasks[256];
	for (int i = 0; i < 256; i++) {
		bitmasks[i] = setbit_256(zero, i);
	}

	for (int l = 0; l < 4; l++) {
		__m256i limit1  = _mm256_set1_epi32 (limits[l]);
		__m256i result1 = _mm256_cmpgt_epi32(values_v, limit1);
		__m256i result2 = _mm256_and_si256(result1, one);
		result          = _mm256_add_epi32  (result, result2);
	}

	dump(result);

	__m256i imprint = _mm256_setzero_si256();

	for (int i = 0; i < 8; i++) {
		imprint = _mm256_or_si256(imprint, bitmasks[_mm256_extract_epi32(result, i)]);

	}
	printBits(sizeof(__m256i), &imprint);

	return 0;
}
