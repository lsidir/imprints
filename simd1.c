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
            printf("%c", byte == 0 ? '_' : 'X');
        }
    }
    puts("");
}

__m256i setbit_256(__m256i x,int k){
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


#define IMPRINT_BITS 256
#define VALUES_PER_IMPRINT 256

#define VALUE_BITS 32


int main(int argc,char** argv) {
	int n = IMPRINT_BITS * 1000; // nicely divisible by 8 and 256

	int* values       = malloc(sizeof(int) * n);
	char* value_ptr   = (char*) values;
	int* limits       = malloc(sizeof(int) * IMPRINT_BITS);
	__m256i* imprints = calloc(sizeof(__m256i) * n/IMPRINT_BITS, 1);

	if (!values || !limits || !imprints) {
		return -1;
	}
	// generate some values
	for (int i = 0; i < n; i++) {
		values[i] = i; // best case for imprints, yes?
	}
	// some even-spaced limits
	for (int i = 0; i < IMPRINT_BITS; i++) {
		limits[i] = i* n/IMPRINT_BITS;
	}

	__m256i zero        = _mm256_setzero_si256();
	__m256i one         = _mm256_set1_epi32(1);

	__m256i bitmasks[256];
	for (int i = 0; i < 256; i++) {
		bitmasks[i] = setbit_256(zero, i);
	}

	for (int chunk = 0; chunk < n/VALUES_PER_IMPRINT; chunk++) {

		for (int chunk2 = 0; chunk2 < 32; chunk2++){

			__m256i values_v    = _mm256_load_si256((__m256i*) value_ptr);
			__m256i result      = _mm256_setzero_si256();

			// TODO: Do this in larger groups
			for (int l = 0; l < IMPRINT_BITS - 1; l++) {
				__m256i limit1  = _mm256_set1_epi32(limits[l]);
				__m256i result1 = _mm256_cmpgt_epi32(values_v, limit1);
				// turn -1 from cmpgt into 1
				result          = _mm256_add_epi32(result, result1);
			}
			result = _mm256_abs_epi32(result);
			
			// todo: can we do better here going from indices to bit patterns?
			for (int i = 0; i < IMPRINT_BITS/VALUE_BITS; i++) {
				imprints[chunk] = _mm256_or_si256(imprints[chunk], bitmasks[_mm256_extract_epi32(result, i)]);
			}
			value_ptr += sizeof(__m256i);

		}
		printBits(sizeof(__m256i), &imprints[chunk]);
	}


	return 0;
}
