all: simd_imprints

simd_imprints: simd_imprints.c helpers.c simd_imprints.h
	clang -O3 -g -Wall -march=native helpers.c simd_imprints.c -lm -o simd_imprints

clean: simd_imprints
	rm simd_imprints
