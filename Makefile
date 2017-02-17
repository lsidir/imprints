all: simd_imprints

simd_imprints: simd_imprints.c utils.c print.c simd_imprints.h
	clang -O3 -g -Wall -march=native utils.c print.c simd_imprints.c -lm -o simd_imprints

clean:
	rm simd_imprints
