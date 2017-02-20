all: simd_imprints

simd_imprints: simd_imprints.c utils.c queries.c print.c imprints.c simd_imprints.h
	clang -g -Wall -march=native utils.c queries.c print.c imprints.c simd_imprints.c -lm -o simd_imprints

clean:
	rm simd_imprints
