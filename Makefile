all: simd_imprints

simd_imprints: main.c utils.c queries.c print.c imprints.c zonemaps.c main.h
	clang -O3 -g -Wall -march=native utils.c queries.c print.c imprints.c main.c zonemaps.c -lm -o simd_imprints

clean:
	rm simd_imprints
