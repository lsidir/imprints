all: simd_imprints

simd_imprints: simd_imprints.c simd_imprints.h
	gcc -g simd_imprints.c -lm -o simd_imprints

fast:
	gcc -O3 simd_imprints.c -lm -o simd_imprints

clean: simd_imprints
	rm simd_imprints
