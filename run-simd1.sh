#!/bin/bash
N=1000000
clang simd1.c -march=native -Wall -O3 -g

bits=10
while (( bits < 257 )); do
	echo -e $bits "\t" `./a.out $N $bits`
    (( bits += 10 ))
done
