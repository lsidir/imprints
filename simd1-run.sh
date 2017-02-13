#!/bin/sh
clang simd1.c -march=native -Wall  -g -O3 && ./a.out 1000000 256
