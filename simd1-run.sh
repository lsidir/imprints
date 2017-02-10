#!/bin/sh
clang simd1.c -march=native -Wall -O3 -g && ./a.out 3000000 256

