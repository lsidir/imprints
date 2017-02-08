all: leanfingerprint

leanfingerprint: leanfingerprint.c leanfingerprint.h
	gcc -g leanfingerprint.c -lm -o leanfingerprint

fast:
	gcc -O3 leanfingerprint.c -lm -o leanfingerprint

clean: leanfingerprint
	rm leanfingerprint
