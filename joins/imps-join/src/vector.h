/*
 * vector.h
 *
 *  Created on: 28 Jul 2017
 *      Author: zeroxwj
 */

#ifndef VECTOR_H_
#define VECTOR_H_

typedef struct {
	uint32_t *array;
	size_t used;
	size_t size;
} Array;

void initArray(Array *a, size_t initialSize) {
	a->array = 0;
	a->array = (uint32_t *) malloc(initialSize * sizeof(uint32_t));

	if (a->array == 0) {
		printf("[Error] malloc failed in initArray()\n");
		exit(EXIT_FAILURE);
	}

	a->used = 0;
	a->size = initialSize;
}

void insertArray(Array *a, uint32_t element) {
	// a->used is the number of used entries, because a->array[a->used++] updates a->used only *after* the array has been accessed.
	// Therefore a->used can go up to a->size
	if (a->used == a->size) {
		a->size *= 2;
		a->array = (uint32_t *) realloc(a->array, a->size * sizeof(uint32_t));

		if (a->array == 0) {
			printf("[Error] realloc failed in insertArray()\n");
			exit(EXIT_FAILURE);
		}
	}
	a->array[a->used++] = element;
}

void freeArray(Array *a) {
	free(a->array);
	a->array = 0;
	a->used = a->size = 0;
}

#endif /* VECTOR_H_ */
