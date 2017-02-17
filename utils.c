#include "simd_imprints.h"

/* check if a column is sorted */
void isSorted(Column *column)
{
	int i;

#define checksorted(T)													\
	for (column->sorted = 1, i = 1; i < column->colcount; i++)			\
		if (((T*)column->col)[i] < ((T*)column->col)[i-1]) {				\
			column->sorted = 0;											\
			break;														\
		}																\
	if (column->sorted == 0)												\
		for (column->sorted=1, i = 1; i<column->colcount; i++) {			\
			if (((T*)column->col)[i] > ((T*)column->col)[i-1]) {			\
				column->sorted = 0;										\
				break;													\
			}															\
		}

	switch (column->coltype) {
	case TYPE_bte:
		checksorted(char);
		break;
	case TYPE_sht:
		checksorted(short);
		break;
	case TYPE_int:
		checksorted(int);
		break;
	case TYPE_lng:
		checksorted(long);
		break;
	case TYPE_oid:
		checksorted(unsigned long);
		break;
	case TYPE_flt:
		checksorted(float);
		break;
	case TYPE_dbl:
		checksorted(double);
	}

	VERBOSE printf("%s sorted property %d\n", column->colname, column->sorted);
}

long usec()
{
	static struct timeval tpbase; /* automatically initialized to 0 */
	struct timeval tp;

	if (tpbase.tv_sec == 0)
		gettimeofday(&tpbase, 0);
	gettimeofday(&tp, 0);
	tp.tv_sec -= tpbase.tv_sec;
	return (long) tp.tv_sec * 1000000 + (long) tp.tv_usec;
}
