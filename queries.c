#include "simd_imprints.h"

unsigned long
simple_scan(Column *column, ValRecord low, ValRecord high, long *timer)
{
	unsigned long i, res_cnt = 0;
	unsigned long colcnt = column->colcount;

	#define simplescan(X,T) {					\
		T  *restrict col = (T *) column->col;	\
		T l = low.X;							\
		T h = high.X;							\
		for (i = 0; i < colcnt; i++) {			\
			if (col[i] > l && col[i] <= h) {	\
				res_cnt++;						\
			}									\
		}										\
	}

	*timer = usec();
	switch(column->coltype){
	case TYPE_bte:
		simplescan(bval, char);
		break;
	case TYPE_sht:
		simplescan(sval, short);
		break;
	case TYPE_int:
		simplescan(ival, int);
		break;
	case TYPE_lng:
		simplescan(lval, long);
		break;
	case TYPE_oid:
		simplescan(ulval, unsigned long);
		break;
	case TYPE_flt:
		simplescan(fval, float);
		break;
	case TYPE_dbl:
		simplescan(dval, double);
		break;
	}
	*timer = usec() - *timer;
	return res_cnt;
}

unsigned long
imprints_scan(Column *column, Imprints *imps, ValRecord low, ValRecord high, long *timer)
{
	unsigned long i, res_cnt = 0;
	unsigned long colcnt = column->colcount;


	#define impsscan(X,T) {					\
		T  *restrict col = (T *) column->col;	\
		T l = low.X;							\
		T h = high.X;							\
		for (i = 0; i < colcnt; i++) {			\
			if (col[i] > l && col[i] <= h) {	\
				res_cnt++;						\
			}									\
		}										\
		\
		for (dcnt = 0; dcnt < dct_cnt; dcnt++) {									\
			if (imps->dct[dcnt].repeated == 0) {									\
				for (k = tf + imps->dct[dcnt].blks; tf < k; n++, tf++) {			\
						if (imprints##B[tf] & mask) {							\
							register T val;										\
							l = n * rpp;										\
							lim = l + rpp;										\
							lim = lim > column->colcount ? column->colcount: lim;				\
							if ((imprints##B[tf] & ~innermask) == 0) {			\
								for (; l < lim; l++) {							\
									oids[oid++] = l;							\
								}												\
							} else {											\
								for (val = ((T*)column->col)[l]; l < lim; l++, val = ((T*)column->col)[l]) {	\
									if (val <= shigh.X && val > slow.X) {		\
										oids[oid++] = l;						\
									}											\
								}												\
							}													\
						}														\
					}															\
			} else { /* repeated mask case */									\
				if (imprints##B[tf] & mask) {								\
						register T val;											\
						l = n * rpp;											\
						lim = l + rpp * imps->dct[dcnt].blks;						\
						lim = lim > column->colcount ? column->colcount : lim;					\
						if ((imprints##B[tf] & ~innermask) == 0) {				\
							for (; l < lim; l++) {								\
								oids[oid++] = l;								\
							}													\
						} else {												\
							for (val = ((T*)column->col)[l]; l < lim; l++, val = ((T*)column->col)[l]) {	\
								if (val <= shigh.X && val > slow.X  ) {			\
									oids[oid++] = l;							\
								}												\
							}													\
						}														\
					}															\
					n += imps->dct[dcnt].blks;										\
					tf++;														\
				}\
			}
	}


	*timer = usec();
	switch(column->coltype){
	case TYPE_bte:
		impsscan(bval, char);
		break;
	case TYPE_sht:
		impsscan(sval, short);
		break;
	case TYPE_int:
		impsscan(ival, int);
		break;
	case TYPE_lng:
		impsscan(lval, long);
		break;
	case TYPE_oid:
		impsscan(ulval, unsigned long);
		break;
	case TYPE_flt:
		impsscan(fval, float);
		break;
	case TYPE_dbl:
		impsscan(dval, double);
		break;
	}
	*timer = usec() - *timer;
	return res_cnt;
}


#if 0
/* zone maps for the future to fix */
		/* zonesmaps scan */
		m   = 0;
		oid = 0;
		/* watch out, the zones are closed (min,max) intervals */
		#define zonequery(X,T) \
			zonetime = usec();					\
			for (j = 0; j < zonemaps->zmaps_cnt; j++) {		\
				STATS zindex[i] += 1;			\
				if (slow.X < zonemaps->zmaps[j].min.X && shigh.X >= zonemaps->zmaps[j].max.X) { /* all qualify */			\
					for (k = j * rpp, lim = k + rpp < column->colcount ? k+rpp : column->colcount; k < lim; k++) {	\
						oids[oid++] = k;															\
					}																				\
				} else if (!(shigh.X < zonemaps->zmaps[j].min.X || slow.X >= zonemaps->zmaps[j].max.X)) {					\
					/* zone maps are inclusive */													\
					for (k = j * rpp, lim = k + rpp < column->colcount ? k+rpp : column->colcount; k < lim; k++) {	\
						STATS zcomparisons[i] += 1;													\
						if (((T*)column->col)[k] <= shigh.X && ((T*)column->col)[k] > slow.X ) {					\
							oids[oid++] = k;														\
						}\
					} \
				}\
			}

		switch(column->coltype){
		case TYPE_bte:
			zonequery(bval,char);
			break;
		case TYPE_sht:
			zonequery(sval,short);
			break;
		case TYPE_int:
			zonequery(ival,int);
			break;
		case TYPE_lng:
			zonequery(lval,long);
			break;
		case TYPE_oid:
			zonequery(ulval,unsigned long);
			break;
		case TYPE_flt:
			zonequery(fval,float);
			break;
		case TYPE_dbl:
			zonequery(dval,double);
		}

		zonetimer[i] = usec() - zonetime;
		if (tuples[i] != oid) /* correct result check */
			printf("%s base %ld zonemap %ld differ\n", column->colname, tuples[i], oid);
		fprintf(devnull," %ld",m); /* to break compiler optimizations */
#endif
