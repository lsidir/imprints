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
imprints_scan(Column *column, Imprints_index *imps, ValRecord low, ValRecord high, long *timer)
{
	unsigned long i, res_cnt = 0;
	unsigned long colcnt = column->colcount;
	unsigned long dcnt, icnt, top_icnt, bcnt, lim;
	int values_per_block = imps->blocksize/column->typesize;

	dcnt = icnt = bcnt = 0;

	#define impsscan(X,T,_T) {							\
		T  *restrict col = (T *) column->col;			\
		_T *restrict imprints = (_T *) imps->imprints;	\
		T l = low.X;									\
		T h = high.X;									\
		_T mask = 0, innermask = 0;								\
		for (i = 0, dcnt = 0; dcnt < imps->dct_cnt; dcnt++) {						\
			if (imps->dct[dcnt].repeated == 0) {									\
				top_icnt = icnt + imps->dct[dcnt].blks;								\
				for (; icnt < top_icnt; bcnt++, icnt++) {							\
					if (imprints[icnt] & mask) {									\
						i = bcnt * values_per_block;								\
						lim = i + values_per_block;									\
						lim = lim > colcnt ? colcnt: lim;							\
						if ((imprints[icnt] & ~innermask) == 0) {					\
							res_cnt += values_per_block;							\
						} else {													\
							for (; i < lim; i++) {									\
								if (col[i] > l && col[i] <= h) {					\
									res_cnt++;										\
								}													\
							}														\
						}															\
					}																\
				}																	\
			} else { /* repeated mask case */										\
				if (imprints[icnt] & mask) {										\
					i = bcnt * values_per_block;									\
					lim = i + values_per_block * imps->dct[dcnt].blks;				\
					lim = lim > colcnt ? colcnt : lim;								\
					if ((imprints[icnt] & ~innermask) == 0) {						\
						res_cnt += values_per_block;								\
					} else {														\
						for (; i < lim; i++) {										\
							if (col[i] > l && col[i] <= h) {						\
								res_cnt++;											\
							}														\
						}															\
					}																\
				}																	\
				bcnt += imps->dct[dcnt].blks;										\
				icnt++;																\
			}																		\
		}																			\
	}

#define COLTYPE_SWITCH(_T)												\
	switch(column->coltype){\
	case TYPE_bte:\
		impsscan(bval, char,_T);\
		break;\
	case TYPE_sht:\
		impsscan(sval, short,_T);\
		break;\
	case TYPE_int:\
		impsscan(ival, int,_T);\
		break;\
	case TYPE_lng:\
		impsscan(lval, long,_T);\
		break;\
	case TYPE_oid:\
		impsscan(ulval, unsigned long,_T);\
		break;\
	case TYPE_flt:\
		impsscan(fval, float,_T);\
		break;\
	case TYPE_dbl:\
		impsscan(dval, double,_T);\
		break;\
	}

	*timer = usec();
	switch (imps->imprintsize) {
		case 1: COLTYPE_SWITCH(unsigned char); break;
		case 2: COLTYPE_SWITCH(unsigned short); break;
		case 4: COLTYPE_SWITCH(unsigned int); break;
		case 8: COLTYPE_SWITCH(unsigned long); break;
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
