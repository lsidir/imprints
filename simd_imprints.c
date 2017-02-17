#include "simd_imprints.h"

#ifdef __APPLE__
#define aligned_alloc(align, n) malloc(n)
#endif

/* auxilary variables used globaly */
FILE *devnull;
unsigned long pages;                            /* total pages in the column */
int rpp;                                        /* rows per page */

/* timing variables */
long zone_create_time;
long imprints_create_time;
long simd_imprints_create_time;

/*****************/
long globalmask;
Dct *dct_scalar;
long *imprints;
long imps_cnt;
long dct_cnt;

/* stat vars for imprints */
long histogram[64]; /* bin filling */
long vectors[65]; /* vector filling distribution */


/* global vars for queries */
ValRecord slow, shigh;
unsigned long mask;
unsigned long innermask;

/* functions (in call order)-ish */
Zonemap_index *create_zonemaps(Column *column);
Imprints_index *scalar_imprints1(Column *column);
Imprints_index *simd_imprints(Column *column, int blocksize, Imprints_index *other);

void stats(Column *column, Imprints_index *);
void queries(Column *column, Zonemap_index *zonemaps, Imprints_index *imps);
void simd_queries(Column *column, Imprints_index *imps, long results);
void genQueryRange(Column *column, Imprints_index *imps, int i);

/* main function for stand alone imprints */
int main(int argc, char **argv)
{
	Column *column;
	FILE *cfile;
	long filesize;
	size_t rd;
	int stride[14] = {0,0,0,1,2,0,4,8,0,0,4,8,8,0};
	Zonemap_index *zonemaps;
	Imprints_index *scalar_imps;
	//Imprints_index *simd_imps;
	Imprints_index *imps2;

	if (argc != 5) {
		printf("usage: %s type count file column\n", argv[0]);
		return -1;
	}

	column = (Column *) malloc(sizeof(Column));

	strcpy(column->colname, argv[4]);
	strcpy(column->filename, argv[3]);
	column->colcount = atoi(argv[2]);
	strcpy(column->typename, argv[1]);

	if (strcmp(column->typename, "tinyint") == 0 || strcmp(argv[1], "boolean") == 0) {
		column->coltype  = TYPE_bte;
		column->min.bval = 127;
		column->max.bval = -127;
	} else if (strcmp(column->typename, "char") == 0 || strcmp(argv[1],"smallint")== 0 || strcmp(argv[1], "short")== 0) {
		column->coltype  = TYPE_sht;
		column->min.sval = 32767;
		column->max.sval = -32767;
	} else if (strcmp(column->typename, "decimal") == 0 || strcmp(argv[1], "int") == 0 || strcmp(argv[1], "date") == 0) {
		column->coltype  = TYPE_int;
		column->min.ival = INT_MAX;
		column->max.ival = INT_MIN;
	} else if (strcmp(column->typename, "long") == 0 || strcmp(argv[1], "bigint") == 0) {
		column->coltype  = TYPE_lng;
		column->min.lval = LONG_MAX;
		column->max.lval = LONG_MIN;
	} else if (strcmp(column->typename, "float") == 0 || strcmp(argv[1], "real") == 0) {
		column->coltype= TYPE_flt;
		column->min.fval = FLT_MAX;
		column->max.fval = FLT_MIN;
	} else if (strcmp(column->typename, "double") == 0 ) {
		column->coltype  = TYPE_dbl;
		column->min.dval = DBL_MAX;
		column->max.dval = -DBL_MAX;
	} else if (strcmp(column->typename, "oid") == 0) {
		column->coltype  = TYPE_oid;
		column->min.lval = ULONG_MAX;
		column->max.lval = 0;
	} else {
		printf("type %s not supported\n", column->typename);
		return -1;
	}
	column->typesize = stride[column->coltype];

	cfile = fopen(column->filename, "r");
	if (cfile == NULL) {
		printf("failed to open column file %s\n", column->filename);
		return -1;
	}
	fseek(cfile, 0, SEEK_END);
	filesize = ftell(cfile);
	if (filesize == 0){
		printf("empty open column file %s\n", column->filename);
		return -1;
	}
	column->col = (char *) aligned_alloc(32, sizeof(column->col)*filesize);
	if (column->col == 0) {
		printf("malloc failed %ld\n", filesize * sizeof(column->col));
		return -1;
	}
	rewind(cfile);
	if ((rd = fread(column->col, 1, filesize, cfile)) != filesize) {
		printf("Could read only %ld of %ld bytes\n", rd, filesize);
		fclose(cfile);
		return -1;
	}
	fclose(cfile);
	devnull = fopen("/dev/null","a");
	if (devnull == NULL){
		printf("can not open /dev/null\n");
		return -1;
	}

	rpp = PAGESIZE/stride[column->coltype];
	if (rpp == 0) {
		printf("rows per pages is 0\n");
		return -1;
	}
	pages = column->colcount/rpp + 1;
	if (pages > MAX_IMPS) {
		printf("there are too many pages %ld\n", pages);
		return -1;
	}

	VERBOSE printf("%s imprint %s "
	             "filesize %ld "
	             "type %s %d "
	             "stride %d "
	             "records %ld "
	             "pagesize %d "
	             "sysconf(pagesize) %ld "
	             "rpp %d "
	             "pages %ld\n",
	             column->colname, column->filename,
	             filesize,
	             column->typename, column->coltype,
	             stride[column->coltype],
	             column->colcount,
	             PAGESIZE,
	             sysconf(_SC_PAGESIZE),
	             rpp,
	             pages);

	/* check if column is sorted and set sorted = 1 if it is */
	isSorted(column);

	/* create zonemaps */
	zonemaps = create_zonemaps(column);
	/*create imprints */
	scalar_imps = scalar_imprints1(column);

	imps2 = create_imprints(column, 64, 64, 0);

	//compareImprintsIndex(column,scalar_imps,imps2);

	//PRINT_IMPRINTS printImprint(column, zonemaps);
	/*create simd_imprints() equal sized as the original imprint */
	//simd_imps = simd_imprints(column, rpp, scalar_imps);
	//compareImprintsIndex(column,scalar_imps,simd_imps);

	VERBOSE printf("%s tuples=%ld size=%ld(bytes), zonemap_sz=%ld(bytes) %ld%% #zones=%ld, imprints_sz=%ld(bytes) %ld%%,",
	             column->colname, column->colcount, filesize,
	             zonemaps->zmaps_cnt * 2 * stride[column->coltype], ((long)zonemaps->zmaps_cnt * 2 * stride[column->coltype] * 100) / filesize, zonemaps->zmaps_cnt,
	             ((long) (imps_cnt / (BITS/scalar_imps->bins)) * sizeof(long) + dct_cnt * sizeof(Dct)),
	             100 * ((long) (imps_cnt / (BITS/scalar_imps->bins)) * sizeof(long) + dct_cnt * sizeof(Dct)) / filesize);
	VERBOSE printf(" #imprints=%ld #dict=%ld\n", imps_cnt, dct_cnt);

	/* run queries */
	queries(column, zonemaps, imps2 );
	//statistics(column);


	VERBOSE printf("end of run\n");

	/* before exiting free memory for cleaness */
	free(column->col);
	free(dct_scalar);
	free(imprints);
	free(zonemaps->zmaps);
	free(zonemaps);
	return 1;
}

Zonemap_index *
create_zonemaps(Column *column)
{
	Zonemap_index *zonemaps;
	ValRecord val;
	long i;
	long t0;
	int new = rpp-1; /*rpp is always power of 2*/

	/* malloc zonemap array */
	zonemaps = (Zonemap_index *) malloc(sizeof(Zonemap_index));
	zonemaps->zmaps = (Zonemap *) malloc (sizeof(Zonemap)*(pages+1));
	memset((char *)zonemaps->zmaps, 0, sizeof(Zonemap)*(pages+1));

#define upd(X) \
		if (val.X < column->min.X) column->min.X = val.X; \
		if (val.X > column->max.X) column->max.X = val.X; \
		if (zonemaps->zmaps[zonemaps->zmaps_cnt].min.X > val.X) \
			zonemaps->zmaps[zonemaps->zmaps_cnt].min.X = val.X; \
		if (zonemaps->zmaps[zonemaps->zmaps_cnt].max.X < val.X) \
			zonemaps->zmaps[zonemaps->zmaps_cnt].max.X = val.X;

	t0 = usec();
	zonemaps->zmaps_cnt = -1;
	for (i=0; i < column->colcount; i++) {
		if (!(i&new)) {
			zonemaps->zmaps_cnt++;
			switch (column->coltype) {
			case TYPE_bte:
				zonemaps->zmaps[zonemaps->zmaps_cnt].min.bval = 127;
				zonemaps->zmaps[zonemaps->zmaps_cnt].max.bval = -127;
				break;
			case TYPE_sht:
				zonemaps->zmaps[zonemaps->zmaps_cnt].min.sval = 32767;
				zonemaps->zmaps[zonemaps->zmaps_cnt].max.sval = -32767;
				break;
			case TYPE_int:
				zonemaps->zmaps[zonemaps->zmaps_cnt].min.ival = INT_MAX;
				zonemaps->zmaps[zonemaps->zmaps_cnt].max.ival = INT_MIN;
				break;
			case TYPE_lng:
				zonemaps->zmaps[zonemaps->zmaps_cnt].min.lval = LONG_MAX;
				zonemaps->zmaps[zonemaps->zmaps_cnt].max.lval = LONG_MIN;
				break;
			case TYPE_oid:
				zonemaps->zmaps[zonemaps->zmaps_cnt].min.ulval = ULONG_MAX;
				zonemaps->zmaps[zonemaps->zmaps_cnt].max.ulval = 0;
				break;
			case TYPE_flt:
				zonemaps->zmaps[zonemaps->zmaps_cnt].min.fval = FLT_MAX;
				zonemaps->zmaps[zonemaps->zmaps_cnt].max.fval = FLT_MIN;
				break;
			case TYPE_dbl:
				zonemaps->zmaps[zonemaps->zmaps_cnt].min.dval = DBL_MAX;
				zonemaps->zmaps[zonemaps->zmaps_cnt].max.dval = -DBL_MAX;
			}
		}

		switch(column->coltype){
			case TYPE_bte: val.bval  = *(char*)   (column->col + i*column->typesize); upd(bval); break;
			case TYPE_sht: val.sval  = *(short *) (column->col + i*column->typesize); upd(sval); break;
			case TYPE_int: val.ival  = *(int*)    (column->col + i*column->typesize); upd(ival); break;
			case TYPE_lng: val.lval  = *(long*)   (column->col + i*column->typesize); upd(lval); break;
			case TYPE_oid: val.ulval = *(unsigned long *) (column->col + i*column->typesize); upd(ulval); break;
			case TYPE_flt: val.fval  = *(float*)  (column->col + i*column->typesize); upd(fval);break;
			case TYPE_dbl: val.dval  = *(double*) (column->col + i*column->typesize); upd(dval); break;
		}
	}
	if ((i-1)%rpp) zonemaps->zmaps_cnt++;
	zonemaps->zmaps_cnt++;
	zone_create_time = usec()-t0;

	VERBOSE printf("%s zonemap  creation time=%ld, %ld usec per thousand values\n", column->colname,  zone_create_time, ((long)zone_create_time*1000)/column->colcount);
	return zonemaps;
}

Imprints_index*
scalar_imprints1(Column *column)
{
	Imprints_index *imps;
	long i, mask, prevmask;
	long t0;
	int bit;
	ValRecord val;
	int fits;
	int new;


	/* sample to create set the number of bins */
	imps = (Imprints_index *) malloc (sizeof(Imprints_index));
	imps->bounds = (ValRecord *) malloc(sizeof(ValRecord) * 64);
	binning(column, imps->bounds, &(imps->bins), 64);

	/* how many mask vectors fit in the imprint */
	fits = 64/imps->bins;
	dct_cnt = 0;
	imps_cnt = 0;
	prevmask = 0;
	mask = 0;
	new = rpp-1; /* rpp is always power of 2 */

	dct_scalar = (Dct *) malloc (sizeof(Dct)*pages);
	imprints = (long *) malloc (sizeof(long)*(pages/fits+1));

	/* init imprints */
	for (i=0, imps_cnt=pages/fits+1; i<imps_cnt; i++) {
		imprints[i]=0;
	}
	imps_cnt=0;

#define GETBIT(Z, X)							\
do {											\
	int _i;										\
	Z = 0;										\
	for (_i = 1; _i < imps->bins; _i++)			\
		Z += ((val.X) > imps->bounds[_i].X);	\
} while (0)

	t0 = usec();
	/* start creation */
	for (i=0; i<column->colcount; i++) {
		if (!(i&new) && i>0) {
			/* compress list */
			if (prevmask == mask && dct_scalar[dct_cnt-1].blks < ((1<<MAXOFFSET)-1)) {
				if (dct_scalar[dct_cnt - 1].repeated == 0) {
					if (dct_scalar[dct_cnt - 1].blks > 1) {
						dct_scalar[dct_cnt - 1].blks--; /* reduce previous by 1 */
						dct_cnt++;
						dct_scalar[dct_cnt-1].blks = 1;   /* the new is a repeat */
					}
					dct_scalar[dct_cnt-1].repeated = 1;
				}
				/* same mask as before */
				dct_scalar[dct_cnt - 1].blks++;
			} else {
				/* new mask */
				prevmask = mask;
				imps->bins == 64? (imprints[imps_cnt] = mask): (imprints[imps_cnt/fits] |= mask<<((imps_cnt%fits)*imps->bins));
				imps_cnt++;

				if (dct_cnt > 0 && dct_scalar[dct_cnt - 1].repeated == 0 && dct_scalar[dct_cnt-1].blks < ((1<<MAXOFFSET)-1)) {
					dct_scalar[dct_cnt - 1].blks++;
				} else {
					dct_scalar[dct_cnt].blks = 1;
					dct_scalar[dct_cnt].repeated = 0;
					dct_cnt++;
				}
			}
			mask = 0;
		}
		switch (column->coltype) {
			case TYPE_bte: val.bval  = *(char*)   (column->col + i*column->typesize); GETBIT(bit, bval); break;
			case TYPE_sht: val.sval  = *(short*)  (column->col + i*column->typesize); GETBIT(bit, sval); break;
			case TYPE_int: val.ival  = *(int*)    (column->col + i*column->typesize); GETBIT(bit, ival); break;
			case TYPE_lng: val.lval  = *(long*)   (column->col + i*column->typesize); GETBIT(bit, lval); break;
			case TYPE_oid: val.ulval = *(unsigned long*) (column->col + i*column->typesize); GETBIT(bit, ulval); break;
			case TYPE_flt: val.fval  = *(float*)  (column->col + i*column->typesize); GETBIT(bit, fval); break;
			case TYPE_dbl: val.dval  = *(double*) (column->col + i*column->typesize); GETBIT(bit, dval); break;
			default: bit =0;
		}
		mask = setBit(mask, bit);
	}
	/* last mask */
	if (prevmask == mask && dct_cnt > 0 && dct_scalar[dct_cnt-1].blks < (1<<MAXOFFSET)) {
		if (dct_scalar[dct_cnt - 1].repeated == 0) {
			if (dct_scalar[dct_cnt - 1].blks == 1) { /* only 1 on previous    */
				dct_scalar[dct_cnt - 1].repeated = 1;
			} else {
				dct_scalar[dct_cnt - 1].blks--;  /* reduce previous by 1 */
				dct_scalar[dct_cnt].blks = 1;    /* the new is a repeat */
				dct_scalar[dct_cnt].repeated = 1;
				imps->bins == 64? (imprints[imps_cnt] = mask): (imprints[imps_cnt/fits] |= mask<<((imps_cnt%fits)*imps->bins));
				imps_cnt++;
				dct_cnt++;
			}
		}
		/* same mask as before */
		dct_scalar[dct_cnt - 1].blks++;
	} else {
		imps->bins == 64? (imprints[imps_cnt] = mask): (imprints[imps_cnt/fits] |= mask<<((imps_cnt%fits)*imps->bins));
		imps_cnt++;
		if (dct_cnt > 0 && dct_scalar[dct_cnt - 1].repeated == 0) {
			dct_scalar[dct_cnt - 1].blks++;
		} else {
			dct_scalar[dct_cnt].blks = 1;
			dct_scalar[dct_cnt].repeated = 0;
			dct_cnt++;
		}
	}

	/* end creation, stop timer */
	imprints_create_time = usec()- t0;

	imps->imprints = (char *)imprints;
	imps->dct_cnt = dct_cnt;
	imps->imps_cnt = imps_cnt;
	imps->blocksize = rpp;  /* blocksize is values per block */
	imps->dct = dct_scalar;
	/* stats gathering */
	stats(column, imps);

	VERBOSE printf("%s imprints creation time=%ld, %ld usec per thousand values\n", column->colname, imprints_create_time, ((long)imprints_create_time*1000)/column->colcount);
	//PRINT_HISTO printHistogram(imps, "Value distribution");

	return imps;
}

void
stats(Column *column, Imprints_index *imps)
{
	long i, j, tf, k;
	int bits;
	unsigned long mask;

	for (i=0; i<imps->bins; i++)
		vectors[i]= histogram[i] = 0;

	tf = 0;
	for (i=0; i<imps->dct_cnt; i++) {
		for (k=0; k < imps->dct[i].blks; k++) {
			if (dct_scalar[i].repeated == 1)
				k = imps->dct[i].blks;
			bits = 0;
			mask = getMask(tf);
			globalmask |= mask;
			for (j=0; j<imps->bins; j++) {
				if (isSet(mask,j)) {
					bits++;
					histogram[j]++;
				}
			}
			vectors[bits]++;
			tf++;
		}
	}
}

Imprints_index* simd_imprints(Column *column, int blocksize, Imprints_index *other)
{
	Imprints_index *imps;
	int bsteps;
	int newmask;
	long i, e, b, max_imprints;

	/* simd stuff */
	__m256i *simd_bitmasks;
	__m256i zero;
	__m256i bitmasks[256];
	__m256i *restrict limits;

	char mask[32], prevmask[32];
	long t0;


	imps = (Imprints_index *) malloc (sizeof(Imprints_index));
	imps->bounds = other->bounds;
	imps->bins = other->bins;
	imps->imprintsize = imps->bins/8;
	imps->blocksize = blocksize;  /* blocksize is bytes per block */

	/* how many imprints we will need worst case */
	max_imprints = (column->colcount/(imps->blocksize/column->typesize))+1;
	imps->dct = (Dct *) malloc (sizeof(Dct) * max_imprints);
	imps->imprints = (char *) malloc (max_imprints*(imps->bins/8));

	imps->dct_cnt = 0;
	imps->imps_cnt = 0;


	simd_bitmasks = (__m256i *) imps->imprints;
	limits = aligned_alloc(32, imps->bins*sizeof(__m256i));

#define MAKE_LIMITS(SIMDTYPE, X)											\
	for (int _i = 0; _i < imps->bins; _i++)								\
		limits[_i] = _mm256_set1_##SIMDTYPE(imps->bounds[_i].X);				\

	switch (column->coltype) {
		case TYPE_bte: MAKE_LIMITS(epi8, bval); break;
		case TYPE_sht: MAKE_LIMITS(epi16, sval); break;
		case TYPE_int: MAKE_LIMITS(epi32, ival); break;
		case TYPE_lng: MAKE_LIMITS(epi64x, lval); break;
		case TYPE_oid: MAKE_LIMITS(epi64x, ulval); break;
		case TYPE_flt: MAKE_LIMITS(ps, fval); break;
		case TYPE_dbl: MAKE_LIMITS(pd, dval); break;
		default: return NULL;
	}

	/* zero simd value */
	zero = _mm256_setzero_si256();
	/* simd bitmasks */
	for (int i = 0; i < 256; i++) {
		bitmasks[i] = setbit_256(zero, i);
	}

#define GETBIT_SIMD(SIMDTYPE)																\
	/* perform 2 bin comparisons per simd instruction until all bins are checked */			\
		for (int bin1 = 0, bin2 = 1; bin1 < imps->bins-2; bin1+=2, bin2+=2) {				\
			result = _mm256_add_##SIMDTYPE(result,											\
						_mm256_add_##SIMDTYPE(												\
							_mm256_cmpgt_##SIMDTYPE(values_v, limits[bin1]),				\
							_mm256_cmpgt_##SIMDTYPE(values_v, limits[bin2])));				\
		}																					\
		result = _mm256_sub_##SIMDTYPE(zero,result);										\
		for (int i1 = 0, i2=1; i1 < bsteps-1; i1+=2, i2+=2) {								\
			simd_mask = _mm256_or_si256(													\
				_mm256_or_si256(simd_mask, bitmasks[_mm256_extract_##SIMDTYPE(result, i1)]),\
				_mm256_or_si256(simd_mask, bitmasks[_mm256_extract_##SIMDTYPE(result, i2)])	\
			);																				\
		}

#define GETBIT_SIMDD(SIMDTYPE) return NULL;
#define GETBIT_SIMDF(SIMDTYPE) return NULL;

	bsteps = 256 / (column->typesize*8);

	t0 = usec();
	/* start creation */
	for (i = 0; i < column->colcount;) {
		__m256i simd_mask = zero;

		for (b = 0; b < imps->blocksize && i < column->colcount; b += bsteps) {
			__m256i values_v = _mm256_load_si256((__m256i*) (column->col+i*column->typesize));
			__m256i result = zero;
			switch (column->coltype) {
				case TYPE_bte: GETBIT_SIMD(epi8); break;
				case TYPE_sht: GETBIT_SIMD(epi16); break;
				case TYPE_int: GETBIT_SIMD(epi32); break;
				case TYPE_lng: GETBIT_SIMD(epi64); break;
				case TYPE_oid: GETBIT_SIMD(epi64); break;
				case TYPE_flt: GETBIT_SIMDF(ps); break;
				case TYPE_dbl: GETBIT_SIMDD(pd); break;
				default: return NULL;;
			}
			i += bsteps;
		}

		/* simd the prevmask == mask before extract */
		newmask = 0;
		for (e = 0; e < imps->imprintsize; e++) {
			mask[e] = _mm256_extract_epi8(simd_mask, e);
			if (mask[e] != prevmask[e]) {
				prevmask[e] = mask[e];
				newmask = 1;
			}
		}


		if (!newmask && imps->dct[imps->dct_cnt-1].blks < ((1<<MAXOFFSET)-1)) {
			if (imps->dct[imps->dct_cnt - 1].repeated == 0) {
				if (imps->dct[imps->dct_cnt - 1].blks > 1) {
					imps->dct[imps->dct_cnt - 1].blks--;
					imps->dct_cnt++;
					imps->dct[imps->dct_cnt-1].blks = 1;
				}
				imps->dct[imps->dct_cnt-1].repeated = 1;
			}
			/* same mask as before */
			imps->dct[imps->dct_cnt - 1].blks++;
		} else {
			unsigned long pos = imps->imps_cnt*imps->imprintsize;

			for (e = 0; e < imps->imprintsize; e++)
				imps->imprints[pos+e] = mask[e];
			imps->imps_cnt++;

			if (imps->dct_cnt > 0 && imps->dct[imps->dct_cnt - 1].repeated == 0 && imps->dct[imps->dct_cnt-1].blks < ((1<<MAXOFFSET)-1)) {
				imps->dct[imps->dct_cnt - 1].blks++;
			} else {
				imps->dct[imps->dct_cnt].blks = 1;
				imps->dct[imps->dct_cnt].repeated = 0;
				imps->dct_cnt++;
			}

			/*
			printMask(imps->imprints+(imps->imps_cnt-1)*8, 8);
			putchar(' ');
			printMask(other->imprints+(imps->imps_cnt-1)*8, 8);
			putchar('\n');*/
		}
	}

	/* end creation, stop timer */
	simd_imprints_create_time = usec() - t0;


	VERBOSE printf("%s simd_imprints creation time=%ld, %ld usec per thousand values\n", column->colname, simd_imprints_create_time, ((long)simd_imprints_create_time*1000)/column->colcount);
	VERBOSE printf("%s %ld=%ld %ld=%ld\n", column->colname, other->imps_cnt, imps->imps_cnt, other->dct_cnt, imps->dct_cnt);

	return imps;

}

void queries(Column *column, Zonemap_index *zonemaps, Imprints_index *imps)
{
	unsigned long *oids, oid, k, j, lim, n, tf, l;
	long m;
	unsigned long tuples[REPETITION];
	long basetime,                 zonetime,                 impstime;
	long basetimer[REPETITION],    zonetimer[REPETITION],    impstimer[REPETITION];
	long bindex[REPETITION],       zindex[REPETITION],       iindex[REPETITION];
	long bcomparisons[REPETITION], zcomparisons[REPETITION], icomparisons[REPETITION];
	int i;

	unsigned char  *imprints8  = (unsigned char *) imprints;
	unsigned short *imprints16 = (unsigned short *)imprints;
	unsigned int   *imprints32 = (unsigned int *)  imprints;
	unsigned long  *imprints64 = (unsigned long *) imprints;

	oids  = (unsigned long *) malloc(column->colcount * sizeof(unsigned long));
	oid   = 0;
	m     = 0;

	for (i = 0; i < REPETITION; i++) {
		tuples[i] = 0;
		basetimer[i] = impstimer[i] = zonetimer[i] = 0;
		bindex[i] = iindex[i] = zindex[i] = 0;
		bcomparisons[i] = icomparisons[i] = zcomparisons[i] = 0;
	}

	for (i = 0; i < REPETITION; i++) {
		/* select a random range from the pool, leads to bias to skew data distribution */
		/* use [slow,shigh) range expression */
		genQueryRange(column, imps, i);

		/* simple scan */
		tuples[i] = simple_scan(column,slow, shigh, &(basetimer[i]));

		/* column imprint filter */
		n   = 0;
		tf  = 0;
		m   = 0;
		oid = 0;
		#define impsquery(X, T, B)											\
			impstime = usec();													\
			for (j = 0; j < dct_cnt; j++) {										\
				if (dct_scalar[j].repeated == 0) {									\
					for (k = tf + dct_scalar[j].blks; tf < k; n++, tf++) {			\
						STATS iindex[i] += 1;									\
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
									STATS icomparisons[i] += 1;					\
									if (val <= shigh.X && val > slow.X) {		\
										oids[oid++] = l;						\
									}											\
								}												\
							}													\
						}														\
					}															\
				} else { /* repeated mask case */								\
					STATS iindex[i] += 1;										\
					if (imprints##B[tf] & mask) {								\
						register T val;											\
						l = n * rpp;											\
						lim = l + rpp * dct_scalar[j].blks;						\
						lim = lim > column->colcount ? column->colcount : lim;					\
						if ((imprints##B[tf] & ~innermask) == 0) {				\
							for (; l < lim; l++) {								\
								oids[oid++] = l;								\
							}													\
						} else {												\
							for (val = ((T*)column->col)[l]; l < lim; l++, val = ((T*)column->col)[l]) {	\
								STATS icomparisons[i] += 1;						\
								if (val <= shigh.X && val > slow.X  ) {			\
									oids[oid++] = l;							\
								}												\
							}													\
						}														\
					}															\
					n += dct_scalar[j].blks;										\
					tf++;														\
				}\
			}

		#define binchoose(X,T)					\
			switch(imps->bins) {						\
			case 8:  impsquery(X,T,8); break;	\
			case 16: impsquery(X,T,16); break;	\
			case 32: impsquery(X,T,32); break;	\
			case 64: impsquery(X,T,64); break;	\
			default: break;						\
			}

		switch(column->coltype){
		case TYPE_bte:
			binchoose(bval,char);
			break;
		case TYPE_sht:
			binchoose(sval,short);
			break;
		case TYPE_int:
			binchoose(ival,int);
			break;
		case TYPE_lng:
			binchoose(lval,long);
			break;
		case TYPE_oid:
			binchoose(ulval,unsigned long);
			break;
		case TYPE_flt:
			binchoose(fval,float);
			break;
			case TYPE_dbl:
			binchoose(dval,double);
		}

		impstimer[i] = usec() - impstime;
		if (tuples[i] != oid)
			printf("%s base %ld imprints %ld differ\n", column->colname, tuples[i], oid);

		simd_queries(column, imps, tuples[i]);
	}


	for (i =0; i< REPETITION; i++) {
		VERBOSE printf("%s query[%d]=%ld\t selectivity=%2.1f%%\t|\tscan = %ld\timprints = %ld\tzone = %ld\t(usec)\n",
		       column->colname, i, tuples[i], tuples[i]* 100.0/column->colcount, basetimer[i],
		       impstimer[i], zonetimer[i]);
		STATS printf ("bindex %ld bcomparisons %ld zindex %ld zcomparisons %ld iindex %ld icomparisons %ld\n",
		       bindex[i], bcomparisons[i], zindex[i], zcomparisons[i], iindex[i], icomparisons[i]);

		if (i) {
			basetimer[0] += basetimer[i];
			impstimer[0] += impstimer[i];
			zonetimer[0] += zonetimer[i];
		}
	}
	free(oids);
}

void
simd_queries(Column *column, Imprints_index *imps, long results) {
	long simd_impstime = 0;
	__m256i low, high;
	unsigned long dcnt; /* dctn is an increment for dct_cnt  */
	unsigned long icnt; /* ictn is an increment for imps_cnt */
	unsigned long bcnt; /* bctn is an increment for blocks   */
	unsigned long top_icnt, b, top_b, bstep;
	unsigned int v_idx;


	unsigned long *restrict oids, oid; /* for materializing the result */

	unsigned char  *imprints8  = (unsigned char *) imprints;
	unsigned short *imprints16 = (unsigned short *)imprints;
	unsigned int   *imprints32 = (unsigned int *)  imprints;
	unsigned long  *imprints64 = (unsigned long *) imprints;

	oids  = (unsigned long *) malloc(column->colcount * sizeof(unsigned long));
	oid   = 0;

	//simd_bitmasks = (__m256i *) imps->imprints;
	//low  = aligned_alloc(32, sizeof(__m256i));
	//high = aligned_alloc(32, sizeof(__m256i));

	bstep = 256 / (column->typesize*8); /* how many values we can fit in one simd register */

#define SIMD_QUERYBOUNDS(SIMDTYPE, X)					\
	low  = _mm256_set1_##SIMDTYPE(slow.X);				\
	high = _mm256_set1_##SIMDTYPE(shigh.X);

	switch (column->coltype) {
		case TYPE_bte: SIMD_QUERYBOUNDS(epi8, bval); break;
		case TYPE_sht: SIMD_QUERYBOUNDS(epi16, sval); break;
		case TYPE_int: SIMD_QUERYBOUNDS(epi32, ival); break;
		case TYPE_lng: SIMD_QUERYBOUNDS(epi64x, lval); break;
		case TYPE_oid: SIMD_QUERYBOUNDS(epi64x, ulval); break;
		case TYPE_flt: SIMD_QUERYBOUNDS(ps, fval); break;
		case TYPE_dbl: SIMD_QUERYBOUNDS(pd, dval); break;
		default: return;
	}

icnt=bcnt=0;

	/* simd imprints queries */
	#define SIMD_IMPSQUERY(X, T, B, SIMDTYPE)													\
		simd_impstime = usec();													\
		for (dcnt = 0; dcnt < imps->dct_cnt; dcnt++) {							\
			if (imps->dct[dcnt].repeated == 0) {								\
				top_icnt = icnt + imps->dct[dcnt].blks;						\
				for (; icnt < top_icnt; bcnt++, icnt++) {						\
					if (imprints##B[icnt] & mask) {								\
						b = bcnt * imps->blocksize;								\
						top_b = b + imps->blocksize;							\
						top_b = top_b > column->colcount ? column->colcount : top_b;	\
						if ((imprints##B[icnt] & ~innermask) == 0) {			\
							for (; b < top_b; b++) {							\
								oids[oid++] = b;								\
							}													\
						} else {												\
							for (; b < top_b; b+=bstep) {						\
								__m256i values_v = _mm256_load_si256((__m256i*) (column->col+b*column->typesize)); \
								v_idx = _mm256_movemask_epi8(\
								_mm256_sub_##SIMDTYPE(\
									_mm256_cmpgt_##SIMDTYPE(values_v, low),	\
									_mm256_cmpgt_##SIMDTYPE(values_v, high)));	\
								for (int i = 0; i < imps->blocksize; i++) { \
									if (v_idx & (1 << (i*column->typesize)) ) oids[oid++] = b + i;\
								}\
							}											\
						}\
					}															\
				}																\
			} else {  /* repeated mask case */									\
					if (imprints##B[icnt] & mask) {								\
						b = bcnt * imps->blocksize;										\
						top_b = b + imps->blocksize * imps->dct[dcnt].blks;					\
						top_b = top_b > column->colcount ? column->colcount : top_b;	\
						if ((imprints##B[icnt] & ~innermask) == 0) {			\
							for (; b < top_b; b++) {								\
								oids[oid++] = b;								\
							}													\
						} else {												\
							for (; b < top_b; b+=bstep) {						\
								__m256i values_v = _mm256_load_si256((__m256i*) (column->col+b*column->typesize)); \
								v_idx = _mm256_movemask_epi8(\
								_mm256_sub_##SIMDTYPE(\
									_mm256_cmpgt_##SIMDTYPE(values_v, low),	\
									_mm256_cmpgt_##SIMDTYPE(values_v, high)));	\
								for (int i = 0; i < imps->blocksize; i++) { \
									if (v_idx & (1 << (i*column->typesize)) ) oids[oid++] = b + i;\
								}\
							}											\
						}														\
					}															\
					bcnt += imps->dct[dcnt].blks;									\
					icnt++;														\
				}\
			}

	#define SIMD_BINCHOOSE(X,T,SIMDTYPE)					\
		switch(imps->bins) {						\
		case 8:  SIMD_IMPSQUERY(X,T,8,SIMDTYPE); break;	\
		case 16: SIMD_IMPSQUERY(X,T,16,SIMDTYPE); break;	\
		case 32: SIMD_IMPSQUERY(X,T,32,SIMDTYPE); break;	\
		case 64: SIMD_IMPSQUERY(X,T,64,SIMDTYPE); break;	\
		default: break;						\
		}

		switch(column->coltype){
		case TYPE_bte:
			SIMD_BINCHOOSE(bval,char,epi8);
			break;
		case TYPE_sht:
			SIMD_BINCHOOSE(sval,short,epi16);
			break;
		case TYPE_int:
			SIMD_BINCHOOSE(ival,int,epi32);
			break;
		case TYPE_lng:
			SIMD_BINCHOOSE(lval,long,epi32);
			break;
		case TYPE_oid:
			SIMD_BINCHOOSE(ulval,unsigned long,epi64);
			break;
		case TYPE_flt:
			SIMD_BINCHOOSE(fval,float,epi32);
			break;
		case TYPE_dbl:
			SIMD_BINCHOOSE(dval,double,epi64);
		}

		simd_impstime = usec() - simd_impstime;
		printf("%s simd_imprints time = %ld \n", column->colname, simd_impstime );
		if (results != oid)
			printf("%s base %ld simd_imprints %ld differ\n", column->colname, results, oid);

}

/* simulate a series of queries */
void genQueryRange(Column *column, Imprints_index *imps, int i)
{
	long low, high;
	int j;
	int lastbit;

	mask = 0;
	innermask = 0;

	/* select a range based on the actual non-empty bins */
	for (lastbit = imps->bins-1; lastbit >0; lastbit--)
		if (isSet(globalmask, lastbit))
			break;

	low = 0 + (int)(rand() * 1.0 / RAND_MAX * lastbit);
	high = low + i * lastbit/ (100.0/REPETITION);
	if (high > lastbit)
		high = lastbit;
	if (low > high)
		low = high;
	
	do {
		/* find at least one non-empty bin */
		for (j = low; j <= high; j++)
			if ( histogram[j]) goto foundrange;
		if ( low > 0){
			low--;
			high--;
		} else {
			high = lastbit + (high-low);
			low= lastbit;
		}
	} while (1);

foundrange:
	for (; low < high; low++)
		if (histogram[low]) break;

	for (; high > low; high--)
		if (histogram[high]) break;

	for (j = low; j <= high; j++) {
		mask = setBit(mask, j);
	}
	/* inner mask should only be set when range bounds call for it */
	/* in this case we use full bins */
	for (j = low+1; j <= high-1; j++) {
		innermask = setBit(innermask, j);
	}

#define setqueryrange(X)						\
	slow.X = imps->bounds[low].X;						\
	shigh.X = imps->bounds[high].X;					\
	/*
	PRINT_QUERIES {								\
		printf("query             ");			\
		printMask(mask,BITS); putchar('\n');	\
		printf("inner msk         ");			\
		printMask(innermask,BITS);				\
		putchar('\n');							\
	}*/

	switch (column->coltype) {
	case TYPE_bte:
		setqueryrange(bval);
		break;
	case TYPE_sht:
		setqueryrange(sval);
		break;
	case TYPE_int:
		setqueryrange(ival);
		break;
	case TYPE_lng:
		setqueryrange(lval);
		break;
	case TYPE_oid:
		setqueryrange(ulval);
		break;
	case TYPE_flt:
		setqueryrange(fval);
		break;
	case TYPE_dbl:
		setqueryrange(dval);
	}
	return;
}
