#include "simd_imprints.h"

/* global bounds */
ValRecord absmin, absmax, mxbins[BITS], mibins[BITS];

/* auxilary variables used globaly */
FILE *devnull;
unsigned long pages;                            /* total pages in the column */
int rpp;                                        /* rows per page */
int sorted;                                     /* column is sorted */

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
int bins;

/* stat vars for imprints */
long histogram[BITS]; /* bin filling */
long vectors[BITS+1]; /* vector filling distribution */


/* global vars for queries */
ValRecord slow, shigh;
unsigned long mask;
unsigned long innermask;

/* functions (in call order)-ish */
void isSorted(Column column);
Zonemap_index *create_zonemaps(Column column);
void scalar_imprints(Column column);
void imps_sample(Column column);
void imps_histogram(Column column, ValRecord *sample, int smp);
void stats();
void queries(Column column, Zonemap_index *zonemaps);
void genQueryRange(Column column, int i);
void printHistogram(long histo[BITS], char *name);
void printBins(Column column);
void printMask(long mask, int limit);
void printImprint(Column column, Zonemap_index *zonemaps);
void statistics();
Imprints_index* simd_imprints(Column column, int blocksize, int imprint_bins);



/* main function for stand alone imprints */
int main(int argc, char **argv)
{
	Column column;
	FILE *cfile;
	long filesize;
	size_t rd;
	Zonemap_index *zonemaps;
	Imprints_index *simd_imps_orig;

	if (argc != 5) {
		printf("usage: %s type count file column\n", argv[0]);
		return -1;
	}

	strcpy(column.colname, argv[4]);
	strcpy(column.filename, argv[3]);
	column.colcount = atoi(argv[2]);
	strcpy(column.typename, argv[1]);

	if (strcmp(column.typename, "tinyint") == 0 || strcmp(argv[1], "boolean") == 0) {
		column.coltype  = TYPE_bte;
		column.min.bval = 127;
		column.max.bval = -127;
	} else if (strcmp(column.typename, "char") == 0 || strcmp(argv[1],"smallint")== 0 || strcmp(argv[1], "short")== 0) {
		column.coltype  = TYPE_sht;
		column.min.sval = 32767;
		column.max.sval = -32767;
	} else if (strcmp(column.typename, "decimal") == 0 || strcmp(argv[1], "int") == 0 || strcmp(argv[1], "date") == 0) {
		column.coltype  = TYPE_int;
		column.min.ival = INT_MAX;
		column.max.ival = INT_MIN;
	} else if (strcmp(column.typename, "long") == 0 || strcmp(argv[1], "bigint") == 0) {
		column.coltype  = TYPE_lng;
		column.min.lval = LONG_MAX;
		column.max.lval = LONG_MIN;
	} else if (strcmp(column.typename, "float") == 0 || strcmp(argv[1], "real") == 0) {
		column.coltype= TYPE_flt;
		column.min.fval = FLT_MAX;
		column.max.fval = FLT_MIN;
	} else if (strcmp(column.typename, "double") == 0 ) {
		column.coltype  = TYPE_dbl;
		column.min.dval = DBL_MAX;
		column.max.dval = -DBL_MAX;
	} else if (strcmp(column.typename, "oid") == 0) {
		column.coltype  = TYPE_oid;
		column.min.lval = ULONG_MAX;
		column.max.lval = 0;
	} else {
		printf("type %s not supported\n", column.typename);
		return -1;
	}
	absmin = column.max;
	absmax = column.min;

	cfile = fopen(column.filename, "r");
	if (cfile == NULL) {
		printf("failed to open column file %s\n", column.filename);
		return -1;
	}
	fseek(cfile, 0, SEEK_END);
	filesize = ftell(cfile);
	if (filesize == 0){
		printf("empty open column file %s\n", column.filename);
		return -1;
	}
	column.col = (char *) aligned_alloc(32, sizeof(column.col)*filesize);
	if (column.col == 0) {
		printf("malloc failed %ld\n", filesize * sizeof(column.col));
		return -1;
	}
	rewind(cfile);
	if ((rd = fread(column.col, 1, filesize, cfile)) != filesize) {
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

	rpp = PAGESIZE/stride[column.coltype];
	if (rpp == 0) {
		printf("rows per pages is 0\n");
		return -1;
	}
	pages = column.colcount/rpp + 1;
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
	             column.colname, column.filename,
	             filesize,
	             column.typename, column.coltype,
	             stride[column.coltype],
	             column.colcount,
	             PAGESIZE,
	             sysconf(_SC_PAGESIZE),
	             rpp,
	             pages);

	/* check if column is sorted and set sorted = 1 if it is */
	isSorted(column);

	/* create zonemaps */
	zonemaps = create_zonemaps(column);
	/*create imprints */
	scalar_imprints(column);
	PRINT_IMPRINTS printImprint(column, zonemaps);
	/*create simd_imprints() equal sized as the original imprint */
	simd_imps_orig = simd_imprints(column, rpp, bins);

	VERBOSE printf("%s tuples=%ld size=%ld(bytes), zonemap_sz=%ld(bytes) %ld%% #zones=%ld, imprints_sz=%ld(bytes) %ld%%,",
	             column.colname, column.colcount, filesize,
	             zonemaps->zmaps_cnt * 2 * stride[column.coltype], ((long)zonemaps->zmaps_cnt * 2 * stride[column.coltype] * 100) / filesize, zonemaps->zmaps_cnt,
	             ((long) (imps_cnt / (BITS/bins)) * sizeof(long) + dct_cnt * sizeof(Dct)),
	             100 * ((long) (imps_cnt / (BITS/bins)) * sizeof(long) + dct_cnt * sizeof(Dct)) / filesize);
	VERBOSE printf(" #imprints=%ld #dict=%ld\n", imps_cnt, dct_cnt);

	/* run queries */
	queries(column, zonemaps);
	statistics(column);


	VERBOSE printf("end of run\n");

	/* before exiting free memory for cleaness */
	free(column.col);
	free(dct_scalar);
	free(imprints);
	free(zonemaps->zmaps);
	free(zonemaps);
	return 1;
}

/* check if a column is sorted */
void isSorted(Column column)
{
	int i;

#define checksorted(T)								\
	for (sorted=1, i = 1; i<column.colcount; i++)			\
		if (((T*)column.col)[i] < ((T*)column.col)[i-1]) {		\
			sorted = 0;								\
			break;									\
		}											\
	if (sorted == 0)								\
		for (sorted=1, i = 1; i<column.colcount; i++)		\
			if (((T*)column.col)[i] > ((T*)column.col)[i-1]) {	\
				sorted = 0;							\
				break;								\
			}

	switch (column.coltype) {
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

	VERBOSE printf("%s sorted property %d\n", column.colname, sorted);
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

Zonemap_index *
create_zonemaps(Column column)
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
		if (val.X < column.min.X) column.min.X = val.X; \
		if (val.X > column.max.X) column.max.X = val.X; \
		if (zonemaps->zmaps[zonemaps->zmaps_cnt].min.X > val.X) \
			zonemaps->zmaps[zonemaps->zmaps_cnt].min.X = val.X; \
		if (zonemaps->zmaps[zonemaps->zmaps_cnt].max.X < val.X) \
			zonemaps->zmaps[zonemaps->zmaps_cnt].max.X = val.X;

	t0 = usec();
	zonemaps->zmaps_cnt = -1;
	for (i=0; i < column.colcount; i++) {
		if (!(i&new)) {
			zonemaps->zmaps_cnt++;
			switch (column.coltype) {
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

		switch(column.coltype){
			case TYPE_bte: val.bval  = *(char*)   (column.col + i*stride[column.coltype]); upd(bval); break;
			case TYPE_sht: val.sval  = *(short *) (column.col + i*stride[column.coltype]); upd(sval); break;
			case TYPE_int: val.ival  = *(int*)    (column.col + i*stride[column.coltype]); upd(ival); break;
			case TYPE_lng: val.lval  = *(long*)   (column.col + i*stride[column.coltype]); upd(lval); break;
			case TYPE_oid: val.ulval = *(unsigned long *) (column.col + i*stride[column.coltype]); upd(ulval); break;
			case TYPE_flt: val.fval  = *(float*)  (column.col + i*stride[column.coltype]); upd(fval);break;
			case TYPE_dbl: val.dval  = *(double*) (column.col + i*stride[column.coltype]); upd(dval); break;
		}
	}
	if ((i-1)%rpp) zonemaps->zmaps_cnt++;
	zonemaps->zmaps_cnt++;
	zone_create_time = usec()-t0;

	VERBOSE printf("%s zonemap  creation time=%ld, %ld usec per thousand values\n", column.colname,  zone_create_time, ((long)zone_create_time*1000)/column.colcount);
	return zonemaps;
}

void
scalar_imprints(Column column)
{
	long i, mask, prevmask;
	long t0;
	int bit;
	ValRecord val;
	int fits;
	int new;


	/* sample to create set the number of bins */
	imps_sample(column);

	/* how many mask vectors fit in the imprint */
	fits = BITS/bins;
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

#define GETBIT(Z, X)					\
do {									\
	int _i;								\
	Z = 0;								\
	for (_i = 1; _i < bins; _i++)		\
		Z += ((val.X) >= mibins[_i].X);	\
} while (0)

	t0 = usec();
	/* start creation */
	for (i=0; i<column.colcount; i++) {
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
				bins == 64? (imprints[imps_cnt] = mask): (imprints[imps_cnt/fits] |= mask<<((imps_cnt%fits)*bins));
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
		switch (column.coltype) {
			case TYPE_bte: val.bval  = *(char*)   (column.col + i*stride[column.coltype]); GETBIT(bit, bval); break;
			case TYPE_sht: val.sval  = *(short*)  (column.col + i*stride[column.coltype]); GETBIT(bit, sval); break;
			case TYPE_int: val.ival  = *(int*)    (column.col + i*stride[column.coltype]); GETBIT(bit, ival); break;
			case TYPE_lng: val.lval  = *(long*)   (column.col + i*stride[column.coltype]); GETBIT(bit, lval); break;
			case TYPE_oid: val.ulval = *(unsigned long*) (column.col + i*stride[column.coltype]); GETBIT(bit, ulval); break;
			case TYPE_flt: val.fval  = *(float*)  (column.col + i*stride[column.coltype]); GETBIT(bit, fval); break;
			case TYPE_dbl: val.dval  = *(double*) (column.col + i*stride[column.coltype]); GETBIT(bit, dval); break;
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
				bins == 64? (imprints[imps_cnt] = mask): (imprints[imps_cnt/fits] |= mask<<((imps_cnt%fits)*bins));
				imps_cnt++;
				dct_cnt++;
			}
		}
		/* same mask as before */
		dct_scalar[dct_cnt - 1].blks++;
	} else {
		bins == 64? (imprints[imps_cnt] = mask): (imprints[imps_cnt/fits] |= mask<<((imps_cnt%fits)*bins));
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

	/* stats gathering */
	stats();

	VERBOSE printf("%s imprints creation time=%ld, %ld usec per thousand values\n", column.colname, imprints_create_time, ((long)imprints_create_time*1000)/column.colcount);
	PRINT_HISTO printHistogram(histogram, "Value distribution");

}

#define cmpvalues(X)	\
static int	\
cmpvalues##X(const void *p1, const void *p2)\
{\
	ValRecord v1 = *(ValRecord *) p1;\
	ValRecord v2 = *(ValRecord *) p2;\
	if ( v1.X < v2.X ) return -1; \
	if ( v1.X > v2.X ) return 1; \
	return 0;\
}

cmpvalues(bval);
cmpvalues(sval);
cmpvalues(ival);
cmpvalues(lval);
cmpvalues(ulval);
cmpvalues(fval);
cmpvalues(dval);



void imps_sample(Column column)
{
	ValRecord sample[SAMPLE_SZ];
	int i,j;

#define sampleDistribution(X,T)												\
	/* draw the sample, */													\
	for (j=0, i=0; i<SAMPLE_SZ; i++) {										\
		j = (rand()*column.colcount)/RAND_MAX;								\
		sample[i].X = ((T*)column.col)[j];									\
	}																		\
	/* sort the sample and make it unique*/									\
	qsort((char *) sample, SAMPLE_SZ, sizeof(ValRecord), cmpvalues##X);		\
	for (j=0, i=1; i<SAMPLE_SZ; i++) {										\
		if (sample[i].X != sample[j].X) {									\
			sample[++j] = sample[i];										\
		}																	\
	}																		\
	j++;

	switch(column.coltype){
	case TYPE_bte:
		sampleDistribution(bval, char);
		break;
	case TYPE_sht:
		sampleDistribution(sval, short);
		break;
	case TYPE_int:
		sampleDistribution(ival, int);
		break;
	case TYPE_lng:
		sampleDistribution(lval, long);
		break;
	case TYPE_oid:
		sampleDistribution(ulval, unsigned long);
		break;
	case TYPE_flt:
		sampleDistribution(fval, float);
		break;
	case TYPE_dbl:
		sampleDistribution(dval, double);
	}
	imps_histogram(column, sample, j);
}

void
imps_histogram(Column column, ValRecord *sample, int smp) {
	int k;

	mibins[0] = absmin;
	mxbins[0] = sample[0];
	mibins[BITS-1] = sample[smp-1];
	mxbins[BITS-1] = absmax;
	bins = 64;

	if (smp < BITS-1) {
		for (k=1; k < smp; k++) {
			mibins[k] = mxbins[k-1];
			mxbins[k] = sample[k];
		}
		if (k<8) bins=8;
		if (8<=k && k<16) bins=16;
		if (16<=k && k<32) bins=32;
		if (32<=k && k<64) bins=64;
		for (; k < BITS-1; k++) {
			mibins[k] = sample[smp-1];
			mxbins[k] = absmax;
		}
	} else {
		double y, ystep = (double)smp/(double)(BITS-2);
		for (k=1, y=ystep; y < smp; y+= ystep, k++) {
			mibins[k] = mxbins[k-1];
			mxbins[k] = sample[(int)y];
		}
		if (k == BITS-2) { /* there is some leftovers */
			assert (y>=smp);
			assert ((y-ystep)<smp);
			mibins[k] = mxbins[k-1];
			mxbins[k] = sample[smp-1];
		}
	}

	VERBOSE printf("%s sample gives %d unique values from %d and %d bins\n", column.colname, smp, SAMPLE_SZ, bins);

	return;
}

void
stats()
{
	long i, j, tf, k;
	int bits;
	unsigned long mask;

	for (i=0; i<bins; i++)
		vectors[i]= histogram[i] = 0;

	tf = 0;
	for (i=0; i<dct_cnt; i++) {
		for (k=0; k < dct_scalar[i].blks; k++) {
			if (dct_scalar[i].repeated == 1)
				k = dct_scalar[i].blks;
			bits = 0;
			mask = getMask(tf);
			globalmask |= mask;
			for (j=0; j<bins; j++) {
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

__m256i setbit_256(__m256i x,int k){
	// constants that will (hopefully) be hoisted out of a loop after inlining  
	__m256i indices = _mm256_set_epi32(224,192,160,128,96,64,32,0);
	__m256i one = _mm256_set1_epi32(-1);
	one = _mm256_srli_epi32(one, 31);    // set1(0x1)
	__m256i kvec = _mm256_set1_epi32(k);  
	// if 0<=k<=255 then kvec-indices has exactly one element with a value between 0 and 31
	__m256i shiftcounts = _mm256_sub_epi32(kvec, indices);
	__m256i kbit        = _mm256_sllv_epi32(one, shiftcounts);   // shift counts outside 0..31 shift the bit out of the element
	return _mm256_or_si256(kbit, x);                             // use _mm256_andnot_si256 to unset the k-th bit
}


Imprints_index* simd_imprints(Column column, int blocksize, int imprints_bins)
{
	Imprints_index *imps;
	int bsteps;
	long i, b, page;

	/* simd stuff */
	__m256i *simd_bitmasks;
	__m256i zero;
	__m256i bitmasks[256];
	__m256i *restrict limits;

	long long mask, prevmask;
	long t0;


	imps = (Imprints_index *) malloc (sizeof(Imprints_index));
	imps->blocksize = blocksize;  /* blocksize is values per block */
	imps->bins = imprints_bins;   /* aka how many bits per imprint */
	page = column.colcount/blocksize + 1; /* aka how many imprints we will need worst case */
	imps->dct = (Dct *) malloc (sizeof(Dct) * page);
	imps->imprints = (char *) malloc (page*(imps->bins/8));
	imps->dct_cnt = 0;
	imps->imps_cnt = 0;


	simd_bitmasks = (__m256i *) imps->imprints;
	limits = aligned_alloc(32, imps->bins*sizeof(__m256i));

#define MAKE_LIMITS(SIMDTYPE, X)											\
	for (int _i = 0; _i < imps->bins; _i++)								\
		limits[_i] = _mm256_set1_##SIMDTYPE(mibins[_i].X);				\

	switch (column.coltype) {
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

	prevmask = 0;



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

	bsteps = 256 / (stride[column.coltype]*8);

	t0 = usec();
	/* start creation */
	for (i = 0; i < column.colcount;) {
		__m256i simd_mask = zero;

		for (b = 0; b < imps->blocksize && i < column.colcount; b += bsteps) {
			__m256i values_v = _mm256_load_si256((__m256i*) (column.col+i*stride[column.coltype]));
			__m256i result = zero;
			switch (column.coltype) {
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

		switch (imps->bins) {
			case 8:
				mask = _mm256_extract_epi8(simd_mask, 0);
				break;
			case 16:
				mask = _mm256_extract_epi16(simd_mask, 0);
				break;
			case 32:
				mask = _mm256_extract_epi32(simd_mask, 0);
				break;
			case 64:
				mask = _mm256_extract_epi64(simd_mask, 0);
				break;
			case 128:
				break;
			case 256:
				break;
			default:
				break;
		}

		if (mask == prevmask && imps->dct[imps->dct_cnt-1].blks < ((1<<MAXOFFSET)-1)) {
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
			/* new mask */
			prevmask = mask;
			imps->bins == 64 ? (imps->imprints[imps->imps_cnt] = mask): (imps->imprints[imps->imps_cnt/(BITS/imps->bins)] |= mask<<((imps->imps_cnt%(BITS/imps->bins))*imps->bins));
			imps->imps_cnt++;

			if (imps->dct_cnt > 0 && imps->dct[imps->dct_cnt - 1].repeated == 0 && imps->dct[imps->dct_cnt-1].blks < ((1<<MAXOFFSET)-1)) {
				imps->dct[imps->dct_cnt - 1].blks++;
			} else {
				imps->dct[imps->dct_cnt].blks = 1;
				imps->dct[imps->dct_cnt].repeated = 0;
				imps->dct_cnt++;
			}
		}
		/* printMask(mask, imps->bins); putchar('\n'); */
	}

	/* end creation, stop timer */
	simd_imprints_create_time = usec() - t0;


	VERBOSE printf("%s simd_imprints creation time=%ld, %ld usec per thousand values\n", column.colname, simd_imprints_create_time, ((long)simd_imprints_create_time*1000)/column.colcount);
	VERBOSE printf("%s %ld=%ld %ld=%ld\n", column.colname, imps_cnt, imps->imps_cnt, dct_cnt, imps->dct_cnt);

	return imps;

}

void queries(Column column, Zonemap_index *zonemaps)
{
	unsigned long *oids, oid, k, j, lim, n, tf, l;
	long m;
	long tuples[REPETITION];
	long basetime,                 zonetime,                 impstime;
	long basetimer[REPETITION],    zonetimer[REPETITION],    impstimer[REPETITION];
	long bindex[REPETITION],       zindex[REPETITION],       iindex[REPETITION];
	long bcomparisons[REPETITION], zcomparisons[REPETITION], icomparisons[REPETITION];
	int i;

	unsigned char  *imprints8  = (unsigned char *) imprints;
	unsigned short *imprints16 = (unsigned short *)imprints;
	unsigned int   *imprints32 = (unsigned int *)  imprints;
	unsigned long  *imprints64 = (unsigned long *) imprints;

	oids  = (unsigned long *) malloc(column.colcount * sizeof(unsigned long));
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
		genQueryRange(column, i);

		/* simple scan */
		m   = 0;
		oid = 0;

		#define simplescan(X,T) \
			basetime = usec();					\
			for (k = 0; k < column.colcount; k++) {	\
				STATS bcomparisons[i] += 1;		\
				if (((T*)column.col)[k] < shigh.X && ((T*)column.col)[k] >= slow.X ) {\
					oids[oid++] = k;\
				}\
			}
		switch(column.coltype){
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
		basetimer[i] = usec() - basetime;
		tuples[i] = oid; /* for result error checking */
		fprintf(devnull, "m = %ld ", m); /* to break compiler optimizations, which may lead to empty loops */

		/* zonesmaps scan */
		m   = 0;
		oid = 0;
		/* watch out, the zones are closed (min,max) intervals */
		#define zonequery(X,T) \
			zonetime = usec();					\
			for (j = 0; j < zonemaps->zmaps_cnt; j++) {		\
				STATS zindex[i] += 1;			\
				if (slow.X <= zonemaps->zmaps[j].min.X && shigh.X > zonemaps->zmaps[j].max.X) { /* all qualify */			\
					for (k = j * rpp, lim = k + rpp < column.colcount ? k+rpp : column.colcount; k < lim; k++) {	\
						oids[oid++] = k;															\
					}																				\
				} else if (!(shigh.X <= zonemaps->zmaps[j].min.X || slow.X > zonemaps->zmaps[j].max.X)) {					\
					/* zone maps are inclusive */													\
					for (k = j * rpp, lim = k + rpp < column.colcount ? k+rpp : column.colcount; k < lim; k++) {	\
						STATS zcomparisons[i] += 1;													\
						if (((T*)column.col)[k] < shigh.X && ((T*)column.col)[k] >= slow.X ) {					\
							oids[oid++] = k;														\
						}\
					} \
				}\
			}

		switch(column.coltype){
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
			printf("%s base %ld zonemap %ld differ\n", column.colname, tuples[i], oid);
		fprintf(devnull," %ld",m); /* to break compiler optimizations */

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
							lim = lim > column.colcount ? column.colcount: lim;				\
							if ((imprints##B[tf] & ~innermask) == 0) {			\
								for (; l < lim; l++) {							\
									oids[oid++] = l;							\
								}												\
							} else {											\
								for (val = ((T*)column.col)[l]; l < lim; l++, val = ((T*)column.col)[l]) {	\
									STATS icomparisons[i] += 1;					\
									if (val < shigh.X && val >= slow.X) {		\
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
						lim = lim > column.colcount ? column.colcount : lim;					\
						if ((imprints##B[tf] & ~innermask) == 0) {				\
							for (; l < lim; l++) {								\
								oids[oid++] = l;								\
							}													\
						} else {												\
							for (val = ((T*)column.col)[l]; l < lim; l++, val = ((T*)column.col)[l]) {	\
								STATS icomparisons[i] += 1;						\
								if (val < shigh.X && val >= slow.X  ) {			\
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
			switch(bins) {						\
			case 8:  impsquery(X,T,8); break;	\
			case 16: impsquery(X,T,16); break;	\
			case 32: impsquery(X,T,32); break;	\
			case 64: impsquery(X,T,64); break;	\
			default: break;						\
			}

		switch(column.coltype){
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
			printf("%s base %ld imprints %ld differ\n", column.colname, tuples[i], oid);
		fprintf(devnull, "m = %ld\n", m); /* to break compiler optimizations */
	}


	for (i =0; i< REPETITION; i++) {
		VERBOSE printf("%s query[%d]=%ld\t selectivity=%2.1f%%\t|\tscan = %ld\timprints = %ld\tzone = %ld\t(usec)\n",
		       column.colname, i, tuples[i], tuples[i]* 100.0/column.colcount, basetimer[i], impstimer[i], zonetimer[i]);
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

/* simulate a series of queries */
void genQueryRange(Column column, int i)
{
	long low, high;
	int j;
	int lastbit;

	mask = 0;
	innermask = 0;

	/* select a range based on the actual non-empty bins */
	for (lastbit = bins-1; lastbit >0; lastbit--)
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
	slow.X = mibins[low].X;						\
	shigh.X = mxbins[high].X;					\
	PRINT_QUERIES {								\
		printf("query             ");			\
		printMask(mask,BITS); putchar('\n');	\
		printf("inner msk         ");			\
		printMask(innermask,BITS);				\
		putchar('\n');							\
	}

	switch (column.coltype) {
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

/*****************************
 * Printing Functions        *
 * ***************************/


/* show the distribution graphically */
void printHistogram(long histo[BITS], char *name)
{
	int i, n;
	double m = 0;

	for (i = 0; i < bins; i++)
		if (histo[i] > m)
			m = histo[i];

	m /= 10.0;
	for (n = 9; n >= 0; n--) {
		printf("                 ");
		for (i = 0; i < bins; i++)
			printf("%c", histo[i] > n * m?'H':' ');
		printf("\n");
	}
}

void printBins(Column column)
{
	int j;
	printf("%s bins ", column.colname);
	for ( j=0; j<bins; j++)
		switch(column.coltype){
		case TYPE_bte:
			printf(" %7d:%d ", mibins[j].bval, mxbins[j].bval);
			break;
		case TYPE_sht:
			printf(" %7d:%d ", mibins[j].sval, mxbins[j].sval);
			break;
		case TYPE_int:
			printf(" %7d:%d ", mibins[j].ival, mxbins[j].ival);
			break;
		case TYPE_lng:
			printf(" %7ld:%ld ", mibins[j].lval, mxbins[j].lval);
			break;
		case TYPE_oid:
			printf(" %7lu:%lu ", mibins[j].ulval, mxbins[j].ulval);
			break;
		case TYPE_flt:
			printf(" %9.8f:%0.8f ", mibins[j].fval, mxbins[j].fval);
			break;
		case TYPE_dbl:
			printf(" %9.8g:%9.8g", mibins[j].dval,  mxbins[j].dval);
		}
	printf("\n");
}

void printMask(long mask, int limit)
{
	int j;
	for ( j =0; j<limit; j++)
		printf("%c", isSet(mask,j)?'x':'.');
}

void printImprint(Column column, Zonemap_index *zonemaps)
{
	int i,j, lzone, blks = 0,tf = 0;
	unsigned long mask;
	ValRecord mx,mi;
	int unique = 0;

	printf("%s rpp=%d imprint cells %ld zone cells %ld \n", column.colname, rpp, dct_cnt, zonemaps->zmaps_cnt);
	printf("                 ");
	for ( j =0; j< bins; j++)
		printf("%c", j% 10 == 0?'0'+ j/10:'.');
	printf("\n");

#define findrange(X) \
	for ( j= lzone+1; j < blks; j++) { \
		if ( zonemaps->zmaps[j].min.X < mi.X) \
			mi.X = zonemaps->zmaps[j].min.X; \
		if ( zonemaps->zmaps[j].max.X > mx.X) \
			mx.X = zonemaps->zmaps[j].max.X; \
	}

	for ( i=0; i< dct_cnt; i++) {
		if (dct_scalar[i].repeated == 0) {
			for (j=0; j<dct_scalar[i].blks;j++) {
				mi = zonemaps->zmaps[blks].min;
				mx = zonemaps->zmaps[blks].max;
				blks ++;
				lzone = blks;
				printf("[ %10d ]   ", blks);
				mask = getMask(tf);
				printMask(mask,bins);
				tf++;
				switch(column.coltype){
					case TYPE_bte:
						printf(" %7d : %d\n", mi.bval,  mx.bval);
					break;
					case TYPE_sht:
						printf(" %7d : %d\n", mi.sval,  mx.sval);
					break;
					case TYPE_int:
						printf(" %7d : %d\n", mi.ival,  mx.ival);
					break;
					case TYPE_lng:
						printf(" %7ld : %ld\n", mi.lval,  mx.lval);
					break;
					case TYPE_oid:
						printf(" %7lu : %lu\n", mi.ulval,  mx.ulval);
					break;
					case TYPE_flt:
						printf(" %9.8f : %9.8f\n", mi.fval,  mx.fval);
					break;
					case TYPE_dbl:
					printf(" %9.8f : %9.8f\n", mi.dval,  mx.dval);
				}
			}
		} else { /* same imprint for imprint[i].blks next blocks */
			unique++;
			lzone = blks;
			mi = zonemaps->zmaps[lzone].min;
			mx = zonemaps->zmaps[lzone].max;
			blks += dct_scalar[i].blks;
			printf("[ %10d ]+  ", blks);
			mask = getMask(tf);
			printMask(mask,bins);
			tf++;
			switch(column.coltype){
			case TYPE_bte:
				findrange(bval);
				printf(" %7d : %d\n", mi.bval,  mx.bval);
				break;
			case TYPE_sht:
				findrange(sval);
				printf(" %7d : %d\n", mi.sval,  mx.sval);
				break;
			case TYPE_int:
				findrange(ival);
				printf(" %7d : %d\n", mi.ival,  mx.ival);
				break;
			case TYPE_lng:
				findrange(lval);
				printf(" %7ld : %ld\n", mi.lval,  mx.lval);
				break;
			case TYPE_oid:
				findrange(ulval);
				printf(" %7lu : %lu\n", mi.ulval,  mx.ulval);
				break;
			case TYPE_flt:
				findrange(fval);
				printf(" %9.8f : %9.8f\n", mi.fval,  mx.fval);
				break;
			case TYPE_dbl:
				findrange(dval);
				printf(" %9.8f : %9.8f\n", mi.dval,  mx.dval);
			}
		}
	}
	printBins(column);

	printf("%s histogram summary ", column.colname);
	i = 0;
	for ( j=0; j< bins; j++)
	if ( histogram[j] == 0)  i++;
		printf(" %d empty cells, bits used %d",i,bins-i);
	printf("\n%s histogram", column.colname);
	for ( j=0; j< bins; j++)
	if( histogram[j])
		printf("[%d] %ld ", j, histogram[j]);
	i=0;
	for ( j=0; j< bins; j++)
		 i += isSet(globalmask,j);
	printf("\n%s vectors dct_cnt %ld imps_cnt %ld unique %d bits %d ", column.colname, dct_cnt, imps_cnt, unique, i);
	for ( j=0; j<= bins; j++)
		if(vectors[j])
			printf("[%d] %ld ", j, vectors[j]);
	printf("\n");
}

void statistics(Column column)
{
	double var = 0, bitvar = 0;
	double delta = 0, bitdelta = 0;
	double mean = 0, bitmean = 0;
	double c;
	long edit, on;
	unsigned long mask;
	int bitcnt, i, j, first, last;

	assert(globalmask);
	assert(imps_cnt);
	for( i= 0; i< imps_cnt; i++){
		first= -1;
		bitcnt = 0;
		mask = getMask(i);
		for( j=0; j< bins; j++)
		if( isSet(mask,j)  && isSet(globalmask,j)){
			if ( first == -1) first= j;
			last = j;
			bitcnt++;	/* number of bits set */
		}
		/* compensate for the rpp */
		assert(first != -1);
		c= (last-first+1.0);
		delta = c-mean;
		mean += c;
		var += (delta*(c - mean/(i+1)));

		bitmean += bitcnt;
		bitdelta += bitcnt-bitmean;
		bitvar += (bitdelta*(bitcnt - bitmean/(i+1)));
	}
	mean /= i;
	bitmean /= i;
	printf("%s bit spread average %5.3f deviation %5.3f bit density average %5.3f dev %5.3f total %d bits rpp %d\n", 
		column.colname, mean, sqrt(var/imps_cnt), bitmean, sqrt(bitvar/imps_cnt), bins, rpp);

	/* edit distance */
	edit = 0; on = 0;
	for (i=0; i< imps_cnt; i++) {
		mask = getMask(i);
		for (j=0; j<bins; j++)
			if (isSet(mask,j)) on++;
		if (i > 0) {
			mask = mask ^ getMask(i-1);
			for (j=0; j<bins; j++)
				if (isSet(mask,j)) edit++;
		}
	}
	printf("%s total bits on %ld total edit distance %ld entropy is %lf\n", 
		column.colname, on, edit, (double)edit/(double)(2*on));
}


