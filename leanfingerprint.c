#include "leanfingerprint.h"

#define STATS if (1)

#define REPETITION  10
#define EXPERIMENTS 1
#define BITRANGE    1    /* search by bin 1 or by sampling the data 0 */
#define WAHTEST  1

char filename[1024];
char column[1024];
char *col, *typename;
int coltype;
unsigned long colcount;
unsigned long pages;
long filesize;
FILE *devnull;
int sorted=1;

int stride[14]= { 0,0,0,1,2,0,4,8,0,0,4,8,8,0};

/* actual number of bins used */
int BINS = BITS;

/* global bounds */
ValRecord min, max, slice; /* 62 bins equidistant in value distribution */
ValRecord absmin, absmax, mxbins[BITS],mibins[BITS];

Zonemap *zmap;
long zonetop;

Fingerprint *cfp;
long *bitmask;
long globalmask;
int cfptop=0;
int masktop=0;

/* each column fingerprint is characterized by a partition vector */

long histogram[BITS]; /* of bin fillers */
/* vector filling distribution */
long vectors[BITS+1];
int HISTO_TYPE=0; /*which type of histogram to build*/

int rpp=1; /* rows per page */

int materialize = 1; /* oid list or aggregate construction */

#if 0
#define setBit(X,Y)  ((((long)1)<< (Y)) | ( ~(((long)1)<< (Y)) & (X)))
#define isSet(X,Y)   (((((long)1)<<(Y)) & (X)) ? 1:0)
#endif

/* WAH encoding over the base table  */
/* this implementation uses the binning derived for cfp */
typedef struct{
	int n;
	int *bits;
} WAHindex ;

long wahtimer;
WAHindex wah[BITS];

#define wahsetBit(X,Y)  ((((int)1)<< (Y)) | ( ~(((int)1)<< (Y)) & (X)))
#define wahisSet(X,Y)   (((((int)1)<<(Y)) & (X)) ? 1:0)
#define wahbit(T,X,Y,T2) val.T2 = ((T*)col)[X]; GETBIT(Y,T2);

int bitcover(long histo[BITS]){
	int i, c=0;
	for( i = 0; i< BINS; i++)
		c += histo[i]>0;
	return c;
}
/* show the distribution graphically */
void printHistogram(long histo[BITS], char *name)
{
	int i,n;
	double m = 0;

	for (i=0;i<BINS; i++)
		if ( histo[i] > m) m = histo[i];

	m /= 10.0;
	for ( n= 9; n >= 0; n --){
		printf("                 ");
		for(i =0; i< BINS ; i++)
		printf("%c", histo[i] > n * m?'H':' ');
		printf("\n");
	}
}

void printBins()
{
	int j;
	printf("%s bins ",column);
	for ( j=0; j<BINS; j++)
		switch(coltype){
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

void printMask(long mask,int limit)
{
	int j, cells =0;
	for ( j =0; j<limit; j++)
		printf("%c", isSet(mask,j)?'x':'.');
}


/* we should calculate the selectivity per cell */
void printFingerprint()
{
	int i,j, lzone, blks = 0,tf = 0;
	unsigned long mask;
	ValRecord mx,mi;
	int unique = 0;

	printf("%s rpp=%d cfp cells %d zone cells %ld \n",column,rpp,cfptop,zonetop);
	printf("                 ");
	for ( j =0; j< BINS; j++)
		printf("%c", j% 10 == 0?'0'+ j/10:'.');
	printf("\n");

#define findrange(X) \
	for ( j= lzone+1; j < blks; j++) { \
		if ( zmap[j].min.X < mi.X) \
			mi.X = zmap[j].min.X; \
		if ( zmap[j].max.X > mx.X) \
			mx.X = zmap[j].max.X; \
	}

	for ( i=0; i< cfptop; i++) {
		if (cfp[i].repeated == 0) {
			for (j=0; j<cfp[i].blks;j++) {
				mi = zmap[blks].min;
				mx = zmap[blks].max;
				blks ++;
				lzone = blks;
				printf("[ %10d ]   ", blks);
				mask = getMask(tf);
				printMask(mask,BINS);
				tf++;
				switch(coltype){
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
		} else { /* same fingerprint for cfp[i].blks next blocks */
			unique++;
			lzone = blks;
			mi = zmap[lzone].min;
			mx = zmap[lzone].max;
			blks += cfp[i].blks;
			printf("[ %10d ]+  ", blks);
			mask = getMask(tf);
			printMask(mask,BINS);
			tf++;
			switch(coltype){
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
	printBins();

	printf("%s histogram summary ",column);
	i = 0;
	for ( j=0; j< BINS; j++)
	if ( histogram[j] == 0)  i++;
		printf(" %d empty cells, bits used %d",i,BINS-i);
	printf("\n%s histogram",column);
	for ( j=0; j< BINS; j++)
	if( histogram[j])
		printf("[%d] %ld ", j, histogram[j]);
	i=0;
	for ( j=0; j< BINS; j++)
		 i += isSet(globalmask,j);
	printf("\n%s vectors cfptop %d masktop %d unique %d bits %d ", column, cfptop, masktop, unique, i);
	for ( j=0; j<= BINS; j++)
		if(vectors[j])
			printf("[%d] %ld ", j, vectors[j]);
	printf("\n");
}

long usec()
{
	static struct timeval tpbase;   /* automatically initialized to 0 */
	struct timeval tp;

	if (tpbase.tv_sec == 0)
		gettimeofday(&tpbase, 0);
	gettimeofday(&tp, 0);
	tp.tv_sec -= tpbase.tv_sec;
	return (long) tp.tv_sec * 1000000 + (long) tp.tv_usec;
}

#define upd(X) \
		if (val.X < min.X) min.X = val.X; \
		if (val.X > max.X) max.X = val.X; \
		if (zmap[zonetop].min.X > val.X) \
			zmap[zonetop].min.X = val.X; \
		if (zmap[zonetop].max.X < val.X) \
			zmap[zonetop].max.X = val.X;

long zone_t0;
void zoneMap()
{
	ValRecord val;
	long i; 
	long t0= usec();
	int new = rpp-1; /*rpp is always power of 2*/

	/* malloc zonemap array */
	zmap = (Zonemap *) malloc (sizeof(Zonemap)*(pages+1));
	memset((char*)zmap, 0, sizeof(Zonemap) * (pages+1));

	zonetop = -1;
	for( i=0; i < colcount ; i ++) {

		if (!(i&new)) {
			zonetop++;
			switch (coltype) {
			case TYPE_bte:
				zmap[zonetop].min.bval = 127;
				zmap[zonetop].max.bval = -127;
				break;
			case TYPE_sht:
				zmap[zonetop].min.sval = 32767;
				zmap[zonetop].max.sval = -32767;
				break;
			case TYPE_int:
				zmap[zonetop].min.ival = INT_MAX;
				zmap[zonetop].max.ival = INT_MIN;
				break;
			case TYPE_lng:
				zmap[zonetop].min.lval = LONG_MAX;
				zmap[zonetop].max.lval = LONG_MIN;
				break;
			case TYPE_oid:
				zmap[zonetop].min.ulval = ULONG_MAX;
				zmap[zonetop].max.ulval = 0;
				break;
			case TYPE_flt:
				zmap[zonetop].min.fval = FLT_MAX;
				zmap[zonetop].max.fval = FLT_MIN;
				break;
			case TYPE_dbl:
				zmap[zonetop].min.dval = DBL_MAX;
				zmap[zonetop].max.dval = -DBL_MAX;
			}
		}

		switch(coltype){
			case TYPE_bte: val.bval = *(char*) (col + i * stride[coltype]);  upd(bval);break;
			case TYPE_sht: val.sval = *(short *) (col + i * stride[coltype]);  upd(sval);break;
			case TYPE_int: val.ival = *(int*) (col + i * stride[coltype]);  upd(ival);break;
			case TYPE_lng: val.lval = *(long*) (col + i * stride[coltype]);  upd(lval);break;
			case TYPE_oid: val.ulval = *(unsigned long*) (col + i * stride[coltype]);  upd(ulval);break;
			case TYPE_flt: val.fval = *(float*) (col + i * stride[coltype]);  upd(fval);break;
			case TYPE_dbl: val.dval = *(double*) (col + i * stride[coltype]); upd(dval); break;
		}
	}
	if ((i-1)% rpp ) zonetop++;
	zonetop++;
	zone_t0=  usec()-t0;
}

/* choose a policy */
/* getBitbin getBitMK */
#define GETBIT getBitbin

#define getBit(Z,X,T) { int z; for( z=0; z<BITS; z++) if ( (T) mibins[z].X <=(T) val.X && (T)mxbins[z].X > (T)val.X) break; Z =z;}

/* martin's take on getbin */
#define SRCH(X, Z, S) \
     S1 = (val.X >= mxbins[Z].X);\
     S2 = -(val.X < mxbins[Z+1].X);\
     Z += (S1+S2) * S;

#define getBitMK(Z,X) \
	{ \
	char S1, S2; \
	Z = 32; \
	/*SRCH(X, Z, 32);*/ \
	SRCH(X, Z, 16); \
	SRCH(X, Z, 8); \
	SRCH(X, Z, 4); \
	SRCH(X, Z, 2); \
	SRCH(X, Z, 1); \
	}

/* assume 64 bits for ultra speed with no divisions:) */
#define getBit2binfor(Z,X,Y) \
		if ((Y)val.X < mxbins[0].X) Z = 0; \
		if ((Y)val.X >= mibins[63].X) Z = 63; \
		if ((Y)val.X >= mibins[31].X && (Y)val.X < mxbins[31].X) Z = 31; \
		if ((Y)val.X >= mibins[32].X) { \
		if ((Y)val.X >= mibins[47].X && (Y)val.X < mxbins[47].X) Z = 47; \
		if ((Y)val.X >= mibins[48].X) {int _s; for(_s=48;_s<63;_s++) if((Y)val.X < mxbins[_s].X) break; Z =_s;} \
		if ((Y)val.X < mxbins[46].X) {int _s; for(_s=46;_s>32;_s--) if((Y)val.X >= mibins[_s].X) break; Z =_s;} \
		} \
		if ((Y)val.X < mxbins[30].X) { \
		if ((Y)val.X >= mibins[15].X && (Y)val.X < mxbins[15].X) Z = 15; \
		if ((Y)val.X >= mibins[16].X) {int _s; for(_s=16;_s<30;_s++) if((Y)val.X < mxbins[_s].X) break; Z =_s;} \
		if ((Y)val.X < mxbins[14].X) {int _s; for(_s=14;_s>0;_s--) if((Y)val.X >= mibins[_s].X) break; Z =_s;} \
		}

#define check(Z,X,W) if (val.X >= mibins[W].X && val.X < mxbins[W].X) Z = W;
#define left(Z,X,W)  if (val.X < mxbins[W].X)
#define right(Z,X,W) if (val.X >= mibins[W].X)

#define getBitbin(Z,X) \
		right(Z,X,32){ \
			right(Z,X,48) { \
				right(Z,X,56) {\
					right(Z,X,60){ \
						right(Z,X,62) {\
							Z = 62;\
							right(Z,X,63) {\
								Z = 63;\
							}\
						}\
						check(Z,X,61); \
						left(Z,X,60) {\
							Z = 60; \
						}\
					}\
					check(Z,X,59);\
					left(Z,X,58) {\
						right(Z,X,58) {\
							Z = 58;\
						}\
						check(Z,X,57);\
						left(Z,X,56) {\
							Z = 56; \
						}\
					}\
				}\
				check(Z,X,55);\
				left(Z,X,54) { \
					right(Z,X,52){ \
						right(Z,X,54) {\
							Z = 54;\
						}\
						check(Z,X,53);\
						left(Z,X,52) {\
							Z = 52; \
						}\
					}\
					check(Z,X,51);\
					left(Z,X,50) {\
						right(Z,X,50) {\
							Z = 50;\
						}\
						check(Z,X,49);\
						left(Z,X,48) {\
							Z = 48; \
						}\
					}\
				}\
			}\
			check(Z,X,47);\
			left(Z,X,46) { \
				right(Z,X,40) {\
					right(Z,X,44){ \
						right(Z,X,46) {\
							Z = 46;\
						}\
						check(Z,X,45); \
						left(Z,X,44) {\
							Z = 44; \
						}\
					}\
					check(Z,X,43);\
					left(Z,X,42) {\
						right(Z,X,42) {\
							Z = 42;\
						}\
						check(Z,X,41);\
						left(Z,X,40) {\
							Z = 40; \
						}\
					}\
				}\
				check(Z,X,39);\
				left(Z,X,38) { \
					right(Z,X,36){ \
						right(Z,X,38) {\
							Z = 38;\
						}\
						check(Z,X,37);\
						left(Z,X,36) {\
							Z = 36; \
						}\
					}\
					check(Z,X,35);\
					left(Z,X,34) {\
						right(Z,X,34) {\
							Z = 34;\
						}\
						check(Z,X,33);\
						left(Z,X,32) {\
							Z = 32; \
						}\
					}\
				}\
			}\
		}\
		check(Z,X,31);\
		left(Z,X,30) { \
			right(Z,X,16) { \
				right(Z,X,24) {\
					right(Z,X,28){ \
						right(Z,X,30) {\
							Z = 30;\
						}\
						check(Z,X,29); \
						left(Z,X,28) {\
							Z = 28; \
						}\
					}\
					check(Z,X,27);\
					left(Z,X,26) {\
						right(Z,X,26) {\
							Z = 26;\
						}\
						check(Z,X,25);\
						left(Z,X,24) {\
							Z = 24; \
						}\
					}\
				}\
				check(Z,X,23);\
				left(Z,X,22) { \
					right(Z,X,20){ \
						right(Z,X,22) {\
							Z = 22;\
						}\
						check(Z,X,21);\
						left(Z,X,20) {\
							Z = 20; \
						}\
					}\
					check(Z,X,19);\
					left(Z,X,18) {\
						right(Z,X,18) {\
							Z = 18;\
						}\
						check(Z,X,17);\
						left(Z,X,16) {\
							Z = 16; \
						}\
					}\
				}\
			}\
			check(Z,X,15);\
			left(Z,X,14) { \
				right(Z,X,8) {\
					right(Z,X,12){ \
						right(Z,X,14) {\
							Z = 14;\
						}\
						check(Z,X,13);\
						left(Z,X,12) {\
							Z = 12; \
						}\
					}\
					check(Z,X,11);\
					left(Z,X,10) {\
						right(Z,X,10) {\
							Z = 10;\
						}\
						check(Z,X,9);\
						left(Z,X,8) {\
							Z = 8; \
						}\
					}\
				}\
				check(Z,X,7);\
				left(Z,X,6) { \
					right(Z,X,4){ \
						right(Z,X,6) {\
							Z = 6;\
						}\
						check(Z,X,5);\
						left(Z,X,4) {\
							Z = 4; \
						}\
					}\
					check(Z,X,3);\
					left(Z,X,2) {\
						right(Z,X,2) {\
							Z = 2;\
						}\
						check(Z,X,1);\
						left(Z,X,0) {\
							Z = 0; \
						}\
					}\
				}\
			}\
		}\

#define getBit2binif(Z,X,Y) \
		if ((Y)val.X < mxbins[0].X) Z = 0; \
		if ((Y)val.X >= mibins[63].X) Z = 63; \
		if ((Y)val.X >= mibins[31].X && (Y)val.X < mxbins[31].X) Z = 31; \
		if ((Y)val.X >= mibins[32].X) { \
			if ((Y)val.X >= mibins[47].X && (Y)val.X < mxbins[47].X) Z = 47; \
			if ((Y)val.X >= mibins[48].X) { \
				if ((Y)val.X < mxbins[48].X) Z = 48; \
				if ((Y)val.X >= mibins[49].X && (Y)val.X < mxbins[49].X) Z = 49; \
				if ((Y)val.X >= mibins[50].X && (Y)val.X < mxbins[50].X) Z = 50; \
				if ((Y)val.X >= mibins[51].X && (Y)val.X < mxbins[51].X) Z = 51; \
				if ((Y)val.X >= mibins[52].X && (Y)val.X < mxbins[52].X) Z = 52; \
				if ((Y)val.X >= mibins[53].X && (Y)val.X < mxbins[53].X) Z = 53; \
				if ((Y)val.X >= mibins[54].X && (Y)val.X < mxbins[54].X) Z = 54; \
				if ((Y)val.X >= mibins[55].X && (Y)val.X < mxbins[55].X) Z = 55; \
				if ((Y)val.X >= mibins[56].X && (Y)val.X < mxbins[56].X) Z = 56; \
				if ((Y)val.X >= mibins[57].X && (Y)val.X < mxbins[57].X) Z = 57; \
				if ((Y)val.X >= mibins[58].X && (Y)val.X < mxbins[58].X) Z = 58; \
				if ((Y)val.X >= mibins[59].X && (Y)val.X < mxbins[59].X) Z = 59; \
				if ((Y)val.X >= mibins[60].X && (Y)val.X < mxbins[60].X) Z = 60; \
				if ((Y)val.X >= mibins[61].X && (Y)val.X < mxbins[61].X) Z = 61; \
				if ((Y)val.X >= mibins[62].X && (Y)val.X < mxbins[62].X) Z = 62; \
			} \
			if ((Y)val.X < mxbins[46].X) {\
				if ((Y)val.X < mxbins[32].X) Z = 32; \
				if ((Y)val.X >= mibins[33].X && (Y)val.X < mxbins[33].X) Z = 33; \
				if ((Y)val.X >= mibins[34].X && (Y)val.X < mxbins[34].X) Z = 34; \
				if ((Y)val.X >= mibins[35].X && (Y)val.X < mxbins[35].X) Z = 35; \
				if ((Y)val.X >= mibins[36].X && (Y)val.X < mxbins[36].X) Z = 36; \
				if ((Y)val.X >= mibins[37].X && (Y)val.X < mxbins[37].X) Z = 37; \
				if ((Y)val.X >= mibins[38].X && (Y)val.X < mxbins[38].X) Z = 38; \
				if ((Y)val.X >= mibins[39].X && (Y)val.X < mxbins[39].X) Z = 39; \
				if ((Y)val.X >= mibins[40].X && (Y)val.X < mxbins[40].X) Z = 40; \
				if ((Y)val.X >= mibins[41].X && (Y)val.X < mxbins[41].X) Z = 41; \
				if ((Y)val.X >= mibins[42].X && (Y)val.X < mxbins[42].X) Z = 42; \
				if ((Y)val.X >= mibins[43].X && (Y)val.X < mxbins[43].X) Z = 43; \
				if ((Y)val.X >= mibins[44].X && (Y)val.X < mxbins[44].X) Z = 44; \
				if ((Y)val.X >= mibins[45].X && (Y)val.X < mxbins[45].X) Z = 45; \
				if ((Y)val.X >= mibins[46].X) Z = 46; \
			} \
		} \
		if ((Y)val.X < mxbins[30].X) { \
			if ((Y)val.X >= mibins[15].X && (Y)val.X < mxbins[15].X) Z = 15; \
			if ((Y)val.X >= mibins[16].X) {\
				if ((Y)val.X < mxbins[16].X) Z = 16; \
				if ((Y)val.X >= mibins[17].X && (Y)val.X < mxbins[17].X) Z = 17; \
				if ((Y)val.X >= mibins[18].X && (Y)val.X < mxbins[18].X) Z = 18; \
				if ((Y)val.X >= mibins[19].X && (Y)val.X < mxbins[19].X) Z = 19; \
				if ((Y)val.X >= mibins[20].X && (Y)val.X < mxbins[20].X) Z = 20; \
				if ((Y)val.X >= mibins[21].X && (Y)val.X < mxbins[21].X) Z = 21; \
				if ((Y)val.X >= mibins[22].X && (Y)val.X < mxbins[22].X) Z = 22; \
				if ((Y)val.X >= mibins[23].X && (Y)val.X < mxbins[23].X) Z = 23; \
				if ((Y)val.X >= mibins[24].X && (Y)val.X < mxbins[24].X) Z = 24; \
				if ((Y)val.X >= mibins[25].X && (Y)val.X < mxbins[25].X) Z = 25; \
				if ((Y)val.X >= mibins[26].X && (Y)val.X < mxbins[26].X) Z = 26; \
				if ((Y)val.X >= mibins[27].X && (Y)val.X < mxbins[27].X) Z = 27; \
				if ((Y)val.X >= mibins[28].X && (Y)val.X < mxbins[28].X) Z = 28; \
				if ((Y)val.X >= mibins[29].X && (Y)val.X < mxbins[29].X) Z = 29; \
				if ((Y)val.X >= mibins[30].X) Z = 30; \
			} \
			if ((Y)val.X < mxbins[14].X) { \
				if ((Y)val.X >= mibins[1].X && (Y)val.X < mxbins[1].X) Z = 1; \
				if ((Y)val.X >= mibins[2].X && (Y)val.X < mxbins[2].X) Z = 2; \
				if ((Y)val.X >= mibins[3].X && (Y)val.X < mxbins[3].X) Z = 3; \
				if ((Y)val.X >= mibins[4].X && (Y)val.X < mxbins[4].X) Z = 4; \
				if ((Y)val.X >= mibins[5].X && (Y)val.X < mxbins[5].X) Z = 5; \
				if ((Y)val.X >= mibins[6].X && (Y)val.X < mxbins[6].X) Z = 6; \
				if ((Y)val.X >= mibins[7].X && (Y)val.X < mxbins[7].X) Z = 7; \
				if ((Y)val.X >= mibins[8].X && (Y)val.X < mxbins[8].X) Z = 8; \
				if ((Y)val.X >= mibins[9].X && (Y)val.X < mxbins[9].X) Z = 9; \
				if ((Y)val.X >= mibins[10].X && (Y)val.X < mxbins[10].X) Z = 10; \
				if ((Y)val.X >= mibins[11].X && (Y)val.X < mxbins[11].X) Z = 11; \
				if ((Y)val.X >= mibins[12].X && (Y)val.X < mxbins[12].X) Z = 12; \
				if ((Y)val.X >= mibins[13].X && (Y)val.X < mxbins[13].X) Z = 13; \
				if ((Y)val.X >= mibins[14].X) Z = 14; \
			} \
		}


#define getBitif(Z,X,Y) \
if ((Y)val.X < mxbins[0].X) Z = 0; \
if ((Y)val.X >= mibins[63].X) Z = 63; \
if ((Y)val.X >= mibins[1].X && (Y)val.X < mxbins[1].X) Z = 1; \
if ((Y)val.X >= mibins[2].X && (Y)val.X < mxbins[2].X) Z = 2; \
if ((Y)val.X >= mibins[3].X && (Y)val.X < mxbins[3].X) Z = 3; \
if ((Y)val.X >= mibins[4].X && (Y)val.X < mxbins[4].X) Z = 4; \
if ((Y)val.X >= mibins[5].X && (Y)val.X < mxbins[5].X) Z = 5; \
if ((Y)val.X >= mibins[6].X && (Y)val.X < mxbins[6].X) Z = 6; \
if ((Y)val.X >= mibins[7].X && (Y)val.X < mxbins[7].X) Z = 7; \
if ((Y)val.X >= mibins[8].X && (Y)val.X < mxbins[8].X) Z = 8; \
if ((Y)val.X >= mibins[9].X && (Y)val.X < mxbins[9].X) Z = 9; \
if ((Y)val.X >= mibins[10].X && (Y)val.X < mxbins[10].X) Z = 10; \
if ((Y)val.X >= mibins[11].X && (Y)val.X < mxbins[11].X) Z = 11; \
if ((Y)val.X >= mibins[12].X && (Y)val.X < mxbins[12].X) Z = 12; \
if ((Y)val.X >= mibins[13].X && (Y)val.X < mxbins[13].X) Z = 13; \
if ((Y)val.X >= mibins[14].X && (Y)val.X < mxbins[14].X) Z = 14; \
if ((Y)val.X >= mibins[15].X && (Y)val.X < mxbins[15].X) Z = 15; \
if ((Y)val.X >= mibins[16].X && (Y)val.X < mxbins[16].X) Z = 16; \
if ((Y)val.X >= mibins[17].X && (Y)val.X < mxbins[17].X) Z = 17; \
if ((Y)val.X >= mibins[18].X && (Y)val.X < mxbins[18].X) Z = 18; \
if ((Y)val.X >= mibins[19].X && (Y)val.X < mxbins[19].X) Z = 19; \
if ((Y)val.X >= mibins[20].X && (Y)val.X < mxbins[20].X) Z = 20; \
if ((Y)val.X >= mibins[21].X && (Y)val.X < mxbins[21].X) Z = 21; \
if ((Y)val.X >= mibins[22].X && (Y)val.X < mxbins[22].X) Z = 22; \
if ((Y)val.X >= mibins[23].X && (Y)val.X < mxbins[23].X) Z = 23; \
if ((Y)val.X >= mibins[24].X && (Y)val.X < mxbins[24].X) Z = 24; \
if ((Y)val.X >= mibins[25].X && (Y)val.X < mxbins[25].X) Z = 25; \
if ((Y)val.X >= mibins[26].X && (Y)val.X < mxbins[26].X) Z = 26; \
if ((Y)val.X >= mibins[27].X && (Y)val.X < mxbins[27].X) Z = 27; \
if ((Y)val.X >= mibins[28].X && (Y)val.X < mxbins[28].X) Z = 28; \
if ((Y)val.X >= mibins[29].X && (Y)val.X < mxbins[29].X) Z = 29; \
if ((Y)val.X >= mibins[30].X && (Y)val.X < mxbins[30].X) Z = 30; \
if ((Y)val.X >= mibins[31].X && (Y)val.X < mxbins[31].X) Z = 31; \
if ((Y)val.X >= mibins[32].X && (Y)val.X < mxbins[32].X) Z = 32; \
if ((Y)val.X >= mibins[33].X && (Y)val.X < mxbins[33].X) Z = 33; \
if ((Y)val.X >= mibins[34].X && (Y)val.X < mxbins[34].X) Z = 34; \
if ((Y)val.X >= mibins[35].X && (Y)val.X < mxbins[35].X) Z = 35; \
if ((Y)val.X >= mibins[36].X && (Y)val.X < mxbins[36].X) Z = 36; \
if ((Y)val.X >= mibins[37].X && (Y)val.X < mxbins[37].X) Z = 37; \
if ((Y)val.X >= mibins[38].X && (Y)val.X < mxbins[38].X) Z = 38; \
if ((Y)val.X >= mibins[39].X && (Y)val.X < mxbins[39].X) Z = 39; \
if ((Y)val.X >= mibins[40].X && (Y)val.X < mxbins[40].X) Z = 40; \
if ((Y)val.X >= mibins[41].X && (Y)val.X < mxbins[41].X) Z = 41; \
if ((Y)val.X >= mibins[42].X && (Y)val.X < mxbins[42].X) Z = 42; \
if ((Y)val.X >= mibins[43].X && (Y)val.X < mxbins[43].X) Z = 43; \
if ((Y)val.X >= mibins[44].X && (Y)val.X < mxbins[44].X) Z = 44; \
if ((Y)val.X >= mibins[45].X && (Y)val.X < mxbins[45].X) Z = 45; \
if ((Y)val.X >= mibins[46].X && (Y)val.X < mxbins[46].X) Z = 46; \
if ((Y)val.X >= mibins[47].X && (Y)val.X < mxbins[47].X) Z = 47; \
if ((Y)val.X >= mibins[48].X && (Y)val.X < mxbins[48].X) Z = 48; \
if ((Y)val.X >= mibins[49].X && (Y)val.X < mxbins[49].X) Z = 49; \
if ((Y)val.X >= mibins[50].X && (Y)val.X < mxbins[50].X) Z = 50; \
if ((Y)val.X >= mibins[51].X && (Y)val.X < mxbins[51].X) Z = 51; \
if ((Y)val.X >= mibins[52].X && (Y)val.X < mxbins[52].X) Z = 52; \
if ((Y)val.X >= mibins[53].X && (Y)val.X < mxbins[53].X) Z = 53; \
if ((Y)val.X >= mibins[54].X && (Y)val.X < mxbins[54].X) Z = 54; \
if ((Y)val.X >= mibins[55].X && (Y)val.X < mxbins[55].X) Z = 55; \
if ((Y)val.X >= mibins[56].X && (Y)val.X < mxbins[56].X) Z = 56; \
if ((Y)val.X >= mibins[57].X && (Y)val.X < mxbins[57].X) Z = 57; \
if ((Y)val.X >= mibins[58].X && (Y)val.X < mxbins[58].X) Z = 58; \
if ((Y)val.X >= mibins[59].X && (Y)val.X < mxbins[59].X) Z = 59; \
if ((Y)val.X >= mibins[60].X && (Y)val.X < mxbins[60].X) Z = 60; \
if ((Y)val.X >= mibins[61].X && (Y)val.X < mxbins[61].X) Z = 61; \
if ((Y)val.X >= mibins[62].X && (Y)val.X < mxbins[62].X) Z = 62; \


#define CMPVAL(X) \
	if ( v1.X < v2.X ) return -1; \
	if ( v1.X > v2.X ) return 1; \
	return 0;

static int
cmpvalues(const void *p1, const void *p2)
{
	ValRecord v1 = *(ValRecord*)p1;
	ValRecord v2 = *(ValRecord*)p2;
	switch(coltype){
	case TYPE_bte: CMPVAL(bval);
	case TYPE_sht: CMPVAL(sval);
	case TYPE_int: CMPVAL(ival);
	case TYPE_lng: CMPVAL(lval);
	case TYPE_oid: CMPVAL(ulval);
	case TYPE_flt: CMPVAL(fval);
	case TYPE_dbl: CMPVAL(dval);
	}
	return 0;
}

/* simulate a series of queries */
/* to make it realistic, we better clean memory before each run, now we measure hot */
/* in each case we are interested in both number of access and the total time. */
ValRecord slow, shigh;
unsigned long mask;
unsigned long innermask;
void genQueryRange(int i, int flag)
{
	long low, high;
	int j;
	int lastbit;
	mask = 0;
	innermask = 0;

/* select a range based on the actual non-empty bins */
	if (flag) {
		for (lastbit = BINS-1; lastbit >0; lastbit--)
			if (isSet(globalmask, lastbit))
				break;

		low = 0 + (int)(rand() * 1.0 / RAND_MAX * lastbit);
		high = low + i * lastbit/ (100.0/REPETITION);
		if (high > lastbit)
			high = lastbit;
		if ( low > high)
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
			mask = setBit(mask,j);
		}
		/* inner mask should only be set when range bounds call for it */
		/* in this case we use full bins */
		for (j = low+1; j <= high-1; j++) {
			innermask = setBit(innermask,j);
		}

#define setqueryrange(X) \
		slow.X = mibins[low].X; \
		shigh.X = mxbins[high].X; \
		printf("query             "); printMask(mask,BITS); putchar('\n'); \
		printf("inner msk         "); printMask(innermask,BITS); putchar('\n');

		switch(coltype){
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

/* random probing the source table leads to skewed searches*/
	low = (int)(rand() * 1.0 / RAND_MAX * colcount);
	high = (int)(rand() * 1.0 / RAND_MAX * colcount);

#define setquerybound(X,T) \
		if ( ((T*)col)[low] > ((T*)col)[high]){ j= low; low=high; high=j; } \
		slow.X = ((T*)col)[low]; \
		shigh.X = ((T*)col)[high]; \
		/* now create the bit mask */ \
		mask=0; for (j = 0; j < BINS; j++) \
		if (!(((T*)col)[low] >= mxbins[j].X|| ((T*)col)[high] < mibins[j].X)) \
			mask = setBit(mask,j);\
		printf("query        "); printMask(mask,BINS); putchar('\n');

	switch(coltype){
	case TYPE_bte:
		setquerybound(bval,char);
		break;
	case TYPE_sht:
		setquerybound(sval,short);
		break;
	case TYPE_int:
		setquerybound(ival,int);
		break;
	case TYPE_lng:
		setquerybound(lval,long);
		break;
	case TYPE_oid:
		setquerybound(ulval,unsigned long);
		break;
	case TYPE_flt:
		setquerybound(fval,float);
		break;
	case TYPE_dbl:
		setquerybound(dval,double);
	}
	/* assume that none of the bounds is included */
	/* go for expensive boundary test if needed */
	for (j = low+1; j <= high-1; j++) {
		innermask = setBit(innermask,j);
	}
}

void queries()
{
	int e,i,k, n;
	long m = 0, tf;
	long tuples[REPETITION];
	long j,l,lim,low = 0, high = 0, total, bit, basetime,cfptime, zonetime, wahtime;
	long basetimer[REPETITION], cfptimer[REPETITION], zonetimer[REPETITION], wahtimer[REPETITION];
	long *oids, oid =0;
	unsigned char *merge;
	unsigned char *bitmask8 = (unsigned char *)bitmask;
	unsigned short *bitmask16 = (unsigned short *)bitmask;
	unsigned int *bitmask32 = (unsigned int *) bitmask;
	unsigned long *bitmask64 = (unsigned long *) bitmask;
	long bindex[REPETITION], bcomparisons[REPETITION];
	long zindex[REPETITION], zcomparisons[REPETITION];
	long findex[REPETITION], fcomparisons[REPETITION];
	long windex[REPETITION], wcomparisons[REPETITION];

	if (materialize){
		oids = (long *) malloc(colcount * sizeof(long));
		merge = (unsigned char *) malloc(colcount/8+1);
	}

	(void) wahtime;
	(void) wahtimer;
	for (e=0; e < EXPERIMENTS; e++) {
		oid = 0; m=0;
		for (i =0; i< REPETITION; i++)
			tuples[i]= basetimer[i] = cfptimer[i] = zonetimer[i] = wahtimer[i] = bindex[i] = bcomparisons[i] = zindex[i] = zcomparisons[i] =  findex[i] = fcomparisons[i] =  windex[i] = wcomparisons[i] = 0;

		for (i= 0; i < REPETITION; i++) {
				/* select a random range from the pool, leads to bias to skew data distribution */
				/* use [slow,shigh) range expression */
				genQueryRange(i,BITRANGE);

				/* simple scan */
				m = 0;
				oid =0;

	#define simplescan(F,X,T) \
					/*printf("zone scan [%"F" -  %"F")\n",slow.X, shigh.X);*/\
					basetime = usec(); /* start timer after printf */ \
					if (materialize) {\
						for (k = 0; k<colcount; k++) { \
							STATS bcomparisons[i] += 1; \
							if (((T*)col)[k] < shigh.X && ((T*)col)[k] >= slow.X )\
								{oids[oid++]=k;}\
						} \
					} else {\
						for (k = 0; k< colcount; k++)\
						if (((T*)col)[k] < shigh.X && ((T*)col)[k] >= slow.X )\
							{oid++; m += ((T*)col)[k] ;}\
					}
				switch(coltype){
				case TYPE_bte:
					simplescan("d",bval,char);
					break;
				case TYPE_sht:
					simplescan("d",sval,short);
					break;
				case TYPE_int:
					simplescan("d",ival,int);
					break;
				case TYPE_lng:
					simplescan("ld",lval,long);
					break;
				case TYPE_oid:
					simplescan("lu",ulval,unsigned long);
					break;
				case TYPE_flt:
					simplescan("f",fval,float);
					break;
				case TYPE_dbl:
					simplescan("g",dval,double);
					break;
				}
				basetime = usec()-basetime;
				basetimer[i] += basetime;
				tuples[i] = oid;
				fprintf(devnull,"m = %ld ", m); /* to break compiler optimizations, which may lead to empty loops */

				/* zone map filter */
				m = 0;
				oid=0;
				/* watch out, the zones are closed (min,max) intervals */
	#define zonequery(X,T)\
			zonetime = usec(); \
			if (materialize) { \
				for (j=0; j< zonetop; j++) { \
					STATS zindex[i] += 1; \
					if (slow.X <= zmap[j].min.X && shigh.X > zmap[j].max.X) { /* all qualify */ \
						for (k = j * rpp, lim = k + rpp < colcount? k+rpp:colcount; k< lim; k++) {\
							oids[oid++] = k;\
						}\
					} else if ( ! (shigh.X <= zmap[j].min.X || slow.X > zmap[j].max.X)) { \
						/* zone maps are inclusive */ \
						for ( k = j * rpp, lim = k + rpp < colcount? k+rpp:colcount; k< lim; k++) { \
							STATS zcomparisons[i] += 1; \
							if (((T*)col)[k] < shigh.X && ((T*)col)[k] >= slow.X ) {\
								oids[oid++]= k; \
							}\
						} \
					}\
				} \
			} else {\
				for (j=0; j< zonetop; j++) { \
					if (slow.X <= zmap[j].min.X && shigh.X > zmap[j].max.X) { /* all qualify */ \
						for ( k = j * rpp, lim = k + rpp < colcount? k+rpp:colcount; k< lim; k++) \
							{m += ((T*)col)[k];} \
					} else if ( ! (shigh.X <= zmap[j].min.X || slow.X > zmap[j].max.X)) { \
						/* zone maps are inclusive */ \
						for (k = j * rpp, lim = k + rpp < colcount? k+rpp:colcount; k< lim; k++) \
							if (((T*)col)[k]< shigh.X && ((T*)col)[k]>= slow.X) \
								{m += ((T*)col)[k];} \
					} \
				}\
			}

				switch(coltype){
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

				zonetime = usec()-zonetime;
				zonetimer[i] += zonetime;
				if (tuples[i] != oid)
					printf("%s base %ld zonemap %ld differ\n", column, tuples[i], oid);
				fprintf(devnull," %ld",m);          /* to break compiler optimizations */

				/* construct filter mask, take care of overflow */
				bit = 1;
				/* column fingerprint filter */
				n = 0;
				tf = 0;
				m = 0;
				oid=0;
	#define cfpquery(X,T,B) \
			cfptime = usec(); \
			if (materialize) {\
				for (j =0; j< cfptop; j++) {\
					if(cfp[j].repeated == 0 ){\
						for( k=tf + cfp[j].blks; tf < k; n++, tf++) { \
							STATS findex[i] += 1; \
							if ( bitmask##B[tf] & mask) { \
								register T val; \
								l= n * rpp; \
								lim = l+ rpp; \
								lim = lim > colcount? colcount: lim;\
								if ( (bitmask##B[tf] & ~innermask) == 0) { \
									for(; l < lim; l++) \
									{oids[oid++]= l;}\
								} else { \
									for (val = ((T*)col)[l] ; l < lim; l++, val = ((T*)col)[l]) { \
										STATS fcomparisons[i] += 1; \
										if ( val< shigh.X && val >= slow.X  ) { \
											oids[oid++]= l; \
										}\
									} \
								} \
							}\
						}\
					} else { /* repeated mask case */\
						STATS findex[i] += 1; \
						if ( bitmask##B[tf] & mask) { \
							register T val;\
							l= n * rpp; \
							lim = l+ rpp * cfp[j].blks; \
							lim = lim > colcount? colcount: lim;\
							if ( (bitmask##B[tf] & ~innermask) == 0) { \
								for (; l < lim; l++) { \
									oids[oid++]= l;\
								}\
							} else { \
								for (val = ((T*)col)[l] ; l < lim; l++, val = ((T*)col)[l]) {\
									STATS fcomparisons[i] += 1; \
									if ( val< shigh.X && val >= slow.X  ) {\
										oids[oid++]= l;\
									}\
								}\
							}\
						}\
						n += cfp[j].blks;\
						tf ++; \
					}\
				} \
			} /* non-materialized case to be done */

	#define binchoose(X,T) \
		switch(BINS) {\
			case 8: cfpquery(X,T,8); break; \
			case 16: cfpquery(X,T,16); break; \
			case 32: cfpquery(X,T,32); break; \
			case 64: cfpquery(X,T,64); break; \
			default: break; \
			}

				switch(coltype){
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
				cfptime = usec()- cfptime;
				cfptimer[i] += cfptime;
				if ( tuples[i] != oid)
					printf("%s base %ld cfp %ld differ\n", column, tuples[i], oid);

#ifdef WAHTEST
	/* the wah compressed bitvectors should be decoded before they can be used */
	/* we assume a limited answer set which is (not?) SORTED in the end */
	/* detection of false positives is the more expensIVE PArt, it requires a scan
	   through the boundaries */
	oid=0;
	#define wahquery(X,T) \
			wahtime = usec(); \
			if (materialize) {\
				int idx, n, m; \
				/* use innermask for direct construction */ \
				/* use an accumulator vector */ \
				memset((char*) merge,0,colcount/8+1);\
				oid = 0;\
				for ( idx = 0; idx < BINS; idx++) {\
					if (isSet(innermask,idx)) { \
						n = 0; \
						for ( k = 0; k < wah[idx].n; k++) { \
							if ( wahisSet(wah[idx].bits[k], 31) ) { /* literal  */\
								j= n;\
								for (m=0 ; j < n+31 && j < colcount; j++, m++) {\
									STATS windex[i] += 1; \
									if ( wahisSet(wah[idx].bits[k],m) ) {\
										merge[j/8] = wahsetBit(merge[j/8],j&7);\
									}\
								}\
								n+= m; \
							} else {  /* fill word */\
								j= n; n+= (~(1<<30)) & wah[idx].bits[k];\
								STATS windex[i] += 1; \
								if ( wahisSet(wah[idx].bits[k],30) ) { \
									for ( ; j < n; j++) {\
										merge[j/8] = wahsetBit(merge[j/8],j&7);\
									}\
								}\
							} \
						}\
					} else if ( isSet(mask,idx)) {\
						n = 0; \
						for ( k = 0; k < wah[idx].n; k++) { \
							if ( wahisSet(wah[idx].bits[k], 31) ) {\
								register T val; /* literal  */\
								j= n;\
								for (m=0 ; j < n+31 && j < colcount; j++, m++) { \
									STATS windex[i] += 1; \
									if ( wahisSet(wah[idx].bits[k],m) ) { \
										val = ((T*)col)[j]; \
										STATS wcomparisons[i] +=1; \
										if ( val < shigh.X && val >= slow.X  ) {\
											merge[j/8] = wahsetBit(merge[j/8],j&7);\
										}\
									}\
								} \
								n += m; \
							} else {\
								register T val; /* fill word */\
								j= n; n+= (~(1<<30)) & wah[idx].bits[k];\
								STATS windex[i] += 1; \
								if ( wahisSet(wah[idx].bits[k],30) ) {\
									for ( ; j < n && j < colcount; j++) {\
										val = ((T*)col)[j]; \
										STATS wcomparisons[i] +=1; \
										if ( val < shigh.X && val >= slow.X  ) {\
											merge[j/8] = wahsetBit(merge[j/8],j&7);\
										}\
									}\
								}\
							}\
						}\
					}\
				}\
				oid = 0; /* materialize the oid list */ \
				for( j = 0; j<(colcount/8)+1; j++) {\
					if (merge[j] != 0)\
						for (m=0; m<8; m++) {oids[oid]=8*j+m; oid += wahisSet(merge[j],m);} \
				}\
			} /* non-materialized case to be done */

	#define wahchoose(X,T) \
		switch(BINS) {\
			case 8: wahquery(X,T); break; \
			case 16: wahquery(X,T); break; \
			case 32: wahquery(X,T); break; \
			case 64: wahquery(X,T); break; \
			default: break; \
			}
				/* Use WAH query */
				switch(coltype){
				case TYPE_bte:
					wahchoose(bval,char);
					break;
				case TYPE_sht:
					wahchoose(sval,short);
					break;
				case TYPE_int:
					wahchoose(ival,int);
					break;
				case TYPE_lng:
					wahchoose(lval,long);
					break;
				case TYPE_oid:
					wahchoose(ulval,unsigned long);
					break;
				case TYPE_flt:
					wahchoose(fval,float);
					break;
				case TYPE_dbl:
					wahchoose(dval,double);
				}
				wahtime = usec()- wahtime;
				wahtimer[i] += wahtime;
				if ( tuples[i] != oid)
					printf("%s base %ld wah %ld differ\n", column, tuples[i], oid);
#endif

				/* we should also check the performance in multithreaded execution */
				fprintf(devnull,"m = %ld\n",m);	/* to break compiler optimizations */
		}

		for (i =0; i< REPETITION; i++) {
				printf("%s %s sorted %d select[%d] cfp %d zones %ld time %ld cfp %ld zone %ld wah %ld tuples %ld %2.1f %%\n",
					column, typename, sorted, i, cfptop, zonetop, 
					basetimer[i],cfptimer[i],zonetimer[i], wahtimer[i], tuples[i], tuples[i]* 100.0/colcount);

				STATS printf ("bindex %ld bcomparisons %ld zindex %ld zcomparisons %ld findex %ld fcomparisons %ld windex %ld wcomparisons %ld\n", bindex[i], bcomparisons[i], zindex[i], zcomparisons[i], findex[i], fcomparisons[i], windex[i], wcomparisons[i] );

				if (i) {
					basetimer[0] += basetimer[i];
					cfptimer[0] += cfptimer[i];
					zonetimer[0] += zonetimer[i];
					wahtimer[0] += wahtimer[i];
				}
		}
	}
	if (materialize)
		free(oids);
}

#define SAMPLE 2048

#define sampleDistribution(X,T)\
	/* draw the sample, */ \
	/* we can use the sample.c from monet - faster */ \
	for ( i = 0; i < SAMPLE; i++){ \
		j = (rand() * colcount)/RAND_MAX ; \
		sample[i].X = ((T*)col)[j]; \
	} \
	/* sort the sample */ \
	qsort((char*) sample, SAMPLE, sizeof(ValRecord), cmpvalues); \
	/* the count[] seems out of sync but it is not used anywhere either */ \
	/* if sample[0]=sample[1] then count[0]=1 and count[1]++ and so on */ \
	j = 0;\
	count[0]=1;\
	for(i=1; i<SAMPLE; i++)\
		if( sample[i].X != sample[j].X) {\
			sample[++j]= sample[i];\
			count[j]=1;\
		} else count[j]++;\
	j += 1; \


int eqHeightHisto(ValRecord *sample, int *count, int smp) {
	int height;
	int k,cnt,j;

	height = SAMPLE/(smp < BITS-2?smp:(BITS-2));

	mibins[0] = absmin;
	mxbins[0] = sample[0];
	mibins[BITS-1] = sample[smp-1];
	mxbins[BITS-1] = absmax;

	if( smp < BITS-1 ) {
		for(k=1; k < smp; k++) {
			mibins[k] = mxbins[k-1];
			mxbins[k] = sample[k];
		}
		for (; k < BITS-1; k++) {
			mibins[k] = sample[smp-1];
			mxbins[k] = absmax;
		}
	} else {
		for (k = 1, j = 0; k < BITS-1; k++) {
			mibins[k] = mxbins[k-1];
			cnt = count[j++];
			while ((j < smp) && (cnt < height)) cnt += count[j++];
			mxbins[k] = sample[j];
			if (j==smp) break;
		}
	}

	return 1;
}

int eqWidthHisto(ValRecord *sample, int smp) {
	int k;

	mibins[0] = absmin;
	mxbins[0] = sample[0];
	mibins[BITS-1] = sample[smp-1];
	mxbins[BITS-1] = absmax;

	if( smp < BITS-1 ) {
		for(k=1; k < smp; k++) {
			mibins[k] = mxbins[k-1];
			mxbins[k] = sample[k];
		}
		if (k<8) BINS=8;
		if (8<=k && k<16) BINS=16;
		if (16<=k && k<32) BINS=32;
		if (32<=k && k<64) BINS=64;
		for (; k < BITS-1; k++) {
			mibins[k] = sample[smp-1];
			mxbins[k] = absmax;
		}
	} else {
		double y , ystep = (double)smp/(double)(BITS-2);
		for (k=1, y=ystep; y < smp; y+= ystep, k++) {
			mibins[k] = mxbins[k-1];
			mxbins[k] = sample[(int)y];
		}
		if (k==BITS-2) { /* there is some leftovers */
			assert (y>=smp);
			assert ((y-ystep)<smp);
			mibins[k] = mxbins[k-1];
			mxbins[k] = sample[smp-1];
		}
	}
	printf("sample gives %d unique values from %d in fingerprints of %d bits\n",smp,SAMPLE,BINS);
	return 1;
}

int cfpSample()
{
	ValRecord sample[SAMPLE];
	int count[SAMPLE];
	int i,j=0,k;

	switch(coltype){
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
		sampleDistribution(lval,long);
		break;
	case TYPE_oid: 
		sampleDistribution(ulval,unsigned long);
		break;
	case TYPE_flt: 
		sampleDistribution(fval,float);
		break;
	case TYPE_dbl: 
		sampleDistribution(dval, double);
	}

	if (HISTO_TYPE == EQWIDTH)
		return eqWidthHisto(sample, j);
	else if (HISTO_TYPE == EQHEIGHT)
		return eqHeightHisto(sample,count,j);
	else
		return -1;
}

void
stats(long timer)
{
	long i, tf=0, k;
	int j,bits;
	unsigned long mask;

	for (i=0; i<BINS; i++)
		vectors[i]= histogram[i] = 0;

	for (i=0; i<cfptop; i++) {
		for (k=0;k<cfp[i].blks;k++) {
			if (cfp[i].repeated == 1) k = cfp[i].blks;
			bits = 0;
			mask = getMask(tf);
			globalmask |= mask;
			for (j=0;j<BINS;j++) {
				if (isSet(mask,j)) {
					bits++;
					histogram[j]++;
				}
			}
			vectors[bits]++;
			tf++;
		}
	}

	//printf("global mask      ");
	//printMask(globalmask,BITS);
	//printf("\n");
	printf("%s cfp size %ld creation time %ld  %ld usec per thousand\n", column, colcount, timer, ((long)timer*1000)/colcount);
	printf("%s zonemap size %ld creation time %ld %ld usec per thousand \n", column, colcount,  zone_t0, ((long)zone_t0*1000)/colcount);
	/* comment uncomment for printing the fingerprints */
	printHistogram(histogram, "Value distribution ");
	printFingerprint();

}

void
cfpMask()
{
	long i,mask,timer,prevmask;
	int j, bit;
	ValRecord val;
	int fits;
	int new;

	/* sample to create the bins in the histogram */
	cfpSample();

	/* how many mask vectors fit in the fingerprint */
	fits = BITS/BINS;
	cfptop = 0;
	masktop = 0;
	prevmask = mask = 0;
	new = rpp-1; /*rpp is always power of 2*/

	cfp = (Fingerprint *) malloc (sizeof(Fingerprint)*pages);
	bitmask = (long *) malloc (sizeof(long)*(pages/fits+1));

	/* init bitmask, not counting on creation time since it can be avoided */
	for (int i=0,masktop=pages/fits+1;i<masktop;i++) {
		bitmask[i]=0;
	}
	masktop=0;

	/* start creation */
	timer = usec();
	for (i=0; i < colcount; i ++) {
		if (!(i&new) && i>0) {
			/* compress list */
			if (prevmask == mask && cfp[cfptop-1].blks < ((1<<MAXOFFSET)-1)) {
				if (cfp[cfptop - 1].repeated == 0) {
					if (cfp[cfptop - 1].blks > 1) {
						cfp[cfptop - 1].blks--; /* reduce previous by 1 */
						cfptop++;
						cfp[cfptop-1].blks = 1;   /* the new is a repeat */
					}
					cfp[cfptop-1].repeated = 1;
				}
				/* same mask as before */
				cfp[cfptop - 1].blks ++;
			} else {
				/* new mask */
				prevmask = mask;
				BINS==64 ? (bitmask[masktop] = mask) : (bitmask[masktop/fits] |= mask<<((masktop%fits)*BINS));
				masktop++;

				if (cfptop> 0 && cfp[cfptop - 1].repeated == 0 && cfp[cfptop-1].blks < ((1<<MAXOFFSET)-1)) {
					cfp[cfptop - 1].blks++;
				} else {
					cfp[cfptop].blks = 1;
					cfp[cfptop].repeated = 0;
					cfptop++;
				}
			}
			mask = 0;
		}
		switch(coltype){
		case TYPE_bte: val.bval = *(char*) (col + i * stride[coltype]);  GETBIT(bit,bval); break;
		case TYPE_sht: val.sval = *(short*) (col + i * stride[coltype]);  GETBIT(bit,sval); break;
		case TYPE_int: val.ival = *(int*) (col + i * stride[coltype]);  GETBIT(bit,ival); break;
		case TYPE_lng: val.lval = *(long*) (col + i * stride[coltype]);  GETBIT(bit,lval); break;
		case TYPE_oid: val.ulval = *(unsigned long*) (col + i * stride[coltype]);  GETBIT(bit,ulval); break;
		case TYPE_flt: val.fval = *(float*) (col + i * stride[coltype]);  GETBIT(bit,fval); break;
		case TYPE_dbl: val.dval = *(double*) (col + i * stride[coltype]); GETBIT(bit,dval); break;
		default: bit =0;
		}
		mask= setBit(mask,bit);
	}
	/* last mask */
	if (prevmask == mask && cfptop > 0 && cfp[cfptop-1].blks < (1<<MAXOFFSET)) {
		if (cfp[cfptop - 1].repeated == 0) {
			if (cfp[cfptop - 1].blks == 1) { /* only 1 on previous    */
				cfp[cfptop - 1].repeated = 1;
			} else {
				cfp[cfptop - 1].blks--;  /* reduce previous by 1 */
				cfp[cfptop].blks = 1;    /* the new is a repeat */
				cfp[cfptop].repeated = 1;
				BINS==64 ? (bitmask[masktop] = mask) : (bitmask[masktop/fits] |= mask<<((masktop%fits)*BINS));
				masktop++;
				cfptop++;
			}
		}
		/* same mask as before */
		cfp[cfptop - 1].blks ++;
	} else {
		BINS==64 ? (bitmask[masktop] = mask) : (bitmask[masktop/fits] |= mask<<((masktop%fits)*BINS));
		masktop++;
		if (cfptop> 0 && cfp[cfptop - 1].repeated == 0) {
			cfp[cfptop - 1].blks++;
		} else {
			cfp[cfptop].blks = 1;
			cfp[cfptop].repeated = 0;
			cfptop++;
		}
	}

	/* end creation, stop timer */
	timer = usec()-timer;

	/* stats gathering and printing */
	stats(timer);
}


/* it makes sense to test the fingerprints for further compression */
/* to be redone */
unsigned long *newmask;
static int
cmp(const void *p1, const void *p2)
{
	int i1 = *(int*)p1;
	int i2 = *(int*)p2;
	return newmask[i1] < newmask[i2]? -1: (newmask[i1]==newmask[i2]? 0: 1);
}

void sortandcount()
{
	int *base, i, j, uniq, limit;
	long old;

	limit= BINS * masktop;
	base= (int*) malloc(sizeof(int) * limit);
	newmask= (unsigned long*) malloc(sizeof(unsigned long) * limit);
	/* should done differently */
	for( i = 0; i < masktop; i++)
	for( j = 0; j < BINS; j++){
		base[i * BINS +j]= i * BINS +j;
		newmask[i * BINS +j] = getMask(i * BINS+j);
	}
	qsort(base, limit, sizeof(long), cmp);
	old= newmask[base[0]];
	uniq=1;
	for ( i=0; i<limit; i++) {
		if ( old != bitmask[base[i]]){
			old = bitmask[base[i]];
			uniq++;
		}
		//printMask(bitmask[base[i]]); putchar('\n');
	}
	printf("%s sorting reduces mask %d to %d uniq with compressed bits %d\n",column, BINS*masktop, uniq, bitcover(histogram));
	free(base);
}

/* scan */

void scan() {
	int i;

	for( i = 1; i < colcount; i++) {

	}
}

/* check if a column is sorted */
void sortedproperty()
{
	int i;
	
#define checksorted(T) \
	for( i = 1; i < colcount; i++) \
		if ( ((T*)col)[i] < ((T*)col)[i-1]){ \
			sorted = 0; \
			break; \
		} \
	if( sorted == 0) \
	for( sorted=1, i = 1; i< colcount; i++) \
		if ( ((T*)col)[i] > ((T*)col)[i-1]){ \
			sorted = 0; \
			break; \
		}

		switch(coltype){
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
	printf("%s sorted property %d\n",column,sorted);
}

/* calculate the statistics of the bits in the fingerprints */
/* calculate the spread of bits over the mask. It determines the filter capability for ranges */
/* calculate the number of bits per mask to determine single bin occupancy */
void
statistics()
{
	double var = 0, bitvar = 0;
	double delta = 0, bitdelta = 0;
	double mean = 0, bitmean = 0;
	double c;
	long d;
	long edit, on;
	unsigned long mask;
	int bitcnt, i, j, first, last;

	assert(globalmask);
	assert(masktop);
	for( i= 0; i< masktop; i++){
		first= -1;
		bitcnt = 0;
		mask = getMask(i);
		for( j=0; j< BINS; j++)
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
		column, mean, sqrt(var/masktop), bitmean, sqrt(bitvar/masktop), BINS, rpp);

	/* edit distance */
	edit = 0; on = 0;
	for (i=0; i< masktop; i++) {
		mask = getMask(i);
		for (j=0; j<BINS; j++)
			if (isSet(mask,j)) on++;
		if (i > 0) {
			mask = mask ^ getMask(i-1);
			for (j=0; j<BINS; j++)
				if (isSet(mask,j)) edit++;
		}
	}
	printf("%s total bits on %ld total edit distance %ld entropy is %lf\n", 
		column, on, edit, (double)edit/(double)(2*on));
}

long WAHencode()
{
	int *encoding;
	int fills=0, literals=0;
	int idx, i,j,n =0, t =0, b2;
	int bin,bit1,bit2;
	ValRecord val;
	long total = 0;


	wahtimer = usec();
	/* stage one, create all binned bitmaps */
	for ( idx = 0; idx < BINS; idx++){
		/* we need 2 bits extra per word at max , take care of small colcount */
		wah[idx].bits = (int*) malloc(sizeof(int) * (colcount/30+2));
		memset(wah[idx].bits, 0, sizeof(int) * (colcount/30+2));
		wah[idx].n = 0;
	}

	switch(coltype){
	case TYPE_bte: 
		for( i = 0; i< colcount; i++ ){
			wahbit(char,i,bin,bval);
			wah[bin].bits[i/32] = wahsetBit(wah[bin].bits[i/32],i & 31);
		}
		break;
	case TYPE_sht: 
		for( i = 0; i< colcount; i++ ){
			wahbit(short,i,bin,sval);
			wah[bin].bits[i/32] = wahsetBit(wah[bin].bits[i/32],i & 31);
		}
		break;
	case TYPE_int: 
		for( i = 0; i< colcount; i++ ){
			wahbit(int,i,bin,ival);
			wah[bin].bits[i/32] = wahsetBit(wah[bin].bits[i/32],i & 31);
		}
		break;
	case TYPE_lng: 
		for( i = 0; i< colcount; i++ ){
			wahbit(long,i,bin,lval);
			wah[bin].bits[i/32] = wahsetBit(wah[bin].bits[i/32],i & 31);
		}
		break;
	case TYPE_oid: 
		for( i = 0; i< colcount; i++ ){
			wahbit(long,i,bin,ulval);
			wah[bin].bits[i/32] = wahsetBit(wah[bin].bits[i/32],i & 31);
		}
		break;
	case TYPE_flt: 
		for( i = 0; i< colcount; i++ ){
			wahbit(float,i,bin,fval);
			wah[bin].bits[i/32] = wahsetBit(wah[bin].bits[i/32],i & 31);
		}
		break;
	case TYPE_dbl: 
		for( i = 0; i< colcount; i++ ){
			wahbit(double,i,bin,dval);
			wah[bin].bits[i/32] = wahsetBit(wah[bin].bits[i/32],i & 31);
		}
		break;
	default:
		printf("Unknown WAH type %d\n",coltype);
	}

	/* stage two, compress all binned bitmaps */
	fills = literals = 0;
	for (idx = 0; idx < BINS; idx++) {
		encoding = (int*) malloc(sizeof(int) * (colcount/30+2));
		memset(encoding, 0, sizeof(int) * (colcount/30+2));
		n = 0; /* number of encoding words */

		for (i = 0; i< colcount; ) {
			bin = wahisSet(wah[idx].bits[i/32], i & 31);
			for (j = i+1; j< colcount; j++) {
				bit2 = wahisSet(wah[idx].bits[j/32], j & 31);
				if (bit2 != bin) break;
			}
			/* see if we can perform runlength compression */
			if (j-i > 30 && j-i < (1<<30)) {
				/* construct fill word */
				fills++;
				if (bin)
					encoding[n] = wahsetBit(encoding[n],30);
				encoding[n] |= (j-i) ;
				i = j;
			} else {
				/* contruct literal word */
				int m;
				literals++;
				encoding[n] = wahsetBit(encoding[n],31);
				for (m=0, j = i; m < 31 && j<colcount; j++,m++) {
					bit2 = wahisSet(wah[idx].bits[j/32], j & 31);
					if (bit2)
						encoding[n] = wahsetBit(encoding[n],m);
				}
				i += m;
			}
			n++;
		}
		total += n;
		wah[idx].n = n;
		if (wah[idx].bits) free(wah[idx].bits);
		wah[idx].bits = encoding;
	}
	wahtimer = usec() - wahtimer;
#if 0
	printf("fills %d literals %d\n",fills, literals);
	for (idx = 0; idx < BINS; idx++) {
		int sum = 0;
		int ones = 0;
		printf (" n = %d ", wah[idx].n);
		for (i=0; i<wah[idx].n; i++) {
			if (wahisSet(wah[idx].bits[i],31)) {
				for (j=0;j<31;j++) if (wahisSet(wah[idx].bits[i],j)) ones++;
				sum+=31;
				if (sum>colcount) {
					printf ("sum is over %d but ", sum);
					sum=colcount;
				}
			} else {
				sum += (~(1<<30)) & wah[idx].bits[i];
				if (wahisSet(wah[idx].bits[i],30)) ones += (~(1<<30)) & wah[idx].bits[i];
			}
		}
		printf ("sum = %d ones = %d\n", sum, ones);
	}
	printf ("\n");
#endif
	return total;
}

int main(int argc, char **argv)
{
	long oid;
	FILE *cfile;
	long ns;
	int i;
	long n,count;
	char *c,*r;

	if ( argc < 5 || argc > 6) {
		printf("usage: %s type count file column [histotype] \n", argv[0]);
		return -1;
	}
	if (argc == 6 && strcmp(argv[5],"eqheight")) {
		HISTO_TYPE = EQHEIGHT;
	} else {
		HISTO_TYPE = EQWIDTH;
	}
	strcpy(column,argv[4]);
	count  = atoi(argv[2]);
	if ( strcmp(argv[1],"tinyint")== 0 || strcmp(argv[1],"boolean") == 0) {
		coltype= TYPE_bte;
		min.bval = 127;
		max.bval = -127;
	}
	if (strcmp(argv[1],"char") == 0 ||  strcmp(argv[1],"smallint")== 0 ||  strcmp(argv[1],"short")== 0) {
		coltype= TYPE_sht;
		min.sval = 32767;
		max.sval = -32767;
	}
	if ( strcmp(argv[1],"decimal")== 0 || strcmp(argv[1],"int")== 0 || strcmp(argv[1],"date")==0){
		coltype= TYPE_int;
		min.ival = INT_MAX;
		max.ival = INT_MIN; 
	}
	if ( strcmp(argv[1],"long")== 0 || strcmp(argv[1],"bigint")== 0){
		coltype= TYPE_lng;
		min.lval = LONG_MAX;
		max.lval = LONG_MIN; 
	}
	if ( strcmp(argv[1],"float")== 0 || strcmp(argv[1],"real") == 0) {
		coltype= TYPE_flt;
		min.fval = FLT_MAX;
		max.fval = FLT_MIN; 
	}
	if ( strcmp(argv[1],"double")== 0 ) {
		coltype= TYPE_dbl;
		min.dval = DBL_MAX;
		max.dval = -DBL_MAX;
	}
	if ( strcmp(argv[1],"oid")== 0) {
		coltype= TYPE_oid;
		min.lval = ULONG_MAX;
		max.lval = 0;
	}

	absmin = max;
	absmax = min;
	if ( strcmp(argv[1],"clob")== 0 || strcmp(argv[1],"string")== 0 || strcmp(argv[1],"varchar")== 0) {
		coltype= TYPE_str;
	}

	strcpy(filename,argv[3]);
	cfile = fopen(argv[3],"r");
	if ( cfile == NULL){
		printf("failed to open column file %s\n", argv[1]);
		return -1;
	}
	fseek(cfile, 0, SEEK_END);
	filesize = ftell(cfile);
	if ( filesize == 0){
		printf("empty open column file %s\n", argv[1]);
		return -1;
	}
	r= strrchr(column,'.');
	if ( r){
		*r = ' ' ;
		typename= r+1;
	}
	if (count && coltype == TYPE_str ) {
		/* the dictionary indices are assumed to be ordered on the underlying strings */
		/* it can only be used for point queries for now */
		if ( filesize/count <= 1) {
			coltype = TYPE_bte;
			min.bval = 127;
			max.bval = -127;
		} else
		if ( filesize/count <= 2){
			coltype= TYPE_sht;
				min.sval = 32767;
				max.sval = - 32767;
		} else
		if ( filesize/count <= 4){
			coltype= TYPE_int;
			min.ival = INT_MAX;
			max.ival = INT_MIN; 
		} else {
			coltype = TYPE_lng;
			min.lval = LONG_MAX;
			max.lval = LONG_MIN; 
		}
		absmin = max;
		absmax = min;
		printf("varchar stride %d size %ld count %ld\n",stride[coltype], filesize,count);
	}
	col= (char*) malloc(sizeof(col) * filesize);
	if( col == 0) {
		printf("malloc failed %ld\n",filesize * sizeof(col));
		return -1;
	}
	rewind(cfile);
	if ( (ns =fread(col, 1, filesize, cfile)) != filesize){
		printf("Could read %ld of %ld bytes\n", ns, filesize);
		return 0;
	}
	fclose(cfile);


	devnull = fopen("/dev/null","a");
	if ( devnull == NULL){
		printf("can not open /dev/null\n");
		return 0;
	}
	rpp  = stride[coltype]? PAGESIZE / stride[coltype]: 1;
	colcount = count;
	pages = colcount/rpp + 1;
	if (pages > MAXCFP) {
		printf("WARNING: there are too many pages\n");
	}

	if ( rpp ==0)
		return 0;
	
	printf("%s fingerprint %s %s "
			"size %ld "
			"type %d "
			"stride %d "
			"records %ld "
			"pagesize %ld "
			"rpp %d "
			"pages %ld \n",
			column, filename, argv[2],
			filesize,
			coltype,
			stride[coltype],
			colcount,
			sysconf(_SC_PAGESIZE),
			rpp,
			pages);

	sortedproperty();
	scan();
	zoneMap();		/* also calculates the domain range in absmin/absmax */
	cfpMask();

	/* sortandcount(); to be adjusted to smaller masks */
	printf("%s storage comparison tuples %ld size %ld zonemap %ld %ld%% cfp %ld %ld%%",
		column, colcount, filesize,
		zonetop*2*stride[coltype], ((long)zonetop * 2 * stride[coltype]*100)/filesize,
		((long)(masktop/(BITS/BINS))* sizeof(long)+ cfptop * sizeof(Fingerprint)),
		100 * ((long)(masktop/(BITS/BINS))* sizeof(long)+ cfptop * sizeof(Fingerprint)) /filesize);
#ifdef WAHTEST
		/* determine the equivalent WAH storage cost */
		n = WAHencode();
		printf(" wah %ld %ld%%", sizeof(int) * n, 100 * n * sizeof(int) /filesize);
		printf(" total fingerprints %d total headers %d\n",masktop,cfptop);
		printf("%s wah creation time %ld %ld usec per thousand\n",column, wahtimer,  ((long)wahtimer * 1000)/colcount);
#else
	printf(" total fingerprints %d total headers %d\n",masktop,cfptop);
#endif
	/* run queries */
	queries();
	statistics();


	printf("end of run\n");
	/* before exiting free memory for cleaness */
	free(cfp);
	free(bitmask);
	free(zmap);
}
