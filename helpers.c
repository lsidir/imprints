#include "simd_imprints.h"

void printMask(char *mask, int bytes)
{
	int i, j;
	for (i = 0; i < bytes; i++)
		for (j = 0; j < 8; j++)
			printf("%c", isSet(mask[i], j) ? 'x' : '.');
}

void
compareImprintsIndex(Column *column, Imprints_index *imps1, Imprints_index *imps2)
{
	unsigned long i, j;

	printf("%s dct_cnt   %lu\t%lu\n", column->colname, imps1->dct_cnt, imps2->dct_cnt);
	printf("%s imps_cnt  %lu\t%lu\n", column->colname, imps1->imps_cnt, imps2->imps_cnt);
	printf("%s bins      %d\t%d\n", column->colname, imps1->bins, imps2->bins);
	printf("%s blocksize %d\t%d\n", column->colname, imps1->blocksize, imps2->blocksize);

	for (i = 0, j = 0; (i <  imps1->dct_cnt) && (j <  imps1->dct_cnt); i++,j++) {
		if ((imps1->dct[i].blks != imps2->dct[j].blks) || (imps1->dct[i].repeated != imps2->dct[i].repeated))
			printf("%lu: imps1.blk = %d.%d imps2.blk = %d.%d\n", i,
					imps1->dct[i].blks, imps1->dct[i].repeated,
					imps2->dct[j].blks,imps2->dct[j].repeated);
	}
}

void
printBounds(Column *column, Imprints_index *imps)
{
	int i;

	printf("%s bounds of %d bins\n", column->colname, imps->bins);
	for ( i = 0; i < imps->bins; i++)
		switch(column->coltype){
		case TYPE_bte:
			printf("%10d:%d", i, imps->bounds[i].bval);
			break;
		case TYPE_sht:
			printf("%10d:%d", i, imps->bounds[i].sval);
			break;
		case TYPE_int:
			printf("%10d:%d", i, imps->bounds[i].ival);
			break;
		case TYPE_lng:
			printf("%10d:%ld", i, imps->bounds[i].lval);
			break;
		case TYPE_oid:
			printf("%10d:%lu", i, imps->bounds[i].ulval);
			break;
		case TYPE_flt:
			printf("%10d:%f", i, imps->bounds[i].fval);
			break;
		case TYPE_dbl:
			printf("%10d:%g", i, imps->bounds[i].dval);
		}
	printf("\n");
}

/*****************************
 * Printing Functions        *
 * ***************************/
#if 0

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

void printBins(Column column, ImprintsIndex *imps)
{
	int j;
	printf("%s bins ", column->colname);
	for ( j = 0; j < imps->bins; j++)
		switch(column->coltype){
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

	printf("%s rpp=%d imprint cells %ld zone cells %ld \n", column->colname, rpp, dct_cnt, zonemaps->zmaps_cnt);
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
				switch(column->coltype){
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
			switch(column->coltype){
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

	printf("%s histogram summary ", column->colname);
	i = 0;
	for ( j=0; j< bins; j++)
	if ( histogram[j] == 0)  i++;
		printf(" %d empty cells, bits used %d",i,bins-i);
	printf("\n%s histogram", column->colname);
	for ( j=0; j< bins; j++)
	if( histogram[j])
		printf("[%d] %ld ", j, histogram[j]);
	i=0;
	for ( j=0; j< bins; j++)
		 i += isSet(globalmask,j);
	printf("\n%s vectors dct_cnt %ld imps_cnt %ld unique %d bits %d ", column->colname, dct_cnt, imps_cnt, unique, i);
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
		column->colname, mean, sqrt(var/imps_cnt), bitmean, sqrt(bitvar/imps_cnt), bins, rpp);

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
		column->colname, on, edit, (double)edit/(double)(2*on));
}
#endif
