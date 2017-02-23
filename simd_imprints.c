#include "simd_imprints.h"

/* auxilary variables used globaly */
unsigned long pages;                            /* total pages in the column */
int rpp;                                        /* rows per page */

/* timing variables */
long zone_create_time;
long imprints_create_time;
long simd_imprints_create_time;

/* functions (in call order)-ish */
Zonemap_index *create_zonemaps(Column *column);

void queries(Column *column, Zonemap_index *zonemaps, Imprints_index *scalar_imps, Imprints_index *simd_imps, Imprints_index **exper_imps);
void simd_queries(Column *column, Imprints_index *imps, ValRecord low, ValRecord high, long results);
void genQueryRange(Column *column, Imprints_index *imps, int selectivity, ValRecord *low, ValRecord *high);

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
	Imprints_index *simd_imps;
	Imprints_index **exper_imps;

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
		column->min.ulval = ULONG_MAX;
		column->max.ulval = 0;
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
	             "typesize %d "
	             "records %ld "
	             "pagesize %d "
	             "sysconf(pagesize) %ld "
	             "rpp %d "
	             "pages %ld\n",
	             column->colname, column->filename,
	             filesize,
	             column->typename, column->coltype,
	             column->typesize,
	             column->colcount,
	             PAGESIZE,
	             sysconf(_SC_PAGESIZE),
	             rpp,
	             pages);

	/* check if column is sorted and set sorted = 1 if it is */
	isSorted(column);
	/* create zonemaps */
	zonemaps = create_zonemaps(column);
	/* create scalar imprints */
	scalar_imps = create_imprints(column, 64, 64, 0);
	/* create equivelant simd imprints */
	simd_imps = create_imprints(column, 64, 64, 1);

	exper_imps = (Imprints_index **)malloc(sizeof(Imprints_index *) * 9);

	exper_imps[0] = create_imprints(column, 64, 64, 1);
	exper_imps[1] = create_imprints(column, 64, 128, 1);
	exper_imps[2] = create_imprints(column, 64, 256, 1);
	exper_imps[3] = create_imprints(column, 128, 64, 1);
	exper_imps[4] = create_imprints(column, 128, 128, 1);
	exper_imps[5] = create_imprints(column, 128, 256, 1);
	exper_imps[6] = create_imprints(column, 256, 64, 1);
	exper_imps[7] = create_imprints(column, 256, 128, 1);
	exper_imps[8] = create_imprints(column, 256, 256, 1);
	/* run queries */
	queries(column, zonemaps, scalar_imps, simd_imps, exper_imps);

	VERBOSE printf("end of run\n");
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
	zonemaps->zonesize = rpp; /* FIX THIS EVERYWHERE */
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

void queries(Column *column, Zonemap_index *zonemaps, Imprints_index *scalar_imps, Imprints_index *simd_imps, Imprints_index **exper_imps)
{
	unsigned long res_cnt;
	unsigned long tuples[REPETITION];
	ValRecord low, high;
	long dummy, basetimer[REPETITION], zonetimer[REPETITION], impstimer[REPETITION], simd_impstimer[REPETITION];
	int i;

	for (i = 0; i < REPETITION; i++) {
		tuples[i] = 0;
		basetimer[i] = impstimer[i] = simd_impstimer[i] = zonetimer[i] = 0;
	}

	for (i = 0; i < REPETITION; i++) {
		genQueryRange(column, scalar_imps, i, &low, &high);

		/* simple scan */
		tuples[i] = simple_scan(column, low, high, &(basetimer[i]));

		res_cnt = zonemaps_scan(column, zonemaps, low, high, &(zonetimer[i]));
		if (tuples[i] != res_cnt) {
			printf("%s expecting %lu results and got %lu results from zonemaps\n", column->colname, tuples[i], res_cnt);
		}

		res_cnt = imprints_scan(column, scalar_imps, low, high, &(impstimer[i]));
		if (tuples[i] != res_cnt) {
			printf("%s expecting %lu results and got %lu results from scalar imprints\n", column->colname, tuples[i], res_cnt);
		}

		res_cnt = imprints_scan(column, simd_imps, low, high, &dummy);
		if (tuples[i] != res_cnt) {
			printf("%s expecting %lu results and got %lu results from simd imprints run on scalar queries (for debuging)\n", column->colname, tuples[i], res_cnt);
		}

		res_cnt = imprints_simd_scan(column, simd_imps, low, high, &(simd_impstimer[i]));
		if (tuples[i] != res_cnt) {
			printf("%s expecting %lu results and got %lu results from simd imprints\n", column->colname, tuples[i], res_cnt);
		}

		for (int k = 0; k < 9; k++) {
			res_cnt = imprints_simd_scan(column, exper_imps[k], low, high, &dummy);
			printf("%s "
					   "query[%d]=%12ld "
					   "selectivity=%2.1f%% \t"
					   "bins = %d \t"
					   "blocksize = %d(bytes)\t"
					   "simd_imprints = %ld"
					   "(usec)\n",
					   column->colname,
					   i, tuples[i],
					   tuples[i]* 100.0/column->colcount,
					   exper_imps[k]->bins,
					   exper_imps[k]->blocksize,
					   dummy
					   );
		}
	}

	for (i = 0; i < REPETITION; i++) {
		VERBOSE printf("%s "
					   "query[%d]=%12ld "
					   "selectivity=%2.1f%% \t"
					   "scan = %ld \t"
					   "zone = %ld \t"
					   "imprints = %ld \t"
					   "simd_imprints = %ld "
					   "(usec)\n",
					   column->colname,
					   i, tuples[i],
					   tuples[i]* 100.0/column->colcount,
					   basetimer[i],
					   zonetimer[i],
					   impstimer[i],
					   simd_impstimer[i]);
	}
}

/* simulate a series of queries */
void genQueryRange(Column *column, Imprints_index *imps, int selectivity, ValRecord *low, ValRecord *high)
{
#define setqueryrange(X)																\
	(*low).X = imps->bounds[1].X;					\
	(*high).X = (*low).X + selectivity * column->max.X/ (100.0/REPETITION);;			\
	if ((*high).X > column->max.X) (*high).X  = column->max.X;							\
	if ((*low).X > (*high).X) (*low).X = (*high).X;										\
	
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
