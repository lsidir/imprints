/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "main.h"

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
	int stride[14] = {0,0,0,1,2,0,4,8,0,0,4,8,8,0}; /* sizeof(column type) according to MonetDB */
	unsigned long pages;                            /* total pages in the column */
	int vpp;                                        /* values per page */
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

	vpp = PAGESIZE/stride[column->coltype];
	if (vpp == 0) {
		printf("rows per pages is 0\n");
		return -1;
	}
	pages = column->colcount/vpp + 1;
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
	             "values per page %d "
	             "pages %ld\n",
	             column->colname, column->filename,
	             filesize,
	             column->typename, column->coltype,
	             column->typesize,
	             column->colcount,
	             PAGESIZE,
	             sysconf(_SC_PAGESIZE),
	             vpp,
	             pages);

	/* check if column is sorted and set sorted = 1 if it is */
	isSorted(column);
	/* create zonemaps */
	zonemaps = create_zonemaps(column, 64);
	/* create scalar imprints */
	scalar_imps = create_imprints(column, 64, 64, 0);
	/* create equivelant simd imprints */
	simd_imps = create_imprints(column, 64, 64, 1);

	/* create more versions of simd imprints with different bitvector and block sizes */
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

	/* run the queries */
	queries(column, zonemaps, scalar_imps, simd_imps, exper_imps);

	VERBOSE printf("end of run\n");
	return 1;
}
