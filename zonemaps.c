/* This Source Code Form is subject to the terms of the Mozilla Public
 * License, v. 2.0. If a copy of the MPL was not distributed with this
 * file, You can obtain one at http://mozilla.org/MPL/2.0/. */

#include "simd_imprints.h"

Zonemap_index *
create_zonemaps(Column *column, int blocksize)
{
	Zonemap_index *zonemaps;
	ValRecord val;
	int values_per_zone;
	unsigned long total_zones;
	long timer;
	long i;

	/* malloc zonemap array */
	zonemaps = (Zonemap_index *) malloc(sizeof(Zonemap_index));
	zonemaps->zonesize = blocksize;
	values_per_zone = zonemaps->zonesize/column->typesize;
	total_zones = column->colcount/values_per_zone + 1;
	zonemaps->zmaps = (Zonemap *) malloc (sizeof(Zonemap)*(total_zones)); /* maybe +1 if core dump hits the fan */
	memset((char *)zonemaps->zmaps, 0, sizeof(Zonemap)*(total_zones));    /* maybe +1 if core dump hits the fan */

#define upd(X) \
		if (val.X < column->min.X) column->min.X = val.X; \
		if (val.X > column->max.X) column->max.X = val.X; \
		if (zonemaps->zmaps[zonemaps->zmaps_cnt].min.X > val.X) \
			zonemaps->zmaps[zonemaps->zmaps_cnt].min.X = val.X; \
		if (zonemaps->zmaps[zonemaps->zmaps_cnt].max.X < val.X) \
			zonemaps->zmaps[zonemaps->zmaps_cnt].max.X = val.X;

	timer = usec();
	zonemaps->zmaps_cnt = -1;
	for (i=0; i < column->colcount; i++) {
		if (!(i&(values_per_zone-1))) {
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
	if ((i-1)%values_per_zone) {
		zonemaps->zmaps_cnt++;
	}
	zonemaps->zmaps_cnt++;
	timer = usec()-timer;

	VERBOSE printf("%s zonemap  creation time=%ld, %ld usec per thousand values\n",
	               column->colname, timer, ((long)timer*1000)/column->colcount);
	return zonemaps;
}

