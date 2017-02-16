#include "simd_imprints.h"

void
compareImprintsIndex(Column column, Imprints_index *imps1, Imprints_index *imps2)
{
	printf("%s\t", column.colname);
	printf("%s\tindex\tIMPS1\tIMPS2\n", column.colname);
	printf("%s\tdct_cnt\t%lu\t%lu\n", column.colname, imps1->dct_cnt, imps2->dct_cnt);
	printf("%s\timps_cnt\t%lu\t%lu\n", column.colname, imps1->imps_cnt, imps2->imps_cnt);
	printf("%s\tbins\t%d\t%d\n", column.colname, imps1->bins, imps2->bins);
	printf("%s\tblocksize\t%d\t%d\n", column.colname, imps1->blocksize, imps2->blocksize);
}
