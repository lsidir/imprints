/**
 *
 */

//#include <cpuid.h>              /* for check_avx() */

#include <getopt.h>             /* getopt */

#include "src/generator.h"
#include "src/types.h"
#include "src/common.h"
#include "src/imprints.h"
#include "src/join.h"

//#include "pcm/system_perf_monitor.h"

#define REPEAT 1

//typedef struct joinalgo_t joinalgo_t;
typedef struct impsalgo_t impsalgo_t;
typedef struct cmdparam_t cmdparam_t;

struct impsalgo_t {
	char name[128];
	int (*impsalgorithm)(Column *, Column *, int, int);
};

struct cmdparam_t {
	//joinalgo_t * joinalgo;
	unsigned int joinalgo;		//0-> default, 1->multipass, 2->visitlist
	impsalgo_t * impsalgo;
	unsigned long l_size;		// outer relation, for probe
	unsigned long r_size;		// inner relation, for build
	unsigned int nthreads;
	unsigned int l_seed;
	unsigned int r_seed;
	int typesize;				//size of type in bytes
	double skew;
	int nonunique_keys;			// non-unique keys allowed?
	int fullrange_keys;			// keys covers full int range?

	unsigned int imps_bits;		// # of bits in the hash value used to encode imprints info; 6 by default

	/* read data from file */
	char r_col_filename[256];	// build column file name
	char l_col_filename[256];	// probe column file name
};



/** All available hashjoin algorithms */
#if 0
struct joinalgo_t {
    char name[128];
    query_result_t (*joinalgorithm)(Column *, Column *, unsigned int);
};

static struct joinalgo_t joinalgos [] =
{
  {"default", hashjoin},
  {"impsjoin", imps_hashjoin},
  {{0}, 0}
};
#endif

/** All available build imprints algorithms */
static struct impsalgo_t impsalgos [] =
{
		{"R2L", buildImpsR2L},
		{"L2R", buildImpsL2R},
		{{0}, 0}
};

//const char * DATAPATH = "/home/zeroxwj/Desktop/AllCode/CWI_imps_join/workspace/scripts/tpch100g_query_data_extract/";
const char * DATAPATH = "/home/zeroxwj/Desktop/AllCode/CWI_imps_join/workspace/microbenchmark/extract_from_monetdb/";

/* command line handling functions */
void parse_args(int argc, char ** argv, cmdparam_t * cmd_params);
//void print_args(cmdparam_t& cmd_params);
void print_help(char * progname);

int
main(int argc, char *argv[])
{
    Column* colL;
    Column* colR;
    //int64_t    results;

	/* Command line parameters */
	cmdparam_t cmd_params;

	/* Default values if not specified on command line */
	//cmd_params.joinalgo 	= &joinalgos[0];
	cmd_params.joinalgo 	= 0;
	cmd_params.impsalgo		= &impsalgos[0];
	cmd_params.nthreads 	= 1;
	cmd_params.l_size 		= 128000000;
	cmd_params.r_size		= 128000000;
	cmd_params.l_seed		= 12345;
	cmd_params.r_seed		= 54321;
	cmd_params.skew			= 0.0;
	cmd_params.nonunique_keys = 0;
	cmd_params.fullrange_keys = 0;
	cmd_params.typesize		= 4;

	cmd_params.imps_bits 	= 3;

	//cmd_params.r_col_filename = "";
	//cmd_params.l_col_filename = "";

	parse_args(argc, argv, &cmd_params);

	joinconf_t joinconf;
	joinconf.joinalgo 		= cmd_params.joinalgo;
	joinconf.imps_bits 		= cmd_params.imps_bits;

    /** first initialize (build) Column R */
    fprintf(stdout,
            "[INFO ] Creating column R with size = %.3lf MiB, #tuples = %ld : \n",
            (double) cmd_params.typesize * cmd_params.r_size/(1024.0*1024.0),
            cmd_params.r_size);
    fflush(stdout);

    colR = (Column *) malloc_aligned(sizeof(Column));
    colR->colcount = cmd_params.r_size;
	strcpy(colR->colname, "colR");
	//strcpy(colR->filename, argv[3]);
	//strcpy(colR->typename, argv[1]);
    colR->coltype  = (cmd_params.typesize == 4) ? TYPE_int : TYPE_lng;	 // use integer by default in microbenchmark
    colR->typesize = cmd_params.typesize;
    //colR->min.ival = INT_MAX;
    //colR->max.ival = INT_MIN;
    colR->min.lval = LONG_MAX;
    colR->max.lval = LONG_MIN;
    colR->imps_idx = NULL;

    colR->col = NULL;
    colR->col = (char *) malloc_aligned(colR->typesize * colR->colcount);

    if (!colR->col) {
    	fprintf(stderr, "Malloc failed for colR\n");
    	return -1;
    }

    /* Column R value generation*/
#if 0
    seed_generator(cmd_params.r_seed);

    if (cmd_params.fullrange_keys) {
    	create_column_nonunique(colR, INT_MAX);
    } else if (cmd_params.nonunique_keys) {
    	create_column_nonunique(colR, colR->colcount);
    } else {
    	create_relation_pk(colR);	// shuffled unique primary key
    }
#endif

    /* Read Column R from file */
    //printf("%s\n", cmd_params.r_col_filename);
    assert(strlen(cmd_params.r_col_filename) > 0);
    ReadColumnFromFile(colR, cmd_params.r_col_filename);

    fprintf(stdout, "[INFO ] Done\n");

	/** next initialize (probe) Column L */
    fprintf(stdout,
            "[INFO ] Creating column L with size = %.3lf MiB, #tuples = %ld : \n",
            (double) cmd_params.typesize * cmd_params.l_size/(1024.0*1024.0),
            cmd_params.l_size);
    fflush(stdout);

    colL = (Column *) malloc_aligned(sizeof(Column));
    colL->colcount = cmd_params.l_size;
	strcpy(colL->colname, "colL");
	//strcpy(colL->filename, argv[3]);
	//strcpy(colL->typename, argv[1]);
    colL->coltype  = (cmd_params.typesize == 4) ? TYPE_int : TYPE_lng;	 // use integer by default in microbenchmark
    colL->typesize = cmd_params.typesize;
    //colL->min.ival = INT_MAX;
    //colL->max.ival = INT_MIN;
    colL->min.lval = LONG_MAX;
    colL->max.lval = LONG_MIN;
    colL->imps_idx = NULL;

    assert(colL->coltype == colR->coltype);

    colL->col = NULL;
    colL->col = (char *) malloc_aligned(colL->typesize * colL->colcount);

    if (!colL->col) {
    	fprintf(stderr, "Malloc failed for colL\n");
    	return -1;
    }

    /* Column L value generation*/
#if 0
    seed_generator(cmd_params.l_seed);

    if (cmd_params.skew > 0) {
    	//gen_zipf()
    } else {
    	create_column_fk_from_pk(colL, colR);
    }
#endif

    /* Read Column L from file */
    assert(strlen(cmd_params.l_col_filename) > 0);
    ReadColumnFromFile(colL, cmd_params.l_col_filename);

    fprintf(stdout, "[INFO ] Done\n");

    /** Building imprints */
    fprintf(stdout, "[INFO ] Building imprints for columns\n");

    colL->imps_idx = (Imprints_index *) malloc_aligned(sizeof(Imprints_index));
    colR->imps_idx = (Imprints_index *) malloc_aligned(sizeof(Imprints_index));

    int max_bins = 64;
    int blocksize = 64;

    int imps_result = cmd_params.impsalgo->impsalgorithm(colL, colR, max_bins, blocksize);

    if (imps_result < 0) {
    	fprintf(stderr, "create imprints failed\n");
    	return -1;
    }

    fprintf(stdout, "[INFO ] Done\n");
#if 1
	query_result_t result;
	for (uint32_t i = 0; i < REPEAT; ++i) {
		result = hashjoin(colL, colR, cmd_params.nthreads, joinconf);
	}
#endif

	/* destroy imprints */
	free(colR->imps_idx->bounds);
	free(colL->imps_idx->bounds);

	if (colR->imps_idx->dct != NULL) {
		free(colR->imps_idx->dct);
	}

	if (colR->imps_idx->imprints != NULL) {
		free(colR->imps_idx->imprints);
	}

	if (colL->imps_idx->dct != NULL) {
		free(colL->imps_idx->dct);
	}

	if (colL->imps_idx->imprints != NULL) {
		free(colL->imps_idx->imprints);
	}

	free(colR->imps_idx);
	free(colL->imps_idx);

	/* destroy columns */
	free(colR->col);
	free(colL->col);
	free(colR);
	free(colL);

	return 0;
}

/* command line handling functions */
void print_help(char * progname) {
	printf("Usage: %s [options]\n", progname);

	printf("\n");

	printf(
			"\
\tData generation options, with default values in []:      		           			\n\
\t\t--gen-data			Just generate data and save as files; do not sort [NO] 		\n\
\t\t-r --nrows=<R>      Number of rows in each column [64M]							\n\
\t\t-c --ncolumns=<C>   Number of columns to be sorted [2]							\n\
\t\t-w --bitwidth=<B>   the bit width of each column; separated by `B' [25B21]		\n\
\t\t-z --zipf=<Z>	    data skew: [0, 20], 0 indicate uniform distribution [0]		\n\
\t\t-C --cardinality=<> # of distinct values, as perc [0,100] of max value [100]	\n\
\t\t-g --ngroups=<G> 	Number of groups in each column [2]							\n\
\t\t-i --inputfile=<F>  read column values from files (--gen-data is not set)		\n\
	   																				\n\
\tGeneral options:																	\n\
\t\t--encode-oid		Encode OID with necessary bits; otherwise use 32|64 [NO]	\n\
\t\t-n --nthreads=<N>   Number of threads to use [1]	                    		\n\
\t\t-a --asc-desc=<A>   Indicate ascending(1)|descending(2) of each column [12]		\n\
\t\t-s --sortalgo=<..>  Implementation of underlying sort algorithms. [ms]			\n\
						   Acceptable sort algorithms are:							\n\
						   ms: Merge sort											\n\
						   rs: Radix sort											\n\
\t\t-p --comptype=<..>  Approaches to compose sorting of multiple columns. [ch]		\n\
	   					   Acceptable compose types are:							\n\
						   ch: chaining-based; baseline as in MonetDB				\n\
						   st: stitching-based; 									\n\
\t\t-P --oid-pack=<..>  Indicate how OID are processed during sorting. [nonpack]	\n\
						   Acceptable pack types are:								\n\
						   pack: 	pack the oid with (stitched) column values		\n\
						   nonpack:	OID stored as separated array	   				\n\
\t\t-b --baseline=<..>	Indicate the versions of baseline	[Disabled...]			\n\
    						Acceptable baseline type are:							\n\
    						aligned32: always stored as 32bit for each column		\n\
    						adaptive:  fit the lowest bank size						\n\
																					\n\
\t\t-x --stitchstyle=<>	Indicate # of bits to stitch left or right [Disabled...]	\n\
    						Acceptable examples are:								\n\
    						l3: stitch to left 3 bits								\n\
    						r4: stitch to right 4 bits								\n\
\t\t-o --ordered=<>		Indicate Group-by case (0) or Order-by case (1) [1]			\n\
																					\n\
\tMerge sort related options:														\n\
\t\t-f --partfanout=<F> Fanout of partition phase in merge sort. [1]				\n\
	   					   (also equivalent to fanin of multiway merge phase)		\n\
																					\n\
\tRadix sort related options:														\n\
                                                                               		\n\
\tBasic user options                                                         		\n\
\t\t-h --help           Show this message                                   		\n");
}

//void print_args(cmdparam_t& cmd_params) {
//
//}

void
parse_args(int argc, char ** argv, cmdparam_t * cmd_params)
{
    int c;

    while(1) {
        static struct option long_options[] =
            {
                /* These options set a flag. */
                {"help",       no_argument,    0, 'h'},
                {"imps-bits",    no_argument,    0, 'm'},
                /* These options don't set a flag.
                   We distinguish them by their indices. */
                {"algo",    required_argument, 0, 'a'},
                {"nthreads",required_argument, 0, 'n'},
                /*{"perfconf",required_argument, 0, 'p'},*/
                /* columns file name */
                {"l-filename",  required_argument, 0, 'p'},	/*l=>probe*/
                {"r-filename",  required_argument, 0, 'b'},	/*r=>build*/
				{"typesize", required_argument, 0, 't'},
                {"l-size",  required_argument, 0, 'l'},
                {"r-size",  required_argument, 0, 'r'},
                /*{"perfout", required_argument, 0, 'o'},*/
                {"r-seed",  required_argument, 0, 'x'},
                {"s-seed",  required_argument, 0, 'y'},
                {"skew",    required_argument, 0, 'z'},
            };
        /* getopt_long stores the option index here. */
        int option_index = 0;

        c = getopt_long (argc, argv, "h:m:a:n:p:b:t:l:r:x:y:z",
                         long_options, &option_index);

        /* Detect the end of the options. */
        if (c == -1)
            break;
        switch (c)
        {
          case 0:
              /* If this option set a flag, do nothing else now. */
              if (long_options[option_index].flag != 0)
                  break;
              printf ("option %s", long_options[option_index].name);
              if (optarg)
                  printf (" with arg %s", optarg);
              printf ("\n");
              break;
          case 'a':
        	  cmd_params->joinalgo = atoi(optarg);	//0-> default, 1->multipass, 2->visitlist
              break;
          case 'm':
        	  cmd_params->imps_bits = atoi(optarg);	// vary # of partitions
              break;
          case 'h':
          case '?':
              /* getopt_long already printed an error message. */
              print_help(argv[0]);
              exit(EXIT_SUCCESS);
              break;

          case 'n':
              cmd_params->nthreads = atoi(optarg);
              break;

          case 'l':	/* the probe column */
              cmd_params->l_size = atol(optarg);
              break;

          case 'r': /* the build column */
              cmd_params->r_size = atol(optarg);
              break;

          case 'x':
              cmd_params->l_seed = atoi(optarg);
              break;

          case 'y':
              cmd_params->r_seed = atoi(optarg);
              break;

          case 'z':
              cmd_params->skew = atof(optarg);
              break;

          case 'p':
        	  //strcpy(cmd_params->l_col_filename, DATAPATH);
        	  //printf("setting l_col_filename\n");
        	  //strcpy(cmd_params->l_col_filename, "/home/zeroxwj/Desktop/AllCode/CWI_imps_join/workspace/scripts/tpch100g_query_data_extract/Q19_l_partkey.dat");
        	  //strcpy(cmd_params->l_col_filename, "/home/zeroxwj/Desktop/AllCode/CWI_imps_join/workspace/scripts/tpch100g_query_data_extract/");
              //strncat(cmd_params->l_col_filename, optarg, 128);
              //printf("%s\n", cmd_params->l_col_filename);
        	  //printf("%s\n", optarg);
        	  sprintf(cmd_params->l_col_filename, "%s%s", DATAPATH, optarg);
        	  //printf("%s\n", cmd_params->l_col_filename);
              break;

          case 'b':
              //strcpy(cmd_params->r_col_filename, DATAPATH);
        	  //printf("setting r_col_filename\n");
        	  //strcpy(cmd_params->r_col_filename, "/home/zeroxwj/Desktop/AllCode/CWI_imps_join/workspace/scripts/tpch100g_query_data_extract/Q19_p_partkey.dat");
        	  //strcpy(cmd_params->r_col_filename, "/home/zeroxwj/Desktop/AllCode/CWI_imps_join/workspace/scripts/tpch100g_query_data_extract/");
              //strncat(cmd_params->r_col_filename, optarg, 128);
              //printf("%s\n", cmd_params->r_col_filename);
        	  //printf("%s\n", optarg);
        	  sprintf(cmd_params->r_col_filename, "%s%s", DATAPATH, optarg);
              break;

          case 't':
              cmd_params->typesize = atoi(optarg);
              printf("typesize: %d\n", cmd_params->typesize);
              break;

          default:
              break;
        }
    }
}
