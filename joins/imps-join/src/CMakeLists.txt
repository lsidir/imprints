list(APPEND impsjoin-core_sources
	generator.c
	imprints.c
	utils.c
	join.c
    )

add_library(impsjoin-core STATIC ${impsjoin-core_sources})
