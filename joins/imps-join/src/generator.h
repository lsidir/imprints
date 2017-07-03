/**
 * @file    generator.h
 * @author  Cagri Balkesen <cagri.balkesen@inf.ethz.ch>
 * @date    Fri May 18 14:05:07 2012
 * @version $Id: generator.h 3017 2012-12-07 10:56:20Z bcagri $
 * 
 * @brief  Provides methods to generate data sets of various types
 * 
 */

#ifndef GENERATOR_H
#define GENERATOR_H

#include "types.h"
#include "common.h"

/** 
 * @defgroup DataGeneration Data Set Generation
 * @{
 */

/**
 * Seed the random number generator before calling create_relation_xx. If not
 * called, then generator will be initialized with the time of the call which
 * produces different random numbers from run to run.
 */
void 
seed_generator(unsigned int seed);

/**
 * Create relation with non-unique keys uniformly distributed between [0, maxid]
 */
int 
create_column_nonunique(Column *col, const int32_t maxid);

/**
 * Create relation with only primary keys (i.e. keys are unique from 1 to
 * col->colcount)
 */
int 
create_relation_pk(Column* col);


/** 
 * Create a foreign-key column using the given primary-key relation.
 * If the keys in pkcol is randomly distributed in
 * the full integer range, then 
 */
int 
create_column_fk_from_pk(Column *fkcol, Column *pkcol);



/** @} */

#endif /* GENERATOR_H */
