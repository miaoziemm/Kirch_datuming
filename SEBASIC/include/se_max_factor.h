#ifndef SE_MAX_FACTOR_H
#define SE_MAX_FACTOR_H

#include <stdint.h>
#include <stddef.h>  // For size_t
#include "se_alloc.h"
#include "se_type.h"

/** If this program is ever modified to factor integers larger
    than 2^128, this constant (and the algorithm) will have to change.  */

#define SE_MAX_N_FACTORS 128

/** Returns a factorization of n0 into the factors array. It is
    assumed that factors can hold up to ZTR_MAX_N_FACTORS numbers */

int se_factor(uint64_t n0, uint64_t *factors);

/** Splits number in nfct aproximatively equal factors */
void se_split_factors(uint64_t number, uint64_t* fct, size_t nfct);

/** Returns the smallest integer bigger or equals to number that has a
    max prime factor smaller or equal than max_factor. If growth_limit
    is > 0, and no integer smaller than number+growth_limit is found,
    then the one with the smallest bigger prime factor between number
    and number+growth_limit is returned. */

uint64_t se_next_int_with_max_factor(uint64_t number, uint64_t max_factor, int64_t growth_limit);



#endif