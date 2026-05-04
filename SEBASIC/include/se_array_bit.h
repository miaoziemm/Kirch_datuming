#ifndef SE_ARRAY_BIT_H
#define SE_ARRAY_BIT_H

#include "se_type.h"
#include "se_byte.h"
#include "se_assert.h"
#include "se_log.h"
#include "se_alloc.h"


typedef void* bit_array;

/**
 * Creates a bit set with initial size large enough to explicitly
 * represent bits with indices in the range 0 through nbits-1. All
 * bits are initially false.
 *
 * If nbits == 0, a default initial size will be used.
 */
bit_array ba_create(size_t nbits);

/** Release used resources. */
void ba_destroy(bit_array _ba);


/**
 * Sets the bit at the specified index to true.
 */
void ba_set(bit_array _ba, size_t bit_index);

/**
 * Sets the bit specified by the index to false.
 */
void ba_clear(bit_array _ba, size_t bit_index);

/**
 * Sets all the bits to false.
 */
void ba_clear_all(bit_array _ba);

/**
 * Returns the value of the bit with the specified index. The value is
 * true if the bit with the index bit_index is currently set;
 * otherwise, the result is false.
 */
int ba_get(const bit_array _ba, size_t bit_index);

/**
 * Performs a logical OR of arguments. The result is stored in the
 * first parameter.
 */
void ba_or(bit_array _ba1, const bit_array _ba2);

/**
 *
 */
void ba_status(bit_array _ba);


byte* ba_serialize(bit_array _ba, int destroy, int64_t* nbytes);

bit_array ba_deserialize(byte* buff, int allocmem);


#endif