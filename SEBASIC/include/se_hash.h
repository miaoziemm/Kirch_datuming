#ifndef SE_HASH_H
#define SE_HASH_H

#include <cstdio>
#include <cstdlib>
/* Avoid C complex in C++ TU to prevent name/macro clashes */
#ifndef __cplusplus
#include <complex.h>
#endif
#include <string.h>
#include <errno.h>
#include <math.h>
#include "se_type.h"
#include "se_alloc.h"
#include "se_log.h"
#include "se_assert.h"
/** A hash table */
typedef void* se_hash;

/**
 * Constructs an empty hash table with the default initial capacity
 * (16) and the default load factor (0.75).
 */
se_hash create_hash();

/**
 * Constructs an empty hash table with the specified initial
 * capacity and load factor.
 */
se_hash create_hash_args(int initial_capacity, float load_factor);

/**
 * Destroys the hash map;
 */
void destroy_hash( se_hash ht );
void destroy_hash_and_entries( se_hash ht, int free_entries );
void destroy_hash_and_keyval( se_hash _ht, int free_keys, int free_values );


/**
 * Returns the value to which the specified key is mapped in the h identity
 * hash map, or NULL if the map contains no mapping for this key.
 * A return value of NULL does not necessarily indicate that the map
 * contains no mapping for the key; it is also possible that the map
 * explicitly maps the key to NULL. The ht_contains_key function may be
 * used to distinguish these two cases.
 */
void* ht_get( se_hash ht, const char* key );
void* ht_get_hash( se_hash _ht, const char* key, size_t keylen, int h );

/**
 * Returns true (1) if this map contains a mapping for the
 * specified key.
 */
int ht_contains_key( se_hash ht, const char* key );
int ht_contains_key_hash( se_hash _ht, const char* key, size_t keylen, int hash );

/**
 * Associates the specified value with the specified key in this map.
 * If the map previously contained a mapping for this key, the old
 * value is replaced.
 *
 * Returns previous value associated with specified key, or NULL
 * if there was no mapping for key. A NULL return can also indicate
 * that the hash map previously associated NULL with the specified key.
 */
char* ht_put( se_hash ht, const char* key, void* value );
char* ht_put_hash( se_hash _ht, const char* key, size_t keylen,
                            int hash, void* value );

/**
 * Removes the mapping for this key from this map if present.
 * Returns previous value associated with specified key, or NULL
 * if there was no mapping for key. A NULL return can also indicate
 * that the map previously associated NULL with the specified key.
 */
void* ht_remove( se_hash ht, const char* key );
void* ht_remove_hash( se_hash _ht, const char* key, size_t keylen, int hash);

/**
 * Returns the number of entries from the hash table
 **/
int ht_size( se_hash ht );

double ht_get_dispersion( se_hash ht );



/** Hash table (read-only) iterator */
typedef void* htiter;

/**
 * Creates a hash table read-only iterator and points it to the first entry
 * from the hash table (if any)
 **/
htiter create_htiter( se_hash ht );

/**
 * Destroys a hash table iterator
 **/
void destroy_htiter( htiter hti );

/**
 * Returns 1 if the hash table iterator has more entries, or 0 otherwise
 **/
int hti_hasnext( htiter hti );

/**
 * Output the current entry from the hash table and move the iterator forward
 **/
void hti_next( htiter hti,
                   const char** key,    /* hash table entry's key, can be NULL */
                   void** value );      /* hash table entry's value, can be NULL */

#endif