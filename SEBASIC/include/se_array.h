#ifndef SE_ARRAY_H
#define SE_ARRAY_H

#include <string.h>
#include "se_byte.h"
#include "se_hash.h"
#include "se_assert.h"


/**
 * The default initial capacity of an array list.
 **/
#define ARRAY_INITIAL_CAPACITY   8

/**
 * The default capacity increment of an array list.
 **/
#define ARRAY_DEF_CAPACITY_INC   8

/**
 * The internal array list struct.
 **/
typedef struct {
    /**
     * The value array.
     */
    void** values;
    /**
     * The value count.
     */
    int size;
    /**
     * The capacity of the value array.
     */
    int capacity;
    /**
     * Array's capacity increment (default \c ARRAY_DEF_CAPACITY_INC).
     */
    int capacity_inc;
} _array;

/** A array is a pointer to the _array structure */
typedef _array* array;


/**
 * Creates an empty array list with the default initial capacity (8).
 **/
array create_array( void );

/**
 * Creates an empty array list with the specified initial capacity.
 **/
array create_array_args( int initial_capacity );

/**
 * Creates an empty array list with the specified initial capacity and
 * capacity increment.
 **/
array create_array_args2( int initial_capacity, int capacity_increment );

/**
 * Destroys the array list.
 * free() the values if 'free_values' is true (non-zero)
 **/
void destroy_array( array a, int free_values );

/**
 * Increases the capacity of the array list, if necessary, to ensure that it can
 * hold at least 'min_capacity' values.
 **/
void array_ensure_capacity( array a, int min_capacity );

/**
 * Returns the number of values from the array list.
 **/
inline int array_size( const array a )
{
    ASSERT(a);

    return a->size;
}

/**
 * Returns the capacity of the array list.
 **/
inline int array_capacity( const array a )
{
    ASSERT(a);

    return a->capacity;
}


/**
 * Returns the value at the specified index.
 **/
inline void* array_get_at( array a,
                                   int index ) {
    ASSERT(a);
    ASSERT(index >= 0 && index < a->size);

    return a->values[index];
}

/**
 * Returns the first value from the array list.
 **/
inline void* array_get_first( array a ) {
    ASSERT(a);
    ASSERT(a->size);

    return a->values[0];
}

/**
 * Returns the last value from the array list.
 **/
inline void* array_get_last( array a ) {
    ASSERT(a);
    ASSERT(a->size);

    return a->values[a->size - 1];
}

/**
 * Sets the value at the specified index.
 * Returns the previous value at the specified index.
 **/
inline void* array_set_at( array a,
                                   int index,
                                   void* value ) {
    void* old_value;

    ASSERT(a);
    ASSERT(index >= 0 && index < a->size);

    old_value = a->values[index];
    a->values[index] = value;

    return old_value;
}

/**
 * Sets all values from within the specified range (inclusive).
 * Does _not_ free() the previous values.
 **/
void array_set_range( array a,
                          int index_from,
                          int index_to,
                          void* value );


/**
 * Adds the specified value at the end of the array list. The capacity of the array
 * list is increased automatically with the default capacity increment (8) if necessary.
 **/
inline void array_add( array a,
                               void* value )
{
    ASSERT(a);

    array_ensure_capacity( a, a->size + 1 );

    a->values[a->size++] = value;
}

/**
 * Swaps the two values at the specified indexes.
 **/
inline void array_swap( array a,
                                int index1,
                                int index2 )
{
    void* temp;
    ASSERT(a);
    ASSERT(index1 >= 0 && index1 < a->size);
    ASSERT(index2 >= 0 && index2 < a->size);

    temp = a->values[index1];
    a->values[index1] = a->values[index2];
    a->values[index2] = temp;
}


/**
 * Adds 'n' times the specified value at the end of the array list. The capacity of
 * the array list is increased automatically if necessary.
 **/
void array_add_n( array a,
                      void* value,
                      int n );

/**
 * Adds all values from the 'a_src' array list at the end of the 'a_dest' array list.
 * The capacity of the 'a_dest' array list is increased automatically if necessary.
 **/
void array_add_all( array a_dest,
                        const array a_src );


/**
 * Inserts the specified value at the specified index in the array list. The capacity
 * of the array list is increased automatically with the default capacity increment (8)
 * if necessary.
 **/
void array_insert( array a,
                       int index,
                       void* value );

/**
 * Inserts 'n' times the specified value at the specified index in the array list. The
 * capacity of the array list is increased automatically with the default capacity
 * increment (8) if necessary.
 **/
void array_insert_n( array a,
                         int index,
                         void* value,
                         int n);

/**
 * Inserts all values from the 'a_src' array list at the specified index in the 'a_dest'
 * array list. The capacity of the 'a_dest' array list is increased automatically with
 * the default capacity increment (8) if * necessary.
 **/
void array_insert_all( array a_dest,
                           int index,
                           const array a_src );


/**
 * Removes from the array list the first value that matches the specified key.
 * Returns the removed value (that can be NULL), or NULL if no match was found.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 * Does _not_ call free() for the removed value.
 **/
void* array_remove( array a,
                        compare_values_fn compare,
                        const void* key );

/**
 * Removes from the array list the value at the specified index.
 * Returns the removed value.
 * Does _not_ call free() for the removed value.
 **/
void* array_remove_at( array a,
                           int index );

/**
 * Remove the specified range (inclusive) of values from the array list.
 * Does _not_ call free() for the removed value.
 **/
void array_remove_range( array a,
                             int index_from,
                             int index_to );

/**
 * Removes from the array all values that matches the specified key.
 * Returns the number of removed values.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 * free() the values if 'free_values' is true (non-zero)
 **/
int array_remove_matches( array a,
                              compare_values_fn compare,
                              const void* key,
                              int free_values );

/**
 * Reverses the order of the values from the array list.
 **/
void array_reverse( array a );

/**
 * Returns the index of the minimum value from the array list, according to the order
 * given by the function referenced by 'compare'.
 * The comparison function must return an integer less than, equal to, or greater
 * than zero if the first argument is considered to be respectively less than, equal
 * to, or greater than the second.
 **/
int array_min( const array a,
                   compare_values_fn compare );

/**
 * Returns the index of the maximum value from the array list, according to the order
 * given by the function referenced by 'compare'.
 * The comparison function must return an integer less than, equal to, or greater
 * than zero if the first argument is considered to be respectively less than, equal
 * to, or greater than the second.
 **/
int array_max( const array a,
                   compare_values_fn compare );

/**
 * Removes all values from the array list.
 * free() the values if 'free_values' is true (non-zero)
 **/
void array_clear( array a, int free_values );



/**
 * Shrinks the capacity of the array list to be equal with it's size.
 **/
void array_shrink( array a );


/**
 * Returns true (non-zero) if the array contains the specified key, or false
 * (zero) no match was found.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 **/
int array_contains( const array a,
                        compare_values_fn compare,
                        const void* key );

/**
 * Returns the index of the first value matching the specified key, or -1 if
 * no match was found.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 **/
int array_indexof( const array a,
                       compare_values_fn compare,
                       const void* key );

/**
 * Returns the index of the last value matching the specified key, or -1 if
 * no match was found.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 **/
int array_last_indexof( const array a,
                            compare_values_fn compare,
                            const void* key );


/**
 * Sorts the values from specified array list in ascending order according to a
 * comparison function referenced to by compare, which is called with two arguments
 * that point to the values being compared.
 * The comparison function must return an integer less than, equal to, or greater
 * than zero if the first argument is considered to be respectively less than, equal
 * to, or greater than the second.  If two members compare as equal, their order in
 * the sorted array is undefined.
 **/
void array_sort( array a,
                     compare_values_fn compare );

/**
 * Returns the index of the value matching the specified key, or -1 if no match
 * was found.
 * The values from the array list should be in ascending sorted order according
 * to the comparison function referenced by 'compare'.
 * The 'compare' function is expected to have two arguments which point to the
 * key object and to an array value, in this order, and should return an integer
 * less than, equal to, or greater than zero if the key is found, respectively,
 * to be less than, to match, or be greater than the array value.
 **/
int array_binary_search( const array a,
                             compare_values_fn compare,
                             const void* key );


/**
 * Process all the values from the array list using the function referenced by
 * 'process'.
 **/
void array_process( array a,
                        process_value_fn process );

/**
 * Process the specified range (inclusive) of values from the array list using
 * the function referenced by 'process'.
 **/
void array_process_range( array a,
                              process_value_fn process,
                              int index_from,
                              int index_to );

/**
 * Returns the number of array values that matches the specified key.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 **/
int array_count_matches( const array _a,
                             compare_values_fn compare,
                             const void* key );

#endif