#include "../include/se_array.h"


/* Internal functions */
static void _array_qsort(void** values, int left, int right, compare_values_fn compare);
static void _array_autoshrink(_array* a);

/**
 * Creates an empty array list with the default initial capacity (8).
 **/
array create_array(void)
{
    return create_array_args2( ARRAY_INITIAL_CAPACITY, ARRAY_DEF_CAPACITY_INC );
}

/**
 * Creates an empty array list with the specified initial capacity.
 **/
array create_array_args(int initial_capacity)
{
    return create_array_args2( initial_capacity, ARRAY_DEF_CAPACITY_INC );
}

/**
 * Creates an empty array list with the specified initial capacity and
 * capacity increment.
 **/
array create_array_args2( int initial_capacity, int capacity_inc )
{
    _array* a;
    ASSERT( initial_capacity >= 0 );
    ASSERT( capacity_inc > 0 );

    a = (_array*) malloc( sizeof(_array) );

    a->values = (void**)malloc( initial_capacity * sizeof(void*) );
    a->size = 0;
    a->capacity = initial_capacity;
    a->capacity_inc = capacity_inc;
    /*DEBUG(("capacity+: %d", a->capacity));*/

    return (array) a;
}

/**
 * Destroys the array list.
 * free() the values if 'free_values' is true (non-zero)
 **/
void destroy_array( array a, int free_values )
{
    ASSERT(a);

    array_clear(a, free_values);
    free(a);
}


/**
 * Sets all values from within the specified range (inclusive).
 * Does _not_ free() the previous values.
 **/
void array_set_range( array _a,
                          int index_from,
                          int index_to,
                          void* value )
{
    _array* a = (_array*) _a;
    ASSERT(a);
    ASSERT(index_from >= 0 && index_from < a->size);
    ASSERT(index_to >= 0 && index_to < a->size);
    ASSERT(index_from <= index_to);

    while( index_from <= index_to ) {
        a->values[index_from++] = value;
    }
}


/**
 * Adds 'n' times the specified value at the end of the array list. The capacity of
 * the array list is increased automatically if necessary.
 **/
void array_add_n( array _a,
                      void* value,
                      int n )
{
    _array* a = (_array*) _a;
    int i;
    ASSERT(a);
    ASSERT(n >= 0);

    array_ensure_capacity( a, a->size + n );

    for (i = a->size; i < a->size + n; ++i) {
        a->values[i] = value;
    }
    a->size += n;
}

/**
 * Adds all values from the 'a_src' array list at the end of the 'a_dest' array list.
 * The capacity of the 'a_dest' array list is increased automatically if necessary.
 **/
void array_add_all( array _a_dest,
                        const array _a_src )
{
    _array* a_dest = (_array*) _a_dest;
    const _array* a_src = (const _array*) _a_src;
    ASSERT(a_dest);
    ASSERT(a_src);

    array_ensure_capacity( a_dest, a_dest->size + a_src->size );

    memcpy( a_dest->values + a_dest->size,
            a_src->values,
            a_src->size * sizeof(void*) );
    a_dest->size += a_src->size;
}


/**
 * Inserts the specified value at the specified index in the array list. The capacity
 * of the array list is increased automatically with the default capacity increment (8)
 * if necessary.
 **/
void array_insert( array _a,
                       int index,
                       void* value )
{
    _array* a = (_array*) _a;
    ASSERT(a);
    ASSERT(index >= 0 && index <= a->size);

    array_ensure_capacity( a, a->size + 1 );

    memmove( a->values + index + 1,
             a->values + index,
             (a->size - index) * sizeof(void*) );
    a->values[index] = value;
    a->size++;
}

/**
 * Inserts 'n' times the specified value at the specified index in the array list. The
 * capacity of the array list is increased automatically with the default capacity
 * increment (8) if necessary.
 **/
void array_insert_n( array _a,
                         int index,
                         void* value,
                         int n )
{
    _array* a = (_array*) _a;
    int i;
    ASSERT(a);
    ASSERT(index >= 0 && index <= a->size);
    ASSERT(n >= 0);

    array_ensure_capacity( a, a->size + n );

    memmove( a->values + index + n,
             a->values + index,
             (a->size - index) * sizeof(void*) );
    for (i = index; i < index + n; ++i) {
        a->values[i] = value;
    }
    a->size += n;
}

/**
 * Inserts all values from the 'a_src' array list at the specified index in the 'a_dest'
 * array list. The capacity of the 'a_dest' array list is increased automatically with
 * the default capacity increment (8) if * necessary.
 **/
void array_insert_all( array _a_dest,
                           int index,
                           const array _a_src )
{
    _array* a_dest = (_array*) _a_dest;
    const _array* a_src = (const _array*) _a_src;
    ASSERT(a_dest);
    ASSERT(a_src);
    ASSERT(index >= 0 && index <= a_dest->size);

    array_ensure_capacity( a_dest, a_dest->size + a_src->size );

    memmove( a_dest->values + index + a_src->size,
             a_dest->values + index,
             (a_dest->size - index) * sizeof(void*) );
    memmove( a_dest->values + index,
             a_src->values,
             a_src->size * sizeof(void*) );
    a_dest->size += a_src->size;
}


/**
 * Returns true (non-zero) if the array contains the specified key, or false
 * (zero) no match was found.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 **/
int array_contains( const array a,
                        compare_values_fn compare,
                        const void* key )
{
    ASSERT(a);
    ASSERT(compare);

    return array_indexof( a, compare, key ) >= 0;
}

/**
 * Removes from the array list the first value that matches the specified key.
 * Returns the removed value (that can be NULL), or NULL if no match was found.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 * Does _not_ call free() for the removed value.
 **/
void* array_remove( array _a,
                        compare_values_fn compare,
                        const void* key )
{
    _array* a = (_array*) _a;
    int index;
    void* old_value;
    ASSERT(a);
    ASSERT(compare);

    index = array_indexof( a, compare, key );

    if( index >= 0 ) {
        old_value = a->values[index];
        array_remove_at( a, index );
        return old_value;
    } else {
        return NULL;
    }
}

/**
 * Removes from the array list the value at the specified index.
 * Returns the removed value.
 * Does _not_ call free() for the removed value.
 **/
void* array_remove_at( array _a,
                           int index )
{
    _array* a = (_array*) _a;
    void* value;
    ASSERT(a);
    ASSERT(index >= 0 && index < a->size);

    value = a->values[index];
    memmove( a->values + index,
             a->values + index + 1,
             (a->size - index - 1) * sizeof(void*) );
    a->size--;

    _array_autoshrink(a);

    return value;
}

/**
 * Remove the specified range (inclusive) of values from the array list.
 * Does _not_ call free() for the removed value.
 **/
void array_remove_range( array _a,
                             int index_from,
                             int index_to )
{
    _array* a = (_array*) _a;
    ASSERT(a);
    ASSERT(index_from >= 0 && index_from < a->size);
    ASSERT(index_to >= 0 && index_to < a->size);
    ASSERT(index_from <= index_to);

    memmove( a->values + index_from,
             a->values + index_to + 1,
             (a->size - index_to - 1) * sizeof(void*) );
    a->size -= index_to - index_from + 1;

    _array_autoshrink(a);
}

/**
 * Removes from the array all values that matches the specified key.
 * Returns the number of removed values.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 * free() the values if 'free_values' is true (non-zero)
 **/
int array_remove_matches( array _a,
                              compare_values_fn compare,
                              const void* key,
                              int free_values )
{
    _array* a = (_array*) _a;
    int i, count;
    void* value;
    ASSERT(a);
    ASSERT(compare);

    for (i = count = 0; i < a->size; ++i) {
        if(! compare( key, a->values[i] )) {
            value = array_remove_at( _a, i-- );
            if( free_values && value ) {
                free(value);
            }
            ++count;
        }
    }
    return count;
}


/**
 * Reverses the order of the values from the array list.
 **/
void array_reverse( array _a )
{
    _array* a = (_array*) _a;
    void **left, **right, *temp;
    ASSERT(a);

    left  = a->values;
    right = a->values + a->size - 1;
    while( left < right ) {
        temp = *left;
        *left = *right;
        *right = temp;
        ++left, --right;
    }
}

/**
 * Returns the index of the minimum value from the array list, according to the order
 * given by the function referenced by 'compare'.
 * The comparison function must return an integer less than, equal to, or greater
 * than zero if the first argument is considered to be respectively less than, equal
 * to, or greater than the second.
 **/
int array_min( const array _a,
                   compare_values_fn compare )
{
    const _array* a = (const _array*) _a;
    int i;
    int min;
    ASSERT(a);
    ASSERT(compare);

    if(! a->size) {
        return -1;
    }

    min = 0;
    for (i = 1; i < a->size; ++i) {
        if( compare( a->values[min], a->values[i] ) > 0 ) {
            min = i;
        }
    }
    return min;
}

/**
 * Returns the index of the maximum value from the array list, according to the order
 * given by the function referenced by 'compare'.
 * The comparison function must return an integer less than, equal to, or greater
 * than zero if the first argument is considered to be respectively less than, equal
 * to, or greater than the second.
 **/
int array_max( array const _a,
                   compare_values_fn compare )
{
    const _array* a = (const _array*) _a;
    int i;
    int max;
    ASSERT(a);
    ASSERT(compare);

    if(! a->size) {
        return -1;
    }

    max = 0;
    for (i = 1; i < a->size; ++i) {
        if( compare( a->values[max], a->values[i] ) < 0 ) {
            max = i;
        }
    }
    return max;
}


/**
 * Removes all values from the array list.
 * free() the values if 'free_values' is true (non-zero)
 **/
void array_clear( array _a, int free_values )
{
    _array* a = (_array*) _a;
    int i;
    ASSERT(a);

    if (free_values) {
        for (i = 0; i < a->size; ++i) {
            if( a->values[i] ) {
                free( a->values[i] );
            }
        }
    }

    if( a->values ) {
        free(a->values);
    }
    a->values = NULL;
    a->size = a->capacity = 0;
}


/**
 * Increases the capacity of the array list, if necessary, to ensure that it can
 * hold at least 'min_capacity' values.
 **/
void array_ensure_capacity( array _a, int min_capacity )
{
    _array* a = (_array*) _a;
    int new_capacity;
    void** new_values;
    ASSERT(a);
    ASSERT(min_capacity >= 0);

    if( min_capacity <= a->capacity ) {
        return;
    }

    /* the new capacity is multiple of ARRAY_CAPACITY_INCREMENT */
    new_capacity = min_capacity + a->capacity_inc - min_capacity % a->capacity_inc;
    /*DEBUG(("capacity+: %d", new_capacity));*/

    new_values = (void**)malloc( new_capacity * sizeof(void*) );
    if( new_values == NULL ) {
        ERROR(("Failed to allocate memory for array expansion"));
        return;
    }
    memcpy( new_values,
            a->values,
            a->size * sizeof(void*) );
    if( a->values ) {
        free( a->values );
    }
    a->values = new_values;
    a->capacity = new_capacity;
}

/**
 * Shrinks the capacity of the array list to be equal with it's size.
 **/
void array_shrink( array _a )
{
    _array* a = (_array*) _a;
    void** new_values;
    ASSERT(a);

    if( a->capacity == a->size ) {
        return;
    }
    /*DEBUG(("capacity-: %d", a->size));*/

    new_values = (void**)malloc( a->size * sizeof(void*) );
    if( new_values == NULL ) {
        ERROR(("Failed to allocate memory for array shrinking"));
        return;
    }
    memcpy( new_values,
            a->values,
            a->size * sizeof(void*) );
    if( a->values ) {
        free( a->values );
    }
    a->values = new_values;
    a->capacity = a->size;
}


/**
 * Returns the index of the first value matching the specified key, or -1 if
 * no match was found.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 **/
int array_indexof( const array _a,
                       compare_values_fn compare,
                       const void* key )
{
    const _array* a = (const _array*) _a;
    int i;
    ASSERT(a);
    ASSERT(compare);

    for (i = 0; i < a->size; ++i) {
        if(! compare(key, a->values[i])) {
             return i; /* found */
        }
    }
    return -1; /* not found */
}

/**
 * Returns the index of the last value matching the specified key, or -1 if
 * no match was found.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 **/
int array_last_indexof( const array _a,
                            compare_values_fn compare,
                            const void* key )
{
    const _array* a = (const _array*) _a;
    int i;
    ASSERT(a);
    ASSERT(compare);

    for (i = a->size - 1; i >= 0; --i) {
        if(! compare(key, a->values[i])) {
             return i; /* found */
        }
    }
    return -1; /* not found */
}

/**
 * Returns the number of array values that matches the specified key.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 **/
int array_count_matches( const array _a,
                             compare_values_fn compare,
                             const void* key )
{
    const _array* a = (_array*) _a;
    int i, count;
    ASSERT(a);
    ASSERT(compare);

    for (i = count = 0; i < a->size; ++i) {
        if(! compare(key, a->values[i])) {
            ++count; /* match fount */
        }
    }
    return count;
}

/**
 * Sorts the values from specified array list in ascending order according to a
 * comparison function referenced to by compare, which is called with two arguments
 * that point to the values being compared.
 * The comparison function must return an integer less than, equal to, or greater
 * than zero if the first argument is considered to be respectively less than, equal
 * to, or greater than the second.  If two members compare as equal, their order in
 * the sorted array is undefined.
 **/
void array_sort( array _a,
                     compare_values_fn compare )
{
    _array* a = (_array*) _a;
    ASSERT(a);
    ASSERT(compare);

    _array_qsort( a->values, 0, a->size - 1, compare );
}

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
int array_binary_search( const array _a,
                             compare_values_fn compare,
                             const void* key )
{
    const _array* a = (const _array*) _a;
    int left, right, mid, comp;
    ASSERT(a);
    ASSERT(compare);

    left = 0;
    right = a->size - 1;
    while( left < right ) {
        mid = left + (right - left) / 2;
        comp = compare( key, a->values[mid] );

        if( comp < 0 ) {
            right = mid - 1;
        } else if( comp > 0 ) {
            left = mid + 1;
        } else {
            return mid; /* found */
        }
    }
    return -1; /* not found */
}


/**
 * Process all the values from the array list using the function referenced by
 * 'process'.
 **/
void array_process( array _a,
                        process_value_fn process )
{
    _array* a = (_array*) _a;
    int i;
    ASSERT(a);
    ASSERT(process);

    for (i = 0; i < a->size; ++i) {
        process( a->values[i] );
    }
}

/**
 * Process the specified range (inclusive) of values from the array list using
 * the function referenced by 'process'.
 **/
void array_process_range( array _a,
                              process_value_fn process,
                              int index_from,
                              int index_to )
{
    _array* a = (_array*) _a;
    int i;
    ASSERT(a);
    ASSERT(process);
    ASSERT(index_from >= 0 && index_from < a->size);
    ASSERT(index_to >= 0 && index_to < a->size);
    ASSERT(index_from <= index_to);

    for (i = index_from; i <= index_to; ++i) {
        process( a->values[i] );
    }
}


/* Internal functions */
static void _array_qsort(void** values, int left, int right, compare_values_fn compare)
{
    int old_l = left, old_r = right, mid;
    void *pivot = values[left], *temp;

    while( left < right )
    {
        while( compare( values[right], pivot ) >= 0 && left < right ) {
            --right;
        }
        if( left != right ) {
            temp = values[left];
            values[left] = values[right];
            values[right] = temp;
            ++left;
        }

        while( compare( values[left], pivot ) <= 0 && left < right ) {
            ++left;
        }
        if( left != right ) {
            temp = values[right];
            values[right] = values[left];
            values[left] = temp;
            --right;
        }
    }
    values[left] = pivot;

    mid = left;
    left = old_l;
    right = old_r;
    if( left < mid ) {
        _array_qsort( values, left, mid - 1, compare );
    }
    if( right > mid ) {
        _array_qsort( values, mid + 1, right, compare );
    }
}

static void _array_autoshrink(_array* a)
{
    int new_capacity;
    void** new_values;
    ASSERT(a);

    if( a->capacity - a->size >= a->capacity_inc * 2 ) {
        new_capacity = a->size + a->capacity_inc - a->size % a->capacity_inc;
        /*DEBUG(("capacity-: %d", new_capacity));*/

        new_values = (void**)malloc( new_capacity * sizeof(void*) );
        memcpy( new_values,
                a->values,
                a->size * sizeof(void*) );
        if( a->values ) {
            free( a->values );
        }
        a->values = new_values;
        a->capacity = new_capacity;
    }
}
