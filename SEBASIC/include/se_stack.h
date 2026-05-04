#ifndef SE_STACK_H
#define SE_STACK_H
#include "se_type.h"
#include "se_assert.h"
#include "se_alloc.h"
#include "se_util.h"


/** stack internal single-linked node */
typedef struct se_stack_node_s
{
    void* value;
    struct se_stack_node_s* next;
} se_stack_node_t;

/** stack internal structure */
typedef struct
{
    se_stack_node_t* head;
    se_stack_node_t* tail;
    int size;
} se_stack_t;

/**
 * Removes all values from the stack.
 * se_free() all values if 'free_values' is true (non-zero)
 **/
void se_stack_clear( se_stack_t* s, int free_values );

/**
 * Creates an empty stack.
 **/
inline se_stack_t* se_create_stack( void )
{
    return (se_stack_t*)alloc1type_zero(se_stack_t,1);
}


/**
 * Destroys the stack.
 * \c se_free() the values if 'free_values' is true (non-zero)
 **/
inline void se_destroy_stack( se_stack_t* s, int free_values )
{
    ASSERT(s);
    se_stack_clear(s, free_values);
    free(s);
}


/**
 * Push a value on top of the stack.
 */
inline void se_stack_push( se_stack_t* s, void* value )
{
    se_stack_node_t* new_node;
    ASSERT(s);

    /* create a new node */
    new_node = (se_stack_node_t*)malloc( sizeof(*new_node) );
    new_node->value = value;
    new_node->next = s->head;

    /* link the new node */
    if (s->size == 0) {
        s->head = s->tail = new_node;
    } else {
        s->head = new_node;
    }

    s->size++;
}


/**
 * Pop the value from the top of the stack.
 * If the stack is empty the program will exit with ERROR.
 **/
inline void* se_stack_pop( se_stack_t* s )
{
    se_stack_node_t* head_node;
    void* value;

    ASSERT(s);

    if( s->size == 0 ) {
        ERROR(("Cannot pop value from an empty stack"));
    }

    /* keep the value */
    value = s->head->value;

    /* remove (unlink) the top node */
    head_node = s->head;
    s->head = head_node->next;
    free(head_node);
    if( ! s->head ) {
        s->tail = NULL;
    }
    s->size--;

    return value;
}

/**
 * Returns the value from the top of the stack.
 * If the stack is empty the program will exit with ERROR.
 **/
inline void* se_stack_top( se_stack_t* s )
{
    ASSERT(s);
    if (s->size == 0) {
        ERROR(("Cannot query top value from an empty stack"));
    }

    return s->head->value;
}


/**
 * Returns the number of values from the stack.
 **/
inline int se_stack_size( const se_stack_t* s )
{
    ASSERT(s);
    return s->size;
}


/**
 * Removes from the stack the first value that matches the specified key.
 * Returns the removed value (that can be NULL), or NULL if no match was found.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 * Does _not_ call free() for the removed value.
 **/
void* se_stack_remove( se_stack_t* s,
                        compare_values_fn compare,
                        const void* key );

/**
 * Removes from the stack all values that matches the specified key.
 * Returns the number of removed values.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 * se_free() the values if 'free_values' is true (non-zero)
 **/
int se_stack_remove_matches( se_stack_t* s,
                              compare_values_fn compare,
                              const void* key,
                              int free_values );


/**
 * Returns true (non-zero) if the stack contains the specified key, or false
 * (zero) no match was found.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 **/
int se_stack_contains( const se_stack_t* s,
                        compare_values_fn compare,
                        const void* key );

/**
 * Counts all the values from the stack that matches the specified key.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 **/
int se_stack_count_matches( const se_stack_t* s,
                             compare_values_fn compare,
                             const void* key );


/**
 * Process all the values from the stack using the function referenced by
 * 'process'.
 **/
void se_stack_process( se_stack_t* s,
                        process_value_fn process );




#endif