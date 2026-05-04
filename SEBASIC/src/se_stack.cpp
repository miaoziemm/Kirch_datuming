#include "../include/se_stack.h"


void se_stack_clear( se_stack_t* s, int free_values )
{
    se_stack_node_t* curr_node;
    se_stack_node_t* next_node;

    ASSERT(s);

    /* for each node */
    curr_node = s->head;
    while(curr_node) {
        /* get next node */
        next_node = curr_node->next;

        /* free the value if requested so */
        if(free_values)
            free(curr_node->value);

        /* free the node */
        free(curr_node);

        /* move to next node */
        curr_node = next_node;
    }
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
                        const void* key )
{
    se_stack_node_t* curr_node;
    se_stack_node_t* prev_node;
    void* value = NULL;

    ASSERT(s);
    ASSERT(compare);

    /* for each node */
    curr_node = prev_node = s->head;
    while (curr_node) {
        /* check for a match */
        if(! compare(key, curr_node->value)) {
            /* keep the value */
            value = curr_node->value;

            /* remove (unlink) the node */
            if (curr_node == s->head) {
                /* removing head node */
                prev_node = curr_node->next;
                free(curr_node);
                curr_node = prev_node;
                s->head = curr_node;
            } else {
                /* removing non-head node */
                prev_node->next = curr_node->next;
                free(curr_node);
                curr_node = prev_node->next;
            }
            if(! curr_node) {
                /* removing tail node */
                s->tail = prev_node;
            }

            s->size--;
            break;
        }
        /* move to next node */
        prev_node = curr_node;
        curr_node = curr_node->next;
    }

    /* return the removed value (or NULL if no value was removed) */
    return value;
}

/**
 * Removes from the stack all values that matches the specified key.
 * Returns the number of removed values.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 * free() the values if 'free_values' is true (non-zero)
 **/
int se_stack_remove_matches( se_stack_t* s,
                              compare_values_fn compare,
                              const void* key,
                              int free_values )
{
    se_stack_node_t* curr_node;
    se_stack_node_t* prev_node;
    int count = 0;

    ASSERT(s);
    ASSERT(compare);

    /* for each node */
    curr_node = prev_node = s->head;
    while (curr_node) {
        /* check for a match */
        if(! compare(key, curr_node->value)) {
            /* free the value if requested so */
            if (free_values) {
                free(curr_node->value);
            }

            /* remove (unlink) the node */
            if (curr_node == s->head) {
                /* removing head node */
                prev_node = curr_node->next;
                free(curr_node);
                curr_node = prev_node;
                s->head = curr_node;
            } else {
                /* removing non-head node */
                prev_node->next = curr_node->next;
                free(curr_node);
                curr_node = prev_node->next;
            }
            if(! curr_node) {
                /* removing tail node */
                s->tail = prev_node;
            }

            /* update counts */
            s->size--;
            ++count;
        } else {
            /* move to next node */
            prev_node = curr_node;
            curr_node = curr_node->next;
        }
    }

    /* return removed node count */
    return count;
}

/**
 * Returns true (non-zero) if the stack contains the specified key, or false
 * (zero) no match was found.
 * The key and values are compared used the comparision function referenced
 * by 'compare'. The comparison function must return zero if it's arguments are
 * equal, or a non-zero value otherwise.
 **/
int se_stack_contains( const se_stack_t* s,
                        compare_values_fn compare,
                        const void* key )
{
    se_stack_node_t* curr_node;

    ASSERT(s);
    ASSERT(compare);

    /* for each node */
    curr_node = s->head;
    while(curr_node) {
        /* check for a match */
        if( ! compare(key, curr_node->value) ) {
            return 1;   /* found */
        }
        curr_node = curr_node->next;
    }
    return 0;   /* not found */
}

int se_stack_count_matches( const se_stack_t* s,
                             compare_values_fn compare,
                             const void* key )
{
    se_stack_node_t* curr_node;
    int count = 0;

    ASSERT(s);
    ASSERT(compare);

    /* for each node */
    curr_node = s->head;
    while (curr_node) {
        /* check for a match */
        if (! compare(key, curr_node->value)) {
            ++count;
        }
        curr_node = curr_node->next;
    }
    return count;
}

/**
 * Process all the values from the stack using the function referenced by
 * 'process'.
 **/
void se_stack_process( se_stack_t* s,
                        process_value_fn process )
{
    se_stack_node_t* curr_node;

    ASSERT(s);
    ASSERT(process);

    /* for each node */
    curr_node = s->head;
    while (curr_node) {
        process(curr_node->value);
        curr_node = curr_node->next;
    }
}
