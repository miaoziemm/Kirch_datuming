#ifndef SE_EXPRESSION_H
#define SE_EXPRESSION_H
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "se_basic_math.h"
#include "se_type.h"
#include "se_alloc.h"
#include "se_log.h"
#include "se_assert.h"


typedef void* se_exp_t;

/** A generic function that takes a double and return a double. It can be
 * sin, cos, tan etc. */
typedef double (*se_exp_function)(double param);

/** Generic function to check if a variable name can be evaluated.
    This function is called durring the parsing; the returned integer 
    would be passed during the eval phase to get the value for a variable */
typedef int (*se_exp_get_var)( void*, const char* str, size_t len );

/** generic function to check if a function name is legal and to
    return a pointer to the actual function */
typedef se_exp_function (*se_exp_get_func)(void*, const char* str, size_t len);

/** Called during the evaluation phase to get the actual values for variables.
    The int is the same as the one return buy se_exp_get_var during parsing. */
typedef double (*se_exp_eval_par)( void*, int );





/**
 * Parse the expression passed in str. The function uses user provided
 * information to get the variables, functions, etc.
 */
se_exp_t se_exp_parse( const char* str, 
                       se_exp_get_var get_var, void* get_var_arg,
                       se_exp_get_func get_func, void* get_func_arg,
                       se_exp_eval_par eval_var );

/**
 * Evaluate an already parsed expression.
 */
double se_exp_eval(se_exp_t exp, void* eval_var_arg);


/**
 * Parse the expression passed in str. The list of available variables
 * is passed in the vars array.
 */
se_exp_t se_exp_parse_simple( const char* str, 
                              const char* vars[],
                              int nvars );

/**
 * Evaluate an already parsed expression. The current values for variables are
 * passed in the vars array. The dimension of this array should be at least
 * equal to the dimension of the vars array passed to the parse_simple function.
 */
inline double se_exp_eval_simle(se_exp_t exp, double vars[])
{
    return se_exp_eval(exp, vars);
}


/**
 * Free memory.
 */
void se_exp_destroy_expression(se_exp_t exp);

/**
 * Test if there is an error on an expression. Call this after calling parse.
 */
int se_exp_have_error(se_exp_t exp);

/**
 * If have_error returns true, a message describing the error can be retreived
 * using this function.
 */
const char* se_exp_error_message(se_exp_t exp);

/**
 * If have_error returns true, the position related to the parsing error can
 * be retreived using this function.
 */
size_t se_exp_error_pos(se_exp_t exp);


/**
 * Output internal structeure of the expression, after parsing.
 * Useful for debugging.
 */
void se_exp_print_expression(FILE* f, se_exp_t expr);

/**
 * Returns the initial string used to parse this expression.
 */
const char* se_exp_initial_str(se_exp_t _exp);


/**
 * Prints a formated version of the error, with position of the error.
 */
void se_exp_print_formated_error(se_exp_t e, int verb, int print_no_error);


#endif