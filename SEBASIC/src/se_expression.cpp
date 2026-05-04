#include "../include/se_expression.h"
#include "se_math_physics.h"


#define EXTRACE(x) TRACE(x)

/** All supported operations */
typedef enum {
    SE_EXP_OP_NOOP,
    SE_EXP_OP_PLUS,
    SE_EXP_OP_MINUS,
    SE_EXP_OP_MINUS_U,
    SE_EXP_OP_MULTIPLY,
    SE_EXP_OP_DIVIDE,
    SE_EXP_OP_MOD,
    SE_EXP_OP_POW,
    SE_EXP_OP_EQ,
    SE_EXP_OP_NE,
    SE_EXP_OP_LE,
    SE_EXP_OP_LT,
    SE_EXP_OP_GE,
    SE_EXP_OP_GT,
    SE_EXP_OP_NEGATE,
    SE_EXP_OP_VALUE
} se_exp_operation;

/** 
 * Internal structure, used to keep a value node. 
 * A value node can be a number (constant), a variable in which case
 * the value is taken form the variable array (see the eval function), or
 * a function (that takes a double and returns a double), in which case the
 * value is the one returned by the function.
 **/
typedef struct se_exp_value_s {
    enum { 
        se_exp_number_type, 
        se_exp_function_type, 
        se_exp_variable_type
    } type;

    union {
        double    n;
        se_exp_function  f;
        int v;
    } val;
} se_exp_value_t;

/** Internal structure used to keep the parsing tree. */
typedef struct se_exp_node_t_s {
    struct se_exp_node_t_s* left;
    struct se_exp_node_t_s* right;

    se_exp_operation op;
    se_exp_value_t val;
} se_exp_node_t;


/** This is created using the parse method.
 * Use the eval method evaluates it to a double value.
 */
typedef struct se_exp_expression_s {
    /** ASCII representation of the expression. */
    char* str;

    /** In case of parsing error it points to an error message.
     * Use have_error(expression*) to check for error and 
     * get_error_message(expression*) to get the error message.
     */
    const char* error;

    /** The current position durring parsing.
     * In case of parsing error it points to the position the error was found.
     * Don't use it directly. Use have_error and get_error_pos instead. */
    size_t pos;

    /** The parsed tree, used to evaluate the expression.
     * Internal use only. */
    se_exp_node_t* tree;


    se_exp_get_var get_var;
    void* get_var_arg;

    se_exp_get_func get_func;
    void* get_func_arg;

    se_exp_eval_par eval_var;

} se_exp_expression_t;




/** All "built-in" functions */
typedef struct se_exp_named_function_s {
    const char* name;
    se_exp_function f;
} se_exp_named_function;

static se_exp_named_function built_in_functions[] = {
    { "sin", sin },
    { "cos", cos },
    { "asin", asin },
    { "acos", acos },
    { "tan", tan },
    { "atan", atan },
    { "exp", exp },
    { "log", log },
    { "log10", log10 },
    { "sqrt", sqrt },
    { "cbrt", cbrt },
    { "fabs", fabs },
    { "ceil", ceil },
    { "floor", floor },
    { "erf", erf },
    { "erfc", erfc },
    { "bessel_j0", se_bessel_j0 },
    { "bessel_j1", se_bessel_j1 },
    { "bessel_y0", se_bessel_y0 },
    { "bessel_y1", se_bessel_y1 },
    { "lgamma", lgamma },
    { NULL, NULL }
};



static se_exp_node_t* create_node( void );
static void destroy_node( se_exp_node_t* n );

static double eval_node(se_exp_expression_t* thisone, se_exp_node_t* n, void* vars);

static se_exp_node_t* parse_expression( se_exp_expression_t* thisone );
static se_exp_node_t* parse_term( se_exp_expression_t* thisone );
static se_exp_node_t* parse_pow( se_exp_expression_t* thisone );
static se_exp_node_t* parse_logical( se_exp_expression_t* thisone );
static se_exp_node_t* parse_operand( se_exp_expression_t* thisone );
static se_exp_node_t* parse_factor( se_exp_expression_t* thisone );
static se_exp_node_t* parse_identifier( se_exp_expression_t* thisone );
static se_exp_node_t* parse_number( se_exp_expression_t* thisone );

static int is_space(char c);
static int is_identifier_part(char c);


typedef struct simple_list_of_vars_s {
    const char** vars;
    int nvars;
} simple_list_of_vars;


static se_exp_function slv_get_func( void* arg, const char* str, size_t len )
{
    se_exp_named_function* functions = (se_exp_named_function*)arg;
    int i;

    for( i = 0; functions[i].name != NULL; ++i ) {
        if( len == strlen(functions[i].name) && 
            strncmp(str, functions[i].name, len) == 0 )
            return functions[i].f;
    }

    return NULL;
}

static int slv_get_var( void* arg,
                        const char* str,
                        size_t len )
{
    simple_list_of_vars* slv = (simple_list_of_vars*)arg;
    int i;

    for(i = 0; i < slv->nvars; i++) {
        if( slv->vars[i] &&
            len == strlen(slv->vars[i]) &&
            strncmp(str, slv->vars[i], len) == 0 )
            return i;
    }

    return -1;
}

static double slv_eval_var( void* arg, int i)
{
    return ((double*)arg)[i];
}



/**
 * Parse the expression passed in str. The list of available variables
 * is passed in the vars array.
 */
se_exp_t se_exp_parse_simple( const char* str, 
                              const char* vars[],
                              int nvars )
{
    simple_list_of_vars slv;
    slv.vars = vars;
    slv.nvars = nvars;
    return se_exp_parse( str, slv_get_var, &slv,
                         slv_get_func, built_in_functions,
                         slv_eval_var );
}

se_exp_t se_exp_parse( const char* str, 
                       se_exp_get_var get_var, void* get_var_arg,
                       se_exp_get_func get_func, void* get_func_arg,
                       se_exp_eval_par eval_var )
{
    se_exp_expression_t* expr;
    expr = (se_exp_expression_t*) malloc( sizeof(se_exp_expression_t) );

    expr->str = strdup(str);

    expr->get_var = get_var;
    expr->get_var_arg = get_var_arg;

    if( get_func != NULL ) {
        expr->get_func = get_func;
        expr->get_func_arg = get_func_arg;
    } else {
        expr->get_func = slv_get_func;
        expr->get_func_arg = built_in_functions;
    }

    expr->eval_var = eval_var;

    expr->error = NULL;
    expr->pos = 0;
    expr->tree = NULL;

    expr->tree = parse_expression( expr );

    return (se_exp_t)expr;
}


double se_exp_eval(se_exp_t _exp, void* eval_var_arg)
{
    se_exp_expression_t* exp = (se_exp_expression_t*)_exp;
    /*TODO: asert(expr != NULL) + exp->tree */
    return eval_node(exp, exp->tree, eval_var_arg);
}

/**
 * Free used memory.
 */
void se_exp_destroy_expression(se_exp_t _exp)
{
    se_exp_expression_t* exp = (se_exp_expression_t*)_exp;
    if(exp == NULL) return;

    if( exp->tree != NULL )
        destroy_node(exp->tree);

    free(exp->str);
    free(exp);
}


/**
 * Test if there is an error on an expression. Call thisone after calling parse.
 */
int se_exp_have_error(se_exp_t _exp)
{
    se_exp_expression_t* exp = (se_exp_expression_t*)_exp;
    return ( exp == NULL || exp->error != NULL );
}


/**
 * If have_error returns true, a message describing the error can be retreived
 * using thisone function.
 */
const char* se_exp_error_message(se_exp_t _exp)
{
    se_exp_expression_t* exp = (se_exp_expression_t*)_exp;
    if( exp != NULL ) {
        return exp->error;
    } else {
        return "Unknown error";
    }
}

/**
 * If have_error returns true, the position related to the parsing error can
 * be retreived using thisone function.
 */
size_t se_exp_error_pos(se_exp_t _exp)
{
    se_exp_expression_t* exp = (se_exp_expression_t*)_exp;
    if( exp != NULL ) {
        return exp->pos;
    } else {
        return 0;
    }
}

/**
 * Returns the initial string used to parse thisone expression.
 */
const char* se_exp_initial_str(se_exp_t _exp)
{
    se_exp_expression_t* exp = (se_exp_expression_t*)_exp;
    if( exp != NULL ) {
        return exp->str;
    } else {
        return "Expression was not parsed";
    }
}

static double eval_node(se_exp_expression_t* thisone, se_exp_node_t* n, void* vars)
{
    switch(n->op) {
    case SE_EXP_OP_PLUS:
        return ( eval_node(thisone, n->left,  vars) + 
                 eval_node(thisone, n->right, vars) );

    case SE_EXP_OP_MINUS:
        return ( eval_node(thisone, n->left,  vars) - 
                 eval_node(thisone, n->right, vars) );

    case SE_EXP_OP_MULTIPLY:
        return ( eval_node(thisone, n->left,  vars) * 
                 eval_node(thisone, n->right, vars) );

    case SE_EXP_OP_DIVIDE:
        return ( eval_node(thisone, n->left,  vars) / 
                 eval_node(thisone, n->right, vars) );

    case SE_EXP_OP_MOD: /* TODO: thisone have to be replaced with something better */
        return (double)( (int64_t)eval_node(thisone, n->left,  vars) % 
                         (int64_t)eval_node(thisone, n->right, vars) );

    case SE_EXP_OP_POW:
        return pow( eval_node(thisone, n->left, vars), 
                    eval_node(thisone, n->right, vars) );

    case SE_EXP_OP_EQ:
        return fabs( eval_node(thisone, n->left, vars) - 
                     eval_node(thisone, n->right, vars) ) < EPSILON ? 
            1.0 : 0.0;

    case SE_EXP_OP_NE:
        return fabs( eval_node(thisone, n->left, vars) - 
                     eval_node(thisone, n->right, vars) ) >= EPSILON ? 
            1.0 : 0.0;

    case SE_EXP_OP_LE:
        return ( eval_node(thisone, n->left, vars) <= 
                 eval_node(thisone, n->right, vars) ) ? 
            1.0 : 0.0;

    case SE_EXP_OP_LT:
        return ( eval_node(thisone, n->left, vars) <  
                 eval_node(thisone, n->right, vars) ) ? 
            1.0 : 0.0;

    case SE_EXP_OP_GE:
        return ( eval_node(thisone, n->left, vars) >= 
                 eval_node(thisone, n->right, vars) ) ? 
            1.0 : 0.0;

    case SE_EXP_OP_GT:
        return ( eval_node(thisone, n->left, vars) >
                 eval_node(thisone, n->right, vars) ) ? 
            1.0 : 0.0;

    case SE_EXP_OP_NEGATE:
        return fabs( eval_node(thisone, n->left, vars) ) < EPSILON ? 1.0 : 0.0 ;
        
    case SE_EXP_OP_MINUS_U:
        return ( - eval_node(thisone, n->left, vars) );

    case SE_EXP_OP_VALUE:
        switch( n->val.type ) {
        case se_exp_value_t::se_exp_number_type:
            return n->val.val.n;

        case se_exp_value_t::se_exp_function_type:
            return n->val.val.f( eval_node(thisone, n->left, vars) );

        case se_exp_value_t::se_exp_variable_type:
            return thisone->eval_var(vars, n->val.val.v);
        }

    case SE_EXP_OP_NOOP:
        break;
    }

    /** TODO: assert(false) */
    fprintf(stderr, "Internal error on eval\n");
    exit(1);
    return 0.0;
}


/**
 * EXPRESSION = TERM[+|-TERM]*
 */
static se_exp_node_t* parse_expression( se_exp_expression_t* thisone )
{
    se_exp_node_t* n;

    EXTRACE(("Parse expression:[%s], pos=%d", thisone->str+thisone->pos, thisone->pos));

    if( ( n = parse_term(thisone) ) == NULL ) goto error;

    /* skip leading spaces */
    while( is_space(thisone->str[thisone->pos]) ) thisone->pos++;

    for(;;) {
        se_exp_operation op;
        se_exp_node_t* nl;
        se_exp_node_t* nr;

        switch( thisone->str[thisone->pos] ) {
        case '+': 
            op = SE_EXP_OP_PLUS; 
            break;
        case '-':
            op = SE_EXP_OP_MINUS; 
            break;

        default:
            return n;
        }

        thisone->pos++;

        nl = n;
        if( ( nr = parse_term(thisone) ) == NULL ) {
            destroy_node(n);
            goto error;
        }

        n = create_node();
        n->left  = nl;
        n->right = nr;
        n->op = op;
    }

 error:
    if(thisone->error == NULL) thisone->error = "Expected term";
    return NULL;
}


/**
 * TERM = POW[*|/POW]
 */
static se_exp_node_t* parse_term( se_exp_expression_t* thisone )
{
    se_exp_node_t* n;

    EXTRACE(("Parse term      :[%s], pos=%d", thisone->str+thisone->pos, thisone->pos));

    if( ( n = parse_pow(thisone) ) == NULL ) goto error;

    /* skip leading spaces */
    while( is_space(thisone->str[thisone->pos]) ) thisone->pos++;

    for(;;) {
        se_exp_operation op;
        se_exp_node_t* nl;
        se_exp_node_t* nr;

        switch( thisone->str[thisone->pos] ) {
        case '*': 
            op = SE_EXP_OP_MULTIPLY; 
            break;
        case '/':
            op = SE_EXP_OP_DIVIDE; 
            break;
        case '%':
            op = SE_EXP_OP_MOD;
            break;

        default:
            return n;
        }

        thisone->pos++;

        nl = n;
        if( ( nr = parse_pow(thisone) ) == NULL ) {
            destroy_node(n);
            goto error;
        }
            
        n = create_node();
        n->left  = nl;
        n->right = nr;
        n->op = op;
    }

 error:
    if(thisone->error == NULL) thisone->error = "Expected pow";
    return NULL;
}


/**
 * POW = LOGICAL[^LOGICAL]
 */
static se_exp_node_t* parse_pow( se_exp_expression_t* thisone )
{
    se_exp_node_t* n;

    EXTRACE(("Parse pow       :[%s], pos=%d", thisone->str+thisone->pos, thisone->pos));

    if( ( n = parse_logical(thisone) ) == NULL ) goto error;

    /* skip leading spaces */
    while( is_space(thisone->str[thisone->pos]) ) thisone->pos++;

    for(;;) {
        se_exp_operation op;
        se_exp_node_t* nl;
        se_exp_node_t* nr;

        if ( thisone->str[thisone->pos] == '^' ) {
            op = SE_EXP_OP_POW; 
        } else {
            return n;
        }

        thisone->pos++;

        nl = n;
        if( ( nr = parse_logical(thisone) ) == NULL ) {
            destroy_node(n);
            goto error;
        }
            
        n = create_node();
        n->left  = nl;
        n->right = nr;
        n->op = op;
    }

 error:
    if(thisone->error == NULL) thisone->error = "Expected logical";
    return NULL;
}

/**
 * LOGICAL = OPERAND[<|>OPERAND]
 */
static se_exp_node_t* parse_logical( se_exp_expression_t* thisone )
{
    se_exp_node_t* n;

    EXTRACE(("Parse logical   :[%s], pos=%d", thisone->str+thisone->pos, thisone->pos));

    if( ( n = parse_operand(thisone) ) == NULL ) goto error;

    /* skip leading spaces */
    while( is_space(thisone->str[thisone->pos]) ) thisone->pos++;

    for(;;) {
        se_exp_operation op;
        se_exp_node_t* nl;
        se_exp_node_t* nr;

        switch( thisone->str[thisone->pos] ) {
        case '>':
            if( thisone->str[thisone->pos + 1] == '=' ) {
                thisone->pos++;
                op = SE_EXP_OP_GE;
            } else {
                op = SE_EXP_OP_GT;
            }
            break;

        case '<':
            if( thisone->str[thisone->pos + 1] == '=' ) {
                thisone->pos++;
                op = SE_EXP_OP_LE;
            } else {
                op = SE_EXP_OP_LT;
            }
            break;

        case '=':
            if( thisone->str[thisone->pos + 1] == '=' ) {
                thisone->pos++;
            }
            op = SE_EXP_OP_EQ;
            break;

        case '!':
            if( thisone->str[thisone->pos + 1] == '=' ) {
                thisone->pos++;
                op = SE_EXP_OP_NE;
            } else {
                return n;
            }
            break;

        default:
            return n;
        }

        thisone->pos++;

        nl = n;
        if( ( nr = parse_operand(thisone) ) == NULL ) {
            destroy_node(n);
            goto error;
        }
            
        n = create_node();
        n->left  = nl;
        n->right = nr;
        n->op = op;
    }

 error:
    if(thisone->error == NULL) thisone->error = "Expected operand";
    return NULL;
}

/**
 * OPERAND = [!|-]OPERAND|FACTOR
 */
static se_exp_node_t* parse_operand( se_exp_expression_t* thisone )
{
    se_exp_node_t* n;
    se_exp_node_t* nl;
    se_exp_operation op;

    EXTRACE(("Parse operand   :[%s], pos=%d", thisone->str+thisone->pos, thisone->pos));

    /* skip leading spaces */
    while( is_space(thisone->str[thisone->pos]) ) thisone->pos++;

    switch( thisone->str[thisone->pos] ) {
    case '!':
        op = SE_EXP_OP_NEGATE;
        break;

    case '-':
        op = SE_EXP_OP_MINUS_U;
        break;

    default:
        return parse_factor(thisone);
    }


    thisone->pos++;
    
    nl = parse_operand(thisone);

    n = create_node();
    n->left  = nl;
    n->right = NULL;
    n->op = op;

    return n;
}

/**
 * FACTOR=IDENTIFIER|(EXPRESSION)
 */
static se_exp_node_t* parse_factor( se_exp_expression_t* thisone )
{
    se_exp_node_t* n;

    EXTRACE(("Parse factor    :[%s], pos=%d", thisone->str+thisone->pos, thisone->pos));

    /* skip leading spaces */
    while( is_space(thisone->str[thisone->pos]) ) thisone->pos++;
    
    if( thisone->str[thisone->pos] == '(' ) {
        thisone->pos++;

        if( ( n = parse_expression(thisone) ) == NULL ) {
            if( thisone->error == NULL ) thisone->error = "Expected expression";
            return NULL;
        }

        if( thisone->str[ thisone->pos] != ')' ) {
            destroy_node(n);
            if( thisone->error == NULL ) thisone->error = "Expected ')'";
            return NULL;
        }

        thisone->pos++;
    }

    else {
        n = parse_identifier(thisone);
    }


    return n;
}

/** Identifier is number, function or variable
 */
static se_exp_node_t* parse_identifier( se_exp_expression_t* thisone )
{
    se_exp_node_t* n;

    EXTRACE(("Parse identifier:[%s], pos=%d", thisone->str+thisone->pos, thisone->pos));

    /* skip leading spaces */
    while( is_space(thisone->str[thisone->pos]) ) thisone->pos++;

    /* first try for a number */
    if( ( n = parse_number(thisone) ) != NULL) {
        return n;
    }

    /* function or variable */
    if( is_identifier_part(thisone->str[thisone->pos]) ) { 
        
        size_t p = thisone->pos; /* save start position */

        n = create_node();
        n->op = SE_EXP_OP_VALUE;

        while( is_identifier_part(thisone->str[thisone->pos]) ) thisone->pos++;
        
        if( ( n->val.val.f = thisone->get_func(thisone->get_func_arg, thisone->str + p, thisone->pos - p) )
            != NULL ) {
            
            if( ( n->left = parse_factor(thisone) ) == NULL ) {
                if(thisone->error == NULL) thisone->error = "Expected expression";
                goto error;
            }

            n->right = NULL;
            n->val.type = se_exp_value_t::se_exp_function_type;
        }

        else if( ( n->val.val.v = thisone->get_var( thisone->get_var_arg,
                                                 thisone->str + p, 
                                                 thisone->pos - p ) ) != -1 ) {
            n->left  = NULL;
            n->right = NULL;
            n->val.type = se_exp_value_t::se_exp_variable_type;
        }

        else {
            if( thisone->error == NULL ) 
                thisone->error = "Unknown variable or function";
            goto error;
        }
        
    }

    else {
        if( thisone->error == NULL ) 
            thisone->error = "Expected variable or function";
        goto error;
    }

    /* skip trailing spaces */
    while( is_space(thisone->str[thisone->pos]) ) thisone->pos++;

    return n;

 error:
    destroy_node(n);
    return NULL;
}

/** Parse a number, if possible. */
static se_exp_node_t* parse_number( se_exp_expression_t* thisone )
{
    const char* ptr1 = thisone->str + thisone->pos;
    char* ptr2;
    double dval = strtod(ptr1, &ptr2);

    EXTRACE(("Parse number    :[%s], pos=%d", thisone->str+thisone->pos, thisone->pos));

    if( ptr1 != ptr2 ) {
        se_exp_node_t* n = create_node();
        n->op = SE_EXP_OP_VALUE;
        n->left  = NULL;
        n->right = NULL;
        n->val.val.n = dval;
        n->val.type = se_exp_value_t::se_exp_number_type;
        thisone->pos += ( ptr2 - ptr1 );

        /* skip trailing spaces */
        while( is_space(thisone->str[thisone->pos]) ) thisone->pos++;

        return n;
    } else {
        /* not a number */
        return NULL;
    }
}


/** Allocate and initialize a node. */
static se_exp_node_t* create_node()
{
    se_exp_node_t* n = (se_exp_node_t*)malloc(sizeof(se_exp_node_t));

    n->left = NULL;
    n->right = NULL;
    n->op = SE_EXP_OP_NOOP;

    return n;
}

/** Deallocate a node. It deallocates recursively all childrens */
static void destroy_node(se_exp_node_t* n)
{
    if( n != NULL ) {
        destroy_node(n->left);
        destroy_node(n->right);
        free(n);
    }
}

/* Next functions ASCII specific ... */
static int is_space(char c)
{
    return ( c == ' '  || 
             c == '\t' ||
             c == '\n' || 
             c == '\r' );
}

static int is_identifier_part(char c)
{ 
    return ( (c >= 'a' && c <= 'z') ||
             (c >= 'A' && c <= 'Z') || 
             (c >= '0' && c <= '9') || 
             c == '_' || c == '$' );
}


/**********************************************************************/
/************************** Testing ***********************************/

static const char* op_to_string(se_exp_operation op)
{
    switch(op) {
    case SE_EXP_OP_NOOP: return "SE_EXP_OP_NOOP";
    case SE_EXP_OP_PLUS: return "SE_EXP_OP_PLUS";
    case SE_EXP_OP_MINUS: return "SE_EXP_OP_MINUS";
    case SE_EXP_OP_MINUS_U: return "SE_EXP_OP_MINUS_U";
    case SE_EXP_OP_MULTIPLY: return "SE_EXP_OP_MULTIPLY";
    case SE_EXP_OP_DIVIDE: return "SE_EXP_OP_DIVIDE";
    case SE_EXP_OP_MOD: return "SE_EXP_OP_MOD";
    case SE_EXP_OP_POW: return "SE_EXP_OP_POW";
    case SE_EXP_OP_EQ: return "SE_EXP_OP_EQ";
    case SE_EXP_OP_NE: return "SE_EXP_OP_NE";
    case SE_EXP_OP_LE: return "SE_EXP_OP_LE";
    case SE_EXP_OP_LT: return "SE_EXP_OP_LT";
    case SE_EXP_OP_GE: return "SE_EXP_OP_GE";
    case SE_EXP_OP_GT: return "SE_EXP_OP_GT";
    case SE_EXP_OP_NEGATE: return "SE_EXP_OP_NEGATE";
    case SE_EXP_OP_VALUE: return "SE_EXP_OP_VALUE";
    }

    return "OPERATION UNKNOWN";
}

static void print_value(FILE* f, se_exp_value_t v)
{
    switch(v.type) {
    case se_exp_value_t::se_exp_number_type:
        fprintf(f, "number:%f", v.val.n);
        break;
    case se_exp_value_t::se_exp_function_type:
        fprintf(f, "function");
        break;
    case se_exp_value_t::se_exp_variable_type:
        fprintf(f, "variable:%d", v.val.v);
        break;
    default:
        fprintf(f, "TYPE unknown");
        break;
    }
}

static void print_node(FILE* f, se_exp_node_t* n, int indent)
{
    int i;
    if(n == NULL) {
        for(i =0; i<indent; i++) fprintf(f, " ");
        fprintf(f, "NULL\n");
    } else {
        for(i =0; i<indent; i++) fprintf(f, " ");
        fprintf(f, "op=%s", op_to_string(n->op));
        if(n->op == SE_EXP_OP_VALUE) {
            fprintf(f, ":");
            print_value(f, n->val);
        }
        fprintf(f, "\n");
        print_node(f, n->left, indent+2);
        print_node(f, n->right, indent+2);
    }
}

void se_exp_print_expression(FILE* f, se_exp_t _exp)
{
    se_exp_expression_t* exp = (se_exp_expression_t*)_exp;

    fprintf(f, "Expression structure\n");
    print_node(f, exp->tree, 0);
    fprintf(f, "\n");
}


void se_exp_print_formated_error(se_exp_t e, int verb, int print_no_error)
{
    if(!se_exp_have_error(e)) {
        if(!print_no_error) return;
        INFOV((verb, "%s has no errors", se_exp_initial_str(e)));
    } else {
        int i;
        char* spaces = alloc1char(se_exp_error_pos(e) + 1);
        INFOV((verb, "%s : %s", se_exp_initial_str(e), se_exp_error_message(e)));
        for(i = 0; i < (int)se_exp_error_pos(e)-1; i++)
            spaces[i] = ' ';
        spaces[i++] = '^';
        spaces[i++] = 0;
        INFOV((verb, "%s", spaces));
        free(spaces);
    }
}
