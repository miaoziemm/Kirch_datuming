#ifndef SE_ASSERT_H
#define SE_ASSERT_H
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <errno.h>
#include "se_type.h"
#include "se_log.h"
#include "se_macro.h"

/**
 * Assertion_handler - called when an error is detected.
 * If format ends with ':' it prints also the error message associated
 * with the current value of 'errno'.
 */
void assertion_handler( const char* condition,
                            const char* file,
                            const int line,
                            const char* format,
                            ... );
#define VERIFY(condition) do {                                        \
        if( ! (condition) ) {                                         \
            assertion_handler(#condition, NULL, 0, NULL);         \
        }                                                             \
   } while(0);


#define VERIFYM(condition, msg) do {                                  \
        if( ! (condition) ) {                                         \
            assertion_handler(#condition, NULL, 0, (msg));        \
        }                                                             \
    } while(0);


#define ASSERT(c)        NOOP
#define ASSERTM(c, msg)  NOOP


#endif
