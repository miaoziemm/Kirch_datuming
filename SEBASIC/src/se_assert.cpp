#include "../include/se_assert.h"


void assertion_handler( const char* condition,
                            const char* file,
                            const int line,
                            const char* format,
                            ... )
{
    va_list args;
    char msg[1024];
    char* p;

    fflush( stdout );
    fflush( stderr );
    p = msg;
    if( condition != NULL )
        p += snprintf(p, sizeof(msg)-(p-msg), "Assertion [%s] failed", condition);
    else
        p += snprintf(p, sizeof(msg)-(p-msg), "Assertion failed");

    if( file != NULL )
        p += snprintf(p, sizeof(msg)-(p-msg), " at %s:%d.\n", file, line);
    else
        p += snprintf(p, sizeof(msg)-(p-msg), ".\n");

    if( format != NULL ) {

        va_start( args, format );

        p += vsnprintf( p, sizeof(msg)-(p-msg), format, args );

        /* if format ends with ':', print system information */
        if( format[0] != '\0' && format[ strlen(format) - 1 ] == ':' ) {
            int err = errno;
            p += snprintf( p, sizeof(msg)-(p-msg), "errno=%d - %s", err, strerror(err) );
        }

        va_end(args);
    }

    ERROR((msg));


    exit(1);
}