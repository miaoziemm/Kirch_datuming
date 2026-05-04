#ifndef _SE_MESSAGE_H
#define _SE_MESSAGE_H
#include <iostream>
#include <cstdio>
#include <cstdarg>
#include <chrono>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <unistd.h>

#ifdef SE_USE_OMP
#include <omp.h>
#endif

#ifndef NO_LOG_TIME
#include <time.h>
/* #include <locale.h> */
/* #include <langinfo.h> */
#endif
#include "se_assert.h"
#include "se_version.h"

enum LOGType
{
    SE_INFO,
    SE_WARNING,
    SE_ERROR,
    SE_FATAL
};

void LOG_CONSOLE(LOGType type = SE_INFO, const char *format = "", ...);


// copy from Z-terra

#ifdef ERROR
#undef ERROR
#endif

#ifdef WARN
#undef WARN
#endif

#ifdef NOTICE
#undef NOTICE
#endif

#ifdef INFO
#undef INFO
#endif

#ifdef INFOV
#undef INFOV
#endif

#ifdef DEBUG
#undef DEBUG
#endif

#ifdef DEBUGV
#undef DEBUGV
#endif

#ifdef TRACE
#undef TRACE
#endif

#ifdef TRACEV
#undef TRACEV
#endif

#ifdef PROGRESS
#undef PROGRESS
#endif

#ifdef _LOG
#undef _LOG
#endif

/** Figure out the log level ... */

/** Set the log level */
#ifndef SE_LOG_LEVEL /* if it is not defined outside, look for SE_DEBUG def */
#  ifdef SE_DEBUG
#    define SE_LOG_LEVEL 40 /* all up to debugging */
#  else
#    define SE_LOG_LEVEL 30 /* all up to info */
#  endif
#endif

/** Enable printing the source file name and line number */
#ifndef SE_LOG_SRC
#  ifdef SE_DEBUG
#    define SE_LOG_SRC
#  endif
#endif


/** ERROR enabled always */
#ifndef SE_ENABLE_ERROR_LOG
#  define SE_ENABLE_ERROR_LOG
#endif

/** Enable warnings */
#ifndef SE_ENABLE_WARN_LOG
#  if SE_LOG_LEVEL >= 10
#    define SE_ENABLE_WARN_LOG
#  endif
#endif

/** Enable notices */
#ifndef SE_ENABLE_NOTICE_LOG
#  if SE_LOG_LEVEL >= 20
#    define SE_ENABLE_NOTICE_LOG
#  endif
#endif

/** Enable progress: same level as NOTICE */
#ifndef SE_ENABLE_PROGRESS_LOG
#  if SE_LOG_LEVEL >= 20
#    define SE_ENABLE_PROGRESS_LOG
#  endif
#endif

/** Enable info */
#ifndef SE_ENABLE_INFO_LOG
#  if SE_LOG_LEVEL >= 30
#    define SE_ENABLE_INFO_LOG
#  endif
#endif

/** Enable debugging messages  */
#ifndef SE_ENABLE_DEBUG_LOG
#  if SE_LOG_LEVEL >= 40
#    define SE_ENABLE_DEBUG_LOG
#  endif
#endif

/** Enable tracing  */
#ifndef SE_ENABLE_TRACE_LOG
#  if SE_LOG_LEVEL >= 50
#    define SE_ENABLE_TRACE_LOG
#  endif
#endif


/** The 5 types (levels) of logs; These are defined to match the names array. */
typedef enum {
    SE_LOG_ERROR    = 0,
    SE_LOG_WARN     = 1,
    SE_LOG_NOTICE   = 2,
    SE_LOG_INFO     = 3,
    SE_LOG_DEBUG    = 4,
    SE_LOG_TRACE    = 5,
    SE_LOG_PROGRESS = 6
} se_log_type;

#define _SE_LOG_MSG_MAX_LEN 2048
#define _SE_LOG_INFO_MAX_LEN 1024

/** please don't use it directly */
extern int _se_private_verb_level_variable;

inline void se_set_verb_level(int v)
{
    _se_private_verb_level_variable = v;
}

inline int se_get_verb_level()
{
    return _se_private_verb_level_variable;
}

inline int se_log_verb_test( int verb, const char* format, ... )
{
    (void)format; // avoid unused parameter warning
    return ( verb <= se_get_verb_level() );
}

/** please don't use it directly */
extern int _se_private_log_id_variable;

inline void se_set_logid(int lid)
{
    _se_private_log_id_variable = lid;
}

inline int se_get_logid()
{
    return _se_private_log_id_variable;
}

/**
 * For the moment there isn't really a need to have more than one
 * additional logger. If more are needed, this function can be
 * implemented in terms of "add_logger" 
 */
void se_log_set_additional_logger(void (*logger)( se_log_type type, int logid, 
                                                   const char* file, int line,
                                                   const char* message, size_t message_len));

/** Check if logging is enabled for the type and file/line */
int se_log_enabled( se_log_type type, const char* file, int line );

/** Called to collect extra log info */
void se_collect_log_info( se_log_type type, const char* file, int line );

/** Called to log something */
void se_log_handler( const char* format, ... );
void se_logv_handler( int verb, const char* format, ... );
void se_progress_handler( double progress, const char* format, ... );


/** Called in case of error */
void se_error_handler( const char* format, ... );

/** Set the output for logging. The default is stderr. */
void se_set_log_file( const char* file_name );

/** Returns the current logging stream. */
FILE* se_get_log_stream();

/** Set the date format for logging */
void se_set_log_date_format( const char* date_format );

int se_get_log_filed();

void _se_log_lock();
void _se_log_unlock();

/**
 * Useful when implementing a different logger.
 *
 * The return, infostr and infostr_size have the same meaning as for snprintf;
 */
int se_default_log_info_format( char* infostr, size_t infostr_size,
                                 se_log_type type, int logid, 
                                 const char* file, int line );

/**
 * This is the default logger: it will format the log info string,
 * append the message and send it to the current log stream.
 *
 * message_len should be the lenght of the message string, not
 * including the trailing '\0'.
 */
void se_log_default_logger( se_log_type type, int logid, 
                             const char* file, int line,
                             const char* message, size_t message_len);


/** _LOG(t,h,f...): call the log handler:
 * t  - type
 * h  - log handler - either log_handler or error_handler
 * f  - message - printf style with parenthesis : ("i=%d\n", i)
 */
#ifdef SE_LOG_SRC
#define _LOG(t, h, f)                                          \
    do {                                                       \
        if( t == SE_LOG_ERROR ||                              \
            se_log_enabled(t, __FILE__, __LINE__) ) {         \
            _se_log_lock();                                   \
            se_collect_log_info(t, __FILE__, __LINE__);       \
            h f;                                               \
            _se_log_unlock();                                 \
        }                                                      \
    } while(0)

#define _LOGV(t, h, f)                                         \
    do {                                                       \
        if( t == SE_LOG_ERROR ||                              \
            ( se_log_verb_test f &&                           \
              se_log_enabled(t, __FILE__, __LINE__) ) ) {     \
            _se_log_lock();                                   \
            se_collect_log_info(t, __FILE__, __LINE__);       \
            h f;                                               \
            _se_log_unlock();                                 \
        }                                                      \
    } while(0)

#else /* no SE_LOG_SRC */

#define _LOG(t, h, f)                               \
    do {                                            \
        _se_log_lock();                            \
        se_collect_log_info(t, NULL, -1);          \
        h f;                                        \
        _se_log_unlock();                          \
    } while(0)

#define _LOGV(t, h, f)                              \
    do {                                            \
        if( se_log_verb_test f ) {                 \
            _se_log_lock();                        \
            se_collect_log_info(t, NULL, -1);      \
            h f;                                    \
            _se_log_unlock();                      \
        }                                           \
    } while(0)

#endif /* defined SE_LOG_SRC */


#ifdef SE_ENABLE_ERROR_LOG
#  define ERROR(f) _LOG(SE_LOG_ERROR, se_error_handler, f)
#else
#  define ERROR(f) NOOP
#endif

#ifdef SE_ENABLE_WARN_LOG
#  define WARN(f) _LOG(SE_LOG_WARN, se_log_handler, f)
#else
#  define WARN(f) NOOP
#endif

#ifdef SE_ENABLE_NOTICE_LOG
#  define NOTICE(f) _LOG(SE_LOG_NOTICE, se_log_handler, f)
#else
#  define NOTICE(f) NOOP
#endif

#ifdef SE_ENABLE_INFO_LOG
#  define INFO(f) _LOG(SE_LOG_INFO, se_log_handler, f)
#  define INFOV(f) _LOGV(SE_LOG_INFO, se_logv_handler, f)
#else
#  define INFO(f) NOOP
#  define INFOV(f) NOOP
#endif

#ifdef SE_ENABLE_DEBUG_LOG
#  define DEBUG(f) _LOG(SE_LOG_DEBUG, se_log_handler, f)
#  define DEBUGV(f) _LOGV(SE_LOG_DEBUG, se_logv_handler, f)
#else
#  define DEBUG(f) NOOP
#  define DEBUGV(f) NOOP
#endif

#ifdef SE_ENABLE_TRACE_LOG
#  define TRACE(f) _LOG(SE_LOG_TRACE, se_log_handler, f)
#  define TRACEV(f) _LOGV(SE_LOG_TRACE, se_logv_handler, f)
#else
#  define TRACE(f) NOOP
#  define TRACEV(f) NOOP
#endif

#ifdef SE_ENABLE_PROGRESS_LOG
#  define PROGRESS(f) _LOG(SE_LOG_PROGRESS, se_progress_handler, f)
#else
#  define PROGRESS(f) NOOP
#endif





#endif