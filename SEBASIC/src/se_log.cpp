#include "../include/se_log.h"
#include <iostream>
#include <iomanip>
#include <chrono>
#include <ctime>
#include <cstdarg>
#include <cstdio>

void LOG_CONSOLE(LOGType type, const char *format, ...)
{
    auto now = std::chrono::system_clock::now();
    std::time_t now_time = std::chrono::system_clock::to_time_t(now);
    char time_buffer[100];
    std::strftime(time_buffer, sizeof(time_buffer), "%Y-%m-%d %H:%M:%S", std::localtime(&now_time));

    FILE *out = stdout;
    const char *type_str = "[INFO] ";
    switch (type)
    {
    case SE_INFO:
        type_str = "[INFO] ";
        out = stdout;
        break;
    case SE_WARNING:
        type_str = "[WARNING] ";
        out = stdout;
        break;
    case SE_ERROR:
        type_str = "[ERROR] ";
        out = stderr;
        break;
    case SE_FATAL:
        type_str = "[FATAL] ";
        out = stderr;
        break;
    }

    fprintf(out, "[%s] %s", time_buffer, type_str);

    va_list args;
    va_start(args, format);
    vfprintf(out, format, args);
    va_end(args);

    fprintf(out, "\n");

    fflush(out);
}


/** Names for the levels of logs.
    The values for se_log_type enum should match exactly this array. */
static const char* _log_names[] = {
    "ERROR",
    "WARN",
    "NOTICE",
    "INFO",
    "DEBUG",
    "TRACE",
    "PROGRESS"
};

#ifndef ZTR_DEFAULT_DATE_FORMAT
#define ZTR_DEFAULT_DATE_FORMAT "%F %H:%M:%S"
#endif

/** The log file descriptor. */
static FILE* _log_file = NULL;

/** Stores the verbosity level. */
int _se_private_verb_level_variable = 0;

/** Stores the log id. */
int _se_private_log_id_variable = 0;

FILE* se_get_log_stream()
{
    if(_log_file) return _log_file;
    return stderr;
}

int se_get_log_filed()
{
    return fileno(se_get_log_stream());
}


void se_set_log_file( const char* file_name )
{
    FILE* f;

    _se_log_lock();
    if(_log_file ) {
        fflush(_log_file);
    }
    _se_log_unlock();

    if(file_name) {
        f = fopen( file_name, "a" );
        if( f == NULL ) {
            int err = errno;
            ERROR(("Cannot open log file %s: %d - %s",
                   file_name, err, strerror(err)));
        }
        setvbuf(f, (char *)NULL, _IOLBF, 0);
    } else {
        f = NULL; /*this will have the effect of using stderr, but we don't close it*/
    }

    _se_log_lock();
    if(_log_file ) { /* we open it, we close it */
        fclose(_log_file);
    }
    _log_file = f;
    _se_log_unlock();
}

static void (*_se_additonal_logger)( se_log_type type, int logid, 
                                      const char* file, int line,
                                      const char* message, size_t message_len) = NULL;

void se_log_set_additional_logger(void (*logger)( se_log_type type, int logid, 
                                                  const char* file, int line,
                                                  const char* message, size_t message_len))
{
    _se_log_lock();
    _se_additonal_logger = logger;
    _se_log_unlock();
}

#ifndef NO_LOG_TIME
static char _loginfo_datestring[256] = ZTR_DEFAULT_DATE_FORMAT;
#endif

void se_set_log_date_format( const char* date_format )
{
#ifndef NO_LOG_TIME
    strncpy(_loginfo_datestring, date_format, sizeof(_loginfo_datestring) - 1);
    _loginfo_datestring[sizeof(_loginfo_datestring) - 1] = '\0';  // 确保以NULL结尾
#endif
}

/**
 * The loginfo is not available in the log handler; the log handler
 * gets only the original format string and the variable list of
 * arguments.
 *
 * This is the reason we collect the extra information in global
 * variable, before calling the log handler.
 *
 * The access to these variables is controlled with a lock (see the se_log.h file).
 */
static se_log_type _loginfo_type;
#ifndef NO_LOG_SRC
static const char* _loginfo_file;
static int _loginfo_line;
#endif

/**
 * Called to save extra info for logging. This function returns fast -
 * it copies at most two integers and a pointer to a string, but
 * usually it copies only an integer - the log type.
 */
void se_collect_log_info( se_log_type type, const char* file, int line )
{
    _loginfo_type = type;
#ifndef NO_LOG_SRC
    _loginfo_file = file;
    _loginfo_line = line;
#else
    (void)file;
    (void)line;
#endif
}


int se_default_log_info_format( char* infostr, size_t infostr_size,
                                 se_log_type type, int logid, 
                                 const char* file, int line )
{
    /* looks ugly because we are trying to do everything with just one
       printf call */
#ifndef NO_LOG_TIME
    char datestr[256];
#else
    const char* datestr = "";
#endif
    int ret;

#ifndef NO_LOG_TIME
    {
        time_t t;
        struct tm *tmp;
        t = time(NULL);
        tmp = localtime(&t);
        if(tmp == NULL) {
            strncpy(datestr, "unknown", sizeof(datestr));
        } else {
            strftime(datestr, sizeof(datestr), _loginfo_datestring, tmp);
            datestr[sizeof(datestr)-1] = 0;
        }
    }
#endif

#ifndef NO_LOG_SRC
    if(line != -1) {
#ifndef NO_LOG_TYPE
#ifndef NO_LOG_LOGID
        ret = snprintf(infostr, infostr_size,
                       "%d %s %s[%s:%d]:", logid, datestr,
                       _log_names[type], file, line);
#else
        ret = snprintf(infostr, infostr_size,
                       "%s %s[%s:%d]:", datestr, _log_names[type], file, line);
#endif
#else
#ifndef NO_LOG_LOGID
        ret = snprintf(infostr, infostr_size,
                       "%d %s [%s:%d]:", logid, datestr, file, line);
#else
        ret = snprintf(infostr, infostr_size,
                       "%s [%s:%d]:", datestr, file, line);
#endif
#endif
    } else
#endif
        {
            (void)file;
            (void)line;
#ifndef NO_LOG_TYPE
#ifndef NO_LOG_LOGID
            ret = snprintf(infostr, infostr_size,
                           "%d %s %s:", logid, datestr, _log_names[type]);
#else
            ret = snprintf(infostr, infostr_size,
                           "%s %s:", datestr, _log_names[type]);
#endif
#else
            ret = snprintf(infostr, infostr_size,
                           "%s %d:", datestr, logid);
#endif
        }

    return ret;
}


void se_log_default_logger( se_log_type type, int logid, 
                             const char* file, int line,
                             const char* message, size_t message_len)
{
     /* Trying not to allocate any new memory in this call.  The call
        to this function is synchronized anyway, so we can use a
        static variable */
    static char infostr[_SE_LOG_INFO_MAX_LEN];
    int infosize = se_default_log_info_format( infostr, sizeof(infostr), type, logid, file, line );
    FILE* log_file = se_get_log_stream();

    /* At this point everything is formated, so lets avoid calling yet another printf */

    if(infosize > 0) {
        /* just making sure the string ends; don't really care if is truncated */
        if(infosize >= (int)sizeof(infostr)) {
            infosize = sizeof(infostr)-1;
            infostr[infosize] = '\0';
        }
        fwrite(infostr, infosize, 1, log_file);
    }

    if(message && message_len > 0) {
        fwrite(message, message_len, 1, log_file);
    }

    fputc('\n', log_file );
}


/* Trying not to allocate any new memory in the next functions calls. The call
   to these functions is synchronized anyway, so we can use a
   static variable */
static char _log_msg[_SE_LOG_MSG_MAX_LEN];

/** Called to log something... */
void se_log_handler( const char* format, ... )
{
    int msg_len;
    va_list args;

    if( format != NULL ) {
        va_start( args, format );
        msg_len = vsnprintf( _log_msg, sizeof(_log_msg), format, args );
        va_end( args );
        if(msg_len > 0) {
            if(msg_len >= (int)sizeof(_log_msg)) {
                msg_len = sizeof(_log_msg) - 1;
                _log_msg[msg_len] = '\0';
            }
        } else {
            msg_len = 0;
            _log_msg[0] = '\0';
        }
    } else {
        msg_len = 0;
        _log_msg[0] = '\0';
    }

#ifndef NO_LOG_SRC
    se_log_default_logger( _loginfo_type, se_get_logid(), 
                            _loginfo_file, _loginfo_line,
                            _log_msg, msg_len);
    if(_se_additonal_logger)
        _se_additonal_logger( _loginfo_type, se_get_logid(), 
                               _loginfo_file, _loginfo_line,
                               _log_msg, msg_len);
#else
    se_log_default_logger( _loginfo_type, se_get_logid(), 
                            NULL, -1,
                            _log_msg, msg_len);
    if(_se_additonal_logger)
        _se_additonal_logger( _loginfo_type, se_get_logid(), 
                               NULL, -1,
                               _log_msg, msg_len);
#endif
}

/** Called to log something for a verbosity level ... */
void se_logv_handler( int verb, const char* format, ... )
{
    int msg_len;
    va_list args;

    if( verb > se_get_verb_level() ) return;

    if( format != NULL ) {
        va_start( args, format );
        msg_len = vsnprintf( _log_msg, sizeof(_log_msg), format, args );
        va_end( args );
        if(msg_len > 0) {
            if(msg_len >= (int)sizeof(_log_msg)) {
                msg_len = sizeof(_log_msg) - 1;
                _log_msg[msg_len] = '\0';
            }
        } else {
            msg_len = 0;
            _log_msg[0] = '\0';
        }
    } else {
        msg_len = 0;
        _log_msg[0] = '\0';
    }

#ifndef NO_LOG_SRC
    se_log_default_logger( _loginfo_type, se_get_logid(), 
                            _loginfo_file, _loginfo_line,
                            _log_msg, msg_len);
    if(_se_additonal_logger)
        _se_additonal_logger( _loginfo_type, se_get_logid(), 
                               _loginfo_file, _loginfo_line,
                               _log_msg, msg_len);
#else
    se_log_default_logger( _loginfo_type, se_get_logid(), 
                            NULL, -1,
                            _log_msg, msg_len);
    if(_se_additonal_logger)
        _se_additonal_logger( _loginfo_type, se_get_logid(), 
                               NULL, -1,
                               _log_msg, msg_len);
#endif
}

void se_progress_handler( double progress, const char* format, ... )
{
    int msg_len;
    int progress_len;
    va_list args;

    if( se_get_verb_level() < 0 ) return;

    progress_len = snprintf(_log_msg, sizeof(_log_msg), "%.2f%% - ", progress);
    if(progress_len < 0) {
        progress_len = 0;
    }

    if( format != NULL ) {
        va_start( args, format );
        msg_len = vsnprintf( _log_msg+progress_len, sizeof(_log_msg)-progress_len, format, args );
        va_end( args );
        if(msg_len > 0) {
            if(msg_len >= (int)(sizeof(_log_msg)-progress_len)) {
                msg_len = sizeof(_log_msg)-progress_len - 1;
                _log_msg[progress_len+msg_len] = '\0';
            }
        } else {
            msg_len = 0;
            _log_msg[progress_len] = '\0';
        }
    } else {
        msg_len = 0;
        _log_msg[progress_len] = '\0';
    }
    msg_len += progress_len;

#ifndef NO_LOG_SRC
    se_log_default_logger( _loginfo_type, se_get_logid(), 
                            _loginfo_file, _loginfo_line,
                            _log_msg, msg_len);
    if(_se_additonal_logger)
        _se_additonal_logger( _loginfo_type, se_get_logid(), 
                               _loginfo_file, _loginfo_line,
                               _log_msg, msg_len);
#else
    se_log_default_logger( _loginfo_type, se_get_logid(), 
                            NULL, -1,
                            _log_msg, msg_len);
    if(_se_additonal_logger)
        _se_additonal_logger( _loginfo_type, se_get_logid(), 
                               NULL, -1,
                               _log_msg, msg_len);
#endif
}


/********************************************************************/
static int _se_error_handler_called_already = 0;
/** Called in case of error ... */
void se_error_handler( const char* format, ... )
{
    int msg_len;
    va_list args;

    if(_se_error_handler_called_already) return;
    _se_error_handler_called_already = 1;

    if( format != NULL ) {
        va_start( args, format );
        msg_len = vsnprintf( _log_msg, sizeof(_log_msg), format, args );
        va_end( args );
        if(msg_len > 0) {
            if(msg_len >= (int)sizeof(_log_msg)) {
                msg_len = sizeof(_log_msg) - 1;
                _log_msg[msg_len] = '\0';
            }
        } else {
            msg_len = 0;
            _log_msg[0] = '\0';
        }
    } else {
        msg_len = 0;
        _log_msg[0] = '\0';
    }

#ifndef NO_LOG_SRC
    se_log_default_logger( _loginfo_type, se_get_logid(), 
                            _loginfo_file, _loginfo_line,
                            _log_msg, msg_len);
    if(_se_additonal_logger)
        _se_additonal_logger( _loginfo_type, se_get_logid(), 
                               _loginfo_file, _loginfo_line,
                               _log_msg, msg_len);
#else
    se_log_default_logger( _loginfo_type, se_get_logid(), 
                            NULL, -1,
                            _log_msg, msg_len);
    if(_se_additonal_logger)
        _se_additonal_logger( _loginfo_type, se_get_logid(), 
                               NULL, -1,
                               _log_msg, msg_len);
#endif

#ifdef ZTR_DEBUG
    se_print_stack_trace();
#endif

    fflush(stdout);
    fflush(stderr);

    if( _log_file ) {
        fclose(_log_file);
        _log_file = NULL;
    }

    _se_log_unlock(); /*we don't expect the next function to return,
                         so we must force unlocking the log lock, just
                         in case there are some calls to log functions
                         inside the exit function */
    exit( EXIT_FAILURE );
}



/********************************************************************/

/** Check if logging is enabled for the type and file/line */
int se_log_enabled( se_log_type type, const char* file, int line )
{
    char* disabled = getenv(_log_names[type]);
    (void)file;
    (void)line;

    if(disabled == NULL) return 1;

    /*return ( ! se_pattern_match(disabled, file) );*/
    return 1;
}

/****************************************************************/
/* synchronzied access to the output file */
#ifdef SE_USE_OMP
static omp_lock_t _log_omp_lock;
static volatile int _log_omp_lock_initialized = 0;
 
static void _se_ensure_log_omp_lock_initialized()
{
    if( _log_omp_lock_initialized ) return;

#pragma omp critical (_log_lock_init_critical_name)
    {
        if( ! _log_omp_lock_initialized ) {
            omp_init_lock(&_log_omp_lock);
            _log_omp_lock_initialized = 1;
        }
    }
}

static void _se_omplock()
{
    _se_ensure_log_omp_lock_initialized();
    omp_set_lock(&_log_omp_lock);
}

static void _se_ompunlock()
{
    ASSERT(_log_omp_lock_initialized);
    omp_unset_lock(&_log_omp_lock);
}

#else
static void _se_omplock()
{
}

static void _se_ompunlock()
{
}
#endif

#ifdef ZTR_USE_PTHREADS
#include <pthread.h>
static pthread_mutex_t log_mutex = PTHREAD_MUTEX_INITIALIZER;

void _se_log_lock()
{
    pthread_mutex_lock(&log_mutex);
    _se_omplock();
}

void _se_log_unlock()
{
    _se_ompunlock();
    pthread_mutex_unlock(&log_mutex);
}

#else /* ! ZTR_USE_PTHREADS */

void _se_log_lock()
{
    _se_omplock();
}

void _se_log_unlock()
{
    _se_ompunlock();
}

#endif /* ZTR_USE_PTHREADS */

