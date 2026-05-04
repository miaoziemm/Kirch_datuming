#ifndef SE_TIMER_H
#define SE_TIMER_H
#include <sys/time.h>
#include "se_type.h"
#include "se_alloc.h"
#include "se_log.h"
#include "se_assert.h"
#include "se_return_code.h"
#include "se_hash.h"
#include "se_util.h"
#include <time.h>    // 用于timespec
#include <stddef.h>  // 用于size_t
#include <stdlib.h>
#include <stdio.h>
#include <stdarg.h>  //可变参数
#include <string.h>
#include <errno.h>  //报错
#include <sys/time.h>       //时间头文件
#include <sys/resource.h>      //资源操作 读写
#include <unistd.h>
#include <limits.h> //用于检测整型数据类型的表达值范围

/**
 * Returns to ut the user time in seconds and to st the system time in
 * secons.
 */
void get_time(real* ut, real* st);

/**
 * Returns the elapsed time in seconds since the beginning of
 * execution of the current process.
 * tarray[0] contains the user time in seconds.
 * tarray[1] contains the system time in seconds.
 *
 * This function is similar to etime from Fortran.
 */
real etime(real tarray[2]);

/**
 * Returns the current time in milliseconds. This is the difference,
 * measured in seconds, between the current time and midnight,
 * January 1, 1970 UTC.
 * The granulaztry of the value depends on the gettimeofday function.
 */
double current_time_sec(void);


#define ZTR_INFINITE_TIMEOUT 0

/**
 * Converts 'millis' timeout from milliseconds to absolute time
 **/
struct timespec millis2abstime(uint64_t millis);


/**
 * A structure that can be used to measure time spent on
 * different operations.
 * Typical use:<pre>
 *  timer* t1 = create_timer("Name1");
 *  timer* t2 = create_timer("Name2");
 *  ...
 *  loop
 *      ...
 *      start_timer(t1);
 *      some_operation1();
 *      stop_timer(t1);
 *      ...
 *      start_timer(t2);
 *      some_operation2();
 *      stop_timer(t2);
 *      ...
 *      report_all_timers();
 *      ...
 *  end loop
 *  report_all_timers();
 *  free_timer(t1);
 *  free_timer(t2);
 * </pre>
 */
typedef struct timer_s {
    char* name;  /* The name of the timer */     //timer的名字

    real utime; /* user time */      //用户时间
    real stime; /* system time */     //系统时间
    real etime; /* elapsed time */      //消逝的时间

    /* instantaneous absolute time */
    struct timeval a_tv;      //瞬时绝对时间  //秒+毫秒

    /* Accumulated values */
    real total_utime;
    real total_stime;
    real total_etime;
    real total_atime;

    int64_t n_calls;
} timer;


/**
 * Creates a named timer. 'name' can be NULL; if not null a copy of
 * the name will be kept (not the pointer).
 */
timer* create_timer(const char* name, ...);

/**
 * Releases resources used by a timer.
 */
void destroy_timer(timer* t);

/**
 * Reset a timer to zero.
 */
void reset_timer(timer* t);

/**
 * Start the watch.
 */
void start_timer(timer* t);

/**
 * Start the watch.  //stop the watch
 */
void stop_timer(timer* t);

/**
 * Log a report (using INFO log level) for an array of timers.
 */
void report_timers(timer** timers, size_t n);

/**
 * Log a report (using INFO log level) for all registered timers.
 */
void report_all_timers(void);

real read_and_stop_timer(timer* t);

real abs_read_timer(timer* t);

void accumulate_timer(timer* dst, timer* src);

/**
 * Save an array of timers to text form, so it can be read in later.
 */
char* timers_to_string(timer** t, size_t n);

timer** timers_from_string(const char* str, timer** t, size_t* n);


#endif