#include "../include/se_timer.h"


/** Functions to keep track of active timers */
static void _add_timer(timer* t);
static void _remove_timer(timer* t);


/**
 * Returns to ut the user time in seconds and to st the system time in
 * secons.
 */
void get_time(real* ut, real* st)
{
    struct rusage rbuff;

    if( getrusage (RUSAGE_SELF, &rbuff) != 0 ) {  //RUSAGE_SELF 获得当前进程的资源使用信息  0返回成功 -1返回失败
        ERROR(("getrusage:%s", strerror(errno)));
    }

    *ut = (real)( (double) rbuff.ru_utime.tv_sec +
                  (double) rbuff.ru_utime.tv_usec / 1000000.0 );  //tv_usec 微秒 所以要除以1000000
    *st = (real)( (double) rbuff.ru_stime.tv_sec +
                  (double) rbuff.ru_stime.tv_usec / 1000000.0 );
}


/**
 * Returns the elapsed time in seconds since the beginning of
 * execution of the current process.
 * tarray[0] contains the user time in seconds.
 * tarray[1] contains the system time in seconds.
 *
 * This function is similar to etime from Fortran.
 */
real etime(real tarray[2])  //用户时间+系统时间
{
    get_time(&tarray[0], &tarray[1]);

    return ( tarray[0] + tarray[1] );
}


/**
 * Returns the current time in milliseconds. This is the difference,
 * measured in seconds, between the current time and midnight,
 * January 1, 1970 UTC.
 * The granulaztry of the value depends on the gettimeofday function.
 */
double current_time_sec(void)  //in seconds not in milliseconds
{
    struct timeval time;

    int status = gettimeofday(&time, NULL);
    if(status == -1) {
        ERROR(("gettimeofday:%s", strerror(errno)));
    }

    return ( (double) time.tv_sec  +
             (double) time.tv_usec / 1000000.0 );
}


/**
 * Converts 'millis' timeout from milliseconds to absolute time
 **/
struct timespec millis2abstime( uint64_t millis )    //nsec纳秒
{
    struct timespec abstime;

    if( millis == ZTR_INFINITE_TIMEOUT ) {
        abstime.tv_sec  = LONG_MAX;
        abstime.tv_nsec = 999999999;
    } else {
        struct timeval now;
        if( gettimeofday(&now, NULL) ) {
            int err = errno;
            ERROR(( "gettimeofday() failed: errno=%d - %s",
                    err, strerror(err) ));
        }

        abstime.tv_sec  = now.tv_sec + millis / 1000;
        abstime.tv_nsec = 1000 * (now.tv_usec + (millis % 1000) * 1000);
        if( abstime.tv_nsec >= 1000000000 ) {
            abstime.tv_sec++;
            abstime.tv_nsec -= 1000000000;   //检测是否超过一个单位
        }  
    }

    return abstime;
}



/**********************************************************************/
/*************************  TIMERS ************************************/


/**
 * Creates a named timer. 'name' can be NULL;
 */
timer* create_timer(const char* name, ...)         //变参数函数
{
    va_list args;
    timer* t;

    t = (timer*)malloc( sizeof(timer) );  //分配内存
    memset(t, 0, sizeof(timer));  //初始化为0

    if( name != NULL ) {
        size_t l = strlen(name);
        if(l < 256) l = 256;
        t->name = alloc1char(l+1);
        va_start( args, name );
        vsnprintf( t->name, l, name, args );
        va_end( args );
        t->name[l] = 0;
    } else {
        t->name = NULL;
    }

    reset_timer(t);  //初始化为0

    _add_timer(t);         //加入到timer_list

    return t;
}

/**
 * Releases resources used by a timer.
 */
void destroy_timer(timer* t)
{
    _remove_timer(t);

    if( t->name != NULL )
        free1char(t->name);

    free(t);
}

/**
 * Reset a timer to zero.
 */
void reset_timer(timer* t)
{
    t->utime = 0.0;
    t->stime = 0.0;
    t->etime = 0.0;

    t->a_tv.tv_sec  = 0;
    t->a_tv.tv_usec  = 0;

    t->total_utime = 0.0;
    t->total_stime = 0.0;
    t->total_etime = 0.0;
    t->total_atime = 0.0;

    t->n_calls = 0;
}


/**
 * Start the watch.
 */
void start_timer(timer* t)
{
    get_time( &(t->utime), &(t->stime) );

    if( gettimeofday(&(t->a_tv), NULL) == -1 ) {
        ERROR(("gettimeofday:%s", strerror(errno)));
    }

    t->n_calls += 1;
}

/**
 * Start the watch.
 */
void stop_timer(timer* t)
{
    real ut, st, at;
    struct timeval tv;

    get_time( &ut, &st);

    if( gettimeofday(&tv, NULL) == -1 ) {
        ERROR(("gettimeofday:%s", strerror(errno)));
    }

    ut -= t->utime;
    st -= t->stime;

    at = (real) ( (double)( tv.tv_sec  - t->a_tv.tv_sec  ) +
                  (double)( tv.tv_usec - t->a_tv.tv_usec ) /
                  (double)1000000.0 );

    if( ut < 0.0 ) ut = 0.0;
    if( st < 0.0 ) st = 0.0;
    if( at < 0.0 ) at = 0.0;

    t->utime = ut;
    t->stime = st;
    t->etime = ut+st;

    t->a_tv.tv_sec  = tv.tv_sec;
    t->a_tv.tv_usec = tv.tv_usec;

    t->total_utime += ut;
    t->total_stime += st;
    t->total_etime += st+ut;
    t->total_atime += at;
}

void accumulate_timer(timer* dst, timer* src)  //计算累计时间
{
    dst->n_calls     += src->n_calls;
    dst->total_utime += src->total_utime;
    dst->total_stime += src->total_stime;
    dst->total_atime += src->total_atime;
    dst->total_etime = dst->total_utime + dst->total_stime;
}

/**
 * Log a report (using INFO log level) for an array of timers.
 */
void report_timers(timer** t, size_t n)
{
    size_t i;
    char format_str[101];
    char format_str_title[101];
    const char* title_name = "Time statistics";
    size_t max_name_len;
    real total;

    max_name_len = strlen(title_name);
    total = 0.0;
    for(i = 0; i < n; i++) {
        if( t[i] == NULL ) continue; /* ignore empty slots */

        if( t[i]->name != NULL ) {
            size_t l = strlen(t[i]->name);
            if( max_name_len < l ) max_name_len = l;
        }

        total += t[i]->total_atime;
    }

    snprintf(format_str_title, sizeof(format_str_title) - 1,
             "%%%lus %%11s%%8s%%11s%%11s%%11s%%11s%%10s",
             (unsigned long)max_name_len);
    snprintf(format_str, sizeof(format_str) - 1,
             "%%%lus:%%11.2f%%8.2f%%11.2f%%11.2f%%11.2f%%9.4f%%10d",
             (unsigned long)max_name_len);

    INFO((format_str_title, title_name,
          "wall", "prct[%]", "user", "system",
          "elps", "MC-Boost", "count"));

    for(i = 0; i < n; i++) {
        const char* name;
        double prct;

        if( t[i] == NULL ) continue; /* ignore empty slots */

        name = t[i]->name != NULL ? t[i]->name : "";
        if( total > 0.0 )
            prct = 100.0 * t[i]->total_atime / total;
        else
            prct = 0.0;

        INFO((format_str, name,
              t[i]->total_atime, prct,
              t[i]->total_utime, t[i]->total_stime, t[i]->total_etime,
              t[i]->total_atime>0?t[i]->total_etime/t[i]->total_atime:0.0,
              t[i]->n_calls));
    }
}


/** a list of active timers, used mostily to print a nice report */
/** TO_DO: the access to this functions must be synchronized */  //同步化

static timer** _timers_list = NULL;       //time list n time time size
static size_t _n_timers = 0;
static size_t _timer_list_size = 0;

#define TIMERS_LIST_GROWS_BY 64

/**
 * Log a report (using INFO log level) for all registered timers.
 */
void report_all_timers()
{
    report_timers(_timers_list, _n_timers);
}


static void _ensure_space(size_t n)
{
    if( _timer_list_size < n ) {
        size_t nsize = n + TIMERS_LIST_GROWS_BY;                //Timers_list_grows_by 64
        _timers_list = (timer**)realloc(_timers_list,
                                                nsize*sizeof(timer*));  //realloc  重新分配
        _timer_list_size = nsize;
    }
}


static void _add_timer(timer* t)
{
#ifdef SE_USE_OMP
#pragma omp critical 
#endif
    {
        size_t n = _n_timers + 1;

        _ensure_space(n);

        _timers_list[n-1] = t;

        _n_timers = n;
    }
}

static void _remove_timer(timer* t)
{
    int found = 0;
#ifdef SE_USE_OMP
#pragma omp critical 
#endif
    {
        size_t i, j;
        size_t n = _n_timers;
        for(i = 0; i < n; i++) {
            if( _timers_list[i] == t ) {
                for(j = i+1; j < n; j++) {
                    _timers_list[j-1] = _timers_list[j];
                }
                if( --_n_timers == 0 ) {       //因为只剩下一个 所以直接移除  很right
                    free(_timers_list);
                    _timers_list = NULL;
                    _timer_list_size = 0;
                }
                found  = 1;
                break;
            }
        }
    }

    if(found) return;
    /* If we got here, we didn't find the timer in the list. */
    WARN(("Couldn't find timer at adress %x in timers list", t));
    if( t->name != NULL ) {
        WARN(("\tThe name of the timer is: %s", t->name));
    } else {
        WARN(("\tThe timer has no name"));
    }
}

/* this stops the timer, also returns the time passed so far (in sec) */
real read_and_stop_timer(timer* t)
{
    real ut, st, at;
    struct timeval tv;

    get_time( &ut, &st);

    if( gettimeofday(&tv, NULL) == -1 ) {
        ERROR(("gettimeofday:%s", strerror(errno)));
    }

    ut -= t->utime;
    st -= t->stime;

    at = (real) ( ( tv.tv_sec  - t->a_tv.tv_sec  ) +
                  (double)( tv.tv_usec - t->a_tv.tv_usec ) /
                  (double)1000000.0 );

    if( ut < 0.0 ) ut = 0.0;
    if( st < 0.0 ) st = 0.0;
    if( at < 0.0 ) at = 0.0;

    t->utime = ut;
    t->stime = st;
    t->etime = ut+st;

    t->a_tv.tv_sec  = tv.tv_sec;
    t->a_tv.tv_usec = tv.tv_usec;

    t->total_utime += ut;
    t->total_stime += st;
    t->total_etime += st+ut;
    t->total_atime += at;

    return at;
}

/* returns the absolute time spent on timer ... */
real abs_read_timer(timer* t)
{
    real at;
    struct timeval tv;

    if( gettimeofday(&tv, NULL) == -1 ) {
        ERROR(("gettimeofday:%s", strerror(errno)));
    }

    at = (real) ( ( tv.tv_sec  - t->a_tv.tv_sec  ) +
                  (real)( tv.tv_usec - t->a_tv.tv_usec ) / (real)1000000.0 );

    if( at < 0.0 ) at = 0.0;

    return at;
}

char* timers_to_string(timer** t, size_t n)
{
    char* str = NULL;
    int len = 0;
    int i;
    ASSERT(t != NULL);
    for(i = 0; i < (int)n; ++i) {
        se_hash ht = create_hash();
        ASSERT(t[i] != NULL);
        ht_put(ht, "n", int64_to_str(t[i]->n_calls));
        ht_put(ht, "u", double_to_str(t[i]->total_utime));
        ht_put(ht, "s", double_to_str(t[i]->total_stime));
        ht_put(ht, "a", double_to_str(t[i]->total_atime));
        ht_put(ht, "o", int64_to_str(i));

        ht_to_str(ht, t[i]->name, &str, &len);
        destroy_hash_and_keyval( ht, 0, 1 );
    }

    return str;
}

typedef struct ordered_timer_s {
    timer* t;
    int o;
}ordered_timer_s;

static int cmp_ordered_timers(const void *p1, const void *p2)
{
    ordered_timer_s* ot1 = *((ordered_timer_s**)p1);
    ordered_timer_s* ot2 = *((ordered_timer_s**)p2);
    if(ot1->o > ot2->o) return 1;
    if(ot1->o < ot2->o) return -1;
    return 0;
}

timer** timers_from_string(const char* str, timer** t, size_t* n)
{
    se_hash htimers;
    se_hash ht;
    htiter hti;
    int have_new_timers = 0;

    ht = ht_from_str( str );
    if(ht == NULL) return t;

    htimers = create_hash();

    ASSERT(n != NULL);
    if(t != NULL) {
        int i;
        for(i = 0; i < (int)*n; ++i) {
            ordered_timer_s* ot = (ordered_timer_s*)malloc(sizeof(*ot));
            ASSERT(t[i] != NULL);
            ot->t = t[i];
            ot->o = i;
            ht_put(htimers, t[i]->name, ot);
        }
    }

    hti = create_htiter(ht);
    while(hti_hasnext(hti)) {
        const char* name;
        se_hash* htm;
        ordered_timer_s* ot;
        char* tmp;
        union {
            se_hash** h;
            void** v;
        }dosomecasting;
        dosomecasting.h = &htm;
        hti_next( hti, &name, dosomecasting.v);

        ot = (ordered_timer_s*)ht_get(htimers, name);
        if(ot == NULL) {
            ot = (ordered_timer_s*)malloc(sizeof(*ot));
            ot->t = create_timer("%s", name);
            ht_put(htimers, ot->t->name, ot);
            have_new_timers = 1;
            tmp = (char *)ht_get(htm, "o");
            if(tmp) ot->o = atol(tmp);
        }
        tmp = (char *)ht_get(htm, "n");
        if(tmp) ot->t->n_calls += atol(tmp);
        tmp = (char *)ht_get(htm, "u");
        if(tmp) ot->t->total_utime += (real)atof(tmp);
        tmp = (char *)ht_get(htm, "s");
        if(tmp) ot->t->total_stime += (real)atof(tmp);
        tmp = (char *)ht_get(htm, "a");
        if(tmp) ot->t->total_atime += (real)atof(tmp);

        ot->t->total_etime = ot->t->total_utime + ot->t->total_stime;

        destroy_hash_and_keyval(htm, 1, 1);
    }
    destroy_htiter( hti );

    if(have_new_timers) {
        int i;
        ordered_timer_s** list;
        *n = ht_size(htimers);
        list = (ordered_timer_s**)malloc(*n*sizeof(*list));
        hti = create_htiter(htimers);
        for(i = 0; hti_hasnext(hti); ++i) {
            const char* name;
            ordered_timer_s* ot;
            union {
                ordered_timer_s** h;
                void** v;
            }dosomecasting1;
            dosomecasting1.h = &ot;
            hti_next( hti, &name, dosomecasting1.v);
            list[i] = ot;
        }
        destroy_htiter( hti );

        qsort(list, *n, sizeof(ordered_timer_s*), cmp_ordered_timers);

        t = (timer**)realloc(t, *n*sizeof(*t));
        for(i = 0; i < (int)*n; ++i) {
            t[i] = list[i]->t;
        }
        free(list);
    }

    destroy_hash_and_keyval(ht, 1, 0);
    destroy_hash_and_keyval(htimers, 0, 1);

    return t;
}
