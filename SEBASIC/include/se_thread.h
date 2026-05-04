#ifndef SE_THREAD_H
#define SE_THREAD_H

#ifdef USE_PTHREADS
#include <pthread.h>
#endif


#ifdef USE_PTHREADS


/*
  mutex

  - simple wrappers for pthread_mutex_*() functions.
*/

typedef struct
{
    pthread_mutex_t mutex;
} se_mutex_t;

#define SE_MUTEX_INIT {.mutex = PTHREAD_MUTEX_INITIALIZER}

void se_mutex_init       (se_mutex_t* m);
void se_mutex_init_shared(se_mutex_t* m);
void se_mutex_destroy    (se_mutex_t* m);

void se_mutex_lock   (se_mutex_t* m);
int  se_mutex_trylock(se_mutex_t* m);
void se_mutex_unlock (se_mutex_t* m);


/*
  condition variable

  - simple wrappers for pthread_cond_*() functions.
*/

typedef struct
{
    pthread_cond_t cond;
} se_cond_t;

#define se_COND_INIT {.cond = PTHREAD_COND_INITIALIZER}

void se_cond_init       (se_cond_t* c);
void se_cond_init_shared(se_cond_t* c);
void se_cond_destroy    (se_cond_t* c);

void se_cond_signal   (se_cond_t* c);
void se_cond_broadcast(se_cond_t* c);

void se_cond_wait        (se_cond_t* c, se_mutex_t* m);
int  se_cond_wait_timeout(se_cond_t* c, se_mutex_t* m, uint64_t millis);


/*
  long lock

  - synchronization object:

  - a long lock have the same interface as a regular mutex (plus the lock_timeout
    function), but it provides better performance when the ammount of execution time
    spent between the lock() and unlock() calls is relatively big.
*/

typedef struct
{
    int is_locked;
    int wanted_count;
    se_mutex_t mutex;
    se_cond_t  cond;
} se_longlock_t;

#define se_LONGLOCK_INIT {.is_locked = 0, .wanted_count = 0, \
                           .mutex = se_MUTEX_INIT, .cond = se_COND_INIT}

void se_longlock_init       (se_longlock_t* ll);
void se_longlock_init_shared(se_longlock_t* ll);
void se_longlock_destroy    (se_longlock_t* ll);

void se_longlock_lock        (se_longlock_t* ll);
int  se_longlock_trylock     (se_longlock_t* ll);
int  se_longlock_lock_timeout(se_longlock_t* ll, uint64_t millis);

void se_longlock_unlock(se_longlock_t* ll);


/*
  event

  - complex synchronization object:

  - an event object can be used in notifying other threads that a particular
    event has occurred.
  - the state of the event can be zero (the reset state) or any other non-zero
    value.
  - threads can wait for particular state(s) of the event by calling any of the
    the event_wait_*() functions.
  - threads can change the state of the event, this will automaticaly notify
    threads waiting for that event state.
  - manual-reset events: the state of the event remains unchanged until it is
      explicitly reset by the event_reset(). seting a manual reset event will
      cause any number of waiting threads can be released.
  - auto-reset events: the state of the event remains unchanged until a single
      waiting thread is released. If more threads are waiting, then a single
      thread is selected arbitrary.
*/

typedef struct
{
    int state;
    int is_auto_reset;
    se_mutex_t mutex;
    se_cond_t  cond;
} se_event_t;

#define se_EVENT_AUTORESET_INIT   {.state = EVENT_STATE_RESET, .is_auto_reset = 1, \
                                    .mutex = se_MUTEX_INIT, .cond = se_COND_INIT}
#define se_EVENT_MANUALRESET_INIT {.state = EVENT_STATE_RESET, .is_auto_reset = 0, \
                                    .mutex = se_MUTEX_INIT, .cond = se_COND_INIT}

#define EVENT_STATE_RESET 0

typedef enum {EVENT_AUTORESET, EVENT_MANUALRESET} se_event_type_t;

void se_event_init       (se_event_t* ev, se_event_type_t type);
void se_event_init_shared(se_event_t* ev, se_event_type_t type);
void se_event_destroy    (se_event_t* ev);

void se_event_set      (se_event_t* ev);
void se_event_reset    (se_event_t* ev);
void se_event_set_state(se_event_t* ev, int state);
void se_event_inc_state(se_event_t* ev);
void se_event_dec_state(se_event_t* ev);

int  se_event_state (se_event_t* ev);
int  se_event_is_set(se_event_t* ev);

int  se_event_wait        (se_event_t* ev);
int  se_event_wait_timeout(se_event_t* ev, uint64_t millis);

void se_event_wait_state        (se_event_t* ev, int state);
int  se_event_wait_state_timeout(se_event_t* ev, int state, uint64_t millis);

int  se_event_wait_multistate        (se_event_t* ev, const int states[], int state_no);
int  se_event_wait_multistate_timeout(se_event_t* ev, const int states[], int state_no,
                                       uint64_t millis);

se_event_type_t se_event_type(const se_event_t* ev);


/*
  multi-readers/single-writer lock

  - synchronization object.

  - when locking the writer thread(s) have priority over the reader threads.
*/

typedef struct
{
    pthread_rwlock_t rwlock;
} se_rwlock_t;

#define se_RWLOCK_INIT {.rwlock = PTHREAD_RWLOCK_INITIALIZER}

void se_rwlock_init       (se_rwlock_t* rwl);
void se_rwlock_init_shared(se_rwlock_t* rwl);
void se_rwlock_destroy    (se_rwlock_t* rwl);

void se_rwlock_rdlock        (se_rwlock_t* rwl);
int  se_rwlock_tryrdlock     (se_rwlock_t* rwl);
int  se_rwlock_rdlock_timeout(se_rwlock_t* rwl, uint64_t millis);

void se_rwlock_wrlock        (se_rwlock_t* rwl);
int  se_rwlock_trywrlock     (se_rwlock_t* rwl);
int  se_rwlock_wrlock_timeout(se_rwlock_t* rwl, uint64_t millis);

void se_rwlock_unlock(se_rwlock_t* rwl);


/*
  barrier

  - synchronization object:

  - a barrier blocks the execution of all threads calling the barrier_wait() function
    until the number of waiting threads reaches the 'thr_count' value specified when
    the barrier was initialized.
  - the barrier_wait() call returns true (non-zero) to a single arbitrary selected
    threads, and zero to the rest of calling threads.
*/

typedef struct
{
    int count;
    int hit_count;
    se_event_t event;
    se_mutex_t mutex;
} se_barrier_t;

#define se_BARRIER_INIT(thr_count) {.count = thr_count, .hit_count = 0,  \
                                     .event = se_EVENT_MANUALRESET_INIT, \
                                     .mutex = se_MUTEX_INIT}

void se_barrier_init       (se_barrier_t* b, int thr_count);
void se_barrier_init_shared(se_barrier_t* b, int thr_count);
void se_barrier_destroy    (se_barrier_t* b);

int  se_barrier_wait(se_barrier_t* b);


/*
  synchronized [priority] queue

  - synchronization object:

  - a synchronized queue is an abstract FIFO data type.
  - if the 'priority' param passed to se_create_syncq() is not NULL, then you will
    have a priority queue.
  - when an value is enqueued using the syncq_signal() function, a single arbitrary
    selected waiting thread is released.
*/

typedef void* se_syncq;

se_syncq se_create_syncq(int max_size, se_value_priority_fn priority);
void se_destroy_syncq( se_syncq _sq, int free_values );

int  se_syncq_is_empty(se_syncq sq);
int  se_syncq_size    (se_syncq sq);

int  se_syncq_signal(se_syncq sq, void* value);

void se_syncq_wait        (se_syncq sq, void** value);
int  se_syncq_wait_timeout(se_syncq sq, void** value, uint64_t millis);

void se_syncq_clear(se_syncq sq, int free_values);

int  se_syncq_remove_matches(se_syncq sq,
                              se_compare_values_fn compare,
                              const void* key, int free_values);

int  se_syncq_contains(se_syncq sq,
                        se_compare_values_fn compare,
                        const void* key);

int  se_syncq_count_matches(se_syncq sq,
                             se_compare_values_fn compare,
                             const void* key);

int  se_syncq_max_size (se_syncq sq);
int  se_syncq_peak_size(se_syncq sq);

void se_syncq_reset_peak_size(se_syncq sq);

void se_syncq_process(se_syncq sq, void(*process)(void* _value));


/*
  join synchronized queue

  - complex object:

  - a join synchronized queue collects threads descriptors of ended joinable threads.
  - other threads can wait or directly join the threads that have their descriptors
    in the queue.
*/

typedef void* se_joinq;

se_joinq se_create_joinq(int max_thr_count);
void se_destroy_joinq(se_joinq jq);

int  se_joinq_size    (se_joinq jq);
int  se_joinq_is_empty(se_joinq jq);

void se_joinq_signal_self(se_joinq jq);
void se_joinq_signal     (se_joinq jq, pthread_t thr);

void se_joinq_wait        (se_joinq jq, pthread_t* thr);
int  se_joinq_wait_timeout(se_joinq jq, pthread_t* thr, uint64_t millis);

void se_joinq_join        (se_joinq jq);
int  se_joinq_join_timeout(se_joinq jq, uint64_t millis);

void se_joinq_joinall        (se_joinq jq, int thr_count);
int  se_joinq_joinall_timeout(se_joinq jq, int thr_count, uint64_t millis);


/*
  worker thread pool

  - complex object:

  - priority-queue based pool of detached (background) threads.
  - each thread pool job have a function, a parameter and a name.
  - min_count, max_count and reserve_idle_count arguments sets the thread
    allocation policy: the poll will try to maintain reserve_idle_count
    threads idle.
  - pool's enqueued jobs can be queried or manipulated.
  - pool's executing jobs can not be queried or manipulated.
  - to wait for the completion of all running jobs, use the thrpool_wait*()
    functions.
*/

#define THRPOOL_JOBNAME_MAXLEN 128
#define THRPOOL_UNNAMEDJOB "(unnamed-job)"

typedef void* se_thrpooljob;

void* se_thrpooljob_get_param(se_thrpooljob job);
const char* se_thrpooljob_get_name(se_thrpooljob job);


typedef void* se_thrpool;

se_thrpool se_create_thrpool(int min_count, int max_count, int reserve_idle_count,
                               se_value_priority_fn priority, int queue_max_size);

void se_thrpool_shutdown             (se_thrpool tp);
void se_thrpool_shutdown_wait        (se_thrpool tp);
int  se_thrpool_shutdown_wait_timeout(se_thrpool tp, uint64_t millis);

int  se_thrpool_queue_put(se_thrpool tp,
                           se_thrpool_job_fn job_fn, void* job_param, const char* job_name,
                           int wait_if_queue_is_full, uint64_t wait_timeout_millis);
int  se_thrpool_queue_contains      (se_thrpool tp, se_compare_values_fn compare, const void* key);
int  se_thrpool_queue_count_matches (se_thrpool tp, se_compare_values_fn compare, const void* key);
int  se_thrpool_queue_remove_matches(se_thrpool tp, se_compare_values_fn compare, const void* key);
void se_thrpool_queue_process       (se_thrpool tp, se_process_value_fn process);
void se_thrpool_queue_clear         (se_thrpool tp);

void se_thrpool_waitall        (se_thrpool tp);
int  se_thrpool_waitall_timeout(se_thrpool tp, uint64_t millis);

int  se_thrpool_queue_size        (se_thrpool tp);
int  se_thrpool_queue_peak_size   (se_thrpool tp);
void se_thrpool_reset_queue_peak  (se_thrpool tp);
int  se_thrpool_queue_max_size    (se_thrpool tp);
void se_thrpool_set_queue_max_size(se_thrpool tp, int queue_max_count);

void se_thrpool_thread_status    (se_thrpool tp, int* thread_count, int* busy_count, int* idle_count);
void se_thrpool_thread_counts    (se_thrpool tp, int* min_count, int* max_count, int* reserve_idle_count);
void se_thrpool_set_thread_counts(se_thrpool tp, int min_count, int max_count, int reserve_idle_count);


/*
  thread utils

  - simple wrappers for pthread_*() functions.
*/

void se_thr_sleep (uint64_t millis);
void se_thr_detach(pthread_t thr);


void se_thr_cancel   (pthread_t thr);
void se_thr_cancelall(const pthread_t* thr, int thr_count);


void se_thr_join(pthread_t thr, void** ret_ptr);


#define THR_ATTR_JOINABLE       0x0001
#define THR_ATTR_DETACHED       0x0002
#define THR_ATTR_NOCANCEL       0x0004
#define THR_ATTR_CANCELDEFERRED 0x0008
#define THR_ATTR_CANCELASYNC    0x0010
#define THR_ATTR_DEFAULT        (THR_ATTR_JOINABLE | THR_ATTR_CANCELDEFERRED)

void se_create_thr( pthread_t* pthr, uint32_t thr_attr, void* (pthr_fn)(void*),
                     void* pthr_param, se_barrier_t* start_b, se_joinq jq );


#else /* se_USE_PTHREADS */



#endif /* se_USE_PTHREADS */


#endif