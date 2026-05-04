#ifndef SE_TYPE_H
#define SE_TYPE_H

#include <cstdint>


/** This is the type used for all computations */
#ifdef USE_DOUBLE_DOUBLE_PRECISION
typedef long double real;
#else
#ifdef USE_DOUBLE_PRECISION
typedef double real;
#else
typedef double real;
#endif
#endif


/** This is the type used to store in arrays - memory footprint is important */
typedef float reala;

typedef uint8_t byte;

typedef float se_complex[2];

/* pointer to function types */
typedef void     (*process_value_fn)  (void* value);
typedef int      (*compare_values_fn) (const void* key, const void* value);
typedef uint32_t (*value_priority_fn) (const void* value);
typedef int      (*process_path_fn)   (const char* path, int is_dir, void* arg);
typedef void     (*thrpool_job_fn)    (void* job_param);
typedef void     (*thr_cleanup_fn)    (void* _param);


inline void swap_real(real* a, real* b)
{
    real tmp = *a;
    *a = *b;
    *b = tmp;
}

inline void swapf(float *a, float *b)
{
    float tmp = *a;
    *a = *b;
    *b = tmp;
}

inline void swap_reala(reala* a, reala* b)
{
    reala tmp = *a;
    *a = *b;
    *b = tmp;
}


#endif