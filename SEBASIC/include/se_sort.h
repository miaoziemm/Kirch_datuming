#ifndef SE_SORT_H
#define SE_SORT_H
#include <stdint.h>
#include <stddef.h>
#include "se_byte.h"
#include "se_type.h"
#include "se_alloc.h"
#include "se_minmax.h"
#include "se_assert.h"
#include "se_basic_math.h"


/** Generic signature for a function that can be used to sort float arrays
 *  \param x   - the float array to be sorted
 *  \param off - the offset in the array to start sorting
 *  \param len - the number of elements to sort
 */
typedef void (*se_sortf_fc)(float* x, size_t off, size_t len);
void se_insert_sortf(float* x, size_t off, size_t len);
void se_quick_sortf(float* x, size_t off, size_t len);
void se_bubble_sortf(float* x, size_t off, size_t len);
ssize_t se_binary_searchf(float key, const float* a, size_t from, size_t to);

typedef void (*se_sortd_fc)(double* x, size_t off, size_t len);
void se_insert_sortd(double* x, size_t off, size_t len);
void se_quick_sortd(double* x, size_t off, size_t len);
void se_bubble_sortd(double* x, size_t off, size_t len);
ssize_t se_binary_searchd(double key, const double* a, size_t from, size_t to);

typedef void (*se_sortr_fc)(real* x, int off, int len);
void se_insert_sortr(real* x, size_t off, size_t len);
void se_quick_sortr(real* x, size_t off, size_t len);
void se_bubble_sortr(real* x, size_t off, size_t len);
ssize_t se_binary_searchr(real key, const real* a, size_t from, size_t to);

typedef void (*se_sortra_fc)(reala* x, size_t off, size_t len);
void se_insert_sortra(reala* x, size_t off, size_t len);
void se_quick_sortra(reala* x, size_t off, size_t len);
void se_bubble_sortra(reala* x, size_t off, size_t len);
ssize_t se_binary_searchra(reala key, const reala* a, size_t from, size_t to);

typedef void (*se_sorti_fc)(int* x, size_t off, size_t len);
void se_insert_sorti(int* x, size_t off, size_t len);
void se_quick_sorti(int* x, size_t off, size_t len);
void se_bubble_sorti(int* x, size_t off, size_t len);
ssize_t se_binary_searchi(int key, const int* a, size_t from, size_t to);

typedef void (*se_sorti64_fc)(int64_t* x, size_t off, size_t len);
void se_insert_sorti64(int64_t* x, size_t off, size_t len);
void se_quick_sorti64(int64_t* x, size_t off, size_t len);
void se_bubble_sorti64(int64_t* x, size_t off, size_t len);
ssize_t se_binary_searchi64(int64_t key, const int64_t* a, size_t from, size_t to);

typedef void (*se_sortsz_fc)(size_t* x, size_t off, size_t len);
void se_insert_sortsz(size_t* x, size_t off, size_t len);
void se_quick_sortsz(size_t* x, size_t off, size_t len);
void se_bubble_sortsz(size_t* x, size_t off, size_t len);
ssize_t se_binary_searchsz(size_t key, const size_t* a, size_t from, size_t to);

typedef void (*se_sortb_fc)(byte* x, size_t off, size_t len);
void se_insert_sortb(byte* x, size_t off, size_t len);
void se_quick_sortb(byte* x, size_t off, size_t len);
void se_bubble_sortb(byte* x, size_t off, size_t len);
ssize_t se_binary_searchb(byte key, const byte* a, size_t from, size_t to);

/*some specialization routines */
typedef struct se_int_pairs_s {
    int key;
    int val;
} se_int_keyval_t;

typedef void (*se_sort_int_keyval_fc)(se_int_keyval_t* x, size_t off, size_t len);
void se_insert_sort_int_keyval(se_int_keyval_t* x, size_t off, size_t len);
void se_quick_sort_int_keyval(se_int_keyval_t* x, size_t off, size_t len);
void se_bubble_sort_int_keyval(se_int_keyval_t* x, size_t off, size_t len);
ssize_t se_binary_search_int_keyval(se_int_keyval_t key, const se_int_keyval_t* a, size_t from, size_t to);

typedef struct se_int64_pairs_s {
    int64_t key;
    int64_t val;
} se_int64_keyval_t;

typedef void (*se_sort_int64_keyval_fc)(se_int64_keyval_t* x, size_t off, size_t len);
void se_insert_sort_int64_keyval(se_int64_keyval_t* x, size_t off, size_t len);
void se_quick_sort_int64_keyval(se_int64_keyval_t* x, size_t off, size_t len);
void se_bubble_sort_int64_keyval(se_int64_keyval_t* x, size_t off, size_t len);
ssize_t se_binary_search_int64_keyval(se_int64_keyval_t key, const se_int64_keyval_t* a, 
                                       size_t from, size_t to);

typedef struct se_int64_double_s {
    double key;
    int64_t val;
} se_int64_double_keyval_t;

typedef void (*se_sort_int64_double_keyval_fc)(se_int64_double_keyval_t* x, size_t off, size_t len);
void se_insert_sort_int64_double_keyval(se_int64_double_keyval_t* x, size_t off, size_t len);
void se_quick_sort_int64_double_keyval(se_int64_double_keyval_t* x, size_t off, size_t len);
void se_bubble_sort_int64_double_keyval(se_int64_double_keyval_t* x, size_t off, size_t len);
ssize_t se_binary_search_int64_double_keyval(se_int64_double_keyval_t key, const se_int64_double_keyval_t* a,
                                              size_t from, size_t to);

typedef struct se_double_int64_s {
    int64_t key;
    double val;
} se_double_int64_keyval_t;

typedef void (*se_sort_double_int64_keyval_fc)(se_double_int64_keyval_t* x, size_t off, size_t len);
void se_insert_sort_double_int64_keyval(se_double_int64_keyval_t* x, size_t off, size_t len);
void se_quick_sort_double_int64_keyval(se_double_int64_keyval_t* x, size_t off, size_t len);
void se_bubble_sort_double_int64_keyval(se_double_int64_keyval_t* x, size_t off, size_t len);
ssize_t se_binary_search_double_int64_keyval(se_double_int64_keyval_t key, const se_double_int64_keyval_t* a,
                                              size_t from, size_t to);




#endif