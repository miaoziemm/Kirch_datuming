#include "../include/se_sort.h"


#define se_sort_lt(a, b) ((a) < (b))
#define se_sort_gt(a, b) ((a) > (b))
#define se_sort_eq(a, b) FEQUAL((a), (b))

#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#define se_sort_type float
#define se_sort_name(name) name ## f
#include "se_sort_impl.cpp"

#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#define se_sort_type double
#define se_sort_name(name) name ## d
#include "se_sort_impl.cpp"

#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#define se_sort_type real
#define se_sort_name(name) name ## r
#include "se_sort_impl.cpp"

#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#define se_sort_type reala
#define se_sort_name(name) name ## ra
#include "se_sort_impl.cpp"


#ifdef se_sort_eq
#undef se_sort_eq
#endif
#define se_sort_eq(a, b) ((a) == (b))


#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#define se_sort_type int
#define se_sort_name(name) name ## i
#include "se_sort_impl.cpp"


#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#define se_sort_type int64_t
#define se_sort_name(name) name ## i64
#include "se_sort_impl.cpp"


#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#define se_sort_type size_t
#define se_sort_name(name) name ## sz
#include "se_sort_impl.cpp"


#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#define se_sort_type byte
#define se_sort_name(name) name ## b
#include "se_sort_impl.cpp"

// se_int64_keyval_t

#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#ifdef se_sort_lt
#undef se_sort_lt
#endif
#ifdef se_sort_gt
#undef se_sort_gt
#endif
#ifdef se_sort_eq
#undef se_sort_eq
#endif

#define se_sort_lt(a, b) ((a).key < (b).key)
#define se_sort_gt(a, b) ((a).key > (b).key)
#define se_sort_eq(a, b) ((a).key == (b).key)

#define se_sort_type se_int_keyval_t
#define se_sort_name(name) name ## _int_keyval
#include "se_sort_impl.cpp"

#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#define se_sort_type se_int64_keyval_t
#define se_sort_name(name) name ## _int64_keyval
#include "se_sort_impl.cpp"

// se_int64_double_keyval_t

#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#ifdef se_sort_lt
#undef se_sort_lt
#endif
#ifdef se_sort_gt
#undef se_sort_gt
#endif
#ifdef se_sort_eq
#undef se_sort_eq
#endif

#define se_sort_lt(a, b) ((a).key < (b).key)
#define se_sort_gt(a, b) ((a).key > (b).key)
#define se_sort_eq(a, b) (FEQUAL((a).key,(b).key))

#define se_sort_type se_int64_double_keyval_t
#define se_sort_name(name) name ## _int64_double_keyval
#include "se_sort_impl.cpp"

// se_double_int64_keyval_t

#ifdef se_sort_type
#undef se_sort_type
#endif
#ifdef se_sort_name
#undef se_sort_name
#endif
#ifdef se_sort_lt
#undef se_sort_lt
#endif
#ifdef se_sort_gt
#undef se_sort_gt
#endif
#ifdef se_sort_eq
#undef se_sort_eq
#endif

#define se_sort_lt(a, b) ((a).key < (b).key)
#define se_sort_gt(a, b) ((a).key > (b).key)
#define se_sort_eq(a, b) ((a).key == (b).key)

#define se_sort_type se_double_int64_keyval_t
#define se_sort_name(name) name ## _double_int64_keyval
#include "se_sort_impl.cpp"
