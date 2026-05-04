#ifndef __SE_MINMAX_H__
#define __SE_MINMAX_H__

#include "se_type.h"

#ifdef ARCH_OSX
#include <sys/types.h>
#endif

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type double
#define minmax_name(name) name ## d
#include "se_minmax_impl.h"


#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type real
#define minmax_name(name) name ## r
#include "se_minmax_impl.h"

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type float
#define minmax_name(name) name ## f
#include "se_minmax_impl.h"

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type int
#define minmax_name(name) name ## i
#include "se_minmax_impl.h"

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type int64_t
#define minmax_name(name) name ## i64
#include "se_minmax_impl.h"

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type size_t
#define minmax_name(name) name ## sz
#include "se_minmax_impl.h"

#ifdef minmax_type
#undef minmax_type
#endif
#ifdef minmax_name
#undef minmax_name
#endif
#define minmax_type ssize_t
#define minmax_name(name) name ## ssz
#include "se_minmax_impl.h"

#endif