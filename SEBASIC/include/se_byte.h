#ifndef SE_BYTE_H
#define SE_BYTE_H
#include <string.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>
#include "se_type.h"
#include "se_log.h"
#include "se_assert.h"

// 运行时字节序检测函数
inline bool is_little_endian() {
    union {
        uint32_t i;
        char c[4];
    } test = {0x01020304};
    return test.c[0] == 4;  // 小端：低位字节在前
}

inline bool is_big_endian() {
    return !is_little_endian();
}

// 编译时字节序检测（优先使用）
#if defined(__BYTE_ORDER__)
    // GCC/Clang 预定义宏
    #if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
        #define SE_IS_LITTLE_ENDIAN 1
        #define SE_IS_BIG_ENDIAN 0
    #elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
        #define SE_IS_LITTLE_ENDIAN 0
        #define SE_IS_BIG_ENDIAN 1
    #else
        #define SE_IS_LITTLE_ENDIAN (is_little_endian())
        #define SE_IS_BIG_ENDIAN (is_big_endian())
    #endif
#elif defined(_WIN32) || defined(_WIN64)
    // Windows 通常是小端
    #define SE_IS_LITTLE_ENDIAN 1
    #define SE_IS_BIG_ENDIAN 0
#elif defined(__i386__) || defined(__x86_64__) || defined(__amd64__) || \
      defined(__arm__) || defined(__aarch64__)
    // 常见的小端架构
    #define SE_IS_LITTLE_ENDIAN 1
    #define SE_IS_BIG_ENDIAN 0
#elif defined(__sparc__) || defined(__powerpc__) || defined(__ppc__)
    // 常见的大端架构
    #define SE_IS_LITTLE_ENDIAN 0
    #define SE_IS_BIG_ENDIAN 1
#else
    // 无法在编译时确定，使用运行时检测
    #define SE_IS_LITTLE_ENDIAN (is_little_endian())
    #define SE_IS_BIG_ENDIAN (is_big_endian())
#endif

// 兼容原有的字节序定义
#define	__LITTLE_ENDIAN	1234
#define	__BIG_ENDIAN	4321

// 根据检测结果设置 __BYTE_ORDER
#if defined(BYTE_ORDER_LE) || SE_IS_LITTLE_ENDIAN
    #define __BYTE_ORDER __LITTLE_ENDIAN
#elif defined(BYTE_ORDER_BE) || SE_IS_BIG_ENDIAN
    #define __BYTE_ORDER __BIG_ENDIAN
#else
    // 运行时检测作为最后的后备方案
    #define __BYTE_ORDER (is_little_endian() ? __LITTLE_ENDIAN : __BIG_ENDIAN)
#endif




/** Enumeration to define endiannes, or byte order. */
typedef enum {
    /** Big endian, or network order: most significant bit has the lowest addresses */
    big_endian,
    /** Little endian: most significant bit has the highest addresses */
    little_endian
} byte_order ;

#if __BYTE_ORDER == __LITTLE_ENDIAN
#   define DEFAULT_BYTE_ORDER little_endian
#   define NEEDS_SWAPPING(bo) ( (bo) != little_endian )

#else
#if __BYTE_ORDER == __BIG_ENDIAN
#   define DEFAULT_BYTE_ORDER big_endian
#   define NEEDS_SWAPPING(bo) ( (bo) != big_endian )

#else
#error "Please fix <endian.h>"
#endif
#endif




/** Swaps bytes in 16 bit value.  */
#define bswap_const_16(x)                           \
    ( (((x) >> 8) & 0xff) | (((x) & 0xff) << 8) )


/** Swaps bytes in 32 bit value.  */
#define bswap_const_32(x) (                                         \
                               ( ( (x) & (uint32_t)0xff000000 ) >> 24 ) | \
                               ( ( (x) & (uint32_t)0x00ff0000 ) >>  8 ) | \
                               ( ( (x) & (uint32_t)0x0000ff00 ) <<  8 ) | \
                               ( ( (x) & (uint32_t)0x000000ff ) << 24 ) )

/** Swaps bytes in 64 bit value.  */
#define bswap_const_64(x) (                                         \
                               ( ( (x) & (uint64_t)0xff00000000000000 ) >> 56 ) | \
                               ( ( (x) & (uint64_t)0x00ff000000000000 ) >> 40 ) | \
                               ( ( (x) & (uint64_t)0x0000ff0000000000 ) >> 24 ) | \
                               ( ( (x) & (uint64_t)0x000000ff00000000 ) >>  8 ) | \
                               ( ( (x) & (uint64_t)0x00000000ff000000 ) <<  8 ) | \
                               ( ( (x) & (uint64_t)0x0000000000ff0000 ) << 24 ) | \
                               ( ( (x) & (uint64_t)0x000000000000ff00 ) << 40 ) | \
                               ( ( (x) & (uint64_t)0x00000000000000ff ) << 56 ) )


/**
 * If necessary, change the byte order of the input array,
 * swaping every 2 bytes.
 *
 * \param[in,out] buff - the buffer to be swapped in place.
 *
 * \param[in] n - the size of the memory block pointed by buff that needs swapping.
 *
 * \param[in] bo_input - the desired order; if the current
 *                       architecture has the same order as the
 *                       desired order, then the function returns
 *                       immediately.
 */
void order_bytes_2(byte* buff, const size_t n, const byte_order bo_input);

/**
 * If necessary, change the byte order of the input array,
 * swaping every 4 bytes: (0,3) and (1,2).
 *
 * \param[in,out] buff - the buffer to be swapped in place.
 *
 * \param[in] n - the size of the memory block pointed by buff that needs swapping.
 *
 * \param[in] bo_input - the desired order; if the current
 *                       architecture has the same order as the
 *                       desired order, then the function returns
 *                       immediately.
 */
void order_bytes_4(byte* buff, const size_t n, const byte_order bo_input);


/**
 * If necessary, change the byte order of the input array,
 * swaping every 8 bytes: (0,7), (1,6), (2,5) and (3,4)
 *
 * \param[in,out] buff - the buffer to be swapped in place.
 *
 * \param[in] n - the size of the memory block pointed by buff that needs swapping.
 *
 * \param[in] bo_input - the desired order; if the current
 *                       architecture has the same order as the
 *                       desired order, then the function returns
 *                       immediately.
 */
void order_bytes_8(byte* buff, const size_t n, const byte_order bo_input);


/**
 * This is the reference implementation for 2-byte swapping. It
 * doesn't use any compiler-specific or library-specific
 * optimizations.
 *
 * \param[in,out] buff - the buffer to be swapped in place.
 *
 * \param[in] n - the size of the memory block pointed by buff that needs swapping.
 *
 * \param[in] bo_input - the desired order; if the current
 *                       architecture has the same order as the
 *                       desired order, then the function returns
 *                       immediately.
 */
void reference_order_bytes_2( byte* buff,
                                  const size_t n,
                                  const byte_order bo_input );

/**
 * This is the reference implementation for 4-byte swapping. It
 * doesn't use any compiler-specific or library-specific
 * optimizations.
 *
 * \param[in,out] buff - the buffer to be swapped in place.
 *
 * \param[in] n - the size of the memory block pointed by buff that needs swapping.
 *
 * \param[in] bo_input - the desired order; if the current
 *                       architecture has the same order as the
 *                       desired order, then the function returns
 *                       immediately.
 */
void reference_order_bytes_4( byte* buff,
                                  const size_t n,
                                  const byte_order bo_input );

/**
 * This is the reference implementation for 8-byte swapping. It
 * doesn't use any compiler-specific or library-specific
 * optimizations.
 *
 * \param[in,out] buff - the buffer to be swapped in place.
 *
 * \param[in] n - the size of the memory block pointed by buff that needs swapping.
 *
 * \param[in] bo_input - the desired order; if the current
 *                       architecture has the same order as the
 *                       desired order, then the function returns
 *                       immediately.
 */
void reference_order_bytes_8( byte* buff,
                                  const size_t n,
                                  const byte_order bo_input );


void ibm_to_ieee_sep(const unsigned int* from,
                         unsigned int* to, int n, int endian);
void ieee_to_ibm_sep(const unsigned int* from,
                         unsigned int* to, int n, int endian);

void ibm_to_ieee(const unsigned int* from,
                     unsigned int* to, int n, int endian);
void ieee_to_ibm(const unsigned int* from,
                     unsigned int* to, int n, int endian);

#endif