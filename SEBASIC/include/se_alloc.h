#ifndef SE_ALLOC_H
#define SE_ALLOC_H

#include <cstdio>
#include <cstdlib>
#include <cstring>
/* Avoid including C complex header here to prevent macro I pollution in C++ builds */

#ifdef __cplusplus
extern "C" {
#endif

/* Function declarations for basic allocation routines */
void *alloc1(size_t n1, size_t size);
void *realloc1(void *v, size_t n1, size_t size);
void free1(void *p);

void **alloc2(size_t n1, size_t n2, size_t size);
void free2(void **p);

void ***alloc3(size_t n1, size_t n2, size_t n3, size_t size);
void free3(void ***p);

void ****alloc4(size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void free4(void ****p);

void *****alloc5(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size);
void free5(void *****p);

void ******alloc6(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6, size_t size);
void free6(void ******p);

/* Zero-initialized allocation function declarations */
void *alloc1_zero(size_t n1, size_t size);
void **alloc2_zero(size_t n1, size_t n2, size_t size);
void ***alloc3_zero(size_t n1, size_t n2, size_t n3, size_t size);
void ****alloc4_zero(size_t n1, size_t n2, size_t n3, size_t n4, size_t size);
void *****alloc5_zero(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t size);
void ******alloc6_zero(size_t n1, size_t n2, size_t n3, size_t n4, size_t n5, size_t n6, size_t size);

#ifdef __cplusplus
}
#endif

/* Generic macros for type-safe memory allocation */


/* 1D array allocation macros */
#define ALLOC1(type, n1) ((type*)alloc1((n1), sizeof(type)))
#define REALLOC1(type, ptr, n1) ((type*)realloc1((ptr), (n1), sizeof(type)))
#define FREE1(ptr) free1(ptr)

/* 2D array allocation macros */
#define ALLOC2(type, n1, n2) ((type**)alloc2((n1), (n2), sizeof(type)))
#define FREE2(ptr) free2((void**)(ptr))

/* 3D array allocation macros */
#define ALLOC3(type, n1, n2, n3) ((type***)alloc3((n1), (n2), (n3), sizeof(type)))
#define FREE3(ptr) free3((void***)(ptr))

/* 4D array allocation macros */
#define ALLOC4(type, n1, n2, n3, n4) ((type****)alloc4((n1), (n2), (n3), (n4), sizeof(type)))
#define FREE4(ptr) free4((void****)(ptr))

/* 5D array allocation macros */
#define ALLOC5(type, n1, n2, n3, n4, n5) ((type*****)alloc5((n1), (n2), (n3), (n4), (n5), sizeof(type)))
#define FREE5(ptr) free5((void*****)(ptr))

/* 6D array allocation macros */
#define ALLOC6(type, n1, n2, n3, n4, n5, n6) ((type******)alloc6((n1), (n2), (n3), (n4), (n5), (n6), sizeof(type)))
#define FREE6(ptr) free6((void******)(ptr))

/* Zero-initialized allocation macros */

/* 1D array zero allocation macros */
#define ALLOC1_ZERO(type, n1) ((type*)alloc1_zero((n1), sizeof(type)))
#define FREE1_ZERO(ptr) free1(ptr)

/* 2D array zero allocation macros */
#define ALLOC2_ZERO(type, n1, n2) ((type**)alloc2_zero((n1), (n2), sizeof(type)))
#define FREE2_ZERO(ptr) free2((void**)(ptr))

/* 3D array zero allocation macros */
#define ALLOC3_ZERO(type, n1, n2, n3) ((type***)alloc3_zero((n1), (n2), (n3), sizeof(type)))
#define FREE3_ZERO(ptr) free3((void***)(ptr))

/* 4D array zero allocation macros */
#define ALLOC4_ZERO(type, n1, n2, n3, n4) ((type****)alloc4_zero((n1), (n2), (n3), (n4), sizeof(type)))
#define FREE4_ZERO(ptr) free4((void****)(ptr))

/* 5D array zero allocation macros */
#define ALLOC5_ZERO(type, n1, n2, n3, n4, n5) ((type*****)alloc5_zero((n1), (n2), (n3), (n4), (n5), sizeof(type)))
#define FREE5_ZERO(ptr) free5((void*****)(ptr))

/* 6D array zero allocation macros */
#define ALLOC6_ZERO(type, n1, n2, n3, n4, n5, n6) ((type******)alloc6_zero((n1), (n2), (n3), (n4), (n5), (n6), sizeof(type)))
#define FREE6_ZERO(ptr) free6((void******)(ptr))


/* Generic type allocation macros - for any type (struct, union, typedef, etc.) */
#define alloc1type(type_name, count) ((type_name*)alloc1((count), sizeof(type_name)))
#define alloc2type(type_name, n1, n2) ((type_name**)alloc2((n1), (n2), sizeof(type_name)))
#define alloc3type(type_name, n1, n2, n3) ((type_name***)alloc3((n1), (n2), (n3), sizeof(type_name)))
#define alloc4type(type_name, n1, n2, n3, n4) ((type_name****)alloc4((n1), (n2), (n3), (n4), sizeof(type_name)))
#define alloc5type(type_name, n1, n2, n3, n4, n5) ((type_name*****)alloc5((n1), (n2), (n3), (n4), (n5), sizeof(type_name)))
#define alloc6type(type_name, n1, n2, n3, n4, n5, n6) ((type_name******)alloc6((n1), (n2), (n3), (n4), (n5), (n6), sizeof(type_name)))

#define realloc1type(type_name, ptr, count) ((type_name*)realloc1((ptr), (count), sizeof(type_name)))
/* Note: realloc for multi-dimensional arrays (2D, 3D, etc.) is not provided 
   as it requires complex memory layout reorganization. Use alloc + copy + free instead. */

#define free1type(ptr) free1(ptr)
#define free2type(ptr) free2((void**)(ptr))
#define free3type(ptr) free3((void***)(ptr))
#define free4type(ptr) free4((void****)(ptr))
#define free5type(ptr) free5((void*****)(ptr))
#define free6type(ptr) free6((void******)(ptr))

/* Zero-initialized type allocation macros - consistent with alloc1type naming style */
#define alloc1type_zero(type_name, count) ((type_name*)alloc1_zero((count), sizeof(type_name)))
#define alloc2type_zero(type_name, n1, n2) ((type_name**)alloc2_zero((n1), (n2), sizeof(type_name)))
#define alloc3type_zero(type_name, n1, n2, n3) ((type_name***)alloc3_zero((n1), (n2), (n3), sizeof(type_name)))
#define alloc4type_zero(type_name, n1, n2, n3, n4) ((type_name****)alloc4_zero((n1), (n2), (n3), (n4), sizeof(type_name)))
#define alloc5type_zero(type_name, n1, n2, n3, n4, n5) ((type_name*****)alloc5_zero((n1), (n2), (n3), (n4), (n5), sizeof(type_name)))
#define alloc6type_zero(type_name, n1, n2, n3, n4, n5, n6) ((type_name******)alloc6_zero((n1), (n2), (n3), (n4), (n5), (n6), sizeof(type_name)))

/* Note: Use the same free functions as regular allocation - freeing is the same regardless of initialization */
#define free1type_zero(ptr) free1type(ptr)
#define free2type_zero(ptr) free2type(ptr)
#define free3type_zero(ptr) free3type(ptr)
#define free4type_zero(ptr) free4type(ptr)
#define free5type_zero(ptr) free5type(ptr)
#define free6type_zero(ptr) free6type(ptr)

/* Zero-initialized struct allocation macros */
#define alloc1struct_zero(struct_type, count) ((struct struct_type*)alloc1_zero((count), sizeof(struct struct_type)))
#define alloc2struct_zero(struct_type, n1, n2) ((struct struct_type**)alloc2_zero((n1), (n2), sizeof(struct struct_type)))
#define alloc3struct_zero(struct_type, n1, n2, n3) ((struct struct_type***)alloc3_zero((n1), (n2), (n3), sizeof(struct struct_type)))
#define alloc4struct_zero(struct_type, n1, n2, n3, n4) ((struct struct_type****)alloc4_zero((n1), (n2), (n3), (n4), sizeof(struct struct_type)))
#define alloc5struct_zero(struct_type, n1, n2, n3, n4, n5) ((struct struct_type*****)alloc5_zero((n1), (n2), (n3), (n4), (n5), sizeof(struct struct_type)))
#define alloc6struct_zero(struct_type, n1, n2, n3, n4, n5, n6) ((struct struct_type******)alloc6_zero((n1), (n2), (n3), (n4), (n5), (n6), sizeof(struct struct_type)))

#define free1struct_zero(ptr) free1struct(ptr)
#define free2struct_zero(ptr) free2struct(ptr)
#define free3struct_zero(ptr) free3struct(ptr)
#define free4struct_zero(ptr) free4struct(ptr)
#define free5struct_zero(ptr) free5struct(ptr)
#define free6struct_zero(ptr) free6struct(ptr)

/* Backward compatibility macros for commonly used types */
#define alloc1int(n1) ALLOC1(int, n1)
#define realloc1int(ptr, n1) REALLOC1(int, ptr, n1)
#define free1int(ptr) FREE1(ptr)
#define alloc2int(n1, n2) ALLOC2(int, n1, n2)
#define free2int(ptr) FREE2(ptr)
#define alloc3int(n1, n2, n3) ALLOC3(int, n1, n2, n3)
#define free3int(ptr) FREE3(ptr)
#define alloc4int(n1, n2, n3, n4) ALLOC4(int, n1, n2, n3, n4)
#define free4int(ptr) FREE4(ptr)
#define alloc5int(n1, n2, n3, n4, n5) ALLOC5(int, n1, n2, n3, n4, n5)
#define free5int(ptr) FREE5(ptr)

#define alloc1float(n1) ALLOC1(float, n1)
#define realloc1float(ptr, n1) REALLOC1(float, ptr, n1)
#define free1float(ptr) FREE1(ptr)
#define alloc2float(n1, n2) ALLOC2(float, n1, n2)
#define free2float(ptr) FREE2(ptr)
#define alloc3float(n1, n2, n3) ALLOC3(float, n1, n2, n3)
#define free3float(ptr) FREE3(ptr)
#define alloc4float(n1, n2, n3, n4) ALLOC4(float, n1, n2, n3, n4)
#define free4float(ptr) FREE4(ptr)
#define alloc5float(n1, n2, n3, n4, n5) ALLOC5(float, n1, n2, n3, n4, n5)
#define free5float(ptr) FREE5(ptr)
#define alloc6float(n1, n2, n3, n4, n5, n6) ALLOC6(float, n1, n2, n3, n4, n5, n6)
#define free6float(ptr) FREE6(ptr)

#define alloc1double(n1) ALLOC1(double, n1)
#define realloc1double(ptr, n1) REALLOC1(double, ptr, n1)
#define free1double(ptr) FREE1(ptr)
#define alloc2double(n1, n2) ALLOC2(double, n1, n2)
#define free2double(ptr) FREE2(ptr)
#define alloc3double(n1, n2, n3) ALLOC3(double, n1, n2, n3)
#define free3double(ptr) FREE3(ptr)

#define alloc1char(n1) ALLOC1(char, n1)
#define realloc1char(ptr, n1) REALLOC1(char, ptr, n1)
#define free1char(ptr) FREE1(ptr)

/* Zero-initialized macros for commonly used types */
#define alloc1int_zero(n1) ALLOC1_ZERO(int, n1)
#define alloc2int_zero(n1, n2) ALLOC2_ZERO(int, n1, n2)
#define alloc3int_zero(n1, n2, n3) ALLOC3_ZERO(int, n1, n2, n3)
#define alloc4int_zero(n1, n2, n3, n4) ALLOC4_ZERO(int, n1, n2, n3, n4)
#define alloc5int_zero(n1, n2, n3, n4, n5) ALLOC5_ZERO(int, n1, n2, n3, n4, n5)
#define alloc6int_zero(n1, n2, n3, n4, n5, n6) ALLOC6_ZERO(int, n1, n2, n3, n4, n5, n6)
#define free1int_zero(ptr) FREE1_ZERO(ptr)
#define free2int_zero(ptr) FREE2_ZERO(ptr)
#define free3int_zero(ptr) FREE3_ZERO(ptr)
#define free4int_zero(ptr) FREE4_ZERO(ptr)
#define free5int_zero(ptr) FREE5_ZERO(ptr)
#define free6int_zero(ptr) FREE6_ZERO(ptr)

#define alloc1float_zero(n1) ALLOC1_ZERO(float, n1)
#define alloc2float_zero(n1, n2) ALLOC2_ZERO(float, n1, n2)
#define alloc3float_zero(n1, n2, n3) ALLOC3_ZERO(float, n1, n2, n3)
#define alloc4float_zero(n1, n2, n3, n4) ALLOC4_ZERO(float, n1, n2, n3, n4)
#define alloc5float_zero(n1, n2, n3, n4, n5) ALLOC5_ZERO(float, n1, n2, n3, n4, n5)
#define alloc6float_zero(n1, n2, n3, n4, n5, n6) ALLOC6_ZERO(float, n1, n2, n3, n4, n5, n6)
#define free1float_zero(ptr) FREE1_ZERO(ptr)
#define free2float_zero(ptr) FREE2_ZERO(ptr)
#define free3float_zero(ptr) FREE3_ZERO(ptr)
#define free4float_zero(ptr) FREE4_ZERO(ptr)
#define free5float_zero(ptr) FREE5_ZERO(ptr)
#define free6float_zero(ptr) FREE6_ZERO(ptr)

#define alloc1double_zero(n1) ALLOC1_ZERO(double, n1)
#define alloc2double_zero(n1, n2) ALLOC2_ZERO(double, n1, n2)
#define alloc3double_zero(n1, n2, n3) ALLOC3_ZERO(double, n1, n2, n3)
#define free1double_zero(ptr) FREE1_ZERO(ptr)
#define free2double_zero(ptr) FREE2_ZERO(ptr)
#define free3double_zero(ptr) FREE3_ZERO(ptr)

#define alloc1char_zero(n1) ALLOC1_ZERO(char, n1)
#define free1char_zero(ptr) FREE1_ZERO(ptr)

/* Complex number macros */
#ifdef __cplusplus
#include <complex>
#define alloc1complexf(n1) ALLOC1(std::complex<float>, n1)
#define realloc1complexf(ptr, n1) REALLOC1(std::complex<float>, ptr, n1)
#define free1complexf(ptr) FREE1(ptr)
#define alloc2complexf(n1, n2) ALLOC2(std::complex<float>, n1, n2)
#define free2complexf(ptr) FREE2(ptr)
#define alloc3complexf(n1, n2, n3) ALLOC3(std::complex<float>, n1, n2, n3)
#define free3complexf(ptr) FREE3(ptr)
#define alloc4complexf(n1, n2, n3, n4) ALLOC4(std::complex<float>, n1, n2, n3, n4)
#define free4complexf(ptr) FREE4(ptr)

#define alloc1complex(n1) ALLOC1(std::complex<double>, n1)
#define realloc1complex(ptr, n1) REALLOC1(std::complex<double>, ptr, n1)
#define free1complex(ptr) FREE1(ptr)
#define alloc2complex(n1, n2) ALLOC2(std::complex<double>, n1, n2)
#define free2complex(ptr) FREE2(ptr)
#define alloc3complex(n1, n2, n3) ALLOC3(std::complex<double>, n1, n2, n3)
#define free3complex(ptr) FREE3(ptr)
#define alloc4complex(n1, n2, n3, n4) ALLOC4(std::complex<double>, n1, n2, n3, n4)
#define free4complex(ptr) FREE4(ptr)

/* Zero-initialized complex number macros for C++ */
#define alloc1complexf_zero(n1) ALLOC1_ZERO(std::complex<float>, n1)
#define alloc2complexf_zero(n1, n2) ALLOC2_ZERO(std::complex<float>, n1, n2)
#define alloc3complexf_zero(n1, n2, n3) ALLOC3_ZERO(std::complex<float>, n1, n2, n3)
#define alloc4complexf_zero(n1, n2, n3, n4) ALLOC4_ZERO(std::complex<float>, n1, n2, n3, n4)
#define free1complexf_zero(ptr) FREE1_ZERO(ptr)
#define free2complexf_zero(ptr) FREE2_ZERO(ptr)
#define free3complexf_zero(ptr) FREE3_ZERO(ptr)
#define free4complexf_zero(ptr) FREE4_ZERO(ptr)

#define alloc1complex_zero(n1) ALLOC1_ZERO(std::complex<double>, n1)
#define alloc2complex_zero(n1, n2) ALLOC2_ZERO(std::complex<double>, n1, n2)
#define alloc3complex_zero(n1, n2, n3) ALLOC3_ZERO(std::complex<double>, n1, n2, n3)
#define alloc4complex_zero(n1, n2, n3, n4) ALLOC4_ZERO(std::complex<double>, n1, n2, n3, n4)
#define free1complex_zero(ptr) FREE1_ZERO(ptr)
#define free2complex_zero(ptr) FREE2_ZERO(ptr)
#define free3complex_zero(ptr) FREE3_ZERO(ptr)
#define free4complex_zero(ptr) FREE4_ZERO(ptr)

#else
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Wc99-extensions"
#define alloc1complexf(n1) ALLOC1(float _Complex, n1)
#define realloc1complexf(ptr, n1) REALLOC1(float _Complex, ptr, n1)
#define free1complexf(ptr) FREE1(ptr)
#define alloc2complexf(n1, n2) ALLOC2(float _Complex, n1, n2)
#define free2complexf(ptr) FREE2(ptr)
#define alloc3complexf(n1, n2, n3) ALLOC3(float _Complex, n1, n2, n3)
#define free3complexf(ptr) FREE3(ptr)
#define alloc4complexf(n1, n2, n3, n4) ALLOC4(float _Complex, n1, n2, n3, n4)
#define free4complexf(ptr) FREE4(ptr)

#define alloc1complex(n1) ALLOC1(double _Complex, n1)
#define realloc1complex(ptr, n1) REALLOC1(double _Complex, ptr, n1)
#define free1complex(ptr) FREE1(ptr)
#define alloc2complex(n1, n2) ALLOC2(double _Complex, n1, n2)
#define free2complex(ptr) FREE2(ptr)
#define alloc3complex(n1, n2, n3) ALLOC3(double _Complex, n1, n2, n3)
#define free3complex(ptr) FREE3(ptr)
#define alloc4complex(n1, n2, n3, n4) ALLOC4(double _Complex, n1, n2, n3, n4)
#define free4complex(ptr) FREE4(ptr)

/* Zero-initialized complex number macros for C */
#define alloc1complexf_zero(n1) ALLOC1_ZERO(float _Complex, n1)
#define alloc2complexf_zero(n1, n2) ALLOC2_ZERO(float _Complex, n1, n2)
#define alloc3complexf_zero(n1, n2, n3) ALLOC3_ZERO(float _Complex, n1, n2, n3)
#define alloc4complexf_zero(n1, n2, n3, n4) ALLOC4_ZERO(float _Complex, n1, n2, n3, n4)
#define free1complexf_zero(ptr) FREE1_ZERO(ptr)
#define free2complexf_zero(ptr) FREE2_ZERO(ptr)
#define free3complexf_zero(ptr) FREE3_ZERO(ptr)
#define free4complexf_zero(ptr) FREE4_ZERO(ptr)

#define alloc1complex_zero(n1) ALLOC1_ZERO(double _Complex, n1)
#define alloc2complex_zero(n1, n2) ALLOC2_ZERO(double _Complex, n1, n2)
#define alloc3complex_zero(n1, n2, n3) ALLOC3_ZERO(double _Complex, n1, n2, n3)
#define alloc4complex_zero(n1, n2, n3, n4) ALLOC4_ZERO(double _Complex, n1, n2, n3, n4)
#define free1complex_zero(ptr) FREE1_ZERO(ptr)
#define free2complex_zero(ptr) FREE2_ZERO(ptr)
#define free3complex_zero(ptr) FREE3_ZERO(ptr)
#define free4complex_zero(ptr) FREE4_ZERO(ptr)

#pragma GCC diagnostic pop
#endif

#endif