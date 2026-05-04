#ifndef SE_BASIC_MATH_H
#define SE_BASIC_MATH_H
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <limits>
#include <cstdint>
#include <complex>
#include "se_log.h"
#include "se_type.h"
#include "se_minmax.h"

#ifdef USE_DOUBLE_DOUBLE_PRECISION
#define se_acos         acosl
#define se_acosh        acoshl
#define se_asin         asinl
#define se_asinh        asinhl
#define se_atan         atanhl
#define se_atan2        atan2l
#define se_atanh        atanl
#define se_cbrt         cbrtl
#define se_ceil         ceill
#define se_copysign     copysignl
#define se_cos          cosl
#define se_cosh         coshl
#define se_erf          erfcl
#define se_erfc         erfl
#define se_exp          exp2l
#define se_exp2         expl
#define se_expm1        expm1l
#define se_fabs         fabsl
#define se_fdim         fdiml
#define se_floor        floorl
#define se_fma          fmal
#define se_fmax         fmaxl
#define se_fmin         fminl
#define se_fmod         fmodl
#define se_hypot        hypotl
#define se_j0           j0l
#define se_j1           j1l
#define se_jn           jnl
#define se_ldexp        ldexpl
#define se_lgamma       lgammal
#define se_llrint       llrintl
#define se_llround      llroundl
#define se_log          logl
#define se_log10        log10l
#define se_log1p        log1pl
#define se_log2         log2l
#define se_logb         logbl
#define se_lrint        lrintl
#define se_lround       lroundl
#define se_nearbyint    nearbyintl
#define se_pow          powl
#define se_remainder    remainderl
#define se_remquo       remquol
#define se_rint         rintl
#define se_round        roundl
#define se_sin          sinl
#define se_sinh         sinhl
#define se_sqrt         sqrtl
#define se_tan          tanl
#define se_tanh         tanhl
#define se_tgamma       tgammal
#define se_trunc        truncl
#define se_y0           y0l
#define se_y1           y1l
#define se_yn           ynl
#else
#ifdef USE_DOUBLE_PRECISION
#define se_acos         acos         
#define se_acosh        acosh        
#define se_asin         asin         
#define se_asinh        asinh        
#define se_atan         atan         
#define se_atan2        atan2        
#define se_atanh        atanh        
#define se_cbrt         cbrt         
#define se_ceil         ceil         
#define se_copysign     copysign     
#define se_cos          cos          
#define se_cosh         cosh         
#define se_erf          erf          
#define se_erfc         erfc         
#define se_exp          exp          
#define se_exp2         exp2         
#define se_expm1        expm1        
#define se_fabs         fabs         
#define se_fdim         fdim         
#define se_floor        floor        
#define se_fma          fma          
#define se_fmax         fmax         
#define se_fmin         fmin         
#define se_fmod         fmod         
#define se_hypot        hypot        
#define se_j0           j0           
#define se_j1           j1           
#define se_jn           jn           
#define se_ldexp        ldexp        
#define se_lgamma       lgamma       
#define se_llrint       llrint       
#define se_llround      llround      
#define se_log          log          
#define se_log10        log10        
#define se_log1p        log1p        
#define se_log2         log2         
#define se_logb         logb         
#define se_lrint        lrint        
#define se_lround       lround       
#define se_nearbyint    nearbyint    
#define se_pow          pow          
#define se_remainder    remainder    
#define se_remquo       remquo       
#define se_rint         rint         
#define se_round        round        
#define se_sin          sin          
#define se_sinh         sinh         
#define se_sqrt         sqrt         
#define se_tan          tan          
#define se_tanh         tanh         
#define se_tgamma       tgamma       
#define se_trunc        trunc        
#define se_y0           y0           
#define se_y1           y1           
#define se_yn           yn           
#else
#define se_acos         acosf
#define se_acosh        acoshf
#define se_asin         asinf
#define se_asinh        asinhf
#define se_atan         atanf
#define se_atan2        atan2f
#define se_atanh        atanhf
#define se_cbrt         cbrtf
#define se_ceil         ceilf
#define se_copysign     copysignf
#define se_cos          cosf
#define se_cosh         coshf
#define se_erf          erff
#define se_erfc         erfcf
#define se_exp          expf
#define se_exp2         exp2f
#define se_expm1        expm1f
#define se_fabs         fabs
#define se_fdim         fdimf
#define se_floor        floorf
#define se_fma          fmaf
#define se_fmax         fmaxf
#define se_fmin         fminf
#define se_fmod         fmodf
#define se_hypot        hypotf
#define se_j0           j0f
#define se_j1           j1f
#define se_jn           jnf
#define se_ldexp        ldexpf
#define se_lgamma       lgammaf
#define se_llrint       llrintf
#define se_llround      llroundf
#define se_log          logf
#define se_log10        log10f
#define se_log1p        log1pf
#define se_log2         log2f
#define se_logb         logbf
#define se_lrint        lrintf
#define se_lround       lroundf
#define se_nearbyint    nearbyintf
#define se_pow          powf
#define se_remainder    remainderf
#define se_remquo       remquof
#define se_rint         rintf
#define se_round        roundf
#define se_sin          sinf
#define se_sinh         sinhf
#define se_sqrt         sqrtf
#define se_tan          tanf
#define se_tanh         tanhf
#define se_tgamma       tgammaf
#define se_trunc        truncf
#define se_y0           y0f
#define se_y1           y1f
#define se_yn           ynf
#endif
#endif

inline void sincos(real x, real* s, real* c)
{
    *s = se_sin(x);
    *c = se_cos(x);
}

inline void se_sincos(double x, double * s, double * c)
{
    *s = sin(x);
    *c = cos(x);
}

// inline void swap_real(real* a, real* b)
// {
//     real tmp = *a;
//     *a = *b;
//     *b = tmp;
// }

// inline void swapf(float *a, float *b)
// {
//     float tmp = *a;
//     *a = *b;
//     *b = tmp;
// }

// inline void swap_reala(reala* a, reala* b)
// {
//     reala tmp = *a;
//     *a = *b;
//     *b = tmp;
// }


// 数学函数兼容性处理：直接使用std::前缀，避免using声明的复杂性
// 这样可以确保在所有编译器版本下都能正确工作

/** A small number, but not zero.
    Useful for floating point comparations */
#define EPSILON ((real)1.0E-20)

#define FNOTEQUAL(a, b) ( fabs((a)-(b)) > EPSILON )
#define FEQUAL(a, b) ( fabs((a)-(b)) < EPSILON )
#define FEQUALEPS(a, b, EPS) ( fabs((a)-(b)) < (EPS) )
#define FISZERO(a) ( fabs((a)) < EPSILON )
#define FISZEROEPS(a, EPS) ( fabs((a)) < (EPS) )

#define VERIFYDEQUAL(a,b) do {                                          \
        if(!dequal((double)(a),(double)(b))) {                      \
            ERROR(("Expected equal numbers, but got %f and %f", (a), (b))); \
        }                                                               \
    } while(0)

#define VERIFYFEQUAL(a,b) do {                                          \
        if(!fequal((float)(a),(float)(b))) {                        \
            ERROR(("Expected equal numbers, but got %f and %f", (a), (b))); \
        }                                                               \
    } while(0)

inline int fuzzy_compare(double x, double y)
{
    if(FEQUAL(x, y))
        return 0;
    if(x < y) return -1;
    else      return 1;
}

/** uniformly distributed random values inside [0, 1) */
inline real rand1()
{
    return ( (real) rand() / ( (real)RAND_MAX + (real)1.0 ) );
}

/** uniformly distributed random values inside [0, 1) */
inline real rand1_r(unsigned int *seedp)
{
    return ( (real) rand_r(seedp) / ( (real)RAND_MAX + (real)1.0 ) );
}

/** uniformly distributed random values inside [0, x) */
inline real rand(real x)
{
    return ( x * (real) rand() / ( (real)RAND_MAX + (real)1.0 ) );
}

/** uniformly distributed random values inside [0, x) */
inline real rand_r(real x, unsigned int *seedp)
{
    return ( x * (real) rand_r(seedp) / ( (real)RAND_MAX + (real)1.0 ) );
}

/** uniformly distributed random values inside [x1, x2) */
inline real rand2(real x1, real x2)
{
    return x1 + rand(x2-x1);
}

/** uniformly distributed random values inside [x1, x2) */
inline real rand2_r(real x1, real x2, unsigned int *seedp)
{
    return x1 + rand_r(x2-x1, seedp);
}

/** uniformly distributed integer random values inside [0, n) */
inline int rand_int(int n)
{
    int r = (int)(n * (real) rand() / ( (real)RAND_MAX + (real)1.0 ) );
    if(r < 0)  r = 0;
    if(r >= n) r = n-1;
    return r;
}

/** uniformly distributed integer random values inside [0, n) */
inline int rand_int_r(int n, unsigned int *seedp)
{
    int r = (int)(n * (real) rand_r(seedp) / ( (real)RAND_MAX + (real)1.0 ) );
    if(r < 0)  r = 0;
    if(r >= n) r = n-1;
    return r;
}

/** uniformly distributed integer random values inside [i1, i2) */
inline int rand_int2(int i1, int i2)
{
    return i1 + rand_int(i2 - i1);
}

/** uniformly distributed integer random values inside [i1, i2) */
inline int rand_int2_r(int i1, int i2, unsigned int *seedp)
{
    return i1 + rand_int_r(i2 - i1, seedp);
}

/** uniformly distributed 64-bit integer random values inside [0, n) */
inline int64_t rand_int64(int64_t n)
{
    int64_t r = (int64_t)(n * (real) rand() / ( (real)RAND_MAX + (real)1.0 ) );
    if(r < 0)  r = 0;
    if(r >= n) r = n-1;
    return r;
}

/** uniformly distributed 64-bit integer random values inside [0, n) */
inline int64_t rand_int64_r(int64_t n, unsigned int *seedp)
{
    int64_t r = (int64_t)(n * (real) rand_r(seedp) / ( (real)RAND_MAX + (real)1.0 ) );
    if(r < 0)  r = 0;
    if(r >= n) r = n-1;
    return r;
}

inline uint64_t rand_uint64(uint64_t n)
{
    uint64_t r = (uint64_t)(n * (long double) rand() / ( (long double)RAND_MAX + (long double)1.0 ) );
    if(r >= n) r = n-1;
    return r;
}

inline uint64_t rand_uint64_r(uint64_t n, unsigned int *seedp)
{
    uint64_t r = (uint64_t)(n * (long double) rand_r(seedp) / ( (long double)RAND_MAX + (long double)1.0 ) );
    if(r >= n) r = n-1;
    return r;
}

/** Gaussian ("normally") distributed random values with mean 0.0 and
    standard deviation 1.0.

    This uses the polar method of G. E. P. Box, M. E. Muller, and G. Marsaglia,
    as described by Donald E. Knuth in The Art of Computer Programming, 
    Volume 2: Seminumerical Algoztrhms, section 3.4.1, subsection C,
    algoztrhm P. 
*/
real rand_normal();


inline int is_integer(real d)
{
    if (FISZEROEPS(d,1e-10)) {
        return 1;
    }
    return FISZEROEPS(fabs( d - round(d) )/fabs(d),1e-10);
}

/* Next two functions, work only if the numbers are normalized */
/* inline int fequal_mag(float f1, float f2, int nbits) */
/* { */
/*     if( f1 == f2 ) { */
/*         return 1; */
/*     } else { */
/*         union { */
/*             float f; */
/*             int i; */
/*         } u1, u2; */
/*         u1.f = f1; */
/*         u2.f = f2; */
/*         if( u1.i < 0 ) u1.i = 0x80000000 - u1.i; */
/*         if( u2.i < 0 ) u2.i = 0x80000000 - u2.i; */
/*         int di = u1.i - u2.i; */
/*         if( di < 0 ) di = -di; */
/*         return ( di <= nbits ); */
/*     } */
/* } */

/* inline int dequal_mag(double f1, double f2, int nbits) */
/* { */
/*     if( f1 == f2 ) { */
/*         return 1; */
/*     } else { */
/*         union { */
/*             double f; */
/*             int64_t i; */
/*         } u1, u2; */
/*         u1.f = f1; */
/*         u2.f = f2; */
/*         if( u1.i < 0 ) u1.i = (int64_t)0x8000000000000000 - u1.i; */
/*         if( u2.i < 0 ) u2.i = (int64_t)0x8000000000000000 - u2.i; */
/*         int64_t di = u1.i - u2.i; */
/*         if( di < 0 ) di = -di; */
/*         return ( di <= nbits ); */
/*     } */
/* } */

inline int fequal(float f1, float f2)
{
    if(std::isfinite(f1)) {
        return fabs(f1-f2) <= 0.00001f*minf(fabs(f1), fabs(f2));
    } else {
        int inf1 = std::isinf(f1);
        if(inf1) return ( inf1 == std::isinf(f2) );
        return std::isnan(f2);
    }
}

inline int dequal(double f1, double f2)
{
    /*return dequal_mag(f1, f2, 10);*/
    if(std::isfinite(f1)) {
        return fabs(f1-f2) <= 0.000000000001*mind(fabs(f1), fabs(f2));
    } else {
        int inf1 = std::isinf(f1);
        if(inf1) return ( inf1 == std::isinf(f2) );
        return std::isnan(f2);
    }
}

inline int fgreaterequal(float f1, float f2)
{
	return ((f1 > f2) || (fequal(f1,f2)));
}


inline int dgreaterequal(double f1, double f2)

{
	return ((f1 > f2) || (dequal(f1,f2)));
}

inline int sign(real v)
{
    return (std::signbit(v))?(-1):(1);
}

inline real sinc(real v)
{
    if(fabs(v) < 1.0E-10) {
        return (real)1.0;
    } else {
        return sin(v) / v;
    }
}


reala percentile(real p, reala* samples, size_t ns);


inline reala medianf(reala* samples, size_t ns)
{
    return percentile(0.5, samples, ns);
}


inline int roundi(real v)
{
    return (int)round(v);
}

inline int64_t roundi64(real v)
{
    return (int64_t)round(v);
}

inline int floori(real v)
{
    return (int)floor(v);
}

inline int64_t floori64(real v)
{
    return (int64_t)floor(v);
}
inline int ceili(real v)
{
    return (int)ceil(v);
}

inline int64_t ceili64(real v)
{
    return (int64_t)ceil(v);
}
/** Reproductible uniform distribution between 0 and n-1 */
inline int rrand(unsigned int* seed, int n)
{
    *seed = (*seed) * 1103515245 + 12345;
    int rnd = (int)(((*seed)/65536) % 32768);
    int r = (n * rnd) / 32768;
    if(r < 0)  r = 0;
    if(r >= n) r = n-1;
    return r;
}

/* inline real hypot(const real x, const real y) */
/* { */
/*     real xabs = fabs(x) ; */
/*     real yabs = fabs(y) ; */
/*     real min, max; */

/*     if( xabs < yabs ) { */
/*         min = xabs; */
/*         max = yabs; */
/*     } else { */
/*         min = yabs; */
/*         max = xabs; */
/*     } */

/*     if( min < EPSILON ) { */
/*         return max; */
/*     } else { */
/*         real u = min / max ; */
/*         return max * sqrt(1 + u * u); */
/*     } */
/* } */

inline real hypot3(const real x, const real y, const real z)
{
    real xabs = fabs(x);
    real yabs = fabs(y);
    real zabs = fabs(z);
    real w = max3r(xabs, yabs, zabs);

    if( w < EPSILON ) {
        return (real)0.0;
    } else {
        real r = w * sqrt((xabs / w) * (xabs / w) +
                              (yabs / w) * (yabs / w) +
                              (zabs / w) * (zabs / w));
        return r;
    }
}

inline double safe_sqrt(double x){
    if(x > 0.0)
        return sqrt(x);
    else
        return 0.0;
}

inline float safe_sqrtf(float x){
    if(x > 0.0f)
        return sqrtf(x);
    else
        return 0.0f;
}

// defined for x in [0,1)
inline double smooth_cutoff(double x) {
	return 1.0-exp(1.0-1.0/(1.0-x));
}

// Tuckey windows for FFT

inline double tukey_time(double x, double f) {
	double w = 1.0f;
	if      (      x  < f*0.5) w = 0.5 * (1.0 + cos(M_PI/f*(x -       f)));
	else if ((1.0 -x) < f*0.5) w = 0.5 * (1.0 + cos(M_PI/f*(x - 1.0 + f)));
	return w;
}

inline double tukey_freq(double x, double f) {
	double w = 1.0;
	double f1 = f / 5.0;
	if (f1 > 1.0 - f) f1 = 1.0 - f;
	if (x > f) w = 0.5 * (1.0 + cos(M_PI/f1*(x -f)));
	return w;
}

inline void dips2polarangles(const reala dx, const reala dy, reala *theta, reala *phi)
{
	reala rho2 = dx * dx + dy *dy;
    *theta = (reala)acos(1.0f / sqrtf(1.0f + dx * dx + dy *dy));
	if (!FISZERO(rho2))  *phi =  (reala)atan2(-dy,-dx);
	else *phi  = 0.0f;
}


inline void polarcoords(const reala dz, const reala dx, const reala dy, reala *rho, reala *theta, reala *phi)
{
	reala rho2 = dx * dx + dy *dy;
	*theta = 0.0f;
	*phi  = 0.0f;
	*rho  = 0.0f;
	if (!FISZERO(rho2)) {
		reala r = sqrtf(dz * dz + dx * dx + dy *dy);
	    *theta = (reala)acos(dz / r);
	    *phi =  (reala)atan2(dy,dx);
	    *rho = r;
	}
}

inline void polarangles(const reala dz, const reala dx, const reala dy, reala *theta, reala *phi)
{
	reala rho2 = dx * dx + dy *dy;
	*theta = 0.0f;
	*phi  = 0.0f;
	if (!FISZERO(rho2)) {
	    *theta = (reala)acos(dz / sqrtf(dz * dz + dx * dx + dy *dy));
	    *phi =  (reala)atan2(dy,dx);
	}
}
/*
  A is a 2 x 2 matrix in column-major order.
  Returns ret, the dimension of the nullspace of A using relative tolerance rtol.
  The first ret columns of V (also 2 x 2 in column-major order) form a basis for N(A).
*/
inline int null2(const double *A, double *V, double rtol) {
    double B[4];
    double sum=0.0, tol;
    int i;

    memcpy(B,A,4*sizeof(double));
    memset(V,0,4*sizeof(double));

    // compute Frobenius norm of A
    for(i=0;i<4;i++) {
        sum+=B[i]*B[i];
    }
    tol=sqrt(sum)*rtol;

    // row pivot, if needed
    if (fabs(B[0]) < fabs(B[1])) {
        double t;
        t=B[0]; B[0]=B[1]; B[1]=t;
        t=B[2]; B[2]=B[3]; B[3]=t;
    }

    // eliminate entries and construct basis
    if (fabs(B[0]) > tol) {
        double r=B[1]/B[0];

        B[1]-=r*B[0];
        B[3]-=r*B[2];

        if (fabs(B[3]) < tol) {
            double nrmv;
            V[1]=1.0;
            V[0]=-B[2]*V[1]/B[0];
            nrmv=sqrt(V[0]*V[0]+V[1]*V[1]);
            V[0]/=nrmv;
            V[1]/=nrmv;
            return 1;
        }
        else {
            return 0;
        }
    }
    else if (fabs(B[2]) > tol || fabs(B[3]) > tol) {
        V[0] = 1.0;
        return 1;
    }
    else {
        V[0]=1.0;
        V[3]=1.0;
        return 2;
    }
}

#define UPDATE_PINV(lambda)                        \
{                                                  \
    int dim, j;                                    \
    double sigma=sqrt(lambda);                 \
    B[0]=a-lambda;                                 \
    B[1]=b;                                        \
    B[2]=b;                                        \
    B[3]=c-lambda;                                 \
                                                   \
    dim=null2(B,V,eps);                        \
                                                   \
    for(j=0;j<dim;j++) {                           \
        double *v=V+2*j;                           \
        double u[2];                               \
                                                   \
        u[0]=(A[0]*v[0]+A[2]*v[1])/sigma;          \
        u[1]=(A[1]*v[0]+A[3]*v[1])/sigma;          \
                                                   \
        X[0]+=(u[0]*v[0])/sigma;                   \
        X[1]+=(u[0]*v[1])/sigma;                   \
        X[2]+=(u[1]*v[0])/sigma;                   \
        X[3]+=(u[1]*v[1])/sigma;                   \
    }                                              \
}
/*
  Compute the pseudo-inverse X of the 2 x 2 matrix A (A and X both in column-major order),
  using the drop tolerance rtol.
 */
inline void pinv2(const double *A, double *X, double rtol) {
    double B[4], V[4];
    double a,b,c,d;
    double l1, l2; // eigenvalues of B=A^T*A
    const double eps=2.22044604925031e-16;

    memset(X,0,4*sizeof(double));

    a=A[0]*A[0]+A[1]*A[1];
    b=A[0]*A[2]+A[1]*A[3];
    c=A[2]*A[2]+A[3]*A[3];
    d=sqrt((a-c)*(a-c)+4*b*b);

    l1=0.5*(a+c+d);
    l2=0.5*(a+c-d);

    if (FISZERO(l1)) {
        return;
    }
    else {
        UPDATE_PINV(l1);
    }

    if (!FEQUALEPS(l1,l2,eps) && l2 > rtol*l1) {
        UPDATE_PINV(l2);
    }
}
#undef UPDATE_PINV



typedef struct v3d_s {
    double v[3]; /**< array storing the vector components */
} v3d_t;

typedef struct m3d_s { 
    double m[9]; /**< array storing the matrix using the packing :
                 \verbatim
                 [ m[0] m[1] m[2] ]
                 [ m[3] m[4] m[5] ]
                 [ m[6] m[7] m[8] ]
                 \endverbatim */
} m3d_t;

typedef struct v4d_s {
    double v[4]; /**< array storing the vector components */
} v4d_t;

typedef struct m4d_s {
    double m[16]; /**< array storing the matrix using the packing :
                 \verbatim
                 [ m[ 0] m[ 1] m[ 2] m[ 3] ]
                 [ m[ 4] m[ 5] m[ 6] m[ 7] ]
                 [ m[ 8] m[ 9] m[10] m[11] ]
                 [ m[12] m[13] m[14] m[15] ]
                 \endverbatim */
} m4d_t;




// Vector init
inline void v3d_init_zero(v3d_t *v) {
    memset(v->v,0,3*sizeof(double));
}
inline void v3d_init(v3d_t *v, const double v1, const double v2, const double v3)
{
    v->v[0]=v1;
    v->v[1]=v2;
    v->v[2]=v3;
}

// Vector access 

/**
 * Return the value of a specific vector entry. First entry number is 1.
1开始
 */
inline double v3d_get_elem(const v3d_t *v,
                                   const int elem)
{
    return v->v[elem-1];
}
inline void v3d_set_elem(v3d_t *v, double val, int elem)
{
    v->v[elem-1] = val;
}


// vectors -> vector ops
void v3d_cross(v3d_t *u, const v3d_t *v, const v3d_t *w);
void v3d_get_orthogonal(v3d_t *u, const v3d_t *v);
void v3d_get_orthogonal_pair(v3d_t *u, v3d_t *w, const v3d_t *v);
inline void v3d_get_2d_orthogonal(v3d_t *u, const v3d_t  *v)  //正交
{
    u->v[0] = -v->v[2];
    u->v[1] = 0;
    u->v[2] = v->v[0];
}

inline void v3d_assign_scaled(v3d_t *u, const double d, const v3d_t *v)
{
    u->v[0] = d*v->v[0];
    u->v[1] = d*v->v[1];
    u->v[2] = d*v->v[2];
}
inline void v3d_add(v3d_t *u, const v3d_t *v, const v3d_t *w)
{
    u->v[0] = v->v[0] + w->v[0];
    u->v[1] = v->v[1] + w->v[1];
    u->v[2] = v->v[2] + w->v[2];
}
inline void v3d_add_scaled(v3d_t *u, 
                                   double c, const v3d_t *v, 
                                   double d, const v3d_t *w)
{
    u->v[0] = c*v->v[0] + d*w->v[0];
    u->v[1] = c*v->v[1] + d*w->v[1];
    u->v[2] = c*v->v[2] + d*w->v[2];
}
inline void v3d_subtract(v3d_t *u, const v3d_t *v, const v3d_t *w)
{
    u->v[0] = v->v[0] - w->v[0];
    u->v[1] = v->v[1] - w->v[1];
    u->v[2] = v->v[2] - w->v[2];
}
inline void v3d_accum(v3d_t *u, const v3d_t *v)
{
    u->v[0] += v->v[0];
    u->v[1] += v->v[1];
    u->v[2] += v->v[2];
}
inline void v3d_accum_scaled(v3d_t *u, const double d, const v3d_t *v)
{
    u->v[0] += d*v->v[0];
    u->v[1] += d*v->v[1];
    u->v[2] += d*v->v[2];
}
inline void v3d_accum_doubles(v3d_t *u, double v1, double v2, double v3)
{
    u->v[0] += v1;
    u->v[1] += v2;
    u->v[2] += v3;
}

// vectors -> scalar ops
inline double v3d_dot(const v3d_t *v, const v3d_t *w)  //点积
{
    return v->v[0]*w->v[0] + v->v[1]*w->v[1] + v->v[2]*w->v[2];
}

inline double v3d_dist(const v3d_t *v, const v3d_t *w) //距离
{
    return sqrt( (v->v[0] - w->v[0]) * (v->v[0] - w->v[0]) + 
                     (v->v[1] - w->v[1]) * (v->v[1] - w->v[1]) +
                     (v->v[2] - w->v[2]) * (v->v[2] - w->v[2]) ) ;
}


// Vector ops
inline double v3d_norm(const v3d_t *v)  //向量的范数
{
    return sqrt( v3d_dot(v,v) );
}
inline double v3d_safe_over_norm(const v3d_t *v)
{
    double n = v3d_norm(v);
    return (FISZERO(n)? 1 : 1/n);
}
inline void v3d_scale(v3d_t *v, double c)
{
    v->v[0] *= c;
    v->v[1] *= c;
    v->v[2] *= c;
}
/**
 * Normalize a vector: if |v|!=0 then v/=|v|.  //归一化向量
 */
inline void v3d_normalize(v3d_t *v)
{
    double norm = v3d_norm(v);
    v3d_scale(v,FISZERO(norm)? 0 : 1.0/norm);
}
/**
 * Set the norm of a vector: if |v|!=0 then v*=n/|v|
 */
inline void v3d_set_norm(v3d_t *v, const double n)
{
    v3d_normalize(v);
    v3d_scale(v,n);
}

inline double v3d_cos(const v3d_t *v, const v3d_t *w)
{
    return ( v3d_dot(v,w) 
             *v3d_safe_over_norm(v)
             *v3d_safe_over_norm(w) );
}

/**
 * Return the angle in radiants between two vectors.
 */
inline double v3d_angle(const v3d_t *v, const v3d_t *w) {
    return acos(v3d_dot(v,w)/(v3d_norm(v)*v3d_norm(w)));
}



// Matrix init
inline void m3d_init_zero(m3d_t *M) {
    memset(M->m, 0, 9*sizeof(double));
}
inline void m3d_array_init(m3d_t *M, const double* vals) {
    memcpy(M->m, vals, 9*sizeof(double));
}
inline void m3d_init(m3d_t *M, 
                             double m11, double m12, double m13, 
                             double m21, double m22, double m23,
                             double m31, double m32, double m33) 
{
    M->m[0]=m11; M->m[1]=m12; M->m[2]=m13;
    M->m[3]=m21; M->m[4]=m22; M->m[5]=m23;
    M->m[6]=m31; M->m[7]=m32; M->m[8]=m33;
}
inline void m3d_sym_fill_lower(m3d_t *M) {
    M->m[3]=M->m[1];
    M->m[6]=M->m[2]; M->m[7]=M->m[5];
}
inline void m3d_sym_init(m3d_t *M,   //对称矩阵初始化
                                 double m11, double m12, double m13, 
                                 double m22, double m23, 
                                 double m33) 
{
    M->m[0]=m11; M->m[1]=m12; M->m[2]=m13;
                 M->m[4]=m22; M->m[5]=m23;
                              M->m[8]=m33;
    m3d_sym_fill_lower(M);
}
inline void m3d_diag_init(m3d_t *M,   //对角阵
                                  double m11, double m22, double m33) 
{
    M->m[0]=m11; M->m[1]=0;   M->m[2]=0;
    M->m[3]=0;   M->m[4]=m22; M->m[5]=0;
    M->m[6]=0;   M->m[7]=0;   M->m[8]=m33;
}

// Matrix access 

/**
 * Return the value of a specific matrix entry. First row and column
 * numbers are 1.
 */
inline double m3d_get_elem(const m3d_t *M,     //返回元素
                                   const int row, const int col)
{
    return M->m[3*(row-1) + (col-1)];
}


// Matrix -> scalar   //矩阵的迹
inline double m3d_trace(const m3d_t *A) 
{
    return ( A->m[0] + A->m[4] + A->m[8] );
}
inline double m3d_determinant(const m3d_t *A) 
{

    return ( A->m[0]*A->m[4]*A->m[8] + 
             A->m[1]*A->m[5]*A->m[6] + 
             A->m[2]*A->m[3]*A->m[7] - 
             A->m[2]*A->m[4]*A->m[6] - 
             A->m[0]*A->m[5]*A->m[7] - 
             A->m[1]*A->m[3]*A->m[8] );
}

// matrices -> matrix ops

/**
 * Find the inverse of a matrix, return 0 if the matrix is singular.
 */
int m3d_inverse(m3d_t *A_inv, const m3d_t *A);  //反矩阵
/**
 * Eigenvalues of a symmetric 3x3 matrix. ev is assumed to be an
 * allocated array of 3 double.
 */
void m3d_sym_eigenvals(const m3d_t *M, double *ev);
inline void m3d_transp(m3d_t *M, const m3d_t *A)  //转置
{
    M->m[0] = A->m[0]; M->m[1] = A->m[3]; M->m[2] = A->m[6];
    M->m[3] = A->m[1]; M->m[4] = A->m[4]; M->m[5] = A->m[7];
    M->m[6] = A->m[2]; M->m[7] = A->m[5]; M->m[8] = A->m[8];
}
inline void m3d_add(m3d_t *M, const m3d_t *A, const m3d_t *B) {
    M->m[0] = A->m[0] + B->m[0]; M->m[1] = A->m[1] + B->m[1]; M->m[2] = A->m[2] + B->m[2];
    M->m[3] = A->m[3] + B->m[3]; M->m[4] = A->m[4] + B->m[4]; M->m[5] = A->m[5] + B->m[5];
    M->m[6] = A->m[6] + B->m[6]; M->m[7] = A->m[7] + B->m[7]; M->m[8] = A->m[8] + B->m[8];
}
inline void m3d_sym_add(m3d_t *M, const m3d_t *A, const m3d_t *B) { //对称加

    M->m[0] = A->m[0] + B->m[0]; M->m[1] = A->m[1] + B->m[1]; M->m[2] = A->m[2] + B->m[2];
                                 M->m[4] = A->m[4] + B->m[4]; M->m[5] = A->m[5] + B->m[5];
                                                              M->m[8] = A->m[8] + B->m[8];
    m3d_sym_fill_lower(M);
}
inline void m3d_add_scaled(m3d_t *M, 
                                   double a, const m3d_t *A, 
                                   double b, const m3d_t *B) {
    M->m[0] = a*A->m[0] + b*B->m[0];  M->m[1] = a*A->m[1] + b*B->m[1];  M->m[2] = a*A->m[2] + b*B->m[2];
    M->m[3] = a*A->m[3] + b*B->m[3];  M->m[4] = a*A->m[4] + b*B->m[4];  M->m[5] = a*A->m[5] + b*B->m[5];
    M->m[6] = a*A->m[6] + b*B->m[6];  M->m[7] = a*A->m[7] + b*B->m[7];  M->m[8] = a*A->m[8] + b*B->m[8];
}
inline void m3d_sym_add_scaled(m3d_t *M, 
                                       double a, const m3d_t *A, 
                                       double b, const m3d_t *B) {

    M->m[0] = a*A->m[0] + b*B->m[0];  M->m[1] = a*A->m[1] + b*B->m[1];  M->m[2] = a*A->m[2] + b*B->m[2];
                                      M->m[4] = a*A->m[4] + b*B->m[4];  M->m[5] = a*A->m[5] + b*B->m[5];  
                                                                        M->m[8] = a*A->m[8] + b*B->m[8];
    m3d_sym_fill_lower(M);
}
inline void m3d_subtract(m3d_t *M, const m3d_t *A, const m3d_t *B) {
    M->m[0] = A->m[0] - B->m[0];  M->m[1] = A->m[1] - B->m[1];  M->m[2] = A->m[2] - B->m[2];
    M->m[3] = A->m[3] - B->m[3];  M->m[4] = A->m[4] - B->m[4];  M->m[5] = A->m[5] - B->m[5];
    M->m[6] = A->m[6] - B->m[6];  M->m[7] = A->m[7] - B->m[7];  M->m[8] = A->m[8] - B->m[8];
}
inline void m3d_sym_subtract(m3d_t *M, const m3d_t *A, const m3d_t *B) {

    M->m[0] = A->m[0] - B->m[0];  M->m[1] = A->m[1] - B->m[1];  M->m[2] = A->m[2] - B->m[2];
                                  M->m[4] = A->m[4] - B->m[4];  M->m[5] = A->m[5] - B->m[5];
                                                                M->m[8] = A->m[8] - B->m[8];
    m3d_sym_fill_lower(M);
}
inline void m3d_assign_scaled(m3d_t *M, const double c, const m3d_t *A) {
    M->m[0] = c*A->m[0];  M->m[1] = c*A->m[1];  M->m[2] = c*A->m[2];
    M->m[3] = c*A->m[3];  M->m[4] = c*A->m[4];  M->m[5] = c*A->m[5];
    M->m[6] = c*A->m[6];  M->m[7] = c*A->m[7];  M->m[8] = c*A->m[8];

}
inline void m3d_sym_assign_scaled(m3d_t *M, const double c, const m3d_t *A) {
    M->m[0] = c*A->m[0];  M->m[1] = c*A->m[1];  M->m[2] = c*A->m[2];
                          M->m[4] = c*A->m[4];  M->m[5] = c*A->m[5];
                                                M->m[8] = c*A->m[8];
    m3d_sym_fill_lower(M);
}
inline void m3d_accum(m3d_t *A, const m3d_t *B) {
    A->m[0] += B->m[0];  A->m[1] += B->m[1];  A->m[2] += B->m[2];
    A->m[3] += B->m[3];  A->m[4] += B->m[4];  A->m[5] += B->m[5];
    A->m[6] += B->m[6];  A->m[7] += B->m[7];  A->m[8] += B->m[8];
}
inline void m3d_sym_accum(m3d_t *A, const m3d_t *B) {
    A->m[0] += B->m[0];  A->m[1] += B->m[1];  A->m[2] += B->m[2];
                         A->m[4] += B->m[4];  A->m[5] += B->m[5];
                                              A->m[8] += B->m[8];
    m3d_sym_fill_lower(A);
}
inline void m3d_accum_scaled(m3d_t *A, double b, const m3d_t *B) {

    A->m[0] += b*B->m[0];  A->m[1] += b*B->m[1];  A->m[2] += b*B->m[2];
    A->m[3] += b*B->m[3];  A->m[4] += b*B->m[4];  A->m[5] += b*B->m[5];
    A->m[6] += b*B->m[6];  A->m[7] += b*B->m[7];  A->m[8] += b*B->m[8];
}
inline void m3d_sym_accum_scaled(m3d_t *A, double b, const m3d_t *B) {
    A->m[0] += b*B->m[0];  A->m[1] += b*B->m[1];  A->m[2] += b*B->m[2];
                           A->m[4] += b*B->m[4];  A->m[5] += b*B->m[5];
                                                  A->m[8] += b*B->m[8];
    m3d_sym_fill_lower(A);
}

inline void m3d_mult(m3d_t *M, const m3d_t *A, const m3d_t *B) { //乘法

    M->m[0] = A->m[0]*B->m[0] + A->m[1]*B->m[3] + A->m[2]*B->m[6];
    M->m[1] = A->m[0]*B->m[1] + A->m[1]*B->m[4] + A->m[2]*B->m[7];
    M->m[2] = A->m[0]*B->m[2] + A->m[1]*B->m[5] + A->m[2]*B->m[8];
    
    M->m[3] = A->m[3]*B->m[0] + A->m[4]*B->m[3] + A->m[5]*B->m[6];
    M->m[4] = A->m[3]*B->m[1] + A->m[4]*B->m[4] + A->m[5]*B->m[7];
    M->m[5] = A->m[3]*B->m[2] + A->m[4]*B->m[5] + A->m[5]*B->m[8];
    
    M->m[6] = A->m[6]*B->m[0] + A->m[7]*B->m[3] + A->m[8]*B->m[6];
    M->m[7] = A->m[6]*B->m[1] + A->m[7]*B->m[4] + A->m[8]*B->m[7];
    M->m[8] = A->m[6]*B->m[2] + A->m[7]*B->m[5] + A->m[8]*B->m[8];

}
inline void m3d_sym_mult(m3d_t *M, const m3d_t *A, const m3d_t *B) {

    M->m[0] = A->m[0]*B->m[0] + A->m[1]*B->m[3] + A->m[2]*B->m[6];
    M->m[1] = A->m[0]*B->m[1] + A->m[1]*B->m[4] + A->m[2]*B->m[7];
    M->m[2] = A->m[0]*B->m[2] + A->m[1]*B->m[5] + A->m[2]*B->m[8];
    
    M->m[3] = A->m[3]*B->m[0] + A->m[4]*B->m[3] + A->m[5]*B->m[6];
    M->m[4] = A->m[1]*B->m[3] + A->m[4]*B->m[4] + A->m[5]*B->m[7];
    M->m[5] = A->m[3]*B->m[2] + A->m[4]*B->m[5] + A->m[5]*B->m[8];
    
    M->m[6] = A->m[6]*B->m[0] + A->m[7]*B->m[3] + A->m[8]*B->m[6];
    M->m[7] = A->m[6]*B->m[1] + A->m[7]*B->m[4] + A->m[8]*B->m[7];
    M->m[8] = A->m[2]*B->m[6] + A->m[5]*B->m[7] + A->m[8]*B->m[8];

}
/** 
 * Product of two matrices that is known a propri to be symmetric
 */
inline void m3d_mult_sym_result(m3d_t *M, const m3d_t *A, const m3d_t *B) {

    M->m[0] = A->m[0]*B->m[0] + A->m[1]*B->m[3] + A->m[2]*B->m[6];
    M->m[1] = A->m[0]*B->m[1] + A->m[1]*B->m[4] + A->m[2]*B->m[7];
    M->m[2] = A->m[0]*B->m[2] + A->m[1]*B->m[5] + A->m[2]*B->m[8];
    
    M->m[4] = A->m[3]*B->m[1] + A->m[4]*B->m[4] + A->m[5]*B->m[7];
    M->m[5] = A->m[3]*B->m[2] + A->m[4]*B->m[5] + A->m[5]*B->m[8];
    
    M->m[8] = A->m[6]*B->m[2] + A->m[7]*B->m[5] + A->m[8]*B->m[8];

    m3d_sym_fill_lower(M);
}
/**
 * Produces the M=H*C*H product of two symmetric matrices H and C. 
 */
inline void m3d_mult_HCH_sym(m3d_t *M, const m3d_t *H, const m3d_t *C) {


    M->m[0] = ( + ( H->m[0]*C->m[0] + 2*H->m[1]*C->m[1] + 2*H->m[2]*C->m[2] )*H->m[0]
                + (                 +   H->m[1]*C->m[4] + 2*H->m[2]*C->m[5] )*H->m[1]
                + (                                     +   H->m[2]*C->m[8] )*H->m[2] );

    M->m[1] = ( + ( H->m[0]*C->m[0] + H->m[1]*C->m[1] + H->m[2]*C->m[2] )*H->m[1]
                + ( H->m[0]*C->m[1] + H->m[1]*C->m[4] + H->m[2]*C->m[5] )*H->m[4]
                + ( H->m[0]*C->m[2] + H->m[1]*C->m[5] + H->m[2]*C->m[8] )*H->m[5] );

    M->m[2] = ( + ( H->m[0]*C->m[0] + H->m[1]*C->m[1] + H->m[2]*C->m[2] )*H->m[2]
                + ( H->m[0]*C->m[1] + H->m[1]*C->m[4] + H->m[2]*C->m[5] )*H->m[5]
                + ( H->m[0]*C->m[2] + H->m[1]*C->m[5] + H->m[2]*C->m[8] )*H->m[8] );

    M->m[4] = ( + ( H->m[1]*C->m[0] + 2*H->m[4]*C->m[1] + 2*H->m[5]*C->m[2] )*H->m[1]
                + (                 +   H->m[4]*C->m[4] + 2*H->m[5]*C->m[5] )*H->m[4]
                + (                                     +   H->m[5]*C->m[8] )*H->m[5] );

    M->m[5] = ( + ( H->m[1]*C->m[0] + H->m[4]*C->m[1] + H->m[5]*C->m[2] )*H->m[2]
                + ( H->m[1]*C->m[1] + H->m[4]*C->m[4] + H->m[5]*C->m[5] )*H->m[5]
                + ( H->m[1]*C->m[2] + H->m[4]*C->m[5] + H->m[5]*C->m[8] )*H->m[8] );


    M->m[8] = ( + ( H->m[2]*C->m[0] + 2*H->m[5]*C->m[1] + 2*H->m[8]*C->m[2] )*H->m[2]
                + (                 +   H->m[5]*C->m[4] + 2*H->m[8]*C->m[5] )*H->m[5]
                + (                                     +   H->m[8]*C->m[8] )*H->m[8] );


    m3d_sym_fill_lower(M);
}

// Matrix ops

/**
 * Test if a matrix is symmetric. Return 0 if not.
 */
int m3d_is_sym(const m3d_t *A);
inline void m3d_scale(m3d_t *M, const double c) {
    M->m[0] *= c;  M->m[1] *= c;  M->m[2] *= c;
    M->m[3] *= c;  M->m[4] *= c;  M->m[5] *= c;
    M->m[6] *= c;  M->m[7] *= c;  M->m[8] *= c;
}
inline void m3d_sym_scale(m3d_t *M, const double c) {
    M->m[0] *= c;  M->m[1] *= c;  M->m[2] *= c;
                   M->m[4] *= c;  M->m[5] *= c;
                                  M->m[8] *= c;
    m3d_sym_fill_lower(M);
}



// Vectors -> Matrix

/**
 * Outer vector product: M = v^t * v .
 */
inline void v3d_outer(m3d_t *M, const v3d_t *v, const v3d_t *w)
{
    M->m[0] = v->v[0] * w->v[0];  M->m[1] = v->v[0] * w->v[1];  M->m[2] = v->v[0] * w->v[2];
    M->m[3] = v->v[1] * w->v[0];  M->m[4] = v->v[1] * w->v[1];  M->m[5] = v->v[1] * w->v[2];
    M->m[6] = v->v[2] * w->v[0];  M->m[7] = v->v[2] * w->v[1];  M->m[8] = v->v[2] * w->v[2];
}


// info
inline void v3d_info(int verb, const v3d_t *v, const char * str) {
    INFOV((verb, "%s -> (%g,%g,%g)",str, v->v[0], v->v[1], v->v[2]));
} 
inline void m3d_info(int verb, const m3d_t *M, char *desc)
{
    INFOV((verb,
           "%s\n"
           "  [%+1.4e, %+1.4e, %+1.4e]\n"
           "  [%+1.4e, %+1.4e, %+1.4e]\n"
           "  [%+1.4e, %+1.4e, %+1.4e]\n",
           desc,
           M->m[0], M->m[1], M->m[2],
           M->m[3], M->m[4], M->m[5],
           M->m[6], M->m[7], M->m[8]));
}

// Vector,Matrix -> Vector
inline void m3d_v3d_mult(v3d_t *u, const m3d_t *M, const v3d_t *v) 
{
    u->v[0] = M->m[0] * v->v[0] + M->m[1] * v->v[1] + M->m[2] * v->v[2];
    u->v[1] = M->m[3] * v->v[0] + M->m[4] * v->v[1] + M->m[5] * v->v[2];
    u->v[2] = M->m[6] * v->v[0] + M->m[7] * v->v[1] + M->m[8] * v->v[2];
}

// Vectors, Matrix -> Scalar

/**
 * Return the the scalar product return <v,Mw>.
 */
inline double v3d_m3d_v3d_inner(const v3d_t *v, const m3d_t *M, const v3d_t *w) 
{
    return (v->v[0] * ( M->m[0] * w->v[0] + M->m[1] * w->v[1] + M->m[2] * w->v[2] ) +
            v->v[1] * ( M->m[3] * w->v[0] + M->m[4] * w->v[1] + M->m[5] * w->v[2] ) +
            v->v[2] * ( M->m[6] * w->v[0] + M->m[7] * w->v[1] + M->m[8] * w->v[2] ) );
}
/*
 * Return the the scalar product return <v,Mv> for a symmetric matrix M.
返回标量积
 */
inline double v3d_m3d_sym_inner(const m3d_t *M, const v3d_t *v) 
{
#ifdef DEBUG_MATRIX_SYMMETRY_CHECK
    if (!m3d_is_sym(M))
        ERROR(("<v,Mv> for symmetric M called with non-symmetric M!"));
#endif
    return (v->v[0] * ( M->m[0] * v->v[0] + 2 * M->m[1] * v->v[1] + 2 * M->m[2] * v->v[2] ) +
            v->v[1] * (                         M->m[4] * v->v[1] + 2 * M->m[5] * v->v[2] ) +
            v->v[2] * (                                                 M->m[8] * v->v[2] ) );


}

inline void v4d_init(v4d_t *v, const double v1, const double v2, const double v3, const double v4)
{
    v->v[0]=v1;
    v->v[1]=v2;
    v->v[2]=v3;
    v->v[3]=v4;
}

inline void m4d_sym_fill_lower(m4d_t *M) {
    M->m[ 4]=M->m[ 1];
    M->m[ 8]=M->m[ 2];  M->m[ 9]=M->m[ 6];
    M->m[12]=M->m[ 3];  M->m[13]=M->m[ 7];  M->m[14]=M->m[11];
}
inline void m4d_sym_init(m4d_t *M, 
                                 double m11, double m12, double m13,  double m14,
                                 double m22, double m23, double m24,
                                 double m33, double m34,
                                 double m44) 
{
    M->m[ 0]=m11;  M->m[ 1]=m12;  M->m[ 2]=m13;  M->m[3]=m14;
                   M->m[ 5]=m22;  M->m[ 6]=m23;  M->m[7]=m24;
                                  M->m[10]=m33;  M->m[11]=m34;
                                                 M->m[15]=m44;
    m4d_sym_fill_lower(M);
}
inline void m4d_diag_init(m4d_t *M, 
                                  double m11, double m22, double m33, double m44) 
{
    memset(M->m, 0, 16*sizeof(double));
    M->m[ 0]=m11;
    M->m[ 5]=m22;
    M->m[10]=m33;
    M->m[15]=m44;
}


inline void v4d_accum_scaled(v4d_t *u, const double d, const v4d_t *v)
{
    u->v[0] += d*v->v[0];
    u->v[1] += d*v->v[1];
    u->v[2] += d*v->v[2];
    u->v[3] += d*v->v[3];
}


/**
 * Get a matrix column as a vector. First column number is 1.
 */
inline void m4d_v4d_get_col(v4d_t* v, const m4d_t *M, const int col)
{
    v->v[0] = M->m[     col-1];
    v->v[1] = M->m[4  + col-1];
    v->v[2] = M->m[8  + col-1];
    v->v[3] = M->m[12 + col-1];
}



inline void m4d_v4d_mult(v4d_t *u, const m4d_t *M, const v4d_t *v) 
{
    u->v[0] = M->m[ 0] * v->v[0] + M->m[ 1] * v->v[1] + M->m[ 2] * v->v[2] + M->m[ 3] * v->v[3];
    u->v[1] = M->m[ 4] * v->v[0] + M->m[ 5] * v->v[1] + M->m[ 6] * v->v[2] + M->m[ 7] * v->v[3];
    u->v[2] = M->m[ 8] * v->v[0] + M->m[ 9] * v->v[1] + M->m[10] * v->v[2] + M->m[11] * v->v[3];
    u->v[3] = M->m[12] * v->v[0] + M->m[13] * v->v[1] + M->m[14] * v->v[2] + M->m[15] * v->v[3];
}


/**
 * Eigenvalues and eigenvectors of a symmetric 4x4 matrix. The matrix
 * E will contain the eigenvectors as columns in the same order as the
 * eigenvalues.  If E is NULL, the eigenvectors will not be
 * computed. ev is asumed to be an array of 4 doubles.
 */
void m4d_sym_eigensystem(const m4d_t *M, double *ev, m4d_t *E);

/**
 * Eigenvalues of a symmetric 4x4 matrix, ev is assumed to be array of 4.
 */
inline void m4d_sym_eigenvals(const m4d_t *M, double *ev)
{
    m4d_sym_eigensystem(M,ev,NULL);
}

inline void m4d_info(int verb, const m4d_t *M, char *desc)
{
   INFOV((verb,
           "%s\n"
           "  [%+1.4e, %+1.4e, %+1.4e, %+1.4e]\n"
           "  [%+1.4e, %+1.4e, %+1.4e, %+1.4e]\n"
           "  [%+1.4e, %+1.4e, %+1.4e, %+1.4e]\n"
           "  [%+1.4e, %+1.4e, %+1.4e, %+1.4e]\n",
           desc,
           M->m[ 0], M->m[ 1], M->m[ 2], M->m[ 3], 
           M->m[ 4], M->m[ 5], M->m[ 6], M->m[ 7], 
           M->m[ 8], M->m[ 9], M->m[10], M->m[11], 
           M->m[12], M->m[13], M->m[14], M->m[15] ));
}


#endif