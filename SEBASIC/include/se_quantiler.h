#ifndef SE_QUANTILER_H
#define SE_QUANTILER_H
#include "se_type.h"
#include "se_alloc.h"
#include "se_assert.h"
#include "se_log.h"
#include "se_basic_math.h"


/**
 * Structure storing the internal state of a quantiler.
 */
typedef struct se_quantiler_s {
    real q;              /**< the desired quantile */
    real fnull;          /**< null sample value, if ignoring nulls */
    real m0,m1,m2,m3,m4; /**< marker positions */
    real q0,q1,q2,q3,q4; /**< marker heights */
    real    f1,f2,f3   ; /**< desired marker positions */
    real    d1,d2,d3   ; /**< desired marker position increments */
    int ignore_null;     /**< true if ignoring fnull samples */
    int initialized;     /**< true if estimator has been initialized */
    int isminmax;        /**< if 1, look for min, if 2 look for max */
} se_quantiler_t;


/**
 * Constructs a quantiler for the specified quantile fraction and an
 * optional null value.  
 *
 * \param[in,out] zqt the quantiler to be initialized.
 *
 * \param[in] q the target quatile fraction; it must be between 0 and
 *              1 inclusive.
 *
 * \param[in] ignore_null flag to specify whether the updates should
 *                        ignore the null value or not.
 *
 * \param[in] fnull the null value, use if ignore_null is true.
 */
inline void se_quantiler_init(se_quantiler_t* zqt, real q, int ignore_null, real fnull)
{
    ASSERT(q >= 0.0);
    ASSERT(q <= 1.0);

    memset(zqt, 0, sizeof(*zqt));

    zqt->q = q;
    zqt->m0 = -1.0;
    if(ignore_null) {
        zqt->ignore_null = 1;
        zqt->q2 = zqt->fnull = fnull;
    } else {
        zqt->ignore_null = 0;
        zqt->q2 = 0.0;
    }

    zqt->initialized = 0;

    if(FEQUAL(zqt->q, 0.0)) {
        zqt->isminmax = 1;
    } else if(FEQUAL(zqt->q, 1.0)) {
        zqt->isminmax = 2;
    } else {
        zqt->isminmax = 0;
    }
}

/**
 * Returns the current quantile estimate.
 * \param[in] zqt the quantiler.
 * \returns the current quantile estimate.
 */
inline real se_quantiler_estimate(const se_quantiler_t* zqt)
{
    return zqt->q2;
}

/**
 * Returns the quantile fraction.
 * \param[in] zqt the quantiler.
 * \returns the quantile fraction.
 */
inline real se_quantiler_fraction(const se_quantiler_t* zqt)
{
    return zqt->q;
}


/************************** PRIVATE ********************************/

/**
 * Helper function used by \link _se_quantiler_update_one \endlink.
 */
inline real _se_quantiler_qp(real mp,
                                  real m0, real m1, real m2,
                                  real q0, real q1, real q2)
{
    real qt = q1 + ( (mp-m0)*(q2-q1)/(m2-m1) + (m2-mp)*(q1-q0)/(m1-m0) ) / (m2-m0);
    if(qt <= q2) {
        return qt;
    } else {
        return q1 + (q2-q1) / (m2-m1);
    }
}

/**
 * Helper function used by \link _se_quantiler_update_one \endlink.
 */
inline real _se_quantiler_qm(real mm,
                                  real m0, real m1, real m2,
                                  real q0, real q1, real q2)
{
    real qt = q1 - ( (mm-m0)*(q2-q1)/(m2-m1) + (m2-mm)*(q1-q0)/(m1-m0) ) / (m2-m0);
    if(q0 <= qt) {
        return qt;
    } else {
        return q1 + (q0-q1)/(m0-m1);
    }
}


/**
 * The function that actually updates the quantiler once it is initialized. 
 * \attention This function is not to be called directly - use \link se_quantiler_update \endlink for this.
 * \see _se_quantiler_init_one, se_quantiler_update
 * \param[in,out] zqt the quantiler to be updated.
 * \param[in] f the value to use to update the quantiler.
 */
inline void _se_quantiler_update_one(se_quantiler_t* zqt, real f)
{
    ASSERT(zqt->initialized);

    /* If ignoring null samples, check for null. */
    if (zqt->ignore_null && FEQUAL(f, zqt->fnull) )
        return;
      
    if(zqt->isminmax == 0) {
        real mm, mp;

        /* Increment marker locations and update min and max. */
        if(f < zqt->q0) {
            zqt->m1 += 1.0;
            zqt->m2 += 1.0;
            zqt->m3 += 1.0;
            zqt->m4 += 1.0;
            zqt->q0 = f;
        } else if (f < zqt->q1) {
            zqt->m1 += 1.0;
            zqt->m2 += 1.0;
            zqt->m3 += 1.0;
            zqt->m4 += 1.0;
        } else if (f < zqt->q2) {
            zqt->m2 += 1.0;
            zqt->m3 += 1.0;
            zqt->m4 += 1.0;
        } else if (f < zqt->q3) {
            zqt->m3 += 1.0;
            zqt->m4 += 1.0;
        } else if (f < zqt->q4) {
            zqt->m4 += 1.0;
        } else {
            zqt->m4 += 1.0;
            zqt->q4 = f;
        }
      
        /* Increment desired marker positions. */
        zqt->f1 += zqt->d1;
        zqt->f2 += zqt->d2;
        zqt->f3 += zqt->d3;

        /* If necessary, adjust height and location of markers 1, 2, and 3. */
        mm = zqt->m1 - 1.0;
        mp = zqt->m1 + 1.0;

        if (zqt->f1 >= mp && zqt->m2 > mp) {
            zqt->q1 = _se_quantiler_qp(mp, zqt->m0, zqt->m1, zqt->m2, zqt->q0, zqt->q1, zqt->q2);
            zqt->m1 = mp;
        } else if (zqt->f1 <= mm && zqt->m0 < mm) {
            zqt->q1 = _se_quantiler_qm(mm, zqt->m0, zqt->m1, zqt->m2, zqt->q0, zqt->q1, zqt->q2);
            zqt->m1 = mm;
        }

        mm = zqt->m2 - 1.0;
        mp = zqt->m2 + 1.0;
        
        if (zqt->f2 >= mp && zqt->m3 > mp) {
            zqt->q2 = _se_quantiler_qp(mp, zqt->m1, zqt->m2, zqt->m3, zqt->q1, zqt->q2, zqt->q3);
            zqt->m2 = mp;
        } else if (zqt->f2 <= mm && zqt->m1 < mm) {
            zqt->q2 = _se_quantiler_qm(mm, zqt->m1, zqt->m2, zqt->m3, zqt->q1, zqt->q2, zqt->q3);
            zqt->m2 = mm;
        }

        mm = zqt->m3 - 1.0;
        mp = zqt->m3 + 1.0;
        
        if (zqt->f3 >= mp && zqt->m4 > mp) {
            zqt->q3 = _se_quantiler_qp(mp, zqt->m2, zqt->m3, zqt->m4, zqt->q2, zqt->q3, zqt->q4);
            zqt->m3 = mp;
        } else if (zqt->f3<=mm && zqt->m2<mm) {
            zqt->q3 = _se_quantiler_qm(mm, zqt->m2, zqt->m3, zqt->m4, zqt->q2, zqt->q3, zqt->q4);
            zqt->m3 = mm;
        }
    } else if(zqt->isminmax == 1) {
        if(f < zqt->q2)
            zqt->q2 = f;
    } else {
        if (f > zqt->q2)
            zqt->q2 = f;
    }
}

/**
 * Helper function use to update the quantiler for the first samples
 * (usualy first 5).
 * \attention This function is not to be called directly - use \link se_quantiler_update \endlink for this.
 * \param[in,out] zqt the quantiler to be updated.
 * \param[in] f the value to use to update the quantiler.
 */
void _se_quantiler_init_one(se_quantiler_t* zqt, real f);

/************************** END OF PRIVATE ********************************/

/**
 * Updates the quantile estimate with the specified sample.
 * \param[in,out] zqt the quantiler to be updated.
 * \param[in] f the value to use to update the quantiler.
 * \return the updated quantile estimate.
 */
inline real se_quantiler_update(se_quantiler_t* zqt, real f)
{
    if (!zqt->initialized) {
        _se_quantiler_init_one(zqt, f);
    } else {
        _se_quantiler_update_one(zqt, f);
    }
    return se_quantiler_estimate(zqt);
}

/**
 * Updates the quantiler from an array of values.
 * \param[in,out] zqt the quantiler to be updated.
 * \param[in] f the array of values to use to update the quantiler.
 * \param[in] n the number of samples from the array.
 * \return the updated quantile estimate.
 */
real se_quantiler_update_array(se_quantiler_t* zqt, const real* f, size_t n);

/**
 * Updates the quantiler from an array of single precision floating point values.
 * \param[in,out] zqt the quantiler to be updated.
 * \param[in] f the array of values to use to update the quantiler.
 * \param[in] n the number of samples from the array.
 * \return the updated quantile estimate.
 */
real se_quantiler_update_farray(se_quantiler_t* zqt, const float* f, size_t n);


#endif