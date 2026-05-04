#include "../include/se_sparse_interpolation.h"


static int64_t idx2int(const rsf3d_t* rsf, int i1, int i2, int i3)
{
    int n1=rsf->g.z.n, n2=rsf->g.x.n;
    return (int64_t)(i3*n2 + i2) * n1 + i1;
}

static void int2idx(const rsf3d_t* rsf, int64_t i, int* k1, int* k2, int* k3)
{
    int i3, i12, i2, i1;
    int n1=rsf->g.z.n, n2=rsf->g.x.n;

    i3 = (int)(i / (n2 * n1));
    i12 = (int)(i - i3 * (n2 * n1));

    i2 = i12 / n1;
    i1 = i12 - i2 * n1;

    *k1 = i1;
    *k2 = i2;
    *k3 = i3;
}

typedef struct intarray_s {
    int64_t* d;
    int64_t size;
    int64_t capacity;
} intarray;

static void ensure_capacity(intarray* ia, int64_t required_size)
{
    if(required_size > ia->capacity) {
        int64_t request = required_size + 1024;

        ia->d = (int64_t *)realloc(ia->d,sizeof(int64_t)*request);
        if (ia->d == NULL) {
            ERROR(("unable to reallocate stack array"));
        }

        ia->capacity = request;
    }
}

static intarray* create_intarray()
{
    intarray* ia;
    ia = (intarray*)malloc(sizeof(intarray));
    if(ia == NULL) exit(1);

    ia->size = 0;
    ia->capacity = 0;
    ia->d = NULL;

    return ia;
}

static void destroy_intarray(intarray* ia)
{
    if(ia->d != NULL) free(ia->d);
    free(ia);
}


static void push(intarray* stack, int64_t i)
{
    ensure_capacity(stack, stack->size+1);
    stack->d[stack->size++] = i;
}

static int pop(intarray* stack)
{
    return stack->d[--stack->size];
}


static void mark_hole(const rsf3d_t* rsf,
                           bit_array this_hole, bit_array all_holes, 
                           intarray* boundary, intarray* stack,
                           int k1, int k2, int k3)
{
    int64_t idx = idx2int(rsf, k1, k2, k3);
    int n1=rsf->g.z.n, n2=rsf->g.x.n, n3=rsf->g.y.n;
    
    stack->size = 0;

    if( ba_get(all_holes, idx) ) 
        push(stack, idx);

    while(stack->size > 0) {
        int i1, i2, i3;
        int j1, j2, j3;
        int64_t crt = pop(stack);
        ba_set(this_hole, crt);
        int2idx(rsf, crt, &i1, &i2, &i3);

        for(j3 = i3-1; j3 <= i3+1; j3++) {
            if( j3 < 0 || j3 >= n3 ) continue;
            for(j2 = i2-1; j2 <= i2+1; j2++) {
                if( j2 < 0 || j2 >= n2 ) continue;
                for(j1 = i1-1; j1 <= i1+1; j1++) {
                    if( j1 < 0 || j1 >= n1 ) continue;

                    idx = idx2int(rsf, j1, j2, j3);
                    if( ba_get(this_hole, idx) ) continue;

                    if( ba_get(all_holes, idx) ) {
                        push(stack, idx);
                    } else {
                        push(boundary, idx);
                    }
                }
            }
        }
    }
}


#define IS_NO_VALUE(v) ( (!std::isfinite(v)) || (fabs( (v) - no_value ) < no_value_tolerance) )

void sparse_fill_holes_interpolator( rsf3d_t* rsf,
                                         float no_value, float no_value_tolerance,
                                         double dpow, int weight_only,
                                         double radius)
{
    int i3, i2, i1;
    int n1=rsf->g.z.n, n2=rsf->g.x.n, n3=rsf->g.y.n;
    double d1=rsf->g.z.d, d2=rsf->g.x.d, d3=rsf->g.y.d;    
    double o1=rsf->g.z.o, o2=rsf->g.x.o, o3=rsf->g.y.o;
    float *data=(float*)rsf->prvt;

    bit_array this_hole = ba_create(1024);
    bit_array all_holes = ba_create(1024);
    bit_array filled_holes = ba_create(1024);
    intarray* boundary = create_intarray();
    intarray* stack = create_intarray();

    (void)radius;

    for(i3 = 0; i3 < n3; i3++) {
        for(i2 = 0; i2 < n2; i2++) {
            for(i1 = 0; i1 < n1; i1++) {
                int64_t idx = idx2int(rsf, i1, i2, i3);
                if( IS_NO_VALUE( data[idx] ) ) {
                    ba_set(all_holes, idx);
                }
            }
        }
    }

    for(i3 = 0; i3 < n3; i3++) {
        for(i2 = 0; i2 < n2; i2++) {
            for(i1 = 0; i1 < n1; i1++) {
                int k3, k2, k1;
                int64_t idx = idx2int(rsf, i1, i2, i3);

                if( ba_get(filled_holes, idx) ) continue;
                if( ! ba_get(all_holes, idx) ) continue;

                INFOV((2, "Found a hole at %d,%d,%d - start marking boundary", i1, i2, i3));

                boundary->size = 0;
                ba_clear_all(this_hole);
                mark_hole(rsf, this_hole, all_holes, boundary, stack, i1, i2, i3);
                ba_or(filled_holes, this_hole);

                INFOV((2, "\tthe boundary has %d points - start filling the hole", boundary->size));

                if(boundary->size == 0) continue;

                for(k3 = 0; k3 < n3; k3++) {
                    for(k2 = 0; k2 < n2; k2++) {
                        for(k1 = 0; k1 < n1; k1++) {
                            int i, hole_idx = idx2int(rsf, k1, k2, k3);
                            float x01, x02, x03;
                            double sum, tot;
                            if( ! ba_get(this_hole, hole_idx) ) continue;

                            sum = 0.0;
                            tot = 0.0;
                            x01 = (float)(o1 + d1 * (k1 + 0.5));
                            x02 = (float)(o2 + d2 * (k2 + 0.5));
                            x03 = (float)(o3 + d3 * (k3 + 0.5));
                            for(i = 0; i < boundary->size; i++) {
                                int j1, j2, j3;
                                float v;
                                double d, x1, x2, x3;
                                int2idx(rsf, boundary->d[i], &j1, &j2, &j3);
                                v = data[boundary->d[i]];
                                x1 = (float)(o1 + d1 * (j1 + 0.5));
                                x2 = (float)(o2 + d2 * (j2 + 0.5));
                                x3 = (float)(o3 + d3 * (j3 + 0.5));
                                d = (x01-x1)*(x01-x1) + 
                                    (x02-x2)*(x02-x2) + 
                                    (x03-x3)*(x03-x3) ;
                                d = pow(d, -0.5*dpow);
                                tot += d;
                                sum += d*v;
                            }
                            if(weight_only) {
                                data[hole_idx] = (float)(tot);
                            } else {
                                data[hole_idx] = (float)(sum/tot);
                            }
                        }
                    }
                }
            }
        }
    }

    if(weight_only) {
        for(i3 = 0; i3 < n3; i3++) {
            for(i2 = 0; i2 < n2; i2++) {
                for(i1 = 0; i1 < n1; i1++) {
                    int64_t idx = idx2int(rsf, i1, i2, i3);
                    if( ! ba_get(all_holes, idx) ) {
                        data[idx] = 0.0;
                    }
                }
            }
        }
    }

    destroy_intarray(boundary);
    destroy_intarray(stack);
    ba_destroy(this_hole);
    ba_destroy(all_holes);
    ba_destroy(filled_holes);
}

/*
static void sparse_onehole_interpolator_with_radius( rsf3d_t* rsf,
                                                              float no_value, float no_value_tolerance,
                                                              double dpow, int weight_only,
                                                              double radius)
{
    int i3, i2, i1;
    int m1, m2, m3;
    int n1=rsf->g.z.n, n2=rsf->g.x.n, n3=rsf->g.y.n;
    double d1=rsf->g.z.d, d2=rsf->g.x.d, d3=rsf->g.y.d;
    double o1=rsf->g.z.o, o2=rsf->g.x.o, o3=rsf->g.y.o;
    float *data = (float*)rsf->prvt;

    bit_array all_alive;
    int n_holes_left;
    int nalives;

    m1 = (int)(0.5 * radius / d1);
    m2 = (int)(0.5 * radius / d2);
    m3 = (int)(0.5 * radius / d3);

    all_alive = ba_create(n3*n2*n1);
    nalives = 0;
    for(i3 = 0; i3 < n3; i3++) {
        for(i2 = 0; i2 < n2; i2++) {
            for(i1 = 0; i1 < n1; i1++) {
                int64_t idx = idx2int(rsf, i1, i2, i3);
                if( ! IS_NO_VALUE( data[idx] ) ) {
                    ba_set(all_alive, idx);
                    nalives += 1;
                }
            }
        }
    }

    if(nalives == 0) {
        INFOV((1, "There are no live points - nothing to do here"));
        return;
    }

    INFOV((2, "There are %d total live points; filling within radius %g", nalives, radius));


    n_holes_left = 0;
    for(i3 = 0; i3 < n3; i3++) {
        for(i2 = 0; i2 < n2; i2++) {
            for(i1 = 0; i1 < n1; i1++) {
                float x01, x02, x03;
                int k;
                double sum, tot;
                int j1, j2, j3;
                int npt;
                int64_t idx = idx2int(rsf, i1, i2, i3);

                if( ba_get(all_alive, idx) ) continue;

                sum = 0.0;
                tot = 0.0;
                x01 = (float)(o1 + d1 * (i1 + 0.5));
                x02 = (float)(o2 + d2 * (i2 + 0.5));
                x03 = (float)(o3 + d3 * (i3 + 0.5));

                npt = 0;
                k = 0;
                while(npt == 0) {
                    for(j3 = i3-(m3+k); j3 < i3+(m3+k); j3++) {
                        if( j3 < 0 || j3 > n3 ) continue;
                        for(j2 = i2-(m2+k); j2 < i2+(m2+k); j2++) {
                            if( j2 < 0 || j2 > n2 ) continue;
                            for(j1 = i1-(m1+k); j1 < i1+(m1+k); j1++) {
                                double w, d, x1, x2, x3;
                                int jdx;

                                if( j1 < 0 || j1 > n1 ) continue;
                                jdx = idx2int(rsf, j1, j2, j3);
                                if( ! ba_get(all_alive, jdx) ) continue;
                                npt++;

                                x1 = o1 + d1 * (j1 + 0.5);
                                x2 = o2 + d2 * (j2 + 0.5);
                                x3 = o3 + d3 * (j3 + 0.5);
                                d = (x01-x1)*(x01-x1) + 
                                    (x02-x2)*(x02-x2) + 
                                    (x03-x3)*(x03-x3) ;

                                w = pow(d, -0.5*dpow);

                                tot += w;
                                sum += w*data[jdx];
                            }
                        }
                    }
                    if(npt == 0) {
                        n_holes_left += 1;
                    } else {
                        data[idx] = (float)(sum/tot);
                    }
                }
            }
        }
    }
    
    ba_destroy(all_alive);
    
    if(n_holes_left > 0) {
        INFOV((2, "There are %d points left to fill", n_holes_left));
        sparse_onehole_interpolator( rsf, no_value, no_value_tolerance,
                                         dpow, weight_only, -1.0 );
    }
}
*/

void sparse_onehole_interpolator( rsf3d_t* rsf,
                                      float no_value, float no_value_tolerance,
                                      double dpow, int weight_only,
                                      double radius)
{
    int i3, i2, i1;
    float min = 1.0E20, max = -1.0E20;
    intarray* live_samples;
    float *data = (float*)rsf->prvt;
    int n1=rsf->g.z.n, n2=rsf->g.x.n, n3=rsf->g.y.n;
    double d1=rsf->g.z.d, d2=rsf->g.x.d, d3=rsf->g.y.d;
    double o1=rsf->g.z.o, o2=rsf->g.x.o, o3=rsf->g.y.o;
    (void)radius;

    /*
    if(radius > 0) {
        sparse_onehole_interpolator_with_radius( rsf, no_value, no_value_tolerance,
                                                     dpow, weight_only, radius );
        return;
    }
    */

    live_samples = create_intarray();
    for(i3 = 0; i3 < n3; i3++) {
        for(i2 = 0; i2 < n2; i2++) {
            for(i1 = 0; i1 < n1; i1++) {
                int64_t idx = idx2int(rsf, i1, i2, i3);
                if( ! IS_NO_VALUE( data[idx] ) ) {
                    push(live_samples, idx);
                }
            }
        }
    }
    
    if(live_samples->size == 0) {
        INFOV((1, "There are no live points - nothing to do here"));
        return;
    }

    INFOV((2, "There are %d live points", live_samples->size));

    for(i3 = 0; i3 < n3; i3++) {
        for(i2 = 0; i2 < n2; i2++) {
            for(i1 = 0; i1 < n1; i1++) {
                float x01, x02, x03;
                double sum, tot;
                int i;
                int64_t idx = idx2int(rsf, i1, i2, i3);

                /* since we are filling the data sequentially, there
                   is no danger in overwriting a 'live' sample */
                if( !weight_only && !IS_NO_VALUE( data[idx] ) )
                    continue;

                sum = 0.0;
                tot = 0.0;
                x01 = (float)(o1 + d1 * (i1 + 0.5));
                x02 = (float)(o2 + d2 * (i2 + 0.5));
                x03 = (float)(o3 + d3 * (i3 + 0.5));
                for(i = 0; i < live_samples->size; i++) {
                    int j1, j2, j3;
                    double w, d, x1, x2, x3;
                    int2idx(rsf, live_samples->d[i], &j1, &j2, &j3);

                    x1 = (float)(o1 + d1 * (j1 + 0.5));
                    x2 = (float)(o2 + d2 * (j2 + 0.5));
                    x3 = (float)(o3 + d3 * (j3 + 0.5));
                    d = (x01-x1)*(x01-x1) + 
                        (x02-x2)*(x02-x2) + 
                        (x03-x3)*(x03-x3) ;

                    if( idx != live_samples->d[i] ) {
                        w = pow(d, -0.5*dpow);
                    } else {
                        w = 0;
                    }

                    if(weight_only) {
                        tot = 1;
                        sum += w;
                    } else {
                        float v = data[live_samples->d[i]];
                        tot += w;
                        sum += w*v;
                    }
                }

                data[idx] = (float)(sum/tot);
                if(weight_only) {
                    if(min > data[idx]) min = data[idx];
                    if(max < data[idx]) max = data[idx];
                }
            }
        }
    }

    if(weight_only) {
        INFOV((2, "Min/max for this section: %f/%f", min, max));
        for(i3 = 0; i3 < n3; i3++) {
            for(i2 = 0; i2 < n2; i2++) {
                for(i1 = 0; i1 < n1; i1++) {
                    int64_t idx = idx2int(rsf, i1, i2, i3);
                    data[idx] = (data[idx] - min)/(max - min);
                }
            }
        }
        for(i1 = 0; i1 < live_samples->size; i1++) {
            data[live_samples->d[i1]] = 1.0;
        }
    }

    destroy_intarray(live_samples);
}


static void internal_linear_interpolator( rsf3d_t* rsf,
                                               float no_value, float no_value_tolerance,
                                               double dpow, int weight_only,
                                               double radius, const int nn)
{
    int i3, i2, i1, i;
    int n1=rsf->g.z.n, n2=rsf->g.x.n, n3=rsf->g.y.n;
    float *data=(float*)rsf->prvt;
    (void)dpow;
    (void)weight_only;
    (void)radius;

    for(i3 = 0; i3 < n3; i3++) {
        for(i2 = 0; i2 < n2; i2++) {
            for(i1 = 0; i1 < n1; i1++) {
                int64_t idx = idx2int(rsf, i1, i2, i3);
                if( IS_NO_VALUE( data[idx] ) ) {
                    data[idx] = no_value;
                }
            }
        }
    }

    for(i3 = 0; i3 < n3; i3++) {
        for(i2 = 0; i2 < n2; i2++) {
            int first_value = 0;
            int iold = 0;
            for(i1 = 0; i1 < n1; i1++) {
                int64_t idx = idx2int(rsf, i1, i2, i3);
                if( ! IS_NO_VALUE( data[idx] ) ) {
                    if( first_value == 0 ) {
                        for(i=0; i < i1; i++) {
                            int64_t idst = (i3 * n2 + i2) * n1 + i;
                            int64_t isrc = idx;
                            data[idst] = data[isrc];
                        }
                        first_value = 1;
                        iold = i1;
                    } else {
                        float dv = (float)(i1 - iold);
                        for(i = iold+1; i < i1; i++) {
                            float dv1 = (i - iold) / dv;
                            float dv2 = 1.0f - dv1;
                            int64_t idst  = (i3 * n2 + i2) * n1 + i;
                            int64_t isrc1 = (i3 * n2 + i2) * n1 + iold;
                            int64_t isrc2 = idx;
                            if(nn) {
                                data[idst] = dv1<0.5?data[isrc1]:data[isrc2];
                            } else {
                                data[idst] = dv2 * data[isrc1] + dv1 * data[isrc2];
                            }
                        }
                        iold = i1;
                    }
                }
            }
            if(iold != n1) {
                for(i = iold+1; i < n1; i++) {
                    int64_t idst = (i3 * n2 + i2) * n1 + i;
                    int64_t isrc = (i3 * n2 + i2) * n1 + iold;
                    data[idst] = data[isrc];
                }
            }
        }
    }

    for(i3 = 0; i3 < n3; i3++) {
	int first_value = 0;
        int iold = 0;

        for(i2 = 0; i2 < n2; i2++) {
            int64_t idx = idx2int(rsf, 0, i2, i3);

            if( ! IS_NO_VALUE( data[idx] ) ) {
                if( first_value == 0 ) {
                    for(i=0; i < i2; i++) {
                        for(i1 = 0; i1 < n1; i1++) {
                            int64_t idst = (i3 * n2 + i ) * n1 + i1;
                            int64_t isrc = idx+i1;
                            data[idst] = data[isrc];
                        }
                    }
                    first_value = 1;
                    iold = i2;
                } else {
                    float dv = (float)(i2 - iold);
                    for(i = iold+1; i < i2; i++) {
                        float dv1 = (i - iold) / dv;
                        float dv2 = 1.0f - dv1;
                        for(i1 = 0; i1 < n1; i1++) {
                            int64_t idst  = (i3 * n2 + i   ) * n1 + i1;
                            int64_t isrc1 = (i3 * n2 + iold) * n1 + i1;
                            int64_t isrc2 = idx+i1;
                            if(nn) {
                                data[idst] = dv1<0.5?data[isrc1]:data[isrc2];
                            } else {
                                data[idst] = dv2 * data[isrc1] + dv1 * data[isrc2];
                            }
                        }
                    }
	       
                    iold = i2;
                }
            }
	}

        if(iold != n2) {
            for(i = iold+1; i < n2; i++) {
                for(i1 = 0; i1 < n1; i1++) {
                    int64_t idst = (i3 * n2 + i   ) * n1 + i1;
                    int64_t isrc = (i3 * n2 + iold) * n1 + i1;
                    data[idst] = data[isrc];
                }
            }
	}
    }

    for(i2 = 0; i2 < n2; i2++) {
	int first_value = 0;
        int iold = 0;

        for(i3 = 0; i3 < n3; i3++) {
            int64_t idx = idx2int(rsf, 0, i2, i3);

            if( ! IS_NO_VALUE( data[idx] ) ) {
                if( first_value == 0 ) {
                    for(i=0; i < i3; i++) {
                        for(i1 = 0; i1 < n1; i1++) {
                            int64_t idst = (i  * n2 + i2) * n1 + i1;
                            int64_t isrc = idx+i1;
                            data[idst] = data[isrc];
                        }
                    }
                    first_value = 1;
                    iold = i3;
                } else {
                    float dv = (float)(i3 - iold);
                    for(i = iold+1; i < i3; i++) {
                        float dv1 = (i - iold) / dv;
                        float dv2 = 1.0f - dv1;
                        for(i1 = 0; i1 < n1; i1++) {
                            int64_t idst  = (i    * n2 + i2) * n1 + i1;
                            int64_t isrc1 = (iold * n2 + i2) * n1 + i1;
                            int64_t isrc2 = idx+i1;
                            if(nn) {
                                data[idst] = dv1<0.5?data[isrc1]:data[isrc2];
                            } else {
                                data[idst] = dv2 * data[isrc1] + dv1 * data[isrc2];
                            }
                        }
                    }
	       
                    iold = i3;
                }
            }
	}

        if(iold != n3) {
            for(i = iold+1; i < n3; i++) {
                for(i1 = 0; i1 < n1; i1++) {
                    int64_t idst = (i    * n2 + i2) * n1 + i1;
                    int64_t isrc = (iold * n2 + i2) * n1 + i1;
                    data[idst] = data[isrc];
                }
            }
	}
    }
}

void sparse_linear_interpolator( rsf3d_t* rsf,
                                     float no_value, float no_value_tolerance,
                                     double dpow, int weight_only,
                                     double radius)
{
    internal_linear_interpolator(rsf, no_value, no_value_tolerance, 
                                 dpow, weight_only, radius, 0);
}

void sparse_nn_interpolator( rsf3d_t* rsf,
                                 float no_value, float no_value_tolerance,
                                 double dpow, int weight_only,
                                 double radius)
{
    internal_linear_interpolator(rsf, no_value, no_value_tolerance, 
                                 dpow, weight_only, radius, 1);
}

void sparse_trace_linear_interpolator( float *tr, const axa_t *ax, float no_value, float no_value_tolerance) {
    int i, i1, n1=ax->n;

    for(i1 = 0; i1 < n1; i1++) {
        if( IS_NO_VALUE( tr[i1] ) ) {
            tr[i1] = no_value;
        }
    }

    int first_value = 0;
    int iold = 0;
    for(i1 = 0; i1 < n1; i1++) {
        if( ! IS_NO_VALUE( tr[i1] ) ) {
            if( first_value == 0 ) {
                for(i=0; i < i1; i++) {
                    tr[i] = tr[i1];
                }
                first_value = 1;
                iold = i1;
            } else {
                float dv = (float)(i1 - iold);
                for(i = iold+1; i < i1; i++) {
                    float dv1 = (i - iold) / dv;
                    float dv2 = 1.0f - dv1;
                    tr[i] = dv2 * tr[iold] + dv1 * tr[i1];
                }
                iold = i1;
            }
        }
    }
    if(iold != n1) {
        for(i = iold+1; i < n1; i++) {
            tr[i] = tr[iold];
        }
    }
}
