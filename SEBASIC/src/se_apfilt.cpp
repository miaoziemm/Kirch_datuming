#include "../include/se_apfilt.h"


typedef struct _se_apfilt_s {
    int n;
    double *b;
} _se_apfilt_t;

se_apfilt_t se_apfilt_init(int nw /* filter order */)
{
    int j, k;
    double bk;
    _se_apfilt_t* ap = (_se_apfilt_t*)malloc(sizeof(*ap));

    ap->n = nw*2;
    ap->b = alloc1double(ap->n+1);

    for (k=0; k <= ap->n; k++) {
        bk = 1.0;
        for (j=0; j < ap->n; j++) {
            if (j < ap->n-k) {
                bk *= (k+j+1.0)/(2*(2*j+1)*(j+1));
            } else {
                bk *= 1.0/(2*(2*j+1));
            }
        }
        ap->b[k] = bk;
    }
    return(se_apfilt_t)ap;
}

void se_apfilt_destroy(se_apfilt_t* _ap)
{
    _se_apfilt_t* ap = (_se_apfilt_t*)_ap;
    free(ap->b);
    // free(ap);
}

/** find filter coefficients */
void se_passfilter (se_apfilt_t _ap,
                      float p  /* slope */, 
                      float* a /* output filter [n+1] */)
{
    _se_apfilt_t* ap = (_se_apfilt_t*)_ap;
    int j, k;
    double ak;
    
    for (k=0; k <= ap->n; k++) {
        ak = ap->b[k];
        for (j=0; j < ap->n; j++) {
            if (j < ap->n-k) {
                ak *= (ap->n-j-p);
            } else {
                ak *= (p+j+1);
            }
        }
        a[k] = (float)ak;
    }
}

/** find coefficients for filter derivative */
void se_aderfilter (se_apfilt_t _ap,
                      float p  /* slope */, 
                      float* a /* output filter [n+1] */)
{
    _se_apfilt_t* ap = (_se_apfilt_t*)_ap;
    int i, j, k;
    float ak, ai;
    
    for (k=0; k <= ap->n; k++) {
        ak = 0.;
        for (i=0; i < ap->n; i++) {
            ai = -1.0;
            for (j=0; j < ap->n; j++) {
                if (j != i) {           
                    if (j < ap->n-k) {
                        ai *= (ap->n-j-p);
                    } else {
                        ai *= (p+j+1);
                    }
                } else if (j < ap->n-k) {
                    ai *= (-1);
                }
            }
            ak += ai;
        }
        a[k] = (float)(ak*ap->b[k]);
    }
}
