#include "../include/se_allpass.h"
#include "../include/se_apfilt.h"

typedef struct _se_allpass_s {
    int nx, ny, nz, nw, nj;
    float *flt, *pp;

    se_allpass_t apf;

    float *flt_tbl;
    int tbl_n;
    double tbl_d;
    double tbl_o;
    int *flt_mini;
    int *flt_maxi;

} _se_allpass_t;


se_allpass_t se_allpass_init(int nw,
                                 int nj,
                                 int nx, int ny, int nz,
                                 float *pp)
{
    _se_allpass_t* ap = (_se_allpass_t*)malloc(sizeof(*ap));
    memset(ap,0,sizeof(*ap));

    ap->nw = nw;
    ap->nj = nj;
    ap->nx = nx;
    ap->ny = ny;
    ap->nz = nz;
    ap->pp = pp;

    ap->flt = alloc1float(2*nw+1);
    ap->apf = se_apfilt_init(nw);

    ap->flt_tbl = NULL;
    ap->tbl_d = 1;
    ap->tbl_o = 0;
    ap->tbl_n = 0;
    ap->flt_mini = NULL;
    ap->flt_maxi = NULL;

    return (se_allpass_t)ap;
}


void se_allpass_destroy(se_allpass_t _ap)
{
    _se_allpass_t* ap = (_se_allpass_t*)_ap;
    se_apfilt_destroy(&ap->apf);
    free(ap->flt);
    if (ap->flt_tbl==NULL) {
        free(ap->flt_tbl);
        ap->flt_tbl = NULL;
        ap->tbl_d = 1;
        ap->tbl_o = 0;
        ap->tbl_n = 0;
    }
    if (ap->flt_mini==NULL) {
        free(ap->flt_mini);
        ap->flt_mini=NULL;
    }
    if (ap->flt_maxi==NULL) {
        free(ap->flt_maxi);
        ap->flt_maxi=NULL;
    }
    free(ap);

}
//初始化查找表
void se_allpass_init_lookup_table(se_allpass_t _ap,
                                    double tbl_min, double tbl_max, int tbl_n) 

{
    _se_allpass_t* ap = (_se_allpass_t*)_ap;
    int i;
    int ns=2*ap->nw+1;

    if ( ap->tbl_n != tbl_n ) {
        ap->flt_tbl = (float *)realloc(ap->flt_tbl, sizeof(float) * tbl_n * ns);
        ap->flt_mini = (int *)realloc(ap->flt_mini, sizeof(int) * tbl_n);
        ap->flt_maxi = (int *)realloc(ap->flt_maxi, sizeof(int) * tbl_n);
    }

    ap->tbl_n = tbl_n;

    if (ap->tbl_n > 1) {
        ap->tbl_d = (tbl_max - tbl_min) /  (ap->tbl_n-1) ;
        ap->tbl_o = tbl_min;
    }
    else {
        ap->tbl_d = 0;
        ap->tbl_o = (tbl_max + tbl_min) /  2.0;
    }

    for (i = 0 ; i < ap->tbl_n; i++) {
        double dip;

        dip = ap->tbl_o + i*ap->tbl_d;

        se_passfilter(ap->apf, (float)dip, ap->flt_tbl + ns*i);

        ap->flt_mini[i] = 0;
        while  ( (ap->flt_mini[i]<2*ap->nw) && 
                 (fabs( ap->flt_tbl[ns*i+ap->flt_mini[i]]) < 1e-5) ) 
            ap->flt_mini[i]++;
            
        ap->flt_maxi[i]=2*ap->nw;
        while  ( (ap->flt_maxi[i]>0) && 
                 (fabs(ap->flt_tbl[ns*i + ap->flt_maxi[i]]) < 1e-5) ) 
            ap->flt_maxi[i]--;

    }

}

void se_allpass_lookup_table_init_antialias(se_allpass_t _ap)
{

     _se_allpass_t* ap = (_se_allpass_t*)_ap;
    int i;
    int ns=2*ap->nw+1;
    float *pwd;
    double *kernel;

    pwd =alloc1float(ns);
    kernel =alloc1double(ns);

    for (i = 0 ; i < ap->tbl_n; i++) {

        /* double dip = ap->tbl_o + i*ap->tbl_d; */
        int j;

        memcpy(pwd,ap->flt_tbl+ns*i,ns*sizeof(float));

        for (j=0; j<ns; j++) {
            ap->flt_tbl[i*ns+j] = 0;
            int k;
            double rho;
            double re_norm;
            rho = 200;//1/(fabs(dip)+1e-8);
            re_norm = 0;
            for (k=0; k<ns; k++) {
                kernel[k] = exp(-rho*(k-j)*(k-j));    //反假频
                re_norm += kernel[k];
            }
            re_norm = 1/re_norm;
            for (k=0; k<ns; k++) {
                ap->flt_tbl[i*ns+j] += (float)(pwd[k]*re_norm*kernel[k]);
            }
        }

        ap->flt_mini[i] = 0;
        while  ( (ap->flt_mini[i]<2*ap->nw) && 
                 (fabs( ap->flt_tbl[ns*i+ap->flt_mini[i]]) < 1e-5) ) 
            ap->flt_mini[i]++;
            
        ap->flt_maxi[i]=2*ap->nw;
        while  ( (ap->flt_maxi[i]>0) && 
                 (fabs(ap->flt_tbl[ns*i + ap->flt_maxi[i]]) < 1e-5) ) 
            ap->flt_maxi[i]--;


    }
    
    free(kernel);
    free(pwd);

}


void se_allpass1(const se_allpass_t _ap, 
                   const int left,
                   const int der, 
                   const float* xx, 
                   float* yy)
{
    _se_allpass_t* ap = (_se_allpass_t*)_ap;
    int nx, ny, nz, i1, i2, ip, dir;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    
    if (left) {
        i1=1; i2=ny;   ip=-nx; dir=-1;
    } else {
        i1=0; i2=ny-1; ip= nx; dir= 1;
    }

    memset(yy, 0, nx*ny*nz*sizeof(float));
 
    int iz;
    for (iz=0; iz < nz; iz++) {
        int iy;           
        for (iy=i1; iy < i2; iy++) {
            int ix;
            for (ix = 0; ix < nx; ix++) {
                int i = ix + nx * (iy + ny * iz);

                float p    = dir*ap->pp[i];
                int pint   = (int) p;
                float pdec = p - pint;

                if (der) {
                    se_aderfilter(ap->apf, pdec, ap->flt);
                } else {
                    se_passfilter(ap->apf, pdec, ap->flt);
                }
          
                int iw;
                for (iw = 0; iw <= 2*ap->nw; iw++) {
                    int is = (iw-ap->nw)*ap->nj;
                    if ( ( (ix+is+pint)<0 ) || ( (ix+is+pint)>(nx-1) ) ) continue;
                    if ( ( (ix-is)<0 ) || ( (ix-is)>(nx-1) ) ) continue;
                    yy[i] += (xx[i+is+ip+pint] - xx[i-is]) * ap->flt[iw]  ;
                }
                yy[i]*= der?dir:1;
            }           
        }
    }
}

void se_allpass2(const se_allpass_t _ap, 
                   const int left,
                   const int der, 
                   const float* xx, 
                   float* yy)
{
    _se_allpass_t* ap = (_se_allpass_t*)_ap;
    int nx, ny, nz, i1, i2, ip, dir;

    nx = ap->nx;
    ny = ap->ny;
    nz = ap->nz;

    if (left) {
        i1=1; i2=nz;   ip=-nx*ny; dir=-1;
    } else {
        i1=0; i2=nz-1; ip= nx*ny; dir= 1;
    }

    memset(yy, 0, nx*ny*nz*sizeof(float));

    int iz;
    for (iz=i1; iz < i2; iz++) {
        int iy;
        for (iy=0; iy < ny; iy++) {
            int ix;
            for (ix = ap->nw*ap->nj; ix < nx-ap->nw*ap->nj; ix++) {
                int i = ix + nx * (iy + ny * iz);

                float p    = dir*ap->pp[i];
                int pint   = (int) p;
                float pdec = p - pint;

                if (der) {
                    se_aderfilter(ap->apf, pdec, ap->flt);
                } else {
                    se_passfilter(ap->apf, pdec, ap->flt);
                }

                int iw;
                for (iw = 0; iw <= 2*ap->nw; iw++) {
                    int is = (iw-ap->nw)*ap->nj;
                    if ( ( (ix+is+pint)<0 ) || ( (ix+is+pint)>(nx-1) ) ) continue;
                    if ( ( (ix-is)<0 ) || ( (ix-is)>(nx-1) ) ) continue;
                    yy[i] += (xx[i+is+ip] - xx[i-is]) * ap->flt[iw];
                }
                yy[i]*= der?dir:1;
            }
        }
    }
}

void se_allpass_filter_lookup(se_allpass_t _ap, double p, float *flt, int *min_i, int *max_i)
{
    _se_allpass_t* ap = (_se_allpass_t*)_ap;  //强制转换类型
    int en;
    
    en = (int)((p-ap->tbl_o)/ap->tbl_d);

    if ( (ap->flt_tbl) && (en>=0) && (en<(ap->tbl_n-1)) ) {
        double w1,w2;
        int i;
        int ns=2*ap->nw+1;
        float *flt_tbl_en=ap->flt_tbl + en*ns;
        float *flt_tbl_en_p=flt_tbl_en + ns;

        w2 = (p - ap->tbl_o)/ap->tbl_d - en;
        w1 = 1 - w2;
         
        (*min_i) = mini(ap->flt_mini[en],ap->flt_mini[en+1]);
        (*max_i) = maxi(ap->flt_maxi[en],ap->flt_maxi[en+1]);
        memset(flt,0,sizeof(float)*ns);
        for (i=(*min_i);i<=(*max_i);i++) {
            flt[i] = (float)(w1*flt_tbl_en[i] + w2*flt_tbl_en_p[i]);
        }
    } else {
        se_passfilter(ap->apf, (float)(p), flt);  //

        (*min_i) = 0;
        while  ( ((*min_i)<2*ap->nw) && 
                 (fabs(flt[(*min_i)])<1e-5) ) (*min_i)++;
        (*max_i) = 2*ap->nw;
        while  ( ((*max_i)>0) && 
                 (fabs(flt[(*max_i)])<1e-5) ) (*max_i)--;
    }

}

void se_allpass_shift_trace(se_allpass_t _ap,
                              double shift,
                              const float *tr, 
                              int n1, 
                              int i1, 
                              int w1, 
                              double *tr_shift)
{    
    _se_allpass_t* ap = (_se_allpass_t*)_ap;  //强制转换
    int st_i;
    int hw1=(w1-1)/2;
    int min_iw_flt, max_iw_flt; //flt filter

    // use the PWD to shift the trace. Given the all pass implementation, the dip
    // that needs to be used is twice the shift.

    /* se_passfilter(ap->apf, (float)(2*shift), ap->flt); */
    /* min_iw_flt=0; */
    /* while  ( (min_iw_flt<2*ap->nw) && (fabs(ap->flt[min_iw_flt])<1e-5) ) min_iw_flt++; */            
    /* max_iw_flt=2*ap->nw; */
    /* while  ( (max_iw_flt>0) && (fabs(ap->flt[max_iw_flt])<1e-5) ) max_iw_flt--; */

    se_allpass_filter_lookup(_ap, 2*shift, ap->flt, &min_iw_flt, &max_iw_flt);

    // the all pass filter is a convolution filter, so compute each sample by a 
    // weighted sum of the trace data  全通滤波是一个卷积，所以通过地震数据的一个加权叠加计算每一个采样点

    memset(tr_shift,0,w1*sizeof(double));

    for (st_i=0; st_i < w1; st_i++) {

        int iw;

        const int j1 = i1 - hw1 + st_i - ap->nw;    //j1=i1-hwl+st_i-ap->nw

        const int min_iw = maxi(min_iw_flt,-j1);
        const int max_iw = mini(max_iw_flt,n1-j1-1);

        double * const tr_shift_sample=tr_shift+st_i;
        const float * const tr_sample = tr + j1;

        for (iw=min_iw; iw<=max_iw; iw++) {            
            (*tr_shift_sample) += tr_sample[iw] * ap->flt[iw];
        }
    }
}

void se_beam_allpass1_L2(const se_allpass_t _ap,
                           const int left,
                           const int der,
                           const float* u,
                           const float* beam_p,
                           const int* beam_crd,
                           const int* beam_halfwid,
                           const int beam_total,
                           float* beam_yy)
{

    _se_allpass_t* ap = (_se_allpass_t*)_ap;
    se_allpass_t _apder = NULL;
    _se_allpass_t* apder = NULL;

    if (der) {
        _apder = se_allpass_init(ap->nw,ap->nj,ap->nx,ap->ny,ap->nz,ap->pp);
        apder = (_se_allpass_t*)_apder;
    }

    memset(beam_yy, 0, beam_total*sizeof(float));
 
    int i;
    for (i=0; i<beam_total; i++) {


        int i1 = beam_crd[3*i];
        int i2 = beam_crd[3*i+1];
        
        float yyi=0;
        
        int dir = left?-1:1;
        float p = dir*beam_p[i];
        int nx = ap->nx;
        int ny = ap->ny;

        int ip = dir*nx;

        se_passfilter(ap->apf, p, ap->flt);
 
        if (der) se_aderfilter(apder->apf, p,apder->flt);
            

        int j2off;
        for (j2off=-beam_halfwid[1]; j2off<=beam_halfwid[1]; j2off++) {
            int j2 = i2+j2off;
            if ( ((j2+dir)<0) || ((j2+dir)>=ny) ) continue;
            if ( (j2<0) || (j2>=ny) ) continue;
                    
            int j1;

            int poff=(int)(beam_p[i]*j2off);

            for (j1=i1-beam_halfwid[0]+poff; j1<=i1+beam_halfwid[0]+poff; j1++) {
                if ((j1<0) || (j1>=nx) ) continue;

                int ind = j1+j2*nx;
                int iw;
                float val=0, valder=0;
                for (iw = 0; iw <= 2*ap->nw; iw++) {
                    int is = (iw-ap->nw);
                    if ( ( (j1+is)<0 ) || ( (j1+is)>(nx-1) ) ) continue;
                    if ( ( (j1-is)<0 ) || ( (j1-is)>(nx-1) ) ) continue;

                    val += (u[ind+is+ip] - u[ind-is])* ap->flt[iw]  ;
                    if (der) {
                        //yyi+= 2*v* (der?dir:1)*(u[ind+is+ip] - u[ind-is])*  ((_se_allpass_t*)apder)->flt[iw]  ;
                        valder += (der?dir:1)*(u[ind+is+ip] - u[ind-is])*  apder->flt[iw]  ;
                    }
                }
                if (der) {
                    yyi+=2*val*valder;
                } else {
                    yyi+=val*val;
                }
            }
        }
        beam_yy[i]= yyi
            / ( (2*beam_halfwid[0]+1) *
                (2*beam_halfwid[1]+1) );
    }

    if (der) se_allpass_destroy(_apder);


}




