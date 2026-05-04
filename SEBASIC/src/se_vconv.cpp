#include "../include/se_vconv.h"

// time to time

void vconv_time_int_to_rms(const axa_t* a, const float *vint, float *vrms) {
    int64_t i;

    double sum_v2=0.0, vm=vint[0];

    vrms[0]=vint[0];
    for(i=1;i<a->n;i++) {
        const double d=vint[i]-vm;
        vm=vint[i];

        // exact integration of piecewise-linear vint
        sum_v2+=vm*vm + d*d/3.0 + vm*d;
        vrms[i]=(float)sqrt(sum_v2/i);
    }
}

void vconv_time_int_to_avg(const axa_t *a, const float *vint, float *vavg) {
    int64_t i;

    double z=0.0, vm=vint[0];

    vavg[0] = vint[0];
    for(i=1;i<a->n;i++) {
        // exact integration of piecewise-linear vint
        z += 0.5*(vint[i]+vm);
        vm=vint[i];
        vavg[i] = (float)(z/i);
    }
}

void vconv_time_rms_to_int(const axa_t* a, const float *vrms, float *vint) {
    int64_t i;

    vint[0]=vrms[0];
    for(i=1;i<a->n;i++) {
        double vrms2m, vrms2p, vint2;

        if (i < a->n-1) { // centered difference
            vrms2m = vrms[i-1]*vrms[i-1];
            vrms2p = vrms[i+1]*vrms[i+1];
            vint2 = i*0.5*(vrms2p-vrms2m)+vrms[i]*vrms[i];
        }
        else { // backward difference
            vrms2m = vrms[i-1]*vrms[i-1];
            vrms2p = vrms[i]*vrms[i];
            vint2 = i*(vrms2p-vrms2m)+vrms[i]*vrms[i];
        }

        if (vint2>0.0) vint[i]=(float)sqrt(vint2);
        else   vint[i]=vint[i-1];
    }
}

void vconv_time_rms_to_avg(const axa_t *a, const float *vrms, float *vavg) {
    int64_t i;
    double z=0.0, vm=vrms[0];

    vavg[0]=vrms[0];
    for(i=1;i<a->n;i++) {
        double vrms2m, vrms2p, vint2, vi;

        if (i < a->n-1) { // centered difference
            vrms2m = vrms[i-1]*vrms[i-1];
            vrms2p = vrms[i+1]*vrms[i+1];
            vint2 = i*0.5*(vrms2p-vrms2m)+vrms[i]*vrms[i];
        }
        else { // backward difference
            vrms2m = vrms[i-1]*vrms[i-1];
            vrms2p = vrms[i]*vrms[i];
            vint2 = i*(vrms2p-vrms2m)+vrms[i]*vrms[i];
        }

        if (vint2>0.0) vi=sqrt(vint2);
        else vi=vm;

        z += 0.5*(vi+vm)*a->d;
        vavg[i] = (float)(z/(i*a->d));

        vm=vi;
    }
}

void vconv_time_avg_to_int(const axa_t *a, const float *vavg, float *vint) {
    int64_t i;

    vint[0]=vavg[0];

    for(i=1;i<a->n;i++) {
        double vi;
        if (i < a->n-1) { // centered difference
            vi=0.5*((i+1)*vavg[i+1]-(i-1)*vavg[i-1]);
        }
        else { // backward difference
            vi=(i*vavg[i]-(i-1)*vavg[i-1]);
        }

        if (vi > 0.0) vint[i]=(float)vi;
        else vint[i]=vint[i-1];
    }
}

void vconv_time_avg_to_rms(const axa_t *a, const float *vavg, float *vrms) {
    int64_t i;

    double sum_v2=0.0, vm=vavg[0];

    vrms[0]=vavg[0];

    for(i=1;i<a->n;i++) {
        double vi, d;
        if (i < a->n-1) { // centered difference
            vi=0.5*((i+1)*vavg[i+1]-(i-1)*vavg[i-1]);
        }
        else { // backward difference
            vi=(i*vavg[i]-(i-1)*vavg[i-1]);
        }

        if (vi < 0.0) vi=vm;
        d = vi-vm;

        sum_v2+=vm*vm + d*d/3.0 + vm*d;
        vrms[i]=(float)sqrt(sum_v2/i);
        vm=vi;
    }
}

// depth to depth

void vconv_depth_int_to_rms(const axa_t* a, const float *vint, float *vrms) {
    int64_t i;

    double num=0.0, denom=0.0, vm=vint[0];

    vrms[0]=vint[0];
    for(i=1;i<a->n;i++) {
        num += vint[i] + vm; // exact integration of piecewise-linear vint
        denom += 1.0/vint[i] + 1.0/vm; // treat 1/v as piecewise-linear
        vm=vint[i];

        vrms[i]=(float)sqrt(num/denom);
    }
}

void vconv_depth_int_to_avg(const axa_t *a, const float *vint, float *vavg) {
    int64_t i;

    // vavg[i]^2 + (i-1)*vint[i]*vavg[i] - i*vint[i]*vavg[i-1] = 0
    vavg[0]=vint[0];
    for(i=1;i<a->n;i++) {
        double v=(1.0-i)*vint[i];
        vavg[i]=(float)(0.5*(v+sqrt(v*v+4*i*vint[i]*vavg[i-1])));
    }
}

void vconv_depth_rms_to_int(const axa_t* a, const float *vrms, float *vint) {
    int64_t i;

    double num=0.0, denom=0.0;
    vint[0]=vrms[0];
    
    for(i=1;i<a->n;i++) {
        const double vrms2 = vrms[i]*vrms[i];
        const double b = vint[i-1]+num-vrms2*(denom+1.0/vint[i-1]);

        vint[i] = (float)(0.5*(-b + sqrt(b*b+4*vrms2)));

        num+=(vint[i]+vint[i-1]);
        denom = num/vrms2;
    }
}

void vconv_depth_rms_to_avg(const axa_t *a, const float *vrms, float *vavg) {
    int64_t i;

    double num=0.0, denom=0.0, vm=vrms[0];
    vavg[0]=vrms[0];

    for(i=1;i<a->n;i++) {
        const double vrms2 = vrms[i]*vrms[i];
        const double b = vm+num-vrms2*(denom+1.0/vm);
        const double vi = 0.5*(-b + sqrt(b*b+4*vrms2));
        const double v=(1.0-i)*vi;

        vavg[i]=(float)(0.5*(v+sqrt(v*v+4*i*vi*vavg[i-1])));

        num+=(vi+vm);
        denom = num/vrms2;
        vm=vi;
    }
}

void vconv_depth_avg_to_int(const axa_t *a, const float *vavg, float *vint) {
    int64_t i;

    // vavg[i]^2 + (i-1)*vint[i]*vavg[i] - i*vint[i]*vavg[i-1] = 0
    vint[0]=vavg[0];
    for(i=1;i<a->n;i++) {
        vint[i]=(float)(vavg[i]*vavg[i]/(vavg[i]-i*(vavg[i]-vavg[i-1])));
    }
}

void vconv_depth_avg_to_rms(const axa_t *a, const float *vavg, float *vrms) {
    int64_t i;

    double num=0.0, denom=0.0, vm=vavg[0];

    vrms[0]=vavg[0];
    for(i=1;i<a->n;i++) {
        const double vi = vavg[i]*vavg[i]/(vavg[i]-i*(vavg[i]-vavg[i-1]));

        num += vi + vm;
        denom += 1.0/vi + 1.0/vm;

        vrms[i]=(float)sqrt(num/denom);
        vm=vi;
    }
}

// axis conversion based on a single trace

void vconv_depth_axis_from_time(const axa_t *at, const float *vint_t, axa_t *az) {
    int64_t i;
    double t_4=0.25*at->d, min_d=(vint_t[0]+vint_t[1])*t_4, d=min_d;

    for(i=1;i<at->n;i++) {
        double di=(vint_t[i-1]+vint_t[i])*t_4;
        if (di<min_d) {
            min_d=di;
        }
        d+=di;
    }

    az->n=(int64_t)ceil(d/min_d);
    az->d=min_d;
    az->o=0.0;
}

void vconv_time_axis_from_depth(const axa_t *az, const float *vint_z, axa_t *at) {
    int64_t i;
    double vm=1.0/vint_z[1], min_t=az->d*(1.0/vint_z[0]+vm), t=min_t;

    for(i=1;i<az->n;i++) {
        double v=1.0/vint_z[i], ti=az->d*(v+vm);

        if (ti<min_t) {
            min_t=ti;
        }
        t+=ti;
        vm=v;
    }

    at->n=(int64_t)ceil(t/min_t);
    at->d=min_t;
    at->o=0.0;
}

// time to depth
int vconv_time_to_depth_int_to_int(const axa_t *at, const float *vint_t, const axa_t *az, float *vint_z, float *work) {

    int64_t i, j;
    double t_2=0.5*at->d, d_i_1=0.0, d_i;
    int64_t jmin=0, jmax=az->n;
    (void)work;

    {
        d_i=d_i_1+t_2*0.5*(vint_t[0]+vint_t[1]);
        vint_z[0]=vint_t[0];
        jmin=1;
        jmax=mini64(az->n,axa_lclosest_indx(az,d_i)+1);

        for(j=jmin;j<jmax;j++) {
            float q=(float)((d_i-axa_x_from_idx(az, j))/(d_i-d_i_1));
            vint_z[j]=vint_t[0]*q+vint_t[1]*(1.0f-q);
        }

        d_i_1=d_i;
    }

    for(i=1;i<at->n-2;i++) {
        d_i=d_i_1+t_2*0.5*(vint_t[i]+vint_t[i+1]);

        // check bounds
        jmin=maxi64(0,axa_lclosest_indx(az,d_i_1)+1);
        jmax=mini64(az->n,axa_lclosest_indx(az,d_i)+1);

        for(j=jmin;j<jmax;j++) {
            float q=(float)((d_i-axa_x_from_idx(az, j))/(d_i-d_i_1));
            vint_z[j]=vint_t[i]*q+vint_t[i+1]*(1.0f-q);
        }
        d_i_1=d_i;
    }

    for(j=jmax;j<az->n;j++) {
        vint_z[j]=vint_t[at->n-1];
    }

    if (jmax < az->n) {
        return 1;
    }
    return 0;
}

int vconv_time_to_depth_int_to_rms(const axa_t *at, const float *vint_t, const axa_t *az, float *vrms_z, float *work) {
    int ret;
    (void)work;

    ret=vconv_time_to_depth_int_to_int(at, vint_t, az, vrms_z, NULL);
    vconv_depth_int_to_rms(az, vrms_z, vrms_z);

    return ret;
}

int vconv_time_to_depth_int_to_avg(const axa_t *at, const float *vint_t, const axa_t *az, float *vavg_z, float *work) {
    int ret;
    (void)work;

    ret=vconv_time_to_depth_int_to_int(at, vint_t, az, vavg_z, NULL);
    vconv_depth_int_to_avg(az, vavg_z, vavg_z);

    return ret;
}

int vconv_time_to_depth_rms_to_int(const axa_t *at, const float *vrms_t, const axa_t *az, float *vint_z, float *work) {
    int ret;

    vconv_time_rms_to_int(at, vrms_t, work);
    ret=vconv_time_to_depth_int_to_int(at, work, az, vint_z, NULL);

    return ret;
}

int vconv_time_to_depth_rms_to_rms(const axa_t *at, const float *vrms_t, const axa_t *az, float *vrms_z, float *work) {
    int ret;

    vconv_time_rms_to_int(at, vrms_t, work);
    ret=vconv_time_to_depth_int_to_rms(at, work, az, vrms_z, NULL);

    return ret;
}

int vconv_time_to_depth_rms_to_avg(const axa_t *at, const float *vrms_t, const axa_t *az, float *vavg_z, float *work) {
    int ret;

    vconv_time_rms_to_int(at, vrms_t, work);
    ret=vconv_time_to_depth_int_to_avg(at, work, az, vavg_z, NULL);

    return ret;
}

int vconv_time_to_depth_avg_to_int(const axa_t *at, const float *vavg_t, const axa_t *az, float *vint_z, float *work) {
    int ret;

    vconv_time_avg_to_int(at, vavg_t, work);
    ret=vconv_time_to_depth_int_to_int(at, work, az, vint_z, NULL);

    return ret;
}

int vconv_time_to_depth_avg_to_rms(const axa_t *at, const float *vavg_t, const axa_t *az, float *vrms_z, float *work) {
    int ret;

    vconv_time_avg_to_int(at, vavg_t, work);
    ret=vconv_time_to_depth_int_to_rms(at, work, az, vrms_z, NULL);

    return ret;
}

int vconv_time_to_depth_avg_to_avg(const axa_t *at, const float *vavg_t, const axa_t *az, float *vavg_z, float *work) {
    int ret;

    vconv_time_avg_to_int(at, vavg_t, work);
    ret=vconv_time_to_depth_int_to_avg(at, work, az, vavg_z, NULL);

    return ret;
}

int vconv_time_to_depth(const axa_t *at, const float *vint_t, float *z) {
    int64_t n=at->n, i;
    double d_2=at->d*0.5;

    z[0] = 0.0f;
    for(i=1;i<n;i++) {
        z[i]=z[i-1]+(float)(0.5*(vint_t[i] + vint_t[i-1])*d_2);
    }

    return 0;
}

int vconv_data_time_to_depth(const axa_t *at, const float *f_t, const float *vint_z, const axa_t *az, float *f_z, float *work) {
    int ret=0;
    int64_t i;

    // compute t(z)
    if (vconv_depth_to_time(az, vint_z, work)) { // work=t(z)
        ret=1;
    }

    // compute output f(z)
    for(i=0;i<az->n;i++) {
        f_z[i]=(float)se_interp_axa(at, f_t, work[i]);
    }
    return ret;
}

// depth to time
int vconv_depth_to_time_int_to_int(const axa_t *az, const float *vint_z, const axa_t *at, float *vint_t, float *work) {
    int64_t i, j;
    double d2=2.0*az->d, t_i_1=0.0, t_i;
    int64_t jmin=0, jmax=at->n;
    (void)work;

    {
        t_i=t_i_1+d2/(0.5*(vint_z[0]+vint_z[1]));
        vint_t[0] = vint_z[0];
        jmin=1;
        jmax=mini64(at->n,axa_lclosest_indx(at,t_i)+1);
        for(j=jmin;j<jmax;j++) {
            float q=(float)((t_i-axa_x_from_idx(at, j))/(t_i-t_i_1));
            vint_t[j]=vint_z[0]*q+vint_z[1]*(1.0f-q);
        }
        t_i_1=t_i;
    }


    for(i=1;i<az->n-2;i++) {
        t_i=t_i_1+d2/(0.5*(vint_z[i]+vint_z[i+1]));

        // check bounds
        jmin=maxi64(0,axa_lclosest_indx(at,t_i_1));
        jmax=mini64(at->n,axa_lclosest_indx(at,t_i)+1);

        for(j=jmin;j<jmax;j++) {
            float q=(float)((t_i-axa_x_from_idx(at, j))/(t_i-t_i_1));
            vint_t[j]=vint_z[i]*q+vint_z[i+1]*(1.0f-q);
        }

        t_i_1=t_i;
    }

    for(j=jmax;j<at->n;j++) {
        vint_t[j]=vint_z[az->n-1];
    }

    if (jmax < at->n) {
        return 1;
    }
    return 0;
}

int vconv_depth_to_time_int_to_rms(const axa_t *az, const float *vint_z, const axa_t *at, float *vrms_t, float *work) {
    int ret;
    (void)work;

    ret=vconv_depth_to_time_int_to_int(az, vint_z, at, vrms_t, NULL);
    vconv_time_int_to_rms(at, vrms_t, vrms_t);

    return ret;
}

int vconv_depth_to_time_int_to_avg(const axa_t *az, const float *vint_z, const axa_t *at, float *vavg_t, float *work) {
    int ret;
    (void)work;

    ret=vconv_depth_to_time_int_to_int(az, vint_z, at, vavg_t, NULL);
    vconv_time_int_to_avg(at, vavg_t, vavg_t);

    return ret;
}

int vconv_depth_to_time_rms_to_int(const axa_t *az, const float *vrms_z, const axa_t *at, float *vint_t, float *work) {
    int ret;

    vconv_depth_rms_to_int(az, vrms_z, work);
    ret=vconv_depth_to_time_int_to_int(az, work, at, vint_t, NULL);

    return ret;
}

int vconv_depth_to_time_rms_to_rms(const axa_t *az, const float *vrms_z, const axa_t *at, float *vrms_t, float *work) {
    int ret;

    vconv_depth_rms_to_int(az, vrms_z, work);
    ret=vconv_depth_to_time_int_to_rms(az, work, at, vrms_t, NULL);

    return ret;
}

int vconv_depth_to_time_rms_to_avg(const axa_t *az, const float *vrms_z, const axa_t *at, float *vavg_t, float *work) {
    int ret;

    vconv_depth_rms_to_int(az, vrms_z, work);
    ret=vconv_depth_to_time_int_to_avg(az, work, at, vavg_t, NULL);

    return ret;
}

int vconv_depth_to_time_avg_to_int(const axa_t *az, const float *vavg_z, const axa_t *at, float *vint_t, float *work) {
    int ret;

    vconv_depth_avg_to_int(az, vavg_z, work);
    ret=vconv_depth_to_time_int_to_int(az, work, at, vint_t, NULL);

    return ret;
}

int vconv_depth_to_time_avg_to_rms(const axa_t *az, const float *vavg_z, const axa_t *at, float *vrms_t, float *work) {
    int ret;

    vconv_depth_avg_to_int(az, vavg_z, work);
    ret=vconv_depth_to_time_int_to_rms(az, work, at, vrms_t, NULL);

    return ret;
}

int vconv_depth_to_time_avg_to_avg(const axa_t *az, const float *vavg_z, const axa_t *at, float *vavg_t, float *work) {
    int ret;

    vconv_depth_avg_to_int(az, vavg_z, work);
    ret=vconv_depth_to_time_int_to_avg(az, work, at, vavg_t, NULL);

    return ret;
}

int vconv_depth_to_time(const axa_t *az, const float *vint_z, float *t) {
    int64_t n=az->n, i;
    double vm=1.0/vint_z[0];

    t[0]=0.0f;
    for(i=1;i<n;i++) {
        double v=1.0/vint_z[i];
        t[i]=t[i-1]+(float)(az->d*(v+vm));
        vm=v;
    }

    return 0;
}

int vconv_data_depth_to_time(const axa_t *az, const float *f_z, const float *vint_t, const axa_t *at, float *f_t, float *work) {
    int ret=0;
    int64_t i;

    // compute z(t)
    if (vconv_time_to_depth(at, vint_t, work)) { // work=z(t)
        ret=1;
    }

    // compute output f(t)
    for(i=0;i<at->n;i++) {
        f_t[i]=(float)se_interp_axa(az, f_z, work[i]);
    }
    return ret;
}
