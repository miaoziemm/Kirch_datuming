#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>

#include <math.h>
#include <stdlib.h>
#include <string.h>

#ifdef SE_USE_OMP
#include <omp.h>
#endif

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define MIN_NUM 1e-10f

#define MAX(a, b) ((a) > (b) ? (a) : (b))
#define MIN(a, b) ((a) < (b) ? (a) : (b))

int main(int argc, char *argv[])
{
    se_par_init(argc, argv);

    int ngx, nx, nz;
    int igx1, igx2;
    float dgx, ogx, dx, dz, ox, oz;
    sep_t *ttabel, *tgreen, *model;
    float datum = 0.0f;
    int idatum = 0;
    int should_datum = 0;
    char *ttabel_f = NULL, *tgreen_f = NULL, *model_f = NULL;
    

    if (!se_have_par("ttabel_file"))ERROR(("Need ttabel_file="));else ttabel_f = se_get_par_str("ttabel_file");
    if (!se_have_par("tgreen_file"))ERROR(("Need tgreen_file="));else tgreen_f = se_get_par_str("tgreen_file");
    if (!se_have_par("model_file"))ERROR(("Need model_file="));else model_f = se_get_par_str("model_file");
    ttabel = sep_open(ttabel_f, SEP_READ, 0);
    tgreen = sep_open(tgreen_f, SEP_WRITE, 0);
    model = sep_open(model_f, SEP_READ, 0);

    if(ttabel->headers->ndim < 3) ERROR(("Need 3D data volume for ttabel= (z,x,gx)"));

    nz = ttabel->headers->n[0];
    nx = ttabel->headers->n[1];
    dz = ttabel->headers->d[0];
    dx = ttabel->headers->d[1];
    oz = ttabel->headers->o[0];
    ox = ttabel->headers->o[1];

    ngx = ttabel->headers->n[2];
    dgx = ttabel->headers->d[2];
    ogx = ttabel->headers->o[2];

    tgreen->headers->ndim = 2;
    tgreen->headers->n[0] = ngx;
    tgreen->headers->d[0] = dgx;
    tgreen->headers->o[0] = ogx;

    if (sep_have_hdr_int(model, "should_datum")) {
        should_datum = sep_get_hdr_int(model, "should_datum", 0);
    } else {
        ERROR(("Need should_datum in model_file header"));
    }
    INFO(("Read should_datum=%d from %s", should_datum, model_f));

    datum = (float)(model->headers->o[0] + should_datum * model->headers->d[0]);
    idatum = (int)((datum-oz)/dz+0.5);
    if (idatum < 0 || idatum >= nz) ERROR(("Datum out of range."));
    // tgreen保存的数据为ngx乘ngx的走时数据,从0平面到datum深度处的走时数据
    float ***ttabel_data = alloc3float(nz, nx, ngx);
    se_fsio_read_float(ttabel->data->io, ttabel_data[0][0], nz*nx*ngx);
    float **tgreen_data = alloc2float(ngx, ngx);
    for (igx1=0; igx1 < ngx; igx1++) {
        for (igx2=0; igx2 < ngx; igx2++) {
            float gx1 = ogx + igx1*dgx;
            float gx2 = ogx + igx2*dgx;
            int ix1 = (int)((gx1-ox)/dx+0.5);
            int ix2 = (int)((gx2-ox)/dx+0.5);
            if (ix1 < 0 || ix1 >= nx || ix2 < 0 || ix2 >= nx) ERROR(("Receiver table too small."));
            tgreen_data[igx1][igx2] = ttabel_data[igx1][ix2][idatum];
            
        }
    }


    se_fsio_write_float(tgreen->data->io, tgreen_data[0], ngx*ngx);

    free2float(tgreen_data);
    free3float(ttabel_data);

    sep_close(model);
    sep_close(ttabel);
    sep_close(tgreen);

    se_par_destroy();
    return 0;
}