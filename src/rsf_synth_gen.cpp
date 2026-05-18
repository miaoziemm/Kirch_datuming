#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <math.h>
#include <stdlib.h>

int main(int argc,char**argv){
 se_par_init(argc,argv);
 char *in_f=se_get_par_str("input_file"), *sg_f=se_get_par_str("sgreen_file"), *rg_f=se_get_par_str("rgreen_file"), *m_f=se_get_par_str("model_file");
 int nt=se_get_par_int("nt"), nh=se_get_par_int("nh"), ns=se_get_par_int("ns"), ng=se_get_par_int("ng");
 float dt=se_get_par_float("dt"), dh=se_get_par_float("dh"), ds=se_get_par_float("ds");
 sep_t *in=sep_open(in_f,SEP_WRITE,0), *sg=sep_open(sg_f,SEP_WRITE,0), *rg=sep_open(rg_f,SEP_WRITE,0), *m=sep_open(m_f,SEP_WRITE,0);
 in->headers->ndim=3; in->headers->n[0]=nt; in->headers->n[1]=nh; in->headers->n[2]=ns; in->headers->d[0]=dt; in->headers->d[1]=dh; in->headers->d[2]=ds; in->headers->o[0]=0; in->headers->o[1]=0; in->headers->o[2]=0;
 float *d=alloc1float(nt*nh*ns);
 srand(7);
 for(int is=0;is<ns;is++) for(int ih=0;ih<nh;ih++) for(int it=0;it<nt;it++){
  float t=it*dt-0.04f, f0=25.0f; float a=M_PI*f0*t; a*=a; float w=(1-2*a)*expf(-a);
  float r=((rand()%1000)/500.0f-1.0f)*((rand()%100)<8?1.0f:0.0f);
  d[(is*nh+ih)*nt+it]=r*w;
 }
 se_fsio_write_float(in->data->io,d,nt*nh*ns);
 sg->headers->ndim=2; sg->headers->n[0]=ng; sg->headers->n[1]=ng; sg->headers->d[0]=ds; sg->headers->o[0]=0;
 rg->headers->ndim=2; rg->headers->n[0]=ng; rg->headers->n[1]=ng; rg->headers->d[0]=dh; rg->headers->o[0]=0;
 float *g=alloc1float(ng*ng);
 for(int i=0;i<ng;i++) for(int j=0;j<ng;j++) g[i*ng+j]=0.015f+0.5f*fabsf(i-j)*ds;
 se_fsio_write_float(sg->data->io,g,ng*ng); se_fsio_write_float(rg->data->io,g,ng*ng);
 m->headers->ndim=1; m->headers->n[0]=2; m->headers->d[0]=1.0f; m->headers->o[0]=1.0f; float mv[2]={0,0}; se_fsio_write_float(m->data->io,mv,2);
 sep_close(in); sep_close(sg); sep_close(rg); sep_close(m); free1float(d); free1float(g); return 0;
}
