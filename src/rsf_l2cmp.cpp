#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <math.h>
int main(int argc,char**argv){
 se_par_init(argc,argv);
 sep_t *a=sep_open(se_get_par_str("a"),SEP_READ,0), *b=sep_open(se_get_par_str("b"),SEP_READ,0);
 int n=a->headers->n[0]*a->headers->n[1]*a->headers->n[2];
 float *x=alloc1float(n), *y=alloc1float(n); se_fsio_read_float(a->data->io,x,n); se_fsio_read_float(b->data->io,y,n);
 double num=0,den=0,mx=0; for(int i=0;i<n;i++){ double d=fabs((double)x[i]-y[i]); if(d>mx)mx=d; num+=d*d; den+=y[i]*y[i]; }
 printf("rel_l2=%e max_abs=%e\n", sqrt(num/(den+1e-30)), mx);
 return 0;
}
