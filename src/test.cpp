#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>



int main(int argc, char *argv[])
{
    se_par_init(argc, argv);

    int a;
    float b;
    char *c;

    sep_t *out_file;

    b=se_get_par_float("b");
    a=se_get_par_int("a");
    c=se_get_par_str("c");

    printf("a=%d, b=%f, c=%s\n", a, b, c);

    float *buff;
    buff=alloc1float(10);
    for(int i=0;i<10;i++) {
        buff[i]=i;
    }

    out_file=sep_open(c,SEP_WRITE,0);

    out_file->headers->ndim=1;
    out_file->headers->n[0]=10;
    out_file->headers->d[0]=1.0;
    out_file->headers->o[0]=0.0;

    se_fsio_write(out_file->data->io,buff,10*sizeof(float));
    sep_close(out_file);



    return 0;
}