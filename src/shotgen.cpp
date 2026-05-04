#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>
#include <SERECKIRCH/include/se_reckirch.h>

int main(int argc, char* argv[])
{
    int nshot = 0;
    float ox=0.0f, dx=0.0f, oy=0.0f, dy=0.0f, zs=0.0f;
    char *shotfile=NULL;
    sep_t *out = nullptr;

    se_par_init(argc, argv);

    if (!se_have_par("shotfile")) ERROR(("Need shotfile=")); else shotfile = se_get_par_str("shotfile");
    if (!se_have_par("nshot")) ERROR(("Need nshot=")); else nshot = se_get_par_int("nshot");
    if (!se_have_par("ox")) ox = 0.0f; else ox = se_get_par_float("ox");
    if (!se_have_par("dx")) dx = 0.0f; else dx = se_get_par_float("dx");
    if (!se_have_par("oy")) oy = 0.0f; else oy = se_get_par_float("oy");
    if (!se_have_par("dy")) dy = 0.0f; else dy = se_get_par_float("dy");
    if (!se_have_par("zs")) zs = 0.0f; else zs = se_get_par_float("zs");

    out=sep_open(shotfile, SEP_WRITE, 0);
    out->headers->ndim = 2;
    out->headers->n[0] = 3;
    out->headers->n[1] = nshot;
    

    for(int is=0; is<nshot; is++)
    {
        float s[3];
        s[2] = ox + is * dx;
        s[1] = oy + is * dy;
        s[0] = zs;
        se_fsio_write_float(out->data->io, s, 3);
    }

    sep_close(out);
    return 0;
}

   