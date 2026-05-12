#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>

static void read_axis_range(const char* parname, int n, int* k, int* e)
{
    if (!se_have_par(parname)) {
        return;
    }

    int size = 0;
    int32_t* vals = se_get_pararray_int(parname, &size);
    if (vals == NULL || size < 2) {
        ERROR(("Need %s=start,end", parname));
    }

    int start = (int)vals[0];
    int end = (int)vals[1];
    free(vals);

    if (start > end) {
        int tmp = start;
        start = end;
        end = tmp;
    }

    if (start < 0 || end >= n) {
        ERROR(("%s out of range [0,%d]", parname, n - 1));
    }

    *k = start;
    *e = end;
}

int main(int argc, char* argv[])
{
    se_par_init(argc, argv);

    int n1 = 0, n2 = 1, n3 = 1;
    int k1 = 0, k2 = 0, k3 = 0;
    int e1 = 0, e2 = 0, e3 = 0;
    char* in_f = NULL;
    char* out_f = NULL;
    sep_t* in = NULL;
    sep_t* out = NULL;

    if (!se_have_par("input_file")) ERROR(("Need input_file=")); else in_f = se_get_par_str("input_file");
    if (!se_have_par("output_file")) ERROR(("Need output_file=")); else out_f = se_get_par_str("output_file");

    in = sep_open(in_f, SEP_READ, 0);
    out = sep_open(out_f, SEP_WRITE, 0);

    int ndim = (int)sep_get_min_ndim(in);
    if (ndim < 1) ndim = 1;
    if (ndim > 3) {
        ERROR(("Only support up to 3D input"));
    }

    n1 = in->headers->n[0];
    n2 = (ndim >= 2) ? in->headers->n[1] : 1;
    n3 = (ndim >= 3) ? in->headers->n[2] : 1;

    e1 = n1 - 1;
    e2 = n2 - 1;
    e3 = n3 - 1;

    read_axis_range("n1_out", n1, &k1, &e1);
    if (ndim >= 2) {
        read_axis_range("n2_out", n2, &k2, &e2);
    } else if (se_have_par("n2_out")) {
        WARN(("Ignore n2_out because input is 1D"));
    }

    if (ndim >= 3) {
        read_axis_range("n3_out", n3, &k3, &e3);
    } else if (se_have_par("n3_out")) {
        WARN(("Ignore n3_out because input is not 3D"));
    }

    int m1 = e1 - k1 + 1;
    int m2 = e2 - k2 + 1;
    int m3 = e3 - k3 + 1;

    sep_copy_headers(out, in);

    char* label1 = sep_get_hdr(in, "label1", NULL);
    char* label2 = sep_get_hdr(in, "label2", NULL);
    char* label3 = sep_get_hdr(in, "label3", NULL);

    sep_set_axis(out, 0, m1, in->headers->o[0] + k1 * in->headers->d[0], in->headers->d[0], label1);
    if (ndim >= 2) {
        sep_set_axis(out, 1, m2, in->headers->o[1] + k2 * in->headers->d[1], in->headers->d[1], label2);
    }
    if (ndim >= 3) {
        sep_set_axis(out, 2, m3, in->headers->o[2] + k3 * in->headers->d[2], in->headers->d[2], label3);
    }

    if (label1) free(label1);
    if (label2) free(label2);
    if (label3) free(label3);

    sep_write_headers(out);

    float* trace = alloc1float(m1);
    for (int i3 = k3; i3 <= e3; ++i3) {
        for (int i2 = k2; i2 <= e2; ++i2) {
            int64_t idx = ((int64_t)i3 * n2 + i2) * n1 + k1;
            off_t byte_off = (off_t)(idx * (int64_t)sizeof(float));
            se_fsio_seek(in->data->io, byte_off);
            se_fsio_read_float(in->data->io, trace, (size_t)m1);
            se_fsio_write_float(out->data->io, trace, (size_t)m1);
        }
    }
    free1float(trace);

    sep_close(in);
    sep_close(out);

    se_par_destroy();
    return 0;
}