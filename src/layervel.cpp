#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static char* build_output_name(const char* base, int layer_index)
{
    const char* ext = strrchr(base, '.');
    if (ext != NULL && strcmp(ext, ".rsf") == 0) {
        int needed = snprintf(NULL, 0, "%.*s_%d%s", (int)(ext - base), base, layer_index, ext);
        char* out = (char*)malloc((size_t)needed + 1);
        if (out == NULL) {
            ERROR(("Out of memory"));
        }
        snprintf(out, (size_t)needed + 1, "%.*s_%d%s", (int)(ext - base), base, layer_index, ext);
        return out;
    }

    int needed = snprintf(NULL, 0, "%s_%d", base, layer_index);
    char* out = (char*)malloc((size_t)needed + 1);
    if (out == NULL) {
        ERROR(("Out of memory"));
    }
    snprintf(out, (size_t)needed + 1, "%s_%d", base, layer_index);
    return out;
}

int main(int argc, char* argv[])
{
    se_par_init(argc, argv);

    int nlayers = 0;
    int overlap = 0;
    char* in_f = NULL;
    char* out_f_base = NULL;
    sep_t* in = NULL;
    sep_t* out = NULL;

    if (!se_have_par("input_file")) ERROR(("Need input_file=")); else in_f = se_get_par_str("input_file");
    if (!se_have_par("nlayers")) ERROR(("Need nlayers=")); else nlayers = se_get_par_int("nlayers");
    if (!se_have_par("overlap")) ERROR(("Need overlap=")); else overlap = se_get_par_int("overlap");
    if (!se_have_par("output_file_base")) ERROR(("Need output_file_base=")); else out_f_base = se_get_par_str("output_file_base"); // 可能会输出很多文件，不能直接指定输出文件名

    if (nlayers <= 0) {
        ERROR(("Invalid nlayers=%d", nlayers));
    }
    if (overlap < 0) {
        ERROR(("Invalid overlap=%d", overlap));
    }

    in = sep_open(in_f, SEP_READ, 0);

    int ndim = (int)sep_get_min_ndim(in);
    if (ndim < 1) ndim = 1;
    if (ndim > 3) {
        ERROR(("Only support up to 3D input"));
    }

    int n1 = in->headers->n[0];
    int n2 = (ndim >= 2) ? in->headers->n[1] : 1;
    int n3 = (ndim >= 3) ? in->headers->n[2] : 1;

    if (n1 <= 0) {
        ERROR(("Invalid n1=%d", n1));
    }

    int64_t total = (int64_t)n1 + (int64_t)(nlayers - 1) * (int64_t)overlap;
    int base = (int)(total / nlayers);
    int rem = (int)(total % nlayers);

    if (base <= 0) {
        ERROR(("nlayers=%d too large for n1=%d", nlayers, n1));
    }
    if (overlap >= base) {
        ERROR(("overlap=%d too large for n1=%d and nlayers=%d", overlap, n1, nlayers));
    }

    char* label1 = sep_get_hdr(in, "label1", NULL);

    int start = 0;
    for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
        int m1 = base + ((ilayer < rem) ? 1 : 0);
        int end = start + m1 - 1;
        if (end >= n1) {
            ERROR(("Layer %d out of range: end=%d n1=%d", ilayer + 1, end, n1));
        }

        char* out_f = build_output_name(out_f_base, ilayer + 1);
        out = sep_open(out_f, SEP_WRITE, 0);
        sep_copy_headers(out, in);
        out->headers->ndim = in->headers->ndim;
        sep_set_axis(out, 0, m1, in->headers->o[0], in->headers->d[0], label1);
        int should_datum = m1 - overlap - 1;
        int should_datum_abs = end - overlap;
        sep_set_header_int(out, "should_datum", should_datum);
        sep_set_header_int(out, "should_datum_abs", should_datum_abs);
        sep_write_headers(out);

        float* trace = alloc1float(m1);
        for (int i3 = 0; i3 < n3; ++i3) {
            for (int i2 = 0; i2 < n2; ++i2) {
                int64_t idx = ((int64_t)i3 * n2 + i2) * n1 + start;
                off_t byte_off = (off_t)(idx * (int64_t)sizeof(float));
                se_fsio_seek(in->data->io, byte_off);
                se_fsio_read_float(in->data->io, trace, (size_t)m1);
                se_fsio_write_float(out->data->io, trace, (size_t)m1);
            }
        }
        free1float(trace);
        sep_close(out);
        free(out_f);

        start = end - overlap + 1;
    }

    if (label1) free(label1);

    sep_close(in);
    se_par_destroy();

    return 0;
}