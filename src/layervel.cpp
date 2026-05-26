#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>

typedef struct {
    int start;
    int end;
    int should_datum_abs;   /* -1 means no datum */
} layer_range_t;

static int is_end_token(const char* s)
{
    return s &&
           (strcmp(s, "END") == 0 ||
            strcmp(s, "End") == 0 ||
            strcmp(s, "end") == 0);
}

static char* build_output_name(const char* base, int layer_index)
{
    const char* ext = strrchr(base, '.');

    if (ext != NULL && strcmp(ext, ".rsf") == 0) {
        int needed = snprintf(NULL, 0, "%.*s_%d%s",
                              (int)(ext - base), base, layer_index, ext);
        char* out = (char*)malloc((size_t)needed + 1);
        if (out == NULL) ERROR(("Out of memory"));

        snprintf(out, (size_t)needed + 1, "%.*s_%d%s",
                 (int)(ext - base), base, layer_index, ext);
        return out;
    }

    int needed = snprintf(NULL, 0, "%s_%d", base, layer_index);
    char* out = (char*)malloc((size_t)needed + 1);
    if (out == NULL) ERROR(("Out of memory"));

    snprintf(out, (size_t)needed + 1, "%s_%d", base, layer_index);
    return out;
}

static int next_token(char** p, char* tok, int tok_size)
{
    char* s = *p;
    int n = 0;

    while (*s) {
        if (isspace((unsigned char)*s) || *s == ',' || *s == ';') {
            s++;
        } else if ((unsigned char)s[0] == 0xEF &&
                   (unsigned char)s[1] == 0xBC &&
                   (unsigned char)s[2] == 0x8C) {
            s += 3; /* Chinese comma: ， */
        } else {
            break;
        }
    }

    if (*s == '\0' || *s == '#') {
        *p = s;
        return 0;
    }

    while (*s) {
        if (isspace((unsigned char)*s) || *s == ',' || *s == ';' || *s == '#') {
            break;
        }
        if ((unsigned char)s[0] == 0xEF &&
            (unsigned char)s[1] == 0xBC &&
            (unsigned char)s[2] == 0x8C) {
            break;
        }

        if (n < tok_size - 1) {
            tok[n++] = *s;
        }
        s++;
    }

    tok[n] = '\0';
    *p = s;
    return 1;
}

static int parse_int_token(const char* tok, const char* name, int line_no)
{
    char* endp = NULL;
    long v = strtol(tok, &endp, 10);

    if (endp == tok || *endp != '\0') {
        ERROR(("Invalid %s at line %d: %s", name, line_no, tok));
    }

    if (v < 0 || v > 2147483647L) {
        ERROR(("Invalid %s at line %d: %s", name, line_no, tok));
    }

    return (int)v;
}

static layer_range_t* read_layer_txt(const char* fname, int n1, int* nlayers_out)
{
    FILE* fp = fopen(fname, "r");
    if (fp == NULL) {
        ERROR(("Cannot open layer_txt=%s", fname));
    }

    int cap = 16;
    int nlayer = 0;
    layer_range_t* layers = (layer_range_t*)malloc((size_t)cap * sizeof(layer_range_t));
    if (layers == NULL) ERROR(("Out of memory"));

    char line[4096];
    int line_no = 0;

    while (fgets(line, sizeof(line), fp) != NULL) {
        line_no++;

        char* p = line;
        char t1[256], t2[256], t3[256];

        if (!next_token(&p, t1, sizeof(t1))) continue;
        if (!next_token(&p, t2, sizeof(t2))) {
            ERROR(("Line %d needs at least start,end", line_no));
        }

        if (nlayer >= cap) {
            cap *= 2;
            layers = (layer_range_t*)realloc(layers, (size_t)cap * sizeof(layer_range_t));
            if (layers == NULL) ERROR(("Out of memory"));
        }

        layers[nlayer].start = parse_int_token(t1, "start", line_no);

        if (is_end_token(t2)) {
            layers[nlayer].end = n1 - 1;
            layers[nlayer].should_datum_abs = -1;
        } else {
            layers[nlayer].end = parse_int_token(t2, "end", line_no);

            if (!next_token(&p, t3, sizeof(t3))) {
                ERROR(("Line %d needs should_datum_abs, except when end is END", line_no));
            }

            if (is_end_token(t3)) {
                layers[nlayer].should_datum_abs = -1;
            } else {
                layers[nlayer].should_datum_abs =
                    parse_int_token(t3, "should_datum_abs", line_no);
            }
        }

        nlayer++;
    }

    fclose(fp);

    if (nlayer <= 0) {
        ERROR(("No valid layers found in %s", fname));
    }

    for (int i = 0; i < nlayer; i++) {
        if (layers[i].start < 0 || layers[i].start >= n1) {
            ERROR(("Layer %d invalid start=%d, n1=%d", i + 1, layers[i].start, n1));
        }

        if (layers[i].end < layers[i].start || layers[i].end >= n1) {
            ERROR(("Layer %d invalid range: start=%d end=%d n1=%d",
                   i + 1, layers[i].start, layers[i].end, n1));
        }

        if (i > 0 && layers[i].start > layers[i - 1].end) {
            ERROR(("Layer %d start=%d must be <= previous end=%d",
                   i + 1, layers[i].start, layers[i - 1].end));
        }

        if (layers[i].should_datum_abs >= 0) {
            if (layers[i].should_datum_abs < layers[i].start ||
                layers[i].should_datum_abs > layers[i].end) {
                ERROR(("Layer %d invalid should_datum_abs=%d, valid range is [%d,%d]",
                       i + 1, layers[i].should_datum_abs,
                       layers[i].start, layers[i].end));
            }
        }
    }

    if (layers[0].start != 0) {
        ERROR(("First layer must start from 0, but start=%d", layers[0].start));
    }

    if (layers[nlayer - 1].end != n1 - 1) {
        ERROR(("Last layer must reach model bottom n1-1=%d. Use END/end/End for the last layer.",
               n1 - 1));
    }

    *nlayers_out = nlayer;
    return layers;
}

int main(int argc, char* argv[])
{
    se_par_init(argc, argv);

    int nlayers = 0;
    int overlap = 0;
    char* in_f = NULL;
    char* out_f_base = NULL;
    char* layer_txt = NULL;

    sep_t* in = NULL;
    sep_t* out = NULL;

    if (!se_have_par("input_file")) ERROR(("Need input_file="));
    else in_f = se_get_par_str("input_file");

    if (!se_have_par("output_file_base")) ERROR(("Need output_file_base="));
    else out_f_base = se_get_par_str("output_file_base");

    if (se_have_par("layer_txt")) {
        layer_txt = se_get_par_str("layer_txt");
    }

    if (!layer_txt) {
        if (!se_have_par("nlayers")) ERROR(("Need nlayers= when layer_txt is not provided"));
        else nlayers = se_get_par_int("nlayers");

        if (!se_have_par("overlap")) ERROR(("Need overlap= when layer_txt is not provided"));
        else overlap = se_get_par_int("overlap");
    }

    in = sep_open(in_f, SEP_READ, 0);

    int ndim = (int)sep_get_min_ndim(in);
    if (ndim < 1) ndim = 1;
    if (ndim > 3) ERROR(("Only support up to 3D input"));

    int n1 = in->headers->n[0];
    int n2 = (ndim >= 2) ? in->headers->n[1] : 1;
    int n3 = (ndim >= 3) ? in->headers->n[2] : 1;

    if (n1 <= 0) ERROR(("Invalid n1=%d", n1));

    layer_range_t* layers = NULL;

    if (layer_txt) {
        layers = read_layer_txt(layer_txt, n1, &nlayers);
    } else {
        if (nlayers <= 0) ERROR(("Invalid nlayers=%d", nlayers));
        if (overlap < 0) ERROR(("Invalid overlap=%d", overlap));

        int64_t total = (int64_t)n1 + (int64_t)(nlayers - 1) * (int64_t)overlap;
        int base = (int)(total / nlayers);
        int rem = (int)(total % nlayers);

        if (base <= 0) ERROR(("nlayers=%d too large for n1=%d", nlayers, n1));
        if (overlap >= base) {
            ERROR(("overlap=%d too large for n1=%d and nlayers=%d", overlap, n1, nlayers));
        }

        layers = (layer_range_t*)malloc((size_t)nlayers * sizeof(layer_range_t));
        if (layers == NULL) ERROR(("Out of memory"));

        int start = 0;
        for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
            int m1 = base + ((ilayer < rem) ? 1 : 0);
            int end = start + m1 - 1;

            if (end >= n1) {
                ERROR(("Layer %d out of range: end=%d n1=%d", ilayer + 1, end, n1));
            }

            layers[ilayer].start = start;
            layers[ilayer].end = end;
            layers[ilayer].should_datum_abs = end - overlap + 1;

            start = end - overlap + 1;
        }
    }

    char* label1 = sep_get_hdr(in, "label1", NULL);

    for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
        int start = layers[ilayer].start;
        int end = layers[ilayer].end;
        int m1 = end - start + 1;

        char* out_f = build_output_name(out_f_base, ilayer + 1);

        out = sep_open(out_f, SEP_WRITE, 0);
        sep_copy_headers(out, in);
        out->headers->ndim = in->headers->ndim;

        sep_set_axis(out, 0, m1,
                     in->headers->o[0] + start * in->headers->d[0],
                     in->headers->d[0],
                     label1);
out->headers->o[0]=0;
        if (layers[ilayer].should_datum_abs >= 0) {
            int should_datum = layers[ilayer].should_datum_abs - start;
            sep_set_header_int(out, "should_datum", should_datum);
            sep_set_header_int(out, "should_datum_abs", layers[ilayer].should_datum_abs);
        } else {
            sep_set_header_int(out, "should_datum", -1);
            sep_set_header_int(out, "should_datum_abs", -1);
        }

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
    }

    if (label1) free(label1);
    if (layers) free(layers);

    sep_close(in);
    se_par_destroy();

    return 0;
}