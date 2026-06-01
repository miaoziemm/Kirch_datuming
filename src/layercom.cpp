#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>
#include <math.h>

typedef struct {
	int start;
	int end;
	int should_datum_abs;   /* -1 means no datum; only used for checking consistency */
} layer_range_t;

static int is_start_token(const char* s)
{
	return s &&
	       (strcmp(s, "START") == 0 ||
	        strcmp(s, "Start") == 0 ||
	        strcmp(s, "start") == 0);
}

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


static const char* get_raw_argv_value(int argc, char* argv[], const char* key)
/*<
 * Read a parameter value directly from argv.
 * This is used for scale= because some SE-style parsers treat commas
 * as separators and se_get_par_str("scale") may return only the
 * first value before the first comma.
 *>*/
{
	size_t key_len = strlen(key);
	for (int i = 1; i < argc; ++i) {
		if (strncmp(argv[i], key, key_len) == 0 && argv[i][key_len] == '=') {
			return argv[i] + key_len + 1;
		}
	}
	return NULL;
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

static float parse_float_token(const char* tok, const char* name, int line_no)
{
	char* endp = NULL;
	double v = strtod(tok, &endp);

	if (endp == tok || *endp != '\0') {
		ERROR(("Invalid %s at line %d: %s", name, line_no, tok));
	}

	if (!isfinite(v)) {
		ERROR(("Invalid %s at line %d: %s", name, line_no, tok));
	}

	return (float)v;
}

static int parse_start_token(const char* tok, int line_no)
{
	if (is_start_token(tok)) {
		return 0;
	}
	return parse_int_token(tok, "start", line_no);
}

static void parse_scale_list(const char* scale_arg, float* scale_layer, int nlayers)
/*<
 * Read layer scale factors from scale=.
 * Examples:
 *   scale=1.0                 -> all layers use 1.0
 *   scale=1.0,1.15,1.25,1.10 -> layer-by-layer scale factors
 *>*/
{
	for (int i = 0; i < nlayers; ++i) {
		scale_layer[i] = 1.0f;
	}

	if (scale_arg == NULL || scale_arg[0] == '\0') {
		return;
	}

	char* buf = (char*)malloc(strlen(scale_arg) + 1);
	if (buf == NULL) {
		ERROR(("Out of memory"));
	}
	strcpy(buf, scale_arg);

	char* p = buf;
	char tok[256];
	int nscale = 0;
	float* tmp = alloc1float(nlayers);

	while (next_token(&p, tok, sizeof(tok))) {
		if (nscale >= nlayers) {
			ERROR(("Too many scale values. nlayers=%d, scale=%s", nlayers, scale_arg));
		}
		tmp[nscale] = parse_float_token(tok, "scale", 0);
		nscale++;
	}

	if (nscale == 0) {
		ERROR(("No valid scale values in scale=%s", scale_arg));
	}

	if (nscale == 1) {
		for (int i = 0; i < nlayers; ++i) {
			scale_layer[i] = tmp[0];
		}
	} else if (nscale == nlayers) {
		for (int i = 0; i < nlayers; ++i) {
			scale_layer[i] = tmp[i];
		}
	} else {
		ERROR(("The number of scale values must be 1 or nlayers. nscale=%d, nlayers=%d, scale=%s",
		       nscale, nlayers, scale_arg));
	}

	free1float(tmp);
	free(buf);
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
	if (layers == NULL) {
		ERROR(("Out of memory"));
	}

	char line[4096];
	int line_no = 0;

	while (fgets(line, sizeof(line), fp) != NULL) {
		line_no++;

		char* p = line;
		char t1[256], t2[256], t3[256];

		if (!next_token(&p, t1, sizeof(t1))) {
			continue;
		}
		if (!next_token(&p, t2, sizeof(t2))) {
			ERROR(("Line %d needs at least start,end", line_no));
		}

		if (nlayer >= cap) {
			cap *= 2;
			layers = (layer_range_t*)realloc(layers, (size_t)cap * sizeof(layer_range_t));
			if (layers == NULL) {
				ERROR(("Out of memory"));
			}
		}

		layers[nlayer].start = parse_start_token(t1, line_no);

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

	for (int i = 0; i < nlayer; ++i) {
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
		ERROR(("First layer must start from 0/start/START, but start=%d", layers[0].start));
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
	int nlayers_from_txt = 0;
	char* full_model_f = NULL;
	char* layer_model_base = NULL;
	char* layer_image_base = NULL;
	char* layer_txt = NULL;
	char* out_f = NULL;
	char* scale_arg = NULL;
	float scalex_max = 1.0f;
	int scalex_aper = 200;

	/* Tail-only stitching parameter.
	 * If overlap=50 and stitch_blend=20:
	 *   k=0..29  : keep the upper-layer image only;
	 *   k=30..49 : weighted blending from upper to lower;
	 *   below overlap: use the lower-layer image.
	 * Set stitch_blend=-1 to blend the whole overlap zone.
	 * Set stitch_blend=0 to keep the upper layer within the whole overlap zone.
	 */
	int stitch_blend = -1;


	if (!se_have_par("full_model")) ERROR(("Need full_model=")); else full_model_f = se_get_par_str("full_model");
	if (!se_have_par("layer_model_base")) ERROR(("Need layer_model_base=")); else layer_model_base = se_get_par_str("layer_model_base");
	if (!se_have_par("layer_image_base")) ERROR(("Need layer_image_base=")); else layer_image_base = se_get_par_str("layer_image_base");
	if (!se_have_par("output_file")) ERROR(("Need output_file=")); else out_f = se_get_par_str("output_file");

	/* Keep layer_txt as the main control of the stitching layers.
	 * nlayers is optional and is used only for consistency checking when layer_txt is given.
	 */
	if (se_have_par("nlayers")) nlayers = se_get_par_int("nlayers");
	if (se_have_par("layer_txt")) layer_txt = se_get_par_str("layer_txt");

	if (!layer_txt && nlayers <= 0) {
		ERROR(("Need layer_txt= or nlayers="));
	}

	/* scale accepts either one value for all layers or nlayers comma-separated values.
	 * Read it directly from argv first, because se_get_par_str("scale") may
	 * keep only the first value when commas are used.
	 */
	scale_arg = (char*)get_raw_argv_value(argc, argv, "scale");
	if (scale_arg == NULL && se_have_par("scale")) scale_arg = se_get_par_str("scale");
	if (se_have_par("scalex_scale")) scalex_max = se_get_par_float("scalex_scale");
	if (se_have_par("scalex_aper")) scalex_aper = se_get_par_int("scalex_aper");

	if (se_have_par("stitch_blend")) stitch_blend = se_get_par_int("stitch_blend");


	if (stitch_blend < -1) {
		ERROR(("Invalid stitch_blend=%d", stitch_blend));
	}
	sep_t* full_model = sep_open(full_model_f, SEP_READ, 0);
	int full_ndim = (int)sep_get_min_ndim(full_model);
	if (full_ndim < 1) full_ndim = 1;
	if (full_ndim > 3) {
		ERROR(("Only support up to 3D full model"));
	}

	int n1_full = full_model->headers->n[0];
	int n2_full = (full_ndim >= 2) ? full_model->headers->n[1] : 1;
	int n3_full = (full_ndim >= 3) ? full_model->headers->n[2] : 1;

	if (n1_full <= 0 || n2_full <= 0 || n3_full <= 0) {
		ERROR(("Invalid full model size n1=%d n2=%d n3=%d", n1_full, n2_full, n3_full));
	}

	layer_range_t* layers = NULL;
	if (layer_txt) {
		layers = read_layer_txt(layer_txt, n1_full, &nlayers_from_txt);
		if (nlayers > 0 && nlayers != nlayers_from_txt) {
			ERROR(("nlayers=%d is inconsistent with layer_txt=%s, which contains %d layers",
			       nlayers, layer_txt, nlayers_from_txt));
		}
		nlayers = nlayers_from_txt;
		INFO(("Using layer_txt=%s. Number of layers = %d", layer_txt, nlayers));
	}

	if (nlayers <= 0) {
		ERROR(("Invalid nlayers=%d", nlayers));
	}

	float* scale_layer = alloc1float(nlayers);
	parse_scale_list(scale_arg, scale_layer, nlayers);
	for (int i = 0; i < nlayers; ++i) {
		INFO(("Layer %d scale factor: %g", i + 1, scale_layer[i]));
	}

	INFO(("Tail blending: stitch_blend=%d (-1 means full-overlap blending)", stitch_blend));

	/* 横向 scale：两侧各按 scalex_aper 个采样点从 scalex_max 线性过渡到 1，
	 * 中间保持 1。
	 */
	float* scalex_layer = alloc1float(n2_full);
	for (int i2 = 0; i2 < n2_full; ++i2) {
		if (scalex_aper > 0) {
			int dist_left = i2;
			int dist_right = n2_full - 1 - i2;
			int dist_edge = dist_left < dist_right ? dist_left : dist_right;
			if (dist_edge < scalex_aper) {
				if (scalex_aper == 1) {
					scalex_layer[i2] = (dist_edge == 0) ? scalex_max : 1.0f;
				} else {
					float t = 1.0f - (float)dist_edge / (float)(scalex_aper - 1);
					scalex_layer[i2] = 1.0f + (scalex_max - 1.0f) * t;
				}
			} else {
				scalex_layer[i2] = 1.0f;
			}
		} else {
			scalex_layer[i2] = 1.0f;
		}
	}

	int64_t total = (int64_t)n1_full * (int64_t)n2_full * (int64_t)n3_full;
	float*** out_data = alloc3float(n1_full, n2_full, n3_full);
	memset(out_data[0][0], 0, (size_t)total * sizeof(float));

	int prev_end = -1;

	for (int ilayer = 0; ilayer < nlayers; ++ilayer) {
		char* model_f = build_output_name(layer_model_base, ilayer + 1);
		char* image_f = build_output_name(layer_image_base, ilayer + 1);

		sep_t* model = sep_open(model_f, SEP_READ, 0);
		sep_t* image = sep_open(image_f, SEP_READ, 0);

		int m1 = model->headers->n[0];
		if (m1 <= 0) {
			ERROR(("Invalid layer model n1=%d for %s", m1, model_f));
		}

		int img_ndim = (int)sep_get_min_ndim(image);
		if (img_ndim < 1) img_ndim = 1;
		if (img_ndim > 3) {
			ERROR(("Only support up to 3D layer image"));
		}

		int img_n1 = image->headers->n[0];
		int img_n2 = (img_ndim >= 2) ? image->headers->n[1] : 1;
		int img_n3 = (img_ndim >= 3) ? image->headers->n[2] : 1;

		if (img_n1 != m1 || img_n2 != n2_full || img_n3 != n3_full) {
			ERROR(("Layer image size mismatch for %s: n1=%d n2=%d n3=%d (expected n1=%d n2=%d n3=%d)",
			       image_f, img_n1, img_n2, img_n3, m1, n2_full, n3_full));
		}

		int start = 0;
		int end = 0;

		if (layers) {
			start = layers[ilayer].start;
			end = layers[ilayer].end;

			int expected_m1 = end - start + 1;
			if (m1 != expected_m1) {
				ERROR(("Layer %d size mismatch with layer_txt: model/image n1=%d, but range [%d,%d] needs n1=%d",
				       ilayer + 1, m1, start, end, expected_m1));
			}
		} else {
			if (!sep_have_hdr_int(model, "should_datum")) {
				ERROR(("Need should_datum in %s", model_f));
			}
			int should_datum = sep_get_hdr_int(model, "should_datum", 0);
			int has_abs = sep_have_hdr_int(model, "should_datum_abs");
			int should_datum_abs = has_abs ? sep_get_hdr_int(model, "should_datum_abs", 0) : 0;

			if (has_abs && should_datum_abs >= 0 && should_datum >= 0) {
				start = should_datum_abs - should_datum;
			} else if (ilayer == 0) {
				start = 0;
			} else if (should_datum >= 0) {
				int overlap_infer = m1 - should_datum - 1;
				start = prev_end - overlap_infer + 1;
			} else {
				ERROR(("Cannot infer start depth for layer %d from %s. Please provide layer_txt=.",
				       ilayer + 1, model_f));
			}

			end = start + m1 - 1;
		}

		if (start < 0 || end >= n1_full) {
			ERROR(("Layer %d out of range: start=%d end=%d n1=%d", ilayer + 1, start, end, n1_full));
		}

		int overlap = 0;
		if (ilayer > 0) {
			overlap = prev_end - start + 1;
			if (overlap < 0) {
				ERROR(("Layer %d starts after previous end: prev_end=%d start=%d",
				       ilayer + 1, prev_end, start));
			}
		}

		int blend_z = overlap;
		int blend_start = 0;
		if (overlap > 0) {
			if (stitch_blend < 0) {
				blend_z = overlap;
			} else {
				blend_z = stitch_blend;
				if (blend_z > overlap) {
					ERROR(("Layer %d: stitch_blend=%d is larger than actual overlap=%d",
					       ilayer + 1, blend_z, overlap));
				}
			}
			blend_start = overlap - blend_z;
		}

		INFO(("Combining layer %d: start=%d end=%d overlap=%d keep_upper=%d blend_tail=%d",
		      ilayer + 1, start, end, overlap, blend_start, blend_z));

		float*** layer_data = alloc3float(m1, n2_full, n3_full);
		se_fsio_read_float(image->data->io, layer_data[0][0],
		                   (size_t)m1 * (size_t)n2_full * (size_t)n3_full);

		float scale_z = scale_layer[ilayer];

		for (int i3 = 0; i3 < n3_full; ++i3) {
			for (int i2 = 0; i2 < n2_full; ++i2) {
				float scale_x = (ilayer > 0) ? scalex_layer[i2] : 1.0f;
				float scale = scale_x * scale_z;

				for (int i1 = 0; i1 < m1; ++i1) {
					int gz = start + i1;
					float val = layer_data[i3][i2][i1] * scale;

					if (overlap > 0 && gz <= prev_end) {
						int k = gz - start;  /* 0 <= k <= overlap-1 */

						if (k < blend_start || blend_z == 0) {
							/* Keep the upper-layer result in the early part of the overlap zone. */
							continue;
						} else {
							int kb = k - blend_start;
							float w = 1.0f;
							if (blend_z > 1) {
								w = (float)kb / (float)(blend_z - 1);
							}
							out_data[i3][i2][gz] =
							    out_data[i3][i2][gz] * (1.0f - w) + val * w;
						}
					} else {
						out_data[i3][i2][gz] = val;
					}
				}
			}
		}

		prev_end = end;

		free3float(layer_data);
		sep_close(model);
		sep_close(image);
		free(model_f);
		free(image_f);
	}

	char* first_image_f = build_output_name(layer_image_base, 1);
	sep_t* first_image = sep_open(first_image_f, SEP_READ, 0);

	sep_t* out = sep_open(out_f, SEP_WRITE, 0);
	sep_copy_headers(out, first_image);
	out->headers->ndim = first_image->headers->ndim;

	char* label1 = sep_get_hdr(first_image, "label1", NULL);
	sep_set_axis(out, 0, n1_full, full_model->headers->o[0], full_model->headers->d[0], label1);
	if (label1) free(label1);

	sep_write_headers(out);
	se_fsio_write_float(out->data->io, out_data[0][0], (size_t)total);

	sep_close(out);
	sep_close(first_image);
	free(first_image_f);
	sep_close(full_model);

	free3float(out_data);
	free1float(scale_layer);
	free1float(scalex_layer);
	if (layers) free(layers);

	se_par_destroy();

	return 0;
}
