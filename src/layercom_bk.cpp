#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>

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
	char* full_model_f = NULL;
	char* layer_model_base = NULL;
	char* layer_image_base = NULL;
	char* out_f = NULL;
	float scalez = 1.0f;
	float scalex_max = 1.0f;
	int scalex_aper=200;
	

	if (!se_have_par("full_model")) ERROR(("Need full_model=")); else full_model_f = se_get_par_str("full_model");
	if (!se_have_par("layer_model_base")) ERROR(("Need layer_model_base=")); else layer_model_base = se_get_par_str("layer_model_base");
	if (!se_have_par("layer_image_base")) ERROR(("Need layer_image_base=")); else layer_image_base = se_get_par_str("layer_image_base");
	if (!se_have_par("nlayers")) ERROR(("Need nlayers=")); else nlayers = se_get_par_int("nlayers");
	if (!se_have_par("output_file")) ERROR(("Need output_file=")); else out_f = se_get_par_str("output_file");

	float scale_layer[nlayers];

	/* scale is used to amplify the layer image from the second layer onward. */
	if (se_have_par("scale")) scalez = se_get_par_float("scale");
	if (se_have_par("scalex_scale")) scalex_max = se_get_par_float("scalex_scale");
	if (se_have_par("scalex_aper")) scalex_aper = se_get_par_int("scalex_aper");

	// 构建每一层的scale，第一层永远时1，最后一层为给定的scale，其他层线性插值
	// TODO 这里可以换成其他的构建方式，比如指数或者分段线性等
	for (int i = 0; i < nlayers; ++i) {
		if (nlayers == 1) {
			scale_layer[i] = 1.0f;
		} else {
			scale_layer[i] = 1.0f + (scalez - 1.0f) * (float)i / (float)(nlayers - 1);
		}
		INFO(("Layer %d scale factor: %g", i + 1, scale_layer[i]));
	}

	if (nlayers <= 0) {
		ERROR(("Invalid nlayers=%d", nlayers));
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

	// 横向 scale：两侧各按 scalex_aper 个采样点从 scalex_max 线性过渡到 1，
	// 中间保持 1。比如 n2_full=200, scalex_aper=80 时，0-79 由 2 过渡到 1，
	// 80-120 保持 1，121-199 由 1 过渡到 2。
	float scalex_layer[n2_full];
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


	


	if (n1_full <= 0 || n2_full <= 0 || n3_full <= 0) {
		ERROR(("Invalid full model size n1=%d n2=%d n3=%d", n1_full, n2_full, n3_full));
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

		if (!sep_have_hdr_int(model, "should_datum")) {
			ERROR(("Need should_datum in %s", model_f));
		}
		int should_datum = sep_get_hdr_int(model, "should_datum", 0);
		int has_abs = sep_have_hdr_int(model, "should_datum_abs");
		int should_datum_abs = has_abs ? sep_get_hdr_int(model, "should_datum_abs", 0) : 0;

		int start = 0;
		if (has_abs) {
			start = should_datum_abs - should_datum;
		} else if (ilayer == 0) {
			start = 0;
		} else {
			int overlap_infer = m1 - should_datum - 1;
			start = prev_end - overlap_infer + 1;
		}

		int end = start + m1 - 1;
		if (start < 0 || end >= n1_full) {
			ERROR(("Layer %d out of range: start=%d end=%d n1=%d", ilayer + 1, start, end, n1_full));
		}

		int overlap = 0;
		if (ilayer > 0) {
			overlap = prev_end - start + 1;
			if (overlap < 0) {
				ERROR(("Layer %d starts after previous end: prev_end=%d start=%d", ilayer + 1, prev_end, start));
			}
		}

		INFO(("Combining layer %d: start=%d end=%d overlap=%d", ilayer + 1, start, end, overlap));

		float*** layer_data = alloc3float(m1, n2_full, n3_full);
		se_fsio_read_float(image->data->io, layer_data[0][0], (size_t)m1 * n2_full * n3_full);

		for (int i3 = 0; i3 < n3_full; ++i3) {
			for (int i2 = 0; i2 < n2_full; ++i2) {
	
				if (ilayer > 0) {
						//在这里应用横向的 scalex_layer，放在最内层循环以避免重复计算
				float scale = scalex_layer[i2];
				if (scale != 1.0f) {
					for (int i1 = 0; i1 < m1; ++i1) {
						layer_data[i3][i2][i1] *= scale;
					}
				}
					}

				


				for (int i1 = 0; i1 < m1; ++i1) {
					int gz = start + i1;
					float val = layer_data[i3][i2][i1];

					/*
					 * Minimal modification:
					 * before stitching the current lower layer with the previous layer,
					 * amplify the current layer image by scale.
					 * The first layer is kept unchanged.
					 */
					if (ilayer > 0) {
						val *= scale_layer[ilayer];
					}

					if (overlap > 0 && gz <= prev_end) {
						int k = gz - start;
						float w = 0.5f;
						if (overlap > 1) {
							w = (float)k / (float)(overlap - 1);
						}
						out_data[i3][i2][gz] = out_data[i3][i2][gz] * (1.0f - w) + val * w;
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
	se_par_destroy();

	return 0;
}