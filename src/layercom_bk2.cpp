#include <SEBASIC/include/se_basic.h>
#include <SEFILESYSTEM/include/se_fs.h>

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdint.h>
#include <ctype.h>

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

static int imax2(int a, int b)
{
	return a > b ? a : b;
}

static int imin2(int a, int b)
{
	return a < b ? a : b;
}

static double dsqrt_safe(double x)
/*< Compute sqrt(x) without requiring libm. >*/
{
	if (x <= 0.0) return 0.0;
	double y = (x >= 1.0) ? x : 1.0;
	for (int i = 0; i < 30; ++i) {
		y = 0.5 * (y + x / y);
	}
	return y;
}

static float clamp_float(float v, float vmin, float vmax)
{
	if (v < vmin) return vmin;
	if (v > vmax) return vmax;
	return v;
}

static float compute_auto_layer_scale(float*** out_data,
							 float*** layer_data,
							 int m1, int n2, int n3,
							 int start, int prev_end,
							 int overlap, int drop_z, int blend_z,
							 int auto_scale_nz,
							 float auto_scale_min,
							 float auto_scale_max,
							 float auto_scale_eps)
/*<
 * Estimate the vertical amplitude scale of the current lower-layer image.
 *
 * If the first effective samples of the lower layer still lie in the overlap
 * zone, the scale is estimated from the same physical z positions:
 *
 *      scale = sqrt( energy(upper result in overlap) /
 *                    energy(lower result in overlap) ).
 *
 * The lower effective part starts after stitch_drop.  Therefore, for the
 * common case overlap=50, stitch_drop=30, stitch_blend=20, the scale is
 * estimated from k=30..49, i.e., the true effective blending part of the
 * lower layer.
 *
 * If stitch_drop removes the whole overlap, no common z positions remain.
 * In this fallback case, the function compares the last nz samples of the
 * existing upper result with the first nz samples of the lower effective
 * result.  This corresponds to matching the amplitude energy across the
 * stitching boundary.
 *>*/
{
	if (overlap <= 0 || prev_end < 0) return 1.0f;
	if (auto_scale_min <= 0.0f) auto_scale_min = 0.1f;
	if (auto_scale_max < auto_scale_min) auto_scale_max = auto_scale_min;
	if (auto_scale_eps <= 0.0f) auto_scale_eps = 1.0e-20f;

	double upper_energy = 0.0;
	double lower_energy = 0.0;
	int nz_use = 0;
	int same_z_match = 0;

	if (drop_z < overlap) {
		/* Preferred case: compare upper and lower images at the same overlapped z positions. */
		int available = overlap - drop_z;
		if (blend_z > 0 && blend_z < available) available = blend_z;
		if (auto_scale_nz > 0 && auto_scale_nz < available) available = auto_scale_nz;

		if (available <= 0) return 1.0f;
		nz_use = available;
		same_z_match = 1;

		for (int iz = 0; iz < nz_use; ++iz) {
			int k = drop_z + iz;
			int gz = start + k;
			if (gz < 0 || gz > prev_end || k < 0 || k >= m1) continue;
			for (int i3 = 0; i3 < n3; ++i3) {
				for (int i2 = 0; i2 < n2; ++i2) {
					double u = (double)out_data[i3][i2][gz];
					double l = (double)layer_data[i3][i2][k];
					upper_energy += u * u;
					lower_energy += l * l;
				}
			}
		}
	} else {
		/* Fallback case: the whole overlap is discarded, so match across the boundary. */
		int lower_start = overlap;
		int upper_available = prev_end + 1;
		int lower_available = m1 - lower_start;
		int available = upper_available < lower_available ? upper_available : lower_available;
		if (auto_scale_nz > 0 && auto_scale_nz < available) available = auto_scale_nz;
		if (auto_scale_nz <= 0 && available > 1) available = 1;

		if (available <= 0) return 1.0f;
		nz_use = available;
		same_z_match = 0;

		for (int iz = 0; iz < nz_use; ++iz) {
			int gz_upper = prev_end - nz_use + 1 + iz;
			int k_lower = lower_start + iz;
			if (gz_upper < 0 || k_lower < 0 || k_lower >= m1) continue;
			for (int i3 = 0; i3 < n3; ++i3) {
				for (int i2 = 0; i2 < n2; ++i2) {
					double u = (double)out_data[i3][i2][gz_upper];
					double l = (double)layer_data[i3][i2][k_lower];
					upper_energy += u * u;
					lower_energy += l * l;
				}
			}
		}
	}

	if (lower_energy <= (double)auto_scale_eps) {
		INFO(("Auto scale skipped: lower-layer energy is too small. Use scale=1.0"));
		return 1.0f;
	}

	double ratio = (upper_energy + (double)auto_scale_eps) /
			   (lower_energy + (double)auto_scale_eps);
	float scale = (float)dsqrt_safe(ratio);
	if (!(scale > 0.0f) || scale != scale) {
		INFO(("Auto scale skipped: invalid scale value. Use scale=1.0"));
		return 1.0f;
	}

	float unclamped = scale;
	scale = clamp_float(scale, auto_scale_min, auto_scale_max);

	INFO(("Auto scale: same_z_match=%d nz=%d upper_energy=%g lower_energy=%g raw_scale=%g clamped_scale=%g",
		  same_z_match, nz_use, upper_energy, lower_energy, unclamped, scale));

	return scale;
}

static void apply_lateral_scale(float*** data, int n1, int n2, int n3, const float* scalex_layer)
/*< Apply the x-dependent lateral amplitude compensation to a layer image. >*/
{
	for (int i3 = 0; i3 < n3; ++i3) {
		for (int i2 = 0; i2 < n2; ++i2) {
			float scale = scalex_layer[i2];
			if (scale != 1.0f) {
				for (int i1 = 0; i1 < n1; ++i1) {
					data[i3][i2][i1] *= scale;
				}
			}
		}
	}
}

static void parse_scale_list(const char* scale_list, float* scale_layer, int nlayers)
/*<
 * Parse a comma-separated per-layer scale list.
 * Example:
 *     scale_list=1.0,0.95,1.0,0.9
 * The number of values must be exactly equal to nlayers.
 *>*/
{
	if (scale_list == NULL || scale_list[0] == '\0') {
		ERROR(("Empty scale_list"));
	}

	char* work = (char*)malloc(strlen(scale_list) + 1);
	if (work == NULL) {
		ERROR(("Out of memory in parse_scale_list"));
	}
	strcpy(work, scale_list);

	int count = 0;
	char* token = strtok(work, ",");

	while (token != NULL) {
		while (*token != '\0' && isspace((unsigned char)*token)) token++;

		char* endptr = NULL;
		float value = strtof(token, &endptr);

		if (endptr == token) {
			ERROR(("Invalid scale_list value near '%s'", token));
		}

		while (*endptr != '\0' && isspace((unsigned char)*endptr)) endptr++;
		if (*endptr != '\0') {
			ERROR(("Invalid scale_list value near '%s'", token));
		}

		if (count >= nlayers) {
			ERROR(("Too many values in scale_list. Expected nlayers=%d", nlayers));
		}

		if (value <= 0.0f) {
			ERROR(("Invalid scale_list value %g at layer %d. Scale must be positive.",
				   value, count + 1));
		}

		scale_layer[count] = value;
		count++;

		token = strtok(NULL, ",");
	}

	if (count != nlayers) {
		ERROR(("Invalid scale_list: got %d values, but nlayers=%d", count, nlayers));
	}

	free(work);
}

static void remove_low_wavenumber(float*** data,
					  int n1, int n2, int n3,
					  int z_aper, int x_aper,
					  float strength)
/*<
 * Remove the low-wavenumber component from each layer image.
 * The implementation uses a separable spatial high-pass approximation:
 *      output = input - strength * smooth_zx(input).
 * z_aper and x_aper are half-window lengths in samples.
 *>*/
{
	if (strength == 0.0f || (z_aper <= 0 && x_aper <= 0)) return;
	if (z_aper < 0) z_aper = 0;
	if (x_aper < 0) x_aper = 0;

	int nmax = imax2(n1, n2);
	float* tmp = (float*)malloc((size_t)n1 * (size_t)n2 * sizeof(float));
	float* smooth = (float*)malloc((size_t)n1 * (size_t)n2 * sizeof(float));
	float* prefix = (float*)malloc(((size_t)nmax + 1) * sizeof(float));

	if (tmp == NULL || smooth == NULL || prefix == NULL) {
		ERROR(("Out of memory in remove_low_wavenumber"));
	}

	for (int i3 = 0; i3 < n3; ++i3) {
		/* 1) smoothing along z/depth direction. */
		if (z_aper > 0) {
			for (int i2 = 0; i2 < n2; ++i2) {
				prefix[0] = 0.0f;
				for (int i1 = 0; i1 < n1; ++i1) {
					prefix[i1 + 1] = prefix[i1] + data[i3][i2][i1];
				}
				for (int i1 = 0; i1 < n1; ++i1) {
					int a = imax2(0, i1 - z_aper);
					int b = imin2(n1 - 1, i1 + z_aper);
					int cnt = b - a + 1;
					tmp[(size_t)i2 * (size_t)n1 + (size_t)i1] = (prefix[b + 1] - prefix[a]) / (float)cnt;
				}
			}
		} else {
			for (int i2 = 0; i2 < n2; ++i2) {
				memcpy(tmp + (size_t)i2 * (size_t)n1, data[i3][i2], (size_t)n1 * sizeof(float));
			}
		}

		/* 2) smoothing along x/lateral direction. */
		if (x_aper > 0) {
			for (int i1 = 0; i1 < n1; ++i1) {
				prefix[0] = 0.0f;
				for (int i2 = 0; i2 < n2; ++i2) {
					prefix[i2 + 1] = prefix[i2] + tmp[(size_t)i2 * (size_t)n1 + (size_t)i1];
				}
				for (int i2 = 0; i2 < n2; ++i2) {
					int a = imax2(0, i2 - x_aper);
					int b = imin2(n2 - 1, i2 + x_aper);
					int cnt = b - a + 1;
					smooth[(size_t)i2 * (size_t)n1 + (size_t)i1] = (prefix[b + 1] - prefix[a]) / (float)cnt;
				}
			}
		} else {
			memcpy(smooth, tmp, (size_t)n1 * (size_t)n2 * sizeof(float));
		}

		/* 3) high-pass output: remove the smooth low-wavenumber component. */
		for (int i2 = 0; i2 < n2; ++i2) {
			for (int i1 = 0; i1 < n1; ++i1) {
				data[i3][i2][i1] -= strength * smooth[(size_t)i2 * (size_t)n1 + (size_t)i1];
			}
		}
	}

	free(prefix);
	free(smooth);
	free(tmp);
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
	char* scale_list = NULL;

	/* Automatic vertical scale estimation. */
	int auto_scale = 1;
	int auto_scale_nz = -1;
	float auto_scale_min = 0.5f;
	float auto_scale_max = 2.0f;
	float auto_scale_eps = 1.0e-20f;

	float scalex_max = 1.0f;
	int scalex_aper = 200;

	/* Low-wavenumber filtering parameters. */
	int lowk_filter = 1;
	int lowk_z_aper = 20;
	int lowk_x_aper = 60;
	float lowk_strength = 1.0f;

	/* Stitching control in the overlap zone.
	 * For example: overlap=50, stitch_drop=30, stitch_blend=20 means:
	 *   k=0..29  : discard the current lower-layer image and keep upper-layer result;
	 *   k=30..49 : weighted blending from upper to lower;
	 *   below the overlap zone: use the current lower-layer image.
	 */
	int stitch_drop = 0;
	int stitch_blend = -1;

	if (!se_have_par("full_model")) ERROR(("Need full_model=")); else full_model_f = se_get_par_str("full_model");
	if (!se_have_par("layer_model_base")) ERROR(("Need layer_model_base=")); else layer_model_base = se_get_par_str("layer_model_base");
	if (!se_have_par("layer_image_base")) ERROR(("Need layer_image_base=")); else layer_image_base = se_get_par_str("layer_image_base");
	if (!se_have_par("nlayers")) ERROR(("Need nlayers=")); else nlayers = se_get_par_int("nlayers");
	if (!se_have_par("output_file")) ERROR(("Need output_file=")); else out_f = se_get_par_str("output_file");

	if (nlayers <= 0) {
		ERROR(("Invalid nlayers=%d", nlayers));
	}

	/* Vertical amplitude scaling.
	 * If scale_list is provided, it overrides the old linearly increasing scale.
	 * Example: scale_list=1.0,0.95,1.0,0.9 for nlayers=4.
	 */
	if (se_have_par("scale")) scalez = se_get_par_float("scale");
	if (se_have_par("scale_list")) scale_list = se_get_par_str("scale_list");

	if (se_have_par("auto_scale")) auto_scale = se_get_par_int("auto_scale");
	if (se_have_par("auto_scale_nz")) auto_scale_nz = se_get_par_int("auto_scale_nz");
	if (se_have_par("auto_scale_min")) auto_scale_min = se_get_par_float("auto_scale_min");
	if (se_have_par("auto_scale_max")) auto_scale_max = se_get_par_float("auto_scale_max");
	if (se_have_par("auto_scale_eps")) auto_scale_eps = se_get_par_float("auto_scale_eps");

	if (se_have_par("scalex_scale")) scalex_max = se_get_par_float("scalex_scale");
	if (se_have_par("scalex_aper")) scalex_aper = se_get_par_int("scalex_aper");

	if (se_have_par("lowk_filter")) lowk_filter = se_get_par_int("lowk_filter");
	if (se_have_par("lowk_z_aper")) lowk_z_aper = se_get_par_int("lowk_z_aper");
	if (se_have_par("lowk_x_aper")) lowk_x_aper = se_get_par_int("lowk_x_aper");
	if (se_have_par("lowk_strength")) lowk_strength = se_get_par_float("lowk_strength");

	if (se_have_par("stitch_drop")) stitch_drop = se_get_par_int("stitch_drop");
	if (se_have_par("stitch_blend")) stitch_blend = se_get_par_int("stitch_blend");

	if (lowk_z_aper < 0 || lowk_x_aper < 0) {
		ERROR(("Invalid lowk_z_aper=%d or lowk_x_aper=%d", lowk_z_aper, lowk_x_aper));
	}
	if (lowk_strength < 0.0f) {
		ERROR(("Invalid lowk_strength=%g", lowk_strength));
	}
	if (stitch_drop < 0 || stitch_blend < -1) {
		ERROR(("Invalid stitch_drop=%d or stitch_blend=%d", stitch_drop, stitch_blend));
	}
	if (auto_scale_nz < -1) {
		ERROR(("Invalid auto_scale_nz=%d. Use -1 for automatic choice, or a positive integer.", auto_scale_nz));
	}
	if (auto_scale_min <= 0.0f || auto_scale_max <= 0.0f || auto_scale_max < auto_scale_min) {
		ERROR(("Invalid auto_scale_min=%g or auto_scale_max=%g", auto_scale_min, auto_scale_max));
	}
	if (auto_scale_eps <= 0.0f) {
		ERROR(("Invalid auto_scale_eps=%g", auto_scale_eps));
	}

	float* scale_layer = (float*)malloc((size_t)nlayers * sizeof(float));
	if (scale_layer == NULL) {
		ERROR(("Out of memory"));
	}

	/* Construct fallback vertical scale factors.
	 * When auto_scale=1, the factors from the second layer onward are
	 * overwritten by the energy-matching estimates during stitching.
	 * When auto_scale=0, scale_list has priority; otherwise the old linear
	 * scale mode is kept for backward compatibility.
	 */
	if (scale_list != NULL) {
		parse_scale_list(scale_list, scale_layer, nlayers);
		INFO(("Using manual per-layer scale_list=%s", scale_list));
	} else {
		for (int i = 0; i < nlayers; ++i) {
			if (nlayers == 1) {
				scale_layer[i] = 1.0f;
			} else {
				scale_layer[i] = 1.0f + (scalez - 1.0f) * (float)i / (float)(nlayers - 1);
			}
		}
		INFO(("No scale_list provided. Use old linear scale mode with scale=%g", scalez));
	}
	for (int i = 0; i < nlayers; ++i) {
		INFO(("Layer %d scale factor: %g", i + 1, scale_layer[i]));
	}

	INFO(("Low-k filter: enable=%d z_aper=%d x_aper=%d strength=%g",
		  lowk_filter, lowk_z_aper, lowk_x_aper, lowk_strength));
	INFO(("Stitch control: stitch_drop=%d stitch_blend=%d (-1 means blending all remaining overlap)",
		  stitch_drop, stitch_blend));
	INFO(("Auto scale: enable=%d nz=%d min=%g max=%g eps=%g",
		  auto_scale, auto_scale_nz, auto_scale_min, auto_scale_max, auto_scale_eps));

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

	// 横向 scale：两侧各按 scalex_aper 个采样点从 scalex_max 线性过渡到 1，
	// 中间保持 1。比如 n2_full=200, scalex_aper=80 时，0-79 由 2 过渡到 1，
	// 80-120 保持 1，121-199 由 1 过渡到 2。
	float* scalex_layer = (float*)malloc((size_t)n2_full * sizeof(float));
	if (scalex_layer == NULL) {
		ERROR(("Out of memory"));
	}

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

		int drop_z = 0;
		int blend_z = overlap;
		if (overlap > 0) {
			drop_z = stitch_drop;
			blend_z = (stitch_blend >= 0) ? stitch_blend : (overlap - drop_z);

			if (drop_z > overlap) {
				ERROR(("Layer %d: stitch_drop=%d is larger than actual overlap=%d",
					   ilayer + 1, drop_z, overlap));
			}
			if (drop_z + blend_z > overlap) {
				ERROR(("Layer %d: stitch_drop + stitch_blend = %d is larger than actual overlap=%d",
					   ilayer + 1, drop_z + blend_z, overlap));
			}
		}

		INFO(("Combining layer %d: start=%d end=%d overlap=%d drop=%d blend=%d direct_lower=%d",
			  ilayer + 1, start, end, overlap, drop_z, blend_z, overlap - drop_z - blend_z));

		float*** layer_data = alloc3float(m1, n2_full, n3_full);
		se_fsio_read_float(image->data->io, layer_data[0][0], (size_t)m1 * (size_t)n2_full * (size_t)n3_full);

		if (lowk_filter) {
			remove_low_wavenumber(layer_data, m1, n2_full, n3_full,
						  lowk_z_aper, lowk_x_aper, lowk_strength);
		}

		/* Apply lateral compensation before estimating the automatic vertical scale. */
		if (ilayer > 0) {
			apply_lateral_scale(layer_data, m1, n2_full, n3_full, scalex_layer);
		}

		if (auto_scale && ilayer > 0) {
			scale_layer[ilayer] = compute_auto_layer_scale(out_data, layer_data,
										 m1, n2_full, n3_full,
										 start, prev_end,
										 overlap, drop_z, blend_z,
										 auto_scale_nz,
										 auto_scale_min, auto_scale_max,
										 auto_scale_eps);
			INFO(("Layer %d final auto scale factor: %g", ilayer + 1, scale_layer[ilayer]));
		}

		for (int i3 = 0; i3 < n3_full; ++i3) {
			for (int i2 = 0; i2 < n2_full; ++i2) {

				for (int i1 = 0; i1 < m1; ++i1) {
					int gz = start + i1;
					float val = layer_data[i3][i2][i1];

					/* Apply the manually specified or automatically constructed scale for this layer. */
					val *= scale_layer[ilayer];

					if (overlap > 0 && gz <= prev_end) {
						int k = gz - start;

						if (k < drop_z) {
							/* Discard the top part of the current lower layer in the overlap zone. */
							continue;
						} else if (k < drop_z + blend_z) {
							int kb = k - drop_z;
							float w = 0.5f;
							if (blend_z > 1) {
								w = (float)kb / (float)(blend_z - 1);
							}
							out_data[i3][i2][gz] = out_data[i3][i2][gz] * (1.0f - w) + val * w;
						} else {
							/* If drop + blend is smaller than the actual overlap, use the current lower layer. */
							out_data[i3][i2][gz] = val;
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
	free(scalex_layer);
	free(scale_layer);
	se_par_destroy();

	return 0;
}
