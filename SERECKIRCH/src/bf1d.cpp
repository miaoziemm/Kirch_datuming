


/*
 * butterfly_apply_1d.c
 *
 * 1D butterfly matrix-vector multiplication for kernels of the form
 *
 *     K(r,s) = Amp(r,s) * exp(i * alpha * tau(r,s)),
 *
 * where, for the phase-amplitude-separated extrapolation used here,
 *
 *     alpha = -omega,
 *     K(r,s) = Amp(r,s) * exp(-i * omega * tau(r,s)).
 *
 * Input convention:
 *   tau_mat[row][col]  : row = output/target index r, col = input/source index s
 *   Amp_mat[row*n+col] : complex amplitude matrix, same convention
 *   Uin_is[col]        : input vector
 *   Uout_is[row]       : output vector
 *
 * Main entry:
 *
 *   butterfly_apply_1d_phase_amp(ns,
 *                                tau_mat,
 *                                Amp_mat,
 *                                Uin_is,
 *                                Uout_is,
 *                                omega,
 *                                p,
 *                                leaf_n);
 *
 * The public function is safe for arbitrary ns: it computes the largest
 * legal leading block by butterfly and computes all tail interactions by
 * exact direct summation.  The internal butterfly core is used only when
 * ns = leaf_n * 2^L.
 */

#include "../include/bf1d.h"
#include <string.h>

/*
 * butterfly_apply_1d.c
 *
 * 1D butterfly matrix-vector multiplication for kernels of the form
 *
 *     K(r,s) = Amp(r,s) * exp(i * alpha * tau(r,s)),
 *
 * where, for the phase-amplitude-separated extrapolation used here,
 *
 *     alpha = -omega,
 *     K(r,s) = Amp(r,s) * exp(-i * omega * tau(r,s)).
 *
 * Input convention:
 *   tau_mat[row][col]  : row = output/target index r, col = input/source index s
 *   Amp_mat[row*n+col] : complex amplitude matrix, same convention
 *   Uin_is[col]        : input vector
 *   Uout_is[row]       : output vector
 *
 * Main entry:
 *
 *   butterfly_apply_1d_phase_amp(ns,
 *                                tau_mat,
 *                                Amp_mat,
 *                                Uin_is,
 *                                Uout_is,
 *                                omega,
 *                                p,
 *                                leaf_n);
 *
 * The public function is safe for arbitrary ns: it computes the largest
 * legal leading block by butterfly and computes all tail interactions by
 * exact direct summation.  The internal butterfly core is used only when
 * ns = leaf_n * 2^L.
 */


#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

/*
 * High-frequency safeguard.
 *
 * For frequencies higher than BF_DIRECT_FREQ_CUTOFF_HZ, this file
 * bypasses the butterfly approximation and uses exact direct summation.
 *
 * Only this macro is used. Its unit is Hz.
 *
 * Examples:
 *
 *     #define BF_DIRECT_FREQ_CUTOFF_HZ 60.0f
 *
 * or compile with:
 *
 *     -DBF_DIRECT_FREQ_CUTOFF_HZ=60.0f
 *
 * Set BF_DIRECT_FREQ_CUTOFF_HZ <= 0.0f to disable the high-frequency
 * direct-summation fallback.
 */
#ifndef BF_DIRECT_FREQ_CUTOFF_HZ
#define BF_DIRECT_FREQ_CUTOFF_HZ 60.0f
#endif

int bf1d_use_direct_for_omega(float omega)
{
    if ((float)BF_DIRECT_FREQ_CUTOFF_HZ <= 0.0f) {
        return 0;
    }

    float omega_cutoff =
        2.0f * (float)M_PI * (float)BF_DIRECT_FREQ_CUTOFF_HZ;

    return fabsf(omega) >= omega_cutoff;
}

/* Forward declarations because the butterfly routines can fall back to
 * direct summation before the direct routines are defined below.
 */
static void direct_add_1d_phase_amp_range(int n,
                                          int row0,
                                          int row1,
                                          int col0,
                                          int col1,
                                          float **tau_mat,
                                          const fftwf_complex *Amp_mat,
                                          const fftwf_complex *Uin_is,
                                          fftwf_complex *Uout_is,
                                          float omega);

void direct_apply_1d_phase_amp(int n,
                               float **tau_mat,
                               const fftwf_complex *Amp_mat,
                               const fftwf_complex *Uin_is,
                               fftwf_complex *Uout_is,
                               float omega);

void butterfly_apply_1d_phase_amp_main_tail(int n,
                                            float **tau_mat,
                                            const fftwf_complex *Amp_mat,
                                            const fftwf_complex *Uin_is,
                                            fftwf_complex *Uout_is,
                                            float omega,
                                            int p,
                                            int leaf_n);

static void butterfly_apply_1d_phase_amp_core(int n_org,
                                              float **tau_mat,
                                              const fftwf_complex *Amp_mat,
                                              const fftwf_complex *Uin_is,
                                              fftwf_complex *Uout_is,
                                              float omega,
                                              int p,
                                              int leaf_n);

typedef struct {
    int i1;       /* inclusive, 0-based */
    int i2;       /* inclusive, 0-based */
    double cidx;
    double x0;
    double width;
    double *cheb; /* length p */
} BFBox1D;

typedef struct {
    int L;
    int p;
    int leaf_n;
    int n_pad;
    int *nbox;       /* length L+1 */
    BFBox1D **level; /* level[lev][ibox] */
} BFTree1D;

static void bf_die(const char *msg)
{
    fprintf(stderr, "butterfly_apply_1d error: %s\n", msg);
    exit(EXIT_FAILURE);
}

static size_t sigma_index(int ibox_s, int ibox_r, int it, int nbox_r, int p)
{
    return ((size_t)ibox_s * (size_t)nbox_r + (size_t)ibox_r) * (size_t)p + (size_t)it;
}

static size_t mat_index(int row, int col, int n)
{
    return (size_t)row * (size_t)n + (size_t)col;
}

static void c_mul_exp_i(fftwf_complex out,
                        float ar, float ai,
                        double phase)
{
    float cp = cosf((float)phase);
    float sp = sinf((float)phase);
    out[0] = ar * cp - ai * sp;
    out[1] = ar * sp + ai * cp;
}

static int next_padded_size(int n, int leaf_n)
{
    if (n <= 0 || leaf_n <= 0) {
        bf_die("n and leaf_n must be positive.");
    }

    int nleaf = (n + leaf_n - 1) / leaf_n;
    int pow2 = 1;

    while (pow2 < nleaf) {
        pow2 <<= 1;
    }

    return pow2 * leaf_n;
}

static int int_log2_exact(int x)
{
    int L = 0;
    if (x <= 0) {
        bf_die("invalid value for log2.");
    }

    while ((1 << L) < x) {
        L++;
    }

    if ((1 << L) != x) {
        bf_die("log2 argument is not an integer power of two.");
    }

    return L;
}

static void build_cheb_std(int p, double *cheb_std)
{
    for (int k = 0; k < p; k++) {
        cheb_std[k] = cos(((2.0 * (double)(k + 1) - 1.0) * M_PI)
                          / (2.0 * (double)p));
    }
}

static BFTree1D *build_binary_tree_1d_c(int n_pad, int leaf_n, int p)
{
    int ratio = n_pad / leaf_n;
    int L = int_log2_exact(ratio);

    double *cheb_std = (double *)malloc((size_t)p * sizeof(double));
    if (cheb_std == NULL) {
        bf_die("out of memory for cheb_std.");
    }
    build_cheb_std(p, cheb_std);

    BFTree1D *T = (BFTree1D *)calloc(1, sizeof(BFTree1D));
    if (T == NULL) {
        bf_die("out of memory for tree.");
    }

    T->L = L;
    T->p = p;
    T->leaf_n = leaf_n;
    T->n_pad = n_pad;

    T->nbox = (int *)malloc((size_t)(L + 1) * sizeof(int));
    T->level = (BFBox1D **)malloc((size_t)(L + 1) * sizeof(BFBox1D *));
    if (T->nbox == NULL || T->level == NULL) {
        bf_die("out of memory for tree levels.");
    }

    for (int lev = 0; lev <= L; lev++) {
        int nbox = 1 << lev;
        int box_n = n_pad / nbox;

        T->nbox[lev] = nbox;
        T->level[lev] = (BFBox1D *)calloc((size_t)nbox, sizeof(BFBox1D));
        if (T->level[lev] == NULL) {
            bf_die("out of memory for tree boxes.");
        }

        for (int b = 0; b < nbox; b++) {
            int i1 = b * box_n;
            int i2 = (b + 1) * box_n - 1;

            double xL = (double)i1;
            double xR = (double)i2;
            double x0 = 0.5 * (xL + xR);
            double width = xR - xL;

            BFBox1D *box = &T->level[lev][b];

            box->i1 = i1;
            box->i2 = i2;
            box->cidx = 0.5 * ((double)i1 + (double)i2);
            box->x0 = x0;
            box->width = width;

            box->cheb = (double *)malloc((size_t)p * sizeof(double));
            if (box->cheb == NULL) {
                bf_die("out of memory for Chebyshev nodes.");
            }

            for (int it = 0; it < p; it++) {
                box->cheb[it] = x0 + 0.5 * width * cheb_std[it];
            }
        }
    }

    free(cheb_std);
    return T;
}

static void free_binary_tree_1d_c(BFTree1D *T)
{
    if (T == NULL) return;

    for (int lev = 0; lev <= T->L; lev++) {
        if (T->level != NULL && T->level[lev] != NULL) {
            for (int b = 0; b < T->nbox[lev]; b++) {
                free(T->level[lev][b].cheb);
            }
            free(T->level[lev]);
        }
    }

    free(T->level);
    free(T->nbox);
    free(T);
}

static void bary_weights_1d_c(const double *nodes, int p, double *w)
{
    for (int j = 0; j < p; j++) {
        w[j] = 1.0;
        for (int k = 0; k < p; k++) {
            if (k != j) {
                w[j] /= (nodes[j] - nodes[k]);
            }
        }
    }
}

static void lagrange_basis_all_1d_c(double x,
                                    const double *nodes,
                                    const double *w,
                                    int p,
                                    double *L)
{
    const double tol = 1e-14;
    int hit = -1;

    for (int j = 0; j < p; j++) {
        L[j] = 0.0;
        if (fabs(x - nodes[j]) < tol) {
            hit = j;
        }
    }

    if (hit >= 0) {
        L[hit] = 1.0;
        return;
    }

    double denom = 0.0;
    for (int j = 0; j < p; j++) {
        L[j] = w[j] / (x - nodes[j]);
        denom += L[j];
    }

    for (int j = 0; j < p; j++) {
        L[j] /= denom;
    }
}

static double *build_lagrange_matrix_1d_c(const double *points,
                                          int m,
                                          const double *nodes,
                                          const double *w,
                                          int p)
{
    double *Lmat = (double *)malloc((size_t)m * (size_t)p * sizeof(double));
    double *basis = (double *)malloc((size_t)p * sizeof(double));

    if (Lmat == NULL || basis == NULL) {
        bf_die("out of memory for Lagrange matrix.");
    }

    for (int k = 0; k < m; k++) {
        lagrange_basis_all_1d_c(points[k], nodes, w, p, basis);
        for (int it = 0; it < p; it++) {
            Lmat[(size_t)k * (size_t)p + (size_t)it] = basis[it];
        }
    }

    free(basis);
    return Lmat;
}

static double interp2_real_linear_rowcol(const float *A,
                                         int n,
                                         double source_coord,
                                         double target_coord)
{
    const double tol = 1e-12;

    double x = source_coord;
    double z = target_coord;

    if (x < -tol || x > (double)(n - 1) + tol ||
        z < -tol || z > (double)(n - 1) + tol) {
        bf_die("real interpolation query is outside the padded matrix.");
    }

    if (x < 0.0) x = 0.0;
    if (x > (double)(n - 1)) x = (double)(n - 1);
    if (z < 0.0) z = 0.0;
    if (z > (double)(n - 1)) z = (double)(n - 1);

    double jf = x;
    double ifv = z;

    int j_round = (int)llround(jf);
    int i_round = (int)llround(ifv);

    if (fabs(jf - (double)j_round) < tol &&
        fabs(ifv - (double)i_round) < tol) {
        return (double)A[mat_index(i_round, j_round, n)];
    }

    int j1 = (int)floor(jf);
    int i1 = (int)floor(ifv);

    if (j1 >= n - 1) j1 = n - 2;
    if (i1 >= n - 1) i1 = n - 2;
    if (j1 < 0) j1 = 0;
    if (i1 < 0) i1 = 0;

    int j2 = j1 + 1;
    int i2 = i1 + 1;

    double tx = jf - (double)j1;
    double tz = ifv - (double)i1;

    double Q11 = (double)A[mat_index(i1, j1, n)];
    double Q21 = (double)A[mat_index(i1, j2, n)];
    double Q12 = (double)A[mat_index(i2, j1, n)];
    double Q22 = (double)A[mat_index(i2, j2, n)];

    return (1.0 - tx) * (1.0 - tz) * Q11
         + tx         * (1.0 - tz) * Q21
         + (1.0 - tx) * tz         * Q12
         + tx         * tz         * Q22;
}

static void interp2_complex_linear_rowcol(const fftwf_complex *A,
                                          int n,
                                          double source_coord,
                                          double target_coord,
                                          float *out_r,
                                          float *out_i)
{
    const double tol = 1e-12;

    double x = source_coord;
    double z = target_coord;

    if (x < -tol || x > (double)(n - 1) + tol ||
        z < -tol || z > (double)(n - 1) + tol) {
        bf_die("complex interpolation query is outside the padded matrix.");
    }

    if (x < 0.0) x = 0.0;
    if (x > (double)(n - 1)) x = (double)(n - 1);
    if (z < 0.0) z = 0.0;
    if (z > (double)(n - 1)) z = (double)(n - 1);

    double jf = x;
    double ifv = z;

    int j_round = (int)llround(jf);
    int i_round = (int)llround(ifv);

    if (fabs(jf - (double)j_round) < tol &&
        fabs(ifv - (double)i_round) < tol) {
        size_t idx = mat_index(i_round, j_round, n);
        *out_r = A[idx][0];
        *out_i = A[idx][1];
        return;
    }

    int j1 = (int)floor(jf);
    int i1 = (int)floor(ifv);

    if (j1 >= n - 1) j1 = n - 2;
    if (i1 >= n - 1) i1 = n - 2;
    if (j1 < 0) j1 = 0;
    if (i1 < 0) i1 = 0;

    int j2 = j1 + 1;
    int i2 = i1 + 1;

    double tx = jf - (double)j1;
    double tz = ifv - (double)i1;

    double w11 = (1.0 - tx) * (1.0 - tz);
    double w21 = tx         * (1.0 - tz);
    double w12 = (1.0 - tx) * tz;
    double w22 = tx         * tz;

    size_t idx11 = mat_index(i1, j1, n);
    size_t idx21 = mat_index(i1, j2, n);
    size_t idx12 = mat_index(i2, j1, n);
    size_t idx22 = mat_index(i2, j2, n);

    *out_r = (float)(w11 * A[idx11][0]
                   + w21 * A[idx21][0]
                   + w12 * A[idx12][0]
                   + w22 * A[idx22][0]);

    *out_i = (float)(w11 * A[idx11][1]
                   + w21 * A[idx21][1]
                   + w12 * A[idx12][1]
                   + w22 * A[idx22][1]);
}

static fftwf_complex *complex_alloc_zero(size_t n)
{
    fftwf_complex *x = (fftwf_complex *)fftwf_malloc(n * sizeof(fftwf_complex));
    if (x == NULL) {
        bf_die("out of memory for complex array.");
    }

    for (size_t i = 0; i < n; i++) {
        x[i][0] = 0.0f;
        x[i][1] = 0.0f;
    }

    return x;
}

static float *real_alloc_zero(size_t n)
{
    float *x = (float *)calloc(n, sizeof(float));
    if (x == NULL) {
        bf_die("out of memory for real array.");
    }
    return x;
}

static void copy_and_pad_input(int n_org,
                               int n_pad,
                               float **tau_mat,
                               const fftwf_complex *Amp_mat,
                               const fftwf_complex *Uin_is,
                               float *tau_pad,
                               fftwf_complex *Amp_pad,
                               fftwf_complex *Uin_pad)
{
    for (int r = 0; r < n_org; r++) {
        for (int s = 0; s < n_org; s++) {
            size_t id_pad = mat_index(r, s, n_pad);
            size_t id_org = mat_index(r, s, n_org);

            tau_pad[id_pad] = tau_mat[r][s];

            Amp_pad[id_pad][0] = Amp_mat[id_org][0];
            Amp_pad[id_pad][1] = Amp_mat[id_org][1];
        }
    }

    for (int s = 0; s < n_org; s++) {
        Uin_pad[s][0] = Uin_is[s][0];
        Uin_pad[s][1] = Uin_is[s][1];
    }
}

/*
 * Main butterfly routine.
 *
 * Arguments:
 *   n_org    : original vector length, e.g., ns or nh before padding
 *   tau_mat  : n_org x n_org, tau_mat[row][col]
 *   Amp_mat  : n_org x n_org complex, row-major, Amp_mat[row*n_org + col]
 *   Uin_is   : n_org complex input vector
 *   Uout_is  : n_org complex output vector, overwritten by this function
 *   omega    : physical angular frequency used in exp(-i*omega*tau)
 *   p        : Chebyshev rank, e.g., 60
 *   leaf_n   : leaf size, e.g., 64
 */
static void butterfly_apply_1d_phase_amp_core(int n_org,
                                              float **tau_mat,
                                              const fftwf_complex *Amp_mat,
                                              const fftwf_complex *Uin_is,
                                              fftwf_complex *Uout_is,
                                              float omega,
                                              int p,
                                              int leaf_n)
{
    if (n_org <= 0) {
        bf_die("n_org must be positive.");
    }
    if (p <= 0 || leaf_n <= 0) {
        bf_die("p and leaf_n must be positive.");
    }

    /*
     * High-frequency safeguard:
     * for high frequencies, use exact direct summation instead of the
     * butterfly approximation to avoid loss of accuracy.
     */
    if (bf1d_use_direct_for_omega(omega)) {
        direct_apply_1d_phase_amp(n_org,
                                  tau_mat,
                                  Amp_mat,
                                  Uin_is,
                                  Uout_is,
                                  omega);
        return;
    }

    int n_pad = next_padded_size(n_org, leaf_n);

    /*
     * This core routine is intentionally restricted to legal butterfly
     * sizes.  The previous version padded arbitrary n_org to n_pad and
     * filled the padded tau/Amp/U entries with zeros.  Those artificial
     * boundary values can be interpolated at Chebyshev nodes and may
     * contaminate the unpadded output near the tail.  Arbitrary sizes
     * should therefore be handled by butterfly_apply_1d_phase_amp_main_tail(),
     * which computes the legal main block by butterfly and the tail blocks
     * exactly by direct summation.
     */
    if (n_pad != n_org) {
        bf_die("internal core called with a non-legal butterfly size.");
    }

    int L = int_log2_exact(n_pad / leaf_n);

    /*
     * The MATLAB code uses exp(+i*omega*tau).
     * Here the separated extrapolation kernel is
     *
     *     K = Amp * exp(-i*omega*tau).
     *
     * Therefore, we use alpha = -omega everywhere in the MATLAB recursion.
     */
    double alpha = -(double)omega;

    float *tau_pad = real_alloc_zero((size_t)n_pad * (size_t)n_pad);
    fftwf_complex *Amp_pad =
        complex_alloc_zero((size_t)n_pad * (size_t)n_pad);
    fftwf_complex *Uin_pad =
        complex_alloc_zero((size_t)n_pad);
    fftwf_complex *u_final =
        complex_alloc_zero((size_t)n_pad);

    copy_and_pad_input(n_org,
                       n_pad,
                       tau_mat,
                       Amp_mat,
                       Uin_is,
                       tau_pad,
                       Amp_pad,
                       Uin_pad);

    for (int i = 0; i < n_org; i++) {
        Uout_is[i][0] = 0.0f;
        Uout_is[i][1] = 0.0f;
    }

    BFTree1D *TR = build_binary_tree_1d_c(n_pad, leaf_n, p);
    BFTree1D *TS = build_binary_tree_1d_c(n_pad, leaf_n, p);

    /*
     * init
     *
     * sigma(ibox_s, ibox_r, it) at:
     *   source level = L, target level = 0.
     */
    int nbox_s = TS->nbox[L];
    int nbox_r = TR->nbox[0];

    fftwf_complex *sigma_prev =
        complex_alloc_zero((size_t)nbox_s * (size_t)nbox_r * (size_t)p);

    for (int ibox_s = 0; ibox_s < nbox_s; ibox_s++) {
        BFBox1D *B = &TS->level[L][ibox_s];
        const double *nodes = B->cheb;

        double *bw = (double *)malloc((size_t)p * sizeof(double));
        if (bw == NULL) {
            bf_die("out of memory for barycentric weights.");
        }
        bary_weights_1d_c(nodes, p, bw);

        int Bi1 = B->i1;
        int Bi2 = B->i2;
        int mB = Bi2 - Bi1 + 1;

        double *points = (double *)malloc((size_t)mB * sizeof(double));
        if (points == NULL) {
            bf_die("out of memory for box points.");
        }

        for (int k = 0; k < mB; k++) {
            points[k] = (double)(Bi1 + k);
        }

        double *Lg = build_lagrange_matrix_1d_c(points, mB, nodes, bw, p);

        for (int ibox_r = 0; ibox_r < nbox_r; ibox_r++) {
            BFBox1D *A = &TR->level[0][ibox_r];
            double r0 = A->x0;

            for (int it = 0; it < p; it++) {
                float acc_r = 0.0f;
                float acc_i = 0.0f;

                for (int ps = Bi1; ps <= Bi2; ps++) {
                    int row_in_box = ps - Bi1;

                    double phi1 =
                        interp2_real_linear_rowcol(tau_pad,
                                                   n_pad,
                                                   (double)ps,
                                                   r0);

                    double ph = alpha * phi1;
                    float cp = cosf((float)ph);
                    float sp = sinf((float)ph);

                    float fr = Uin_pad[ps][0];
                    float fi = Uin_pad[ps][1];

                    float er = fr * cp - fi * sp;
                    float ei = fr * sp + fi * cp;

                    float Lval = (float)Lg[(size_t)row_in_box * (size_t)p
                                           + (size_t)it];

                    acc_r += Lval * er;
                    acc_i += Lval * ei;
                }

                double phi2 =
                    interp2_real_linear_rowcol(tau_pad,
                                               n_pad,
                                               nodes[it],
                                               r0);

                /*
                 * sigma *= exp(-i * alpha * phi2)
                 */
                double ph2 = -alpha * phi2;
                float cp2 = cosf((float)ph2);
                float sp2 = sinf((float)ph2);

                size_t idx = sigma_index(ibox_s, ibox_r, it, nbox_r, p);

                sigma_prev[idx][0] = acc_r * cp2 - acc_i * sp2;
                sigma_prev[idx][1] = acc_r * sp2 + acc_i * cp2;
            }
        }

        free(Lg);
        free(points);
        free(bw);
    }

    /*
     * r1: first-half recursion
     */
    int lv_max = L / 2;

    for (int lv = 1; lv <= lv_max; lv++) {
        BFBox1D *boxes_r_now = TR->level[lv];
        BFBox1D *boxes_s_now = TS->level[L - lv];
        BFBox1D *boxes_s_prev = TS->level[L - lv + 1];

        int nbox_r_now = TR->nbox[lv];
        int nbox_s_now = TS->nbox[L - lv];

        int nbox_r_prev = TR->nbox[lv - 1];

        fftwf_complex *sigma_curr =
            complex_alloc_zero((size_t)nbox_s_now
                               * (size_t)nbox_r_now
                               * (size_t)p);

        for (int ibox_r = 0; ibox_r < nbox_r_now; ibox_r++) {
            BFBox1D *A = &boxes_r_now[ibox_r];
            double r0 = A->x0;
            int ibox_r_parent = ibox_r / 2;

            for (int ibox_s = 0; ibox_s < nbox_s_now; ibox_s++) {
                BFBox1D *B = &boxes_s_now[ibox_s];
                const double *nodesB = B->cheb;

                double *wB = (double *)malloc((size_t)p * sizeof(double));
                if (wB == NULL) {
                    bf_die("out of memory for wB.");
                }
                bary_weights_1d_c(nodesB, p, wB);

                fftwf_complex *acc = complex_alloc_zero((size_t)p);

                int child_ids[2] = {2 * ibox_s, 2 * ibox_s + 1};

                for (int child = 0; child < 2; child++) {
                    int ibox_s_child = child_ids[child];
                    BFBox1D *Bc = &boxes_s_prev[ibox_s_child];
                    const double *child_nodes = Bc->cheb;

                    double *Lmat =
                        build_lagrange_matrix_1d_c(child_nodes,
                                                   p,
                                                   nodesB,
                                                   wB,
                                                   p);

                    fftwf_complex *tmp = complex_alloc_zero((size_t)p);

                    for (int jt = 0; jt < p; jt++) {
                        size_t idx_child =
                            sigma_index(ibox_s_child,
                                        ibox_r_parent,
                                        jt,
                                        nbox_r_prev,
                                        p);

                        double tau_val =
                            interp2_real_linear_rowcol(tau_pad,
                                                       n_pad,
                                                       child_nodes[jt],
                                                       r0);

                        double ph = alpha * tau_val;
                        c_mul_exp_i(tmp[jt],
                                    sigma_prev[idx_child][0],
                                    sigma_prev[idx_child][1],
                                    ph);
                    }

                    /*
                     * acc(t) += sum_{t'} Lmat(t',t) * tmp(t')
                     */
                    for (int it = 0; it < p; it++) {
                        for (int jt = 0; jt < p; jt++) {
                            float Lval =
                                (float)Lmat[(size_t)jt * (size_t)p
                                            + (size_t)it];

                            acc[it][0] += Lval * tmp[jt][0];
                            acc[it][1] += Lval * tmp[jt][1];
                        }
                    }

                    fftwf_free(tmp);
                    free(Lmat);
                }

                for (int it = 0; it < p; it++) {
                    double tau_val =
                        interp2_real_linear_rowcol(tau_pad,
                                                   n_pad,
                                                   nodesB[it],
                                                   r0);

                    double ph = -alpha * tau_val;
                    float cp = cosf((float)ph);
                    float sp = sinf((float)ph);

                    size_t idx_curr =
                        sigma_index(ibox_s,
                                    ibox_r,
                                    it,
                                    nbox_r_now,
                                    p);

                    float ar = acc[it][0];
                    float ai = acc[it][1];

                    sigma_curr[idx_curr][0] = ar * cp - ai * sp;
                    sigma_curr[idx_curr][1] = ar * sp + ai * cp;
                }

                fftwf_free(acc);
                free(wB);
            }
        }

        fftwf_free(sigma_prev);
        sigma_prev = sigma_curr;
    }

    /*
     * switch
     */
    int mid_lv = lv_max;

    BFBox1D *boxes_r_mid = TR->level[mid_lv];
    BFBox1D *boxes_s_mid = TS->level[L - mid_lv];

    int nbox_r_mid = TR->nbox[mid_lv];
    int nbox_s_mid = TS->nbox[L - mid_lv];

    fftwf_complex *sigma_sw =
        complex_alloc_zero((size_t)nbox_s_mid
                           * (size_t)nbox_r_mid
                           * (size_t)p);

    for (int ibox_r = 0; ibox_r < nbox_r_mid; ibox_r++) {
        BFBox1D *A = &boxes_r_mid[ibox_r];
        const double *r_nodes = A->cheb;

        for (int ibox_s = 0; ibox_s < nbox_s_mid; ibox_s++) {
            BFBox1D *B = &boxes_s_mid[ibox_s];
            const double *s_nodes = B->cheb;

            for (int it = 0; it < p; it++) {
                float sum_r = 0.0f;
                float sum_i = 0.0f;

                double rt = r_nodes[it];

                for (int jt = 0; jt < p; jt++) {
                    double sj = s_nodes[jt];

                    double tau_val =
                        interp2_real_linear_rowcol(tau_pad,
                                                   n_pad,
                                                   sj,
                                                   rt);

                    float amp_r, amp_i;
                    interp2_complex_linear_rowcol(Amp_pad,
                                                  n_pad,
                                                  sj,
                                                  rt,
                                                  &amp_r,
                                                  &amp_i);

                    double ph = alpha * tau_val;
                    float cp = cosf((float)ph);
                    float sp = sinf((float)ph);

                    /*
                     * K = Amp * exp(i * alpha * tau).
                     */
                    float kr = amp_r * cp - amp_i * sp;
                    float ki = amp_r * sp + amp_i * cp;

                    size_t idx_old =
                        sigma_index(ibox_s,
                                    ibox_r,
                                    jt,
                                    nbox_r_mid,
                                    p);

                    float br = sigma_prev[idx_old][0];
                    float bi = sigma_prev[idx_old][1];

                    sum_r += kr * br - ki * bi;
                    sum_i += kr * bi + ki * br;
                }

                size_t idx_new =
                    sigma_index(ibox_s,
                                ibox_r,
                                it,
                                nbox_r_mid,
                                p);

                sigma_sw[idx_new][0] = sum_r;
                sigma_sw[idx_new][1] = sum_i;
            }
        }
    }

    fftwf_free(sigma_prev);
    sigma_prev = sigma_sw;

    /*
     * r2: second-half recursion
     */
    for (int lv = mid_lv + 1; lv <= L; lv++) {
        BFBox1D *boxes_r_now = TR->level[lv];

        BFBox1D *boxes_r_prev = TR->level[lv - 1];
        BFBox1D *boxes_s_prev = TS->level[L - lv + 1];

        int nbox_r_now = TR->nbox[lv];
        int nbox_s_now = TS->nbox[L - lv];

        int nbox_r_prev = TR->nbox[lv - 1];

        fftwf_complex *sigma_curr =
            complex_alloc_zero((size_t)nbox_s_now
                               * (size_t)nbox_r_now
                               * (size_t)p);

        for (int ibox_r = 0; ibox_r < nbox_r_now; ibox_r++) {
            BFBox1D *A = &boxes_r_now[ibox_r];
            const double *r_nodes = A->cheb;

            int ibox_r_parent = ibox_r / 2;
            BFBox1D *Ap = &boxes_r_prev[ibox_r_parent];
            const double *rp_nodes = Ap->cheb;

            double *wAp = (double *)malloc((size_t)p * sizeof(double));
            if (wAp == NULL) {
                bf_die("out of memory for wAp.");
            }
            bary_weights_1d_c(rp_nodes, p, wAp);

            double *Lmat =
                build_lagrange_matrix_1d_c(r_nodes,
                                           p,
                                           rp_nodes,
                                           wAp,
                                           p);

            for (int ibox_s = 0; ibox_s < nbox_s_now; ibox_s++) {
                fftwf_complex *acc = complex_alloc_zero((size_t)p);
                int child_ids[2] = {2 * ibox_s, 2 * ibox_s + 1};

                for (int child = 0; child < 2; child++) {
                    int ibox_s_child = child_ids[child];
                    BFBox1D *Bc = &boxes_s_prev[ibox_s_child];
                    double s0Bc = Bc->x0;

                    fftwf_complex *tmp = complex_alloc_zero((size_t)p);

                    for (int jt = 0; jt < p; jt++) {
                        size_t idx_child =
                            sigma_index(ibox_s_child,
                                        ibox_r_parent,
                                        jt,
                                        nbox_r_prev,
                                        p);

                        double tau_rp =
                            interp2_real_linear_rowcol(tau_pad,
                                                       n_pad,
                                                       s0Bc,
                                                       rp_nodes[jt]);

                        double ph_rp = -alpha * tau_rp;

                        c_mul_exp_i(tmp[jt],
                                    sigma_prev[idx_child][0],
                                    sigma_prev[idx_child][1],
                                    ph_rp);
                    }

                    /*
                     * y(it) = sum_jt Lmat(it,jt) * tmp(jt)
                     */
                    fftwf_complex *y = complex_alloc_zero((size_t)p);

                    for (int it = 0; it < p; it++) {
                        for (int jt = 0; jt < p; jt++) {
                            float Lval =
                                (float)Lmat[(size_t)it * (size_t)p
                                            + (size_t)jt];

                            y[it][0] += Lval * tmp[jt][0];
                            y[it][1] += Lval * tmp[jt][1];
                        }
                    }

                    for (int it = 0; it < p; it++) {
                        double tau_r =
                            interp2_real_linear_rowcol(tau_pad,
                                                       n_pad,
                                                       s0Bc,
                                                       r_nodes[it]);

                        double ph_r = alpha * tau_r;
                        float cp = cosf((float)ph_r);
                        float sp = sinf((float)ph_r);

                        float yr = y[it][0];
                        float yi = y[it][1];

                        acc[it][0] += yr * cp - yi * sp;
                        acc[it][1] += yr * sp + yi * cp;
                    }

                    fftwf_free(y);
                    fftwf_free(tmp);
                }

                for (int it = 0; it < p; it++) {
                    size_t idx_curr =
                        sigma_index(ibox_s,
                                    ibox_r,
                                    it,
                                    nbox_r_now,
                                    p);

                    sigma_curr[idx_curr][0] = acc[it][0];
                    sigma_curr[idx_curr][1] = acc[it][1];
                }

                fftwf_free(acc);
            }

            free(Lmat);
            free(wAp);
        }

        fftwf_free(sigma_prev);
        sigma_prev = sigma_curr;
    }

    /*
     * termination
     */
    BFBox1D *boxes_r_leaf = TR->level[L];
    BFBox1D *Broot = &TS->level[0][0];
    double s0B = Broot->x0;

    int nbox_r_leaf = TR->nbox[L];

    for (int ibox_r = 0; ibox_r < nbox_r_leaf; ibox_r++) {
        BFBox1D *A = &boxes_r_leaf[ibox_r];

        int Ai1 = A->i1;
        int Ai2 = A->i2;
        int nA = Ai2 - Ai1 + 1;

        double *Arange = (double *)malloc((size_t)nA * sizeof(double));
        if (Arange == NULL) {
            bf_die("out of memory for Arange.");
        }

        for (int k = 0; k < nA; k++) {
            Arange[k] = (double)(Ai1 + k);
        }

        const double *r_nodes = A->cheb;

        double *wA = (double *)malloc((size_t)p * sizeof(double));
        if (wA == NULL) {
            bf_die("out of memory for wA.");
        }
        bary_weights_1d_c(r_nodes, p, wA);

        double *LgA =
            build_lagrange_matrix_1d_c(Arange,
                                       nA,
                                       r_nodes,
                                       wA,
                                       p);

        fftwf_complex *tmp = complex_alloc_zero((size_t)p);

        for (int it = 0; it < p; it++) {
            size_t idx_leaf =
                sigma_index(0,
                            ibox_r,
                            it,
                            nbox_r_leaf,
                            p);

            double tau_rt =
                interp2_real_linear_rowcol(tau_pad,
                                           n_pad,
                                           s0B,
                                           r_nodes[it]);

            double ph_rt = -alpha * tau_rt;

            c_mul_exp_i(tmp[it],
                        sigma_prev[idx_leaf][0],
                        sigma_prev[idx_leaf][1],
                        ph_rt);
        }

        for (int ir = 0; ir < nA; ir++) {
            int global_r = Ai1 + ir;

            float sum_r = 0.0f;
            float sum_i = 0.0f;

            for (int it = 0; it < p; it++) {
                float Lval =
                    (float)LgA[(size_t)ir * (size_t)p + (size_t)it];

                sum_r += Lval * tmp[it][0];
                sum_i += Lval * tmp[it][1];
            }

            double tau_r =
                interp2_real_linear_rowcol(tau_pad,
                                           n_pad,
                                           s0B,
                                           Arange[ir]);

            double ph_r = alpha * tau_r;
            float cp = cosf((float)ph_r);
            float sp = sinf((float)ph_r);

            u_final[global_r][0] = sum_r * cp - sum_i * sp;
            u_final[global_r][1] = sum_r * sp + sum_i * cp;
        }

        fftwf_free(tmp);
        free(LgA);
        free(wA);
        free(Arange);
    }

    /*
     * Copy back only the original, unpadded part.
     */
    for (int i = 0; i < n_org; i++) {
        Uout_is[i][0] = u_final[i][0];
        Uout_is[i][1] = u_final[i][1];
    }

    fftwf_free(sigma_prev);

    free_binary_tree_1d_c(TR);
    free_binary_tree_1d_c(TS);

    fftwf_free(u_final);
    fftwf_free(Uin_pad);
    fftwf_free(Amp_pad);
    free(tau_pad);
}



/*
 * Fast block adaptive-cross-approximation (ACA) version.
 *
 * This version keeps the previous high-accuracy blockwise low-rank idea, but
 * removes the main speed bottlenecks of the first ACA implementation:
 *
 *   1. no malloc/free inside each matrix block;
 *   2. no preliminary nonzero-block scan;
 *   3. residual, pivot column, and pivot row are stored in compact real/imag
 *      work arrays;
 *   4. zero blocks are skipped automatically when the initial residual norm is
 *      zero.
 *
 * The algorithm is still a true low-rank butterfly-style compression; it does
 * not use accuracy-test fallback to direct summation.  Each leaf_n-by-leaf_n
 * block is approximated by a rank-p cross expansion,
 *
 *     K_block ~= sum_k col_k * row_k,
 *
 * and then applied to Uin in factored form.
 */
#ifndef BF_ACA_PIVOT_REL_TOL
#define BF_ACA_PIVOT_REL_TOL 1.0e-7f
#endif

static inline float bf_cabs2(float ar, float ai)
{
    return ar * ar + ai * ai;
}

static inline void bf_cmul(float ar, float ai,
                           float br, float bi,
                           float *cr, float *ci)
{
    *cr = ar * br - ai * bi;
    *ci = ar * bi + ai * br;
}

static inline void bf_cdiv(float ar, float ai,
                           float br, float bi,
                           float *cr, float *ci)
{
    float den = br * br + bi * bi;
    if (den <= 0.0f) {
        *cr = 0.0f;
        *ci = 0.0f;
        return;
    }
    *cr = (ar * br + ai * bi) / den;
    *ci = (ai * br - ar * bi) / den;
}

typedef struct {
    int leaf_cap;
    float *Rr;
    float *Ri;
    float *colr;
    float *coli;
    float *rowr;
    float *rowi;
} BFAcaWork;

static void bf_aca_work_init(BFAcaWork *W, int leaf_n)
{
    memset(W, 0, sizeof(*W));
    W->leaf_cap = leaf_n;
    size_t m2 = (size_t)leaf_n * (size_t)leaf_n;
    W->Rr = (float *)malloc(m2 * sizeof(float));
    W->Ri = (float *)malloc(m2 * sizeof(float));
    W->colr = (float *)malloc((size_t)leaf_n * sizeof(float));
    W->coli = (float *)malloc((size_t)leaf_n * sizeof(float));
    W->rowr = (float *)malloc((size_t)leaf_n * sizeof(float));
    W->rowi = (float *)malloc((size_t)leaf_n * sizeof(float));
    if (W->Rr == NULL || W->Ri == NULL || W->colr == NULL || W->coli == NULL ||
        W->rowr == NULL || W->rowi == NULL) {
        bf_die("out of memory for fast ACA workspace.");
    }
}

static void bf_aca_work_free(BFAcaWork *W)
{
    if (W == NULL) return;
    free(W->Rr);
    free(W->Ri);
    free(W->colr);
    free(W->coli);
    free(W->rowr);
    free(W->rowi);
    memset(W, 0, sizeof(*W));
}

static inline void build_kernel_entry_phase_amp_fast(int n,
                                                     int row,
                                                     int col,
                                                     float **tau_mat,
                                                     const fftwf_complex *Amp_mat,
                                                     float omega,
                                                     float *kr,
                                                     float *ki)
{
    size_t imat = (size_t)row * (size_t)n + (size_t)col;
    float amp_r = Amp_mat[imat][0];
    float amp_i = Amp_mat[imat][1];

    if (fabsf(amp_r) < 1e-20f && fabsf(amp_i) < 1e-20f) {
        *kr = 0.0f;
        *ki = 0.0f;
        return;
    }

    float phase = omega * tau_mat[row][col];
    float cosp = cosf(phase);
    float sinp = sinf(phase);

    /* Amp * exp(-i*phase). */
    *kr = amp_r * cosp + amp_i * sinp;
    *ki = amp_i * cosp - amp_r * sinp;
}

static void aca_add_block_phase_amp_fast(int n,
                                         int row0,
                                         int row1,
                                         int col0,
                                         int col1,
                                         float **tau_mat,
                                         const fftwf_complex *Amp_mat,
                                         const fftwf_complex *Uin_is,
                                         fftwf_complex *Uout_is,
                                         float omega,
                                         int rank,
                                         BFAcaWork *W)
{
    int nr = row1 - row0;
    int nc = col1 - col0;

    if (nr <= 0 || nc <= 0 || rank <= 0) return;
    if (nr > W->leaf_cap || nc > W->leaf_cap) {
        bf_die("ACA block is larger than workspace capacity.");
    }
    if (rank > nr) rank = nr;
    if (rank > nc) rank = nc;

    float *Rr = W->Rr;
    float *Ri = W->Ri;
    float *colr = W->colr;
    float *coli = W->coli;
    float *rowr = W->rowr;
    float *rowi = W->rowi;

    float max0 = 0.0f;
    for (int ir = 0; ir < nr; ir++) {
        int row = row0 + ir;
        for (int jc = 0; jc < nc; jc++) {
            float kr, ki;
            build_kernel_entry_phase_amp_fast(n,
                                              row,
                                              col0 + jc,
                                              tau_mat,
                                              Amp_mat,
                                              omega,
                                              &kr,
                                              &ki);
            size_t id = (size_t)ir * (size_t)nc + (size_t)jc;
            Rr[id] = kr;
            Ri[id] = ki;
            float a2 = bf_cabs2(kr, ki);
            if (a2 > max0) max0 = a2;
        }
    }

    if (max0 <= 0.0f) return;

    float pivot_stop2 = (float)BF_ACA_PIVOT_REL_TOL *
                        (float)BF_ACA_PIVOT_REL_TOL * max0;

    for (int k = 0; k < rank; k++) {
        int piv_i = -1;
        int piv_j = -1;
        float piv_abs2 = 0.0f;

        for (int ir = 0; ir < nr; ir++) {
            size_t base = (size_t)ir * (size_t)nc;
            for (int jc = 0; jc < nc; jc++) {
                float a2 = bf_cabs2(Rr[base + jc], Ri[base + jc]);
                if (a2 > piv_abs2) {
                    piv_abs2 = a2;
                    piv_i = ir;
                    piv_j = jc;
                }
            }
        }

        if (piv_i < 0 || piv_j < 0 || piv_abs2 <= pivot_stop2) {
            break;
        }

        size_t piv_id = (size_t)piv_i * (size_t)nc + (size_t)piv_j;
        float piv_r = Rr[piv_id];
        float piv_iq = Ri[piv_id];

        for (int ir = 0; ir < nr; ir++) {
            size_t id = (size_t)ir * (size_t)nc + (size_t)piv_j;
            colr[ir] = Rr[id];
            coli[ir] = Ri[id];
        }

        size_t piv_row_base = (size_t)piv_i * (size_t)nc;
        for (int jc = 0; jc < nc; jc++) {
            bf_cdiv(Rr[piv_row_base + jc],
                    Ri[piv_row_base + jc],
                    piv_r,
                    piv_iq,
                    &rowr[jc],
                    &rowi[jc]);
        }

        float coeff_r = 0.0f;
        float coeff_i = 0.0f;
        for (int jc = 0; jc < nc; jc++) {
            int col = col0 + jc;
            float tr, ti;
            bf_cmul(rowr[jc], rowi[jc],
                    Uin_is[col][0], Uin_is[col][1],
                    &tr, &ti);
            coeff_r += tr;
            coeff_i += ti;
        }

        for (int ir = 0; ir < nr; ir++) {
            int row = row0 + ir;
            float tr, ti;
            bf_cmul(colr[ir], coli[ir], coeff_r, coeff_i, &tr, &ti);
            Uout_is[row][0] += tr;
            Uout_is[row][1] += ti;
        }

        /* R <- R - col * row. */
        for (int ir = 0; ir < nr; ir++) {
            float cr = colr[ir];
            float ci = coli[ir];
            size_t base = (size_t)ir * (size_t)nc;
            for (int jc = 0; jc < nc; jc++) {
                float pr = cr * rowr[jc] - ci * rowi[jc];
                float pi = cr * rowi[jc] + ci * rowr[jc];
                Rr[base + jc] -= pr;
                Ri[base + jc] -= pi;
            }
        }
    }
}

void butterfly_apply_1d_phase_amp_aca(int n,
                                      float **tau_mat,
                                      const fftwf_complex *Amp_mat,
                                      const fftwf_complex *Uin_is,
                                      fftwf_complex *Uout_is,
                                      float omega,
                                      int p,
                                      int leaf_n)
{
    if (n <= 0) {
        bf_die("n must be positive.");
    }
    if (p <= 0 || leaf_n <= 1) {
        bf_die("p must be positive and leaf_n must be larger than one.");
    }
    if (p >= leaf_n) {
        bf_die("for meaningful interpolation/compression, require p < leaf_n.");
    }

    for (int i = 0; i < n; i++) {
        Uout_is[i][0] = 0.0f;
        Uout_is[i][1] = 0.0f;
    }

    BFAcaWork W;
    bf_aca_work_init(&W, leaf_n);

    for (int row0 = 0; row0 < n; row0 += leaf_n) {
        int row1 = row0 + leaf_n;
        if (row1 > n) row1 = n;

        for (int col0 = 0; col0 < n; col0 += leaf_n) {
            int col1 = col0 + leaf_n;
            if (col1 > n) col1 = n;

            aca_add_block_phase_amp_fast(n,
                                         row0,
                                         row1,
                                         col0,
                                         col1,
                                         tau_mat,
                                         Amp_mat,
                                         Uin_is,
                                         Uout_is,
                                         omega,
                                         p,
                                         &W);
        }
    }

    bf_aca_work_free(&W);
}



typedef struct {
    int row0, row1;
    int col0, col1;
    int nr, nc;
    int rank;
    float *colr; /* rank x nr */
    float *coli;
    float *rowr; /* rank x nc */
    float *rowi;
} BFAcaBlockFactor;

typedef struct BFAcaFactor {
    int n;
    int p;
    int leaf_n;
    int nblocks;
    BFAcaBlockFactor *blocks;
} BFAcaFactor;

static void bf_aca_block_factor_free(BFAcaBlockFactor *B)
{
    if (B == NULL) return;
    free(B->colr);
    free(B->coli);
    free(B->rowr);
    free(B->rowi);
    memset(B, 0, sizeof(*B));
}

void bf1d_aca_factor_destroy(BFAcaFactor *F)
{
    if (F == NULL) return;
    for (int i = 0; i < F->nblocks; i++) {
        bf_aca_block_factor_free(&F->blocks[i]);
    }
    free(F->blocks);
    free(F);
}

static int bf_aca_factorize_one_block(int n,
                                      int row0,
                                      int row1,
                                      int col0,
                                      int col1,
                                      float **tau_mat,
                                      const fftwf_complex *Amp_mat,
                                      float omega,
                                      int rank_max,
                                      BFAcaWork *W,
                                      BFAcaBlockFactor *BF)
{
    int nr = row1 - row0;
    int nc = col1 - col0;
    if (nr <= 0 || nc <= 0 || rank_max <= 0) return 0;
    if (rank_max > nr) rank_max = nr;
    if (rank_max > nc) rank_max = nc;

    float *Rr = W->Rr;
    float *Ri = W->Ri;
    float *colr = W->colr;
    float *coli = W->coli;
    float *rowr = W->rowr;
    float *rowi = W->rowi;

    float max0 = 0.0f;
    for (int ir = 0; ir < nr; ir++) {
        int row = row0 + ir;
        for (int jc = 0; jc < nc; jc++) {
            float kr, ki;
            build_kernel_entry_phase_amp_fast(n, row, col0 + jc,
                                              tau_mat, Amp_mat, omega, &kr, &ki);
            size_t id = (size_t)ir * (size_t)nc + (size_t)jc;
            Rr[id] = kr;
            Ri[id] = ki;
            float a2 = bf_cabs2(kr, ki);
            if (a2 > max0) max0 = a2;
        }
    }
    if (max0 <= 0.0f) return 0;

    BF->row0 = row0;
    BF->row1 = row1;
    BF->col0 = col0;
    BF->col1 = col1;
    BF->nr = nr;
    BF->nc = nc;
    BF->rank = 0;
    BF->colr = (float *)malloc((size_t)rank_max * (size_t)nr * sizeof(float));
    BF->coli = (float *)malloc((size_t)rank_max * (size_t)nr * sizeof(float));
    BF->rowr = (float *)malloc((size_t)rank_max * (size_t)nc * sizeof(float));
    BF->rowi = (float *)malloc((size_t)rank_max * (size_t)nc * sizeof(float));
    if (BF->colr == NULL || BF->coli == NULL || BF->rowr == NULL || BF->rowi == NULL) {
        bf_die("out of memory for ACA block factors.");
    }

    float pivot_stop2 = (float)BF_ACA_PIVOT_REL_TOL *
                        (float)BF_ACA_PIVOT_REL_TOL * max0;

    for (int k = 0; k < rank_max; k++) {
        int piv_i = -1;
        int piv_j = -1;
        float piv_abs2 = 0.0f;

        for (int ir = 0; ir < nr; ir++) {
            size_t base = (size_t)ir * (size_t)nc;
            for (int jc = 0; jc < nc; jc++) {
                float a2 = bf_cabs2(Rr[base + jc], Ri[base + jc]);
                if (a2 > piv_abs2) {
                    piv_abs2 = a2;
                    piv_i = ir;
                    piv_j = jc;
                }
            }
        }

        if (piv_i < 0 || piv_j < 0 || piv_abs2 <= pivot_stop2) break;

        size_t piv_id = (size_t)piv_i * (size_t)nc + (size_t)piv_j;
        float piv_r = Rr[piv_id];
        float piv_iq = Ri[piv_id];

        for (int ir = 0; ir < nr; ir++) {
            size_t id = (size_t)ir * (size_t)nc + (size_t)piv_j;
            colr[ir] = Rr[id];
            coli[ir] = Ri[id];
        }

        size_t piv_row_base = (size_t)piv_i * (size_t)nc;
        for (int jc = 0; jc < nc; jc++) {
            bf_cdiv(Rr[piv_row_base + jc], Ri[piv_row_base + jc],
                    piv_r, piv_iq, &rowr[jc], &rowi[jc]);
        }

        size_t koff_col = (size_t)k * (size_t)nr;
        for (int ir = 0; ir < nr; ir++) {
            BF->colr[koff_col + ir] = colr[ir];
            BF->coli[koff_col + ir] = coli[ir];
        }
        size_t koff_row = (size_t)k * (size_t)nc;
        for (int jc = 0; jc < nc; jc++) {
            BF->rowr[koff_row + jc] = rowr[jc];
            BF->rowi[koff_row + jc] = rowi[jc];
        }
        BF->rank++;

        for (int ir = 0; ir < nr; ir++) {
            float cr = colr[ir];
            float ci = coli[ir];
            size_t base = (size_t)ir * (size_t)nc;
            for (int jc = 0; jc < nc; jc++) {
                float pr = cr * rowr[jc] - ci * rowi[jc];
                float pi = cr * rowi[jc] + ci * rowr[jc];
                Rr[base + jc] -= pr;
                Ri[base + jc] -= pi;
            }
        }
    }

    if (BF->rank <= 0) {
        bf_aca_block_factor_free(BF);
        return 0;
    }
    return 1;
}

BFAcaFactor *bf1d_aca_factor_create_phase_amp(int n,
                                              float **tau_mat,
                                              const fftwf_complex *Amp_mat,
                                              float omega,
                                              int p,
                                              int leaf_n)
{
    if (n <= 0) bf_die("n must be positive.");
    if (p <= 0 || leaf_n <= 1) bf_die("p must be positive and leaf_n must be larger than one.");
    if (p >= leaf_n) bf_die("for meaningful interpolation/compression, require p < leaf_n.");

    int max_blocks_per_dim = (n + leaf_n - 1) / leaf_n;
    int max_blocks = max_blocks_per_dim * max_blocks_per_dim;

    BFAcaFactor *F = (BFAcaFactor *)calloc(1, sizeof(BFAcaFactor));
    if (F == NULL) bf_die("out of memory for ACA factor.");
    F->n = n;
    F->p = p;
    F->leaf_n = leaf_n;
    F->blocks = (BFAcaBlockFactor *)calloc((size_t)max_blocks, sizeof(BFAcaBlockFactor));
    if (F->blocks == NULL) bf_die("out of memory for ACA block list.");

    BFAcaWork W;
    bf_aca_work_init(&W, leaf_n);

    for (int row0 = 0; row0 < n; row0 += leaf_n) {
        int row1 = row0 + leaf_n;
        if (row1 > n) row1 = n;
        for (int col0 = 0; col0 < n; col0 += leaf_n) {
            int col1 = col0 + leaf_n;
            if (col1 > n) col1 = n;
            BFAcaBlockFactor B;
            memset(&B, 0, sizeof(B));
            if (bf_aca_factorize_one_block(n, row0, row1, col0, col1,
                                           tau_mat, Amp_mat, omega, p, &W, &B)) {
                F->blocks[F->nblocks++] = B;
            }
        }
    }

    bf_aca_work_free(&W);
    return F;
}

void bf1d_aca_factor_apply(const BFAcaFactor *F,
                           const fftwf_complex *Uin_is,
                           fftwf_complex *Uout_is)
{
    if (F == NULL) bf_die("null ACA factor.");
    int n = F->n;
    for (int i = 0; i < n; i++) {
        Uout_is[i][0] = 0.0f;
        Uout_is[i][1] = 0.0f;
    }

    for (int ib = 0; ib < F->nblocks; ib++) {
        const BFAcaBlockFactor *B = &F->blocks[ib];
        int nr = B->nr;
        int nc = B->nc;
        int row0 = B->row0;
        int col0 = B->col0;

        for (int k = 0; k < B->rank; k++) {
            size_t roff = (size_t)k * (size_t)nc;
            float coeff_r = 0.0f;
            float coeff_i = 0.0f;
            for (int jc = 0; jc < nc; jc++) {
                int col = col0 + jc;
                float rr = B->rowr[roff + jc];
                float ri = B->rowi[roff + jc];
                float ur = Uin_is[col][0];
                float ui = Uin_is[col][1];
                coeff_r += rr * ur - ri * ui;
                coeff_i += rr * ui + ri * ur;
            }

            size_t coff = (size_t)k * (size_t)nr;
            for (int ir = 0; ir < nr; ir++) {
                int row = row0 + ir;
                float cr = B->colr[coff + ir];
                float ci = B->coli[coff + ir];
                Uout_is[row][0] += cr * coeff_r - ci * coeff_i;
                Uout_is[row][1] += cr * coeff_i + ci * coeff_r;
            }
        }
    }
}

/*
 * Public entry.
 *
 * This wrapper is safe for arbitrary n.  It avoids the loss of accuracy
 * caused by padding tau/Amp/U with zeros: a legal leading block is evaluated
 * by the butterfly core, while all tail interactions are evaluated exactly by
 * direct summation.
 */
void butterfly_apply_1d_phase_amp(int n_org,
                                  float **tau_mat,
                                  const fftwf_complex *Amp_mat,
                                  const fftwf_complex *Uin_is,
                                  fftwf_complex *Uout_is,
                                  float omega,
                                  int p,
                                  int leaf_n)
{
    butterfly_apply_1d_phase_amp_aca(n_org,
                                     tau_mat,
                                     Amp_mat,
                                     Uin_is,
                                     Uout_is,
                                     omega,
                                     p,
                                     leaf_n);
}

/*
 * Direct summation for rows [row0,row1) and cols [col0,col1).
 * The result is ADDED into Uout_is[row].
 *
 * This uses the same kernel as the phase-amplitude-separated extrapolation:
 *
 *     K(row,col) = Amp_mat[row,col] * exp(-i * omega * tau_mat[row][col]).
 */
static void direct_add_1d_phase_amp_range(int n,
                                          int row0,
                                          int row1,
                                          int col0,
                                          int col1,
                                          float **tau_mat,
                                          const fftwf_complex *Amp_mat,
                                          const fftwf_complex *Uin_is,
                                          fftwf_complex *Uout_is,
                                          float omega)
{
    if (row0 < 0) row0 = 0;
    if (col0 < 0) col0 = 0;
    if (row1 > n) row1 = n;
    if (col1 > n) col1 = n;

    for (int row = row0; row < row1; row++) {
        for (int col = col0; col < col1; col++) {

            size_t imat = (size_t)row * (size_t)n + (size_t)col;

            float amp_r = Amp_mat[imat][0];
            float amp_i = Amp_mat[imat][1];

            if (fabsf(amp_r) < 1e-20f &&
                fabsf(amp_i) < 1e-20f) {
                continue;
            }

            float tau = tau_mat[row][col];

            float phase = omega * tau;
            float cosp = cosf(phase);
            float sinp = sinf(phase);

            /*
             * K = Amp * exp(-i*phase).
             * If Amp = amp_r + i amp_i, then
             *
             *     K_r = amp_r*cos + amp_i*sin
             *     K_i = amp_i*cos - amp_r*sin.
             */
            float kr = amp_r * cosp + amp_i * sinp;
            float ki = amp_i * cosp - amp_r * sinp;

            float ar = Uin_is[col][0];
            float ai = Uin_is[col][1];

            Uout_is[row][0] += kr * ar - ki * ai;
            Uout_is[row][1] += kr * ai + ki * ar;
        }
    }
}

/*
 * Full direct summation. This is mainly for verification.
 */
void direct_apply_1d_phase_amp(int n,
                               float **tau_mat,
                               const fftwf_complex *Amp_mat,
                               const fftwf_complex *Uin_is,
                               fftwf_complex *Uout_is,
                               float omega)
{
    for (int i = 0; i < n; i++) {
        Uout_is[i][0] = 0.0f;
        Uout_is[i][1] = 0.0f;
    }

    direct_add_1d_phase_amp_range(n,
                                  0,
                                  n,
                                  0,
                                  n,
                                  tau_mat,
                                  Amp_mat,
                                  Uin_is,
                                  Uout_is,
                                  omega);
}

static int largest_legal_main_size(int n, int leaf_n)
{
    if (n <= 0 || leaf_n <= 0) {
        bf_die("n and leaf_n must be positive.");
    }

    if (n < leaf_n) {
        return 0;
    }

    int n0 = leaf_n;

    while (n0 <= n / 2) {
        n0 *= 2;
    }

    return n0;
}

/*
 * Arbitrary-size version using:
 *
 *     main block butterfly + tail direct summation.
 *
 * Let n0 be the largest valid butterfly size:
 *
 *     n0 = leaf_n * 2^L <= n.
 *
 * Then the matrix-vector product is decomposed as:
 *
 *     [U0]   [K00 K01] [V0]
 *     [U1] = [K10 K11] [V1]
 *
 * where:
 *
 *     K00 is n0 x n0 and is computed by butterfly;
 *     K01, K10, and K11 are computed by direct summation.
 *
 * No padding is used. Therefore, no artificial boundary values are introduced.
 *
 * Input convention:
 *
 *     tau_mat[row][col]
 *     Amp_mat[row*n+col]
 *     Uin_is[col]
 *     Uout_is[row]
 */
void butterfly_apply_1d_phase_amp_main_tail(int n,
                                            float **tau_mat,
                                            const fftwf_complex *Amp_mat,
                                            const fftwf_complex *Uin_is,
                                            fftwf_complex *Uout_is,
                                            float omega,
                                            int p,
                                            int leaf_n)
{
    if (n <= 0) {
        bf_die("n must be positive.");
    }
    if (p <= 0 || leaf_n <= 0) {
        bf_die("p and leaf_n must be positive.");
    }

    /*
     * High-frequency safeguard:
     * at high frequency, compute the full matrix-vector product exactly.
     * This bypasses both the main-block butterfly and the tail decomposition.
     */
    if (bf1d_use_direct_for_omega(omega)) {
        direct_apply_1d_phase_amp(n,
                                  tau_mat,
                                  Amp_mat,
                                  Uin_is,
                                  Uout_is,
                                  omega);
        return;
    }

    for (int i = 0; i < n; i++) {
        Uout_is[i][0] = 0.0f;
        Uout_is[i][1] = 0.0f;
    }

    int n0 = largest_legal_main_size(n, leaf_n);

    /*
     * If the vector is shorter than leaf_n, using butterfly is not meaningful.
     * The exact direct summation is used for the whole vector.
     */
    if (n0 <= 0) {
        direct_add_1d_phase_amp_range(n,
                                      0,
                                      n,
                                      0,
                                      n,
                                      tau_mat,
                                      Amp_mat,
                                      Uin_is,
                                      Uout_is,
                                      omega);
        return;
    }

    /*
     * 1. Main block K00: butterfly on rows 0:n0-1 and cols 0:n0-1.
     *
     * tau_mat can be passed by row pointers because only tau_mat[row][col]
     * with row,col < n0 is accessed. Amp_mat must be copied because the
     * original row stride is n, whereas butterfly_apply_1d_phase_amp expects
     * a row stride of n0.
     */
    float **tau00 = (float **)malloc((size_t)n0 * sizeof(float *));
    fftwf_complex *Amp00 =
        (fftwf_complex *)fftwf_malloc((size_t)n0
                                      * (size_t)n0
                                      * sizeof(fftwf_complex));
    fftwf_complex *Uin0 =
        (fftwf_complex *)fftwf_malloc((size_t)n0 * sizeof(fftwf_complex));
    fftwf_complex *Uout0 =
        (fftwf_complex *)fftwf_malloc((size_t)n0 * sizeof(fftwf_complex));

    if (tau00 == NULL || Amp00 == NULL || Uin0 == NULL || Uout0 == NULL) {
        bf_die("out of memory in main-tail butterfly wrapper.");
    }

    for (int row = 0; row < n0; row++) {
        tau00[row] = tau_mat[row];

        for (int col = 0; col < n0; col++) {
            size_t id0 = (size_t)row * (size_t)n0 + (size_t)col;
            size_t id  = (size_t)row * (size_t)n  + (size_t)col;

            Amp00[id0][0] = Amp_mat[id][0];
            Amp00[id0][1] = Amp_mat[id][1];
        }
    }

    for (int col = 0; col < n0; col++) {
        Uin0[col][0] = Uin_is[col][0];
        Uin0[col][1] = Uin_is[col][1];

        Uout0[col][0] = 0.0f;
        Uout0[col][1] = 0.0f;
    }

    butterfly_apply_1d_phase_amp_core(n0,
                                      tau00,
                                 Amp00,
                                 Uin0,
                                 Uout0,
                                 omega,
                                 p,
                                 leaf_n);

    for (int row = 0; row < n0; row++) {
        Uout_is[row][0] = Uout0[row][0];
        Uout_is[row][1] = Uout0[row][1];
    }

    fftwf_free(Uout0);
    fftwf_free(Uin0);
    fftwf_free(Amp00);
    free(tau00);

    /*
     * If n is already legal, K00 covers the whole matrix.
     */
    if (n0 == n) {
        return;
    }

    /*
     * 2. Top-right block K01: rows 0:n0-1, cols n0:n-1.
     *
     * This contribution is added to the butterfly result in Uout_is[0:n0-1].
     */
    direct_add_1d_phase_amp_range(n,
                                  0,
                                  n0,
                                  n0,
                                  n,
                                  tau_mat,
                                  Amp_mat,
                                  Uin_is,
                                  Uout_is,
                                  omega);

    /*
     * 3. Bottom blocks K10 and K11: rows n0:n-1, cols 0:n-1.
     *
     * These rows are computed exactly by direct summation.
     */
    direct_add_1d_phase_amp_range(n,
                                  n0,
                                  n,
                                  0,
                                  n,
                                  tau_mat,
                                  Amp_mat,
                                  Uin_is,
                                  Uout_is,
                                  omega);
}
