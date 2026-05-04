#include "../include/se_coh_func.h"

float se_coh_func_semblance(const float *data, int n_tr, int n_samp)
{
    float num = 0.0;
    float den = 0.0;
    float temp;

    for (int i = 0; i < n_samp; i++)
    {
        temp = 0.0;
        for (int j = 0; j < n_tr; j++)
        {
            temp += data[j * n_samp + i];
        }
        num += temp * temp;
    }

    for (int i = 0; i < n_samp; i++)
    {
        for (int j = 0; j < n_tr; j++)
        {
            den += data[j * n_samp + i] * data[j * n_samp + i];
        }
    }

    if (FISZERO(den))
        return 0.0f;

    return 1.0f - num / (den * n_tr);


}


float se_coh_func_variation(const float *data, int n_tr, int n_samp)
{
    // NOTE: This is a "variation/residual"-style coherence metric.
    // Input `data` is assumed to be an already-aligned gather of size [n_tr x n_samp]
    // stored trace-major (data[n*n_samp + i]).

    if (!data || n_tr <= 0 || n_samp <= 0)
        return 1.0f;

    // Use double for accumulation stability.
    double *stack = (double *)calloc((size_t)n_samp, sizeof(double));
    double *tr_ene = (double *)calloc((size_t)n_tr, sizeof(double));
    double *diff_norm = (double *)calloc((size_t)n_tr, sizeof(double));

    if (!stack || !tr_ene || !diff_norm)
    {
        free(stack);
        free(tr_ene);
        free(diff_norm);
        return 1.0f;
    }

    // Compute per-trace energy and stack (sum) trace.
    for (int n = 0; n < n_tr; n++)
    {
        const float *tr = data + n * n_samp;
        for (int i = 0; i < n_samp; i++)
        {
            const double v = (double)tr[i];
            tr_ene[n] += v * v;
            stack[i] += v;
        }
    }

    // Convert stack to mean stack and compute its energy.
    double st_ene = 0.0;
    for (int i = 0; i < n_samp; i++)
    {
        stack[i] /= (double)n_tr;
        st_ene += stack[i] * stack[i];
    }

    if (FISZERO(st_ene))
    {
        free(stack);
        free(tr_ene);
        free(diff_norm);
        return 1.0f;
    }

    // For each trace, compute a normalized residual to the stack after best scaling (alpha).
    double res = 0.0;
    for (int n = 0; n < n_tr; n++)
    {
        if (FISZERO(tr_ene[n]))
        {
            diff_norm[n] = 1.0;
        }
        else
        {
            const float *tr = data + n * n_samp;

            double alpha = 0.0;
            for (int i = 0; i < n_samp; i++)
                alpha += (double)tr[i] * stack[i];

            alpha = fabs(alpha / st_ene);

            double dn = 0.0;
            for (int i = 0; i < n_samp; i++)
            {
                const double val = (double)tr[i] - alpha * stack[i];
                dn += val * val;
            }

            // Normalize residual to [0,1] (same form as se_trace_bundle_dip_residual_var).
            dn = 0.5 * dn / (alpha * alpha * st_ene + tr_ene[n]);
            if (dn > 1.0)
                dn = 1.0;

            diff_norm[n] = dn;
        }

        res += diff_norm[n];
    }

    res /= (double)n_tr;
    if (res > 1.0)
        res = 1.0;

    free(stack);
    free(tr_ene);
    free(diff_norm);

    return (float)res;
}
