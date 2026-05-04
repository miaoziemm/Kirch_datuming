#include "../include/se_agc.h"


struct se_agc_s {
    int window;
    int dwind;
    int detect;
    real threshold;

    reala* buff;
    int crt_buff_size;
};

static void se_agc_detect_signal_range( const se_agc_t* agci, 
                                              const float* data, 
                                              int nsamples, 
                                              int* beg, int* end );

se_agc_t* se_agc_init(int window, int dwind, int detect, real threshold)
{
    se_agc_t* agci = alloc1type_zero(se_agc_t,1);


    agci->window    = window;
    agci->dwind     = dwind;
    agci->detect    = detect;
    agci->threshold = threshold;
    
    if( 0.0 > agci->threshold || 
        1.0 <= agci->threshold ) {
        WARN(("AGC detect threshold=%g is out of range (0,1]. Will set it to 0.025",
              agci->threshold ));
        agci->threshold = 0.025;
    }

    agci->buff = NULL;
    agci->crt_buff_size = 0;

    return agci;
}

se_agc_t* se_agc_clone(const se_agc_t* src)
{
    return se_agc_init(src->window, src->dwind, src->detect, src->threshold);
}

void se_agc_destroy( se_agc_t* agci )
{
    if (NULL != agci) {
        if( NULL != agci->buff )
            free1float(agci->buff);

        free(agci);
    }
}

void se_agc_get_pars( const se_agc_t* agc, int *window, int *dwind, int *detect, real *threshold)
{
    if (window)    (*window)     = agc->window;
    if (dwind)     (*dwind)      = agc->dwind;
    if (detect)    (*detect)     = agc->detect;
    if (threshold) (*threshold)  = agc->threshold;

}


void se_apply_agc( se_agc_t* agci, reala* data, int nsamples )
{
    int i, beg, end;
    int wnd = agci->window;
    real norm;
    reala* x;
    int half_wnd;

    remove_nans(data, nsamples, 1, 0.0f);

    if( agci->detect ) {
        se_agc_detect_signal_range(agci, data, nsamples, &beg, &end);
    } else {
        beg = 0;
        end = nsamples - 1;
    }
    if(wnd >= nsamples) wnd = nsamples-1;

    if(end - wnd + 1 <= beg) return;

    /* make sure we have enough memory for this trace */
    if( agci->crt_buff_size < nsamples ) {
        agci->crt_buff_size = nsamples;
        reala *new_buff = (reala *)realloc(agci->buff, agci->crt_buff_size * sizeof(reala));
        if (new_buff == NULL) {
            ERROR(("Failed to allocate memory for AGC buffer"));
            return;
        }
        agci->buff = new_buff;
    }
    x = agci->buff; /* makes the code cleaner */

    /* save input data here */
    memcpy( x, data, nsamples*sizeof(reala) );

    norm = 0.0;
    for(i = beg; i < wnd; i++) {
        float v = x[i];
        if(v < 0.0) v = -v;
        norm += v;
    }
    half_wnd = wnd / 2;

    /* use same norm for the wnd/2 samples from the beginning of trace */
    for(i = beg; i < half_wnd; i++) {
        data[i] *= (reala)( norm / ( norm*norm + EPSILON ) );
    }

    /* compute the norm as we move to the end of trace */
    for(i = half_wnd; i < end - half_wnd; i++) {
        reala v1 = x[i-half_wnd];
        reala v2 = x[i+half_wnd];
        if( v1 < 0.0 ) v1 = -v1;
        if( v2 < 0.0 ) v2 = -v2;
        norm += (v2 - v1);

        data[i] *= (reala)( norm / ( norm*norm + EPSILON ) );
    }

    /* use the last computed norm for the last wnd/2 samples */
    for(i = end - half_wnd; i <= end; i++) {
        data[i] *= (reala)( norm / ( norm*norm + EPSILON ) );
    }
}

static void se_agc_detect_signal_range( const se_agc_t* agci,
                                              const float* data,
                                              int nsamples,
                                              int* beg, int* end )
{
    int i;
    int dw = agci->dwind;
    real init, probe;

    if(dw > nsamples) dw = nsamples;

    init = 0.0;
    for( i = 0; i < dw; ++i ) {
        reala v = data[i];
        if( v < 0.0 ) v = -v;
        init += v;
    }
    probe = init;
    for( i = dw/2; i < nsamples && probe <= init * agci->threshold ; ++i ) {
        reala v = data[i];
        if( v < 0.0 ) v = -v;
        probe += v;
    }

    *beg = i - dw/2;
    
    init = 0.0;
    for( i = nsamples; i >= nsamples - dw; --i ) {
        reala v = data[i];
        if( v < 0.0 ) v = -v;
        init += v;
    }

    probe = init;
    for( i = nsamples - dw/2;
         i >= 0 && probe <= init * agci->threshold ; --i ) {
        reala v = data[i];
        if( v < 0.0 ) v = -v;
        probe += v;
    }
    
    *end = i + dw/2;

    /* peace of mind boundary check */
    if( *beg >= nsamples ) *beg = nsamples - 1;
    if( *end >= nsamples ) *end = nsamples - 1;
    if( *beg < 0 ) *beg = 0;
    if( *end < 0 ) *end = 0;
}
