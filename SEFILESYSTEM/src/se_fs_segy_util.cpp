#include "../include/se_par_sep.h"
#include "../include/se_fs_segy.h"
#include "../include/se_fs_sep.h"
#include "../include/se_fs_segy_trchead_eval.h"


char* segy_get_first_last_trace_parameters(int ntraces,
                                              int* _first_trace, int* _last_trace)
{
    int first_trace, last_trace;

    if(se_have_par("first_trace")) {
        first_trace = se_get_par_int32("first_trace");
        if(first_trace < 1) {
            return asprintf("first_trace parameter needs to be bigger than 1; "
                                "current value is %d", first_trace);
        }
        if(first_trace > ntraces) {
            return asprintf("first_trace=%d parameter needs to be less than "
                                "total traces found in file: %d", first_trace, ntraces);
        }

        INFOV((1, "Will read starting with trace %d", first_trace));
    } else {
        first_trace = 1;
    }

    first_trace -= 1; /* 0 indexing */

    if(se_have_par("ntraces")) {
        int nt = se_get_par_int32("ntraces");
        if(nt < 0) {
            return asprintf("ntraces=%d cannot be negative", nt);
        }

        if(first_trace + nt > ntraces) {
            WARN(("Specified number of traces (%d) is too big. Will adjust to %d",
                  nt, ntraces-first_trace));
            nt = ntraces-first_trace;
        } else {
            INFOV((1, "Will read %d traces", nt));
        }
        last_trace = first_trace + nt;
    } else if(se_have_par("last_trace")) {
        last_trace = se_get_par_int32("last_trace");
        if(last_trace < first_trace+1) {
            return asprintf("last_trace=%d cannot be less than the first trace", last_trace);
        }
        if(last_trace > ntraces) {
            WARN(("Specified last_trace=%d is too big. Will adjust to %d",
                  last_trace, ntraces));
            last_trace = ntraces;
        } else {
            INFOV((1, "Last trace to read will be %d", last_trace));
        }
    } else {
        last_trace = ntraces;
    }
    last_trace -= 1; /* 0 indexing */


    *_first_trace = first_trace;
    *_last_trace = last_trace;

    return NULL;
}

static int _have_par(const char* pname)
{
    if(se_have_par(pname)) {
        char* str = se_get_par_str(pname);
        int have = (str[0] != 0);
        free(str);
        return have;
    }
    return 0;
}

static int32_t _get_defpar_int32(const char* pname, int def)
{
    if(!_have_par(pname)) {
        return def;
    }
    return se_get_par_int32(pname);
}
static char* _get_defpar_str(const char* pname, const char* def)
{
    if(!_have_par(pname)) {
        return se_strdup(def);
    }
    return se_get_par_str(pname);
}

static int _have_parg(const char* gname, const char* pname)
{
    int ret = 0;
    char* name = asprintf("%s.%s", gname, pname);
    if( _have_par(name) ) {
        ret = 1;
    } else {
        ret = _have_par(pname);
    }

    free(name);

    return ret;
}

static double _get_defparg_double(const char* gname, const char* pname, double dval)
{
    double ret = dval;
    char* name = asprintf("%s.%s", gname, pname);
    if( _have_par(name) ) {
        ret = se_get_par_double(name);
    } else {
        ret = se_get_defpar_double(pname, dval);
    }
    free(name);

    return ret;
}

static double _get_parg_double(const char* gname, const char* pname)
{
    return _get_defparg_double(gname, pname, 0);
}

void segy_bin_check_and_convert_friendly_pars(void)
{
    int ilineidx = 3, xlineidx = 2;
    if( se_get_par_int32("gathers") && _have_par("offset_header_byte") && !_have_par("f2")) {
        char* fp = _get_defpar_str("offset_header_type", "i");
        char* f = asprintf("%s%d", fp, _get_defpar_int32("offset_header_byte", 37));
        se_set_par_str("f2", f);
        INFOV((2, "Adding pars: f2=%s", f));
        free(f);
        free(fp);
        if(_have_par("offset_increment") && !_have_par("d2")) {
            double v = se_get_par_double("offset_increment");
            se_set_par_double("d2", v);
            INFOV((2, "Adding pars: d2=%f", v));
        }
        if(_have_par("min_offset") && !_have_par("o2")) {
            double v = se_get_par_double("min_offset");
            se_set_par_double("o2", v);
            INFOV((2, "Adding pars: o2=%f", v));
        }
        if(_have_par("max_offset") && !_have_par("n2")) {
            double v1 = se_get_par_double("max_offset");
            double v0 = se_get_par_double("o2");
            double d  = se_get_par_double("d2");
            int n = (int)ceil(fabs(v1-v0)/fabs(d))+1;
            se_set_par_int32("n2", n);
            INFOV((2, "Adding pars: n2=%d", n));
        }

        xlineidx = 3;
        ilineidx = 4;
        if(!_have_par("label2")) {
            se_set_par_str("label2", "offset");
            INFOV((2, "Adding pars: label2=offset"));
        }
    }

    if( _have_par("crossline_header_byte") ) {
        char* par = asprintf("f%d", xlineidx);
        if(!_have_par(par)) {
            double atorig;
            double d, o, inc = _get_defparg_double("segy_bin", "crossline_increment", 1);
            char* fp = _get_defpar_str("crossline_header_type", "i");
            char* f = asprintf("%s%d", fp, _get_defpar_int32("crossline_header_byte", 181));
            se_set_par_str(par, f);
            INFOV((2, "Adding pars: %s=%s", par, f));
            free(f);
            free(fp);

            par[0] = 'd';
            if(_have_parg("segy_bin", "crossline_spacing") && !_have_par(par)) {
                char* pout = asprintf("dout%d", xlineidx);
                d = inc*_get_parg_double("segy_bin", "crossline_spacing");
                se_set_par_double(par, inc);
                se_set_par_double(pout, d);
                INFOV((2, "Adding pars: %s=%f, %s=%f", par, inc, pout, d));
                free(pout);
            } else {
                if(_have_par(par)) {
                    d = se_get_par_double(par);
                } else {
                    d = 1;
                    ERROR(("crossline spacing must be specified in this case"));
                }
            }

            par[0] = 'o';
            if(_have_parg("segy_bin", "first_crossline") && !_have_par(par)) {
                char* pout = asprintf("oout%d", xlineidx);
                double f1 = _get_parg_double("segy_bin", "first_crossline");
                atorig = _get_defparg_double("segy_bin", "crossline_at_origin", f1);
                o = (f1-atorig)*d/inc;
                se_set_par_double(par, f1);
                se_set_par_double(pout, o);
                INFOV((2, "Adding pars: %s=%f, %s=%f", par, f1, pout, o));
                free(pout);
            } else {
                o = se_get_par_double(par);
                atorig = o;
            }

            par[0] = 'n';
            if(_have_parg("segy_bin", "last_crossline") && !_have_par(par)) {
                double l = _get_parg_double("segy_bin", "last_crossline");
                int n = (int)ceil( fabs((l-atorig)/inc - o/d) )+1;
                se_set_par_int32(par, n);
                INFOV((2, "Adding pars: %s=%d", par, n));
            }

            {
                char* plabel = asprintf("label%d", xlineidx);
                if(!_have_par(plabel)) {
                    se_set_par_str(plabel, "crosslines");
                    INFOV((2, "Adding pars: %s=crosslines", plabel));
                }
                free(plabel);
            }
        }
        free(par);
    }

    if( _have_par("inline_header_byte") && _have_parg("segy_bin", "first_inline")) {
        char* par = asprintf("f%d", ilineidx);
        if(!_have_par(par)) {
            double atorig;
            double d, o, inc = _get_defparg_double("segy_bin", "inline_increment", 1);
            char* fp = _get_defpar_str("inline_header_type", "i");
            char* f = asprintf("%s%d", fp, _get_defpar_int32("inline_header_byte", 181));
            se_set_par_str(par, f);
            INFOV((2, "Adding pars: %s=%s", par, f));
            free(f);
            free(fp);

            par[0] = 'd';
            if(_have_parg("segy_bin", "inline_spacing") && !_have_par(par)) {
                char* pout = asprintf("dout%d", ilineidx);
                d = inc*_get_parg_double("segy_bin", "inline_spacing");
                se_set_par_double(par, inc);
                se_set_par_double(pout, d);
                INFOV((2, "Adding pars: %s=%f, %s=%f", par, inc, pout, d));
                free(pout);
            } else {
                if(_have_par(par)) 
                    d = se_get_par_double(par);
                else {
                    d = 1;
                    ERROR(("inline spacing must be specified in this case"));
                }
            }

            par[0] = 'o';
            if(_have_parg("segy_bin", "first_inline") && !_have_par(par)) {
                char* pout = asprintf("oout%d", ilineidx);
                double f1 = _get_parg_double("segy_bin", "first_inline");
                atorig = _get_defparg_double("segy_bin", "inline_at_origin", f1);
                o = (f1-atorig)*d/inc;
                se_set_par_double(par, f1);
                se_set_par_double(pout, o);
                INFOV((2, "Adding pars: %s=%f, %s=%f", par, f1, pout, o));
                free(pout);
            } else {
                o = se_get_par_double(par);
                atorig = o;
            }

            par[0] = 'n';
            if(_have_parg("segy_bin", "last_inline") && !_have_par(par)) {
                double l = _get_parg_double("segy_bin", "last_inline");
                int n = (int)ceil( fabs((l-atorig)/inc - o/d) )+1;
                se_set_par_int32(par, n);
                INFOV((2, "Adding pars: %s=%d", par, n));
            }

            {
                char* plabel = asprintf("label%d", ilineidx);
                if(!_have_par(plabel)) {
                    se_set_par_str(plabel, "inlines");
                    INFOV((2, "Adding pars: %s=inlines", plabel));
                }
                free(plabel);
            }
        }
        free(par);
    }
}

