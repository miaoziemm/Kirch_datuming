#include "../include/se_par_getset.h"


int se_have_par_int(const sep_t *sep, const sesi_map_t map, const char *grp, const char *p) {

    int have_par=0;
    
    if (sep) have_par = have_par || sep_have_hdr_int(sep,p);
    if (map) have_par = have_par || sesi_map_have_hdr_int(map,p);
    have_par = have_par || se_have_par(p);

    if ((grp)&&(!have_par)) {
        char *grp_par;
        grp_par = asprintf("%s:%s",grp,p);
        have_par = se_have_par_int(sep, map, NULL, grp_par);
        free(grp_par);
    }

    return have_par;

}

int se_have_par_float(const sep_t *sep, const sesi_map_t map, const char *grp, const char *p) {

    int have_par=0;
    
    if (sep) have_par = have_par || sep_have_hdr_float(sep,p);
    if (map) have_par = have_par || sesi_map_have_hdr_float(map,p);
    have_par = have_par || se_have_par(p);

    if ((grp)&&(!have_par)) {
        char *grp_par;
        grp_par = asprintf("%s:%s",grp,p);
        have_par = se_have_par_float(sep, map, NULL, grp_par);
        free(grp_par);
    }

    return have_par;

}


int se_have_par_str(const sep_t *sep, const sesi_map_t map, const char *grp, const char *p) {

    int have_par=0;
    
    if (sep) have_par = have_par || sep_have_hdr(sep,p);
    if (map) have_par = have_par || sesi_map_have_hdr(map,p);
    have_par = have_par || se_have_par(p);

    if ((grp)&&(!have_par)) {
        char *grp_par;
        grp_par = asprintf("%s:%s",grp,p);
        have_par = se_have_par_str(sep, map, NULL, grp_par);
        free(grp_par);
    }

    return have_par;

}

int se_get_defpar_int32(const sep_t *sep, const sesi_map_t map, 
                          const char *grp, const char *p, int v)
{
    int ret;
    
    if (sep && sep_have_hdr_int(sep,p)) {
        ret = sep_get_hdr_int(sep,p,v);
        return ret;
    } else if (map && sesi_map_have_hdr_int(map,p)) {
        ret = sesi_map_get_hdr_int(map,p,v);
        return ret;
    } else if (se_have_par(p)) {
        ret = se_get_par_int32(p);
        return ret;
    }

    if (grp) {
        char *grp_par;
        grp_par = asprintf("%s:%s",grp,p);
        ret = se_get_defpar_int32(sep, map, NULL, grp_par, v);
        free(grp_par);
    } else {
        ret = v;
    }

    return ret;        

}

int64_t se_get_defpar_int64(const sep_t *sep, const sesi_map_t map, 
                              const char *grp, const char *p, int64_t v)
{
    int64_t ret;
    
    if (sep && sep_have_hdr_int(sep,p)) {
        ret = sep_get_hdr_int(sep,p,v);
        return ret;
    } else if (map && sesi_map_have_hdr_int(map,p)) {
        ret = sesi_map_get_hdr_int(map,p,v);
        return ret;
    } else if (se_have_par(p)) {
        ret = se_get_par_int64(p);
        return ret;
    }

    if (grp) {
        char *grp_par;
        grp_par = asprintf("%s:%s",grp,p);
        ret = se_get_defpar_int64(sep, map, NULL, grp_par, v);
        free(grp_par);
    } else {
        ret = v;
    }

    return ret;        

}


double se_get_defpar_float(const sep_t *sep, const sesi_map_t map, 
                             const char *grp, const char *p, double v)
{
    double ret;
        
    if (sep && sep_have_hdr_float(sep,p)) {
        ret = sep_get_hdr_float(sep,p,v);
        return ret;
    } else if (map && sesi_map_have_hdr_float(map,p)) {
        ret = sesi_map_get_hdr_float(map,p,v);
        return ret;
    } else if (se_have_par(p)) {
        ret = se_get_par_float(p);
        return ret;
    }

    if (grp) {
        char *grp_par;
        grp_par = asprintf("%s:%s",grp,p);
        ret = se_get_defpar_float(sep, map, NULL, grp_par, v);
        free(grp_par);
    } else {
        ret = v;
    }

    return ret;        

}

char* se_get_defpar_str(const sep_t *sep, const sesi_map_t map, 
                          const char *grp, const char *p, const char *v)
{
    char *ret=NULL;

    if (sep && sep_have_hdr(sep,p)) {
        ret = sep_get_hdr(sep,p,v);
        return ret;
    } else if (map && sesi_map_have_hdr(map,p)) {
        ret = sesi_map_get_hdr(map,p,v);
        return ret;
    } else if (se_have_par(p)) {
        ret = se_get_par_str(p);
        return ret;
    }

    if (grp) {
        char *grp_par;
        grp_par = asprintf("%s:%s",grp,p);
        ret = se_get_defpar_str(sep, map, NULL, grp_par, v);
        free(grp_par);
    } else {
        if (v) {
            ret = strdup(v);
        } else {
            ret = NULL;
        }
    }

    return ret;        

}


void se_get_survey_description(const sep_t *sep, const sesi_map_t map,
                                 const char* grp, const char* prefix, 
                                 survey_description_t* sd) 
{

    char *grp_prefix=NULL;
    survey_description_t tmp_sd;

    init_survey_description(sd);
    se_get_survey_description_from_pars(sd, NULL);

    if (grp) { 
        grp_prefix = asprintf("%s:%s",grp,prefix);
    }
    
    if (map) {
        if (grp_prefix) {
            tmp_sd = *sd;
            sesi_map_get_header_survey_description(map, sd, grp_prefix, &tmp_sd);
        }
        tmp_sd = *sd;
        sesi_map_get_header_survey_description(map, sd, prefix, &tmp_sd);
    }

    if (sep) {
        if (grp_prefix) {
            tmp_sd = *sd;
            sep_get_header_survey_description(sep, sd, grp_prefix, sd, 0);
        }
        tmp_sd = *sd;
        sep_get_header_survey_description(sep, sd, prefix, sd, 0);
    }

    if (grp_prefix) {
        se_get_named_survey_description_from_pars(sd, grp_prefix, sd, NULL, 0);
    }
    se_get_named_survey_description_from_pars(sd, prefix, sd, NULL, 0);

    if (grp_prefix) free(grp_prefix);

}


void se_set_header_int(sep_t *sep, sesi_map_t map, 
                         const char* grp, const char *p, int v) {

    char *par;

    if (grp) {
        par = asprintf("%s:%s",grp,p);
    } else {
        par = asprintf("%s",p);
    }

    if (sep) sep_set_header_int(sep,par,v);
    if (map) sesi_map_set_header_int(map,par,v);

    free(par);
}

void se_set_header_float(sep_t *sep, sesi_map_t map, 
                           const char* grp, const char *p, double v) {

    char *par;

    if (grp) {
        par = asprintf("%s:%s",grp,p);
    } else {
        par = asprintf("%s",p);
    }

    if (sep) sep_set_header_float(sep, par, v);
    if (map) sesi_map_set_header_float(map, par, v);

    free(par);
}

void se_set_header_str(sep_t *sep, sesi_map_t map, 
                         const char* grp, const char *p, const char* v) {

    char *par;

    if (!v) return;

    if (grp) {
        par = asprintf("%s:%s",grp,p);
    } else {
        par = asprintf("%s",p);
    }

    if (sep) sep_set_header(sep,par,v);
    if (map) sesi_map_set_header(map,par,v);

    free(par);
        
}

void se_set_header_survey_description(sep_t *sep, sesi_map_t map,
                                        const char* grp,
                                        const char* prefix, 
                                        const survey_description_t* sd)
{
    char *prf;

    if (grp) {
        prf = asprintf("%s:%s",grp,prefix);
    } else {
        prf = asprintf("%s",prefix);
    }

    if (sep) sep_set_header_survey_description(sep, sd, prf, 0);
    if (map) sesi_map_set_header_survey_description(map, sd, prf, 0);
    
    free(prf);

}


void se_strip_par_group(const char *grp, const char *par)
{
    char *parname;
    if (se_have_par(par)) return;
    parname = asprintf("%s:%s",grp,par);
    if (se_have_par(parname)) {
        char *val = se_get_par_str(parname);
        se_set_par_str(par,val);
        free(val);
    }
    free(parname);
}

void se_set_version_from_selfdoc(const char *par, const char *selfdoc)
{
    char *info;
    char *ver=NULL, *tag=NULL, *blt=NULL, *str=NULL;
    int len;
   
    info = (char *)strstr(selfdoc," Version: ")+10;
    if (info) {
        len=0;
        while ((info[len]!=0)&&(info[len]!=0x0A)) {
            len++;
        }
        ver = strndup(info, len );
    }

    info = (char *)strstr(selfdoc," Tag: ")+6;
    if (info) {
        len=0;
        while ((info[len]!=0)&&(info[len]!=0x0A)) {
            len++;
        }
        tag = strndup(info, len );
    }

    info = (char *)strstr(selfdoc," Built: ")+8;
    if (info) {
        len=0;
        while ((info[len]!=0)&&(info[len]!=0x0A)) {
            len++;
        }
        blt = strndup(info, len );
    }
        
    str = asprintf("Version: %s -- Tag: %s -- Built: %s",ver,tag,blt);
    se_set_par_str(par, str);

    free(ver);
    free(tag);
    free(blt);
    free(str);

    

}
