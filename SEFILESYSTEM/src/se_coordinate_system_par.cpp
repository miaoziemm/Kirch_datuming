#include "../include/se_coordinate_system_par.h"


static int _se_get_project_survey_description(survey_description_t* sd)
{
    char* bdir;
    int ret;

    /**\todo get the job_dir from JCS; solve dependency problems */
    if( se_have_par("job_dir") ) {
        bdir =  se_get_par_str("job_dir");
    } else {
        bdir =  get_cwd();
    }

    ret = se_get_associated_survey_description( sd, bdir, NULL, NULL);
    free(bdir);
    return ret;
}

void se_get_survey_description(survey_description_t* sd)
{
    se_get_survey_description_from_pars(sd, NULL);
}

void se_get_survey_description_from_pars(survey_description_t* sd, se_hash pars)
{

    int have_prj = _se_get_project_survey_description(sd);
    se_get_named_survey_description_from_pars(sd, NULL, NULL, pars, !have_prj);
}

void
se_get_named_survey_description( survey_description_t* sd,
                                  const char* name,
                                  const survey_description_t* parent_sd )
{
    se_get_named_survey_description_from_pars(sd, name, parent_sd, NULL, 1);
}

void
se_get_named_survey_description_from_pars( survey_description_t* sd,
                                            const char* name,
                                            const survey_description_t* parent_sd,
                                            se_hash pars,
                                            int overwrite)
{
    char parname[256];
    char* dot;
    survey_description_t default_parent;

    if(overwrite) {
        init_survey_description(sd);
    }

    if( name == NULL ||
        name[0] == 0 ||
        0 == strcmp("default", name) ) {
        name = "";
        dot = const_cast<char*>("");
        snprintf(sd->name, sizeof(sd->name)-1, "default");
    } else {
        dot = const_cast<char*>(".");
        if( parent_sd == NULL ) {
            se_get_survey_description_from_pars(&default_parent, pars);
            parent_sd = &default_parent;
        }

        snprintf(parname, sizeof(parname), "%s%s%s", name, dot, "survey_name");
        if( ( pars && se_have_htpar(parname, pars) ) ||
            se_have_par(parname) ) {
            char* sname;
            if(pars) {
                sname = se_get_htpar_str(parname, pars);
            } else {
                sname = se_get_par_str(parname);
            }
            snprintf(sd->name, sizeof(sd->name)-1, "%s", sname);
            sd->has_name = 1;
            free(sname);
        } else if( ( pars && se_have_htpar("survey_name", pars) ) ||
                   se_have_par("survey_name") ) {
            char* sname;
            if(pars) {
                sname = se_get_htpar_str("survey_name", pars);
            } else {
                sname = se_get_par_str("survey_name");
            }
            snprintf(sd->name, sizeof(sd->name)-1, "%s", sname);
            sd->has_name = 1;
            free(sname);
        } else {
            // 确保目标缓冲区有足够的空间防止截断
            size_t available_space = sizeof(sd->name) - 1;
            size_t parent_name_len = strlen(parent_sd->name);
            size_t name_len = strlen(name);
            
            // 检查是否有足够的空间存储合并的名称
            if (parent_name_len + 1 + name_len < available_space) {
                // 确保不会截断
                char combined_name[512]; // 足够大的缓冲区
                snprintf(combined_name, sizeof(combined_name), "%s.%s", parent_sd->name, name);
                strncpy(sd->name, combined_name, sizeof(sd->name)-1);
                sd->name[sizeof(sd->name)-1] = '\0';
            } else {
                // 如果空间不足，则截断名称
                size_t max_parent_len = available_space > (name_len + 1) ? 
                                        available_space - (name_len + 1) : 0;
                char truncated_parent[256];
                strncpy(truncated_parent, parent_sd->name, max_parent_len);
                truncated_parent[max_parent_len] = '\0';
                
                char combined_name[512]; // 足够大的缓冲区
                snprintf(combined_name, sizeof(combined_name), "%s.%s", truncated_parent, name);
                strncpy(sd->name, combined_name, sizeof(sd->name)-1);
                sd->name[sizeof(sd->name)-1] = '\0';
            }
            sd->has_name = 0;
        }
    }


#define _SET_SD_PAR(p, def, type) do {                                  \
        int needs_default = 0;                                          \
        snprintf(parname, sizeof(parname), "%s%s%s", name, dot, #p);    \
        if(pars) {                                                      \
            if( se_have_htpar(parname, pars) ) {                       \
                sd->p = se_get_htpar_##type(parname, pars);            \
                sd->has_##p = 1;                                        \
            } else if( parent_sd ) {                                    \
                sd->p = parent_sd->p;                                   \
                sd->has_##p = 0;                                        \
            } else if( se_have_htpar(#p, pars) ) {                     \
                sd->p = se_get_htpar_##type(#p, pars);                 \
                sd->has_##p = 0;                                        \
            } else {                                                    \
                needs_default = 1;                                      \
            }                                                           \
        } else {                                                        \
            if( se_have_par(parname) ) {                               \
                sd->p = se_get_par_##type(parname);                    \
                sd->has_##p = 1;                                        \
            } else if( parent_sd ) {                                    \
                sd->p = parent_sd->p;                                   \
                sd->has_##p = 0;                                        \
            } else if( se_have_par(#p) ) {                             \
                sd->p = se_get_par_##type(#p);                         \
                sd->has_##p = 0;                                        \
            } else {                                                    \
                needs_default = 1;                                      \
            }                                                           \
        }                                                               \
        if(needs_default) {                                             \
            sd->p = (def);                                              \
            sd->has_##p = 0;                                            \
        }                                                               \
    } while(0)

    _SET_SD_PAR(ox_survey     , sd->ox_survey     , double);
    _SET_SD_PAR(oy_survey     , sd->oy_survey     , double);
    _SET_SD_PAR(survey_azimuth, sd->survey_azimuth, double);
    _SET_SD_PAR(left_handed   , sd->left_handed   , int32);
    _SET_SD_PAR(inline_azimuth, sd->inline_azimuth, int32);

    _SET_SD_PAR(first_inline    , sd->first_inline, double);
    _SET_SD_PAR(inline_at_origin, overwrite?sd->first_inline:sd->inline_at_origin, double);
    _SET_SD_PAR(last_inline     , overwrite?sd->first_inline:sd->last_inline, double);
    _SET_SD_PAR(inline_increment, overwrite?(sd->first_inline <= sd->last_inline ? 1.0 : -1.0):sd->inline_increment, double);
    _SET_SD_PAR(inline_spacing  , sd->inline_spacing, double);

    _SET_SD_PAR(first_crossline    , sd->first_crossline, double);
    _SET_SD_PAR(crossline_at_origin, overwrite?sd->first_crossline:sd->crossline_at_origin, double);
    _SET_SD_PAR(last_crossline     , overwrite?sd->first_crossline:sd->last_crossline, double);
    _SET_SD_PAR(crossline_increment, 
                overwrite?(sd->first_crossline <= sd->last_crossline ? 1.0 : -1.0):sd->crossline_increment, double);
    _SET_SD_PAR(crossline_spacing  , sd->crossline_spacing, double);

    _SET_SD_PAR(survey_units, sd->survey_units, se_survey_units_t);

#undef _SET_SD_PAR
    survey_description_update_cached_values(sd);

    if( fabs(sd->inline_increment) < 1.0e-20 ) {
        ERROR(("%s%sinline_increment cannot be zero: %f",
               name, dot, sd->inline_increment));
    }
    if( fabs(sd->crossline_increment) < 1.0e-20 ) {
        ERROR(("%s%scrossline_increment cannot be zero: %f",
               name, dot, sd->crossline_increment));
    }
    if( sd->inline_spacing < 0.0 ) {
        ERROR(("%s%sinline_spacing cannot be negative: %f",
               name, dot, sd->inline_spacing));
    }
    if( sd->inline_spacing < 1.0E-20 ) {
        ERROR(("%s%sinline_spacing cannot be zero: %f",
               name, dot, sd->inline_spacing));
    }
    if( sd->crossline_spacing < 0.0 ) {
        ERROR(("%s%scrossline_spacing cannot be negative: %f",
               name, dot, sd->crossline_spacing));
    }
    if( sd->crossline_spacing < 1.0E-20 ) {
        ERROR(("%s%scrossline_spacing cannot be zero: %f",
               name, dot, sd->crossline_spacing));
    }
    if( ! is_integer
        ( (sd->last_inline - sd->first_inline) / sd->inline_increment ) ) {
        INFOV((10, "%s%slast_inline(%g) doesn't match the increment(%g) "
               "from %s%sfirst_inline(%g)",
               name, dot, sd->last_inline, sd->inline_increment, name, dot,
               sd->first_inline));
    }
    if( ! is_integer
        ((sd->last_crossline - sd->first_crossline)/sd->crossline_increment)) {
        INFOV((10, "%s%slast_crossline(%g) doesn't match the increment(%g) "
               "from %s%sfirst_crossline(%g)",
               name, dot, sd->last_crossline, sd->crossline_increment,
               name, dot, sd->first_crossline));
    }

    switch(sd->survey_units) {
    case SURVEY_UNITS_METERS:
    case SURVEY_UNITS_FEET:
        break;
    default:
        if(se_get_verb_level() > 0) {
            WARN(("Unknown value for %s%ssurvey_units: %d. Will set to %d (%s)", name, dot,
                  sd->survey_units, 
                  SURVEY_UNITS_METERS, survey_get_units_name(SURVEY_UNITS_METERS)));
        }
        sd->survey_units = SURVEY_UNITS_METERS;
    }
}


int
se_get_associated_survey_description_from_file( survey_description_t* sd,
                                                 const char* sd_file,
                                                 const char* filename,
                                                 const char* name,
                                                 const survey_description_t* parent_sd)
{
    array sdfiles;
    int ret = 1;
    ASSERT(sd);

    /* by default it will return whatever values are specified as parameters */
    se_get_named_survey_description_from_pars( sd, name, parent_sd, NULL, 1);

    if(sd_file == NULL) sd_file = SE_DEFAULT_SURVEY_DESCRIPTION_FILE;
    sdfiles = se_find_file_upwards( sd_file, filename );
    if(array_size(sdfiles) > 0) {
        const char* sdf = (const char*)array_get_at(sdfiles, 0);
        se_hash pars = se_parse_pars_file(sdf, NULL, NULL);
        se_get_named_survey_description_from_pars( sd, name, parent_sd, pars, 0);
        se_destroy_parameters(pars);
        ret = 1;
    } else {
        ret = 0;
    }

    destroy_array(sdfiles, 1);

    return ret;
}
