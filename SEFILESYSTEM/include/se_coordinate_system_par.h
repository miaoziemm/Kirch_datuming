#ifndef SE_COORDINATE_SYSTEM_PAR_H
#define SE_COORDINATE_SYSTEM_PAR_H
#include <SEBASIC/include/se_basic.h>
#include "se_module_par_desc_sep.h"


#define SE_DEFAULT_SURVEY_DESCRIPTION_FILE "se_survey_description"

/**
 * Fills the sd structure with the default survey description
 * available from parameters.
 *
 * The result is the same as calling se_get_named_survey_description
 * with NULL or "default" for name, and NULL for the parent.
 */
void se_get_survey_description( survey_description_t* sd );

/**
 * Fills the sd structure with the survey description found for the
 * "name" argument. If name is NULL or "default", the parameters are
 * read for the "main" or "default" survey description available at
 * run-time.
 *
 * The missing parameters are filled from the parent_sd, if this is
 * different from NULL; otherwise default values will be used.
 *
 * If no origin (ox_survey,oy_survey) is found available for this
 * survey, and no in/crossline_at_origin is found, then the
 * in/crossline_at_origin from the parent is used.
 *
 * If an origin is specified, then in/crossline_at_origin defaults to
 * the first_in/crossline available for this (sub)survey.
 *
 */
void
se_get_named_survey_description( survey_description_t* sd,
                                  const char* name,
                                  const survey_description_t* parent_sd);


void
se_get_named_survey_description_from_pars( survey_description_t* sd,
                                            const char* name,
                                            const survey_description_t* parent_sd,
                                            se_hash pars, int overwrite);

void se_get_survey_description_from_pars(survey_description_t* sd, se_hash pars);


/**
 * Returns true if it found a survey description in something named after sd_file.
 */
int
se_get_associated_survey_description_from_file( survey_description_t* sd,
                                                 const char* sd_file,
                                                 const char* filename,
                                                 const char* name,
                                                 const survey_description_t* parent_sd);

inline int
se_get_associated_survey_description( survey_description_t* sd,
                                       const char* filename,
                                       const char* name,
                                       const survey_description_t* parent_sd)
{
    return se_get_associated_survey_description_from_file(sd, 
                                                           SE_DEFAULT_SURVEY_DESCRIPTION_FILE,
                                                           filename, name, parent_sd);
}



inline survey_units_t se_get_htpar_se_survey_units_t( const char* parname, const se_hash ht )
{
    return survey_get_units_from_int(se_get_htpar_int(parname, ht));
}

inline survey_units_t se_get_par_se_survey_units_t(const char* parname)
{
    return survey_get_units_from_int(se_get_par_int(parname));
}


#endif