#ifndef SE_PAR_GETSET_H
#define SE_PAR_GETSET_H
#include <SEBASIC/include/se_basic.h>
#include "se_fs_sep.h"
#include "se_par_sep.h"
#include "sesi_map.h"
#include "se_coordinate_system_par.h"



/**
 * Return a logical as a result of testing if a parameter is defined.
 * Try to get the parameter from sep. If sep=NULL, try to get the
 * parameter from map.  If map=NULL, then use the command line. If 
 * the parameter group is specified, also check if grp:p is a defined.
 *
 * \param[in] sep SEP file
 * \param[in] map map file
 * \param[in] grp parameter group name.
 * \param[in] p parameter name
 *
 * \return non-zero if the parameter is defined; 0 otherwise
*/
int se_have_par_int(const sep_t *sep, const sesi_map_t map,            //测试某一参数是否被定义
                      const char *grp, const char *p);

/**
 * Return a logical as a result of testing if a parameter is defined.
 * Try to get the parameter from sep. If sep=NULL, try to get the
 * parameter from map.  If map=NULL, then use the command line. If 
 * the parameter group is specified, also check if grp:p is a defined.
 *
 * \param[in] sep SEP file
 * \param[in] map map file
 * \param[in] grp parameter group name.
 * \param[in] p parameter name
 *
 * \return non-zero if the parameter is defined; 0 otherwise
*/
int se_have_par_float(const sep_t *sep, const sesi_map_t map,          //同上 测试单精度
                        const char *grp, const char *p);

/**
 * Return a logical as a result of testing if a parameter is defined.
 * Try to get the parameter from sep. If sep=NULL, try to get the
 * parameter from map.  If map=NULL, then use the command line. If 
 * the parameter group is specified, also check if grp:p is a defined.
 *
 * \param[in] sep SEP file
 * \param[in] map map file
 * \param[in] grp parameter group name.
 * \param[in] p parameter name
 *
 * \return non-zero if the parameter is defined; 0 otherwise
*/
int se_have_par_str(const sep_t *sep, const sesi_map_t map,            //测试字符串
                      const char *grp, const char *p);

/**
 * Return an integer value of a parameter.  Try to get the parameter
 * from sep. If sep=NULL, try to get the parameter from map.  If
 * map=NULL, then use the command line.  If the parameter group is
 * specified and p is not a defined parameter then also try grp:p as
 * the parameter name. If the parameter is not defined, return the
 * default value.
 *
 * \param[in] sep SEP file
 * \param[in] map map file
 * \param[in] grp parameter group name
 * \param[in] p parameter name
 * \param[in] v default value for the parameter
 *
 * \return the parameter value
*/
int se_get_defpar_int32(const sep_t *sep, const sesi_map_t map, 
                          const char *grp, const char *p, int v);

/**
 * Return a 64 bit integer value of a parameter.  Try to get the
 * parameter from sep. If sep=NULL, try to get the parameter from map.
 * If map=NULL, then use the command line.  If the parameter group is
 * specified and p is not a defined parameter then also try grp:p as
 * the parameter name. If the parameter is not defined, return the
 * default value.
 *
 * \param[in] sep SEP file
 * \param[in] map map file
 * \param[in] grp parameter group name
 * \param[in] p parameter name
 * \param[in] v default value for the parameter
 *
 * \return the parameter value
*/
int64_t se_get_defpar_int64(const sep_t *sep, const sesi_map_t map, 
                              const char *grp, const char *p, int64_t v);

/**
 * Return the value of a parameter.  Try to get the parameter from
 * sep. If sep=NULL, try to get the parameter from map.  If map=NULL,
 * then use the command line.  If the parameter group is specified and
 * p is not a defined parameter then also try grp:p as the parameter
 * name.  If the parameter is not defined, return the default value.
 *
 * \param[in] sep SEP file
 * \param[in] map map file
 * \param[in] grp parameter group name
 * \param[in] p parameter name
 * \param[in] v default value for the parameter
 *
 * \return the parameter value
*/
double se_get_defpar_float(const sep_t *sep, const sesi_map_t map, 
                             const char *grp, const char *p, double v);

/**
 * Return a string value of a parameter.  Try to get the parameter
 * from sep. If sep=NULL, try to get the parameter from map.  If
 * map=NULL, then use the command line.  If the parameter group is
 * specified and p is not a defined parameter then also try grp:p as
 * the parameter name. If the parameter is not defined, return the
 * default value.
 *
 * \param[in] sep SEP file
 * \param[in] map map file
 * \param[in] grp parameter group name
 * \param[in] p parameter name
 * \param[in] v default value for the parameter
 *
 * \return the parameter value
*/
char * se_get_defpar_str(const sep_t *sep, const sesi_map_t map, 
                           const char *grp, const char *p, const char *v);

/**
 * Get the survey description. If sep=NULL, try to get the survey from
 * map.  If map=NULL, then use the command line. If group is specified,
 * fill parameters of the form "grp:prefix.ox_survey" from that as well. 
 * 
 * \param[in] sep SEP file
 * \param[in] map map file
 * \param[in] grp parameter group name
 * \param[in] prefix the prefix for the survey parameters
 * \param[out] sd the survey description to fill
 *
 */
void se_get_survey_description(const sep_t *sep, const sesi_map_t map,         //获得survey description
                                 const char* grp, const char* prefix, 
                                 survey_description_t* sd);

/**
 * Add an integer parameter to a SEP and/or map header, if their
 * pointers are non NULL. If group is non-NULL use grp:p for the
 * parameter name.
 *
 * \param[in] sep SEP file
 * \param[in] map map file
 * \param[in] grp parameter group name
 * \param[in] p parameter name
 * \param[in] v value for the parameter
 */
void se_set_header_int(sep_t *sep, sesi_map_t map, 
                         const char *grp, const char *p, int v);

/**
 * Add a parameter to a SEP and/or map header, if their pointers are
 * non NULL. If group is non-NULL use grp:p for the parameter name.
 *
 * \param[out] sep SEP file
 * \param[out] map map file
 * \param[in] grp parameter group name
 * \param[in] p parameter name
 * \param[in] v value for the parameter
 */
void se_set_header_float(sep_t *sep, sesi_map_t map, 
                           const char *grp, const char *p, double v);

/**
 * Add a string parameter to a SEP and/or map header, if their
 * pointers are non NULL. If group is non-NULL use grp:p for the
 * parameter name.
 *
 * \param[out] sep SEP file
 * \param[out] map map file
 * \param[in] grp parameter group name
 * \param[in] p parameter name
 * \param[in] v value for the parameter
 */
void se_set_header_str(sep_t *sep, sesi_map_t map,
                         const char *grp, const char *p, const char *v);

/**
 * Add a survey description to a SEP and/or map header, if their
 * pointers are non NULL. If group is non-NULL use grp:prefix for the 
 * parameter name.
 *
 * \param[out] sep SEP file
 * \param[out] map map file
 * \param[in] grp parameter group name
 * \param[in] prefix the prefix for the survey parameters
 * \param[in] sd the survey description
 *
 */
void se_set_header_survey_description(sep_t *sep, sesi_map_t map,
                                        const char* grp,
                                        const char* prefix, 
                                        const survey_description_t* sd);

/**
 * Strip the group of a parameter. If the parameter is already defined,
 * do nothing otherwise check if grp:par is defined and assign its value
 * to par.
 *
 * \param[in] grp the group name of the parameter
 * \param[in] par the name of the parameter
 *
 */
void se_strip_par_group(const char *grp, const char *par);


/**
 * Find the version info from the selfdoc string and set a parameter
 * to it. 
 *
 * \param[in] par the name of the parameter that will contain the version
 * \param[in] selfdoc the selfdoc string
 *
 */
void se_set_version_from_selfdoc(const char *par, const char *selfdoc);  //查找版本信息 设置参数

#endif