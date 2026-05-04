#ifndef SE_PAR_SEP_H
#define SE_PAR_SEP_H

#include <sys/time.h>
#include <SEBASIC/include/se_basic.h>
#include "se_uri.h"
#include "se_par_sep.h"

typedef struct
{
    char*     original;     /* original par value (may include variables) */
    array val_array;    /* parsed par values */
    char*     host;         /* par hostname */
    int       src_type;     /* one of PAR_SRCTYPE_* */
    char*     src_info;     /* par src details (line no, cmd-line index, etc) */
    int       update_count; /* par update count */
} par_t;

typedef struct
{
    array par_array;    /* per host param entries */
} parinfo_t;

extern se_hash _ht_params;  //外部变量

extern se_hash _ht_sys_params; //外部变量



int se_par_initialized(void);

/** Initialize the default parameter list. If check_parafile is set,
    it prints a warning if no parameter file is found in command
    line */
int se_par_init_check_parfile(int argn, char** argv, int check_parfile);

/** Initialize the default parameter list. Equivalent to calling
    se_par_init_check_parfile with check_parfile=0. It also checks out
    the default license, if automatic checking is enabled */
inline int se_par_init(int argn, char** argv)
{
    int ret = se_par_init_check_parfile(argn, argv, 0);

    return ret;
}


void se_par_destroy_internal(void);
/** Releases any resources taken by the parameter lib. */
inline void se_par_destroy(void)
{
    se_par_destroy_internal();
}

void se_destroy_parameters( se_hash ht );

/** Read and parse a parameter file. */
se_hash se_parse_pars_file( const char* file_name,
                              se_hash ht, const char* prefix );

se_hash se_parse_pars_str(const char* str,
                            se_hash ht,
                            const char* prefix);

void se_set_pars_from_file(const char* fn, int set_user, int set_system);

void se_set_pars_from_strdump(const char* str, int set_user, int set_system);

char* se_dump_htpars( const se_hash pars, const char* prefix );

char* se_dump_pars( int dump_app_pars, int dump_se_pars );

void se_str_to_pars( const char* str, int erase_pars, int overwrite );
char* se_pars_to_str( int add_user_pars, int add_sys_pars );

int se_have_namedhtpar( const char* parset,
                         const char* parname,
                         se_hash ht,
                         int strict );

int se_have_namedse_par( const char* parset, const char* parname, int strict );

int se_have_namedpar( const char* parset, const char* parname, int strict );

int se_have_htpar( const char* parname, se_hash ht );

int se_have_se_par( const char* parname );
int se_have_ztrpar( const char* parname );
int se_have_par( const char* parname );

void se_set_namedhtpar_str( const char* parset,
                             const char* parname,
                             se_hash ht,
                             const char* val );

void se_set_namedhtpar_double( const char* parset,
                                const char* parname,
                                se_hash ht,
                                double val );
void se_set_namedhtpar_float( const char* parset,
                               const char* parname,
                               se_hash ht,
                               float val );

void se_set_namedhtpar_int64( const char* parset,
                               const char* parname,
                               se_hash ht,
                               int64_t val );

void se_set_namedhtpar_int32( const char* parset,
                               const char* parname,
                               se_hash ht,
                               int32_t val );

inline void se_set_namedhtpar_int( const char* parset,
                                        const char* parname,
                                        se_hash ht,
                                        int val )
{
    se_set_namedhtpar_int32( parset, parname, ht, (int32_t)val );
}

void se_set_namedpar_str( const char* parset,
                              const char* parname,
                              const char* val );

void se_set_namedpar_double( const char* parset,
                                 const char* parname,
                                 double val );

void se_set_namedpar_float( const char* parset,
                                const char* parname,
                                float val );

void se_set_namedpar_int64( const char* parset,
                                const char* parname,
                                int64_t val );

void se_set_namedpar_int32( const char* parset,
                                const char* parname,
                                int32_t val );

inline void se_set_namedse_par_int( const char* parset,
                                         const char* parname,
                                         int val )
{
    se_set_namedpar_int32( parset, parname, (int32_t)val );
}

void se_set_namedpar_str( const char* parset,
                           const char* parname,
                           const char* val );

void se_set_namedpar_double( const char* parset,
                              const char* parname,
                              double val );

void se_set_namedpar_float( const char* parset,
                             const char* parname,
                             float val );

void se_set_namedpar_int64( const char* parset,
                             const char* parname,
                             int64_t val );

void se_set_namedpar_int32( const char* parset,
                             const char* parname,
                             int32_t val );

inline void se_set_namedpar_int( const char* parset,
                                      const char* parname,
                                      int val )
{
    se_set_namedpar_int32( parset, parname, (int32_t)val);
}

void se_set_htpararray_str( const char* parname,
                             se_hash ht,
                             const char** data,
                             int size );

void se_set_htpararray_double( const char* parname,
                                se_hash ht,
                                const double* data,
                                int size );

void se_set_htpararray_real( const char* parname,
                              se_hash ht,
                              const real* data,
                              int size );

void se_set_htpararray_float( const char* parname,
                               se_hash ht,
                               const float* data,
                               int size );

void se_set_htpararray_int64( const char* parname,
                               se_hash ht,
                               const int64_t* data,
                               int size );

void se_set_htpararray_int32( const char* parname,
                               se_hash ht,
                               const int32_t* data,
                               int size );

inline void se_set_htpararray_int( const char* parname,
                                        se_hash ht,
                                        const int32_t* data,
                                        int size )
{
    se_set_htpararray_int32( parname, ht, data, size );
}

void se_set_se_pararray_double( const char* parname,
                                 const double* data,
                                 int size );

void se_set_pararray_str( const char* parname,
                              const char** data,
                              int size );

void se_set_pararray_float( const char* parname,
                                const float* data,
                                int size );

void se_set_pararray_int32( const char* parname,
                                const int32_t* data,
                                int size );

void se_set_pararray_int64( const char* parname,
                                const int64_t* data,
                                int size );
void se_set_pararray_int32( const char* parname,
                             const int32_t* data,
                             int size );
inline void se_set_pararray_int( const char* parname,
                                      const int32_t* data,
                                      int size )
{
    se_set_pararray_int32( parname, data, size );
}
void se_set_pararray_int64( const char* parname,
                             const int64_t* data,
                             int size );
void se_set_pararray_float( const char* parname,
                             const float* data,
                             int size );
void se_set_pararray_double( const char* parname,
                              const double* data,
                              int size );
void se_set_pararray_str( const char* parname,
                           const char** data,
                           int size );

void se_set_htpar_int32( const char* parname, se_hash ht, int32_t val );
inline void se_set_htpar_int( const char* parname, se_hash ht, int val )
{
    se_set_htpar_int32(parname, ht, (int32_t)val );
}
void se_set_htpar_int64( const char* parname, se_hash ht, int64_t val );
void se_set_htpar_float( const char* parname, se_hash ht, float val );
void se_set_htpar_double( const char* parname, se_hash ht, double val );
void se_set_htpar_real( const char* parname, se_hash ht, real val );
void se_set_htpar_reala( const char* parname, se_hash ht, reala val );
void se_set_htpar_str( const char* parname, se_hash ht, const char* val );

// void se_set_par_int32( const char* parname, int32_t val );
// inline void se_set_par_int( const char* parname, int val )
// {
//     se_set_par_int32( parname, (int32_t)val );
// }

void se_set_par_int64( const char* parname, int64_t val );
void se_set_par_float( const char* parname, float val );
void se_set_par_double( const char* parname, double val );
void se_set_par_str( const char* parname, const char* val );

void se_set_par_int32( const char* parname, int32_t val );
inline void se_set_par_int( const char* parname, int val )
{
    se_set_par_int32(parname, (int32_t)val);
}

// void se_set_par_int64( const char* parname, int64_t val );
// void se_set_par_float( const char* parname, float val );
// void se_set_par_double( const char* parname, double val );
// void se_set_par_str( const char* parname, const char* val );
void se_set_par_real( const char* parname, real val );
void se_set_par_reala( const char* parname, reala val );


char* se_get_defnamedhtpar_str( const char* parset,
                                 const char* parname,
                                 const se_hash ht,
                                 const char* def,
                                 int strict );
float se_get_defnamedhtpar_float( const char* parset,
                                   const char* parname,
                                   const se_hash ht,
                                   float def,
                                   int strict );
double se_get_defnamedhtpar_double( const char* parset,
                                     const char* parname,
                                     const se_hash ht,
                                     double def,
                                     int strict );

int32_t se_get_defnamedhtpar_int32( const char* parset,
                                     const char* parname,
                                     const se_hash ht,
                                     int32_t def,
                                     int strict );
inline int se_get_defnamedhtpar_int( const char* parset,
                                          const char* parname,
                                          const se_hash ht,
                                          int def,
                                          int strict )
{
    return (int)se_get_defnamedhtpar_int32( parset, parname, ht,  (int32_t) def, strict );
}

int64_t se_get_defnamedhtpar_int64( const char* parset,
                                     const char* parname,
                                     const se_hash ht,
                                     int64_t def,
                                     int strict );

char* se_get_defnamedse_par_str( const char* parset,
                                  const char* parname,
                                  const char* def,
                                  int strict );

float se_get_defnamedse_par_float( const char* parset,
                                    const char* parname,
                                    float def,
                                    int strict );
double se_get_defnamedse_par_double( const char* parset,
                                      const char* parname,
                                      double def,
                                      int strict );

int32_t se_get_defnamedse_par_int32( const char* parset,
                                      const char* parname,
                                      int32_t def,
                                      int strict );
inline int se_get_defnamedse_par_int( const char* parset,
                                           const char* parname,
                                           int def,
                                           int strict )
{
    return (int)se_get_defnamedse_par_int32( parset, parname, (int32_t)def, strict);
}

int64_t se_get_defnamedse_par_int64( const char* parset,
                                      const char* parname,
                                      int64_t def,
                                      int strict );

char* se_get_defnamedpar_str( const char* parset,
                               const char* parname,
                               const char* def,
                               int strict );

float se_get_defnamedpar_float( const char* parset,
                                 const char* parname,
                                 float def,
                                 int strict );
double se_get_defnamedpar_double( const char* parset,
                                   const char* parname,
                                   double def,
                                   int strict );

int32_t se_get_defnamedpar_int32( const char* parset,
                                   const char* parname,
                                   int32_t def,
                                   int strict );
inline int se_get_defnamedpar_int( const char* parset,
                                        const char* parname,
                                        int def,
                                        int strict )
{
    return (int)se_get_defnamedpar_int32( parset, parname, (int32_t)def, strict );
}

int64_t se_get_defnamedpar_int64( const char* parset,
                                   const char* parname,
                                   int64_t def,
                                   int strict );

char* se_get_namedhtpar_str( const char* parset,
                              const char* parname,
                              const se_hash ht );

float se_get_namedhtpar_float( const char* parset,
                                const char* parname,
                                const se_hash ht );
double se_get_namedhtpar_double( const char* parset,
                                  const char* parname,
                                  const se_hash ht );

int32_t se_get_namedhtpar_int32( const char* parset,
                                  const char* parname,
                                  const se_hash ht );
inline int se_get_namedhtpar_int( const char* parset,
                            const char* parname,
                            const se_hash ht )
{
    return (int)se_get_namedhtpar_int32( parset, parname, ht );
}

int64_t se_get_namedhtpar_int64( const char* parset,
                                  const char* parname,
                                  const se_hash ht );

int32_t se_get_namedse_par_int32( const char* parset, const char* parname );
inline int se_get_namedse_par_int( const char* parset, const char* parname )
{
    return (int)se_get_namedse_par_int32(parset, parname );
}
int64_t se_get_namedse_par_int64( const char* parset, const char* parname );
float se_get_namedse_par_float( const char* parset, const char* parname );
double se_get_namedse_par_double( const char* parset, const char* parname );
char* se_get_namedse_par_str( const char* parset, const char* parname );

int32_t se_get_namedpar_int32( const char* parset, const char* parname );
inline int se_get_namedpar_int( const char* parset, const char* parname )
{
    return (int)se_get_namedpar_int32( parset, parname );
}
int64_t se_get_namedpar_int64( const char* parset, const char* parname );
float se_get_namedpar_float( const char* parset, const char* parname );
double se_get_namedpar_double( const char* parset, const char* parname );
char* se_get_namedpar_str( const char* parset, const char* parname );

int32_t se_get_defhtpar_int32( const char* parname,
                                const se_hash ht,
                                int32_t def );
inline int se_get_defhtpar_int( const char* parname,
                                     const se_hash ht,
                                     int def )
{
    return (int)se_get_defhtpar_int32( parname, ht, (int32_t) def );
}
int64_t se_get_defhtpar_int64( const char* parname,
                                const se_hash ht,
                                int64_t def );
float se_get_defhtpar_float( const char* parname,
                              const se_hash ht,
                              float def );
double se_get_defhtpar_double( const char* parname,
                                const se_hash ht,
                                double def );
real se_get_defhtpar_real( const char* parname,
                            const se_hash ht,
                            real def );
char* se_get_defhtpar_str( const char* parname,
                            const se_hash ht,
                            const char* def );

int32_t se_get_defse_par_int32(const char* parname, int32_t def);
inline int se_get_defse_par_int(const char* parname, int def)
{
    return (int)se_get_defse_par_int32(parname, (int32_t)def);
}

int64_t se_get_defse_par_int64( const char* parname, int64_t def );
float se_get_defse_par_float( const char* parname, float def );
double se_get_defse_par_double( const char* parname, double def );
char* se_get_defse_par_str( const char* parname, const char* def );

int32_t se_get_defpar_int32( const char* parname, int32_t def );
inline int se_get_defpar_int( const char* parname, int def )
{
    return (int)se_get_defpar_int32(parname, (int32_t)def );
}

int64_t se_get_defpar_int64( const char* parname, int64_t def );
float se_get_defpar_float( const char* parname, float def );
double se_get_defpar_double( const char* parname, double def );
real se_get_defpar_real( const char* parname, real def );
char* se_get_defpar_str( const char* parname, const char* def );

char** se_get_htpararray_str( const char* parname, se_hash ht, int* size );
double* se_get_htpararray_double( const char* parname, se_hash ht,
                                   int* size );
real* se_get_htpararray_real( const char* parname, se_hash ht,
                               int* size );
float* se_get_htpararray_float( const char* parname, se_hash ht, int* size );
int64_t* se_get_htpararray_int64( const char* parname, se_hash ht,
                                   int* size );
int32_t* se_get_htpararray_int32( const char* parname, se_hash ht,
                                   int* size );
inline int32_t* se_get_htpararray_int( const char* parname, se_hash ht,
                                              int* size )
{
    return se_get_htpararray_int32( parname, ht, size );
}

int32_t* se_get_se_pararray_int32( const char* parname, int* size );
inline int32_t* se_get_se_pararray_int( const char* parname, int* size )
{
    return se_get_se_pararray_int32( parname, size );
}

int64_t* se_get_se_pararray_int64( const char* parname, int* size );
float* se_get_se_pararray_float( const char* parname, int* size );
double* se_get_se_pararray_double( const char* parname, int* size );
char** se_get_se_pararray_str( const char* parname, int* size );

int32_t* se_get_pararray_int32( const char* parname, int* size );
inline int32_t* se_get_pararray_int( const char* parname, int* size )
{
    return se_get_pararray_int32(parname, size );
}

int64_t* se_get_pararray_int64( const char* parname, int* size );
float* se_get_pararray_float( const char* parname, int* size );
double* se_get_pararray_double( const char* parname, int* size );
real* se_get_pararray_real( const char* parname, int* size );
char** se_get_pararray_str( const char* parname, int* size );

char* se_get_htpar_str( const char* parname, se_hash ht );
double se_get_htpar_double( const char* parname, se_hash ht );
real se_get_htpar_real( const char* parname, se_hash ht );
float se_get_htpar_float( const char* parname, const se_hash ht );
reala se_get_htpar_reala( const char* parname, se_hash ht );
int32_t se_get_htpar_int32( const char* parname, const se_hash ht );
inline int se_get_htpar_int( const char* parname, const se_hash ht )
{
    return (int)se_get_htpar_int32(parname, ht );
}
int64_t se_get_htpar_int64( const char* parname, se_hash ht );

int32_t se_get_separ_int32(const char* parname);
inline int se_get_separ_int(const char* parname)
{
    return (int)se_get_separ_int32(parname);
}

int64_t se_get_separ_int64( const char* parname );
float se_get_separ_float( const char* parname );
double se_get_separ_double( const char* parname );
char* se_get_separ_str( const char* parname );

int32_t se_get_par_int32(const char* parname);
inline int se_get_par_int(const char* parname)
{
    return (int)se_get_par_int32(parname);
}

int64_t se_get_par_int64( const char* parname );
float se_get_par_float( const char* parname );
double se_get_par_double( const char* parname );
real se_get_par_real( const char* parname );
reala se_get_par_reala( const char* parname );
char* se_get_par_str( const char* parname );

void se_write_pars(FILE* f, int write_regular_pars, int write_system_pars);
void se_write_htpars(FILE* f, const se_hash htpars, const char* prefix);

void se_dbg_dump_params(se_hash ht_params);
void se_dbg_dump_app_params(void);
void se_dbg_dump_sys_params(void);


int se_compare_user_params( se_hash pars );
int se_compare_system_params( se_hash pars );
int se_compare_ht_params( se_hash ht_params_actual,
                           se_hash ht_params_saved );



#endif