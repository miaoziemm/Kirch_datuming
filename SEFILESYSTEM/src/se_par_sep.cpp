#include "../include/se_par_sep.h"


se_hash _ht_params;
se_hash _ht_sys_params;

#define PAR_SRCTYPE_NOTSET  0
#define PAR_SRCTYPE_CMDLINE 1
#define PAR_SRCTYPE_FILE    2
#define PAR_SRCTYPE_STRING  3
#define PAR_SRCTYPE_SET     4

#define _PARSE_PARS_BYTE_LIMIT          (100 * 1024 * 1024)  /* 100MB */
#define _PARSE_PARS_NONPRINT_CHAR_LIMIT 1
#define _PARSE_PARS_INCLUDE_MAX_DEPTH   16
#define _PARSE_PARS_LINE_MAXLEN         (24 * 1024)          /* 24KB */

static const char _se_PARS_TO_STR_USR_NAME[] = "se_user_pars";
static const char _se_PARS_TO_STR_SYS_NAME[] = "se_sys_pars";
static const char _PARS_TO_STR_ENCODE_CHARS[] = "{}=,@ ";

static int app_argn = 0;
static char** app_args = NULL;

static int _se_par_lib_initialized = 0;

int se_par_initialized(void)
{
    return _se_par_lib_initialized;
}

static void init_random(void)
{
    /*
      better srandom() initialization then srandom(time(NULL));
      needed for se_generate_uuid()
    */
    unsigned int seed, sec, usec;
    struct timeval now;

    VERIFY(! gettimeofday( &now, NULL ));
    sec  = (unsigned int) now.tv_sec;
    usec = (unsigned int) now.tv_usec;

    seed = ((sec << 16) ^ (usec << 16)) ^ (sec ^ usec);
    srandom( seed );
}

static void _se_runtime_update_from_pars(void)
{
    if( se_have_ztrpar( "logfile" ) ) {
        char* logfile = se_get_separ_str("logfile");
        se_set_log_file( logfile );
        free(logfile);
    }

    if( se_have_ztrpar( "logid" ) ) {
        se_set_logid(se_get_separ_int32("logid"));
    } else {
        se_set_logid((int)getpid());
    }

    if( se_have_ztrpar( "verb" ) ) {
        se_set_verb_level(se_get_separ_int32("verb"));
    } else if( se_have_ztrpar( "verbose" ) ) {
        se_set_verb_level(se_get_separ_int32("verbose"));
    } else if( se_have_par( "verb" ) ) {
        se_set_verb_level(se_get_par_int32("verb"));
    } else if( se_have_par( "verbose" ) ) {
        se_set_verb_level(se_get_par_int32("verbose"));
    }

    if( se_have_ztrpar("date_format") ) {
        char* df = se_get_separ_str("date_format");
        se_set_log_date_format( df );
        free(df);
    }
}

static par_t* par_create( const char* host )
{
    par_t* p;

    p = (par_t*)malloc( sizeof(par_t) );
    p->original = se_strdup("");
    p->val_array = create_array_args(1);
    array_add( p->val_array, se_strdup("") );
    p->host = se_strdup(host);
    p->src_type = PAR_SRCTYPE_NOTSET;
    p->src_info = se_strdup("");
    p->update_count = 0;

    return p;
}

static void par_destroy(par_t* p)
{
    ASSERT(p);

    if (p->original) free1char(p->original);
    destroy_array( p->val_array, 1 );
    if (p->host) free1char(p->host);
    if (p->src_info) free1char(p->src_info);
    free(p);
}

static parinfo_t* parinfo_create(void)
{
    parinfo_t* pi;

    pi = (parinfo_t*)malloc( sizeof(parinfo_t) );
    pi->par_array = create_array_args(1);
    return pi;
}

static void parinfo_destroy(parinfo_t* pi)
{
    ASSERT(pi);

    array_process( pi->par_array, (process_value_fn) par_destroy );
    destroy_array( pi->par_array, 0 );
    free(pi);
}

/* Splits the specified param value in comma-separated values */
static array parval_split( const char* parval );

/* Parses the 'parval' and updates all fields of the specified par_t */
static void par_update( par_t* p, const char* parval,
                             int src_type, const char* src_info );

/* Create a deep clone of the specified par_t */
static par_t* par_clone(const par_t* p);

static void par_assign_str( se_hash ht,
                                 const char* parnamehost,
                                 const char* parval,
                                 int src_type,
                                 const char* src_info );

static array parinfo_lookup( const parinfo_t* pi );

static array par_lookup( se_hash ht, const char* parname );

static se_hash se_parse_pars_file_internal( const char* file_name,
                                                   se_hash ht,
                                                   const char* prefix,
                                                   se_stack_t* file_stack );

static se_hash se_parse_pars_single(const char* par,
                                           se_hash ht,
                                           const char* prefix,
                                           int src_type,
                                           const char* src_info );

static void se_make_par_subst(se_hash ht);

static char* _get_namedpar( const char* parset, const char* parname );

static se_hash _se_parse_pars_str_internal( const char* str,
                                                   se_hash ht,
                                                   const char* prefix,
                                                   se_stack_t* file_stack,
                                                   int src_type,
                                                   const char* src_info );
static void _pars_to_str_internal( se_hash ht,
                                        const char* name,
                                        char** str,
                                        int* len );

static void par_assign_array( se_hash ht,
                                   const char* parnamehost,
                                   array val_array,
                                   int src_type,
                                   const char* src_info );

static void _par_parse_keyval( char* param, se_hash ht_inner );
static se_hash _ht_from_str( const char* str );

static void _str_to_pars_internal( se_hash ht_in,
                                        se_hash* ht_out,
                                        int erase_pars,
                                        int overwrite );

int se_par_init_check_parfile(int argn, char** argv, int check_parfile)
{
    char *par, *filename, *cmdline_index;
    int i;
    const char* se_force_exe;
    int haveparfile = 0;
    int   exe_from_path = 0;
    char* exe_name = NULL;
    char* exe_path = NULL;
    char* exe_absolute_path = NULL;
    int show_info = 0, show_detailed_info = 0, show_version = 0;
    int have_prlaction = 0;


    app_argn = argn;
    app_args = (char**)malloc(app_argn * sizeof(char*));
    memset(app_args, 0, app_argn * sizeof(char*));
    for(i = 0; i < argn; ++i) {
        app_args[i] = se_strdup(argv[i]);
    }

    exe_from_path = 0;
    se_force_exe = getenv("SE_FORCE_EXE");
    if(se_force_exe != NULL) {
        if(!se_is_path_absolute(se_force_exe) && *se_force_exe != '.') {
            exe_from_path = 1;
        }

        exe_absolute_path = se_abs_path_from_peer(se_force_exe, ".");
        exe_path = se_get_base_path(exe_absolute_path);
        exe_name = se_get_filename(exe_absolute_path);
    } else {
        /* if path is not absolute, there is no . or .. specification,
           then we consider that the executable is found from PATH */
        if(!se_is_path_absolute(argv[0]) && argv[0][0] != '.') {
            exe_from_path = 1;
        }

        char* exe = se_exe_path(argv[0]);
        if(se_is_path_absolute(exe)) {
            exe_absolute_path = exe;
        } else {
            /* normaly it should be absolute, but if not ... I don't
               know, first try to solve it relative to the current
               directory; later on maybe we can parse PATH ??*/
            exe_absolute_path = se_abs_path_from_peer( exe, "./somename" );
            free1char(exe);
        }
        exe_path = se_get_base_path(exe_absolute_path);
        exe_name = se_get_filename(exe_absolute_path);
    }

    se_set_exe(exe_from_path, exe_name, exe_path, exe_absolute_path);

    init_random();

    if(_ht_params)
        se_destroy_parameters(_ht_params);
    if(_ht_sys_params)
        se_destroy_parameters(_ht_sys_params); //设置_ht_params _ht_sys_params为空

    /* create the params hashes */
    _ht_params = create_hash();
    _ht_sys_params = create_hash();

    /* parse each command-line param */ //分析每一个命令行参数
    for(i = 1; i < argn; ++i) {
        par = argv[i];
        if(!par)
            continue;

        if( str_startswith(par, "par=") ) {
            haveparfile = 1;
            filename = se_strdup(par + strlen("par="));
            se_parse_pars_file( filename, _ht_params, NULL );
            se_parse_pars_file( filename, _ht_sys_params, "ztr" );
            free1char(filename);
        } else if( str_startswith(par, "-ztr:par=") ) {
            haveparfile = 1;
            filename = se_strdup(par + strlen("-ztr:par="));
            se_parse_pars_file( filename, _ht_sys_params, NULL );
            free1char(filename);
        } else if( strcmp(par, "--show-info") == 0 ) {
            show_info = 1;
        } else if( strcmp(par, "--show-version") == 0 ) {
            show_version = 1;
        } else if( strcmp(par, "--show-detailed-info") == 0 ) {
            show_detailed_info = 1;
        } else {
            cmdline_index = int32_to_str(i);
            se_parse_pars_single( par, _ht_params, NULL,
                                    PAR_SRCTYPE_CMDLINE, cmdline_index );
            se_parse_pars_single( par, _ht_sys_params, "ztr",
                                   PAR_SRCTYPE_CMDLINE, cmdline_index );
            free1char(cmdline_index);
        }

        if( str_startswith(par, "-ztr:prlaction=") ) { /*ugly hack to check if started by libZTRDC*/
            have_prlaction = 1;
        }
    }

    se_make_par_subst( _ht_params );
    se_make_par_subst( _ht_sys_params );

    _se_runtime_update_from_pars();


    if(show_detailed_info) {
        INFO(("Application path          :[%s]", se_get_exe_path()));
        INFO(("Application name          :[%s]", se_get_exe_name()));
        INFO(("Application absolute path :[%s]", se_get_exe_absolute_path()));
        INFO(("Started from PATH         :[%s]", se_using_exe_from_path()?"true":"false"));
    }

    if(show_info || show_detailed_info || show_version) {
        exit(EXIT_SUCCESS);
    }

    if( !haveparfile && check_parfile && !have_prlaction)
        WARN(("No parameter file was entered. "
              "Please verify that all necessary parameters have been provided in command line."));

    if(se_get_verb_level() >= 1000) {
        INFO(("Using executable path         :[%s]", se_get_exe_path()));
        INFO(("Using executable name         :[%s]", se_get_exe_name()));
        INFO(("Using executable absolute path:[%s]", se_get_exe_absolute_path()));
        INFO(("Using from path               :[%s]", se_using_exe_from_path()?"true":"false"));
    }

    _se_par_lib_initialized = 1;

    if(exe_name) {
        free(exe_name);
        exe_name = NULL;
    }

    if(exe_path) {
        free(exe_path);
        exe_path = NULL;
    }

    if(exe_absolute_path) {
        free(exe_absolute_path);
        exe_absolute_path = NULL;
    }


    INFO(("%s %Ld", se_get_exe_name(),(long long int)getpid()));

    return 0;
}

se_hash se_parse_pars_file( const char* file_name,
                              se_hash ht, const char* prefix )
{
    return se_parse_pars_file_internal( file_name, ht, prefix, NULL );
}

void se_par_destroy_internal(void)
{
    se_destroy_parameters(_ht_params);
    se_destroy_parameters(_ht_sys_params);
    _ht_params = _ht_sys_params = NULL;

    se_exe_reset();

    if(app_argn > 0) {
        int i;
        for(i = 0; i < app_argn; ++i) {
            free(app_args[i]);
        }
        free(app_args);
    }

    _se_par_lib_initialized = 0;
}

void se_destroy_parameters( se_hash ht )
{
    htiter hti;
    const char* parname;
    void* ht_val;
    parinfo_t* pi;
    ASSERT(ht);

    /* destroy all parinfo structures */
    hti = create_htiter( ht );
    while( hti_hasnext( hti )) {
        hti_next( hti, &parname, &ht_val );
        pi = (parinfo_t*) ht_val;

        if(parname) {
            free1char((void*)parname);
        }
        parinfo_destroy( pi );
    }
    destroy_htiter(hti);

    /* destroy the hash table */
    destroy_hash_and_entries( ht, 0 );
}


void se_dbg_dump_app_params(void)
{
    se_dbg_dump_params(_ht_params);
}

void se_dbg_dump_sys_params(void)
{
    se_dbg_dump_params(_ht_sys_params);
}

void se_dbg_dump_params(se_hash ht_params)
{
    htiter hti;
    const char* key;
    void* ht_val;
    parinfo_t* pi;
    par_t* p;
    int i, j;

    hti = create_htiter( ht_params );
    while( hti_hasnext( hti )) {
        hti_next( hti, &key, &ht_val );
        pi = (parinfo_t*) ht_val;

        for (i = 0; i < array_size( pi->par_array ); ++i) {
            p = (par_t*) array_get_at(pi->par_array, i);

            if( array_size( p->val_array ) == 1 ) {
                INFO(("%s%s%s='%s'", key, strlen(p->host) ? "@" : "", p->host,
                      array_get_at( p->val_array, 0 )));
            } else {
                for (j = 0; j < array_size( p->val_array ); ++j) {
                    INFO(("%s%s%s='%s' [#%d]", key,
                          strlen(p->host) ? "@" : "", p->host,
                          array_get_at( p->val_array, j ), j));
                }
            }
        }
    }
    destroy_htiter(hti);
}


int32_t se_get_par_int32(const char* parname)
{
    return (int32_t) se_get_htpar_double( parname, _ht_params );
}

int64_t se_get_par_int64( const char* parname )
{
    return se_get_htpar_int64( parname, _ht_params );
}

float se_get_par_float( const char* parname )
{
    return (float) se_get_htpar_double( parname, _ht_params );
}

double se_get_par_double( const char* parname )
{
    return se_get_htpar_double( parname, _ht_params );
}
real se_get_par_real( const char* parname )
{
    return (real)se_get_htpar_double( parname, _ht_params );
}
reala se_get_par_reala( const char* parname )
{
    return (reala)se_get_htpar_double( parname, _ht_params );
}

char* se_get_par_str( const char* parname )
{
    return se_get_htpar_str( parname, _ht_params );
}

/* get par [ztr] */


int32_t se_get_separ_int32(const char* parname)
{
    return (int32_t) se_get_htpar_double( parname, _ht_sys_params );
}

int64_t se_get_separ_int64( const char* parname )
{
    return se_get_htpar_int64( parname, _ht_sys_params );
}

float se_get_separ_float( const char* parname )
{
    return (float) se_get_htpar_double( parname, _ht_sys_params );
}


char* se_get_separ_str( const char* parname )
{
    return se_get_htpar_str( parname, _ht_sys_params );
}



double se_get_separ_double( const char* parname )
{
    return se_get_htpar_double( parname, _ht_sys_params );
}


/* get par [hash] */

int32_t se_get_htpar_int32( const char* parname, const se_hash ht )
{
    return (int32_t) se_get_htpar_double( parname, ht );
}

int64_t se_get_htpar_int64( const char* parname, se_hash ht )
{
    array val_array;
    const char* val_str;
    int64_t val_int64;

    val_array = par_lookup( ht, parname );
    if(! val_array) {
        return 0;
    }
    val_str = (const char*) array_get_first( val_array );
    val_int64 = (int64_t) atoll(val_str);

    return val_int64;
}

float se_get_htpar_float( const char* parname, const se_hash ht )
{
    return (float) se_get_htpar_double( parname, ht );
}
real se_get_htpar_real( const char* parname, const se_hash ht )
{
    return (real)se_get_htpar_double( parname, ht );
}
reala se_get_htpar_reala( const char* parname, const se_hash ht )
{
    return (reala)se_get_htpar_double( parname, ht );
}

double se_get_htpar_double( const char* parname, se_hash ht )
{
    array val_array;
    const char* val_str;
    double val_d;

    val_array = par_lookup( ht, parname );
    if(! val_array) {
        return 0.0;
    }
    val_str = (const char*) array_get_first( val_array );
    val_d = atof(val_str);

    return val_d;
}

char* se_get_htpar_str( const char* parname, se_hash ht )
{
    array val_array;
    const char* val_str;

    val_array = par_lookup( ht, parname );
    if(! val_array) {
        return se_strdup("");
    }
    val_str = (const char*) array_get_first( val_array );
    return se_strdup(val_str);
}


/* get par array [app] */

int32_t* se_get_pararray_int32( const char* parname, int* size )
{
    return se_get_htpararray_int32( parname, _ht_params, size );
}

int64_t* se_get_pararray_int64( const char* parname, int* size )
{
    return se_get_htpararray_int64( parname, _ht_params, size);
}

float* se_get_pararray_float( const char* parname, int* size )
{
    return se_get_htpararray_float( parname, _ht_params, size);
}

double* se_get_pararray_double( const char* parname, int* size )
{
    return se_get_htpararray_double( parname, _ht_params, size );
}

real* se_get_pararray_real( const char* parname, int* size )
{
    return se_get_htpararray_real( parname, _ht_params, size );
}

char** se_get_pararray_str( const char* parname, int* size )
{
    return se_get_htpararray_str( parname, _ht_params, size );
}

/* get par array [ztr] */

int32_t* se_get_ztrpararray_int32( const char* parname, int* size )
{
    return se_get_htpararray_int32( parname, _ht_sys_params, size );
}

int64_t* se_get_ztrpararray_int64( const char* parname, int* size )
{
    return se_get_htpararray_int64( parname, _ht_sys_params, size);
}

float* se_get_ztrpararray_float( const char* parname, int* size )
{
    return se_get_htpararray_float( parname, _ht_sys_params, size);
}

double* se_get_ztrpararray_double( const char* parname, int* size )
{
    return se_get_htpararray_double( parname, _ht_sys_params, size );
}

char** se_get_ztrpararray_str( const char* parname, int* size )
{
    return se_get_htpararray_str( parname, _ht_sys_params, size );
}

/* get par array [hash] */

int32_t* se_get_htpararray_int32( const char* parname, se_hash ht, int* size )
{
    array val_array;
    const char* val_str;
    int32_t* values;
    int i;
    ASSERT(size);

    val_array = par_lookup( ht, parname );
    if(! val_array) {
        *size = 0;
        return NULL;
    }
    *size = array_size( val_array );
    values = (int32_t*)malloc( *size * sizeof(*values) );
    for (i = 0; i < *size; ++i) {
        val_str = (const char*) array_get_at( val_array, i );
        values[i] = (int32_t) atof(val_str);
    }
    return values;
}

int64_t* se_get_htpararray_int64( const char* parname, se_hash ht, int* size )
{
    array val_array;
    const char* val_str;
    int64_t* values;
    int i;
    ASSERT(size);

    val_array = par_lookup( ht, parname );
    if(! val_array) {
        *size = 0;
        return NULL;
    }
    *size = array_size( val_array );
    values = (int64_t*)malloc( *size * sizeof(*values) );
    for (i = 0; i < *size; ++i) {
        val_str = (const char*) array_get_at( val_array, i );
        values[i] = (int64_t) atoll(val_str);
    }
    return values;
}

float* se_get_htpararray_float( const char* parname, se_hash ht, int* size )
{
    array val_array;
    const char* val_str;
    float* values;
    int i;
    ASSERT(size);

    val_array = par_lookup( ht, parname );
    if(! val_array) {
        *size = 0;
        return NULL;
    }
    *size = array_size( val_array );
    values = (float*)malloc( *size * sizeof(*values) );
    for (i = 0; i < *size; ++i) {
        val_str = (const char*) array_get_at( val_array, i );
        values[i] = (float) atof(val_str);
    }
    return values;
}

double* se_get_htpararray_double( const char* parname, se_hash ht, int* size )
{
    array val_array;
    const char* val_str;
    double* values;
    int i;

    ASSERT(size);

    val_array = par_lookup( ht, parname );
    if(! val_array) {
        *size = 0;
        return NULL;
    }
    *size = array_size( val_array );
    values = (double*)malloc( *size * sizeof(*values) );
    for (i = 0; i < *size; ++i) {
        val_str = (const char*) array_get_at( val_array, i );
        values[i] = atof(val_str);
    }
    return values;
}

real* se_get_htpararray_real( const char* parname, se_hash ht, int* size )
{
    array val_array;
    const char* val_str;
    real* values;
    int i;

    ASSERT(size);

    val_array = par_lookup( ht, parname );
    if(! val_array) {
        *size = 0;
        return NULL;
    }
    *size = array_size( val_array );
    values = (real*)malloc( *size * sizeof(*values) );
    for (i = 0; i < *size; ++i) {
        val_str = (const char*) array_get_at( val_array, i );
        values[i] = (real)atof(val_str);
    }
    return values;
}

char** se_get_htpararray_str( const char* parname, se_hash ht, int* size )
{
    array val_array;
    const char* val_str;
    char** values;
    int i;

    ASSERT(size);

    val_array = par_lookup( ht, parname );
    if(! val_array) {
        *size = 0;
        return NULL;
    }
    *size = array_size( val_array );
    values = (char**)malloc( *size * sizeof(*values) );
    for (i = 0; i < *size; ++i) {
        val_str = (const char*) array_get_at( val_array, i );
        values[i] = se_strdup(val_str);
    }
    return values;
}


/* get default par (app) */

int32_t se_get_defpar_int32( const char* parname, int32_t def )
{
    return (int32_t) se_get_defhtpar_int64( parname, _ht_params,
                                             (int64_t) def );
}

int64_t se_get_defpar_int64( const char* parname, int64_t def )
{
    return se_get_defhtpar_int64( parname, _ht_params, def );
}

float se_get_defpar_float( const char* parname, float def )
{
    return (float) se_get_defhtpar_double( parname, _ht_params, (double) def );
}

double se_get_defpar_double( const char* parname, double def )
{
    return se_get_defhtpar_double( parname, _ht_params, def );
}

real se_get_defpar_real( const char* parname, real def )
{
    return se_get_defhtpar_real( parname, _ht_params, def );
}

char* se_get_defpar_str( const char* parname, const char* def )
{
    return se_get_defhtpar_str( parname, _ht_params, def );
}

/* get default par (ztr) */

int32_t se_get_defsepar_int32(const char* parname, int32_t def)
{
    return (int32_t) se_get_defhtpar_int64( parname, _ht_sys_params,
                                             (int64_t) def );
}

int64_t se_get_defsepar_int64( const char* parname, int64_t def )
{
    return se_get_defhtpar_int64( parname, _ht_sys_params, def );
}

float se_get_defsepar_float( const char* parname, float def )
{
    return (float) se_get_defhtpar_double( parname, _ht_sys_params,
                                            (double) def );
}

double se_get_defsepar_double( const char* parname, double def )
{
    return se_get_defhtpar_double( parname, _ht_sys_params, def );
}

char* se_get_defsepar_str( const char* parname, const char* def )
{
    return se_get_defhtpar_str( parname, _ht_sys_params, def );
}

/* get default par (hash) */

int32_t se_get_defhtpar_int32( const char* parname,
                                const se_hash ht,
                                int32_t def )
{
    return (int32_t) se_get_defhtpar_int64( parname, ht, (int64_t) def );
}

int64_t se_get_defhtpar_int64( const char* parname,
                                const se_hash ht,
                                int64_t def )
{
    if(! se_have_htpar( parname, ht )) {
        return def;
    }
    return se_get_htpar_int64( parname, ht );
}

float se_get_defhtpar_float( const char* parname,
                              const se_hash ht,
                              float def )
{
    return (float) se_get_defhtpar_double( parname, ht, (double) def );
}

double se_get_defhtpar_double( const char* parname,
                                const se_hash ht,
                                double def )
{
    if(! se_have_htpar( parname, ht )) {
        return def;
    }
    return (double) se_get_htpar_double( parname, ht );
}

real se_get_defhtpar_real( const char* parname,
                            const se_hash ht,
                            real def )
{
    if(! se_have_htpar( parname, ht )) {
        return def;
    }
    return se_get_htpar_real( parname, ht );
}

char* se_get_defhtpar_str( const char* parname,
                            const se_hash ht,
                            const char* def )
{
    if(! se_have_htpar( parname, ht ) ) {
        return se_strdup(def);
    }
    return se_get_htpar_str( parname, ht );
}


/* get named par (app) */

int32_t se_get_namedpar_int32( const char* parset, const char* parname )
{
    return se_get_namedhtpar_int32( parset, parname, _ht_params );
}

int64_t se_get_namedpar_int64( const char* parset, const char* parname )
{
    return se_get_namedhtpar_int64( parset, parname, _ht_params );
}

float se_get_namedpar_float( const char* parset, const char* parname )
{
    return se_get_namedhtpar_float( parset, parname, _ht_params );
}

double se_get_namedpar_double( const char* parset, const char* parname )
{
    return se_get_namedhtpar_double( parset, parname, _ht_params );
}

char* se_get_namedpar_str( const char* parset, const char* parname )
{
    return se_get_namedhtpar_str( parset, parname, _ht_params );
}

/* get named par (ztr) */

int32_t se_get_namedsepar_int32( const char* parset, const char* parname )
{
    return (int32_t) se_get_namedhtpar_int64( parset, parname,
                                               _ht_sys_params );
}

int64_t se_get_namedsepar_int64( const char* parset, const char* parname )
{
    return se_get_namedhtpar_int64( parset, parname, _ht_sys_params );
}

float se_get_namedsepar_float( const char* parset, const char* parname )
{
    return (float) se_get_namedhtpar_double( parset, parname, _ht_sys_params );
}

double se_get_namedsepar_double( const char* parset, const char* parname )
{
    return se_get_namedhtpar_double( parset, parname, _ht_sys_params );
}

char* se_get_namedsepar_str( const char* parset, const char* parname )
{
    return se_get_namedhtpar_str( parset, parname, _ht_sys_params );
}


int32_t se_get_namedhtpar_int32( const char* parset,
                                  const char* parname,
                                  const se_hash ht )
{
    return (int32_t) se_get_namedhtpar_int64( parset, parname, ht );
}

int64_t se_get_namedhtpar_int64( const char* parset,
                                  const char* parname,
                                  const se_hash ht )
{
    char* par;
    int64_t val;

    par = _get_namedpar( parset, parname );
    if( se_have_htpar( par, ht )) {
        val = se_get_par_int64( par );
    } else {
        val = se_get_par_int64( parname );
    }
    free1char( par );

    return val;
}

float se_get_namedhtpar_float( const char* parset,
                                const char* parname,
                                const se_hash ht )
{
    return (float) se_get_namedhtpar_double( parset, parname, ht );
}

double se_get_namedhtpar_double( const char* parset,
                                  const char* parname,
                                  const se_hash ht )
{
    char* par;
    double val;

    par = _get_namedpar( parset, parname );
    if( se_have_htpar( par, ht ) ) {
        val = se_get_par_double( par );
    } else {
        val = se_get_par_double( parname );
    }
    free1char( par );

    return val;
}

char* se_get_namedhtpar_str( const char* parset,
                              const char* parname,
                              const se_hash ht )
{
    char* par;
    char* val;

    par = _get_namedpar( parset, parname );
    if( se_have_htpar( par, ht )) {
        val = se_get_par_str( par );
    } else {
        val = se_get_par_str( parname );
    }
    free1char( par );

    return val;
}



int32_t se_get_defnamedpar_int32( const char* parset,
                                   const char* parname,
                                   int32_t def,
                                   int strict)
{
    return (int32_t) se_get_defnamedhtpar_int64( parset, parname, _ht_params,
                                                  (int64_t) def, strict );
}

int64_t se_get_defnamedpar_int64( const char* parset,
                                   const char* parname,
                                   int64_t def,
                                   int strict)
{
    return se_get_defnamedhtpar_int64( parset, parname, _ht_params, def, strict );
}

float se_get_defnamedpar_float( const char* parset,
                                 const char* parname,
                                 float def,
                                 int strict )
{
    return (float) se_get_defnamedhtpar_double( parset, parname, _ht_params,
                                                 (double) def, strict );
}

double se_get_defnamedpar_double( const char* parset,
                                   const char* parname,
                                   double def,
                                   int strict )
{
    return se_get_defnamedhtpar_double( parset, parname, _ht_params, def, strict );
}

char* se_get_defnamedpar_str( const char* parset,
                               const char* parname,
                               const char* def,
                               int strict )
{
    return se_get_defnamedhtpar_str( parset, parname, _ht_params, def, strict );
}

/* get default named par (ztr) */

int32_t se_get_defnamedsepar_int32( const char* parset,
                                      const char* parname,
                                      int32_t def,
                                      int strict )
{
    return (int32_t) se_get_defnamedhtpar_int64( parset, parname,
                                                  _ht_sys_params,
                                                  (int64_t) def, strict );
}

int64_t se_get_defnamedsepar_int64( const char* parset,
                                      const char* parname,
                                      int64_t def,
                                      int strict)
{
    return se_get_defnamedhtpar_int64( parset, parname, _ht_sys_params, def, strict );
}

float se_get_defnamedsepar_float( const char* parset,
                                    const char* parname,
                                    float def,
                                    int strict )
{
    return (float) se_get_defnamedhtpar_double( parset, parname,
                                                 _ht_sys_params, (double) def, strict );
}

double se_get_defnamedsepar_double( const char* parset,
                                      const char* parname,
                                      double def,
                                      int strict )
{
    return se_get_defnamedhtpar_double( parset, parname, _ht_sys_params, def, strict );
}

char* se_get_defnamedsepar_str( const char* parset,
                                  const char* parname,
                                  const char* def,
                                  int strict )
{
    return se_get_defnamedhtpar_str( parset, parname, _ht_sys_params, def, strict );
}


int32_t se_get_defnamedhtpar_int32( const char* parset,
                                     const char* parname,
                                     const se_hash ht,
                                     int32_t def,
                                     int strict )
{
    return (int32_t) se_get_defnamedhtpar_int64( parset, parname, ht,
                                                  (int64_t) def, strict );
}

int64_t se_get_defnamedhtpar_int64( const char* parset,
                                     const char* parname,
                                     const se_hash ht,
                                     int64_t def,
                                     int strict )
{
    if(! se_have_namedhtpar( parset, parname, ht, strict )) {
        return def;
    }
    return se_get_namedhtpar_int64( parset, parname, ht );
}

float se_get_defnamedhtpar_float( const char* parset,
                                   const char* parname,
                                   const se_hash ht,
                                   float def,
                                   int strict )
{
    return (float) se_get_defnamedhtpar_double( parset,
                                                 parname, ht,
                                                 (double) def, 
                                                 strict );
}

double se_get_defnamedhtpar_double( const char* parset,
                                     const char* parname,
                                     const se_hash ht,
                                     double def,
                                     int strict )
{
    if(! se_have_namedhtpar( parset, parname, ht, strict )) {
        return def;
    }
    return se_get_namedhtpar_double( parset, parname, ht );
}

char* se_get_defnamedhtpar_str( const char* parset,
                                 const char* parname,
                                 const se_hash ht,
                                 const char* def,
                                 int strict )
{
    if(! se_have_namedhtpar( parset, parname, ht, strict )) {
        return se_strdup(def);
    }
    return se_get_namedhtpar_str( parset, parname, ht );
}


/* set par (app) */

void se_set_par_int32( const char* parname, int32_t val )
{
    se_set_htpar_int64( parname, _ht_params, (int64_t) val );
}

void se_set_par_int64( const char* parname, int64_t val )
{
    se_set_htpar_int64( parname, _ht_params, val );
}

void se_set_par_float( const char* parname, float val )
{
    se_set_htpar_double( parname, _ht_params, (double) val );
}

void se_set_par_double( const char* parname, double val )
{
    se_set_htpar_double( parname, _ht_params, val);
}
void se_set_par_real( const char* parname, real val )
{
    se_set_htpar_double( parname, _ht_params, (double)val);
}

void se_set_par_reala( const char* parname, reala val )
{
    se_set_htpar_double( parname, _ht_params, (double)val);
}

void se_set_par_str( const char* parname, const char* val )
{
    se_set_htpar_str( parname, _ht_params, val);
}


void se_set_separ_int32( const char* parname, int32_t val )
{
    se_set_htpar_int64( parname, _ht_sys_params, (int64_t) val );
}

void se_set_separ_int64( const char* parname, int64_t val )
{
    se_set_htpar_int64( parname, _ht_sys_params, val );
}

void se_set_separ_float( const char* parname, float val )
{
    se_set_htpar_double( parname, _ht_sys_params, (double) val);
}

void se_set_separ_double( const char* parname, double val )
{
    se_set_htpar_double( parname, _ht_sys_params, val);
}

void se_set_separ_str( const char* parname, const char* val )
{
    se_set_htpar_str( parname, _ht_sys_params, val);
}


void se_set_htpar_int32( const char* parname, se_hash ht, int32_t val )
{
    se_set_htpar_int64( parname, ht, (int64_t) val );
}

void se_set_htpar_int64( const char* parname, se_hash ht, int64_t val )
{
    char* val_str = int64_to_str(val);
    par_assign_str( ht, parname, val_str, PAR_SRCTYPE_SET, "" );
    se_make_par_subst( ht );
    free1char(val_str);
}

void se_set_htpar_float( const char* parname, se_hash ht, float val )
{
    se_set_htpar_double( parname, ht, (double) val );
}
void se_set_htpar_real( const char* parname, se_hash ht, real val )
{
    se_set_htpar_double( parname, ht, (double) val );
}
void se_set_htpar_reala( const char* parname, se_hash ht, reala val )
{
    se_set_htpar_double( parname, ht, (double) val );
}

void se_set_htpar_double( const char* parname, se_hash ht, double val )
{
    char* val_str = double_to_str(val);
    par_assign_str( ht, parname, val_str, PAR_SRCTYPE_SET, "" );
    se_make_par_subst( ht );
    free1char(val_str);
}

void se_set_htpar_str( const char* parname, se_hash ht, const char* val )
{
    par_assign_str( ht, parname, val, PAR_SRCTYPE_SET, "" );
    se_make_par_subst( ht );
}



void se_set_pararray_int32( const char* parname,
                             const int32_t* data,
                             int size )
{
    se_set_htpararray_int32( parname, _ht_params, data, size );
}

void se_set_pararray_int64( const char* parname,
                             const int64_t* data,
                             int size )
{
    se_set_htpararray_int64( parname, _ht_params, data, size);
}

void se_set_pararray_float( const char* parname,
                             const float* data,
                             int size )
{
    se_set_htpararray_float( parname, _ht_params, data, size );
}

void se_set_pararray_double( const char* parname,
                              const double* data,
                              int size )
{
    se_set_htpararray_double( parname, _ht_params, data, size );
}

void se_set_pararray_str( const char* parname,
                           const char** data,
                           int size )
{
    se_set_htpararray_str( parname, _ht_params, data, size );
}

/* set par array (ztr) */

void se_set_ztrpararray_int32( const char* parname,
                                const int32_t* data,
                                int size )
{
    se_set_htpararray_int32( parname, _ht_sys_params, data, size );
}

void se_set_ztrpararray_int64( const char* parname,
                                const int64_t* data,
                                int size )
{
    se_set_htpararray_int64( parname, _ht_sys_params, data, size);
}

void se_set_ztrpararray_float( const char* parname,
                                const float* data,
                                int size )
{
    se_set_htpararray_float( parname, _ht_sys_params, data, size );
}

void se_set_ztrpararray_double( const char* parname,
                                 const double* data,
                                 int size )
{
    se_set_htpararray_double( parname, _ht_sys_params, data, size );
}

void se_set_ztrpararray_str( const char* parname,
                              const char** data,
                              int size )
{
    se_set_htpararray_str( parname, _ht_sys_params, data, size );
}


void se_set_htpararray_int32( const char* parname,
                               se_hash ht,
                               const int32_t* data,
                               int size )
{
    array val_array;
    char* val_str;
    int i;

    if(! data || ! size) {
        return;
    }
    val_array = create_array_args(size);
    for (i = 0; i < size; ++i) {
        val_str = int32_to_str( data[i] );
        array_add( val_array, val_str );
    }
    par_assign_array( ht, parname, val_array, PAR_SRCTYPE_SET, "" );
    destroy_array( val_array, 1 );
}

void se_set_htpararray_int64( const char* parname,
                               se_hash ht,
                               const int64_t* data,
                               int size )
{
    array val_array;
    char* val_str;
    int i;

    if(! data || ! size) {
        return;
    }
    val_array = create_array_args(size);
    for (i = 0; i < size; ++i) {
        val_str = int64_to_str( data[i] );
        array_add( val_array, val_str );
    }
    par_assign_array( ht, parname, val_array, PAR_SRCTYPE_SET, "" );
    destroy_array( val_array, 1 );
}

void se_set_htpararray_float( const char* parname,
                               se_hash ht,
                               const float* data,
                               int size )
{
    array val_array;
    char* val_str;
    int i;

    if(! data || ! size) {
        return;
    }
    val_array = create_array_args(size);
    for (i = 0; i < size; ++i) {
        val_str = double_to_str( (double) data[i] );
        array_add( val_array, val_str );
    }
    par_assign_array( ht, parname, val_array, PAR_SRCTYPE_SET, "" );
    destroy_array( val_array, 1 );
}

void se_set_htpararray_double( const char* parname,
                                se_hash ht,
                                const double* data,
                                int size )
{
    array val_array;
    char* val_str;
    int i;

    if(! data || ! size) {
        return;
    }

    val_array = create_array_args(size);
    for (i = 0; i < size; ++i) {
        val_str = double_to_str( data[i] );
        array_add( val_array, val_str );
    }

    par_assign_array( ht, parname, val_array, PAR_SRCTYPE_SET, "" );
    destroy_array( val_array, 1 );
}

void se_set_htpararray_str( const char* parname,
                             se_hash ht,
                             const char** data,
                             int size )
{
    array val_array;
    char* val_str;
    int i;

    if(! data || ! size) {
        return;
    }

    val_array = create_array_args(size);

    for (i = 0; i < size; ++i) {
        val_str = se_strdup( data[i] );
        array_add( val_array, val_str );
    }

    par_assign_array( ht, parname, val_array, PAR_SRCTYPE_SET, "" );

    destroy_array( val_array, 1 );
}



void se_set_namedpar_int32( const char* parset,
                             const char* parname,
                             int32_t val )
{
    se_set_namedhtpar_int64( parset, parname, _ht_params, (int64_t) val );
}

void se_set_namedpar_int64( const char* parset,
                             const char* parname,
                             int64_t val )
{
    se_set_namedhtpar_int64( parset, parname, _ht_params, val );
}

void se_set_namedpar_float( const char* parset,
                             const char* parname,
                             float val )
{
    se_set_namedhtpar_double( parset, parname, _ht_params, (double) val );
}

void se_set_namedpar_double( const char* parset,
                              const char* parname,
                              double val )
{
    se_set_namedhtpar_double( parset, parname, _ht_params, val );
}

void se_set_namedpar_str( const char* parset,
                           const char* parname,
                           const char* val )
{
    se_set_namedhtpar_str( parset, parname, _ht_params, val );
}

void se_set_namedsepar_int32( const char* parset,
                                const char* parname,
                                int32_t val )
{
    se_set_namedhtpar_int64( parset, parname, _ht_sys_params, (int64_t) val );
}

void se_set_namedsepar_int64( const char* parset,
                                const char* parname,
                                int64_t val )
{
    se_set_namedhtpar_int64( parset, parname, _ht_sys_params, val );
}

void se_set_namedsepar_float( const char* parset,
                                const char* parname,
                                float val )
{
    se_set_namedhtpar_double( parset, parname, _ht_sys_params, (double) val );
}

void se_set_namedsepar_double( const char* parset,
                                 const char* parname,
                                 double val )
{
    se_set_namedhtpar_double( parset, parname, _ht_sys_params, val );
}

void se_set_namedsepar_str( const char* parset,
                              const char* parname,
                              const char* val )
{
    se_set_namedhtpar_str( parset, parname, _ht_sys_params, val );
}


void se_set_namedhtpar_int32( const char* parset,
                               const char* parname,
                               se_hash ht,
                               int32_t val )
{
    char* par = _get_namedpar( parset, parname );
    se_set_htpar_int32( par, ht, val );
    free1char( par );
}

void se_set_namedhtpar_int64( const char* parset,
                               const char* parname,
                               se_hash ht,
                               int64_t val )
{
    char* par = _get_namedpar( parset, parname );
    se_set_htpar_int64( par, ht, val );
    free1char( par );
}

void se_set_namedhtpar_float( const char* parset,
                               const char* parname,
                               se_hash ht,
                               float val )
{
    char* par = _get_namedpar( parset, parname );
    se_set_htpar_float( par, ht, val );
    free1char( par );
}

void se_set_namedhtpar_double( const char* parset,
                                const char* parname,
                                se_hash ht,
                                double val )
{
    char* par = _get_namedpar( parset, parname );
    se_set_htpar_double( par, ht, val );
    free1char( par );
}

void se_set_namedhtpar_str( const char* parset,
                             const char* parname,
                             se_hash ht,
                             const char* val )
{
    char* par = _get_namedpar( parset, parname );
    se_set_htpar_str( par, ht, val );
    free1char( par );
}


/* have par */

int se_have_par( const char* parname )
{
    if (_ht_params == NULL) {
        ERROR(("Parameters are not initialized."));
    }
    return se_have_htpar( parname, _ht_params );
}

int se_have_ztrpar( const char* parname )
{
    if (_ht_sys_params == NULL) {
        ERROR(("Parameters are not initialized."));
    }
    return se_have_htpar( parname, _ht_sys_params );
}

int se_have_se_par( const char* parname )
{
     if (_ht_sys_params == NULL) {
        ERROR(("Parameters are not initialized."));
    }
    return se_have_htpar( parname, _ht_sys_params );
}

int se_have_htpar( const char* parname, se_hash ht )
{
    return par_lookup( ht, parname ) != NULL;
}


int se_have_namedpar( const char* parset, const char* parname, int strict )
{
    if (_ht_params == NULL) {
        ERROR(("Parameters are not initialized."));
    }
    return se_have_namedhtpar( parset, parname, _ht_params, strict );
}

int se_have_namedztrpar( const char* parset, const char* parname, int strict )
{
    if (_ht_sys_params == NULL) {
        ERROR(("Parameters are not initialized."));
    }
    return se_have_namedhtpar( parset, parname, _ht_sys_params, strict );
}

int se_have_namedhtpar( const char* parset,
                         const char* parname,
                         se_hash ht,
                         int strict )
{
    char* par;
    int have;

    par = _get_namedpar( parset, parname );
    have = se_have_htpar( par, ht );
    if( (!have) && (!strict) ) {
        have = se_have_htpar( parname, ht );
    }
    free1char( par );

    return have;
}


char* se_pars_to_str( int add_user_pars, int add_sys_pars )
{
    char* str = NULL;
    int len = 0;

    if (add_user_pars) {
        _pars_to_str_internal( _ht_params, _se_PARS_TO_STR_USR_NAME,
                               &str, &len );
    }
    if (add_sys_pars) {
        _pars_to_str_internal( _ht_sys_params, _se_PARS_TO_STR_SYS_NAME,
                               &str, &len );
    }

    if (str == NULL) {
        str = se_strdup("");
    }

    return str;
}

void se_str_to_pars( const char* str, int erase_pars, int overwrite )
{
    se_hash hts;
    se_hash ht;
    htiter htis;
    const char *key;

    VERIFY(str);

    hts = _ht_from_str(str);
    VERIFY(hts);

    /* parse app params */
    ht = ht_get( hts, _se_PARS_TO_STR_USR_NAME );
    if(ht) {
        _str_to_pars_internal( ht, &_ht_params, erase_pars, overwrite );
    }

    /* parse ztr params */
    ht = ht_get( hts, _se_PARS_TO_STR_SYS_NAME );
    if(ht) {
        _str_to_pars_internal( ht, &_ht_sys_params, erase_pars, overwrite );
    }

    /* cleanup */
    htis = create_htiter( hts );
    while( hti_hasnext(htis) ) {
        hti_next( htis, &key, &ht );

        free1char( (void *)key );
        se_destroy_parameters( (se_hash) ht );
    }
    destroy_htiter( htis );
    destroy_hash( hts );

    /* update ztr params dependent behaviour */
    _se_runtime_update_from_pars();
}


char* se_dump_pars( int dump_app_pars, int dump_se_pars )
{
    char *result, *tmp;
    strbuf sb;

    ASSERT( dump_app_pars || dump_se_pars );

    sb = create_strbuf();

    /* append user (app) params */
    if (dump_app_pars) {
        tmp = se_dump_htpars( _ht_params, NULL );
        sb_append( sb, tmp );
        free1char( tmp );
    }

    /* append system (ztr) params */
    if (dump_se_pars) {
        tmp = se_dump_htpars( _ht_sys_params, "ztr" );
        sb_append( sb, tmp );
        free1char( tmp );
    }

    result = sb_release( sb );
    destroy_strbuf( sb );

    return result;
}

char* se_dump_htpars( const se_hash _pars, const char* prefix )
{
    htiter hti;
    const char* key;
    void* ht_val;
    const parinfo_t* pi;
    const par_t* p;
    int i, have_prefix;
    strbuf sb;
    char* result;

    se_hash pars = _pars?_pars:_ht_params;
    ASSERT(pars);

    sb = create_strbuf_args( 2048 );
    if (prefix == NULL)
        prefix = "";
    have_prefix = (strlen(prefix) > 0);

    /* for each param */
    hti = create_htiter( pars );
    while( hti_hasnext( hti )) {
        hti_next( hti, &key, &ht_val );
        pi = (parinfo_t*) ht_val;

        /* for each host entry */
        for (i = 0; i < array_size( pi->par_array ); ++i) {
            p = (par_t*) array_get_at(pi->par_array, i);

            /* append the param: [-prefix:]name[@host]="original_val"\n */
            if (have_prefix) {
                sb_append( sb, "-" );
                sb_append( sb, prefix );
                sb_append( sb, ":" );
            }
            sb_append( sb, key );
            if (p->host[0] != '\0') {
                sb_append( sb, "@" );
                sb_append( sb, p->host );
            }
            sb_append( sb, "=" );
            sb_appendf(sb, "\"%s\"", p->original);
            sb_append( sb, "\n" );
        }
    }
    destroy_htiter( hti );

    result = sb_release( sb );

    destroy_strbuf( sb );

    return result;
}


void se_set_pars_from_strdump(const char* str, int set_user, int set_system)
{
    if(set_user) {
        se_parse_pars_str( str, _ht_params, NULL );
        se_make_par_subst( _ht_params );
    }

    if(set_system) {
        se_parse_pars_str( str, _ht_sys_params, "ztr" );
        se_make_par_subst( _ht_sys_params );
    }
    _se_runtime_update_from_pars();
}

void se_set_pars_from_file(const char* fn, int set_user, int set_system)
{
    if(set_user) {
        se_parse_pars_file( fn, _ht_params, NULL );
        se_make_par_subst( _ht_params );
    }

    if(set_system) {
        se_parse_pars_file( fn, _ht_sys_params, "ztr" );
        se_make_par_subst( _ht_sys_params );
        _se_runtime_update_from_pars();
    }
}


se_hash se_parse_pars_str(const char* str,
                            se_hash ht,
                            const char* prefix)
{
    return _se_parse_pars_str_internal( str, ht, prefix,
                                         NULL, PAR_SRCTYPE_STRING, "" );
}





static se_hash _se_parse_pars_fstream_internal( FILE* file,
                                                       se_hash ht,
                                                       const char* prefix,
                                                       se_stack_t* file_stack )
{
    char line[_PARSE_PARS_LINE_MAXLEN], *p, *lineno_str;
    uint64_t byte_count = 0, nonprint_char_count = 0;
    int limit_reached = 0, lineno = 0, line_len;

    ASSERT(file);

    /* create the hash if necessary */
    if (ht == NULL) {
        ht = create_hash();
    }

    do {
        /* read next line */
        ++lineno;
        if(fgets( line, _PARSE_PARS_LINE_MAXLEN, file ) == NULL) {
            if( ! feof(file)) {
                ERROR(("Error reading line %d from file", lineno));
            }
            break;
        }

        /* check for total bytes limit */
        line_len = (uint64_t)strlen(line);
        byte_count += line_len;
        if( byte_count > _PARSE_PARS_BYTE_LIMIT ) {
            WARN(("Aborting parameter parsing: Parsed more than %d bytes "
                  "from a single stream",
                  _PARSE_PARS_BYTE_LIMIT));
            line[line_len - (byte_count - _PARSE_PARS_BYTE_LIMIT)] = '\0';
            limit_reached = 1;
        }

        /* check for the non-printable char limit */
        for(p = line; *p; ++p) {
            if( ! isprint(*p) && *p != '\n' && *p != '\r' && *p != '\t') {
                nonprint_char_count++;
                if(nonprint_char_count >= _PARSE_PARS_NONPRINT_CHAR_LIMIT) {
                    WARN(("Aborting param parsing: %d non-printable "
                          "char(s) found.", nonprint_char_count));
                    *p = '\0';
                    limit_reached = 1;
                    break;
                }
            }
        }

        /* strip comments and trim spaces */
        strip_comment( line, '#' );
        str_trim( line );

        /* prepare source info: "/path/filename:line_no" */
        lineno_str = NULL;
        int stack_size = se_stack_size(file_stack);
        if (stack_size > 0) {
            const char* filename = (const char*) se_stack_top(file_stack);
            if (filename != NULL) {
                lineno_str = asprintf("%s:line %d", filename, lineno);
            }
        }
        if (lineno_str == NULL) {
            lineno_str = int32_to_str(lineno);
        }

        /* parse the line */
        _se_parse_pars_str_internal( line, ht, prefix, file_stack,
                                      PAR_SRCTYPE_FILE, lineno_str );

        free1char(lineno_str);

    } while(!limit_reached);

    return ht;
}

static int _file_stack_compare(const char* filename1,
                                    const char* filename2)
{
    return ! se_is_same_file( filename1, filename2 );
}

static se_hash se_parse_pars_file_internal( const char* file_name,
                                                   se_hash ht,
                                                   const char* prefix,
                                                   se_stack_t* file_stack )
{
    FILE* f;
    char *file_name_abs, *peer, *temp;

    ASSERT(file_name != NULL);

    if(ht == NULL)
        ht = create_hash();

    /* open the file */
    file_name_abs = se_abs_path_from_peer( file_name, "." );
    f = se_fopen( file_name_abs, "r" );
    if(f == NULL) {
        if( se_is_path_relative(file_name) &&
            file_stack &&
            se_stack_size(file_stack) ) {
            /* try again using the base path from the last included file */
            peer = (char*)se_stack_top( file_stack );
            free1char(file_name_abs);
            file_name_abs = se_abs_path_from_peer( file_name, peer );
            f = se_fopen( file_name_abs, "r" );
        }
        if(f == NULL) {
            int err = errno;
            ERROR(("Cannot open file '%s' (current dir is '%s'): errno %d - %s",
                   file_name, get_cwd(), err, strerror(err)));
        }
    }
    TRACE(("Parameter file opened: '%s'", file_name_abs));

    /* ignore already included parameter files */
    if( file_stack &&
        se_stack_contains( file_stack,
                            (compare_values_fn)_file_stack_compare,
                            file_name_abs ) ) {
        WARN(("Already included parameter file ignored: '%s'", file_name_abs));
        fclose(f);
        free1char(file_name_abs);
        return ht;
    }

    /* push the file_name in the stack */
    if( ! file_stack ) {
        file_stack = se_create_stack();
    }

    se_stack_push( file_stack, file_name_abs );

    /* check for max include depth */
    if( se_stack_size(file_stack) > _PARSE_PARS_INCLUDE_MAX_DEPTH) {
        ERROR(("Parameter include maximum depth limit (%d) exceeded in "
               "file '%s' (current dir is '%s')",
               _PARSE_PARS_INCLUDE_MAX_DEPTH, file_name, get_cwd()));
    }

    /* parse the stream */
    ht = _se_parse_pars_fstream_internal( f, ht, prefix, file_stack );

    /* close the stream */
    if( fclose(f) ) {
        int err = errno;
        WARN(("Cannot close stream for file '%s' (current dir is '%s'): "
              "errno %d - %s",
              file_name, get_cwd(), err, strerror(err)));
    }

    /* pop the file_name from the stack */
    temp = (char*)se_stack_pop( file_stack );
    if(temp) free1char(temp);
    if( se_stack_size( file_stack ) == 0 ) {
        se_destroy_stack( file_stack, 0 );
    }

    return ht;
}



static int _is_namedset( const char* par )
{
    const char *eq, *ob;

    if( ! str_endswith( par, "}" ) ) {
        return 0;
    }

    eq = strchr( par, '=' );
    ob = strchr( par, '{' );
    if( eq && ob ) {
        return ob < eq;
    }
    if( ! ob ) {
        return 0;
    }
    if( ! eq ) {
        return ob != NULL;
    }
    return 0;
}

static int _have_prefix(const char* str, int str_len,
                             const char* prefix, int prefix_len)
{
    ASSERT(str);
    ASSERT(prefix);
    ASSERT(str_len >= 0);
    ASSERT(prefix_len >= 0);

    return
        str_len >= prefix_len + 2 &&
        str[0] == '-' &&
        ! strncmp( str + 1, prefix, prefix_len ) &&
        str[prefix_len + 1] == ':';
}

static void _remove_prefix(char* str, int str_len, int prefix_len)
{
    ASSERT(str);
    ASSERT(str_len >= 0);
    ASSERT(prefix_len >= 0);
    ASSERT(prefix_len <= str_len - 2);

    memmove( str, str + prefix_len + 2,
             (str_len - prefix_len - 1) * sizeof(char) );
}


static se_hash se_parse_pars_namedset( const char* _par,
                                              se_hash ht,
                                              const char* prefix,
                                              int src_type,
                                              const char* src_info )
{
    char *par, *ob, *setname, *temp, *temp2;
    array array_values;
    int i;

    par = se_strdup( _par );

    ob = strchr( par, '{' );
    ASSERT(ob);
    *ob = ' ';

    ASSERT( strlen(par) && str_endswith( par, "}" ));
    par[strlen(par) - 1] = '\0';

    array_values = split2( par, ' ',
                               SPLIT_TRIMFIELDS |
                               SPLIT_QUOTES |
                               SPLIT_NO_EMPTYFIELDS );
    ASSERT(array_size( array_values ));

    setname = (char*)array_remove_at( array_values, 0 );
    for (i = 0; i < array_size( array_values ); ++i) {
        temp = (char*)array_get_at( array_values, i );

        temp2 = alloc1char( strlen(setname) + strlen(temp) + 2 );
        temp2[0] = '\0';
        strcat( temp2, setname );
        strcat( temp2, "." );
        strcat( temp2, temp );

        free1char( temp );
        array_set_at( array_values, i, temp2 );
    }
    free1char( setname );

    for (i = 0; i < array_size( array_values ); ++i) {
        temp = (char*)array_get_at( array_values, i );

        se_parse_pars_single( temp, ht, prefix, src_type, src_info );
    }

    destroy_array( array_values, 1 );
    free1char( par );

    return ht;
}

static se_hash se_parse_pars_single( const char* par,
                                            se_hash ht,
                                            const char* prefix,
                                            int src_type,
                                            const char* src_info )
{
    char *par_key, *par_val, *eq, *temp;
    int ignore_par, prefix_len = -1, par_key_len;

    if( _is_namedset(par) ) {
        return se_parse_pars_namedset( par, ht, prefix,
                                        src_type, src_info );
    }


    if( ht == NULL ) {
        ht = create_hash();
    }

    /* split param in key and value */
    eq = (char *)strchr( par, '=' );
    if( ! eq ) {
        /* ignore params if '=' is missing */
        return ht;
    } else {
        par_key = se_strndup( par, (size_t) (eq - par) );
        par_val = se_strdup( eq + 1 );
    }
    str_trim( par_key );
    str_trim( par_val );
    str_unquote( par_val );

    ignore_par = 0;
    do {
        /* check for prefix */
        if (prefix) {
            if( prefix_len == -1 ) {
                prefix_len = strlen(prefix);
            }
            par_key_len = strlen(par_key);
            if(! _have_prefix( par_key, par_key_len, prefix, prefix_len)) {
                /* ignore params that does not start with '-prefix:' */
                ignore_par = 1;
                break;
            } else {
                /* remove the prefix: '-ztr:some_key' > 'some_key' */
            _remove_prefix( par_key, par_key_len, prefix_len );
            }
        } else {
            if (par_key[0] == '-' && strchr(par_key, ':') > par_key) {
                /* ignore params that does start with a prefix ('-prefix:') */
                ignore_par = 1;
                break;
            }
        }

        /* ignore params with empty key ("=value") */
        if (par_key[0] == '\0') {
            if (src_info != NULL && *src_info != '\0') {
                WARN(("%s: Empty parameter ignored (extra spaces before '=' ?)", src_info));
            } else {
                WARN(("Empty parameter ignored (extra spaces before '=' ?)"));
            }
            ignore_par = 1;
            break;
        }
    } while(0);

    if (ignore_par) {
        free1char(par_key);
        free1char(par_val);
        return ht;
    }

    /* decode the value of the param (eg: "a\\=b" > "a=b") */
    temp = backslash_decode( par_val, "\"=", 0 );
    free1char(par_val);
    par_val = temp;

    par_assign_str( ht, par_key, par_val, src_type, src_info );
    free1char(par_key);
    free1char(par_val);

    return ht;
}


static void parnamehost_split( const char* parnamehost,
                                    const char* parval,
                                    char** parname, char** host )
{
    const char* at;
    ASSERT(parnamehost);

    at = strchr( parnamehost, '@' );
    if (at) {
        if (parname) {
            *parname = se_strndup( parnamehost, (size_t) (at - parnamehost) );
        }
        if(host) {
            *host = se_strdup( at + 1 );
            if (**host == '\0') {
                WARN(("Empty hostname string for parameter '%s=%s'",
                      parnamehost, parval));
            }
        }
    } else {
        if (parname) {
            *parname = se_strdup( parnamehost );
        }
        if (host) {
            *host = se_strdup("");
        }
    }
}

/* Reverse hostname alphabetical order comparator for par_t structs */
static int par_revhost_compare(const par_t* p1, const par_t* p2)
{
    return strcmp(p1->host, p2->host) * -1;
}

static void par_assign_str( se_hash ht,
                                 const char* parnamehost,
                                 const char* parval,
                                 int src_type,
                                 const char* src_info )
{
    parinfo_t* pi;
    par_t* p;
    char *parname = NULL, *host = NULL;
    int i, found;

    ht = ht?ht:_ht_params;
    VERIFY(ht);
    ASSERT(parnamehost);
    ASSERT(parval);

    /* split 'parname@host' in 'parname' and 'host' */
    parnamehost_split( parnamehost, parval, &parname, &host );

    /* get the existing parinfo for the current parname */
    pi = (parinfo_t*)ht_get( ht, parname );

    if(pi == NULL) {
        /* new param, add it */
        p = par_create( host );
        par_update( p, parval, src_type, src_info );
        pi = parinfo_create();
        array_add( pi->par_array, p );

        ht_put( ht, parname, pi );

    } else {

        /* existing param. for each host entry */
        for(i = 0, found = 0; i < array_size(pi->par_array); ++i) {
            p = (par_t*)array_get_at( pi->par_array, i );

            if( strcmp(p->host, host) == 0 ) {
                /* update the par when the host matches exactly */
                par_update( p, parval, src_type, src_info );
                found = 1;
            }
        }

        if(! found) {
            /* add a new entry if no exact match has been found */
            p = par_create( host );
            par_update( p, parval, src_type, src_info );
            array_add( pi->par_array, p );

            /* sort the entries in reverse alphabetical order */
            array_sort( pi->par_array,
                            (compare_values_fn) par_revhost_compare );
        }
        free1char(parname);
    }

    free1char(host);
}

/* Splits the specified param value in comma-separated values */
static array parval_split( const char* parval )
{
    array val_array;
    ASSERT(parval);

    /* split the param's value: "1, 2, 3" > "1" "2" "3" */
    val_array = split2( parval, ',', SPLIT_TRIMFIELDS | SPLIT_QUOTES );
    ASSERT(array_size(val_array) > 0);
    return val_array;
}

/* Parses the 'parval' and updates all fields of the specified par_t */
static void par_update( par_t* p, const char* parval,
                             int src_type, const char* src_info )
{
    ASSERT(p);

    free1char(p->original);
    p->original = se_strdup(parval);

    destroy_array(p->val_array, 1);
    p->val_array = parval_split(parval);

    p->src_type = src_type;
    free1char(p->src_info);
    p->src_info = se_strdup(src_info);

    p->update_count++;
}



static int _par_have_var( const char* par_val,
                               char** var_name,
                               int* idx_from,
                               int* idx_to )
{
    char *dollar, *ob, *cb, b;

    dollar =(char *) strchr(par_val, '$');
    if(! dollar)
        return 0;

    if (*(dollar + 1) == '{') {
        ob = dollar + 1;
        b = '}';
    }
    else if (*(dollar + 1) == '(') {
        ob = dollar + 1;
        b = ')';
    }
    else {
        ob = dollar;
        b = ' ';
    }

    if (b == '}' || b == ')') {
        cb = strchr( ob + 1, b );
        if(! cb) {
            WARN(("Missing close parantesis: '%s'", par_val));
            return 0;
        }
    } else {
        for (cb = ob + 1; *cb &&
                 (isalnum(*cb) || *cb=='.' || *cb=='-' || *cb == '_'); ++cb);
        --cb;
    }

    if(cb <= ob + 1) {
        /* empty variable "${}", ignore it */
        return 0;
    }

    if( var_name ) {
        if (b == '}' || b == ')') {
            *var_name = se_strndup( ob + 1, cb - ob - 1 );
        } else {
            *var_name = se_strndup( ob + 1, cb - ob );
        }
    }
    if( idx_from ) {
        *idx_from = dollar - par_val;
    }
    if( idx_to ) {
        *idx_to = cb - par_val;
    }

    return 1;
}

static char* _parval_subst( const char* par_val,
                                 int idx_from,
                                 int idx_to,
                                 const char* var_val )
{
    char* new_par_val;
    strbuf sb;

    sb = create_strbuf();
    sb_appendn( sb, par_val, idx_from );
    sb_append ( sb, var_val );
    sb_append ( sb, par_val + idx_to + 1 );
    new_par_val = sb_release( sb );
    destroy_strbuf( sb );

    return new_par_val;
}

static char* _par_get_fullval( const char* parname,
                                    const char* var_name,
                                    se_hash ht )
{
    array val_array;
    char *val, *fullval;
    int i, val_count;

    val_array = par_lookup( ht, var_name );
    if(! val_array) {
        /* no such param, try for an env variable */
        val = getenv(var_name);
        if(! val) {
            ERROR(("Param '%s' references undefined param or "
                   "environment variable '%s'", parname, var_name));
        }
        return se_strdup(val);
    }

    val_count = array_size(val_array);
    if( val_count == 1 ) {
        val = (char*) array_get_first( val_array );
        fullval = se_strdup( val );
    } else {
        strbuf sb = create_strbuf();
        for (i = 0; i < val_count; ++i ) {
            val = (char*) array_get_at( val_array, i );
            sb_append( sb, val );
            if( i < val_count - 1 ) {
                sb_append( sb, ", " );
            }
        }
        fullval = sb_release( sb );
        destroy_strbuf( sb );
    }
    return fullval;
}

/*
 * Substitutes the param variables from the specified hash table.
 * Eg: "p1=val p2=${p1}" => "p1=val p2=val".
 */
static void se_make_par_subst(se_hash ht)
{
    parinfo_t* pi;
    par_t* p;
    htiter hti;
    const char* parname;
    char *subst_val, *new_subst_val, *var_name, *var_val;
    void* ht_val;
    int i, idx_from, idx_to;
    ASSERT(ht);

    /* for each param info */
    hti = create_htiter( ht );
    while( hti_hasnext( hti )) {
        hti_next( hti, &parname, &ht_val );
        pi = (parinfo_t*) ht_val;

        /* for each par */
        for (i = 0; i < array_size(pi->par_array); ++i) {
            p = (par_t*) array_get_at( pi->par_array, i );

            subst_val = se_strdup(p->original);
            while( _par_have_var( subst_val, &var_name, &idx_from, &idx_to )) {
                /* check for circular deps */
                if(! strcmp(parname, var_name) ) {
                    ERROR(("Circular dependency detected for param '%s'",
                           parname));
                }

                /* substitute the value of var_name in subst_val */
                var_val = _par_get_fullval( parname, var_name, ht );
                free1char( var_name );
                new_subst_val = _parval_subst( subst_val, idx_from,
                                               idx_to, var_val );
                free1char( var_val );
                free1char( subst_val );
                subst_val = new_subst_val;

                /* split subst_val and update the val_array */
                destroy_array(p->val_array, 1);
                p->val_array = parval_split(subst_val);
            }
            free1char(subst_val);
        }
    }
    destroy_htiter( hti );
}


/*
 * Returns the values array for the specified parameter, or NULL if the
 * param is not found for the current hostname.
 */
static array par_lookup( se_hash ht, const char* parname )
{
    parinfo_t* pi;
    ht = ht?ht:_ht_params;
    VERIFY(ht);
    ASSERT(parname);

    /* lookup in the hash for parname */
    pi = (parinfo_t*) ht_get( ht, parname );
    if (pi) {
        /* found, look for the right entry for the current hostname */
        return parinfo_lookup(pi);
    } else {
        /* not found */
        return NULL;
    }
}

/*
 * Returns the values array of the specified parameter if a match
 * for the current hostname if found, or NULL otherwise.
 */
static array parinfo_lookup( const parinfo_t* pi )
{
    char* hostname;
    const par_t* p;
    array val_array;
    int i;
    ASSERT(pi);
    ASSERT(array_size(pi->par_array) > 0);

    /* get the hostname */
    hostname = get_hostname();
    VERIFYM( strcmp(hostname, "localhost") != 0, "Cannot obtain hostname!" );

    /* search for the value array from first partial matching hostname */
    val_array = NULL;
    for (i = 0; i < array_size(pi->par_array); ++i) {
        p = (par_t*) array_get_at( pi->par_array, i );

        if (str_startswith( hostname, p->host )) {
            /* found */
            ASSERT(array_size(p->val_array) > 0);
            val_array = p->val_array;
            break;
        }
    }

    free1char( hostname );
    return val_array;
}

/* Create a deep clone of the specified parinfo_t */
static parinfo_t* parinfo_clone( const parinfo_t* pi )
{
    parinfo_t* pi_clone;
    const par_t *p;
    par_t *p_clone;
    int i;
    ASSERT(pi);

    pi_clone = parinfo_create();
    for (i = 0; i < array_size(pi->par_array); ++i) {
        p = (const par_t*) array_get_at( pi->par_array, i );

        p_clone = par_clone(p);
        array_add( pi_clone->par_array, p_clone );
    }

    return pi_clone;
}

/* Create a deep clone of the specified par_t */
static par_t* par_clone(const par_t* p)
{
    par_t* p_clone;
    const char* value;
    int i;
    ASSERT(p);

    p_clone = (par_t*)malloc( sizeof(*p_clone) );

    p_clone->original = se_strdup(p->original);
    p_clone->val_array = create_array_args(1);
    for (i = 0; i < array_size(p->val_array); ++i) {
        value = (const char*) array_get_at( p->val_array, i );
        array_add( p_clone->val_array, se_strdup(value) );
    }
    p_clone->host = se_strdup(p->host);
    p_clone->src_type = p->src_type;
    p_clone->src_info = se_strdup(p->src_info);
    p_clone->update_count = p->update_count;

    return p_clone;
}

static se_hash _se_parse_pars_str_internal( const char* str,
                                                   se_hash ht,
                                                   const char* prefix,
                                                   se_stack_t* file_stack,
                                                   int src_type,
                                                   const char* src_info )
{
    array array_pars;
    char *par, *temp, *_str, *p;
    int i, done;
    strbuf sb;

    ASSERT(str);

    if (ht == NULL) {
        ht = create_hash();
    }

    /* dupe str and repace tabs with spaces */
    _str = se_strdup(str);
    for (p = _str; *p; ++p) {
        if (*p == '\t') {
            *p = ' ';
        }
    }

    /* split str in space separated params and remove empty params */
    array_pars = split2( _str, ' ', SPLIT_USE_ISSPACE|SPLIT_TRIMFIELDS | SPLIT_QUOTES );
    free1char(_str);
    for (i = 0; i < array_size(array_pars); ++i) {
        par = (char*)array_get_at( array_pars, i );

        /* remove empty params */
        if (par[0] == '\0') {
            temp = (char*)array_remove_at( array_pars, i-- );
            free1char(temp);
        }
    }

    /* merge named sets */
    for (i = 0; i < array_size(array_pars); ++i) {
        par = (char*)array_get_at( array_pars, i );

        if( str_startswith( par, "{" ) ||
            str_count_char( par, '{' ) > str_count_char( par, '}' )) {

            sb = create_strbuf_args( 64 );
            if (str_startswith( par, "{" )) {
                if( i == 0 ) {
                    ERROR(("Cannot parse params: unnamed param set found."));
                }
                temp = (char*)array_get_at( array_pars, i - 1 );
                sb_append( sb, temp );
                sb_append( sb, " " );
                temp = (char*)array_remove_at( array_pars, --i );
                free1char( temp );
            }
            sb_append( sb, par );

            ++i;
            done = 0;
            while (! done && i < array_size(array_pars)) {
                par = (char*)array_remove_at( array_pars, i );
                sb_append( sb, " " );
                sb_append( sb, par );
                if( str_endswith( par, "}" ) &&
                    str_count_char( par, '}' ) >
                    str_count_char( par, '{' ) ) {
                    done = 1;
                }
                free1char( par );
            }
            if(! done) {
                ERROR(("Cannot parse params: unterminated named param "
                       "set (missing '}')."));
            }
            --i;

            temp = (char*)array_get_at( array_pars, i );
            free1char( temp );
            temp = sb_release( sb );

            temp = (char*)array_set_at( array_pars, i, temp );
            destroy_strbuf( sb );
        }
    }

    /* process each param */
    for (i = 0; i < array_size(array_pars); ++i) {
        par = (char*)array_get_at( array_pars, i );

        /* check for an include directive */
        if(! strcmp(par, "include") && i < array_size(array_pars) - 1) {
            par = (char*)array_get_at( array_pars, ++i );
            se_parse_pars_file_internal( par, ht, prefix, file_stack );
            continue;
        }

        se_parse_pars_single( par, ht, prefix, src_type, src_info );
    }

    destroy_array( array_pars, 1 );

    return ht;
}

static char* _get_namedpar( const char* parset, const char* parname )
{
    char* par;
    ASSERT(parname);

    if( parset == NULL ) {
        par = se_strdup( parname );
    } else {
        par = alloc1char( strlen(parset) + strlen(parname) + 2 );
        par[0] = '\0';
        strcat( par, parset );
        strcat( par, "." );
        strcat( par, parname );
    }

    return par;
}

/* Sets the value of the specified parameter. */
static void par_assign_array( se_hash ht,
                                   const char* parnamehost,
                                   array val_array,
                                   int src_type,
                                   const char* src_info )
{
    char *parval, *val;
    int i;

    ASSERT(parnamehost);
    ASSERT(val_array);
    ASSERT(src_info);

    if( array_size(val_array) == 1 ) {
        /* we have a single value, assign the param directly */
        parval = (char*)array_get_first( val_array );
        ASSERT(parval);
        par_assign_str( ht, parnamehost, parval, src_type, src_info );

    } else {

        /* merge all values from the val_array into a single string */
        strbuf sb = create_strbuf();
        for (i = 0; i < array_size(val_array); ++i) {
            val = (char*)array_get_at( val_array, i );
            ASSERT(val);

            sb_appendf(sb, "\"%s\"", val);
            if( i < array_size(val_array) - 1 ) {
                sb_append(sb, ", ");
            }
        }
        parval = sb_release( sb );
        destroy_strbuf(sb);

        /* assign the param */
        par_assign_str( ht, parnamehost, parval, src_type, src_info );
        free1char( parval );
    }
}


static void _pars_to_str_internal( se_hash ht,
                                        const char* name,
                                        char** str,
                                        int* len )
{
    strbuf sb;
    htiter hti;
    const char* parname;
    void* ht_val;
    char *ht_str, *new_str, *value;
    parinfo_t* pi;
    par_t* p;
    int ht_len, i, j, str_len;

    ASSERT(ht);
    ASSERT(name);
    ASSERT(str);
    ASSERT(len);

    sb = create_strbuf_args( 128 );

    /* append a space if the string is not empty */
    if (*len > 0) {
        sb_append(sb, " ");
    }

    /* append encoded hash_name */
    value = backslash_encode( name, _PARS_TO_STR_ENCODE_CHARS, 0 );
    sb_append( sb, value );
    free1char(value);

    sb_append(sb, "{");

    /* for each param */
    hti = create_htiter( ht );
    i = 0;
    while( hti_hasnext(hti) ) {
        hti_next( hti, &parname, &ht_val );
        pi = (parinfo_t*) ht_val;
        ASSERT(pi);

        /* for each host entry */
        for (j = 0; j < array_size(pi->par_array); ++j) {
            p = (par_t*) array_get_at( pi->par_array, j );
            ASSERT(p);

            /* append " " to separate par if more than one host */
            if( j > 0 ) {
                sb_append( sb, " " );
            }

            /* append encoded par name */
            value = backslash_encode(parname, _PARS_TO_STR_ENCODE_CHARS, 1);
            sb_append( sb, value );
            free1char( value );

            /* append encoded par hostname */
            if( strlen(p->host) > 0) {
                value = backslash_encode( p->host,
                                              _PARS_TO_STR_ENCODE_CHARS, 1 );
                sb_appendf( sb, "@%s", value );
                free1char( value );
            }

            /* append encoded par value (the original param value,
               not the parsed values) */
            sb_append( sb, "=" );
            value = backslash_encode( p->original,
                                          _PARS_TO_STR_ENCODE_CHARS, 1 );
            sb_append ( sb, value );
            free1char( value );

            /* append update_count */
            sb_append( sb, "\t" );
            value = int32_to_str( (int32_t) p->update_count );
            sb_append ( sb, value );
            free1char( value );

            /* append src_type */
            sb_append( sb, "\t" );
            value = int32_to_str( (int32_t) p->src_type );
            sb_append ( sb, value );
            free1char( value );

            /* append encoded src_info */
            sb_append( sb, "\t" );
            value = backslash_encode( p->src_info,
                                          _PARS_TO_STR_ENCODE_CHARS, 1 );
            sb_append ( sb, value );
            free1char(value);
        }

        if( ++i < ht_size(ht) ) {
            sb_append( sb, " " );
        }
    }
    destroy_htiter(hti);

    sb_append( sb, "}" );
    ht_len = sb_length( sb );
    ht_str = sb_release( sb );
    destroy_strbuf( sb );

    /* append the string to str and update len */
    str_len = *len > 0 ? strlen(*str) : 0;
    if( *len - str_len >= ht_len ) {
        strcpy( *str + str_len, ht_str );
        *len += ht_len;
    } else {
        *len = ht_len + str_len + 1;
        new_str = alloc1char( *len );
        memcpy( new_str, *str, str_len * sizeof(char) );
        memcpy( new_str + str_len, ht_str, (ht_len + 1) * sizeof(char) );
        if (*str)
            free1char(*str);
        *str = new_str;
    }
    free1char(ht_str);
}


static se_hash _ht_from_str( const char* str )
{
    const char *s1, *s2, *ob, *cb, *sp;
    char *ht_name, *ht_name_decoded, *ht_str, *param;
    se_hash ht, ht_inner, ht_old;
    ASSERT(str);

    ht = create_hash();

    /* split by '}' */
    s1 = str;
    while ((cb = strchr_encoded(s1, '}'))) {
        ob = strchr_encoded(s1, '{');
        if(! ob) {
            ERROR(("se_str_to_pars(): Missing '{'"));
            return NULL;
        }
        ht_inner = create_hash();

        /* parse and decode hash name */
        ht_name = se_strndup( s1, ob - s1 );
        ht_name_decoded = backslash_decode( ht_name,
                                                _PARS_TO_STR_ENCODE_CHARS, 0 );
        free1char( ht_name );

        /* parse hash value */
        ht_str = se_strndup( ob + 1, cb - ob - 1 );

        /* split by ' ' */
        s2 = ht_str;
        while ((sp = strchr_encoded(s2, ' '))) {
            /* parse param */
            param = se_strndup( s2, sp - s2 );
            _par_parse_keyval( param, ht_inner );
            free1char(param);

            s2 = sp + 1;
        }
        /* parse the last key & value */
        param = se_strdup( s2 );
        _par_parse_keyval( param, ht_inner );
        free1char(param);

        /* add to hash */
        free1char(ht_str);
        ht_old = ht_put( ht, ht_name_decoded, ht_inner );
        if (ht_old != NULL) {
            se_destroy_parameters( ht_old );
        }
        s1 = cb + 1;
        while (*s1 == ' ') {
            ++s1;
        }

        /* substitute param variables */
        se_make_par_subst( ht_inner );
    }

    return ht;
}

static void _par_parse_keyval( char* param, se_hash ht_inner )
{
    char *value, *src_info;
    array fields_array;
    int src_type;

    while (*param == ' ') {
        ++param;
    }

    if( strlen(param) == 0 ) {
        return;
    }

    /* split by '\t' */
    fields_array = split2( param, '\t', SPLIT_TRIMFIELDS );

    /* check field count */
    if( array_size( fields_array ) != 4) {
        ERROR(("se_str_to_pars(): Missing '\\t'"));
    }

    value = (char*)array_get_at( fields_array, 1 );

    value = (char*)array_get_at( fields_array, 2 );
    src_type = atoi( value );

    src_info = (char*)array_get_at( fields_array, 3 );

    value = (char*)array_get_at( fields_array, 0 );
    value = backslash_decode( value, _PARS_TO_STR_ENCODE_CHARS, 1 );
    se_parse_pars_single( value, ht_inner, NULL, src_type, src_info );
    free1char(value);

    destroy_array( fields_array, 1 );
}

static void _str_to_pars_internal( se_hash ht_in,
                                        se_hash* ht_out,
                                        int erase_pars,
                                        int overwrite )
{
    htiter hti;
    const char *key, *key_out;
    void* val;
    parinfo_t *pi_in, *pi_out, *old_pi;

    if (erase_pars) {
        se_destroy_parameters( *ht_out );
        *ht_out = create_hash();
    }

    hti = create_htiter( ht_in );
    while( hti_hasnext( hti )) {
        hti_next( hti, &key, &val );
        pi_in = (parinfo_t*) val;

        if( ht_contains_key( *ht_out, key ) &&
            ! overwrite && ! erase_pars ) {
            continue;
        }

        key_out = se_strdup(key);
        pi_out = parinfo_clone(pi_in);
        old_pi = (parinfo_t*) ht_put( *ht_out, key_out, pi_out );
        if (old_pi) {
            free1char((void *)key_out);
            parinfo_destroy(old_pi);
        }
    }
    destroy_htiter(hti);
}


void se_write_pars(FILE* f, int write_regular_pars, int write_system_pars)
{
    if(write_regular_pars) se_write_htpars(f, _ht_params, NULL);
    if(write_system_pars)  se_write_htpars(f, _ht_sys_params, "-se:");
}

void se_write_htpars(FILE* f, const se_hash _htpars, const char* prefix)
{
    htiter hti;
    const char* parname;
    void* ht_val;
    parinfo_t* pi;
    par_t* p;
    int i, j, n;

    se_hash htpars = _htpars?_htpars:_ht_params;

    if (prefix == NULL) {
        prefix = "";
    }

    hti = create_htiter( htpars );
    while( hti_hasnext( hti )) {
        hti_next( hti, &parname, &ht_val );
        pi = (parinfo_t*) ht_val;
        ASSERT(parname);
        ASSERT(pi);

        for (i = 0; i < array_size(pi->par_array); ++i) {
            p = (par_t*) array_get_at( pi->par_array, i );
            ASSERT(p);

            n = array_size(p->val_array);
            ASSERT(n > 0);
            fprintf( f, "%s%s%s%s=", prefix, parname, strlen(p->host) > 0 ? "@" : "", p->host );
            if (n >= 2) {
                fputc( '\"', f );
            }
            for (j = 0; j < n; ++j) {
                const char* val_str = (const char*) array_get_at( p->val_array, j );
                fprintf( f, "%s%s", val_str, (j == n-1) ? "" : "," );
            }
            if (n >= 2) {
                fputc( '\"', f );
            }
            fputc( '\n', f );
        }
    }
    destroy_htiter(hti);
}

int se_compare_user_params( se_hash pars )
{
    return se_compare_ht_params(_ht_params, pars);
}

int se_compare_system_params( se_hash pars )
{
    return se_compare_ht_params(_ht_sys_params, pars);
}

int se_compare_ht_params( se_hash ht_params_actual,
                           se_hash ht_params_saved )
{
    htiter hti;
    const char* key;
    int diff = 0;

    /* for each actual parameter */
    hti = create_htiter( ht_params_actual );
    while( hti_hasnext( hti )) {
        hti_next( hti, &key, NULL );
        if( se_have_htpar(key, ht_params_saved) ) {
            char* val_old = se_get_htpar_str(key, ht_params_saved);
            char* val_new = se_get_htpar_str(key, ht_params_actual);
            if( 0 != strcmp(val_old, val_new) ) {
                diff = 1;
            }
            free(val_new);
            free(val_old);
        } else {
            diff = 1;
        }
    }
    destroy_htiter(hti);

    /* for each old parameter, check if it is in the new list  */
    hti = create_htiter( ht_params_saved );
    while( hti_hasnext( hti )) {
        hti_next( hti, &key, NULL );
        if( ! se_have_htpar(key, ht_params_actual) ) {
            diff = 1;
        }
    }
    destroy_htiter(hti);

    return diff;
}
