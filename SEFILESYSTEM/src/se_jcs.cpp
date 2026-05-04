#include "../include/se_jcs.h"


 /* if null, then JCS was not initialized */
static se_hash _ht_jcs_params = NULL;
static se_hash _ht_saved_params = NULL;
static char* _se_jcs_id = NULL;
static char* _se_jcs_job_dir = NULL;
static char* _se_jcs_dir = NULL;
static char* _se_jcs_param_file = NULL;
static char* _se_jcs_id_file = NULL;
static int _se_jcs_status = SE_JCS_STATUS_NEW;
static int _se_jcs_params_modified = 0;


static void _make_path_abs(char** path)
{
    ASSERT(path);
    ASSERT(*path);

    if( se_is_path_relative(*path)) {
        char* rel_path = *path;
        *path = se_abs_path_from_peer( rel_path, "./somename" );
        free1char(rel_path);
    }
}

static int _jcs_initialized()
{
    return (_ht_jcs_params != NULL);
}

/**
 * Returns the unique identifier of the current job
 */
const char* se_jcs_get_id(void)
{
    if(!_jcs_initialized()) return "";

    ASSERT(_se_jcs_id);
    return _se_jcs_id;
}

char* se_jcs_read_id(void)
{
    char *jcs_dir, *jcs_id_file;
    char* jcs_id;
    char line[1024];
    size_t line_len;
    strbuf sb;
    FILE *f;

    /* create JCS dir name: ".ztr-jcs/<exename>" */
    sb = create_strbuf();
    sb_append( sb, ".ztr-jcs/" );
    sb_append( sb, se_get_exe_name() );
    jcs_dir = sb_release( sb );
    _make_path_abs(&jcs_dir);

    /* create JCS id filename: "<jcs-dir>/id" */
    sb_append( sb, jcs_dir );
    sb_append( sb, "/" );
    sb_append( sb, "id" );
    jcs_id_file = sb_release( sb );
    destroy_strbuf(sb);

    /* read the job id from the file */
    f = se_fopen(jcs_id_file, "r");
    if(! f) {
        ERROR(("Cannot open id file '%s': errno=%d - %s",
                jcs_id_file, errno, strerror(errno)));
    }
    if( NULL == fgets( line, sizeof(line)-1, f ) ) {
        ERROR(("Cannot read the id from file '%s'", jcs_id_file));
    }
    line[sizeof(line)-1] = 0;
    line_len = strlen(line);
    while( line_len > 0 && 
           (line[line_len-1] == '\n' || line[line_len-1] == '\r') ) {
        line[line_len-1] = '\0';
        --line_len;
    }

    if(fclose(f)) {
        ERROR(("Cannot close id file '%s': errno=%d - %s",
                jcs_id_file, errno, strerror(errno)));
    }
    jcs_id = se_strdup(line);

    return jcs_id;
}

const char* se_jcs_get_jobdir(void)
{
    if(_jcs_initialized()) {
        ASSERT(_se_jcs_job_dir);
        return _se_jcs_job_dir;
    } else {
        if( ! _se_jcs_job_dir ) {
            _se_jcs_job_dir = get_cwd();
            _make_path_abs(&_se_jcs_job_dir);
        }
        return _se_jcs_job_dir;
    }
}


/** Use it, esp. if JCS was not initialized */
char* se_jcs_get_jobdir_copy(void)
{
    if(_jcs_initialized()) {
        VERIFY(_se_jcs_job_dir);
    }

    if( ! _se_jcs_job_dir ) {
        char* tmp = get_cwd();
        _make_path_abs(&tmp);
        return tmp;
    } else {
        return se_strdup(_se_jcs_job_dir);
    }
}

const char* se_jcs_get_jcs_dir(void)
{
    if(_jcs_initialized()) {
        ASSERT(_se_jcs_dir);
        return _se_jcs_dir;
    }
    return NULL;
}

/**
 * Returns the JCS status of the current program (one of the
 * SE_JCS_STATUS_* values)
 */
int se_jcs_get_status(void)
{
    ASSERT(_se_jcs_status >= 0);
    return _se_jcs_status;
}

/**
 * Returns true (non-zero) if the actual params (both regular and "-ztr:")
 * are different from the params previously saved when the job was started
 * or then the job has been reset.
 */
int se_jcs_params_modified(void)
{
    if( _jcs_initialized() && se_jcs_get_status() == SE_JCS_STATUS_RESTART ) {
        return _se_jcs_params_modified;
    } else {
        return 0;
    }
}


/**
 * For the current job, it associate a value with a name - these are
 * run-time parameters the program may set; ex. check points.
 * These are not user provided parameters; these are values produced
 * by the program for various internal needs.
 */
void se_set_jcspar_int32(const char* parname, int32_t val)
{
    if(!_jcs_initialized()) return;
    se_set_htpar_int32(parname, _ht_jcs_params, val);
}

void se_set_jcspar_int64(const char* parname, int64_t val)
{
    if(!_jcs_initialized()) return;
    se_set_htpar_int64(parname, _ht_jcs_params, val);
}

void se_set_jcspar_float(const char* parname, float val)
{
    if(!_jcs_initialized()) return;
    se_set_htpar_double(parname, _ht_jcs_params, (double) val);
}

void se_set_jcspar_double(const char* parname, double val)
{
    if(!_jcs_initialized()) return;
    se_set_htpar_double(parname, _ht_jcs_params, val);
}

void se_set_jcspar_real(const char* parname, real val)
{
    if(!_jcs_initialized()) return;
    se_set_htpar_real(parname, _ht_jcs_params, val);
}

void se_set_jcspar_reala(const char* parname, reala val)
{
    if(!_jcs_initialized()) return;
    se_set_htpar_reala(parname, _ht_jcs_params, val);
}

void se_set_jcspar_str(const char* parname, const char* val)
{
    if(!_jcs_initialized()) return;
    se_set_htpar_str(parname, _ht_jcs_params, val);
}


/**
 * The pairs for the set* functions - this is how you retrieve
 * a previously saved value.
 */
int32_t se_get_jcspar_int32(const char* parname)
{
    if(!_jcs_initialized()) return 0;
    return (int32_t) se_get_htpar_int32(parname, _ht_jcs_params);
}

int64_t se_get_jcspar_int64(const char* parname)
{
    if(!_jcs_initialized()) return 0;
    return se_get_htpar_int64(parname, _ht_jcs_params);
}

float se_get_jcspar_float(const char* parname)
{
    if(!_jcs_initialized()) return 0.0;
    return (float) se_get_htpar_double(parname, _ht_jcs_params);
}

double se_get_jcspar_double(const char* parname)
{
    if(!_jcs_initialized()) return 0.0;
    return se_get_htpar_double(parname, _ht_jcs_params);
}

real se_get_jcspar_real(const char* parname)
{
    if(!_jcs_initialized()) return 0.0;
    return se_get_htpar_real(parname, _ht_jcs_params);
}

reala se_get_jcspar_reala(const char* parname)
{
    if(!_jcs_initialized()) return 0.0;
    return se_get_htpar_reala(parname, _ht_jcs_params);
}

char* se_get_jcspar_str(const char* parname)
{
    if(!_jcs_initialized()) return se_strdup("");
    return se_get_htpar_str(parname, _ht_jcs_params);
}


/**
 * Queries the existance of the jcs parameter with the specified name.
 * Returns true (non-zero) if the param is present, or zero otherwise.
 */
int se_have_jcspar(const char* parname)
{
    ASSERT(parname);
    return _jcs_initialized() &&
        ht_contains_key( _ht_jcs_params, parname );
}

/**
 * Returns true (non-zero) if the value or the presence of the specified
 * parameter was modified between after a job restart.
 */
int se_par_modified(const char* parname)
{
    char *value, *old_value;
    int ret;

    ASSERT(parname);
    if(!_jcs_initialized()) return 0;
    if(se_jcs_get_status() != SE_JCS_STATUS_RESTART) return 0;

    value = se_get_par_str(parname);
    old_value = se_get_htpar_str(parname, _ht_saved_params);
    /* get par returns allways a non-NULL value */
    ASSERT(value);
    ASSERT(old_value);

    ret = ( strcmp(old_value, old_value) != 0 );
    free(value);
    free(old_value);
    return ret;
}


/**
 * JCS internal functions
 */
void se_jcs_init(void)
{
    se_hash ht_saved_se_params;
    strbuf sb;
    FILE* f;

    
    if( _se_jcs_job_dir ) {
        free(_se_jcs_job_dir);
        _se_jcs_job_dir = NULL;
    }

    /* create the jcs params hash */
    _ht_jcs_params = create_hash();

    /* ignore all JCS stuff if the "-ztr:jcs_notajob" param is present */
    if( se_have_ztrpar(SEJCS_PARAM_NOTAJOB) ) {
        _se_jcs_status = SE_JCS_STATUS_NOTAJOB;
        return;
    }

    /* ugly like hell, but this is the spec for now ... */
    if(! se_have_par(SEJCS_PARAM_JOBDIR) ) {
        _se_jcs_job_dir = get_cwd();
        _make_path_abs(&_se_jcs_job_dir);
        /* force the job_dir parameter to be set anyway, for pieces of
           code that calls se_get_par */
        se_set_par_str(SEJCS_PARAM_JOBDIR, _se_jcs_job_dir);
    } else {
        _se_jcs_job_dir = se_get_par_str(SEJCS_PARAM_JOBDIR);
        _make_path_abs(&_se_jcs_job_dir);
    }

    /* create and chdir to _se_jcs_job_dir */
    if( se_mkdir(_se_jcs_job_dir,
                  S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH) == CODE_SUCCESS ) {
        int ret = chdir(_se_jcs_job_dir);
        if (ret != 0) {
            WARN(("Failed to change directory to %s: %s", 
                 _se_jcs_job_dir, strerror(errno)));
        }
    }

    /* create JCS dir name: ".ztr-jcs/<exename>" */
    sb = create_strbuf();
    sb_append( sb, ".ztr-jcs/" );
    sb_append( sb, se_get_exe_name() );
    _se_jcs_dir = sb_release( sb );
    _make_path_abs(&_se_jcs_dir);

    /* create JCS param filename: "<jcs-dir>/params" */
    sb_append( sb, _se_jcs_dir );
    sb_append( sb, "/" );
    sb_append( sb, "params" );
    _se_jcs_param_file = sb_release( sb );

    /* create JCS id filename: "<jcs-dir>/id" */
    sb_append( sb, _se_jcs_dir );
    sb_append( sb, "/" );
    sb_append( sb, "id" );
    _se_jcs_id_file = sb_release( sb );
    destroy_strbuf(sb);

    /* determine the JCS status: new, reset or restart */
    if( se_get_separ_int32(SEJCS_PARAM_RESET) ) {
        _se_jcs_status = SE_JCS_STATUS_RESET;
        INFOV((1, "Job Control System: Job has been reset"));
    } else {
        if(! se_io_is_directory(_se_jcs_dir) ) {
            _se_jcs_status = SE_JCS_STATUS_NEW;
            INFOV((1, "Job Control System: New job has been started"));
        } else {
            _se_jcs_status = SE_JCS_STATUS_RESTART;
            INFOV((1, "Job Control System: Job has been restarted"));
            INFOV((1, " jcs directory [%s]",_se_jcs_dir));
        }
    }

    if (_se_jcs_status == SE_JCS_STATUS_NEW ||
        _se_jcs_status == SE_JCS_STATUS_RESET) {
        /* create the JCS directory: ".ztr-jcs/<exename>" */
        if( CODE_SUCCESS != se_mkdir(_se_jcs_dir, 
                                S_IRWXU|S_IRWXG|S_IROTH|S_IXOTH) ) {
            ERROR(("Cannot create directory '%s': errno=%d - %s",
                   _se_jcs_dir, errno, strerror(errno)));
        }

        /* create and fill the JCS param file with regular and ztr params */
        f = se_fopen(_se_jcs_param_file, "w+");
        if(! f) {
            ERROR(("Cannot create file '%s': errno=%d - %s",
                   _se_jcs_param_file, errno, strerror(errno)));
        }
        se_write_pars(f, 1, 1);
        if(fclose(f)) {
            ERROR(("Cannot close file '%s': errno=%d - %s",
                _se_jcs_param_file, errno, strerror(errno)));
        }

        /* generate an unique JCS id */
        _se_jcs_id = generate_uuid();

        /* create and fill the JCS id file */
        f = se_fopen(_se_jcs_id_file, "w+");
        if(! f) {
            ERROR(("Cannot create file '%s': errno=%d - %s",
                   _se_jcs_id_file, errno, strerror(errno)));
        }
        fprintf(f, "%s\n", _se_jcs_id);
        if(fclose(f)) {
            ERROR(("Cannot close file '%s': errno=%d - %s",
                   _se_jcs_id_file, errno, strerror(errno)));
        }

    } else if (_se_jcs_status == SE_JCS_STATUS_RESTART) {
        char line[1024];
        size_t line_len;
        /* check the actual params against saved params */
        _ht_saved_params = se_parse_pars_file( _se_jcs_param_file, NULL,
                                                NULL );
        _se_jcs_params_modified |= se_compare_user_params( _ht_saved_params );

        ht_saved_se_params = se_parse_pars_file( _se_jcs_param_file, NULL,
                                                   "ztr" );
        _se_jcs_params_modified |= se_compare_system_params( ht_saved_se_params );
        se_destroy_parameters( ht_saved_se_params );

        /* read the job id from the file */
        f = se_fopen(_se_jcs_id_file, "r");
        if(! f) {
            ERROR(("Cannot open id file '%s': errno=%d - %s",
                   _se_jcs_id_file, errno, strerror(errno)));
        }
        if( NULL == fgets( line, sizeof(line)-1, f ) ) {
            ERROR(("Cannot read the id from file '%s'", _se_jcs_id_file));
        }
        line[sizeof(line)-1] = 0;
        line_len = strlen(line);
        while( line_len > 0 && 
               (line[line_len-1] == '\n' || line[line_len-1] == '\r') ) {
            line[line_len-1] = '\0';
            --line_len;
        }

        if(fclose(f)) {
            ERROR(("Cannot close id file '%s': errno=%d - %s",
                   _se_jcs_id_file, errno, strerror(errno)));
        }
        _se_jcs_id = se_strdup(line);
    }

    /* TODO: ugly hack */
    setenv("PVM_VMID", _se_jcs_id, 1);
}

void se_jcs_cleanup(void)
{
    char* temp;
    int is_empty;
    if(!_jcs_initialized()) {
        if( _se_jcs_job_dir ) {
            free(_se_jcs_job_dir);
            _se_jcs_job_dir = NULL;
        }
        return;
    }

    /* delete the JCS files and dirs on clean job exit */
    if (se_jcs_get_status() != SE_JCS_STATUS_NOTAJOB &&
        !se_jcs_is_enforcing_no_cleanup() ) {
        INFOV((2, "Cleaning up JCS info"));
        if( unlink(_se_jcs_param_file)) {
            WARN(("Cannot remove file '%s': errno=%d - %s",
                  _se_jcs_param_file, errno, strerror(errno)));
        }
        if( unlink(_se_jcs_id_file)) {
            WARN(("Cannot remove file '%s': errno=%d - %s",
                  _se_jcs_id_file, errno, strerror(errno)));
        }
        if( se_rmrf(_se_jcs_dir) != CODE_SUCCESS ) {
            WARN(("Cannot remove directory '%s': errno=%d - %s",
                  _se_jcs_dir, errno, strerror(errno)));
        }
        temp = se_get_base_path( _se_jcs_dir );
        if( se_dir_is_empty( temp, &is_empty ) == CODE_SUCCESS) {
            if( is_empty ) {
                if( rmdir(temp)) {
                    WARN(("Cannot remove directory '%s': errno=%d - %s",
                          temp, errno, strerror(errno)));
                }
            }
        }
        free1char(temp);
    } else {
        INFOV((2, "JCS info will not be removed"));
    }

    /* destroy the jcs params hash */
    se_destroy_parameters(_ht_jcs_params);
    _ht_jcs_params = NULL;

    /* free some JCS strings */
    if (_se_jcs_id) {
        free1char(_se_jcs_id);
        _se_jcs_id = NULL;
    }
    if (_se_jcs_param_file) {
        free1char(_se_jcs_param_file);
        _se_jcs_param_file = NULL;
    }
    if (_se_jcs_id_file) {
        free1char(_se_jcs_id_file);
        _se_jcs_id_file = NULL;
    }
    if (_se_jcs_dir) {
        free1char(_se_jcs_dir);
        _se_jcs_dir = NULL;
    }
    if (_se_jcs_job_dir) {
        free1char(_se_jcs_job_dir);
        _se_jcs_job_dir = NULL;
    }

    if (_ht_saved_params) {
        se_destroy_parameters( _ht_saved_params );
        _ht_saved_params = NULL;
    }

    INFOV((10, "Finished cleaning up JCS info"));
}


/**
 * This is useful for parameters holding file paths: the convention is
 * to solve the relative path to the job directory.
 * 
 * If the parameter is set, and the value is a relative path, then it
 * solves the path relative to the job directory and sets the new
 * value to the parameter list, so it becomes available for subsequent
 * calls to se_get_par_str().
 *
 * This function doesn't do anything if a parameter is not set or if
 * the value for the parameter is an absolute path.
 *
 * It is safe to call this function repeatedly for the same parameter:
 * once the path is made absolute, the function doesn't do anything.
 */
void se_jcs_set_file_parameter_to_jobdir(const char* parname)
{
    if( ! se_have_par(parname) ) {
        return;
    } else {
        char* file = se_get_par_str(parname);
        if( ! se_is_path_absolute(file) ) {
            char* newfile;
            char* jobdir = se_jcs_get_jobdir_copy();
            newfile = se_make_path_absolute(file, jobdir);
            se_set_par_str(parname, newfile);
            INFOV((2,"Set absolute path %s for parameter %s=%s",newfile,parname,file));
            free(newfile);
            free(jobdir);
        }
        free(file);
    }
}


int se_jcs_is_enforcing_no_cleanup()
{
    if( se_get_separ_int32(SEJCS_NO_CLEANUP) ) {
        return 1;
    }
    return 0;
}
