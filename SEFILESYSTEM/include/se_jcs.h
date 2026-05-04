#ifndef SE_JCS_H
#define SE_JCS_H

#include <SEBASIC/include/se_basic.h>
#include "se_uri.h"
#include "se_par_sep.h"


/**
 * Returns the unique identifier of the current job, or NULL if the
 * current program is not a job.
 */
const char* se_jcs_get_id(void);

/**
 * Returns the unique identifier of the current job, read from the
 * .se_jcs dir.
 */
char* se_jcs_read_id(void);

/**
 * Returns the value of the job_dir parameter. It defaults to the
 * current working directory. If JCS was not initialized, returns
 * simple the "." (current directory).
 */
const char* se_jcs_get_jobdir(void);
/** As above, but you have to free the memory.
    Use it, esp. if JCS was not initialized */
char* se_jcs_get_jobdir_copy(void);

const char* se_jcs_get_jcs_dir(void);


/* Indicates that the current program is a new run of a job */
#define SE_JCS_STATUS_NEW     0
/* Indicates that the current program is a reset of a job */
#define SE_JCS_STATUS_RESET   1
/* Indicates that the current program is a restart of a job */
#define SE_JCS_STATUS_RESTART 2
/* Indicates that the current program is not a job */
#define SE_JCS_STATUS_NOTAJOB 3


/**
 * If a job program launches other programs that are not jobs,
 * then the launched programs should receive this ztr parameter
 */
#define SEJCS_PARAM_NOTAJOB "jcs_notajob"

/**
 * If a job is reset then it will receive this ztr parameter
 */
#define SEJCS_PARAM_RESET "jcs_reset"

/**
 * Name of the parameter pointing to th job directory
 */
#define SEJCS_PARAM_JOBDIR "job_dir"

/**
 * Force JCS not to cleanup the job directory.
 */
#define SEJCS_NO_CLEANUP "jcs_no_cleanup"

/**
 * Standard way of testing if JCS should not perform the final cleanup.
 */
int se_jcs_is_enforcing_no_cleanup();

/**
 * Returns the JCS status of the current program (one of the
 * SE_JCS_STATUS_* values)
 */
int se_jcs_get_status(void);

/**
 * Returns true (non-zero) if the actual params (both regular and "-ztr:")
 * are different from the params previously saved when the job was started
 * or then the job has been reset.
 * Call this function only for restarted jobs! (when se_jcs_get_status()
 * returns SE_JCS_STATUS_RESTART)
 */
int se_jcs_params_modified(void);


/**
 * For the current job, it associate a value with a name - these are
 * run-time parameters the program may set; ex. check points.
 * These are not user provided parameters; these are values produced
 * by the program for various internal needs.
 */
void se_set_jcspar_int32(const char* parname, int32_t val);
void se_set_jcspar_int64(const char* parname, int64_t val);
void se_set_jcspar_float(const char* parname, float val);
void se_set_jcspar_double(const char* parname, double val);
void se_set_jcspar_real(const char* parname, real val);
void se_set_jcspar_reala(const char* parname, reala val);
void se_set_jcspar_str(const char* parname, const char* val);

/**
 * The pairs for the set* functions - this is how you retrieve
 * a previously saved value.
 */
int32_t se_get_jcspar_int32(const char* parname);
int64_t se_get_jcspar_int64(const char* parname);
float se_get_jcspar_float(const char* parname);
double se_get_jcspar_double(const char* parname);
real se_get_jcspar_real(const char* parname);
reala se_get_jcspar_reala(const char* parname);
char* se_get_jcspar_str(const char* parname);

/**
 * Queries the existance of the jcs parameter with the specified name.
 * Returns true (non-zero) if the param is present, or zero otherwise.
 */
int se_have_jcspar(const char* parname);

/**
 * Returns true (non-zero) if the value or the presence of the specified
 * parameter was modified between after a job restart.
 */
int se_par_modified(const char* parname);


void se_jcs_init(void);
void se_jcs_cleanup(void);


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
void se_jcs_set_file_parameter_to_jobdir(const char* parname);


#endif