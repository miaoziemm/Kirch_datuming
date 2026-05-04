#include "../include/se_io.h"

/**
 * Remove the files after closing it.
 */
int se_io_remove( const char* name )
{
    LOG_CONSOLE(SE_INFO, "Removing file '%s'", name);

    if( remove(name) < 0 ) {
        int err = errno;
        if(err != ENOENT) {
            LOG_CONSOLE(SE_ERROR, "Cannot remove file %s: errno %d: %s", 
                            name, err, strerror(err));
            return CODE_ERROR;
        }
    }
    return CODE_SUCCESS;
}

/**
 * Returns the current _total_ length of the file(s).
 */
int64_t se_io_length( const char* name )
{
    struct stat st;

    if( 0 != stat(name, &st) ) {
        LOG_CONSOLE(SE_WARNING,"Cannot stat file %s", name);
        return -1;
    }

    return st.st_size;
}


/**
 * Tests whether the file or directory denoted by the path argument
 * exists.
 */
int se_io_exists( const char* path )
{
    struct stat st;

    if(path == NULL) {
        LOG_CONSOLE(SE_WARNING,"se_io_exists: NULL path");
        ASSERT((0));
        return 0;
    }

    if( 0 != stat(path, &st) ) {
        int err = errno;
        if(err == ENOENT) {
            return 0;
        }
        LOG_CONSOLE(SE_WARNING, "Cannot stat %s: %d - %s", path, err, strerror(err));
        return 0;
    }

    return 1;
}

/**
 * Tests whether the file or directory denoted by the path argument
 * exists.
 */
int se_io_is_directory( const char* path )
{
    struct stat st;

    if(path == NULL) {
        LOG_CONSOLE(SE_WARNING,"se_io_is_directory: NULL path");
        ASSERT((0));
        return 0;
    }

    if( 0 != stat(path, &st) ) {
        int err = errno;
        if(err == ENOENT || err == ENAMETOOLONG) {
            return 0;
        }
        LOG_CONSOLE(SE_WARNING, "Cannot stat %s: %d - %s", path, err, strerror(err));
        return 0;
    }

    return S_ISDIR(st.st_mode);
}

/**
 */
int se_io_readable_file( const char* path )
{
    struct stat st;

    if(path == NULL) {
        LOG_CONSOLE(SE_WARNING,"se_io_readable_file: NULL path");
        ASSERT((0));
        return 0;
    }

    if( 0 != stat(path, &st) ) {
        int err = errno;
        if(err == ENOENT) {
            return 0;
        }
        LOG_CONSOLE(SE_WARNING, "Cannot stat %s: %d - %s", path, err, strerror(err));
        return 0;
    }

    if(S_ISDIR(st.st_mode)) return 0;

    if( access(path, R_OK) == 0) return 1;

    return 0;
}


time_t se_io_mtime(const char* path)
{
    struct stat st;

    if(path == NULL) {
        LOG_CONSOLE(SE_WARNING,"se_io_mtime_file: NULL path");
        ASSERT((0));
        return 0;
    }

    if( 0 != stat(path, &st) ) {
        int err = errno;
        if(err == ENOENT) {
            return 0;
        }
        LOG_CONSOLE(SE_WARNING, "Cannot stat %s: %d - %s", path, err, strerror(err));
        return 0;
    }

    return st.st_mtime;
}

char* se_replace_extension(const char* file, const char* rep)
{
    size_t replen = strlen(rep);
    size_t filelen = strlen(file);
    char* result = alloc1char(filelen + replen + 1);
    char* p;

    if (result == NULL) {
        ERROR(("Failed to allocate memory for file extension replacement"));
        return NULL;
    }

    result[0] = 0;
    strcat(result, file);

    p = strrchr(result, '.');
    if( p != NULL ) *p = 0;

    strcat(result, rep);

    return result;
}

const char* se_get_file_extension(const char* file)
{
    const char* p;

    p = strrchr(file, '.');
    if( p != NULL ) {
        size_t len = strlen(file);
        if((size_t)(p - file) < len-1) return p+1;
    }
    return NULL;
}

FILE *se_fopen(const char *path, const char *mode)
{
    int retry;
    FILE* f;
    int err;

    for(retry = 0; retry < 60; ++retry) {
        f = fopen(path, mode);
        if(f) break;
        err = errno;
        if(!se_retry_file_open_error(err)) break;
        
        LOG_CONSOLE(SE_WARNING, "Cannot open file '%s' in mode %s: current dir is %s, "
                            "error code is %d - %s",
                            path, mode, get_cwd(), err, strerror(err) );

        sleep(1);
    }

    return f;
}
