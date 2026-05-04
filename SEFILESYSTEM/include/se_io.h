#ifndef SE_IO_H
#define SE_IO_H
#include <SEBASIC/include/se_basic.h>
#include <stdio.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <fcntl.h>
#include <unistd.h>
#include <time.h>
#include <errno.h>


/**
 * Remove the file.
 */
int se_io_remove( const char* name );

/**
 * Returns the current _total_ length of the file(s).
 */
int64_t se_io_length( const char* name );


/**
 * Tests whether the file or directory denoted by the path argument
 * exists.
 */
int se_io_exists( const char* path );

/**
 * Tests whether the file or directory denoted by the path argument
 * exists.
 */
int se_io_is_directory( const char* path );


int se_io_readable_file( const char* path );


time_t se_io_mtime(const char* path);

/**
 * Replacement for fopen to do some retrying.
 */
FILE *se_fopen(const char *path, const char *mode);

/** returns true is it is ok to retry to open a file given the errno */
inline int se_retry_file_open_error(int err)
{
    switch(err) {
    case EACCES:
    case EEXIST:
    case EFAULT: 
    case EFBIG:
    case EISDIR:
    case ELOOP:
    case EMFILE:
    case ENAMETOOLONG:
    case ENFILE:
    case ENOMEM:
    case ENOTDIR:
    case EOVERFLOW:
    case EPERM:
    case EROFS:
        return 0;
    }
    return 1;
}



#endif