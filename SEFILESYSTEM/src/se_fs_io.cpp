#include "../include/se_fs_io.h"

#include <limits>


#ifdef ARCH_OSX
#define   LOCK_SH   1    /* shared lock */
#define   LOCK_EX   2    /* exclusive lock */
#define   LOCK_NB   4    /* don't block when locking */
#define   LOCK_UN   8    /* unlock */
#endif

/* few helpers */

#define CHECK_VALID(io, op) do {                             \
        if( io == NULL || io->file_name == NULL ) {          \
            LOG_CONSOLE(SE_ERROR,"I/O structure not valid for %s ", op);  \
        }                                                    \
    } while(0)


/* for printing off_t values */
/* SGI doesn't know about %Ld */
#ifndef IRIX64
#define HAVE_LONG_LONG_PRINTF
#endif

static int se_mul_size_checked(size_t a, size_t b, size_t* out,
                               const char* what, const char* file_name)
{
    if( a != 0 && b > std::numeric_limits<size_t>::max() / a ) {
        LOG_CONSOLE(SE_ERROR, "Size overflow computing %s for file %s", what, file_name);
        return CODE_ERROR;
    }
    *out = a * b;
    return CODE_SUCCESS;
}


se_fsio* se_fsio_init( const char* file_name,
                       const char* mode,
                       const char* data_format )
{
    se_fsio* io;
    int have_plus;
    int retry;


    ASSERT( file_name != NULL );

    /* prepare the defaults */
    if( data_format == NULL ) data_format = "xdr_float";
    if( mode == NULL ) mode = "r";

    io = (se_fsio*)malloc(sizeof(se_fsio));
    if( io == NULL ) {
        LOG_CONSOLE(SE_ERROR, "Out of memory creating se_fsio for file %s", file_name);
        return NULL;
    }
    memset( io, 0, sizeof(se_fsio) );

    /* save the file name and data_format */
    io->file_name = strdup(file_name);
    io->data_format = strdup(data_format);
    io->mode_str = strdup(mode);
    if( io->file_name == NULL || io->data_format == NULL || io->mode_str == NULL ) {
        LOG_CONSOLE(SE_ERROR, "Out of memory duplicating io strings for file %s", file_name);
        if(io->file_name) free(io->file_name);
        if(io->data_format) free(io->data_format);
        if(io->mode_str) free(io->mode_str);
        free(io);
        return NULL;
    }

    /* parse open mode */
    switch( *(mode+1) ) {
    case '+' : have_plus = 1; break;
    case '\0':   have_plus = 0; break;
    default:
        LOG_CONSOLE(SE_ERROR, "Invalid I/O mode: [%s]", mode);
        have_plus = 0;
    }

    if( have_plus && *(mode+2) != '\0' ) {
        LOG_CONSOLE(SE_ERROR, "Invalid I/O mode: [%s]", mode);
    }

    switch( *mode ) {
    case 'r':
        if( have_plus ) {
            /* read/write - preserve data; create if necessary */
            io->open_flags = O_RDWR | O_CREAT;
        } else {
            /* read only */
            io->open_flags = O_RDONLY;
        }
        break;

    case 'w':
        if( have_plus ) {
            /* read/write - truncate; create if necessary */
            io->open_flags = O_RDWR | O_CREAT | O_TRUNC;
        } else {
            /* write only  - truncate; create if necessary */
            io->open_flags = O_WRONLY | O_CREAT | O_TRUNC;
        }
        break;

    case 'a':
        if( have_plus ) {
            /* read/write - append; create if necessary */
            io->open_flags = O_RDWR | O_CREAT | O_APPEND;
        } else {
            /* write only  - append; create if necessary */
            io->open_flags = O_WRONLY | O_CREAT | O_APPEND;
        }
        break;

    default:
        LOG_CONSOLE(SE_ERROR, "Invalid I/O mode: [%s]", mode);
    }


    /* parsing data_format */
    if( strcmp("xdr_float", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_XDR_FLOAT;
    } else if( strcmp("xdr_integer", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_XDR_INTEGER;
    } else if( strcmp("xdr_short", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_XDR_SHORT;
    } else if( strcmp("native_float", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_NATIVE_FLOAT;
    } else if( strcmp("native_integer", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_NATIVE_INTEGER;
    } else if( strcmp("native_short", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_NATIVE_SHORT;
    } else if( strcmp("ibm_float", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_IBM_FLOAT;
    } else {
        ERROR(("Invalid I/O data format: [%s]", io->data_format));
    }


    /* for now open it here - later we may delay
       the open moment -> multi0file IO */
    TRACE(("Open file '%s'. Mode=[%s]", io->file_name, mode));
    
    for(retry = 0;;++retry) {
        int err;
        io->fd = open( io->file_name, io->open_flags,
                       S_IRUSR|S_IWUSR | S_IRGRP|S_IWGRP | S_IROTH|S_IWOTH );
        if( io->fd >= 0 ) break;
        err = errno;
        if(!se_retry_file_open_error(err) || retry > 60 ) {
             ERROR(( "Cannot open file '%s' in mode %s: current dir is %s, "
                    "error code is %d - %s",
                    io->file_name, mode, get_cwd(), err, strerror(err) ));
        }

         WARN(( "Cannot open file '%s' in mode %s: current dir is %s, "
               "error code is %d - %s; Will retry",
               io->file_name, mode, get_cwd(), err, strerror(err) ));
        sleep(10);
    }

    io->read_count = 0;
    io->write_count = 0;

    return io;
}

int se_fsio_parse_mmap_flags(const char* flags, int def)
{
    int flg;
    array a;
    int i;
    if(flags == NULL || *flags == 0) {
        return def;
    }

    a = split2(flags, '|', SPLIT_NO_EMPTYFIELDS | SPLIT_TRIMFIELDS | SPLIT_USE_ISSPACE);

    flg = 0;
    for (i = 0; i < array_size(a); ++i) {
        const char* s = (const char*) array_get_at(a, i);
        if(strcmp(s, "MAP_SHARED") == 0) {
            flg |= MAP_SHARED;
        } else if(strcmp(s, "MAP_PRIVATE") == 0) {
            flg |= MAP_PRIVATE;
#ifdef HAVE_MAP_HUGETLB
        } else if(strcmp(s, "MAP_HUGETLB") == 0) {
            flg |= MAP_HUGETLB;
#endif
#ifdef HAVE_MAP_NORESERVE
        } else if(strcmp(s, "MAP_NORESERVE") == 0) {
            flg |= MAP_NORESERVE;
#endif
#ifdef HAVE_MAP_POPULATE
        } else if(strcmp(s, "MAP_POPULATE") == 0) {
            flg |= MAP_POPULATE;
#endif
#ifdef HAVE_MAP_UNINITIALIZED
        } else if(strcmp(s, "MAP_UNINITIALIZED") == 0) {
            flg |= MAP_UNINITIALIZED;
#endif
#ifdef HAVE_MAP_LOCKED
        } else if(strcmp(s, "MAP_LOCKED") == 0) {
            flg |= MAP_LOCKED;
#endif
        } else {
            char* ptr;
            int f = (int)strtol(s, &ptr, 0);
            if(*s && (*ptr) == 0) { /*the entire string is a valid number */
                flg |= f;
            } else {
                 WARN(("Flag %s is not supported in this executable", s));
            }
        }
    }

    destroy_array(a, 1);
    return flg;
}

int se_fsio_set_mmap( se_fsio* io, const char* flags )
{
    off_t len;
    int flg;
    long psize;
    CHECK_VALID(io, "set mmap");

    if(io->open_flags != O_RDONLY) {
         WARN(("Memory mapping mode for a file that is not opened read-only is not supported: %s",
              io->file_name));
        return CODE_ERROR;
    }

    if(io->mmap_ptr) {
         WARN(("The file %s seem to be already memory mapped", io->file_name));
        if(munmap(io->mmap_ptr, io->mmap_len)) {
            int err = errno;
             WARN(("Problems un-mapping memmory-mapped file %s: errno=%d - %s",
                  io->file_name, err, strerror(err) ));
            return CODE_ERROR;
        }
        io->mmap_ptr = NULL;
        io->mmap_offset = 0;
        io->mmap_len = 0;
    }

    len = se_fsio_length(io);
    psize = sysconf(_SC_PAGE_SIZE);
    if((len%psize) != 0) {
        len = psize*(len/psize) + psize;
         INFOV((100, "Memory mapping file %s: mmap_len changed to %Ld", io->file_name, (long long int)len));
    }
    if((uint64_t)len > ((uint64_t)1<<(8*sizeof(size_t) - 1)) - 1) {
        WARN(("File %s too big for memory mapping. File size is %Ld, but should be maximum %Ld",
              io->file_name, (long long int)len, 
              (long long int)(((uint64_t)1<<(8*sizeof(size_t) - 1)) - 1)));
        return CODE_ERROR;
    }
    
    flg = se_fsio_parse_mmap_flags(flags, MAP_PRIVATE);
    io->mmap_offset = 0;
    io->mmap_len = (size_t)len;
    io->mmap_ptr = mmap(NULL, io->mmap_len, PROT_READ, flg, io->fd, io->mmap_offset);
    if(io->mmap_ptr == MAP_FAILED) {
        int err = errno;
        WARN(("Problems memmory mapping file %s: errno=%d - %s",
              io->file_name, err, strerror(err) ));
        io->mmap_ptr = NULL;
        io->mmap_offset = 0;
        io->mmap_len = 0;
        return CODE_ERROR;
    }
    INFOV((100, "Successful memory mapping file %s with flags %04X", io->file_name, flg));


    return CODE_SUCCESS;
}


/**
 * Close the possible open files and release any resource.  The io
 * structure cannot be used anymore.
 */
int se_fsio_close( se_fsio* io )
{
    int ret = CODE_SUCCESS;

    CHECK_VALID(io, "closing");

    if(io->mmap_ptr) {
        int res = munmap(io->mmap_ptr, io->mmap_len);
        if(res) {
            int err = errno;
            WARN(("Problems un-mapping memmory-mapped file %s: errno=%d - %s",
                  io->file_name, err, strerror(err) ));
            
            ret = CODE_ERROR;
        }
    }

    if( io->fd >= 0 ) {
        int res;
        do {
            res = close( io->fd );
        } while (res < 0 && errno == EINTR);
        if( res < 0 ) {
            int err = errno;
            ERROR(("Problems closing %s: errno=%d - %s",
                   io->file_name, err, strerror(err) ));
            ret = CODE_ERROR;
        }
    }

    se_fsio_detach( io );
    return ret;
}


/**
 * Remove the files after closing it.
 */
void se_fsio_remove( se_fsio* io )
{
    char* file_name;

    CHECK_VALID(io, "removing");
    TRACE(("Removing file '%s'", io->file_name));

    /* first save the file name */
    file_name = strdup(io->file_name);
    if(file_name == NULL) {
        WARN(("Cannot duplicate file name for remove operation."));
        se_fsio_close(io);
        return;
    }

    /* call close */
    se_fsio_close(io);

    if( remove(file_name) < 0 )
        WARN(("Cannot remove file %s.", file_name));

    free( file_name );
}

/**
 * Returns the current _total_ length of the file(s).
 */
off_t se_fsio_length( se_fsio* io )
{
    struct stat st;

    CHECK_VALID(io, "getting length");

    if( 0 != fstat(io->fd, &st) ) {
        ERROR(("Cannot stat file %s", io->file_name));
    }

    return st.st_size;
}

static void _se_fsio_update_pos_from_file_pointer(se_fsio* io)
{
    io->crt_offset = lseek(io->fd, 0, SEEK_CUR);
    if(io->crt_offset < 0) {
        int err = errno;
         WARN(("Problems retreiving the current position for file descriptor %d: errno=%d - %s; will force zero",
              io->fd, err, strerror(err) ));
        lseek(io->fd, 0, SEEK_SET);
        io->crt_offset = 0;
    }
}

/**
 * Attaches an se_fsio structure to the specified file descriptor.
 * The file descriptor should be already opened.
 */
se_fsio* se_fsio_attach( const char* file_name,
                         int fd,
                         const char* data_format )
{
    se_fsio* io;

    ASSERT( file_name != NULL );
    ASSERT( fd >= 0 );

    /* prepare the defaults */
    if( data_format == NULL ) data_format = "xdr_float";

    io = (se_fsio*)malloc(sizeof(se_fsio));
    if( io == NULL ) {
        LOG_CONSOLE(SE_ERROR, "Out of memory attaching se_fsio for file %s", file_name);
        return NULL;
    }
    memset( io, 0, sizeof(se_fsio) );

    /* save the file name and data_format */
    io->file_name = strdup(file_name);
    io->data_format = strdup(data_format);
    io->mode_str = strdup("");
    if( io->file_name == NULL || io->data_format == NULL || io->mode_str == NULL ) {
        LOG_CONSOLE(SE_ERROR, "Out of memory duplicating attached io strings for file %s", file_name);
        if(io->file_name) free(io->file_name);
        if(io->data_format) free(io->data_format);
        if(io->mode_str) free(io->mode_str);
        free(io);
        return NULL;
    }

    /* parse the data_format */
    if( strcmp("xdr_float", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_XDR_FLOAT;
    } else if( strcmp("xdr_integer", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_XDR_INTEGER;
    } else if( strcmp("xdr_short", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_XDR_SHORT;
    } else if( strcmp("native_float", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_NATIVE_FLOAT;
    } else if( strcmp("native_integer", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_NATIVE_INTEGER;
    } else if( strcmp("native_short", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_NATIVE_SHORT;
    } else if( strcmp("ibm_float", io->data_format) == 0 ) {
        io->io_data_type = FIO_DATA_TYPE_IBM_FLOAT;
    } else {
        ERROR(("Invalid I/O data format: [%s]", io->data_format));
    }

    /* initialize the counters */
    io->read_count = 0;
    io->write_count = 0;

    /* copy the file descriptor */
    io->fd = fd;

    _se_fsio_update_pos_from_file_pointer(io);

    return io;
}

/**
 * Dettaches the specified se_fsio structure from it's file descriptor
 * (without closing it).
 */
void se_fsio_detach( se_fsio* io )
{
    CHECK_VALID(io, "detaching");

    free( io->file_name );
    free( io->data_format );
    free( io->mode_str );
    free( io );
}

/**
 * Associates a stream with the io.
 */
FILE* se_fsio_fopen( se_fsio* io )
{
    FILE* f;
    CHECK_VALID(io, "fopen");

    f = fdopen(io->fd, io->mode_str);

    return f;
}

/**
 * Position the current pointer to offset, counted from the beginning
 * of file.  If the file pointer cannot be positioned to the desired
 * offset, the function terminate the process with an CODE_ERROR.
 * Attention when open read-only and seeking behind the file limit!
 */
int se_fsio_seek( se_fsio* io, off_t offset )
{
    off_t pos;

    CHECK_VALID(io, "seeking");

    pos = lseek(io->fd, offset, SEEK_SET );

    if( pos < 0 ) {
        int err = errno;
        ERROR(("Problems seeking to file %s at %Ld: errno=%d - %s",
               io->file_name, offset, err, strerror(err) ));
    }

    else if( pos != offset ) {
         ERROR(( "Cannot seek at position %Ld on file %s.",
                offset, io->file_name ));
    }

    io->crt_offset = pos;

    return CODE_SUCCESS;
}

off_t se_fsio_seek_relative( se_fsio* io, off_t offset, int whence )
{
    off_t pos;

    CHECK_VALID(io, "seeking relative");

    pos = lseek(io->fd, offset, whence );

    if( pos < 0 ) {
        int err = errno;
        ERROR(("Problems relative seeking to file %s at %Ld from %d: errno=%d - %s",
               io->file_name, offset, whence, err, strerror(err) ));
    }

    io->crt_offset = pos;

    return pos;
}

int se_fsio_skip( se_fsio* io, off_t size )
{
    off_t pos;

    CHECK_VALID(io, "skipping");

    ASSERT(size >= 0);

    if(size == 0) return CODE_SUCCESS;

    pos = lseek(io->fd, size, SEEK_CUR );

    if( pos < 0 ) {
        int err = errno;
        if(err == ESPIPE) {
            byte buf[4096];
            off_t crt = 0;
            while(crt < size) {
                size_t bsize;
                if((unsigned long)(size - crt) >= sizeof(buf)) bsize = sizeof(buf);
                else                          bsize = (size_t)(size - crt);
                if( se_fsio_raw_read(io, buf, bsize) != CODE_SUCCESS ) {
                    err = errno;
                    ERROR(( "Problems skipping on pipe %s %LdB: errno=%d - %s",
                            io->file_name, size, err, strerror(err) ));
                }
                crt += bsize;
            }
        } else {
            ERROR(( "Problems skipping on file %s %LdB: errno=%d - %s",
                    io->file_name, size, err, strerror(err) ));
        }
    } else {
        io->crt_offset = pos;
    }

    return CODE_SUCCESS;
}


off_t se_fsio_pos( se_fsio* io )
{
    off_t pos;

    CHECK_VALID(io, "getting pos");

    pos = lseek(io->fd, 0, SEEK_CUR );
    if( pos < 0 ) {
        int err = errno;
        ERROR(("Problems getting position on file %s: errno=%d - %s",
               io->file_name, err, strerror(err) ));
    }
    ASSERTM(pos == io->crt_offset, 
            "The current offset in the file is different than the one from the internal structure");

    return pos;
}
int se_fsio_sync( se_fsio* io )
{
    CHECK_VALID(io, "sync'ing");

#if defined(_POSIX_SYNCHRONIZED_IO) && !defined(__APPLE__)
    if( fdatasync(io->fd) < 0 ) {
        ERROR(( "Problems synchronizing file %s: errno=%d - %s",
                io->file_name, errno, strerror(errno) ));
    }
#else
    if( fsync(io->fd) < 0 ) {
        ERROR(( "Problems synchronizing file %s: errno=%d - %s",
                io->file_name, errno, strerror(errno) ));
    }
#endif

    return CODE_SUCCESS;
}


/**
 * Read 'count' bytes to buff. If this is not possible, the routine
 * terminates the program with an CODE_ERROR.  No conversion is performed
 * on the data.
 */
int se_fsio_raw_read ( se_fsio* io, void *buf, size_t count )
{
    byte* bytes = (byte*)buf;

    CHECK_VALID(io, "reading");

    if(io->mmap_ptr) {
        size_t nbytes;
        if(io->crt_offset < 0) {
            ERROR(("Negative current offset while reading memory mapped file %s", io->file_name));
        }

        if((size_t)io->crt_offset >= io->mmap_len) {
            WARN(("Trying to read past the END of file %s while memory mapped; will return zeros. "
                  "Requested %Ld bytes, starting at %Ld, but only %Ld bytes are mapped from the beginning",
                  io->file_name, (long long int)count, (long long int)io->crt_offset,
                  (long long int)io->mmap_len));
            nbytes = 0;
        } else {
            const size_t avail = io->mmap_len - (size_t)io->crt_offset;
            nbytes = (count < avail) ? count : avail;
        }

        if(nbytes < count) {
            WARN(("Trying to read past the END of file %s while memory mapped; will return zeros. "
                  "Requested %Ld bytes, starting at %Ld, but only %Ld bytes are mapped from the beginning",
                  io->file_name, (long long int)count, (long long int)io->crt_offset,
                  (long long int)io->mmap_len));
        }

        if(nbytes > 0) {
            memcpy(bytes, (byte*)(io->mmap_ptr) + io->crt_offset, nbytes);
        }
        io->read_count += nbytes;
        io->crt_offset += nbytes;
        if(nbytes < count) {
            memset(bytes + nbytes, 0, count-nbytes);
        }
    } else {
        size_t pos = 0;
        while( count > 0 ) {
            ssize_t res;

            do { res = read( io->fd, bytes + pos, count); }
            while( res < 0 && errno == EINTR ); /* retry if just interrupted */

            if( res < 0 ) {
                int err = errno;
                ERROR(( "Problems reading from file %s: errno=%d - %s",
                        io->file_name, err, strerror(err) ));
            }
            if( res == 0 ) {
                ERROR(( "Unexpected end-of-file reading from file %s: "
                        "req %lu, read %lu",
                        io->file_name,
                        (unsigned long)count+pos,
                        (unsigned long)pos));
            }

            count -= res;
            pos   += res;
            io->read_count += res;
            io->crt_offset += res;
        }
    }

    return CODE_SUCCESS;
}

/**
 * Read up to count. Return the number of bytes read.
 */
ssize_t se_fsio_raw_read_available ( se_fsio* io, void *buf, size_t count )
{
    ssize_t res;

    CHECK_VALID(io, "reading");
    if(io->mmap_ptr) {
        size_t nbytes;
        if(io->crt_offset < 0) {
            LOG_CONSOLE(SE_ERROR, "Negative current offset while reading memory mapped file %s", io->file_name);
            return CODE_ERROR;
        }
        if((size_t)io->crt_offset >= io->mmap_len) {
            nbytes = 0;
        } else {
            const size_t avail = io->mmap_len - (size_t)io->crt_offset;
            nbytes = (count < avail) ? count : avail;
        }

        if(nbytes > 0) {
            memcpy(buf, (byte*)(io->mmap_ptr) + io->crt_offset, nbytes);
        }
        res = nbytes;
    } else {
        do { res = read( io->fd, buf, count); }
        while( res < 0 && errno == EINTR ); /* retry if just interrupted */

        if( res < 0 ) {
            int err = errno;
            LOG_CONSOLE(SE_ERROR, "Problems reading from file %s: errno=%d - %s",
                   io->file_name, err, strerror(err) );
        }
    }

    if( res > 0 ) {
        io->read_count += res;
        io->crt_offset += res;
    }

    return res;
}


/**
 * Write 'count' bytes from buff. If this is not possible, the routine
 * terminates the process with an CODE_ERROR.
 * No conversion performed.
 */
int se_fsio_raw_write( se_fsio* io, const void *buf, size_t count )
{
    int retry;
    size_t pos;
    const byte* bytes = (byte*)buf;

    CHECK_VALID(io, "writing");

    pos = 0;
    retry = 0;
    while( count > 0 ) {
        ssize_t res;

        do { res = write( io->fd, bytes + pos, count); }
        while( res < 0 && errno == EINTR ); /* retry if just interrupted */

        if( res < 0 ) {
            int err = errno;
            LOG_CONSOLE(SE_ERROR, "Problems writing to file %s: errno=%d - %s",
                   io->file_name, err, strerror(err) );
            return CODE_ERROR;
        }

        if( res == 0 ) {
            if( retry > 10 ) {
                LOG_CONSOLE(SE_ERROR, "Writing zero bytes to %s after %d trials",
                      io->file_name, retry);
                return CODE_ERROR;
            } else {
                LOG_CONSOLE(SE_WARNING, "Writing zero bytes to %s. Trying again.",
                      io->file_name);
                retry++;
                if( retry%4 == 0 )
                    sleep(1); /* sleep from time to time ... */
                continue;
            }
        }
        retry = 0;

        count -= res;
        pos   += res;
        io->write_count += res;
        io->crt_offset += res;
    }

    return CODE_SUCCESS;
}


/**
 * Fills a file with zeros and set the size of the file to 'size'.  If
 * the file is bigger it will be truncated. If shorter, it will be
 * extended to the new size.  The function write zeros to the file,
 * ensuring this way that the space on disk is allocated.  The file
 * have to be opened in write mode.
 *
 * At the end, the file pointer is positioned at the beginning of the
 * file.
 */
int se_fsio_zero_file( se_fsio* io, off_t size )
{
    byte* buff;
    size_t buff_size;
    off_t old_size;

    CHECK_VALID(io, "filling with zero");

    if( size <= 0 ) {
        if( se_fsio_seek( io, (off_t)0 ) != CODE_SUCCESS )
            LOG_CONSOLE(SE_ERROR, "Cannot seek at the beginning file %s", io->file_name);
        return CODE_SUCCESS;
    }

    if( (size / 10) > 100*1024*1024 )
        buff_size = 100*1024*1024;
    else
        buff_size = (size_t)size / 10; /* size is smaller than 1G */

    if( buff_size == 0 ) buff_size = 1;

    old_size = se_fsio_length(io);

    LOG_CONSOLE(SE_INFO, "Zero file %s. Size:%.2fMB. Old size:%.2fMB. Buffer size:%.2fMB",
          io->file_name,
          (double)size/1024.0/1024.0,
          (double)old_size/1024.0/1024.0,
          (double)buff_size/1024.0/1024.0);

    if( size < old_size ) {
        if( ftruncate( io->fd, size ) < 0 ) {
            int err = errno;
           
            LOG_CONSOLE(SE_ERROR, "Cannot truncate the file %s from %Ld to %Ld:errno=%d - %s",
                  io->file_name, old_size, size, err, strerror(err));
        }
    }

    if( se_fsio_seek( io, (off_t)0 ) != CODE_SUCCESS )
        LOG_CONSOLE(SE_ERROR, "Cannot seek at the beginning file %s", io->file_name);

    buff = (byte*)malloc(buff_size);
    if( buff == NULL ) {
        LOG_CONSOLE(SE_ERROR, "Out of memory allocating zero buffer for file %s", io->file_name);
        return CODE_ERROR;
    }
    memset(buff, 0, buff_size);

    while( size > 0 ) {
        size_t l;
        if( buff_size < (size_t)size) l = buff_size;
        else                  l = (size_t)size;

        if( se_fsio_raw_write( io, buff, l ) != CODE_SUCCESS ) {
            LOG_CONSOLE(SE_ERROR, "Cannot write zeros to file %s", io->file_name);
            free(buff);
            return CODE_ERROR;
        }
        size -= l;
    }

    free(buff);

    if( se_fsio_seek( io, (off_t)0 ) != CODE_SUCCESS )
        LOG_CONSOLE(SE_ERROR, "Cannot seek at the beginning file %s", io->file_name);

    return CODE_SUCCESS;
}

int se_fsio_fake_zero_file( se_fsio* io, off_t size )
{
    byte zero = 0;

    CHECK_VALID(io, "filling with zero");

    if(size <= 0) {
        if( se_fsio_seek( io, (off_t)0 ) != CODE_SUCCESS )
            LOG_CONSOLE(SE_ERROR, "Cannot seek at the beginning file %s", io->file_name);
        return CODE_SUCCESS;
    }

    if( ftruncate( io->fd, 0 ) < 0 ) {
        int err = errno;
        LOG_CONSOLE(SE_ERROR, "Cannot truncate the file %s to zero:errno=%d - %s",
              io->file_name, err, strerror(err));
    }

    if( se_fsio_seek( io, (off_t)0 ) != CODE_SUCCESS )
        LOG_CONSOLE(SE_ERROR, "Cannot seek at the beginning file %s", io->file_name);

    if( se_fsio_raw_write( io, &zero, 1 ) != CODE_SUCCESS )
        LOG_CONSOLE(SE_ERROR, "Cannot write zero at the end file %s", io->file_name);

    if( se_fsio_seek( io, size-1 ) != CODE_SUCCESS )
        LOG_CONSOLE(SE_ERROR, "Cannot seek at the end file %s", io->file_name);

    if( se_fsio_raw_write( io, &zero, 1 ) != CODE_SUCCESS )
        LOG_CONSOLE(SE_ERROR, "Cannot write final zero to file %s", io->file_name);

    if( se_fsio_seek( io, (off_t)0 ) != CODE_SUCCESS )
        LOG_CONSOLE(SE_ERROR, "Cannot seek at the beginning file %s", io->file_name);

    return CODE_SUCCESS;
}

/**
 * Truncates a file.  If the file previously was larger than
 * 'new_size', the extra data is lost.  If the file previously was
 * shorter, it is extended, and the extended part reads as zero bytes.
 * The file pointer is not changed.
 */
int se_fsio_truncate( se_fsio* io, off_t new_size )
{
    CHECK_VALID(io, "truncating");

    if( ftruncate( io->fd, new_size ) != 0 ) {
        int err = errno;
        LOG_CONSOLE(SE_ERROR, "Cannot truncate the file '%s' to %Ld bytes: errno=%d - %s",
              io->file_name, new_size, err, strerror(err) );
    }
    _se_fsio_update_pos_from_file_pointer(io);

    return CODE_SUCCESS;
}


/** Copies the content of iosrc to iodst */
int se_fsio_copy_file( se_fsio* iosrc, se_fsio* iodst )
{
    byte* buff;
    size_t buff_size;
    off_t size, old_size;
    int ret = CODE_SUCCESS;

    CHECK_VALID(iosrc, "copy from");
    CHECK_VALID(iodst, "copy to");

    size = se_fsio_length(iosrc);
    old_size = se_fsio_length(iodst);
    if( size < 0 || old_size < 0 ) {
        LOG_CONSOLE(SE_ERROR, "Invalid file sizes while copying from %s to %s",
              iosrc->file_name, iodst->file_name);
        return CODE_ERROR;
    }

    if( size > 40*1024*1024 )
        buff_size = 40*1024*1024;
    else
        buff_size = (size_t)size; /* size is smaller than 1G */

    if( buff_size == 0 ) buff_size = 1;

    if( size < old_size ) {
        if( ftruncate( iodst->fd, size ) < 0 ) {
            int err = errno;
            LOG_CONSOLE(SE_ERROR, "Cannot truncate the file %s from %Ld to %Ld:errno=%d - %s",
                  iodst->file_name, old_size, size, err, strerror(err));
        }
    }

    if( se_fsio_seek( iosrc, (off_t)0 ) != CODE_SUCCESS )
        LOG_CONSOLE(SE_ERROR, "Cannot seek at the beginning file %s", iosrc->file_name);
    if( se_fsio_seek( iodst, (off_t)0 ) != CODE_SUCCESS )
        LOG_CONSOLE(SE_ERROR, "Cannot seek at the beginning file %s", iodst->file_name);

    buff = (byte*)malloc(buff_size);
    if( buff == NULL ) {
        LOG_CONSOLE(SE_ERROR, "Out of memory allocating copy buffer for file %s", iodst->file_name);
        return CODE_ERROR;
    }
    memset(buff, 0, buff_size);

    while( size > 0 ) {
        size_t l;
        if( buff_size < (size_t)size) l = buff_size;
        else                  l = (size_t)size;

        if( se_fsio_raw_read ( iosrc, buff, l ) != CODE_SUCCESS ) {
            LOG_CONSOLE(SE_ERROR, "Cannot read from file %s", iosrc->file_name);
            ret = CODE_ERROR;
            break;
        }
        if( se_fsio_raw_write( iodst, buff, l ) != CODE_SUCCESS ) {
            LOG_CONSOLE(SE_ERROR, "Cannot write to file %s", iodst->file_name);
            ret = CODE_ERROR;
            break;
        }
        size -= l;
    }

    free(buff);

    if( se_fsio_seek( iosrc, (off_t)0 ) != CODE_SUCCESS )
        LOG_CONSOLE(SE_ERROR, "Cannot seek at the beginning file %s", iosrc->file_name);
    if( se_fsio_seek( iodst, (off_t)0 ) != CODE_SUCCESS )
        LOG_CONSOLE(SE_ERROR, "Cannot seek at the beginning file %s", iodst->file_name);

    return ret;
}

/**
 * Read 'count' bytes to buff.
 * Based on the data_format, it may apply a conversion on raw data.
 * The recommended way is to use a typed routine.
 */
int se_fsio_read ( se_fsio* io, void *buf, size_t count )
{
    int ret = se_fsio_raw_read( io, buf, count );
    if( ret != CODE_SUCCESS ) return ret;

    switch( io->io_data_type ) {
    case FIO_DATA_TYPE_XDR_FLOAT:
    case FIO_DATA_TYPE_XDR_INTEGER:
    case FIO_DATA_TYPE_NATIVE_FLOAT:
    case FIO_DATA_TYPE_NATIVE_INTEGER:
    case FIO_DATA_TYPE_IBM_FLOAT:

        if( count%4 != 0 ) {
            LOG_CONSOLE(SE_WARNING, "Read count is not multiple of 4. File: %s", io->file_name);
            count = count - count%4;
        }

        order_bytes_4( (byte*)buf, count, IO_BYTE_ORDER(io->io_data_type) );
        break;

    case FIO_DATA_TYPE_XDR_SHORT:
    case FIO_DATA_TYPE_NATIVE_SHORT:

        if( count%2 != 0 ) {
            LOG_CONSOLE(SE_WARNING, "Read count is not multiple of 2. File: %s", io->file_name);
            count = count - count%2;
        }

        order_bytes_2( (byte*)buf, count, IO_BYTE_ORDER(io->io_data_type) );
        break;

    case FIO_DATA_TYPE_UNKNOWN:
    default:
        /* do nothing in this case */
        break;
    }

    return CODE_SUCCESS;
}


/**
 * Read 'count' 4B float to buff.
 * Based on the data_format, it tries to convert the raw data to native floats.
 */
int se_fsio_read_float ( se_fsio* io, float *buf, size_t count )
{
    int ret;
    size_t nbytes;
    if( ! IS_FLOAT( io->io_data_type) ) {
        /*WARN(( "Trying to read float while data_format is not float. File %s",
          io->file_name ));*/
    }

    ret = se_mul_size_checked(count, 4, &nbytes, "read_float bytes", io->file_name);
    if(ret != CODE_SUCCESS) return ret;

    ret = se_fsio_raw_read( io, buf, nbytes );
    if( ret != CODE_SUCCESS ) return ret;
    order_bytes_4( (byte*)buf, nbytes, IO_BYTE_ORDER(io->io_data_type) );
    return CODE_SUCCESS;
}

/**
 * Read 'count' 8B float to buff.
 * Based on the data_format, it tries to convert the raw data to native floats.
 */
int se_fsio_read_complex ( se_fsio* io, complex *buf, size_t count )
{
    size_t n;
    if(se_mul_size_checked(count, 2, &n, "read_complex float count", io->file_name) != CODE_SUCCESS)
        return CODE_ERROR;
    return se_fsio_read_float ( io, (float *)buf, n );
}

static int se_fsio_read_complex_part ( se_fsio* io,
                                       se_complex *buf,
                                       size_t count, int part )
{
    while(count > 0) {
        float tmp[2048]; /* assuming we have at least 8KB stack... */
        size_t i;
        size_t nbuf = 2048;
        int ret;

        if(nbuf > count) nbuf = count;

        ret = se_fsio_read_float ( io, tmp, nbuf );
        if( ret != CODE_SUCCESS ) return ret;
        for(i = 0; i < nbuf; ++i) {
            buf[i][part] = tmp[i];
        }
        buf += nbuf;
        count -= nbuf;
    }
    return CODE_SUCCESS;
}

/**
 * Read 'count' 4B float to buff and assign them to the real part.
 * Based on the data_format, it tries to convert the raw data to
 * native floats.
 */
int se_fsio_read_complex_real ( se_fsio* io, se_complex *buf, size_t count )
{
    return se_fsio_read_complex_part ( io, buf, count, 0 );
}

/**
 * Read 'count' 4B float to buff and assign them to the imaginary
 * part.  Based on the data_format, it tries to convert the raw data
 * to native floats.
 */
int se_fsio_read_complex_imag ( se_fsio* io, se_complex *buf, size_t count )
{
    return se_fsio_read_complex_part ( io, buf, count, 1 );
}


/**
 * Read 'count' 2B integers to buff.  Based on the data_format, it
 * tries to convert the raw data to 2B integers.
 */
int se_fsio_read_int16 ( se_fsio* io, int16_t *buf, size_t count )
{
    int ret;
    size_t nbytes;
    if( ! IS_INT16( io->io_data_type) ) {
        /*WARN(( "Trying to read int16 while data_format is not int16. File %s",
          io->file_name ));*/
    }

    ret = se_mul_size_checked(count, 2, &nbytes, "read_int16 bytes", io->file_name);
    if(ret != CODE_SUCCESS) return ret;

    ret = se_fsio_raw_read( io, buf, nbytes );
    if(ret != CODE_SUCCESS) return ret;
    order_bytes_2( (byte*)buf, nbytes, IO_BYTE_ORDER(io->io_data_type) );
    return CODE_SUCCESS;
}


/**
 * Read 'count' 4B integers to buff.  Based on the data_format, it
 * tries to convert the raw data to 4B integers.
 */
int se_fsio_read_int32 ( se_fsio* io, int32_t *buf, size_t count )
{
    int ret;
    size_t nbytes;
    if( ! IS_INT32( io->io_data_type) ) {
        /*WARN(( "Trying to read int32 while data_format is not int32. File %s",
          io->file_name ));*/
    }

    ret = se_mul_size_checked(count, 4, &nbytes, "read_int32 bytes", io->file_name);
    if(ret != CODE_SUCCESS) return ret;

    ret = se_fsio_raw_read( io, buf, nbytes );
    if(ret != CODE_SUCCESS) return ret;
    order_bytes_4( (byte*)buf, nbytes, IO_BYTE_ORDER(io->io_data_type) );
    return CODE_SUCCESS;
}

/**
 * Read 'count' 8B integers to buff.  Based on the data_format, it
 * tries to convert the raw data to 8B integers.
 */
int se_fsio_read_int64 ( se_fsio* io, int64_t *buf, size_t count )
{
    int ret;
    size_t nbytes;
    ret = se_mul_size_checked(count, 8, &nbytes, "read_int64 bytes", io->file_name);
    if(ret != CODE_SUCCESS) return ret;

    ret = se_fsio_raw_read( io, buf, nbytes );
    if(ret != CODE_SUCCESS) return ret;
    order_bytes_8( (byte*)buf, nbytes, IO_BYTE_ORDER(io->io_data_type) );
    return CODE_SUCCESS;
}



/**
 * Write 'count' bytes from buff.  Based on the data_format, it may
 * apply a conversion on raw data.  The recommended way is to use a
 * typed routine.
 */
int se_fsio_write ( se_fsio* io, void *buf, size_t count )
{
    switch( io->io_data_type ) {
    case FIO_DATA_TYPE_XDR_FLOAT:
    case FIO_DATA_TYPE_XDR_INTEGER:
    case FIO_DATA_TYPE_NATIVE_FLOAT:
    case FIO_DATA_TYPE_NATIVE_INTEGER:
    case FIO_DATA_TYPE_IBM_FLOAT:

        if( count%4 != 0 ) {
            LOG_CONSOLE(SE_WARNING, "Write count is not multiple of 4. File: %s", io->file_name);
            count = count - count%4;
        }

        order_bytes_4( (byte*)buf, count, IO_BYTE_ORDER(io->io_data_type) );
        break;

    case FIO_DATA_TYPE_XDR_SHORT:
    case FIO_DATA_TYPE_NATIVE_SHORT:

        if( count%2 != 0 ) {
            LOG_CONSOLE(SE_WARNING, "Write count is not multiple of 2. File: %s", io->file_name);
            count = count - count%2;
        }

        order_bytes_2( (byte*)buf, count, IO_BYTE_ORDER(io->io_data_type) );
        break;

    case FIO_DATA_TYPE_UNKNOWN:
    default:
        /* do nothing in this case */
        break;
    }

    return se_fsio_raw_write( io, buf, count );
}

/**
 * Write 'count' 4B floats from buff.  Based on the data_format, it
 * tries to convert the raw data to specified format.
 */
int se_fsio_write_float ( se_fsio* io, const float *buf, size_t count )
{
    size_t nbytes;
    if( ! IS_FLOAT( io->io_data_type) ) {
        /*WARN(( "Trying to write float while data_format is not float. File %s",
          io->file_name ));*/
    }

    if( NEEDS_SWAPPING(IO_BYTE_ORDER(io->io_data_type)) ) {
        while(count > 0) {
            float tmp[2048];
            int ret, n = 2048;
            if((size_t)n > count) n = (int)count;
            memcpy(tmp, buf, n*4);
            order_bytes_4( (byte*)tmp, n*4,
                               IO_BYTE_ORDER(io->io_data_type) );
            ret = se_fsio_raw_write( io, tmp, 4*n );
            if(ret != CODE_SUCCESS) return ret;
            buf += n;
            count -= n;
        }
        return CODE_SUCCESS;
    } else {
        if(se_mul_size_checked(count, 4, &nbytes, "write_float bytes", io->file_name) != CODE_SUCCESS)
            return CODE_ERROR;
        return se_fsio_raw_write( io, buf, nbytes );
    }
}

/**
 * Write 'count' 8B floats from buff.  Based on the data_format, it
 * tries to convert the raw data to specified format.
 */
int se_fsio_write_complex ( se_fsio* io, se_complex *buf, size_t count )
{
    size_t n;
    if(se_mul_size_checked(count, 2, &n, "write_complex float count", io->file_name) != CODE_SUCCESS)
        return CODE_ERROR;
    return se_fsio_write_float ( io, (float *)buf, n );
}

/**
 * Helper function for writing complex parts
 */
static int se_fsio_write_complex_part ( se_fsio* io, const se_complex *buf,
                                       size_t count, const int part )
{
    while(count > 0) {
        float tmp[2048]; /* assuming we have at least 8KB stack... */
        size_t i;
        size_t nbuf = 2048;
        int ret;
        if(nbuf > count) nbuf = count;
        for(i = 0; i < nbuf; ++i) {
            tmp[i] = buf[i][part];
        }
        ret = se_fsio_write_float ( io, tmp, nbuf );
        if(ret != CODE_SUCCESS) return ret;
        buf += nbuf;
        count -= nbuf;
    }
    return CODE_SUCCESS;
}

/**
 * Write 'count' 4B floats from buff, the real part of the complex
 * numbers Based on the data_format, it tries to convert the raw data
 * to specified format.
 */
int se_fsio_write_complex_real ( se_fsio* io, const se_complex *buf,
                                 size_t count )
{
    return se_fsio_write_complex_part ( io, buf, count, 0 );
}

/**
 * Write 'count' 4B floats from buff, the imaginary part of the
 * complex numbers Based on the data_format, it tries to convert the
 * raw data to specified format.
 */
int se_fsio_write_complex_imag ( se_fsio* io, const se_complex *buf,
                                 size_t count )
{
    return se_fsio_write_complex_part ( io, buf, count, 1 );
}

/**
 * Write 'count' 2B integers from buff.  Based on the data_format, it
 * tries to convert the raw data to the specified format.
 */
int se_fsio_write_int16 ( se_fsio* io, int16_t *buf, size_t count )
{
    size_t nbytes;
    if( ! IS_INT16( io->io_data_type) ) {
        /*WARN(( "Trying to write int16 while data_format is not int16. File %s",
          io->file_name ));*/
    }

    if(se_mul_size_checked(count, 2, &nbytes, "write_int16 bytes", io->file_name) != CODE_SUCCESS)
        return CODE_ERROR;

    order_bytes_2( (byte*)buf, nbytes, IO_BYTE_ORDER(io->io_data_type) );
    return se_fsio_raw_write( io, buf, nbytes );
}


/**
 * Write 'count' 4B integers from buff.  Based on the data_format, it
 * tries to convert the raw data to the specified format.
 */
int se_fsio_write_int32 ( se_fsio* io, int32_t *buf, size_t count )
{
    size_t nbytes;
    if( ! IS_INT32( io->io_data_type) ) {
        /*WARN(( "Trying to write int32 while data_format is not int32. File %s",
          io->file_name ));*/
    }

    if(se_mul_size_checked(count, 4, &nbytes, "write_int32 bytes", io->file_name) != CODE_SUCCESS)
        return CODE_ERROR;

    order_bytes_4( (byte*)buf, nbytes, IO_BYTE_ORDER(io->io_data_type) );
    return se_fsio_raw_write( io, buf, nbytes );
}

/**
 * Write 'count' 8B integers from buff.  Based on the data_format, it
 * tries to convert the raw data to the specified format.
 */
int se_fsio_write_int64 ( se_fsio* io, int64_t *buf, size_t count )
{
    size_t nbytes;
    if(se_mul_size_checked(count, 8, &nbytes, "write_int64 bytes", io->file_name) != CODE_SUCCESS)
        return CODE_ERROR;

    order_bytes_8( (byte*)buf, nbytes, IO_BYTE_ORDER(io->io_data_type) );
    return se_fsio_raw_write( io, buf, nbytes );
}



static void io_lock_internal( se_fsio* io, int op )
{
    int res;

    CHECK_VALID(io, "locking");

    do { res = flock( io->fd, op ); }
    while( res < 0 && errno == EINTR ); /* retry if just interrupted */
    if( res < 0 ) {
        int err = errno;
        LOG_CONSOLE(SE_ERROR, "Problems locking file %s: errno=%d - %s",
              io->file_name, err, strerror(err) );
    }
}

/** file locking */
void se_fsio_lock( se_fsio* io )
{
    io_lock_internal(io, LOCK_EX);
}

/** file un-locking */
void se_fsio_unlock( se_fsio* io )
{
    io_lock_internal(io, LOCK_UN);
}

static int io_inteval_lock_internal( se_fsio* io, int cmd, short type, off_t start, off_t len )
{
    struct flock lock;
    int res;

    CHECK_VALID(io, "locking");

    lock.l_type = type;
    lock.l_start = start;
    lock.l_whence = SEEK_SET;
    lock.l_len = len;

    do { res = fcntl(io->fd, cmd, &lock); }
    while( res < 0 && errno == EINTR ); /* retry if just interrupted */

    return res;
}


void se_fsio_interval_read_lock(se_fsio* io, off_t start, off_t len)
{
    int res = io_inteval_lock_internal(io, F_SETLKW, F_RDLCK, start, len);
    if( res < 0 ) {
        int err = errno;
        LOG_CONSOLE(SE_WARNING, "Problems region locking for reading file %s: "
              "errno=%d - %s; will try flock", io->file_name, err, strerror(err));
        io_lock_internal( io, LOCK_SH );
    }

    if(res < 0) {
        int err = errno;
        LOG_CONSOLE(SE_ERROR, "Problems region locking for reading file %s: errno=%d - %s",
              io->file_name, err, strerror(err) );
    }
}

void se_fsio_interval_write_lock(se_fsio* io, off_t start, off_t len)
{
    int res = io_inteval_lock_internal(io, F_SETLKW, F_WRLCK, start, len);
    if( res < 0 ) {
        int err = errno;
        LOG_CONSOLE(SE_WARNING, "Problems region locking for writing file %s: "
              "errno=%d - %s; will try flock", io->file_name, err, strerror(err));
        io_lock_internal( io, LOCK_EX );
    }

    if(res < 0) {
        int err = errno;
        LOG_CONSOLE(SE_ERROR, "Problems region locking for writing file %s: errno=%d - %s",
              io->file_name, err, strerror(err) );
    }
}

void se_fsio_interval_unlock(se_fsio* io, off_t start, off_t len)
{
    int res = io_inteval_lock_internal(io, F_SETLKW, F_UNLCK, start, len);
    if( res < 0 ) {
        int err = errno;
        LOG_CONSOLE(SE_WARNING, "Problems unlocking region in file %s: "
              "errno=%d - %s; will try flock", io->file_name, err, strerror(err));
        io_lock_internal( io, LOCK_UN );
    }

    if(res < 0) {
        int err = errno;
        LOG_CONSOLE(SE_ERROR, "Problems unlocking region in file %s: errno=%d - %s",
              io->file_name, err, strerror(err) );
    }
}

byte* file_to_memory(const char* fn, int fd, size_t* len)
{
    int need_to_close = (fd == -1);
    if(fd == -1 && (fd = open(fn, O_RDONLY)) == -1) {
        int err = errno;
        LOG_CONSOLE(SE_WARNING, "Cannot read file %s: errno=%d - %s",
              fn, err, strerror(err) );
        return NULL;
    } else {
        size_t buf_increment = 2048;
        size_t buf_len = buf_increment;
        byte* buf = (byte *)malloc(buf_len);
        if( buf == NULL ) {
            LOG_CONSOLE(SE_WARNING, "Out of memory while reading file %s", fn);
            if(need_to_close) close(fd);
            return NULL;
        }
        size_t count;
        lseek(fd, 0L, SEEK_SET);
        for(count = 0;;) {
            ssize_t res;
            do { res = read( fd, buf + count, buf_len-1-count); }
            while( res < 0 && errno == EINTR ); /* retry if just interrupted */
            if( res < 0 ) {
                int err = errno;
                LOG_CONSOLE(SE_WARNING, "Problems reading from file %s: errno=%d - %s",
                      fn, err, strerror(err) );
                free(buf);
                if(need_to_close) close(fd);
                return NULL;
            }
            if(res == 0) break;
            count += res;
            if(count >= buf_len-1) {
                byte* tmp;
                if(buf_len + buf_increment < buf_len) {
                    LOG_CONSOLE(SE_WARNING, "Buffer size overflow while reading file %s", fn);
                    free(buf);
                    if(need_to_close) close(fd);
                    return NULL;
                }
                buf_len += buf_increment;
                tmp = (byte *)realloc(buf, buf_len);
                if(tmp == NULL) {
                    LOG_CONSOLE(SE_WARNING, "Out of memory while expanding buffer for file %s", fn);
                    free(buf);
                    if(need_to_close) close(fd);
                    return NULL;
                }
                buf = tmp;
            }
        }
        buf[count] = 0;
        *len = count;
        if(need_to_close) close(fd);
        return buf;
    }
}

void file_from_memory(const char* fn, int fd, byte* buf, size_t count)
{
    int need_to_close = (fd == -1);
    if(fd == -1 && (fd = open(fn, O_WRONLY | O_CREAT | O_TRUNC,
                              S_IRUSR|S_IWUSR | S_IRGRP|S_IWGRP | S_IROTH|S_IWOTH)) == -1) {
        int err = errno;
        LOG_CONSOLE(SE_ERROR, "Cannot create/overwrite file %s: errno=%d - %s",
              fn, err, strerror(err) );
    } else {
        int retry;
        size_t pos;
        const byte* bytes = (byte*)buf;

        pos = 0;
        retry = 0;
        while( count > 0 ) {
            ssize_t res;

            do { res = write( fd, bytes + pos, count); }
            while( res < 0 && errno == EINTR ); /* retry if just interrupted */

            if( res < 0 ) {
                int err = errno;
                LOG_CONSOLE(SE_ERROR, "Problems writing to file %s: errno=%d - %s",
                      fn, err, strerror(err) );
                break;
            }

            if( res == 0 ) {
                if( retry > 10 ) {
                    LOG_CONSOLE(SE_ERROR, "Writing zero bytes to %s after %d trials", fn, retry);
                    break;
                } else {
                    LOG_CONSOLE(SE_WARNING, "Writing zero bytes to %s. Trying again.", fn);
                    retry++;
                    if( retry%4 == 0 )
                        sleep(1); /* sleep from time to time ... */
                    continue;
                }
            }
            retry = 0;

            count -= res;
            pos   += res;
        }
        if(need_to_close) close(fd);
    }
}

