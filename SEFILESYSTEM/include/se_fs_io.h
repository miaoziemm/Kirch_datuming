#ifndef SE_FS_IO_H
#define SE_FS_IO_H
#include <SEBASIC/include/se_basic.h>
#ifdef ARCH_LINUX
#include <features.h>
#endif

#ifdef ARCH_OSX
#include <sys/syscall.h>
#endif

#include <stdio.h>
#include <string.h>
#include <errno.h>
#include <sys/file.h>
#include <sys/mman.h>
#include <fcntl.h>
#include "se_io.h"


typedef enum {
    FIO_DATA_TYPE_UNKNOWN,        /* unknown data type - no swapping */
    FIO_DATA_TYPE_XDR_FLOAT,      /* IEEE 4B big endian float */
    FIO_DATA_TYPE_XDR_INTEGER,    /* 4B integer, big endian */
    FIO_DATA_TYPE_XDR_SHORT,      /* 2B integer, big endian */
    FIO_DATA_TYPE_NATIVE_FLOAT,   /* IEEE 4B float, little endian */
    FIO_DATA_TYPE_NATIVE_INTEGER, /* 4B integer, little endian */
    FIO_DATA_TYPE_NATIVE_SHORT,   /* 2B integer, little endian*/
    FIO_DATA_TYPE_IBM_FLOAT       /* IBM 4B float, big endian */
} fio_data_type_t;

typedef struct se_fsio_s
{
    char* file_name;         //文件名
    char* data_format;       //数据格式
    char* mode_str;      //模式字符串

    /** The file descriptor */
    int fd;  //文件描述符

    /** type size - so it will correctly swap the data */  //枚举量 大头 小头 各种格式
    fio_data_type_t io_data_type;

    /** flags to be used when open the file */
    int open_flags;

    /** related to memory-mapped files */
    void*  mmap_ptr;    /** if not NULL, then we have a memory mapped file, and this is the pointer */
    off_t  mmap_offset; /** offset in the file where the memory mapping starts */
    size_t mmap_len;    /** lenght of the memory mapped buffer */
    off_t  crt_offset;  /** this is to update the current file pointer, maintainig the API interface */

    /* for statistics */
    off_t read_count;
    off_t write_count;

}se_fsio;


#define IS_FLOAT(dt) ( (dt) == FIO_DATA_TYPE_XDR_FLOAT    ||     \
                       (dt) == FIO_DATA_TYPE_NATIVE_FLOAT ||     \
                       (dt) == FIO_DATA_TYPE_IBM_FLOAT )

#define IS_INT16(dt) ( (dt) == FIO_DATA_TYPE_XDR_SHORT    ||     \
                       (dt) == FIO_DATA_TYPE_NATIVE_SHORT )

#define IS_INT32(dt) ( (dt) == FIO_DATA_TYPE_XDR_INTEGER    ||   \
                       (dt) == FIO_DATA_TYPE_NATIVE_INTEGER )

#define IO_BYTE_ORDER(dt) ( ( (dt) == FIO_DATA_TYPE_NATIVE_FLOAT    ||   \
                              (dt) == FIO_DATA_TYPE_NATIVE_INTEGER  ||   \
                              (dt) == FIO_DATA_TYPE_NATIVE_SHORT ) ?     \
                            little_endian : big_endian )

/* All function will exit in case of error - no need to check the result */

/**
 * Initialize the IO structures. It may or may not open the file(s).
 * Since one IO structure can control the access to more than one file,
 * there is no open function - this is done as needed.
 *
 * \param mode specifies the access mode in which the file(s) are to be opened.
 *        The permitted values and their meanings are similar, but not
 *        exactly the same to fopen. If \c NULL is passed, it is the same as \c "r" :
 *        - \b "r"  Open for reading only. Invoking any of the write functions
 *                  will terminate the application with error.
 *                  If the file(s) does not exist, the program will terminate with
 *                  error.
 *        - \b "r+" Open for updating - reading and writing - all existing data is
 *                  preserved. If the file(s) does not already exist then an attempt
 *                  will be made to create them.
 *        - \b "w"  Open for write only. Invoking any of the read functions
 *                  will terminate the application with error.
 *                  Truncate to zero length, if file(s) exist.
 *                  If the file(s) does not already exist then an attempt will be
 *                  made to create them.
 *        - \b "w+" Open for read and write.
 *                  Truncate to zero length, if file(s) exist.
 *                  If the file(s) does not already exist then an attempt will be
 *                  made to create them.
 *        - \b "a"  Open for append: write only at the end-of-file. Invoking any of
 *                  the read functions will terminate the application with error.
 *                  If the file(s) does not already exist then an attempt will be
 *                  made to create them.
 *        - \b "a+" Open for update: read any place, write only at the end-of-file
 *                  If the file(s) does not already exist then an attempt will be
 *                  made to create them.
 *
 * \param data_format specifies the type of data in the file.
 *   The typed read functions (fio_read_\<type\>) will convert to the native
 *   representation of \<type\>.
 *   The \c read routine will also swap if necessary - the number of bytes to
 *   swap depends on the value of data_format.
 *   The \c raw_read routine does not perform any swapping.
 *   The \c write routine will convert from native to the appropriate format
 *   before writing to disk.
 *   The \c write_raw does not perform any conversion.
 *   Please note that \c NULL is equivalent with \c "xdr_float".
 *   Possible values:
 *     - \c xdr_float      - IEEE 4B float, big endian
 *     - \c xdr_integer    - 4B integer, big endian
 *     - \c xdr_short      - 2B integer, big endian
 *     - \c native_float   - IEEE 4B float, little endian
 *     - \c native_integer - 4B integer, little endian
 *     - \c native_short   - 2B integer, little endian
 *     - \c ibm_float      - IBM 4B float, big endian
 *
 * \param file_name the path to the file.
 *
 * Any created files will have mode
 * S_IRUSR|S_IWUSR|S_IRGRP|S_IWGRP|S_IROTH|S_IWOTH (0666), as modified by
 * the process' umask value. A typical example is when umask is 022, in which
 * case the file will have the permissions 0666 & ~022 =0644 = rw-r--r--.
 */
se_fsio* se_fsio_init( const char* file_name,
                       const char* mode,
                       const char* data_format );

int se_fsio_set_mmap( se_fsio* io, const char* flags );

/**
 * Close the possible open files and release any resource.
 * The io structure cannot be used anymore.
 */
int se_fsio_close( se_fsio* io );

/**
 * Remove the files after closing it.
 */
void se_fsio_remove( se_fsio* io );

/**
 * Returns the current _total_ length of the file(s).
 */
off_t se_fsio_length( se_fsio* io );

/**
 * Attaches an fio structure to the specified file descriptor.
 * The file descriptor should be already opened.
 */
se_fsio* se_fsio_attach( const char* file_name,
                       int fd,
                       const char* data_format );

/**
 * Dettaches the specified fio structure from it's file descriptor
 * (without closing it).
 */
void se_fsio_detach( se_fsio* io );

/**
 * Associates a stream with the io.
 */
FILE* se_fsio_fopen( se_fsio* io );


/**
 * Position the current pointer to offset, counted from the beginning of file.
 * If the file pointer cannot be positioned to the desired offset, the function
 * terminate the process with an error.
 * Attention when open read-only and seeking behind the file limit!
 */
int se_fsio_seek( se_fsio* io, off_t offset );

off_t se_fsio_seek_relative( se_fsio* io, off_t offset, int whence );

inline off_t se_fsio_current_pos( se_fsio* io )
{
    return se_fsio_seek_relative( io, 0, SEEK_CUR );
}


int se_fsio_skip( se_fsio* io, off_t size );

off_t se_fsio_pos( se_fsio* io );

int se_fsio_sync( se_fsio* io );

/**
 * Fills a file with zeros and set the size of the file to 'size'.  If
 * the file is bigger it will be truncated. If shorter, it will be
 * extended to the new size.  The function write zeros to the file,
 * ensuring this way that the space on disk is allocated.  The file
 * have to be opened in write mode.
 *
 * At the end, the file pointer is positioned at the beginning of the file.
 */
int se_fsio_zero_file( se_fsio* io, off_t size );
int se_fsio_fake_zero_file( se_fsio* io, off_t size );

/**
 * Truncates a file.  If the file previously was larger than
 * 'new_size', the extra data is lost.  If the file previously was
 * shorter, it is extended, and the extended part reads as zero bytes.
 * The file pointer is not changed.
 */
int se_fsio_truncate( se_fsio* io, off_t new_size );

/** Copies the content of iosrc to iodst */
int se_fsio_copy_file( se_fsio* iosrc, se_fsio* iodst );

/**
 * Read 'count' bytes to buff. If this is not possible, the routine
 * terminates the program with an error.  No conversion is performed
 * on the data.
 */
int se_fsio_raw_read ( se_fsio* io, void *buf, size_t count );

/**
 * Read up to count. Return the number of bytes read.
 */
ssize_t se_fsio_raw_read_available ( se_fsio* io, void *buf, size_t count );

/**
 * Read 'count' bytes to buff.
 * Based on the data_format, it may apply a conversion on raw data.
 * The recommended way is to use a typed routine.
 */
int se_fsio_read ( se_fsio* io, void *buf, size_t count );

/**
 * Read 'count' 4B float to buff.
 * Based on the data_format, it tries to convert the raw data to native floats.
 */
int se_fsio_read_float ( se_fsio* io, float *buf, size_t count );

/**
 * Read 'count' 8B float to buff.
 * Based on the data_format, it tries to convert the raw data to native floats.
 */
int se_fsio_read_complex ( se_fsio* io, complex *buf, size_t count );

/**
 * Read 'count' 4B float to buff and assign them to the real part.
 * Based on the data_format, it tries to convert the raw data to native floats.
 */
int se_fsio_read_complex_real ( se_fsio* io, complex *buf, size_t count );

/**
 * Read 'count' 4B float to buff and assign them to the imaginary
 * part.  Based on the data_format, it tries to convert the raw data
 * to native floats.
 */
int se_fsio_read_complex_imag ( se_fsio* io, complex *buf, size_t count );

/**
 * Read 'count' 2B integers to buff.
 * Based on the data_format, it tries to convert the raw data to 2B integers.
 */
int se_fsio_read_int16 ( se_fsio* io, int16_t *buf, size_t count );

/**
 * Read 'count' 4B integers to buff.
 * Based on the data_format, it tries to convert the raw data to 4B integers.
 */
int se_fsio_read_int32 ( se_fsio* io, int32_t *buf, size_t count );

/**
 * Read 'count' 8B integers to buff.
 * Based on the data_format, it tries to convert the raw data to 8B integers.
 */
int se_fsio_read_int64 ( se_fsio* io, int64_t *buf, size_t count );


/**
 * Write 'count' bytes from buff. If this is not possible, the routine
 * terminates the process with an error.  No conversion performed.
 */
int se_fsio_raw_write( se_fsio* io, const void *buf, size_t count );

/**
 * Write 'count' bytes from buff.
 * Based on the data_format, it may apply a conversion on raw data.
 * The recommended way is to use a typed routine.
 */
int se_fsio_write ( se_fsio* io, void *buf, size_t count );

/**
 * Write 'count' 4B floats from buff.
 * Based on the data_format, it tries to convert the raw data to
 * specified format.
 */
int se_fsio_write_float ( se_fsio* io, const float *buf, size_t count );

/**
 * Write 'count' 8B floats from buff.
 * Based on the data_format, it tries to convert the raw data to
 * specified format.
 */
int se_fsio_write_complex ( se_fsio* io, complex *buf, size_t count );

/**
 * Write 'count' 4B floats from buff, the real part of the complex
 * numbers Based on the data_format, it tries to convert the raw data
 * to specified format.
 */
int se_fsio_write_complex_real ( se_fsio* io, const complex *buf, size_t count );

/**
 * Write 'count' 4B floats from buff, the imaginary part of the
 * complex numbers Based on the data_format, it tries to convert the
 * raw data to specified format.
 */
int se_fsio_write_complex_imag ( se_fsio* io, 
                                 const complex *buf, size_t count );

/**
 * Write 'count' 2B integers from buff.  Based on the data_format, it
 * tries to convert the raw data to the specified format.
 */
int se_fsio_write_int16 ( se_fsio* io, int16_t *buf, size_t count );

/**
 * Write 'count' 4B integers from buff.  Based on the data_format, it
 * tries to convert the raw data to the specified format.
 */
int se_fsio_write_int32 ( se_fsio* io, int32_t *buf, size_t count );

/**
 * Write 'count' 8B integers from buff.  Based on the data_format, it
 * tries to convert the raw data to the specified format.
 */
int se_fsio_write_int64 ( se_fsio* io, int64_t *buf, size_t count );


/** file locking */
void se_fsio_lock( se_fsio* io );

/** file un-locking */
void se_fsio_unlock( se_fsio* io );

void se_fsio_interval_read_lock(se_fsio* io, off_t start, off_t len);

void se_fsio_interval_write_lock(se_fsio* io, off_t start, off_t len);

void se_fsio_interval_unlock(se_fsio* io, off_t start, off_t len);


byte* file_to_memory(const char* fn, int fd, size_t* len);

void file_from_memory(const char* fn, int fd, byte* buf, size_t count);

/** Helper function to parse MAP_ flags from a string. Example "MAP_PRIVATE|MAP_POPULATE" */
int se_fsio_parse_mmap_flags(const char* flags, int def);


#endif