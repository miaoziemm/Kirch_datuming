#ifndef SE_FS_SEP_FAKE_H
#define SE_FS_SEP_FAKE_H
#include <SEBASIC/include/se_basic.h>

#define SEPLIB_IN_NAME "in"
#define SEPLIB_OUT_NAME "out"

/** true is tag is out or stdout */
inline int seplib_is_out_tag(const char* tag)
{
    return (
            tag &&
            ( 
             ( tag[0] == 'o' && tag[1] == 'u' && tag[2] == 't' && tag[3] == 0
               ) /*out*/ || 
             ( tag[0] == 's' && tag[1] == 't' && tag[2] == 'd' &&
               tag[3] == 'o' && tag[4] == 'u' && tag[5] == 't' && tag[6] == 0 
               ) /*stdout*/
              )
            );
}

/** true is tag is in or stdin */
inline int seplib_is_in_tag(const char* tag)
{
    return (
            tag &&
            ( 
             ( tag[0] == 'i' && tag[1] == 'n' && tag[2] == 0
               ) /*in*/ || 
             ( tag[0] == 's' && tag[1] == 't' && tag[2] == 'd' &&
               tag[3] == 'i' && tag[4] == 'n' && tag[5] == 0 
               ) /*stdout*/
              )
            );
}

/** Call this at the end of the program, before destroying parameters */
void seplib_cleanup(void);


/** 
 * From SEPLib: "get parameter from auxilary file"
 */
int seplib_auxpar( const char *name, const char *type, void* ptr, const char* tag );


/** 
 * From SEPLib: "Grab a parameter from the command line"
 * The return value is not the number of matches found, but rather a
 * boolean if the parameter was found or not.
 */ 
int seplib_getch( const char *name, const char *type, void *ptr );


/** 
 * From SEPLib: "Grab a parameter from the command line or history file"
 * The return value is not the number of matches found, but rather a
 * boolean if the parameter was found or not.
 */ 
inline int seplib_fetch( const char *name, const char *type, void* ptr )
{
    if(seplib_getch( name, type, ptr ))
        return 1;

    return seplib_auxpar( name, type, ptr, SEPLIB_IN_NAME );
}

inline int seplib_hetch(const char* name, const char* type, void* ptr)
{
    return seplib_auxpar( name, type, ptr, SEPLIB_IN_NAME );
}

/**
 * From SEPLib: put parameter into auxilary file.
 */
int seplib_auxputch( const char *name, const char *type, const void* ptr, const char* tag );

/**
 * From SEPLib: put argument in output history file.
 */
inline int seplib_putch ( const char *name, const char *type, const void* val )
{
    return seplib_auxputch(name, type, val, SEPLIB_OUT_NAME);
}

/**
 * From SEPLib: seek to a position in a SEPlib dataset
 */
int64_t seplib_sseek_block(const char *tag, int64_t offset, int64_t block_size, int whence);

/**
 * From SEPLib: seek to a position in a SEPlib dataset
 */
inline int64_t seplib_sseek(const char *tag, int64_t offset, int whence)
{
    return seplib_sseek_block(tag, offset, 1, whence);
}

/**
 * From SEPLib: read in an array from a seplib file
 */
ssize_t seplib_sreed(const char *tag, void *buf, size_t nbytes);

/**
 * From SEPLib: write an array to seplib tag
 */
ssize_t seplib_srite(const char *tag, void *buf, size_t nbytes);

/**
 * From SEPLib: obtain the size of a seplib file
 */
int64_t seplib_ssize_block (const char* tag, int block_size);

/**
 * From SEPLib: Close a SEPlib history file
 */
int seplib_auxclose( const char *tag);


/**
 * From SEPLib: close the seplib output history file
 */
inline int seplib_hclose(void)
{
    return seplib_auxclose(SEPLIB_OUT_NAME);
}

/**
 * From SEPLib: Close a SEPlib history file
 */
inline int seplib_auxhclose( const char *tag)
{
    return seplib_auxclose(tag);
}

/**
 * From SEPLib: Opens auxilary output
 *
 * SEP version returns a FILE, but it doesn't seem to be used, so the
 * safest way is to return nothing.
 */
void seplib_auxout( const char* tag );

/**
 * From SEPLib: Opens auxilary input
 *
 * SEP version returns a FILE, but it doesn't seem to be used, so the
 * safest way is to return nothing.
 */
void seplib_auxin( const char* tag );

/**
 * From SEPLib: Opens file for auxilary input/output
 *
 * SEP version returns a FILE, but it doesn't seem to be used, so the
 * safest way is to return nothing.
 */
void seplib_auxinout( const char* tag );

/**
 * Clears all history, from memory and from the file.
 */
void seplib_aux_clear_history( const char* tag );

/**
 * For the moment this function does nothing.
 */
inline void seplib_set_format( const char *tag, const char *format )
{
    (void)tag;
    (void)format;
}

/**
 * Checks if the file pointed by the par value of tag exists and ce be read.
 */
int seplib_tag_exists( const char* tag );

FILE* seplib_input( void );

int seplib_copy_history( const char *from, const char *to );

#ifdef FAKE_SEPLIB


/** 
 * From SEPLib: "get parameter from auxilary file"
 */
inline int auxpar( const char *name, const char *type, void* ptr, const char* tag )
{
    return seplib_auxpar(name, type, ptr, tag);
}


/** 
 * From SEPLib: "Grab a parameter from the command line"
 * The return value is not the number of matches found, but rather a
 * boolean if the parameter was found or not.
 */ 
inline int getch( const char *name, const char *type, void *ptr )
{
    return seplib_getch(name, type, ptr);
}


/** 
 * From SEPLib: "Grab a parameter from the command line or history file"
 * The return value is not the number of matches found, but rather a
 * boolean if the parameter was found or not.
 */ 
inline int fetch( const char *name, const char *type, void* ptr )
{
    return seplib_fetch(name, type, ptr);
}

inline int hetch(const char* name, const char* type, void* ptr)
{
    return seplib_hetch( name, type, ptr );
}

/**
 * From SEPLib: put parameter into auxilary file.
 */
inline int auxputch( const char *name, const char *type, const void* ptr, const char* tag )
{
    return seplib_auxputch( name, type, ptr, tag );
}

/**
 * From SEPLib: put argument in output history file.
 */
inline int putch ( const char *name, const char *type, const void* val )
{
    return seplib_putch(name, type, val);
}

/**
 * From SEPLib: seek to a position in a SEPlib dataset
 */
inline int64_t sseek_block(const char *tag, int64_t offset, int64_t block_size, int whence)
{
    return seplib_sseek_block(tag, offset, block_size, whence);
}

/**
 * From SEPLib: seek to a position in a SEPlib dataset
 */
inline int64_t sseek(const char *tag, int64_t offset, int whence)
{
    return seplib_sseek(tag, offset, whence);
}

/**
 * From SEPLib: read in an array from a seplib file
 */
inline ssize_t sreed(const char *tag, void *buf, size_t nbytes)
{
    return seplib_sreed(tag, buf, nbytes);
}

/**
 * From SEPLib: write an array to seplib tag
 */
inline ssize_t srite(const char *tag, void *buf, size_t nbytes)
{
    return seplib_srite(tag, buf, nbytes);
}

/**
 * From SEPLib: obtain the size of a seplib file
 */
inline int64_t ssize_block (const char* tag, int block_size)
{
    return seplib_ssize_block(tag, block_size);
}

/**
 * From SEPLib: Close a SEPlib history file
 */
inline int auxclose( const char *tag)
{
    return seplib_auxclose(tag);
}


/**
 * From SEPLib: close the seplib output history file
 */
inline int hclose(void)
{
    return seplib_hclose();
}

inline int auxhclose( const char *tag)
{
    return seplib_auxhclose(tag);
}

/**
 * From SEPLib: Opens auxilary output
 *
 * SEP version returns a FILE, but it doesn't seem to be used, so the
 * safest way is to return nothing.
 */
inline void auxout( const char* tag )
{
    return seplib_auxout(tag);
}

/**
 * From SEPLib: Opens auxilary input
 *
 * SEP version returns a FILE, but it doesn't seem to be used, so the
 * safest way is to return nothing.
 */
inline void auxin( const char* tag )
{
    return seplib_auxin(tag);
}

/**
 * From SEPLib: Opens file for auxilary input/output
 *
 * SEP version returns a FILE, but it doesn't seem to be used, so the
 * safest way is to return nothing.
 */
inline void auxinout( const char* tag )
{
    return seplib_auxinout(tag);
}

/**
 * For the moment this function does nothing.
 */
inline void set_format( const char *tag, const char *format )
{
    return seplib_set_format(tag, format);
}

inline int tag_exists( const char* tag )
{
    return seplib_tag_exists(tag);
}

inline FILE* input( void )
{
    return seplib_input();
}

inline int copy_history( const char *from, const char *to )
{
    return seplib_copy_history( from, to );
}

#else
#ifdef USE_SEP_LIB
#include <seplib.h>
int tag_exists( const char* tag );
#else

/** 
 * From SEPLib: "get parameter from auxilary file"
 */
int auxpar( const char *name, const char *type, void* ptr, const char* tag );


/**
 * From SEPLib: "Grab a parameter from the command line"
 */ 
int getch( const char *name, const char *type, void *ptr );


/** 
 * From SEPLib: "Grab a parameter from the command line or history file"
 */ 
int fetch( const char *name, const char *type, void* ptr );

int hetch(const char* name, const char* type, void* ptr);

/**
 * From SEPLib: put parameter into auxilary file.
 */
int auxputch( const char *name, const char *type, const void* ptr, const char* tag );

/**
 * From SEPLib: put argument in output history file.
 */
int putch ( const char *name, const char *type, const void* val );

/**
 * From SEPLib: seek to a position in a SEPlib dataset
 */
int sseek_block(const char *tag, int offset, int block_size, int whence);

/**
 * From SEPLib: seek to a position in a SEPlib dataset
 */
int sseek(const char *tag, int offset, int whence);

/**
 * From SEPLib: read in an array from a seplib file
 */
int sreed(const char *tag, void *buf, int nbytes);

/**
 * From SEPLib: write an array to seplib tag
 */
int srite(const char *tag, void *buf, int nbytes);

/**
 * From SEPLib: obtain the size of a seplib file
 */
int ssize_block (const char* tag, int block_size);

/**
 * From SEPLib: Close a SEPlib history file
 */
int auxclose( const char *tag);

/**
 * From SEPLib: close the seplib output history file
 */
int hclose(void);

int auxhclose( const char *tag);

/**
 * From SEPLib: Opens auxilary output
 */
void auxout( const char* tag );

/**
 * From SEPLib: Opens auxilary input
 */
void auxin( const char* tag );

/**
 * From SEPLib: Opens file for auxilary input/output
 */
void auxinout( const char* tag );

void set_format( const char *tag, const char *format );

int tag_exists( const char* tag );

FILE* input( void );

int copy_history( const char *from, const char *to );

#endif
#endif

#endif