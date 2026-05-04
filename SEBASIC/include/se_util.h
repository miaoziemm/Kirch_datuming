#ifndef SE_UTIL_H
#define SE_UTIL_H
#include "se_hash.h"
#include "se_array.h"
#include "se_minmax.h"
#include "se_assert.h"
#include "se_return_code.h"
#include "sqlite3.h"
#include <math.h>
#include <stdlib.h>
#include <errno.h>
#include <stdarg.h>
#include <stdio.h>
#include <string.h>
#include <ctype.h>
#include <limits.h>
#include <unistd.h>
#include <time.h>
#include <signal.h>
#include <locale.h>
#include <sys/wait.h>
/* few usefull definitions */

#define ARRAY_SIZE(x) (sizeof(x)/sizeof(x[0]))

#define POWER_OF_TWO(n) ( ( (n) > 0 ) && ( ( (n) & ((n) - 1) ) == 0 ) )

#define UNUSED(v) ((void)(v))

/**
 * Checks whether a string matches a given wildcard pattern. Only does '?'
 * and '*', and multiple patterns separated by '|'.
 * '*' are parsed recursively.
 * '*' does not stands for 'nothing' - it is more like + from reg. exp. :
 *  For example test* matches test1, test2 but doesn't matches test.
 * Later version should take care of this.
 **/
int pattern_match( const char* pattern, const char* string );


/**
 * Splits a string into fields using a char separator.
 *
 * Flags: SPLIT_NO_EMPTYFIELDS - don't return empty fields
 *        SPLIT_TRIMFIELDS - trim fields (removes control and white spaces)
 *        SPLIT_QUOTES - don't split the text inside double-quotes \",
 *                       escape sequences (like \") are not supported.
 *        SPLIT_COLLAPSESEP - 1 or more separators would be treated as one
 *
 *
 * For example splitting " aa: bb::c:" with the ':' separator results in:
 * - " aa" " bb" "" "c" "" with SPLIT_NOFLAGS flag,
 * - "aa" "bb" "" "c" "" with SPLIT_TRIMFIELDS flag,
 * - " aa" " bb" "c" with SPLIT_NO_EMPTYFIELDS flags.
 * - "aa" "bb" "c" with SPLIT_TRIMFIELDS and SPLIT_NO_EMPTYFIELDS flags.
 **/

#define SPLIT_NOFLAGS         (0x00000000u)
#define SPLIT_NO_EMPTYFIELDS  (0x00000001u)
#define SPLIT_TRIMFIELDS      (0x00000002u)
#define SPLIT_QUOTES          (0x00000004u)
#define SPLIT_COLLAPSESEP     (0x00000008u)
#define SPLIT_USE_ISSPACE     (0x00000010u)

void split( const char* text,     /* text to be splitted in fields */
                char sep,             /* field separator */
                char*** fields,       /* array of fields (out) */
                int* field_nr,        /* size of the fields array (out) */
                unsigned int flags ); /* flags: bitmask of SPLIT_* macros */

/**
 * Splits a string into fields using a char separator.
 * Same as split() function from above, but returns an array with the fields.
 **/
array split2( const char* text,     /* text to be splitted in fields */
                      char sep,             /* field separator */
                      unsigned int flags ); /* flags: bitmask of SPLIT_* macros */


/**
 * Strips the comment from the given line of text. 'ch' is beginning-of-
 * comment character.
 **/
void strip_comment( char* line,         /* line to strip the comment from (in-out)*/
                        char comment_ch );  /* beginning of comment char */

/**
 * Strips the prefix (if found) from the specified string.
 * Returns ture (non-zero) if the prefix is present, or zero otherwise.
 * eg: strip_prefix("error: blah", "error: ") => "blah"
 **/
int strip_prefix( char* str, const char* prefix );

/**
 * Remove control and white spaces from the beginning and from the end
 * of the specified string.
 **/
void str_trim( char* line );          /* line to trim (in-out) */


/**
 * Parse the string argument as a signed decimal integer.
 * All chars in the string must be decimal digits, except for the sign char.
 **/
int parse_int( const char* str,   /* string to be parsed */
                   int* val );        /* parsed int value (out) */

/**
 * Parse the string argument as a double.
 * All chars in the string must be part of the number.
 **/
int parse_double( const char* str,   /* string to be parsed */
                      double* val );     /* parsed double value (out) */


/**
 * Same as the standard strdup(), but uses alloc_char() to allocate the memory.
 * The caller should call free on the returned string.
 **/
char* se_strdup( const char* str );

/**
 * Same as strndup(), but uses alloc_char() to allocate the memory.
 * The caller should call free on the returned string.
 **/
char* se_strndup( const char* str, size_t n );

/**
 * Returns the minimum between the len of the 'str' string and 'n'.
 * Renamed to avoid conflicts with system strnlen
 */
size_t se_strnlen( const char* str, size_t n );

/**
 * Remove quotes (") from the beginning and end of a string.
 * Eg: ["aa bb"] > [aa bb], [aa bb"] > [aa bb"]
 **/
void str_unquote( char* str );

/**
 * Returns true (non-zero) if the specified string starts with the specified prefix,
 * or zero otherwise.
 **/
int str_startswith( const char* str, const char* prefix );

/**
 * Returns true (non-zero) if the specified string ends with the specified suffix,
 * or zero otherwise.
 **/
int str_endswith( const char* str, const char* suffix );

/**
 * Counts how many 'c' chars appears in the 'str' string.
 */
int str_count_char(const char* str, char c);

/**
 * Returns the index of the first occurence of 'c' in 's' - if any, or -1 otherwise.
 *
 * \param[in] s the string to be searched for c
 * \param[in] c the characther to search for
 *
 * \returns the index of the first occurence of \c c in \c s, or -1 if
 *          \c c is not contained in \c s.
 */
int str_indexof(const char* s, char c);

/**
 * Returns the index of the last occurence of 'c' in 's' - if any, or -1 otherwise
 */
int str_rindexof(const char* s, char c);

/**
 * Same as the standard sprintf(), except that it will allocate the right size buffer
 * automatically. The caller should call free_char() on the returned string. The
 * current implementation requires C99 compiler support.
 */
char* asprintf( const char* format, ... );

/**
 * Frees a string array (free all strings then free the array).
 **/
void free_str_array( char** array, int array_size );

/** Creates a vector of strings from a array. The list is NULL ended */
char** array_to_null_ended_str_list(array list);

/** Frees a NULL ended vector of strings */
void free_null_ended_str_list(char** list);

/**
 * Converts a int32 value to string.
 * The caller should call free on the returned string.
 */
char* int32_to_str(int32_t v);

/**
 * Converts a int64 value to string.
 * The caller should call free on the returned string.
 */
char* int64_to_str(int64_t v);

/**
 * Converts a double value to string.
 * The caller should call free on the returned string.
 */
char* double_to_str(double v);


/* simple object pool */
typedef void* objpool;

objpool create_objpool( int max_count, process_value_fn destroy_obj_fn );
void destroy_objpool( objpool opool );

void* objpool_get( objpool opool );
void  objpool_put( objpool opool, void* obj );

int   objpool_count    ( const objpool opool );
int   objpool_max_count( const objpool opool );

void  objpool_dump_stats( const objpool opool, const char* msg );

int sqlite3_table_exists(sqlite3 *db, const char* name);
/* string buffer */
typedef void* strbuf;

/**
 * Creates a new string buffer.
 */
strbuf create_strbuf( void );

/**
 * Creates a new string buffer.
 * Uses the specified buffer capacity increment.
 */
strbuf create_strbuf_args( int grow_by );

/**
 * Destroys the specified string buffer.
 */
void destroy_strbuf( strbuf sb );

/**
 * Appends the 'str' string to the string buffer.
 */
void sb_append( strbuf sb, const char* str );

void sb_append_char( strbuf _sb, char c );

/**
 * Appends to the string buffer according to the specified format.
 * The format complies with the format of the *printf functions.
 * The current implementation requires C99 compiler support.
 */
void sb_appendf( strbuf sb, const char* format, ... );

/**
 * Appends at most 'n' chars from the 'str' string to the string buffer.
 */
void sb_appendn( strbuf sb, const char* str, int n );

/**
 * Inserts he 'str' string to the string buffer at the specified index.
 */
void sb_insert( strbuf sb, int index, const char* str );

/**
 * Inserts to the string buffer at the specified index according to the
 * specified format. The format complies with the format of the *printf
 * functions. The current implementation requires C99 compiler support.
 */
void sb_insertf( strbuf sb, int index, const char* format, ... );

/**
 * Insert at most 'n' chars from the 'str' string to the string buffer
 * at the specified index.
 */
void sb_insertn( strbuf sb, int index, const char* str, int n );

/**
 * Removed from the string buffer 'length' chars, starting from 'index' index.
 */
void sb_remove( strbuf sb, int index, int length );

/**
 * Returns the length of the string buffer.
 */
int sb_length( const strbuf sb );

/**
 * Gets a copy of the string from the string buffer if 'str' is not NULL.
 * Gets the length of the string from the string buffer if 'len' is not NULL.
 * This call does not modify the content of the string buffer.
 */
char* sb_copy( const strbuf sb );

/**
 * Returns a pointer to the (const) string buffer. Use this pointer only till
 * the first modification of the string buffer.
 */
const char* sb_get( const strbuf sb );

/**
 * Release the string from the string buffer.
 * Gets the length of the string from the string buffer if 'len' is not NULL.
 * This call will clear the content of the string buffer.
 */
char* sb_release( strbuf sb );

/**
 * Clears the content of the specified string buffer.
 */
void sb_clear( strbuf sb );

/**
 * Returns the capacity of the string buffer.
 */
int sb_cap( const strbuf sb );

/**
 * Ensures that the capacity of the string buffer is at least the specified value.
 */
void sb_ensure_cap( strbuf sb, int capacity );


/**
 * Returns the backslash-style encoded version of the 'str' string.
 * The characters encoded are '\\' and all chars from the 'encode_chars' string.
 * If 'newlines' flag is non-zero then '\\n' and '\\r' chars are also encoded to "\\n" and "\\r".
 * Eg:
 *   backslash_encode("\\ xy", "x", 0 ) returns "\\\\ \\xy"
 *   backslash_encode("a\\b 0 1 \r\n", "01", 1 ) returns "a\\\\b \\0 \\1 \\r\\n"
 *   backslash_encode("a\\b 0 1 \r\n", "01", 0 ) returns "a\\\\b \\0 \\1 \\r\\n"
 * The caller should call free on the returned string.
 */
char* backslash_encode( const char* str, const char* encode_chars, int newlines );

/**
 * Returns the backslash-style decoded version of the 'str' string.
 * The characters decoded are '\\' and all chars from the 'encode_chars' string.
 * If 'newlines' flag is non-zero then \\n" and "\\r" chars are also decoded to '\\n' and '\\r'.
 * Eg:
 *   backslash_decode("\\\\ \\xy", "x", 0 ) returns "\\ xy"
 *   backslash_decode("a\\\\b \\0 \\1 \\r\\n", "01", 1 ) returns "a\\b 0 1 \\r\\n"
 *   backslash_decode("a\\\\b \\0 \\1 \r\n", "01", 0 ) returns "a\\b 0 1 \\r\\n"
 * The caller should call free on the returned string.
 */
char* backslash_decode( const char* str, const char* encode_chars, int newlines );

/**
 * Returns a pointer to the first non backslash-encoded occurrence of the 'ch' character
 * in the 'str' string. If no match if found the function returns NULL.
 * Eg:
 *   strchr_encoded("\\aa", 'a') returns a pointer to the second 'a' character.
 *   strchr_encoded("\\a ", 'a') returns NULL.
 */
char* strchr_encoded( const char* str, char ch );


/**
 * Creates a string representing the content of a named string hash and append it
 * to the string pointed by 'str'.
 * eg: "hash-name {key1=val1 key2=val2}"
 */
int ht_to_str( se_hash ht, const char* name, char** str, int* len );

/**
 * Parses a string obtained by calling the ht_to_str() function and
 * returns a hash with the coresponding values.
 * On error the function returns NULL.
 * eg: "hash-name-1 {key1=val1 key2=val2} hash-name-2 {key1=val1 key2=val2}"
 */
se_hash ht_from_str( const char* str );

/**
 * Destroys an hash obtained by calling the ht_from_str() function.
 */
void destroy_htht( se_hash ht );


/**
 * Returns the name of the current UNIX user.
 * The user name is obtained by querying the "USER" environment variable.
 * If the "USER" variable is not defined the "<unkn-user>" string is returned.
 * The caller should call free on the returned string.
 */
char* get_user(void);

/**
 * Returns the name of the current host.
 * On error the "localhost" string is returned.
 * The caller should call free on the returned string.
 */
char* get_hostname(void);

/** force a hostname - just for this app - useful for debugging and testing */
void set_hostname(const char* hn);

/**
 * Returns a string representing the current local datetime.
 * The caller should call free on the returned string.
 */
char* get_datetime(void);

char* get_adatetime(time_t t);

/**
 * Returns the current working directory.
 * The caller should call free on the returned string.
 */
char* get_cwd(void);

/**
 * Returns the user or system defined temp directory.
 *
 * It never returns NULL, and the return value needs to be
 * de-allocated with free.
 *
 * Rules:
 * 1. If the env. variable TMPDIR is defined, then this is used.
 * 2. If the env. variable TEMPDIR is defined, then this is used.
 * 3. If P_tmpdir is defined (usually in stdio.h), then this is used.
 * 4. /tmp is used.
 */
char* get_tmp_dir();


/**
 * Generates a random (type 4) Universal Unique IDentifier.
 * Use this function only after init() was called.
 * The caller should call free on the returned string.
 */
char* generate_uuid(void);


/**
 * Determines if the process with the specified pid exists.
 */
int process_exists(int pid, int* exists);

/**
 * Convenient wrapper for the "system" function.
 */
int se_system(const char* command, ...);

/**
 * Conveniently formats a number of bytes (eg 2*1024 => "2KB" )
 * The caller should call free_char on the returned string.
 */
char* format_nbytes( off_t bytes );

/**
 * Conveniently formats a number of seconds (eg 70.2 => "1m10.2s" )
 * The caller should call free_char on the returned string.
 */
char* format_secs( double total_sec );


void ebcdic2ascii( char *s, size_t len );
void ascii2ebcdic( char *s, size_t len );

char* center_string(const char* str, char fill, int n, int truncate);


/**
 * Some systems have stricmp().  Others have strcasecmp().  Because
 * there is no consistency, we will define our own.
 */
int stricmp(const char *a, const char *b);
int strnicmp(const char *a, const char *b, size_t _n);


char* human_readable_size(int64_t size);
void print_human_readable_size(FILE* f, int64_t size);



size_t remove_nans(float* a, size_t n, int copy_last, float replacement);

#endif