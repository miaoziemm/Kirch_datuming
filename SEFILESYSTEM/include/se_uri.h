#ifndef SE_URI_H
#define SE_URI_H
#include <SEBASIC/include/se_basic.h>
#include "se_io.h"
#include "se_fs_io.h"



/**
 * Returns the base path (ending in a slash) of the specified path
 * eg: "./home/je/ss.conf" > "./home/je/"
 **/
char* se_get_base_path( const char* path );


/**
 * Converts the specified path to an absolute path.
 * If 'path' is a relative path (does not start with a slash), then it is
 * considered as relative to 'base_path'.
 *   eg: "../d3/./g.xml" + "/d1/d2/" => "/d1/d3/g.xml"
 * If 'path' is an absolute path (starts with a slash), then the 'base_path'
 * parameter is ignored.
 *   eg: "/d1/d2/../d3/./" + "/" => "/d1/d3/"
 * The 'base_path' parameter must point to an absolute path.
 **/
char* se_make_path_absolute( const char* path, const char* base_path );


/**
 * Returns true (non-zero) if the specified path is an absolute path,
 * or zero otherwise
 **/
inline int se_is_path_relative( const char* path )
{
    ASSERT(path);
    return path[0] != '/' && path[0] != '~';
}


/**
 * Returns true (non-zero) if the specified path is a relative path,
 * or zero otherwise
 **/
inline int se_is_path_absolute( const char* path )
{
    ASSERT(path);
    return path[0] == '/' || path[0] == '~';
}


/**
 * Returns the filename of the specified path
 * eg: "./home/je/ss.conf" > "ss.conf"
 **/
inline char* se_get_filename( const char* path )
{
    const char* last_slash;

    ASSERT(path);

    last_slash = strrchr( path, '/' );
    if( ! last_slash )
        return se_strdup( path );
    else
        return se_strdup( last_slash + 1 );
}


/**
 * Returns the normalized path of the specified path.
 * eg: "./d1/d2/./../file" > "d1/file"
 **/
char* se_normalize_path( const char* path );


/**
 * Read all data from the specified stream and return it as string
 **/
char* se_read_text_stream(FILE* f);

/**
 * Read an entire line from the specifed stream. The callers needs to
 * call se_free on the returned pointer.
 *
 * If max_len is not zero and after max_len characters no EOL or EOF
 * is found, then NULL is returned (it means that there is no real
 * text line in the file).
 *
 * Return NULL if EOF or error - an warning will be printed if error.
 *
 */
char* se_read_text_stream_line(FILE* f, size_t max_len);

/**
 * Read an entire line from the specifed stream. The callers needs to
 * call se_free on the returned pointer.
 *
 * Return NULL if EOF or error - an warning will be printed if error.
 */
inline char* se_read_text_stream_full_line(FILE* f)
{
    return se_read_text_stream_line(f, 0);
}

/**
 * Read all data from the specified file and return it as string
 **/
char* se_read_text_file( const char* filename );

/**
 * Read all input from stdin and return it as string
 **/
inline char* se_read_stdin_file( void )
{
    return se_read_text_stream(stdin);
}



/**
 * Determines if the specified directory is empty or not.
 * Returns se_OK on success, on fail a WARN is generated and se_ERROR
 * is returned.
 */
int se_dir_is_empty( const char* dir, int* is_empty );


/**
 * Returns the absolute path version of filename, resolved relative to
 * a peer.
 *
 * If the peer path is relative, then it is resolved first against the
 * current working directory.
 *
 * This function is useful for example to get the absolute path for
 * the "in" parameter in a .H file.
 *
 * It always returns a newely allocated array, even if the path is
 * absolute - so the caller is always responsible to free the returned
 * value.
 */
char* se_abs_path_from_peer( const char* filename, const char* peer );

/**
 * Returns true (non-zero) if the two filenames point to the same file,
 * of false (zero) otherwise
 */
int se_is_same_file( const char* filename1, const char* filename2 );


/**
 * Recursive mkdir-like function
 * eg: se_mkdir("./dir1/dir2", 0755) or se_mkdir("/dir1/dir2", 0755)
 * Returns se_OK on success, or se_ERROR on error.
 */
int se_mkdir( const char* path, mode_t mode );

/**
 * Lists the contents of the specified dir (without "." and "..").
 * Returns se_OK on success, or se_ERROR on error.
 */
int se_ls( const char* dir, array entries );

#define se_PATHTREE_DEFAULT            (0x00000001)
#define SE_PATHTREE_TOP_BOTTOM         (0x00000001)
#define SE_PATHTREE_BOTTOM_TOP         (0x00000002)
#define SE_PATHTREE_IGNORE_LS_ERRORS   (0x00000004)

/**
 * Recursively traverse a directory sub-tree, top-to-bottom or bottom-to-top,
 * The 'process' function is called for each file or dir from the sub-tree,
 * and is expected to return se_OK on success or se_ERROR if there is no need
 * to traverse the rest of the sub-tree.
 * Returns se_OK on success, or se_ERROR on error.
 */
int se_traverse_pathtree( const char* dir, process_path_fn process, void* arg, uint32_t flags );

/**
 * This function "rm -rf" the specified sub-tree.
 * Paths that cannot be removed will cause a WARN.
 * Returns se_OK only.
 */
int se_rmrf( const char* dir );

/**
 * Finds all the occurences of a file starting from a directory and
 * moving upwards the file hierarchy.
 *
 * \param[in] file_name - the name of the file to be searched.
 *
 * \param[in] start_point - start point for the search, which can
 *                          either be a directory or a file; if the
 *                          \c start_point is a relative path, then 
 *                          it is resolved to the current working directory.
 * 
 * \return a list of full paths to the occurences of the \c file_name,
 *         the closest to the \c start_point being at the beginnning
 *         of the list ant the closest to the root of the filesystem
 *         to the end of the list.
 */
array se_find_file_upwards( const char* file_name, const char* start_point );


/**
 * Returns the home directory for the current user.
 * On error a WARN is generated and the "/" value is returned.
 * The caller must free the returned string.
 *
 * The function first checks for the environment variable HOME, which
 * is usually set by the login programs. If this is not found, then it
 * tries to use the record from the password database.
 */
char* se_get_home_directory( void );


/**
 * Reads the content of a symbolic link.
 *
 * If something goes wrong, the function prints a warning message and
 * returns a duplicate of the path argument.
 *
 * The function allways returns a non-NULL value, and the caller is
 * responsible to call se_free on it.
 */
char* se_readlink(const char* path);

/**
 * Returns the path to the current executable. The argv0 argument can
 * be NULL or the first element of the command line arguments. Any
 * other value may cause undexpected results.
 *
 * If the program was invoked using full path to the executable, then
 * that path is returned, without any other processing. Otherwise, on
 * Linux it reads the value of the /proc/self/exe link.
 *
 * The function allways returns a non-NULL value, and the caller is
 * responsible to call se_free on it.
 */
char* se_exe_path(const char* argv0);


char* se_replace_extension(const char* file, const char* rep);

const char* se_get_file_extension(const char* file);

int se_file_exists( const char* path );

int64_t se_file_length( const char* path );





/** Public structure / handle. */
typedef void* se_uri;


/**
 * Creates an empty URI
 **/
se_uri se_create_uri( void );

/**
 * Creates a copy of the specified URI
 **/
se_uri se_create_uri_copy( const se_uri uri );

/**
 * Destroy an the specified URI
 **/
void se_destroy_uri( se_uri uri );


/**
 * Parse the specified URI.
 * Returns a se_uri structure on success or NULL on failure.
 **/
se_uri se_parse_uri( const char* str, const char* default_scheme );

/**
 * Format the string corresponding to the specified URI
 * Use the normalize flag to get an normalized URI.
 **/
char* se_format_uri( const se_uri uri, int normalize );

int se_uri_is_absolute( const se_uri uri );
int se_uri_is_relative( const se_uri uri );
int se_uri_is_url( const se_uri uri );
int se_uri_is_urn( const se_uri uri );
const char* se_uri_get_scheme( const se_uri uri );
const char* se_uri_get_scheme_specific( const se_uri uri );
const char* se_uri_get_user( const se_uri uri );
const char* se_uri_get_pass( const se_uri uri );
const char* se_uri_get_host( const se_uri uri );
int se_uri_get_port( const se_uri uri );
int se_uri_get_strong_port( const se_uri uri );
const char* se_uri_get_path( const se_uri uri );
const char* se_uri_get_obj( const se_uri uri );
const char* se_uri_get_ext( const se_uri uri );
const char* se_uri_get_query( const se_uri uri );
const char* se_uri_get_frag( const se_uri uri );
const char* se_uri_get_user_info( const se_uri uri );
const char* se_uri_get_host_port( const se_uri uri );
const char* se_uri_get_host_info( const se_uri uri );
const char* se_uri_get_path_obj( const se_uri uri );
const char* se_uri_get_path_obj_query( const se_uri uri );
se_hash se_uri_get_query_params( const se_uri uri );


/**
 * Returns the default port for the specified scheme, or -1 if the scheme
 * is unknown. The scheme should be lower case
 **/
int se_uri_default_port( const char* scheme );


/**
 * Compares the specified normalized URIs.
 * The path & object (/dir/x.xml) are compared case-sensitive.
 * Return true (non-zero) if the URIs are equals, or zero otherwise
 **/
int se_uri_equals( const se_uri norm_uri1, const se_uri norm_uri2 );


/**
 * Normalizes the specified URI
 * eg: "http://site.com/dir/dir1/../index.php" > "http://site.com/dir/index.php"
 **/
void se_normalize_uri( se_uri uri );

/**
 * Creates a base URI from the specified absolute URI
 * eg: "http://site.com/dir/index.php" > "http://site.com/dir/"
 **/
se_uri se_base_uri( const se_uri abs_uri );


/**
 * Creates a absolute URI from the specified relative & base URIs
 * eg: "../dir2/index.cgi" + "http://site.com/dir/dir1/" > "http://site.com/dir/dir2/index.cgi"
 **/
se_uri se_uri_rel2abs( const se_uri rel_uri, const se_uri base_uri );

/**
 * Creates a relative URI from the specified absolute & base URIs
 * eg: "http://site.com/dir/dir2/index.cgi" + "http://site.com/dir/dir1/" > "../dir2/index.cgi"
 **/
se_uri se_uri_abs2rel( const se_uri abs_uri, const se_uri base_uri );


#endif