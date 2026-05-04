#include "../include/se_uri.h"



/**
 * Returns the normalized path of the specified path.
 * eg: "./d1/d2/./../file" > "d1/file"
 **/
char* se_normalize_path( const char* path )
{
    char* norm_path;
    int len;

    ASSERT( path );

    norm_path = se_make_path_absolute( path, "/" );

    /* remove the leading '/' for relative paths */
    if( se_is_path_relative(path) ) {
        len = strlen(norm_path);
        memmove( norm_path, norm_path + 1, len * sizeof(char) );
    }

    return norm_path;
}


/**
 * Returns the base path (ending in a slash) of the specified path
 * eg: "./home/je/ss.conf" > "./home/je/"
 **/
char* se_get_base_path( const char* path )
{
    const char* last_slash;

    ASSERT(path);

    last_slash = strrchr( path, '/' );
    if( ! last_slash ) {
        return se_strdup("./");
    } else {
        if( last_slash[1] == 0 && last_slash != path) {
            do {
                --last_slash;
            } while(last_slash > path && *last_slash == '/');

            while( last_slash > path && *last_slash != '/' )
                --last_slash;
        }

        return se_strndup( path, last_slash - path + 1 );
    }
}


/**
 * Converts the specified path to an absolute path.
 * If 'path' is a relative path (does not start with a slash), then it is
 * considered as relative to 'base_path'.
 *   eg: "../d3/./g.xml" + "/d1/d2/" => "/d1/d3/g.xml"
 * If 'path' is an absolute path (starts with a slash), then the 'base_path'
 * parameter is ignored.
 *   eg: "/d1/d2/../d3/./" + "/" => "/d1/d3/"
 *
 * The base_path parameter is assumed to be a directory.
 **/
char* se_make_path_absolute( const char* path, const char* _base_path )
{
    array path_array, abs_array;
    char *abs_path, *item;
    int i, size, home_dir = 0;
    char* base_path;

    if( _base_path == NULL || (_base_path[0] == '.' && _base_path[1] == 0) ) {
        base_path = get_cwd();
    } else {
        if( se_is_path_relative( _base_path ) ) {
            base_path = se_make_path_absolute( _base_path, NULL );
        } else {
            base_path = se_strdup(_base_path);
        }
    }

    ASSERT(path);

    /* split path in '/'-separated items */
    path_array = split2( path, '/', SPLIT_NO_EMPTYFIELDS );

    if( se_is_path_relative( path )) {
        /* path is relative to base_path */
        /* split base_path in '/'-separated items */
        abs_array = split2( base_path, '/', SPLIT_NO_EMPTYFIELDS );
    } else {
        /* path is not relative to base_path, ignore base_path parameter */
        abs_array = create_array();
    }

    /* create the abs_array items */
    for (i = 0; i < array_size(path_array); ++i) {
        item = (char*) array_get_at( path_array, i );

        if(i == 0 && 0 == strcmp( item, "~" )) {
            /* substitute "~" with the home dir */
            array_add( abs_array, se_get_home_directory() );
            home_dir = 1;
            continue;
        }

        if(0 == strcmp( item, "." )) {
            /* ignore single dots */
            continue;
        }

        if(0 == strcmp( item, ".." )) {
            /* remove last item from abs_array (if any) */
            size = array_size( abs_array );
            if( size ) {
                item = (char*)array_remove_at( abs_array, size - 1 );
                free1char( item );
            }
        } else {
            /* add the item at the end of abs_array */
            array_add( abs_array, se_strdup(item) );
        }
    }

    /* now create abs_path by merging all items from abs_array */

    /* compute the size of abs_path */
    size = 0;
    for (i = 0; i < array_size(abs_array); ++i) {
        item = (char*) array_get_at( abs_array, i );
        size += strlen(item);
    }
    size += array_size( abs_array ) + 1;

    /* alloc abs_path */
    abs_path = alloc1char( size + 1 );

    /* merge all items from abs_array */
    if(! home_dir) {
        strcpy( abs_path, "/" );
    }
    size = array_size(abs_array);
    for (i = 0; i < size; ++i) {
        item = (char*) array_get_at( abs_array, i );
        strcat( abs_path, item );
        if( i < size - 1 )
            strcat( abs_path, "/" );
    }
    if( str_endswith( path, "/" ) && ! str_endswith( abs_path, "/" ))
        strcat( abs_path, "/" );

    /* cleanup time */
    destroy_array( path_array, 1 );
    destroy_array( abs_array, 1 );

    free(base_path);

    /* all done */
    return abs_path;
}

/**
 * Read all data from the specified file and return it as string
 **/
char* se_read_text_file(const char* filename)
{
    char* text = NULL;
    FILE* f = fopen(filename, "r");
    if(f == NULL) {

        LOG_CONSOLE(SE_ERROR, "Cannot open file %s", filename);
        return NULL;
    }
    text = se_read_text_stream(f);
    fclose(f);
    return text;
}

/**
 * Read all data from the specified stream and return it as string
 **/
char* se_read_text_stream(FILE* f)
{
    char *text = NULL, *buf;
    int pos, len, cnt, newlen;
    const int BUFF_GROW_BY = 1024;

    for (pos = len = 0; ! feof(f); ) {

        if(pos >= len-1) {
            newlen = pos + BUFF_GROW_BY + 1;
            buf = alloc1char(newlen);
            if(text) {
                memcpy(buf, text, len);
                free(text);
            }
            text = buf;
            len = newlen;
        }

        cnt = fread( text+pos, 1, len-pos-1, f );
        if( ferror(f) ) {
         
            LOG_CONSOLE(SE_ERROR, "Problems reading from a stream");
            return NULL;
        }

        if(cnt > 0) {
            pos += cnt;
            *(text+pos) = 0;
        }
    }

    text[len-1] = 0;
    return text;
}

char* se_read_text_stream_line(FILE* f, size_t max_len)
{
    char* line = NULL;
    size_t n, len;

    for(len = 0, n = 200;; n+= 64) {
        line = realloc1char(line, n+1);
        line[n] = 0;

        if( fgets(line+len, n-len, f) == NULL ) {
            int err = ferror(f);
            if( feof(f) ) {
                free(line);
                return NULL;
            }
            if( err ) {
                LOG_CONSOLE(SE_WARNING, "I/O error reading a full line: code=%d", err);
                free(line);
                return NULL;
            }
        }

        if( feof(f) ) break;

        len = strlen(line);
        if( len == 0 ) { /* non text file ? */
            free(line);
            return NULL;
        }

        if(line[len - 1] == '\n')
            break;

        if(max_len > 0 && len >= max_len) {
            free(line);
            return NULL;
        }
    }

    return line;
}

int se_file_exists( const char* path )
{
    int fd;
    ASSERT(path);

    if( se_io_is_directory( path ) ) {
        return 0;
    }

    fd = open(path, O_RDONLY);
    if (fd >= 0) {
        VERIFY(!close(fd));
        return 1;
    }
    return 0;
}

int64_t se_file_length( const char* path )
{
    ASSERT(path);

    if(path) {
        struct stat st;
        if( 0 == stat(path, &st) ) {
            return (int64_t)st.st_size;
        } else {
            int err = errno;
            if(err != ENOENT && err != ENAMETOOLONG) {
                LOG_CONSOLE(SE_WARNING, "Cannot stat %s: %d - %s", path, err, strerror(err));
            }
        }
    }

    return 0;
}

/**
 * Determines if the specified directory is empty or not.
 * Returns CODE_SUCCESS on success, on fail a WARN is generated and CODE_ERROR
 * is returned.
 */
// int se_dir_is_empty( const char* dir, int* is_empty )
// {
//     DIR *d;
//     struct dirent dentry, *pdentry;
//     int err, empty;

//     ASSERT(dir);
//     ASSERT(is_empty);

//     *is_empty = 0;

//     /* open the specified directory */
//     d = opendir(dir);
//     if( !d ) {
//         err = errno;
//         LOG_CONSOLE(SE_WARNING, "Cannot open directory %s: %d - %s", dir, err, strerror(err));
//         return CODE_ERROR;
//     }

//     /* for each directory entry */
//     empty = 1;
//     for(;;) {
//         /* get next directory entry */
//         err = readdir_r( d, &dentry, &pdentry );

//         /* check for an error */
//         if (err != 0) {
//             LOG_CONSOLE(SE_WARNING, "readdir_r() failed for dir '%s': errno=%d - %s", dir, err, strerror(err));
//             closedir(d);
//             return CODE_ERROR;
//         }

//         /* check for no more dir entries */
//         if (pdentry == NULL) {
//             break;
//         }

//         /* check for a directory entry other than "." and ".." */
//         if(dentry.d_name[0] == '.') {
//             if( ( dentry.d_name[1] == 0 ) ||
//                 ( dentry.d_name[1] == '.' && dentry.d_name[2] == 0 ) ) {
//                 continue;
//             }
//         }
//         empty = 0;
//         break;
//     }

//     *is_empty = empty;

//     /* close the directory */
//     if(closedir(d)) {
//         err = errno;
//         LOG_CONSOLE(SE_WARNING, "closedir() failed for dir '%s': errno=%d - %s", dir, err, strerror(err));
//         return CODE_ERROR;
//     }

//     return CODE_SUCCESS;
// }

// 修改 se_dir_is_empty 函数 (大约在326行附近)
int se_dir_is_empty( const char* dir, int* is_empty )
{
    DIR *d;
    struct dirent *dentry;
    int err, empty;

    ASSERT(dir);
    ASSERT(is_empty);

    *is_empty = 0;

    /* open the specified directory */
    d = opendir(dir);
    if( !d ) {
        err = errno;
        LOG_CONSOLE(SE_WARNING, "Cannot open directory %s: %d - %s", dir, err, strerror(err));
        return CODE_ERROR;
    }

    /* for each directory entry */
    empty = 1;
    errno = 0;
    while ((dentry = readdir(d)) != NULL) {
        /* check for a directory entry other than "." and ".." */
        if(dentry->d_name[0] == '.') {
            if( ( dentry->d_name[1] == 0 ) ||
                ( dentry->d_name[1] == '.' && dentry->d_name[2] == 0 ) ) {
                continue;
            }
        }
        empty = 0;
        break;
    }

    /* check for an error */
    if (errno != 0) {
        err = errno;
        LOG_CONSOLE(SE_WARNING, "readdir() failed for dir '%s': errno=%d - %s", dir, err, strerror(err));
        closedir(d);
        return CODE_ERROR;
    }

    *is_empty = empty;

    /* close the directory */
    if(closedir(d)) {
        err = errno;
        LOG_CONSOLE(SE_WARNING, "closedir() failed for dir '%s': errno=%d - %s", dir, err, strerror(err));
        return CODE_ERROR;
    }

    return CODE_SUCCESS;
}



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
char* se_abs_path_from_peer(const char* filename, const char* peer)
{
    ASSERT(filename);
    ASSERT(peer);

    if( se_is_path_relative( filename ) ) {
        char* ret;
        char* base = se_get_base_path(peer);
        if( se_is_path_relative( base ) ) {
            char* cwd = get_cwd();
            char* abase = se_make_path_absolute( base, cwd );
            ret = se_make_path_absolute( filename, abase );
            free(abase);
            free(cwd);
        } else {
            ret = se_make_path_absolute( filename, base );
        }

        free(base);
        return ret;
    }
    else {
        return se_strdup(filename);
    }
}

/**
 * Returns true (non-zero) if the two filenames point to the same file,
 * of false (zero) otherwise
 */
int se_is_same_file(const char* filename1, const char* filename2)
{
    ASSERT(filename1);
    ASSERT(filename2);
    /* needs to be replaced .... */
    return ! strcmp( filename1, filename2 );
}

/**
 * Recursive mkdir-like function
 * eg: se_mkdir("./dir1/dir2", 0755) or se_mkdir("/dir1/dir2", 0755)
 * Returns CODE_SUCCESS on success, or CODE_ERROR on error.
 */
int se_mkdir( const char* path, mode_t mode )
{
    char *path_norm, *p;
    int rc;
    ASSERT(path);

    path_norm = se_normalize_path(path);
    rc = CODE_SUCCESS;

    for( p = path_norm; p; ) {
        p = strchr(p + 1, '/');
        if(p) {
            *p = '\0';
        }
        if( mkdir(path_norm, mode) < 0 && errno != EEXIST ) {
            LOG_CONSOLE(SE_WARNING, "Cannot create directory '%s': errno=%d - %s", path_norm, errno, strerror(errno));
            rc = CODE_ERROR;
            break;
        }
        if (p) {
            *p = '/';
        }
    }

    free1char(path_norm);
    return rc;
}

/**
 * Lists the contents of the specified dir (without "." and "..").
 * Returns CODE_SUCCESS on success, or CODE_ERROR on error.
 */
// int se_ls( const char* dir, array entries )
// {
//     DIR *d;
//     struct dirent dentry, *pdentry;
//     int err;

//     ASSERT(dir);
//     ASSERT(entries);

//     if (array_size(entries)) {
//         array_clear(entries, 1);
//     }

//     /* open the specified directory */
//     d = opendir(dir);
//     if( !d ) {
//         err = errno;
     
//         LOG_CONSOLE(SE_WARNING, "Cannot open directory %s: %d - %s", dir, err, strerror(err));
//         return CODE_ERROR;
//     }

//     /* for each directory entry */
//     for(;;) {
//         char* entry = NULL;

//         /* get next directory entry */
//         err = readdir_r( d, &dentry, &pdentry );

//         /* check for an error */
//         if (err != 0) {
//             LOG_CONSOLE(SE_WARNING, "readdir_r() failed for dir '%s': errno=%d - %s", dir, err, strerror(err));
//             closedir(d);
//             array_clear(entries, 1);
//             return CODE_ERROR;
//         }

//         /* check for no more dir entries */
//         if (pdentry == NULL) {
//             break;
//         }

//         /* check for a directory entry other than "." and ".." */
//         if (strcmp(dentry.d_name, ".") == 0 || strcmp(dentry.d_name, "..") == 0) {
//             continue;
//         }

//         entry = se_strdup(dentry.d_name);
//         array_add(entries, entry);
//     }

//     /* close the directory */
//     if(closedir(d)) {
//         err = errno;
//         LOG_CONSOLE(SE_WARNING, "closedir() failed for dir '%s': errno=%d - %s", dir, err, strerror(err));
//         return CODE_ERROR;
//     }

//     return CODE_SUCCESS;
// }

// 修改 se_ls 函数 (大约在479行附近)
int se_ls( const char* dir, array entries )
{
    DIR *d;
    struct dirent *dentry;
    int err;

    ASSERT(dir);
    ASSERT(entries);

    if (array_size(entries)) {
        array_clear(entries, 1);
    }

    /* open the specified directory */
    d = opendir(dir);
    if( !d ) {
        err = errno;
        LOG_CONSOLE(SE_WARNING, "Cannot open directory %s: %d - %s", dir, err, strerror(err));
        return CODE_ERROR;
    }

    /* for each directory entry */
    errno = 0;
    while ((dentry = readdir(d)) != NULL) {
        char* entry = NULL;

        /* check for a directory entry other than "." and ".." */
        if (strcmp(dentry->d_name, ".") == 0 || strcmp(dentry->d_name, "..") == 0) {
            continue;
        }

        entry = se_strdup(dentry->d_name);
        array_add(entries, entry);
    }

    /* check for an error */
    if (errno != 0) {
        err = errno;
        LOG_CONSOLE(SE_WARNING, "readdir() failed for dir '%s': errno=%d - %s", dir, err, strerror(err));
        closedir(d);
        array_clear(entries, 1);
        return CODE_ERROR;
    }

    /* close the directory */
    if(closedir(d)) {
        err = errno;
        LOG_CONSOLE(SE_WARNING, "closedir() failed for dir '%s': errno=%d - %s", dir, err, strerror(err));
        return CODE_ERROR;
    }

    return CODE_SUCCESS;
}


static int se_traverse_pathtree_internal( const char* path, int is_dir,
                                           process_path_fn process, void* arg,
                                           uint32_t flags )
{
    int rc = CODE_SUCCESS;

    /* process the path */
    if (flags & SE_PATHTREE_TOP_BOTTOM) {
        rc = process(path, is_dir, arg);
    }

    if( rc == CODE_SUCCESS && is_dir ) {
        array entries = create_array();
        int i;

        /* list the entries of the directory */
        if( se_ls(path, entries) != CODE_SUCCESS ) {
            if (! (flags & SE_PATHTREE_IGNORE_LS_ERRORS)) {
                rc = CODE_ERROR;
            }
        }

        /* process all entries */
        for (i = 0; rc == CODE_SUCCESS && i < array_size(entries); ++i) {
            const char* name = (const char*) array_get_at(entries, i);
            const char* entry = asprintf("%s/%s", path, name);

            int entry_is_dir = se_io_is_directory(entry);
            rc = se_traverse_pathtree_internal( entry, entry_is_dir, process, arg, flags );

            free1char((void*)entry);
        }

        destroy_array(entries, 1);
    }

    /* process the path */
    if ((flags & SE_PATHTREE_BOTTOM_TOP) && rc == CODE_SUCCESS) {
        rc = process(path, is_dir, arg);
    }

    return rc;
}

/**
 * Recursively process a directory sub-tree.
 * The 'process' function is called for each file or dir from the sub-tree,
 * and is expected to return CODE_SUCCESS on success or CODE_ERROR if there is no need
 * to process the rest of the sub-tree.
 * Returns CODE_SUCCESS on success, or CODE_ERROR on error.
 */
int se_traverse_pathtree( const char* path, process_path_fn process, void* arg, uint32_t flags )
{
    /* check the flags */
    if (flags == 0) {
        flags = se_PATHTREE_DEFAULT;
    }
    if (! (flags & (SE_PATHTREE_TOP_BOTTOM | SE_PATHTREE_BOTTOM_TOP))) {
        ERROR(("se_traverse_pathtree(): invalid flags: SE_PATHTREE_TOP_BOTTOM or SE_PATHTREE_BOTTOM_TOP is required"));
    }
    if ((flags & SE_PATHTREE_TOP_BOTTOM) && (flags & SE_PATHTREE_BOTTOM_TOP)) {
        ERROR(("se_traverse_pathtree(): invalid flags: SE_PATHTREE_TOP_BOTTOM and SE_PATHTREE_BOTTOM_TOP are exclusive"));
    }

    int is_dir = se_io_is_directory(path);
    return se_traverse_pathtree_internal(path, is_dir, process, arg, flags);
}

static int se_rmrf_internal( const char* path, int is_dir, void* arg )
{
    int err;
    UNUSED(arg);

/*    WARN(("%s: %s", path, is_dir ? "dir" : "file")); */
    if (is_dir) {
        if (rmdir(path) != 0) {
            err = errno;
            LOG_CONSOLE(SE_WARNING, "se_rmrf_internal(): cannot remove dir '%s': errno=%d - %s",
                  path, err, strerror(err));
        }
    } else {
        if (unlink(path) != 0) {
            err = errno;
        
            LOG_CONSOLE(SE_WARNING, "se_rmrf_internal(): cannot remove file '%s': errno=%d - %s", 
                  path, err, strerror(err));
        }
    }

    return CODE_SUCCESS;
}

/**
 * This function "rm -rf" the specified sub-tree.
 * Paths that cannot be removed will cause a WARN.
 * Returns CODE_SUCCESS only.
 */
int se_rmrf( const char* dir )
{
    return se_traverse_pathtree( dir, se_rmrf_internal, NULL,
        SE_PATHTREE_BOTTOM_TOP | SE_PATHTREE_IGNORE_LS_ERRORS );
}

array se_find_file_upwards( const char* file, const char* start_point )
{
    array files;
    char *start_dir, *dir;
    int len;

    /* get the path if 'start_point' is not a dir */
    if (se_io_is_directory(start_point)) {
        start_dir = se_strdup(start_point);
    } else {
        start_dir = se_get_base_path(start_point);
    }

    /* work with absolute path */
    if( se_is_path_relative(start_dir) ) {
        dir = se_make_path_absolute( start_dir, NULL);
    } else {
        dir = se_normalize_path(start_dir);
    }

    /* remove trailing '/' */
    len = strlen(dir);
    if (len > 0 && dir[len - 1] == '/') {
        dir[len - 1] = '\0';
    }

    files = create_array();
    while (1) {
        /* format file's full path */
        char* f = asprintf("%s/%s", dir, file);
        int idx;

        /* check if the file exists */
        if (se_io_exists(f) && ! se_io_is_directory(f)) {
            array_add(files, f);
        } else {
            free(f);
        }

        /* move one dir up */
        len = strlen(dir);
        if (len == 0) {
            break;
        }
        idx = str_rindexof(dir, '/');
        if (idx < 0) {
            break;
        }
        dir[idx] = '\0';
    }

    free(start_dir);
    free(dir);

    return files;
}

/**
 * Returns the home directory for the current user.
 * On error a WARN is generated and the "/" value is returned.
 * The caller must free the returned string.
 *
 * The function first checks for the environment variable HOME, which
 * is usually set by the login programs. If this is not found, then it
 * tries to use the record from the password database.
 */
char* se_get_home_directory( void )
{
    char *home = getenv("HOME");
    if( home == NULL ) {
        struct passwd *pwent;
        errno = 0;
        pwent = getpwuid(getuid());
        if( pwent == NULL ) {
            int err = errno;
            LOG_CONSOLE(SE_WARNING, " Cannot obtain home path - getenv(\"HOME\") "
                  "and getpwuid failed: %d - %s",
                  err, strerror(err) );
        } else {
            home = pwent->pw_dir;
        }
    }

    if(home == NULL) return se_strdup("/");

    return se_strdup( home );
}

/**
 * Reads the content of a symbolic link.
 *
 * If something goes wrong, the function prints a warning message and
 * returns a duplicate of the path argument.
 *
 * The function allways returns a non-NULL value, and the caller is
 * responsible to call free on it.
 */
char* se_readlink(const char* path)
{
    size_t size = 256;
    char* buf = (char *)malloc(size * sizeof(char));
    memset(buf, 0, size * sizeof(char));
    ASSERT((path));
    for(;;) {
        ssize_t len = readlink(path, buf, size);
        if(len < 0) {
            int err = errno;
            free(buf);
            LOG_CONSOLE(SE_WARNING, "Cannot read the link %s: errno=%d - %s",
                  path, err, strerror(err));
            return se_strdup(path);
        } else if(len == (ssize_t)size) {
            free(buf);
            size += 64;
            buf =(char *)malloc(size * sizeof(char));
            memset(buf, 0, size * sizeof(char));
        } else {
            buf[len] = 0;
            return buf;
        }
    }

    /* never reaching this point */
    /*return se_strdup(path);*/
}

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
 * responsible to call free on it.
 */
char* se_exe_path(const char* argv0)
{
    /* if absolute, then we are lucky - just return the base dir */
    if( argv0 != NULL && se_is_path_absolute( argv0 ) ) {
        return se_strdup(argv0);
    } else {
#ifdef ARCH_LINUX
        return se_readlink("/proc/self/exe");
#else
#ifdef ARCH_OSX
	int ret;
	pid_t pid; 
	char pathbuf[PROC_PIDPATHINFO_MAXSIZE];

	pid = getpid();
	ret = proc_pidpath (pid, pathbuf, sizeof(pathbuf));
	if ( ret <= 0 ) {
	  WARN(("unable to determine exe path for %s",argv0));
	}
	return se_strdup(pathbuf);
#else
        /* Fallback: just return argv0 or an empty string if not available */
        if (argv0 != NULL) {
            return se_strdup(argv0);
        } else {
            return se_strdup("");
        }
#endif
#endif
    }
}




static void se_uri_split_in_two_left( const char* text,
                                           char delim,
                                           char** left,
                                           char** right )
{
    const char* s;
    ASSERT(text);
    ASSERT(left);
    ASSERT(right);

    s = strchr( text, delim );
    if (s) {
        *left  = se_strndup( text, s - text );
        *right = se_strdup( s + 1 );
    } else {
        *left  = se_strdup(text);
        *right = se_strdup("");
    }
}

static void se_uri_split_in_two_right( const char* text,
                                            char delim,
                                            char** left,
                                            char** right )
{
    const char* s;
    ASSERT(text);
    ASSERT(left);
    ASSERT(right);

    s = strchr( text, delim );
    if (s) {
        *left  = se_strndup( text, s - text );
        *right = se_strdup( s + 1 );
    } else {
        *left  = se_strdup("");
        *right = se_strdup(text);
    }
}

static char* se_uri_merge_two( const char* left,
                                    const char* right,
                                    char delim )
{
    char* s;
    int len_left;
    ASSERT(left);
    ASSERT(right);

    len_left = strlen(left);
    s = alloc1char( len_left + strlen(right) + 2 );

    strcpy( s, left );
    if( left[0] != '\0' && right[0] != '\0' && left[len_left - 1] != delim ) {
        strncat( s, &delim, 1 );
    }
    strcat( s, right );

    return s;
}


/** Internal URI structure */
typedef struct
{
    int   is_abs;          /* absolute/relative uri flag */
    int   is_url;          /* is URL flag */
    char* scheme;          /* protocol ("http", "file", etc) */
    char* scheme_specific; /* scheme-specific part */
    char* user;            /* username */
    char* pass;            /* password */
    char* host;            /* hostname or IP */
    int   port;            /* port no if present, or -1 */
    int   strong_port;     /* port no if present, default port for
                              scheme, or -1 */
    char* path;            /* path ending in slash, default "/" ("/faq/dir/") */
    char* obj;             /* object ("index.php") */
    char* ext;             /* object's extension ("php") */
    char* query;           /* everything between '?' and '#' */
    char* frag;            /* everything after '#' */
    char* user_info;       /* combined user[:pass] */
    char* host_port;       /* combined host[:port] */
    char* host_info;       /* combined [user[:pass]@]host[:port] */
    char* path_obj;        /* combined path[obj] */
    char* path_obj_query;  /* combined path[obj]?query */
    se_hash ht_query;     /* query params (key/value hashtable) */
} _uri_t;


/**
 * Creates an empty URI
 **/
se_uri se_create_uri( void )
{
    _uri_t* uri = (_uri_t*)malloc( sizeof(_uri_t) );

    uri->is_abs    = 1;
    uri->is_url    = 1;
    uri->scheme    = se_strdup("");
    uri->scheme_specific = se_strdup("");
    uri->user      = se_strdup("");
    uri->pass      = se_strdup("");
    uri->host      = se_strdup("");
    uri->port        = -1;
    uri->strong_port = -1;
    uri->path      = se_strdup("");
    uri->obj       = se_strdup("");
    uri->ext       = se_strdup("");
    uri->query     = se_strdup("");
    uri->frag      = se_strdup("");
    uri->user_info = se_strdup("");
    uri->host_port = se_strdup("");
    uri->host_info = se_strdup("");
    uri->path_obj  = se_strdup("");
    uri->path_obj_query = se_strdup("");
    uri->ht_query  = create_hash();

    return uri;
}

/**
 * Creates a copy of the specified URI
 **/
se_uri se_create_uri_copy( const se_uri _src_uri )
{
    const _uri_t* src_uri = (const _uri_t* ) _src_uri;
    _uri_t* uri;
    htiter hti;
    const char* q_key;
    void* q_val;
    ASSERT(src_uri);

    uri = (_uri_t*)malloc( sizeof(_uri_t) );

    uri->is_abs    = src_uri->is_abs;
    uri->is_url    = src_uri->is_url;
    uri->scheme    = se_strdup(src_uri->scheme);
    uri->scheme_specific = se_strdup(src_uri->scheme_specific);
    uri->user      = se_strdup(src_uri->user);
    uri->pass      = se_strdup(src_uri->pass);
    uri->host      = se_strdup(src_uri->host);
    uri->port        = src_uri->port;
    uri->strong_port = src_uri->strong_port;
    uri->path      = se_strdup(src_uri->path);
    uri->obj       = se_strdup(src_uri->obj);
    uri->ext       = se_strdup(src_uri->ext);
    uri->query     = se_strdup(src_uri->query);
    uri->frag      = se_strdup(src_uri->frag);
    uri->user_info = se_strdup(src_uri->user_info);
    uri->host_port = se_strdup(src_uri->host_port);
    uri->host_info = se_strdup(src_uri->host_info);
    uri->path_obj  = se_strdup(src_uri->path_obj);
    uri->path_obj_query = se_strdup(src_uri->path_obj_query);
    uri->ht_query = create_hash();
    hti = create_htiter( src_uri->ht_query );
    while( hti_hasnext(hti) ) {
        hti_next( hti, &q_key, &q_val );
        ht_put( uri->ht_query, se_strdup(q_key), 
                    se_strdup((const char*)q_val) );
    }
    destroy_htiter(hti);

    return uri;
}

/**
 * Destroy an the specified URI
 **/
void se_destroy_uri( se_uri _uri )
{
    _uri_t* uri = (_uri_t*) _uri;
    if (! _uri) {
        return;
    }

    if (uri->scheme)    free1char(uri->scheme);
    if (uri->scheme_specific) free1char(uri->scheme_specific);
    if (uri->user)      free1char(uri->user);
    if (uri->pass)      free1char(uri->pass);
    if (uri->host)      free1char(uri->host);
    if (uri->path)      free1char(uri->path);
    if (uri->obj)       free1char(uri->obj);
    if (uri->ext)       free1char(uri->ext);
    if (uri->query)     free1char(uri->query);
    if (uri->frag)      free1char(uri->frag);
    if (uri->user_info) free1char(uri->user_info);
    if (uri->host_port) free1char(uri->host_port);
    if (uri->host_info) free1char(uri->host_info);
    if (uri->path_obj)  free1char(uri->path_obj);
    if (uri->path_obj_query) free1char(uri->path_obj_query);
    if (uri->ht_query)  destroy_hash_and_entries(uri->ht_query, 1);

    free( uri );
}


/**
 * Parse the specified URI.
 * Returns a se_uri_t structure on success or NULL on failure.
 **/
se_uri se_parse_uri( const char* str, const char* default_scheme )
{
    _uri_t* uri;
    char *port, *param_key, *param_val, *old_val;
    const char *s, *s1, *s2, *s3;
    array params_array;
    int err = 1, i;
    ASSERT(str);

    uri = (_uri_t*)malloc( sizeof(_uri_t) );
    memset( uri, 0, sizeof(_uri_t) );
    s = str;

    do {
        /* parse scheme & scheme_specific */
        s1 = strstr( s, ":" );
        if (s1 == NULL) {
            /* relative uri */
            if (default_scheme != NULL) {
                uri->scheme = se_strdup( default_scheme );
            } else {
                uri->scheme = se_strdup("");
            }
            uri->is_abs = 0;
        } else {
            /* absolute uri */
            uri->scheme = se_strndup( s, s1 - str );
            if( default_scheme != NULL && 
                strcmp(uri->scheme, default_scheme) ) {
                LOG_CONSOLE(SE_WARNING, "Invalid URI '%s'. Expected '%s' scheme, got '%s'",
                      str, default_scheme, uri->scheme);
                break;
            }
            s = s1 + strlen(":");
            uri->is_abs = 1;
        }
        uri->scheme_specific = se_strdup(s);
        s1 = strchr( uri->scheme_specific, '#' );
        if (s1) {
            *((char*)(s1)) = '\0';
        }

        /* parse is_url */
        if (uri->is_abs) {
            if( str_startswith( s, "//" )) {
                uri->is_url = 1;
                s += strlen("//");
            } else {
                uri->is_url = 0;
            }
        } else {
            uri->is_url = 1;
        }

        /* parse host_info, path_obj, query & fragment */
        if (uri->is_abs) {
            s1 = strchr(s, '/');
            s2 = strchr(s, '?');
            s3 = strchr(s, '#');
            if (s1 > s2 && s2) s1 = s2;
            if (s2 > s3 && s3) s2 = s3;
            if (s1 > s3 && s3) s1 = s3;
            if(! s1) s1 = s2 ? s2 : s3 ? s3 : s + strlen(s);
            if(! s2) s2 = s + strlen(s);
            if(! s3) s3 = s + strlen(s);
            uri->host_info = se_strndup(s, s1 - s);
            if (*s1 == '/')
                uri->path_obj = se_strndup(s1 + 1, s2 - s1 - 1);
            else
                uri->path_obj = (uri->is_url) ? 
                    se_strdup("/") :
                    se_strdup("");
            if (*s2 == '?')
                uri->query = se_strndup(s2 + 1, s3 - s2 - 1);
            else
                uri->query = se_strdup("");
            if (*s3 == '#')
                uri->frag = se_strdup(s3 + 1);
            else
                uri->frag = se_strdup("");
        } else {
            s1 = s;
            s2 = strchr(s, '?');
            s3 = strchr(s, '#');
            if (s1 > s2 && s2) s1 = s2;
            if (s2 > s3 && s3) s2 = s3;
            if (s1 > s3 && s3) s1 = s3;
            if(! s1) s1 = s2 ? s2 : s3 ? s3 : s + strlen(s);
            if(! s2) s2 = s + strlen(s);
            if(! s3) s3 = s + strlen(s);
            uri->host_info = se_strndup(s, s1 - s);
            uri->path_obj = se_strndup(s1, s2 - s1);
            if (*s2 == '?')
                uri->query = se_strndup(s2 + 1, s3 - s2 - 1);
            else
                uri->query = se_strdup("");
            if (*s3 == '#')
                uri->frag = se_strdup(s3 + 1);
            else
                uri->frag = se_strdup("");
        }

        /* split host_info in user_info & host_port */
        se_uri_split_in_two_right( uri->host_info, '@',
                                   &uri->user_info, &uri->host_port );
        /* split user_info in user & pass */
        se_uri_split_in_two_left( uri->user_info, ':', &uri->user, &uri->pass );
        /* split host_port in host & port */
        se_uri_split_in_two_left( uri->host_port, ':', &uri->host, &port );

        /* check for empty host */
        if( strlen(uri->host) == 0 ) {
            if(! strcmp(uri->scheme, "file")) {
                /* default empty host to 'localhost' for 'file://' protocol
                   "file:///dir/f1.txt" > "file://localhost/dir/f1.txt" */
                free1char(uri->host);
                uri->host = se_strdup("localhost");

                /* also update host_port & host_info */
                free1char(uri->host_port);
                free1char(uri->host_info);
                uri->host_port = se_uri_merge_two( uri->host, port, ':' );
                uri->host_info = se_uri_merge_two( uri->user_info,
                                                   uri->host_port, '@' );
            } else if (uri->is_abs) {
                LOG_CONSOLE(SE_WARNING, "Invalid URI: Empty host: '%s'", str);
                break;
            }
        }

        /* parse the port */
        uri->port = -1;
        uri->strong_port = -1;
        if( strlen(port) ) {
            if( parse_int( port, &uri->port )) {
                LOG_CONSOLE(SE_WARNING, "Invalid URI: Cannot parse port: '%s'", str);
                free1char(port);
                break;
            }
            if( uri->port <= 0 ) {
                LOG_CONSOLE(SE_WARNING, "Invalid URI: Invalid port no: '%s'", str);
                free1char(port);
                break;
            }
        }
        free1char(port);
        if (uri->port == -1)
            uri->strong_port = se_uri_default_port(uri->scheme);

        /* split path_obj in path & 
           obj ('/dir/index.php' > '/dir/' & 'index.php') */
        if (uri->is_url) {
            s1 = strrchr( uri->path_obj, '/' );
            if(s1 == NULL) {
                uri->path = se_strdup("");
                uri->obj = se_strdup(uri->path_obj);
            } else {
                uri->path = se_strndup( uri->path_obj,
                                         s1 - uri->path_obj + 1 );
                uri->obj = se_strdup(s1 + 1);
            }
        } else {
            uri->path = se_strdup("");
            uri->obj = se_strdup("");
        }

        /* set path_obj_query as path_obj ['?' + query] */
        uri->path_obj_query = se_uri_merge_two( uri->path_obj,
                                                uri->query, '?' );

        /* extract ext from obj ('index.cgi' > 'cgi') */
        s1 = strrchr(uri->obj, '.');
        uri->ext = se_strdup(s1 ? s1 + 1 : "");

        /* parse query params and add them in ht_query */
        uri->ht_query = create_hash();
        params_array = split2( uri->query, '&',
                                   SPLIT_NO_EMPTYFIELDS | SPLIT_TRIMFIELDS );
        for (i = 0; i < array_size(params_array); ++i) {
            s = (char*) array_get_at( params_array, i);

            se_uri_split_in_two_left( s, '=', &param_key, &param_val );
            old_val = ht_put( uri->ht_query, param_key, param_val );
            if (old_val) {
                free1char(old_val);
                free1char(param_key);
            }
        }
        destroy_array( params_array, 1 );

        err = 0; /* success */
    }
    while (0);

    if (err) {
        se_destroy_uri( uri );
        return NULL;
    }

    return uri;
}

/**
 * Format the string corresponding to the specified URI
 * Use the normalize flag to get an normalized URI.
 **/
char* se_format_uri( const se_uri _uri, int normalize )
{
    _uri_t* uri = (_uri_t*) _uri;
    se_uri norm_uri = 0;
    strbuf sb;
    char* str;

    if(! uri) {
        return se_strdup("(null)");
    }

    /* create a normalized copy if necessary */
    if (normalize) {
        norm_uri = se_create_uri_copy( _uri );
        se_normalize_uri( norm_uri );
        uri = (_uri_t*) norm_uri;
    } else {
        uri = (_uri_t*) _uri;
    }

    sb = create_strbuf();
    if (uri->is_abs) {
        sb_append( sb, uri->scheme );
        sb_append( sb, ":" );
        if (uri->is_url) {
            sb_append( sb, "//" );
        }
        if( strcmp(uri->scheme, "file") || 
            strcmp(uri->host_info, "localhost") ) {
            sb_append( sb, uri->host_info );
        }
        if (uri->is_url) {
            sb_append( sb, "/" );
        }
    }
    sb_append( sb, uri->path_obj_query );
    if (uri->frag[0] != '\0' ) {
        sb_append( sb, "#" );
        sb_append( sb, uri->frag );
    }
    str = sb_release( sb );
    destroy_strbuf( sb );

    /* cleanup */
    if (normalize) {
        se_destroy_uri( uri );
    }

    return str;
}


/* absolute uri ? */
int se_uri_is_absolute( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->is_abs;
}

/* relative uri ? */
int se_uri_is_relative( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return ! uri->is_abs;
}

/* is URL ? */
int se_uri_is_url( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->is_url;
}

/* is URN ? */
int se_uri_is_urn( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return ! uri->is_url;
}

/* protocol ("http", "file", etc) */
const char* se_uri_get_scheme( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->scheme;
}

/* scheme-specific part */
const char* se_uri_get_scheme_specific( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->scheme_specific;
}

/* username */
const char* se_uri_get_user( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->user;
}

/* password */
const char* se_uri_get_pass( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->pass;
}

/* hostname or IP */
const char* se_uri_get_host( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->host;
}

/* port no if present, or -1 */
int se_uri_get_port( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->port;
}

/* port no if present, default port for scheme, or -1 */
int se_uri_get_strong_port( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->strong_port;
}

/* path ending in slash, default "/" ("/faq/dir/") */
const char* se_uri_get_path( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->path;
}

/* object ("index.php") */
const char* se_uri_get_obj( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->obj;
}

/* object's extension ("php") */
const char* se_uri_get_ext( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->ext;
}

/* everything between '?' and '#' */
const char* se_uri_get_query( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->query;
}

/* everything after '#' */
const char* se_uri_get_frag( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->frag;
}

/* combined user[:pass] */
const char* se_uri_get_user_info( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->user_info;
}

/* combined host[:port] */
const char* se_uri_get_host_port( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->host_port;
}

/* combined [user[:pass]@]host[:port] */
const char* se_uri_get_host_info( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->host_info;
}

/* combined path[obj] */
const char* se_uri_get_path_obj( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->path_obj;
}

/* combined path[obj]?query */
const char* se_uri_get_path_obj_query( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->path_obj_query;
}

/* query params (key/value hashtable) */
se_hash se_uri_get_query_params( const se_uri _uri )
{
    const _uri_t* uri = (const _uri_t*) _uri;
    ASSERT(uri);

    return uri->ht_query;
}


/**
 * Returns the default port for the specified scheme, or -1 if the scheme
 * is unknown. The scheme should be lower case
 **/
int se_uri_default_port( const char* scheme )
{
    ASSERT(scheme);

    if(! strcmp(scheme, "ftp"))      return   21;
    if(! strcmp(scheme, "ssh"))      return   22;
    if(! strcmp(scheme, "telnet"))   return   23;
    if(! strcmp(scheme, "smtp"))     return   25;
    if(! strcmp(scheme, "gopher"))   return   70;
    if(! strcmp(scheme, "http"))     return   80;
    if(! strcmp(scheme, "pop"))      return  110;
    if(! strcmp(scheme, "sftp"))     return  115;
    if(! strcmp(scheme, "nntp"))     return  119;
    if(! strcmp(scheme, "imap"))     return  143;
    if(! strcmp(scheme, "prospero")) return  191;
    if(! strcmp(scheme, "wais"))     return  210;
    if(! strcmp(scheme, "ldap"))     return  389;
    if(! strcmp(scheme, "https"))    return  443;
    if(! strcmp(scheme, "smtps"))    return  465;
    if(! strcmp(scheme, "rtsp"))     return  554;
    if(! strcmp(scheme, "snews"))    return  563;
    if(! strcmp(scheme, "acap"))     return  674;
    if(! strcmp(scheme, "ftps"))     return  990;
    if(! strcmp(scheme, "imaps"))    return  993;
    if(! strcmp(scheme, "pop3s"))    return  995;
    if(! strcmp(scheme, "nfs"))      return 2049;
    if(! strcmp(scheme, "tip"))      return 3372;
    if(! strcmp(scheme, "sip"))      return 5060;
    if(! strcmp(scheme, "sips"))     return 5061;

    return -1; /* unknown scheme */
}


/**
 * Compares the specified normalized URIs.
 * The path & object (/dir/x.xml) are compared case-sensitive.
 * Return true (non-zero) if the URIs are equals, or zero otherwise
 **/
int se_uri_equals( const se_uri _norm_uri1, const se_uri _norm_uri2 )
{
    const _uri_t* norm_uri1 = (const _uri_t*) _norm_uri1;
    const _uri_t* norm_uri2 = (const _uri_t*) _norm_uri2;
    ASSERT(norm_uri1);
    ASSERT(norm_uri2);
    ASSERT(norm_uri1->is_abs);
    ASSERT(norm_uri2->is_abs);

    if( norm_uri1->is_url != norm_uri2->is_url )
        return 0;
    if( strcmp( norm_uri1->scheme, norm_uri2->scheme ))
        return 0;
    if( norm_uri1->is_url ) {
        if( strcmp( norm_uri1->host_info, norm_uri2->host_info ))
            return 0;
        if( strcmp( norm_uri1->path_obj_query, norm_uri2->path_obj_query ))
            return 0;
    }
    else {
        if( strcmp( norm_uri1->scheme_specific, norm_uri2->scheme_specific ))
            return 0;
    }
    if( strcmp( norm_uri1->frag, norm_uri2->frag ))
        return 0;

    return 1;
}


/**
 * Normalizes the specified URI
 * eg: "http://site.com/dir/dir1/../index.php" > "http://site.com/dir/index.php"
 **/
void se_normalize_uri( se_uri _uri )
{
    _uri_t* uri = (_uri_t*) _uri;
    char* path;

    ASSERT (uri);
    ASSERT (uri->is_url);
    VERIFYM(uri->is_abs, "Normalizing relative URIs is not implemented");

    /* normalize the path */
    path = se_normalize_path( uri->path );
    if(! strcmp( path, uri->path ))
    {
        free1char(path);
        return;
    }
    if (uri->path) free1char(uri->path);
    uri->path = path;

    /* update path_obj field */
    if (uri->path_obj) free1char(uri->path_obj);
    uri->path_obj = se_uri_merge_two( uri->path, uri->obj, '/' );

    /* update path_obj_query field */
    if (uri->path_obj_query) free1char(uri->path_obj_query);
    uri->path_obj_query = se_uri_merge_two( uri->path_obj, uri->query, '?' );

    /* TODO: update scheme specific */
}

/**
 * Creates a base URI from the specified absolute URI
 * eg: "ftp://site.com/dir/index.php" > "ftp://site.com/dir/"
 **/
se_uri se_base_uri( const se_uri _abs_uri )
{
    const _uri_t* abs_uri = (_uri_t*) _abs_uri;
    _uri_t* base_uri;

    VERIFY(abs_uri);
    VERIFY(abs_uri->is_abs);
    VERIFY(abs_uri->is_url);

    base_uri = (_uri_t*) se_create_uri_copy( _abs_uri );

    /* clear obj and ext */
    if (base_uri->obj) free1char(base_uri->obj);
    base_uri->obj = se_strdup("");
    if (base_uri->ext) free1char(base_uri->ext);
    base_uri->ext = se_strdup("");

    /* update path_obj: clear obj */
    if (base_uri->path_obj) free1char(base_uri->path_obj);
    base_uri->path_obj = se_strdup(base_uri->path);

    /* update path_obj_query: clear obj */
    if (base_uri->path_obj_query) free1char(base_uri->path_obj_query);
    base_uri->path_obj_query = se_uri_merge_two( base_uri->path_obj,
                                                 base_uri->query, '?' );

    /* TODO: update scheme specific */

    return base_uri;
}


/**
 * Creates a absolute URI from the specified relative & base URIs
 * eg: "../dir2/index.cgi" + "http://site.com/dir/dir1/faq.cgi" > "http://site.com/dir/dir2/index.cgi"
 **/
se_uri se_uri_rel2abs( const se_uri _rel_uri, const se_uri _base_uri )
{
    const _uri_t* rel_uri  = (const _uri_t*) _rel_uri;
    const _uri_t* base_uri = (const _uri_t*) _base_uri;
    const char* s;
    _uri_t* abs_uri;
    char *abs_path, *base_path;
    ASSERT(rel_uri);
    ASSERT(! rel_uri->is_abs);
    ASSERT(base_uri);
    ASSERT(base_uri->is_abs);

    /* get the absolute path */
    if( str_startswith(base_uri->path, "/") ) {
        base_path = se_strdup(base_uri->path);
    } else {
        base_path = se_uri_merge_two("/", base_uri->path, '/');
    }
    abs_path = se_make_path_absolute( rel_uri->path, base_path );
    if (str_startswith(abs_path, "/") ) {
        memmove( abs_path, abs_path + 1, strlen(abs_path) );
    }
    free1char(base_path);

    /* create a copy of base uri */
    abs_uri = (_uri_t*) se_create_uri_copy( _base_uri );

    /* change it's path & object fields */
    free1char( abs_uri->path );
    abs_uri->path = abs_path;
    free1char( abs_uri->obj );
    abs_uri->obj = se_strdup( rel_uri->obj );

    /* update it's scheme_specific, ext, path_obj & path_obj_query fields */
    if (abs_uri->ext) free1char(abs_uri->ext);
    s = strrchr(abs_uri->obj, '.');
    abs_uri->ext = se_strdup(s ? s + 1 : "");
    if (abs_uri->path_obj) free1char(abs_uri->path_obj);
    abs_uri->path_obj = se_uri_merge_two( abs_uri->path, abs_uri->obj, '/' );
    if (abs_uri->path_obj_query) free1char(abs_uri->path_obj_query);
    abs_uri->path_obj_query = se_uri_merge_two( abs_uri->path_obj,
                                                abs_uri->query, '?' );
    /* TODO: update scheme specific */

    return abs_uri;
}

/**
 * Creates a relative URI from the specified absolute & base URIs
 * eg: "http://site.com/dir/dir2/index.cgi" + "http://site.com/dir/dir1/" > "../dir2/index.cgi"
 **/
se_uri se_uri_abs2rel( const se_uri _abs_uri, const se_uri _base_uri )
{
    _uri_t* abs_uri = (_uri_t*) _abs_uri;
    _uri_t* base_uri = (_uri_t*) _base_uri;
    VERIFY(abs_uri);
    VERIFY(abs_uri->is_abs);
    VERIFY(base_uri);
    VERIFY(base_uri->is_abs);

    /* TODO */
VERIFYM(0, "not implemented");
return NULL;
}


