#include "../include/se_fs_sep.h"
#include "../include/se_module_par_desc_sep.h"
#include "../include/se_coordinate_system_par.h"
// 增加行缓冲区大小以确保足够容纳所有数据
#define SEP_H_LINESIZE 16384

const char* const SEP_LABELS_X_AXES[] = {
    SEP_DEFAULT_LABEL_X,         //“crossline"
    "x",
    "crossline",
    "x-lines",
    "x-line",
    "xlines",
    "xline",
    "Crosslines (X)",
    NULL
};

const char* const SEP_LABELS_Y_AXES[] = {
    SEP_DEFAULT_LABEL_Y,
    "y",
    "inline",
    "i-lines",
    "i-line",
    "ilines",
    "iline",
    "Inlines (Y)",
    NULL
};


static int sep_is_some_axis(const char* label, const char* const names[])
{
    if(label)
        for(; *names; ++names) {
            if(strcasecmp(label, *names) == 0)
                return 1;
        }
    return 0;
}

int sep_is_x_axis(const char* label)
{
    return sep_is_some_axis(label, SEP_LABELS_X_AXES);
}

int sep_is_y_axis(const char* label)
{
    return sep_is_some_axis(label, SEP_LABELS_Y_AXES);
}


static int* sep_clone_int( int value )
{
    int* ptr = (int*) malloc( sizeof(*ptr) );
    *ptr = value;
    return ptr;
}

static double* sep_clone_double( double value )
{
    double* ptr = alloc1double( 1 );
    *ptr = value;
    return ptr;
}

static char* _sep_print_a_double(double d, char* s, size_t maxsize)
{
    if(is_integer(d)) {
        snprintf(s, maxsize, "%lld", (long long int)d );
    } else {
        char* p;

        snprintf(s, maxsize, "%.12f", d);
        p = strrchr(s, '.');
        if(p) {
            char* z = NULL;
            ++p;
            while(*p) {
                if(*p == '0') {
                    if(z == NULL) z = p;
                } else {
                    z = NULL;
                }
                ++p;
            }
            if(z) {
                *z = 0;
            }
        }
    }
    return s;
}


static void sep_parse_ds_path( const char* path_uri,
                                       char** ds_path,
                                       int* data_only )
{
    se_uri uri;
    char* tmp;

    ASSERT(path_uri);
    ASSERT(ds_path);
    ASSERT(data_only);

    /* check for 'stdin' 'stdout' and '-' first */
    if( strcmp(path_uri, "stdin") == 0 ||
        strcmp(path_uri, "stdout") == 0 ||
        strcmp(path_uri, "-") == 0) {
        *ds_path = strdup( path_uri );
        *data_only = 0;
        return;
    }

    /* parse the uri */
    uri = se_parse_uri( path_uri, "file" );
    if(!uri) {
        ERROR(("Cannot parse file URI: '%s'", path_uri));
    }

    /* extract normalized absolute sqlite file path */
    *ds_path = se_abs_path_from_peer( se_uri_get_path_obj(uri), "." );
    tmp = se_normalize_path( *ds_path );
    free1char(*ds_path);
    *ds_path = tmp;

    /* set data_only flag */
    *data_only = ht_contains_key( se_uri_get_query_params(uri),
                                      "dataonly" );

    /* cleanup */
    se_destroy_uri( uri );
}

static sep_headers_t* sep_headers_create( void )
{
    sep_headers_t* headers = (sep_headers_t*)malloc( sizeof(*headers) );
    int i;

    headers->filename = NULL;
    headers->fd = -1;
    headers->ndim = 0;
    for(i = 0; i < SEP_NDIM_MAX; ++i) {
        headers->n[i] = 0;
        headers->o[i] = 0.0;
        headers->d[i] = 0.0;
    }
    headers->in = strdup("");
    headers->esize = 0;
    headers->le = 0;
    headers->ht_int = create_hash();
    headers->ht_float = create_hash();
    headers->ht_str = create_hash();
    headers->raw_size = 0;
    headers->history = NULL;
    headers->has_inline_data = 0;
    headers->is_modified = 0;
    headers->is_commited = 0;

    return headers;
}

static void sep_headers_destroy( sep_headers_t* headers )
{
    if (headers->fd >= 0) {
        int rc;
        do {
            rc = close( headers->fd );
        } while (rc < 0 && errno == EINTR);
        if (rc < 0) {
            int err = errno;
        
           LOG_CONSOLE(SE_WARNING, "Cannot close SEP file '%s': current dir is %s, "
                   "error code is %d - %s",
                   headers->filename, get_cwd(), err, strerror(err) );
        }
    }

    if(headers->filename) free1char( headers->filename );
    if(headers->in) free1char( headers->in );
    if(headers->expanded_in) free1char( headers->expanded_in );
    if(headers->history) free1char( headers->history );
    destroy_hash_and_entries( headers->ht_int, 1 );
    destroy_hash_and_entries( headers->ht_float, 1 );
    destroy_hash_and_entries( headers->ht_str, 1 );

    free( headers );
}


/** Reads the next byte from 'fd', returns the byte or EOF */
static int sep_headers_readbyte( int fd, const char* filename )
{
    byte b;
    ssize_t rc;

    do {
        rc = read( fd, &b, 1 );
    } while (rc < 0 && errno == EINTR);

    if (rc < 0) {
        int err = errno;
        ERROR(( "Cannot read from SEP file '%s': current dir is %s, "
                "error code is %d - %s",
                filename, get_cwd(), err, strerror(err) ));
    }

    return (rc == 0) ? EOF : (int) b;
}

/**
 * Reads a line from 'fd' into the buffer pointed by 'buffer'.
 *
 * Reading stops when:
 * - An '\n' is found,
 * - An EOT (ASCII 0x04) is found (the '*eot_found' flag is set),
 * - End of file is reached.
 *
 * The '\n' and EOT chars are copied into the buffer.
 *
 * Returns 0 on success or EOF (-1) on end of file.
 */
static int sep_headers_readline( int fd, const char* filename,
                                         char* buffer, int buffer_size,
                                         int* nread, int* eot_found )
{
    char* s = buffer;
    int c = 0;

    ASSERT(fd >= 0);
    ASSERT(filename);
    ASSERT(buffer);
    ASSERT(buffer_size > 0);
    ASSERT(eot_found);

    if (nread) *nread = 0;

    while( --buffer_size > 0 ) {
        c = sep_headers_readbyte( fd, filename );
        if( c == EOF )
            break;
        if( nread )
            (*nread)++;

        *s++ = (char)c;
        if( c == '\n' || c == SEP_EOT )
            break;
    }

    *s = '\0';
    *eot_found = (c == SEP_EOT);

    return (c == EOF) ? EOF : 0;
}

static int sep_headers_syncdata( sep_headers_t* headers,
                                         char* line )
{
    int pid, len = strlen(line);
    ASSERT(headers);
    ASSERT(line);

    /* check for Simple Method (EOL EOL EOT separator) */
    if( len >= 3 && line[len - 2] == SEP_EOL &&
        line[len - 3] == SEP_EOL ) {
        /* EOL EOL EOT found, SEP .H file has inline data */
        headers->has_inline_data = 1;
        /* rewind the two EOLs */
        line[len - 3] = '\0';
        return 0;
    }

    /* check for Pipe Method (INET socket data) */
    if( str_startswith(line, "PIPE ") ) {
        char pipe[16], host[256], ack[16];
        int port;
        se_socket socket;

        /* parse "PIPE <host> <port> <pid>" */
        if( 4 != sscanf(line, "%s %s %d %d", pipe, host, &port, &pid) ) {
            ERROR(("Invalid SEP headers/data synchronization info for "
                   "file '%s': expected 'PIPE <host> <port> <pid>', got '%s'",
                   headers->filename, line));
        }

        /* open socket to <host>:<port> */
        socket = se_socket_connect(host, port);
        if (! socket) {
            ERROR(("Cannot synchronize SEP headers/data for file '%s': "
                   "Cannot connect socket to %s:%d",
                   headers->filename, host, port));
        }

        /* write "GOTIT\0" */
        if( se_socket_write(socket, "GOTIT", 6) == CODE_ERROR ) {
            ERROR(("Cannot synchronize SEP headers/data for file '%s': "
                   "Cannot write to socket %s:%d",
                   headers->filename, host, port));
        }

        /* read "ACK\0" */
        if( se_socket_read(socket, ack, 4) == CODE_ERROR ) {
            ERROR(("Cannot synchronize SEP headers/data for file '%s': "
                   "Cannot read SEP data through socket %s:%d",
                   headers->filename, host, port));
        }
        if( strcmp(ack, "ACK") != 0 ) {
            ERROR(("Cannot synchronize SEP headers/data for file '%s': "
                   "expected 'ACK', got '%s'",
                   headers->filename, ack));
        }
        se_socket_close(socket);

        headers->has_inline_data = 1;
        return 1;
    }

    /* send SIGALRM signal */
    if( sscanf(line, "%d", &pid) == 1 ) {
        if( kill( pid, SIGALRM ) < 0) {
            int err = errno;
            ERROR(("Cannot synchronize SEP headers/data for file '%s': "
                   "Failed to send SIGALRM signal to process %d: errno=%d - %s",
                   headers->filename, pid, err, strerror(err)));
        }
        return 1;
    }

    ERROR(("Unknown SEP headers/data synchronization method for file '%s': "
           "sync. line is '%s'", headers->filename, line));
    return 1;
}

static void sep_set_sepheader(sep_headers_t* headers, const char* name, const char* value)
{
    se_hash ht;
    int ival;
    double dval;
    void* val;
    void* old_val;

    if (! parse_int(value, &ival)) {
        ht = headers->ht_int;
        val = sep_clone_int(ival);
    } else if (! parse_double(value, &dval)) {
        ht = headers->ht_float;
        val = sep_clone_double(dval);
    } else {
        ht = headers->ht_str;
        val = strdup(value);
    }

    if(! ht_contains_key( ht, name )) {
        name = strdup(name);
    }

    old_val = (void*) ht_put( ht, name, val );
    if (old_val) free(old_val);

    headers->is_modified = 1;
}

void sep_set_header(sep_t* sep, const char* name, const char* value)
{
    sep_set_sepheader(sep->headers, name, value);
}

static int sep_axis_from_char(const char c)
{
    switch (c) {
    case '1': return 0;
    case '2': return 1;
    case '3': return 2;
    case '4': return 3;
    case '5': return 4;
    case '6': return 5;
    case '7': return 6;
    case '8': return 7;
    case '9': return 8;
    default: return -1;
    }
}
static int sep_set_header_intercept(sep_t* sep, const char* name, double value)  //si_map_include_all
{
    int i, ret = 0;
    if(name && name[0] && name[1] && name[2] == 0) {
        i = sep_axis_from_char(name[1]);
        if(i >= 0) {
            switch (name[0]) {
            case 'n': sep->headers->n[i] = (int)value; ret = 1; break;
            case 'd': sep->headers->d[i] = value; ret = 1; break;
            case 'o': sep->headers->o[i] = value; ret = 1; break;
            }
        }
    }
    if(ret) {
        if(sep->headers->ndim < i+1) 
            sep->headers->ndim = i+1;
        sep->headers->is_modified = 1;
    }

    return ret;
}

void sep_set_header_int(sep_t* sep, const char* name, int value)
{
    void* val;
    void* old_val;

    if(sep_set_header_intercept(sep, name, (double)value)) return;   //判断是否中断

    if(! ht_contains_key( sep->headers->ht_int, name )) {         //ht_int
        name = strdup(name);
    }

    val = sep_clone_int(value);
    old_val = (void*) ht_put( sep->headers->ht_int, name, val );
    if (old_val) free(old_val);
    sep->headers->is_modified = 1;
}

void sep_set_header_float(sep_t* sep, const char* name, double value)
{
    void* val;
    void* old_val;
 
    if(sep_set_header_intercept(sep, name, (double)value)) return;
    if(! ht_contains_key( sep->headers->ht_float, name )) {
        name = strdup(name);
    }

    val = sep_clone_double(value);         //克隆alue
    old_val = (void*) ht_put( sep->headers->ht_float, name, val );
    if (old_val) free(old_val);
    sep->headers->is_modified = 1;
}


static void sep_headers_parse( sep_headers_t* headers,
                                       int enable_warnings )
{
    char line[SEP_H_LINESIZE], **fields = NULL, **h_nameval = NULL;
    const char *h_name, *h_val;
    int field_nr = 0, h_nameval_nr = 0;
    int hidx, ival, i, line_idx, done, nontext_counter = 0;
    double dval;
    strbuf sb_history;

    ASSERT(headers);

    /* set the default values for the optionsl headers
       ('n1' and 'in' headers are required) */
    for(i = 0; i < SEP_NDIM_MAX; ++i) {
        headers->n[i] = (i == 0) ? 0 : 1;
        headers->o[i] = 0.0;
        headers->d[i] = 1.0;
    }
    headers->esize = 4;
    headers->le = 0;

    sb_history = create_strbuf_args( 256 );

    /* parse the SEP headers from the stream */
    done = 0;
    line_idx = 0;
    nontext_counter = 0;
    while (! done) {
        int nread = 0, eot_found = 0;

        /* read next line */
        if( EOF == sep_headers_readline( headers->fd, headers->filename,
                                            line, SEP_H_LINESIZE,
                                            &nread, &eot_found)) {
            done = 1;   /* stop parsing after this line */
        }
        line_idx++;

        /* check for non-printable chars (other than EOL & EOT) */
        for (i = 0; i < nread; ++i) {
            char c = line[i];
            if(! isprint(c) && ! isspace(c) &&
               c != SEP_EOL && c != SEP_EOT) {
                ++nontext_counter;
                if (nontext_counter >= SEP_MAX_NONTEXT_CHARS) {
                    WARN(("SEP header parser: too many non-text characters "
                           "found in file %s. Probably this is a binary file.",
                           headers->filename));
                }
            }
        }

        /* check for lines longer than our buffer */
        if (nread == SEP_H_LINESIZE - 1 && line[nread - 1] != '\n') {
            ERROR(("SEP header parser: line too long found in file %s:%d",
                   headers->filename, line_idx));
        }

        /* update headers raw size */
        headers->raw_size += nread;

        /* check for data synchronization method */
        if (eot_found) {
            done = 1;
            if (nontext_counter > 0) {
                WARN(("SEP header parser: %d non-text character(s) "
                      "found in file %s ",
                      nontext_counter, headers->filename));
            }
            if( sep_headers_syncdata( headers, line )) {
                /* don't parse the current line anymore */
                break;
            }
        }

        /* add the line to history buffer */
        sb_append( sb_history, line );

        /* strip #-comments and trim white spaces */
        strip_comment(line, '#');
        str_trim(line);

        /* split the line in space separated fields */
        free_str_array(fields, field_nr);
        split( line, ' ', &fields, &field_nr,
                   SPLIT_QUOTES | SPLIT_TRIMFIELDS | SPLIT_NO_EMPTYFIELDS | SPLIT_USE_ISSPACE );

        for (i = 0; i < field_nr; ++i) {
            /* split the line in '='-separated header name and value */
            free_str_array(h_nameval, h_nameval_nr);
            split( fields[i], '=', &h_nameval, &h_nameval_nr,
                       SPLIT_QUOTES );
            if( h_nameval_nr != 2 )
                continue;
            str_trim(h_nameval[0]);

            /* get header name and value */
            h_name = h_nameval[0];
            h_val = h_nameval[1];
            if( strlen(h_name) == 0 )
                continue;

            /* let's see what header did we found ... */

            /* n1 to n9 headers */
            if( str_startswith( h_name, "n" )) {
                if(! parse_int(h_name + 1, &hidx) &&
                   hidx >= 1 && hidx <= SEP_NDIM_MAX ) {
                    if(! parse_double(h_val, &dval) && dval >= 1.0) {
                        headers->n[hidx - 1] = (int) dval;
                        headers->ndim = hidx > headers->ndim ?
                            hidx : headers->ndim;
                    } else if (enable_warnings) {
                        WARN(("Unexpected SEP global header value in "
                              "file '%s': '%s', header ignored",
                              headers->filename, fields[i]));
                    }
                    continue;
                }
            }

            /* o1 to o9 headers */
            if( str_startswith( h_name, "o" )) {
                if(! parse_int(h_name + 1, &hidx) &&
                   hidx >= 1 && hidx <= SEP_NDIM_MAX ) {
                    if(! parse_double(h_val, &dval)) {
                        headers->o[hidx - 1] = dval;
                        headers->ndim = hidx > headers->ndim ?
                            hidx : headers->ndim;
                    } else if (enable_warnings) {
                        WARN(("Unexpected SEP global header value in "
                              "file '%s': '%s', header ignored",
                              headers->filename, fields[i]));
                    }
                    continue;
                }
            }

            /* d1 to d9 headers */
            if( str_startswith( h_name, "d" )) {
                if(! parse_int(h_name + 1, &hidx) &&
                   hidx >= 1 && hidx <= SEP_NDIM_MAX ) {
                    if(! parse_double(h_val, &dval) &&
                       ! FEQUAL(dval, 0.0) ) {
                        headers->d[hidx - 1] = dval;
                        headers->ndim = hidx > headers->ndim ?
                            hidx : headers->ndim;
                    } else if (enable_warnings) {
                        WARN(("Unexpected SEP global header value in "
                              "file '%s': '%s', header ignored",
                              headers->filename, fields[i]));
                    }
                    continue;
                }
            }

            /* 'in' header */
            if(! strcmp(h_name, "in")) {
                if (strlen(h_val) > 0) {
                    free1char(headers->in);
                    headers->in = alloc1char(strlen(h_val) + 1);
                    memcpy( headers->in, h_val, (strlen(h_val) + 1) * sizeof(char) );
                } else if (enable_warnings) {
                    WARN(("Unexpected SEP global header value in "
                          "file '%s': '%s', header ignored",
                          headers->filename, fields[i]));
                }

                continue;
            }

            /* 'esize' header */
            if(! strcmp(h_name, "esize")) {
                if (! parse_int(h_val, &ival) &&
                    (ival == 1 || ival == 2 || ival == 4 || ival == 8 || ival == 12)) {
                    headers->esize = ival;
                }
                else if (enable_warnings) {
                    WARN(("Unexpected SEP global header value in "
                          "file '%s': '%s' "
                          "(expected 1, 2, 4, 8, 12), header ignored",
                          headers->filename, fields[i]));
                }
                continue;
            }

            /* 'data_format' header */
            if(! strcmp(h_name, "data_format")) {
                if(! strcmp(h_val, "native_float"))
                    headers->le = 1;
                else if(! strcmp(h_val, "xdr_float"))
                    headers->le = 0;
                else if (enable_warnings) {
                    WARN(("Unexpected SEP global header value in "
                          "file '%s': '%s' (expected 'xdr_float' or "
                          "'native_float'), header ignored",
                          headers->filename, fields[i]));
                }
                continue;
            }

            /* other header, add it to one of the hash tables */
            sep_set_sepheader(headers, h_name, h_val);
        }
    }

    /* cleanup */
    free_str_array( fields, field_nr );
    free_str_array( h_nameval, h_nameval_nr );

    if (nontext_counter > 0) {
        WARN(("SEP header parser: %d non-text character(s) found in file %s ",
              nontext_counter, headers->filename));
    }

    /* check for the required headers ('n1' & 'in') */
    // if (headers->n[0] == 0) {
    //     ERROR(("Missing required SEP global header 'n1' in file '%s'",
    //            headers->filename));
    // }
    if (! strlen(headers->in)) {
        ERROR(("Missing required SEP global header 'in' in file '%s'",
               headers->filename));
    }

    /* replace '-' with 'stdin' */
    if( strcmp(headers->in, "-") == 0 ) {
        free(headers->in);
        headers->in = strdup("stdin");
    }

    /* set history */
    headers->history = sb_release( sb_history );
    destroy_strbuf( sb_history );

    /* reset the modified flag */
    headers->is_modified = 0;
}


static char*
sep_expand_binary_path( const sep_t* sep,
                           const char* ds_path,
                           const sep_headers_t* headers)
{
    char* filename;
    if( sep->data_only ) {
        filename = strdup( ds_path );
    } else {
        if( strcmp(headers->in, "-") == 0 ||
            strcmp(headers->in, "stdin") == 0 ||
            strcmp(headers->in, "stdout") == 0) {
            filename = strdup( headers->in );
        } else {
            filename = se_abs_path_from_peer( headers->in,
                                               headers->filename );
        }
    }
    return filename;
}


static sep_headers_t*
sep_headers_open( const sep_t* sep, const char* ds_path )
{
    sep_headers_t* headers;
    int i;

    headers = sep_headers_create();

    if (sep->data_only) {
        /* data only */
        headers->ndim = 1;
        for (i = 0; i < SEP_NDIM_MAX; ++i) {
            headers->n[i] = 1;
            headers->o[i] = 0.0;
            headers->d[i] = 1.0;
        }
        headers->esize = 4;
        headers->le = 1;
    } else {
        /* replace '-' with 'stdin' for read or 'stdout' for write */
        if( strcmp(ds_path, "-") == 0 ) {
            ds_path = (sep->mode & SEP_READ) ? "stdin":"stdout";
        }
        headers->filename = strdup(ds_path);

        if( sep->mode & SEP_READ ) {
            /* open and read SEP header file (.H) */
            if( strcmp(headers->filename, "stdin") == 0 ) {
                headers->fd = dup(STDIN_FILENO);
            } else {
                headers->fd = open( headers->filename,
                                    (sep->mode & SEP_WRITE) ?
                                    O_RDWR | O_CREAT : O_RDONLY,
                                    S_IRUSR|S_IWUSR |
                                    S_IRGRP|S_IWGRP |
                                    S_IROTH|S_IWOTH );
            }
            if( headers->fd < 0 ) {
                int err = errno;
                ERROR(( "Cannot open SEP file '%s': current dir is %s, "
                        "error code is %d - %s",
                        headers->filename, get_cwd(), err, strerror(err) ));
            }

            /* parse SEP global headers */
            sep_headers_parse( headers, 1 /* warns on */);

        } else { /* completely new header file */

            /* open the SEP header file (.H) */
            if( strcmp(headers->filename, "stdout") == 0 ) {
                headers->fd = dup(STDOUT_FILENO);
            } else {
                headers->fd = open( headers->filename,
                                    O_RDWR | O_CREAT | O_TRUNC,
                                    S_IRUSR|S_IWUSR |
                                    S_IRGRP|S_IWGRP |
                                    S_IROTH|S_IWOTH );
            }
            if( headers->fd < 0 ) {
                int err = errno;
                ERROR(( "Cannot open SEP file '%s': current dir is %s, "
                        "error code is %d - %s",
                        headers->filename, get_cwd(), err, strerror(err) ));
            }

            headers->ndim = 1;
            for (i = 0; i < SEP_NDIM_MAX; ++i) {
                headers->n[i] = 1;
                headers->o[i] = 0.0;
                headers->d[i] = 1.0;
            }
            headers->esize = 4;
            headers->le = 1;
            headers->is_modified = 1;

            if( strcmp(headers->filename, "stdout") != 0 ) {
                /* set headers->in to "<filename>@" */
                if (headers->in) free1char(headers->in);
                headers->in = asprintf( "%s@", headers->filename );
            } else {
                if (headers->in) free1char(headers->in);
                headers->in = strdup("stdin");
                headers->has_inline_data = 1;
            }
        }
    }

    headers->expanded_in = sep_expand_binary_path( sep, ds_path, headers );

    return headers;
}

static void sep_headers_append( sep_headers_t* headers )
{
    se_fsio* io;
    char line[SEP_H_LINESIZE], *cwd, *user, *host, *tmp;
    char fltstring[SEP_H_LINESIZE / 2]; // 使用一半的SEP_H_LINESIZE确保格式化后的字符串不会超过line
    htiter hti;
    const char *h_name, *exe;
    void* h_vval;
    time_t t;
    int i;

    /* save global headers only if modified */
    if(! headers->is_modified) {
        return;  /* success */
    }

    /* append begin */

    io = se_fsio_attach( headers->filename, headers->fd, NULL );

    /* attached */

    /* write the history first */
    if (headers->history) {
        se_fsio_raw_write( io, headers->history, strlen(headers->history) );
    }

    /* write exe name (robust against NULL) */
    tmp = strdup("\n");
    se_fsio_raw_write( io, tmp, strlen(tmp) );
    free(tmp);
    exe = se_get_exe_name();
    if (exe == NULL || *exe == '\0') {
        exe = "unknown_exe";
    }
    se_fsio_raw_write( io, exe, strlen(exe) );

    /* write current dir */
    cwd = get_cwd();
    if (cwd == NULL) {
        cwd = strdup(".");
    }
    tmp = strdup("\t\"");
    se_fsio_raw_write( io, tmp, strlen(tmp) );
    free(tmp);
    se_fsio_raw_write( io, cwd, strlen(cwd) );
    tmp = strdup("\"");
    se_fsio_raw_write( io, tmp, strlen(tmp) );
    free(tmp);
    free1char(cwd);

    /* wrote exe and cwd */

    /* write user@host */
    host = get_hostname();
    if (host == NULL) host = strdup("unknown_host");
    user = get_user();
    if (user == NULL) user = strdup("unknown_user");
    snprintf( line, SEP_H_LINESIZE, "\t%s@%s\t", user, host );
    se_fsio_raw_write( io, line, strlen(line) );
    free(user);
    free(host);

    /* write date & time */
    t = time(NULL);
    ctime_r( &t, line );
    se_fsio_raw_write( io, line, strlen(line) );

    /* write in header */
    {
        char* dirh  = se_get_base_path( headers->filename );
        char* dirin = se_get_base_path( headers->in );
        char* in;
        if( se_is_same_file(dirh, dirin) ) {
            in = se_get_filename( headers->in );
        } else {
            in = strdup(headers->in);
        }
        snprintf( line, SEP_H_LINESIZE, "\tin=\"%s\"\n", in );
        se_fsio_raw_write( io, line, strlen(line) );
        free(in);
        free(dirin);
        free(dirh);
    }

    /* wrote in header */

    /* write data_format header */
    snprintf( line, SEP_H_LINESIZE, "\tdata_format=\"%s\"\n",
              (headers->le == 1) ? "native_float" : "xdr_float" );
    se_fsio_raw_write( io, line, strlen(line) );

    /* write n1 n2... headers */
    for(i = 0; i < headers->ndim; ++i) {
        snprintf( line, SEP_H_LINESIZE, "\tn%d=%d\n", i + 1, headers->n[i] );
        se_fsio_raw_write( io, line, strlen(line) );
    }

    /* wrote n headers */

    /* write d1 d2... headers */
    for(i = 0; i < headers->ndim; ++i) {
        // 限制fltstring的长度，确保不会超过缓冲区
        _sep_print_a_double(headers->d[i], fltstring, sizeof(fltstring) - 1);
        
        // 直接格式化到最终的行缓冲区
        snprintf(line, sizeof(line) - 1, "\td%d=%s\n", i + 1, fltstring);
        se_fsio_raw_write(io, line, strlen(line));
    }

    /* wrote d headers */

    /* write o1 o2... headers */
    for(i = 0; i < headers->ndim; ++i) {
        // 限制fltstring的长度，确保不会超过缓冲区
        _sep_print_a_double(headers->o[i], fltstring, sizeof(fltstring) - 1);
        
        // 直接格式化到最终的行缓冲区
        snprintf(line, sizeof(line) - 1, "\to%d=%s\n", i + 1, fltstring);
        se_fsio_raw_write(io, line, strlen(line));
    }

    /* wrote o headers */

    /* write esize header */
    snprintf( line, SEP_H_LINESIZE, "\tesize=%d\n", headers->esize );
    se_fsio_raw_write( io, line, strlen(line) );

    /* write int global headers */
    hti = create_htiter( headers->ht_int );
    while( hti_hasnext(hti) ) {
        hti_next( hti, &h_name, (void**) &h_vval );
        snprintf( line, SEP_H_LINESIZE, "\t%s=%d\n",
                  h_name, ((int*) h_vval)[0] );
        se_fsio_raw_write( io, line, strlen(line) );
    }
    destroy_htiter( hti );

    /* wrote int headers */

    /* write float global headers */
    hti = create_htiter( headers->ht_float );
    while( hti_hasnext(hti) ) {
        double fval;
        hti_next( hti, &h_name, (void**) &h_vval );
        fval = ((double*)h_vval)[0];
        
        // 限制fltstring的长度，确保不会超过缓冲区
        _sep_print_a_double(fval, fltstring, sizeof(fltstring) - 1);
        
        // 直接格式化到最终的行缓冲区
        snprintf(line, sizeof(line) - 1, "\t%s=%s\n", h_name, fltstring);
        se_fsio_raw_write(io, line, strlen(line));
    }
    destroy_htiter( hti );

    /* wrote float headers */

    /* write string global headers */
    hti = create_htiter( headers->ht_str );
    while( hti_hasnext(hti) ) {
        hti_next( hti, &h_name, (void**) &h_vval );
        snprintf( line, SEP_H_LINESIZE, "\t%s=\"%s\"\n",
                  h_name, (char*) h_vval );
        se_fsio_raw_write( io, line, strlen(line) );
    }
    destroy_htiter( hti );

    /* wrote string headers */

    /* write EOL EOL EOT if necessary */
    if( headers->has_inline_data ) {
        char buff[] = {'\n', SEP_EOL, SEP_EOL, SEP_EOT};
        se_fsio_raw_write( io, buff, 4 * sizeof(char) );
    }

    se_fsio_detach( io );

    headers->is_modified = 0;
    headers->is_commited = 1;

    /* append end */
}

static sep_data_t* sep_data_create( void )
{
    sep_data_t* data = (sep_data_t*)malloc( sizeof(*data) );
    memset( data, 0, sizeof(*data) );

    data->can_seek = 1;

    return data;
}

static void sep_data_destroy( sep_data_t* data,
                                      int has_inline_data )
{
    if(data->filename)
        free1char( data->filename );
    if(data->io) {
        if (! has_inline_data)
            se_fsio_close( data->io );
        else
            se_fsio_detach( data->io );
    }

    free(data);
}

static sep_data_t* sep_data_open( const sep_t* sep)
{
    sep_data_t* data = sep_data_create();
    const char* df;
    const char* smode;
    int i, *val;

    /* set data filename */
    if( ! sep->data_only ) {
        if( strcmp(sep->headers->in, "-") == 0 ||
            strcmp(sep->headers->in, "stdin") == 0 ||
            strcmp(sep->headers->in, "stdout") == 0) {
            data->can_seek = 0;
        }
    }
    data->filename = strdup(sep->headers->expanded_in);

    /* set data format */
    df = sep->headers->le ? "native_float" : "xdr_float";

    /* init the fio object used for accessing data */
    if(! sep->data_only && sep->headers->has_inline_data ) {
        data->io = se_fsio_attach( data->filename, sep->headers->fd, df );
        if( data->can_seek ) {
            data->data_offset = sep->headers->raw_size;
        }
    } else {
        if( sep->mode & SEP_READ ) {
            smode = (sep->mode & SEP_WRITE) ? "r+" : "r";
        } else {
            smode = "w+";
        }
        data->io =se_fsio_init( data->filename, smode, df );
    }

    /* set ntraces */
    data->ntraces = 1;
    for (i = 1; i < sep->headers->ndim; ++i) {
        data->ntraces *= sep->headers->n[i];
    }

    /* set trace_size */
    data->trace_size = sep->headers->n[0];

    if(! sep->data_only) {
        /* apply segy_data_offset param */
        val = (int*)ht_get( sep->headers->ht_int, "segy_data_offset" );
        if( val && *val > 0) {
            data->data_offset += *val;
        }

        /* apply segy_header_size param */
        val = (int*)ht_get( sep->headers->ht_int, "segy_header_size" );
        if( val && *val > 0) {
            data->header_size = *val;
        }
    }

    return data;
}

int sep_have_hdr_int(const sep_t* sep, const char* name)
{
    int* val = (int*)ht_get(sep->headers->ht_int, name);
    return (val != NULL);
}

int sep_get_hdr_int(const sep_t* sep, const char* name, int def)
{
    int* val = (int*)ht_get(sep->headers->ht_int, name);
    if(val) return *val;
    return def;
}

int sep_have_hdr_float(const sep_t* sep, const char* name)
{
    double* val = (double*)ht_get(sep->headers->ht_float, name);
    if(val) return 1;
    return sep_have_hdr_int(sep, name);
}

double sep_get_hdr_float(const sep_t* sep, const char* name, double def)
{
    double* val = (double*)ht_get(sep->headers->ht_float, name);
    if(val) return *val;
    if(sep_have_hdr_int(sep, name))
        return (double)sep_get_hdr_int(sep, name, 0);
    return def;
}

int sep_have_hdr(const sep_t* sep, const char* name)
{
    char* sval = (char*)ht_get(sep->headers->ht_str, name);
    if(sval) return 1;
    return sep_have_hdr_float(sep, name);
}

char* sep_get_hdr(const sep_t* sep, const char* name, const char* def)
{
    char*  sval;

    sval = (char*)ht_get(sep->headers->ht_str, name);
    if(sval) {
        return strdup(sval);
    } else {
        float* fval = (float*)ht_get(sep->headers->ht_float, name);
        if(fval) {
            return _sep_print_a_double(*fval, alloc1char(SEP_H_LINESIZE), SEP_H_LINESIZE);
        } else {
            int* ival = (int*)ht_get(sep->headers->ht_int, name);
            if(ival) {
                return asprintf("%d", *ival);
            }
        }
    }

    if(def)
        return strdup(def);
    else
        return NULL;
}

void sep_set_axis(sep_t* sep, int aindex, int n, double o, double d, const char* label)
{
    if(aindex < 0 || aindex >= SEP_NDIM_MAX) {
        ERROR(("Setting SEP file axis: axis index out of range: %d", aindex));
    }
    if(n <= 0) {
        ERROR(("Setting SEP file axis: number of samples is not positive: %d", n));
    }

    sep->headers->n[aindex] = n;
    sep->headers->o[aindex] = o;
    sep->headers->d[aindex] = d;

    if(label != NULL) {
        char* lname = asprintf("label%d", aindex+1);
        sep_set_header(sep, lname, label);
        free(lname);
    }

    if(aindex >= sep->headers->ndim) {
        sep->headers->ndim = aindex + 1;
    }
    sep->headers->is_modified = 1;
}

void sep_insert_axis(sep_t* sep, int aindex, int n, double o, double d, const char* label)
{
    int i;
    if(aindex < 0 || aindex >= SEP_NDIM_MAX) {
        ERROR(("Inserting SEP file axis: axis index out of range: %d", aindex));
    }

    if(sep->headers->ndim < SEP_NDIM_MAX) {
        sep->headers->ndim += 1;
    }

    for(i = sep->headers->ndim - 2; i >= aindex; --i) {
        char* lname = asprintf("label%d", i+1);
        char* lval = sep_get_hdr(sep, lname, NULL);
        sep_set_axis(sep, i+1, sep->headers->n[i], sep->headers->o[i], sep->headers->d[i], lval);
        if(lval) free(lval);
        free(lname);
    }

    sep_set_axis(sep, aindex, n, o , d, label);
}

sep_t* sep_open( const char* path_uri, int mode, int headers_only )
{
    sep_t* sep;
    char* ds_path;

    /* check for invalid arguments */
    if( (mode & SEP_WRITE) &&
        ( strcmp(path_uri, "stdin") == 0 ||
          strcmp(path_uri, "-") == 0) ) {
        ERROR(("Cannot open 'stdin' SEP file for updating."));
    }

    if( (mode & SEP_READ) &&
        (strcmp(path_uri, "stdout") == 0 || strcmp(path_uri, "-") == 0) ) {
        ERROR(("Cannot open 'stdout' SEP file in read mode."));
    }

    /* create the sep_t structure */
    sep = (sep_t*)malloc( sizeof(*sep) );

    memset( sep, 0, sizeof(*sep) );

    sep->mode = mode;
    sep->headers_only = headers_only;

    /* parse path into file and options "dir1/file?noghdr" */
    sep_parse_ds_path( path_uri, &ds_path, &sep->data_only );

    /* open SEP .H file */
    sep->headers = sep_headers_open( sep, ds_path );

    if(! sep->headers_only ) {
        /* open SEP .H@ file */
        sep->data = sep_data_open( sep );
    }

    free1char( ds_path );

    return sep;
}

void sep_close( sep_t* sep )
{
    int has_inline_data = sep->headers->has_inline_data;

    /* close global headers */
    if( (sep->mode & SEP_WRITE) && (! sep->data_only)) {
        /* does something only when is_modified is set */
        sep_headers_append( sep->headers );
    }
    sep_headers_destroy( sep->headers );


    if(! sep->headers_only ) {
        /* close data */
        sep_data_destroy( sep->data, has_inline_data );
    }

    free(sep);
}

void sep_write_headers( sep_t* sep )
{
    /* close global headers */
    if( (sep->mode & SEP_WRITE) && (! sep->data_only)) {
        /* does something only when is_modified is set */
        sep_headers_append( sep->headers );
    }
}

void sep_clear_history( sep_t* sep )
{
    if(sep->headers->history) {
        free(sep->headers->history);
        sep->headers->history = NULL;
    }
}

void sep_copy_headers(sep_t* dst, const sep_t* src)
{
    int i;
    dst->headers->ndim = src->headers->ndim;
    dst->headers->esize = src->headers->esize;
    /*dst->headers->le = src->headers->le;*/
    for(i = 0; i < dst->headers->ndim; ++i) {
        dst->headers->n[i] = src->headers->n[i];
        dst->headers->d[i] = src->headers->d[i];
        dst->headers->o[i] = src->headers->o[i];
    }
    if(src->headers->history != NULL) {
        if(dst->headers->history != NULL)
            free1char(dst->headers->history);
        dst->headers->history = strdup(src->headers->history);
    }
    dst->headers->is_modified = 1;
}

void sep_fill_headers_from_par(sep_t* sep, const char* pargroup)
{
    int i, last_axis_set;
    for(last_axis_set = i = 0; i < SEP_NDIM_MAX; ++i) {
        char pname[10];

        snprintf(pname, sizeof(pname), "n%d", i+1);pname[sizeof(pname)-1]=0;
        if( se_have_namedpar( pargroup, pname, 0 ) ) {
            int n = se_get_namedpar_int32(pargroup, pname);
            if(sep->headers->n[i] != n) {
                sep->headers->n[i] = n;
                sep->headers->is_modified = 1;
            }
            last_axis_set = i;
        }

        snprintf(pname, sizeof(pname), "d%d", i+1);pname[sizeof(pname)-1]=0;
        if( se_have_namedpar( pargroup, pname, 0 ) ) {
            double d = se_get_namedpar_double(pargroup, pname);
            if(!FEQUAL(sep->headers->d[i], d)) {
                sep->headers->d[i] = d;
                sep->headers->is_modified = 1;
            }
            last_axis_set = i;
        }

        snprintf(pname, sizeof(pname), "o%d", i+1);pname[sizeof(pname)-1]=0;
        if( se_have_namedpar( pargroup, pname, 0 ) ) {
            double o = se_get_namedpar_double(pargroup, pname);
            if(!FEQUAL(sep->headers->o[i], o)) {
                sep->headers->o[i] = o;
                sep->headers->is_modified = 1;
            }
            last_axis_set = i;
        }
    }
    if(sep->headers->ndim < last_axis_set+1) {
        sep->headers->ndim = last_axis_set+1;
        sep->headers->is_modified = 1;
    }
}

static void sep_read_part_fromfile(se_fsio*  io, float * buff,
                                      const int n1, const int n2,
                                      const int m1, const int m2, const int m3,
                                      const int k1, const int k2, const int k3)
{
    int64_t skip, offset;
    int64_t esize = sizeof(float);
    int i2, i3;

    skip = k3 *(int64_t)n2 * n1 + k2 * (int64_t)n1 + (int64_t)k1;
    skip *= esize;
    offset = 0;
   se_fsio_skip(io, skip);
    for(i3 = k3; i3 < k3 + m3; ++i3) {
        for(i2 = k2; i2 < k2 + m2; ++i2) {
            se_fsio_read_float(io, buff + offset, (int64_t)m1);
            offset += (int64_t)m1;
            se_fsio_skip(io, esize *(int64_t)(n1 - m1));
        }
    	se_fsio_skip(io, esize *(int64_t)n1 * (n2 - m2));
    }
}



static float* sep_read_fromsep_internal(sep_t* sep, int64_t size, int64_t offset)
{
    float* buff = alloc1float((size_t)size);

    se_fsio_seek( sep->data->io, offset );
   se_fsio_read_float(sep->data->io, buff, (size_t)size);

    return buff;
}

static float* sep_read_fromfile_internal(const char* file, int64_t size, int64_t offset)
{
    sep_t* sep;
    float* buff;

    sep = sep_open( file, SEP_READ, 0 );

    buff = sep_read_fromsep_internal(sep, size, offset);

    sep_close(sep);

    return buff;
}

static float* sep_read_whole_file_internal(const char* file, size_t* n_floats)
{
    sep_t* sep;
    int64_t size;
    size_t bs;
    float* buff;

    sep = sep_open( file, SEP_READ, 0 );
    size =  sep_get_total_size(sep);

    if( sizeof(size_t) < sizeof(int64_t) ) {
        if( size > 512*(int64_t)1024*(int64_t)1024*(int64_t)1024 ) {
            ERROR(("SEP file too big to fit in memory"));
        }
    }

    bs = (size_t)size;

    buff = alloc1float(bs);

    se_fsio_read_float(sep->data->io, buff, bs);

    sep_close(sep);

    if(n_floats) *n_floats = bs;
    return buff;
}
static float* sep_read_part_sep_internal(const sep_t* sep,
                                            size_t* n_floats,
                                            const int m1, const int m2, const int m3,
                                            const int k1, const int k2, const int k3)
{
    int64_t size;
    size_t bs;
    float* buff;
    int n1, n2;


    size = m1 * (int64_t)m2 * m3;

    if( sizeof(size_t) < sizeof(int64_t) ) {
        if( size > 512*(int64_t)1024*(int64_t)1024*(int64_t)1024 ) {
            ERROR(("SEP file too big to fit in memory"));
        }
    }

    bs = (size_t)size;

    buff = alloc1float(bs);
    n1 = sep->headers->n[0];
    n2 = sep->headers->n[1];
    sep_read_part_fromfile(sep->data->io, buff, n1, n2, m1, m2, m3, k1, k2, k3);
    if(n_floats) *n_floats = bs;
    return buff;
}



float* sep_read_fromsep(sep_t* sep, int64_t size, int64_t offset, float** buff, size_t* crt_size)
{
    if(*crt_size < (size_t)size) {
        *crt_size = (size_t)size;
        *buff = realloc1float(*buff, *crt_size);
    }

    se_fsio_seek( sep->data->io, offset );
    se_fsio_read_float(sep->data->io, *buff, (size_t)size);

    return *buff;
}

float* sep_read_fromfile(const char* file, int64_t size, int64_t offset)
{
    return sep_read_fromfile_internal(file, size, offset);
}

float* sep_read_whole_file(const char* file)
{
    return sep_read_whole_file_internal(file, NULL);
}

int interp_rsf3d_from_sep(sep_t *sep, rsf3d_t *rsf) {
    grid3d_t *gridsep = init_grid3d_from_sepheaders(sep);
    grid3d_t subgrid;
    rsf3d_t *subsep;
    int64_t j;
    int ret=0;

    // find covering subgrid
    if (grid3d_find_minimal_covering_subgrid(&rsf->g, gridsep, &subgrid)) {
        ret=1; // no covering subgrid exists
    }

    // read this portion of sep
    subsep=init_rsf3d_part_from_sep_offset(sep, 0,
                                               subgrid.z.n, subgrid.x.n, subgrid.y.n,
                                               axa_closest_indx(&gridsep->z, subgrid.z.o),
                                               axa_closest_indx(&gridsep->x, subgrid.x.o),
                                               axa_closest_indx(&gridsep->y, subgrid.y.o));

    if (grid3d_equal(&subgrid, &rsf->g)) { // grids match: copy the values
        memcpy(rsf->prvt, subsep->prvt, subgrid.z.n*subgrid.x.n*subgrid.y.n*sizeof(float));
    }
    else { // interpolate
        for(j=0;j<rsf->g.y.n;j++) {
            int64_t i;
            double y=axa_x_from_idx(&rsf->g.y,j);

            for(i=0;i<rsf->g.x.n;i++) {
                int64_t k;
                double x=axa_x_from_idx(&rsf->g.x,i);

                for(k=0;k<rsf->g.z.n;k++) {
                    double z=axa_x_from_idx(&rsf->g.z,k);
                    *SEI3(rsf,i,j,k)=(float)se_rsf3d_value_binterp(x,y,z,subsep);
                }
            }
        }
    }

    free(gridsep);
    rsf_destroy3d_and_data_defmem(subsep);
    return ret;
}

float* sep_read_part_sep(const sep_t* sep,
                            const int m1, const int m2, const int m3,
                            const int k1, const int k2, const int k3)
{
    return sep_read_part_sep_internal(sep, NULL, m1, m2, m3, k1, k2, k3);
}

float* sep_read_whole_file_check_values(const char* file,
                                           float min, float max, float replace,
                                           int exit_on_error)
{
    size_t i, n;
    float* vals = sep_read_whole_file_internal(file, &n);
    int warn_printed = 0;
    if(n == 0 || vals == NULL) return vals;
    for(i = 0; i < n; ++i) {
        float v = vals[i];
        char* msg = NULL;
        int err = 0;
        if(!std::isfinite(v)) {
            msg = asprintf("File %s contains non-finite values (first at index %Ld)",
                               file, (long long int)i);
            err = 1;
        } else if(v < min || v > max) {
            msg = asprintf("File %s contains values out of range [%g, %g] (first at index %Ld)",
                               file, min, max, (long long int)i);
            err = 1;
        }

        if(err) {
            if(exit_on_error) {
                ERROR(("%s", msg));
                break;
            } else {
                if(!warn_printed) {
                    WARN(("%s", msg));
                    warn_printed = 1;
                }
                vals[i] = replace;
            }
        }
    }

    return vals;
}

rsf3d_t* init_rsf3d_part_from_sep(const sep_t* sep,
                                          const int m1, const int m2, const int m3,
                                          const int k1, const int k2, const int k3)
{
    grid3d_t* grid;
    rsf3d_t* rsf = NULL;
    float* values;
    int inside_grid;
    int n1 = sep->headers->n[0];
    int n2 = sep->headers->n[1];
    int n3 = sep->headers->n[2];
    int64_t esize =  sizeof(float);
    float* buff;
    int i2,i3;
    int k1r, k2r, k3r, m1r, m2r, m3r;
    int64_t offin, offout;

    if (!((m1 > 0)&&(m2 > 0)&&(m3 > 0))) ERROR(("Dimensions of the extracted box should be positive. "
                                                "Here m1=%d m2=%d m3=%d",m1,m2,m3));
    grid = init_grid3d_part_from_sepheaders(sep, m1, m2, m3, k1, k2, k3);
    inside_grid = !((k1 < 0) || (k1 + m1 > n1) ||
                    (k2 < 0) || (k2 + m2 > n2) ||
                    (k3 < 0) || (k3 + m3 > n3) );

    k1r = maxi(k1, 0);
    m1r = mini(m1 + k1, n1) - k1r;
    if (m1r <= 0) inside_grid = -1;
    k2r = maxi(k2, 0);
    m2r = mini(m2 + k2, n2) - k2r;
    if (m2r <= 0) inside_grid = -1;
    k3r = maxi(k3, 0);
    m3r = mini(m3 + k3, n3) - k3r;
    if (m3r <= 0) inside_grid = -1;

    switch(inside_grid) {
    case -1:
    	values = alloc1float_zero(m1*(int64_t)m2*m3);
        break;
    case 0:
    	values = alloc1float_zero(m1*(int64_t)m2*m3);
        buff = alloc1float(m1r * (int64_t)m2r * m3r);
    	sep_read_part_fromfile(sep->data->io, buff, n1, n2, m1r, m2r, m3r, k1r, k2r, k3r);

    	offin = 0;
    	offout = (int64_t)(k1r - k1) + (k2r - k2) * (int64_t)m1 + (k3r - k3) * (int64_t)m1 * m2;
        for (i3=k3r; i3<k3r + m3r; i3++) {
            for (i2=k2r; i2<k2r + m2r; i2++) {
                ASSERT((offout < m1*(int64_t)m2*m3));
                ASSERT((offin < m1r*(int64_t)m2r*m3r));
                memcpy(values + offout, buff + offin, esize * m1r);
                offin += m1r;
                offout += m1;
            }
            offout += m1 * (int64_t)(m2 - m2r);
        }
        free(buff);
        break;
    case 1:
        values = sep_read_part_sep(sep, m1, m2, m3, k1, k2, k3);
        break;
    default:
    	values = NULL;
        ERROR((" Unrecognized case in function init_rsf3d_part_from_sep"));
    }
    rsf = rsf_create3d_defmem_from_grid(grid, values);
    destroy_grid3d(grid);
    return rsf;
}

float* sep_read_part_sep_offset(const sep_t* sep, off_t off,
                                   const int m1, const int m2, const int m3,
                                   const int k1, const int k2, const int k3) {
    se_fsio_seek(sep->data->io, off);
    return sep_read_part_sep(sep, m1, m2, m3, k1, k2, k3);
}

rsf3d_t* init_rsf3d_part_from_sep_offset(const sep_t* sep, off_t off,
                                                 const int m1, const int m2, const int m3,
                                                 const int k1, const int k2, const int k3) {
    se_fsio_seek(sep->data->io, off);
    return init_rsf3d_part_from_sep(sep, m1, m2, m3, k1, k2, k3);
}

rsf3d_t* init_rsf3d_from_sep(const char* file)
{
    sep_t* sep;
    grid3d_t* grid;
    rsf3d_t* rsf;
    float* values;

    sep = sep_open( file, SEP_READ, 0 );
    grid = init_grid3d_from_sepheaders(sep);
    values = sep_read_whole_file(file);
    rsf = rsf_create3d_defmem_from_grid(grid, values);
    destroy_grid3d(grid);
    sep_close(sep);

    return rsf;
}

rsf2d_t* init_rsf2d_from_sep(const char* file)
{
    sep_t* sep;
    grid2d_t* grid;
    rsf2d_t* rsf;
    float* values;

    sep = sep_open( file, SEP_READ, 0 );
    grid = init_grid2d_from_sepheaders(sep);
    values = sep_read_whole_file(file);
    rsf = rsf_create2d_defmem_from_grid(grid, values);
    destroy_grid2d(grid);
    sep_close(sep);

    return rsf;
}

void sep_set_header_survey_description(sep_t* sep,
                                          const survey_description_t* sd,
                                          const char* prefix,
                                          int use_only_set_fields)
{
    char name[256];
#define _SET_SD_FIELD(p) do { \
        if(!use_only_set_fields || sd->has_##p ) {                      \
            if(prefix) {                                                \
                snprintf(name, sizeof(name), "%s.%s", prefix, #p);      \
            } else {                                                    \
                snprintf(name, sizeof(name), "%s", #p);                 \
            }                                                           \
            sep_set_header_float(sep, name, (float)sd->p);           \
        }                                                               \
    } while(0)

    _SET_SD_FIELD(ox_survey);
    _SET_SD_FIELD(oy_survey);
    _SET_SD_FIELD(survey_azimuth);
    _SET_SD_FIELD(left_handed);
    _SET_SD_FIELD(inline_azimuth);

    _SET_SD_FIELD(first_inline);
    _SET_SD_FIELD(inline_at_origin);
    _SET_SD_FIELD(last_inline);
    _SET_SD_FIELD(inline_increment);
    _SET_SD_FIELD(inline_spacing);

    _SET_SD_FIELD(first_crossline);
    _SET_SD_FIELD(crossline_at_origin);
    _SET_SD_FIELD(last_crossline);
    _SET_SD_FIELD(crossline_increment);
    _SET_SD_FIELD(crossline_spacing);

#undef _SET_SD_FIELD
}


void sep_adjust_survey_description(const sep_t* sep,
                                       survey_description_t* sd)
{
    int i, xaxis = -1, yaxis = -1;
    for(i = 0; i < sep->headers->ndim && (xaxis < 0 || yaxis < 0); ++i) {
        char* lname = asprintf("label%d", i+1);
        char* label = sep_get_hdr(sep, lname, "");

        /* Sometimes there are axis with the same label, or labels
           that we recognize as being X or Y; for the moment we will
           keep only the first axis found to match. This fixes files
           produced by SEP Stack, which when collapsing an axis,
           doesn't remove the label from the higher order axes. */
        if(xaxis < 0 && sep_is_x_axis(label)) {
            xaxis = i;
        } else if(yaxis < 0 && sep_is_y_axis(label)) {
            yaxis = i;
        }
        free(label);
        free(lname);
    }

    if(xaxis < 0) {
        int ndim = sep_get_min_ndim(sep);
        switch(ndim) {
        case 2:
        case 3: xaxis = 1; break;
        case 4: xaxis = 2; break;
        }
    }
    if(yaxis < 0) {
        int ndim = sep_get_min_ndim(sep);
        switch(ndim) {
        case 3: yaxis = 2; break;
        case 4: yaxis = 3; break;
        }
    }

    if(xaxis >= 0) {
        int n    = sep->headers->n[xaxis];
        double d = sep->headers->d[xaxis];
        double o = sep->headers->o[xaxis];
        double ato = sd->crossline_at_origin;
        double sp  = sd->crossline_spacing;
        sd->first_crossline = ato + o / sp;
        sd->crossline_increment = d / sp;
        sd->last_crossline = sd->first_crossline + n*sd->crossline_increment - 1;
    }
    if(yaxis >= 0) {
        int n    = sep->headers->n[yaxis];
        double d = sep->headers->d[yaxis];
        double o = sep->headers->o[yaxis];
        double ato = sd->inline_at_origin;
        double sp  = sd->inline_spacing;
        sd->first_inline = ato + o / sp;
        sd->inline_increment = d / sp;
        sd->last_inline = sd->first_inline + n*sd->inline_increment - 1;
    }
}

void sep_get_header_survey_description(const sep_t* sep,
                                          survey_description_t* sd,
                                          const char* prefix,
                                          const survey_description_t* project_sd,
                                          int adjust_values)
{
    char name[256];

    /* the initial value is comming from project sd or from a project
       file associated with this file */
    if(project_sd) {
        *sd = *project_sd;
    } else {
        if( ! se_get_associated_survey_description( sd, sep->headers->filename, prefix, NULL ) ) {
            se_get_named_survey_description(sd, prefix, NULL);
        }
    }

    /* allow for any value to be overwritten inside this file */

#define _GET_SD_FIELD(p, type) do {                                     \
        if(prefix) {                                                    \
            snprintf(name, sizeof(name), "%s.%s", prefix, #p);          \
        } else {                                                        \
            snprintf(name, sizeof(name), "%s", #p);                     \
        }                                                               \
        if(sep_have_hdr_float(sep, name)) {                          \
            sd->has_##p = 1;                                            \
            sd->p = (type)sep_get_hdr_float(sep, name, sd->p);       \
        }                                                               \
    } while(0)

    _GET_SD_FIELD(ox_survey, double);
    _GET_SD_FIELD(oy_survey, double);
    _GET_SD_FIELD(survey_azimuth, double);
    _GET_SD_FIELD(left_handed, int);
    _GET_SD_FIELD(inline_azimuth, int);

    _GET_SD_FIELD(first_inline, double);
    _GET_SD_FIELD(inline_at_origin, double);
    _GET_SD_FIELD(last_inline, double);
    _GET_SD_FIELD(inline_increment, double);
    _GET_SD_FIELD(inline_spacing, double);

    _GET_SD_FIELD(first_crossline, double);
    _GET_SD_FIELD(crossline_at_origin, double);
    _GET_SD_FIELD(last_crossline, double);
    _GET_SD_FIELD(crossline_increment, double);
    _GET_SD_FIELD(crossline_spacing, double);

#undef _GET_SD_FIELD

    /* finally adjust the survey description to match this file */
    if(adjust_values)
        sep_adjust_survey_description(sep, sd);

    survey_description_update_cached_values(sd);
}

/** At the moment this is not 100% consistent with sep_open */
int sep_is_good_sep_header(const char* file)
{
    int ret = 0;
    if(file && se_file_exists(file) ) {
        se_hash ht = se_parse_pars_file( file, NULL, NULL);
        if(se_have_htpar("n1", ht)) {
            if(se_get_htpar_int32("n1", ht) > 0) {
                if(se_have_htpar("in", ht)) {
                    ret = 1;
                }
            }
        }
        se_destroy_parameters( ht );
    }

    return ret;
}
