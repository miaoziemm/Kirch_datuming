#include "../include/se_util.h"

/* internal helpers */

static void _split_handle_field( const char** field_start,
                                      const char* ch,
                                      char*** fields,
                                      int* field_cap,
                                      int* field_nr,
                                      char sep,
                                      unsigned int flags );

static void _split_strip_quotes( char* field );

static void _util_parse_keyval( char* keyval, se_hash ht_inner );

inline int _split_test_separator(char c, char sep, int flag)
{
    if(flag & SPLIT_USE_ISSPACE) {
        return isspace(c);
    } else {
        return (c == sep);
    }
}


int sqlite3_table_exists(sqlite3 *db, const char* name)
{
      sqlite3_stmt* stmt;
    char* sql = sqlite3_mprintf("SELECT name FROM sqlite_master WHERE type='table' AND name='%q'",
                                name);

    int ret = sqlite3_prepare_v2( db, sql, -1, &stmt, NULL );

    sqlite3_free(sql);
    
    if(ret != SQLITE_OK) return 0;

    ret = sqlite3_step(stmt);
    
    sqlite3_finalize(stmt);

    return (ret == SQLITE_ROW)?1:0;
}

/**
 * Checks whether a string matches a given wildcard pattern. Only does '?'
 * and '*', and multiple patterns separated by '|'.
 * '*' are parsed recursively.
 * '*' does not stands for 'nothing' - it is more like + from reg. exp. :
 *  For example test* matches test1, test2 but doesn't matches test.
 * Later version should take care of this.
 **/
int pattern_match( const char* pattern, const char* string )
{
    ssize_t p, s;
    ssize_t plen, slen;

    plen = strlen(pattern);
    slen = strlen(string);

    for( p = 0; ; ++p ) {
        for( s = 0; ; ++p, ++s ) {
            int s_end = ( s >= slen );
            int p_end = ( p >= plen || pattern[p] == '|' );

            if ( s_end && p_end )
                return 1;

            if ( s_end || p_end )
                break;

            if ( pattern[p] == '?' )
                continue;

            if ( pattern[p] == '*' ) {
                ssize_t i;
                ++p;
                for( i = slen; i >= s; --i ) {
                    if( pattern_match(pattern + p, string + i) )
                        return 1;
                }
                break;
            }

            if ( pattern[p] != string[s] )
                break;
        }
        while( p < plen && pattern[p] != '|') p++;
        if(p >= plen) break;
    }

    return 0;
} /* end of matches */


/**
 * Splits a string into fields using a char separator.
 *
 * Flags: SPLIT_NO_EMPTYFIELDS - don't return empty fields
 *        SPLIT_TRIMFIELDS - trim fields (removes control and white spaces)
 *        SPLIT_QUOTES - don't split the text inside double-quotes \",
 *                       escape sequences (like \") are not supported.
 *
 * For example splitting " aa: bb::c:" with the ':' separator results in:
 * - " aa" " bb" "" "c" "" with SPLIT_NOFLAGS flag,
 * - "aa" "bb" "" "c" "" with SPLIT_TRIMFIELDS flag,
 * - " aa" " bb" "c" with SPLIT_NO_EMPTYFIELDS flags.
 * - "aa" "bb" "c" with SPLIT_TRIMFIELDS and SPLIT_NO_EMPTYFIELDS flags.
 **/

#define _SPLIT_ARRAY_GROWS_BY 8

void split( const char* text,       /* text to be splitted in fields */
                char sep,               /* field separator */
                char*** fields,         /* array of fields (out) */
                int* field_nr,          /* size of the fields array (out) */
                unsigned int flags )    /* flags: see SPLIT_* macros */
{
    const char *field_start = text, *ch;
    int in_quote = 0, field_cap = 0;

    ASSERT(text);
    ASSERT(fields);
    ASSERT(field_nr);

    *fields = NULL;
    *field_nr = 0;

    if (flags & SPLIT_QUOTES) {
        for( ch = text; *ch; ++ch ) {
            if( _split_test_separator(*ch, sep, flags) && ! in_quote ) {
                _split_handle_field( &field_start, ch, fields, &field_cap, field_nr, sep, flags );
            }
            else if( *ch == '"' ) {
                in_quote = ! in_quote;
            }
        }
        _split_handle_field( &field_start, ch, fields, &field_cap, field_nr, sep, flags );
    } else {
        for( ch = text; *ch; ++ch ) {
            if( _split_test_separator(*ch, sep, flags) )
                _split_handle_field( &field_start, ch, fields, &field_cap, field_nr, sep, flags );
        }
        _split_handle_field( &field_start, ch, fields, &field_cap, field_nr, sep, flags );
    }
}

/**
 * Splits a string into fields using a char separator.
 * Same as split() function from above, but returns an array with the fields.
 **/
array split2( const char* text,     /* text to be splitted in fields */
                      char sep,             /* field separator */
                      unsigned int flags )  /* flags: bitmask of SPLIT_* macros */
{
    array array;
    char** fields = NULL;
    int field_nr, i;

    ASSERT(text);

    /* split the text */
    split( text, sep, &fields, &field_nr, flags );

    /* move the fields in the array  */
    array = create_array_args( field_nr );
    for (i = 0; i < field_nr; ++i) {
        array_add( array, fields[i] );
    }

    /* cleanup */
    if( fields )
        free( fields );

    return array;
}

/* split() helper function */
static void _split_strip_quotes( char* field )
{
    int len_old = strlen(field), len_new;

    if( field[0] != '"' || field[0] == '\0' )
        return;
    if( len_old == 1 ) {
        field[0] = '\0';
        return;
    }

    len_new = len_old - 1;
    if( field[len_old - 1] == '"' )
        --len_new;

    memmove( field, field + 1, len_new * sizeof(char) );
    field[len_new] = '\0';
}

/* split() helper function */
static void _split_handle_field( const char** field_start,
                                      const char* ch,
                                      char*** fields,
                                      int* field_cap,
                                      int* field_nr,
                                      char sep,
                                      unsigned int flags )
{
    char *field, **new_fields;
    size_t field_len;

    ASSERT(field_start);
    ASSERT(ch);
    ASSERT(fields);
    ASSERT(field_cap);
    ASSERT(field_nr);

    /* extract the new field */
    field_len = (size_t) (ch - *field_start);
    field = (char *)malloc(field_len + 1);
    memcpy( field, *field_start, field_len * sizeof(char) );
    field[field_len] = '\0';

    *field_start = ch + 1;

    /* apply the flags */
    if( flags & SPLIT_COLLAPSESEP ) {
        const char *l, *r;
        l = field;
        while( *l && _split_test_separator(*l, sep, flags) ) ++l;
        r = field + strlen(field) - 1;
        while( r > l && _split_test_separator(*l, sep, flags)) --r;
        if( l > r ) {
            free(field);
            return;
        } else {
            if(field != l)
                memmove( field, l, (r - l + 1) * sizeof(char) );
            field[r - l + 1] = '\0';
        }
    }
    if( flags & SPLIT_TRIMFIELDS ) {
        str_trim(field);
    }
    if( flags & SPLIT_QUOTES ) {
        _split_strip_quotes(field);
    }
    if( (flags & SPLIT_NO_EMPTYFIELDS) && field[0] == '\0' ) {
        free(field);
        return;
    }

    /* increase the capacity of the fields array if needed */
    if( *field_nr == *field_cap ) {
        new_fields = (char**) malloc( (*field_cap + _SPLIT_ARRAY_GROWS_BY) * sizeof(char*) );
        if( new_fields == NULL ) {
            free(field);
            return;
        }
        memcpy( new_fields, *fields, *field_nr * sizeof(char*) );
        if( *fields )
            free(*fields);
        *fields = new_fields;
        *field_cap += _SPLIT_ARRAY_GROWS_BY;
    }

    /* add the new field into the fields array */
    (*fields)[*field_nr] = field;
    ++(*field_nr);
}


/**
 * Strips the comment from the given line of text. 'ch' is beginning-of-
 * comment character.
 **/
void strip_comment( char* line,          /* line to strip the comment from (in-out) */
                        char comment_ch )    /* beginning of comment char */
{
    char* ch = line;

    ASSERT(line);

    for( ch = line; *ch; ++ch ) {
        if( *ch == comment_ch ) {
            *ch = '\0';
            return;
        }
    }
}

/**
 * Strips the prefix (if found) from the specified string.
 * Returns ture (non-zero) if the prefix is present, or zero otherwise.
 * eg: strip_prefix("error: blah", "error: ") => "blah"
 **/
int strip_prefix( char* str, const char* prefix )
{
    int len, len_prefix;
    ASSERT(str);
    ASSERT(prefix);

    len = strlen(str);
    len_prefix = strlen(prefix);
    if (len < len_prefix)
        return 0;

    if(! memcmp( str, prefix, len_prefix )) {
        memmove( str, str + len_prefix, (len - len_prefix + 1) * sizeof(char) );
        return 1;
    }

    return 0;
}

/**
 * Remove control and white spaces from the beginning and from the end
 * of the specified string.
 **/
void str_trim( char* line )          /* line to trim (in-out) */
{
    const char *l, *r;

    ASSERT(line);

    l = line;
    while( *l ) {
        if( iscntrl(*l) || isspace(*l) )
            ++l;
        else
            break;
    }

    r = line + strlen(line) - 1;
    while( r > l ) {
        if( iscntrl(*r) || isspace(*r) )
            --r;
        else
            break;
    }

    if( l > r )
        line[0] = '\0';
    else {
        if (line != l)
            memmove( line, l, (r - l + 1) * sizeof(char) );
        line[r - l + 1] = '\0';
    }
}


/**
 * Parse the string argument as a signed decimal integer.
 * All chars in the string must be decimal digits, except for the sign char.
 **/
int parse_int( const char* str,  /* string to be parsed */
                   int* val )        /* parsed int value (out) */
{
    long lval = -1;
    char* val_end = NULL;

    ASSERT(str);
    ASSERT(val);

    *val = 0;

    errno = 0;
    lval = strtol( str, &val_end, 10 );

    if( (errno == ERANGE && (lval == LONG_MAX || lval == LONG_MIN)) ||
        (errno != 0 && lval == 0) ||
        ((errno == 0) && (lval < INT_MIN || lval > INT_MAX)) ) {
        return 1; /* error: out of range */
    }
    if( val_end == str ) {
        return 2; /* error: no digits found */
    }
    if( (*str != '+' && *str != '-' && (*str < '0' || *str > '9')) ||
        *val_end != '\0' ) {
        return 3; /* error: other char(s) before or after the number */
    }

    *val = (int) lval;
    return 0; /* success */
}

/**
 * Parse the string argument as a double.
 * All chars in the string must be part of the number.
 **/
int parse_double( const char* str,   /* string to be parsed */
                      double* val )      /* parsed double value (out) */
{
    char *val_end = NULL, *str2 = NULL;
    double dval;
    int rc, len;

    ASSERT(str);
    ASSERT(val);

    *val = 0.0;
    len = strlen(str);
    rc = 0;

    /* transform ".45" > "0.45" */
    if( len >= 2 && str[0] == '.' && isdigit(str[1]) ) {
        str2 = (char *)malloc(len + 2);
        str2[0] = '0';
        strcpy(str2 + 1, str);
        str = str2;
    }

    /* transform "-.32" > "-0.32" */
    if( len >= 3 && (str[0] == '+' || str[0] == '-') && str[1] == '.' && isdigit(str[2]) ) {
        str2 = alloc1char(len + 2);
        str2[0] = str[0];
        str2[1] = '0';
        strcpy(str2 + 2, str + 1);
        str = str2;
    }

    /* now use strtod() to parse the thing */
    errno = 0;
    dval = strtod(str, &val_end);

    if( errno != 0 ) {
        rc = 1; /* error: out of range */
    } else if( (*str != '+' && *str != '-' && (*str < '0' || *str > '9')) ||
             *val_end != '\0' ) {
        rc = 2; /* error: other char(s) before or after the number */
    } else {
        *val = dval;
    }

    if (str2)
        free1char(str2);

    return rc;
}


/**
 * Same as the standard strdup(), but uses alloc1char() to allocate the memory.
 * The caller should call free on the returned string.
 **/
char* se_strdup( const char* str )
{
    char* dupe;
    size_t len;
    ASSERT(str != NULL);

    len = strlen(str);
    dupe = alloc1char( len + 1 );
    memcpy( dupe, str, (len + 1) * sizeof(char) );
    return dupe;
}

/**
 * Same as strndup(), but uses alloc1char() to allocate the memory.
 * The caller should call free on the returned string.
 **/
char* se_strndup( const char* str, size_t n )
{
    char* dupe;
    size_t len;
    ASSERT(str != NULL);
    ASSERT((int)n >= 0);

    len = se_strnlen( str, n );
    dupe = alloc1char( len + 1 );
    memcpy( dupe, str, len * sizeof(char) );
    dupe[len] = '\0';
    return dupe;
}

/**
 * Returns the minimum between the len of the 'str' string and 'n'.
 * Renamed to avoid conflicts with system strnlen
 */
size_t se_strnlen( const char* str, size_t n )
{
    const char* zero;
    ASSERT(str != NULL);
    ASSERT((int)n >= 0);

    zero = (const char*)memchr( str, '\0', n );
    if (zero != NULL)
        return (size_t) (zero - str);
    else
        return n;
}


/**
 * Remove quotes (") from the beginning and end of a string.
 * Eg: ["aa bb"] > [aa bb], [aa bb"] > [aa bb"]
 **/
void str_unquote( char* str )
{
    int len;
    ASSERT(str);

    len = strlen(str);
    if( len >= 2 && str[0] == '"' && str[len - 1] == '"') {
        memmove( str, str + 1, (len - 2) * sizeof(char) );
        str[len - 2] = '\0';
    }
}

/**
 * Returns true (non-zero) if the specified string starts with the specified prefix,
 * or zero otherwise.
 **/
int str_startswith( const char* str, const char* prefix )
{
    ASSERT(str);
    ASSERT(prefix);

    return ! strncmp( str, prefix, strlen(prefix) );
}

/**
 * Returns true (non-zero) if the specified string ends with the specified suffix,
 * or zero otherwise.
 **/
int str_endswith( const char* str, const char* suffix )
{
    int len_str, len_suffix;
    ASSERT(str);
    ASSERT(suffix);

    len_str = strlen(str);
    len_suffix = strlen(suffix);

    if( len_suffix > len_str ) {
        return 0;
    }
    return ! strcmp( str + len_str - len_suffix, suffix );
}

/**
 * Counts how many 'c' chars appears in the 'str' string.
 */
int str_count_char( const char* str, char c )
{
    int count;
    ASSERT(str);

    for (count = 0; *str; ++str) {
        if (*str == c) {
            ++count;
        }
    }
    return count;
}

/**
 * Same as the standard sprintf(), except that it will allocate the
 * right size buffer automatically. The caller should call
 * free1char() on the returned string. The current implementation
 * requires C99 compiler support.
 */
char* asprintf( const char* format, ... )
{
    char* str;

    if( format != NULL ) {
        va_list args;
        int len;

        va_start( args, format );
        len = vsnprintf( NULL, 0, format, args  );
        va_end( args );

        str = alloc1char( len + 1 );

        va_start( args, format );
        vsnprintf( str, len + 1, format, args );
        va_end( args );
    } else {
        str = strdup("(null)");
    }

    return str;
}

int str_indexof(const char* s, char c)
{
    const char* p = s;
    ASSERT(s);

    for (; *p; ++p) {
        if (*p == c) {
            return p - s;
        }
    }
    return -1;
}

/**
 * Returns the index of the last occurence of 'c' in 's' - if any, or -1 otherwise
 */
int str_rindexof(const char* s, char c)
{
    int idx;
    ASSERT(s);

    for (idx = strlen(s) - 1; idx >= 0; --idx) {
        if (s[idx] == c) {
            return idx;
        }
    }
    return -1;
}


/**
 * Free a string array (free all strings then free the array).
 **/
void free_str_array( char** array,
                         int array_size )
{
    int i;
    ASSERT(array_size >= 0);

    if( array ) {
        for (i = 0; i < array_size; ++i) {
            if( array[i] ) {
                free1char(array[i]);
            }
        }
        free(array);
    }
}

char** array_to_null_ended_str_list(array list)
{
    int i, n = array_size(list);
    char** ret;
    if(n < 0) n = 0;
    ret =(char **) malloc((n+1)*sizeof(*ret));
    for(i = 0; i < n; ++i) {
        ret[i] = (char*)array_get_at(list, i);
    }
    ret[n] = NULL;
    return ret;
}

void free_null_ended_str_list(char** list)
{
    if(list) {
        char** ptr = list;
        while(*ptr) {
            free1char(*ptr);
            ++ptr;
        }
        free(list);
    }
}


/**
 * Converts a int32 value to string.
 * The caller should call free on the returned string.
 */
char* int32_to_str(int32_t v)
{
    return asprintf("%d", v);
}

/**
 * Converts a int64 value to string.
 * The caller should call free on the returned string.
 */
char* int64_to_str(int64_t v)
{
    return asprintf("%Ld", v);
}

/**
 * Converts a double value to string.
 * The caller should call free on the returned string.
 */
char* double_to_str(double v)
{
    return asprintf("%f", v);
}



#define _OBJPOOL_GROW_BY 128

typedef struct {
    array            obj_array;      /* object array */
    int                  max_count;      /* max pool capacity */
    process_value_fn destroy_obj;    /* object "destructor" */
    int                  stat_reused;    /* reused objects count */
    int                  stat_max_count; /* max object count */
} _objpool_t;

objpool create_objpool( int max_count, process_value_fn destroy_obj_fn )
{
    _objpool_t* opool = (_objpool_t*)malloc( sizeof(_objpool_t) );
    ASSERT(max_count >= 0);
    ASSERT(destroy_obj_fn);

    opool->obj_array   = create_array_args2( 0, _OBJPOOL_GROW_BY );
    opool->max_count   = max_count;
    opool->destroy_obj = destroy_obj_fn;
    opool->stat_reused    = 0;
    opool->stat_max_count = 0;

    return (objpool) opool;
}

void destroy_objpool( objpool _opool )
{
    _objpool_t* opool = (_objpool_t*) _opool;

    if (opool) {
        array_process( opool->obj_array, opool->destroy_obj );
        destroy_array( opool->obj_array, 0 );
        free( opool );
    }
}

void* objpool_get( objpool _opool )
{
    _objpool_t* opool = (_objpool_t*) _opool;
    int obj_count;
    ASSERT(opool);

    obj_count = array_size(opool->obj_array);
    if (obj_count > 0) {
        opool->stat_reused++;
        return array_remove_at( opool->obj_array, obj_count - 1 );
    } else {
        return NULL;
    }
}

void objpool_put( objpool _opool, void* obj )
{
    _objpool_t* opool = (_objpool_t*) _opool;
    int obj_count;
    ASSERT(opool);

    if (obj == NULL) return;

    obj_count = array_size(opool->obj_array);
    if (obj_count < opool->max_count) {
        array_add( opool->obj_array, obj );
        if (obj_count + 1 > opool->stat_max_count ) {
            opool->stat_max_count = obj_count + 1;
        }
    } else {
        opool->destroy_obj( obj );
    }
}

int objpool_count( const objpool _opool )
{
    _objpool_t* opool = (_objpool_t*) _opool;
    ASSERT(opool);

    return array_size(opool->obj_array);
}

int objpool_max_count( const objpool _opool )
{
    _objpool_t* opool = (_objpool_t*) _opool;
    ASSERT(opool);

    return opool->max_count;
}

void objpool_dump_stats( const objpool _opool, const char* msg )
{
    _objpool_t* opool = (_objpool_t*) _opool;
    ASSERT(opool);
    INFO(("%s object pool statistics: reused %d obj(s), max count: %d obj(s).",
          msg, opool->stat_reused, opool->stat_max_count));
    
}


#define _STRBUF_DEFAULT_GROW_BY 16

typedef struct {
    char* buf;
    int   len, cap, grow_by;
    char  zero;
} _strbuf_t;

/**
 * Creates a new string buffer.
 */
strbuf create_strbuf( void )
{
    return create_strbuf_args( _STRBUF_DEFAULT_GROW_BY );
}

/**
 * Creates a new string buffer.
 * Uses the specified buffer capacity increment.
 */
strbuf create_strbuf_args( int grow_by )
{
    _strbuf_t* sb = (_strbuf_t*)malloc( sizeof(_strbuf_t) );
    ASSERT(grow_by > 0);

    sb->buf = NULL;
    sb->len = 0;
    sb->cap = 0;
    sb->grow_by = grow_by;
    sb->zero = '\0';

    return (strbuf) sb;
}

/**
 * Destroys the specified string buffer.
 */
void destroy_strbuf( strbuf sb )
{
    sb_clear(sb);
    free(sb);
}

static void _sb_append_exactly( _strbuf_t* sb, const char* str, size_t len )
{
    ASSERT(sb);
    ASSERT(str);

    sb_ensure_cap( (strbuf) sb, sb->len + len + 1 );

    memcpy( sb->buf + sb->len, str, len * sizeof(char) );
    sb->len += len;
    sb->buf[sb->len] = '\0';
}

void sb_append_char( strbuf _sb, char c )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    ASSERT(sb);
    sb_ensure_cap( _sb, sb->len + 1 );
    sb->buf[sb->len] = c;
    sb->len += 1;
    sb->buf[sb->len] = '\0';
}

/**
 * Appends the 'str' string to the string buffer.
 */
void sb_append( strbuf _sb, const char* str )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    ASSERT(sb);
    ASSERT(str);

    _sb_append_exactly( sb, str, strlen(str) );
}

/**
 * Appends to the string buffer according to the specified format.
 * The format complies with the format of the *printf functions.
 * The current implementation requires C99 compiler support.
 */
void sb_appendf( strbuf _sb, const char* format, ... )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    va_list args;
    char* str;
    int len;
    ASSERT(sb);
    ASSERT(format);

    va_start( args, format );
    len = vsnprintf( NULL, 0, format, args  );
    va_end( args );

    str = alloc1char( len + 1 );

    va_start( args, format );
    vsnprintf( str, len + 1, format, args );
    va_end( args );

    _sb_append_exactly( sb, str, len );

    free1char(str);
}

/**
 * Appends at most 'n' chars from the 'str' string to the string buffer.
 */
void sb_appendn( strbuf _sb, const char* str, int n )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    ASSERT(sb);
    ASSERT(str);
    ASSERT(n >= 0);

    _sb_append_exactly( sb, str, mini(n, strlen(str)) );
}

static void _sb_insert_exactly( _strbuf_t* sb, size_t index, const char* str, size_t len )
{
    int capacity;
    ASSERT(sb);
    ASSERT(index <= sb->len);
    ASSERT(str);

    capacity = sb->len + len + 1;
    if( capacity > sb->cap ) {
        char* new_buf;

        sb->cap = sb->grow_by * (capacity / sb->grow_by + 1);
        new_buf = alloc1char( sb->cap );

        if (sb->buf != NULL) {
            memcpy( new_buf, sb->buf, index * sizeof(char) );
            memcpy( new_buf + index, str, len * sizeof(char) );
            memcpy( new_buf + index + len, sb->buf + index, (sb->len - index) * sizeof(char) );
        } else {
            memcpy( new_buf, str, len * sizeof(char) );
        }

        if (sb->buf)
            free1char(sb->buf);

        sb->buf = new_buf;
    } else {
        memmove( sb->buf + index + len, sb->buf + index, (sb->len - index) * sizeof(char) );
        memcpy( sb->buf + index, str, len * sizeof(char) );
    }
    sb->len += len;
    sb->buf[sb->len] = '\0';
}

/**
 * Inserts he 'str' string to the string buffer at the specified index.
 */
void sb_insert( strbuf _sb, int index, const char* str )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    ASSERT(sb);
    ASSERT(index >= 0 && index <= sb->len);
    ASSERT(str);

    _sb_insert_exactly( sb, index, str, strlen(str) );
}

/**
 * Inserts to the string buffer at the specified index according to the
 * specified format. The format complies with the format of the *printf
 * functions. The current implementation requires C99 compiler support.
 */
void sb_insertf( strbuf _sb, int index, const char* format, ... )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    va_list args;
    char* str;
    int len;
    ASSERT(sb);
    ASSERT(index >= 0 && index <= sb->len);

    va_start( args, format );
    len = vsnprintf( NULL, 0, format, args  );
    va_end( args );

    str = alloc1char( len + 1 );

    va_start( args, format );
    vsnprintf( str, len + 1, format, args );
    va_end( args );

    _sb_insert_exactly( sb, index, str, len );

    free1char(str);
}

/**
 * Insert at most 'n' chars from the 'str' string to the string buffer
 * at the specified index.
 */
void sb_insertn( strbuf _sb, int index, const char* str, int n )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    ASSERT(sb);
    ASSERT(index >= 0 && index <= sb->len);
    ASSERT(str);
    ASSERT(n >= 0);

    _sb_insert_exactly( sb, index, str, mini(n, strlen(str)) );
}

/**
 * Removed from the string buffer 'length' chars, starting from 'index' index.
 */
void sb_remove( strbuf _sb, int index, int length )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    ASSERT(_sb);
    ASSERT(index >= 0 && index <= sb->len);
    ASSERT(length >= 0 && index + length <= sb->len);

    if (sb->buf != NULL) {
        memmove( sb->buf + index, sb->buf + index + length, (sb->len - index - length) * sizeof(char) );
        sb->len -= length;
        sb->buf[sb->len] = '\0';
    }
}

/**
 * Returns the length of the string buffer.
 */
int sb_length( const strbuf _sb )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    ASSERT(_sb);

    return sb->len;
}

/**
 * Gets a copy of the string from the string buffer if 'str' is not NULL.
 * Gets the length of the string from the string buffer if 'len' is not NULL.
 * This call does not modify the content of the string buffer.
 */
char* sb_copy( const strbuf _sb )
{
    const _strbuf_t* sb = (const _strbuf_t*) _sb;
    ASSERT(sb);

    if (sb->buf != NULL)
        return strndup(sb->buf, sb->len);
    else
        return strdup("");
}

/**
 * Returns a pointer to the (const) string buffer. Use this pointer only till
 * the first modification of the string buffer.
 */
const char* sb_get( const strbuf _sb )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    ASSERT(sb);

    if (sb->buf != NULL)
        return sb->buf;
    else
        return &(sb->zero);
}

/**
 * Release the string from the string buffer.
 * Gets the length of the string from the string buffer if 'len' is not NULL.
 * This call will clear the content of the string buffer.
 */
char* sb_release( strbuf _sb )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    char* str;
    ASSERT(sb);

    if (sb->buf != NULL)
        str = sb->buf;
    else
        str = strdup("");

    sb->buf = NULL;
    sb->len = 0;
    sb->cap = 0;

    return str;
}

/**
 * Clears the content of the specified string buffer.
 */
void sb_clear( strbuf _sb )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    ASSERT(sb);

    if (sb->buf) free1char(sb->buf);

    sb->buf = NULL;
    sb->len = 0;
    sb->cap = 0;
}

/**
 * Returns the capacity of the string buffer.
 */
int sb_cap( const strbuf _sb )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    ASSERT(_sb);

    return sb->cap;
}

/**
 * Ensures that the capacity of the string buffer is at least the specified value.
 */
void sb_ensure_cap( strbuf _sb, int capacity )
{
    _strbuf_t* sb = (_strbuf_t*) _sb;
    ASSERT(_sb);
    ASSERT(capacity >= 0);
    ++capacity;
    if( capacity > sb->cap ) {
        char* new_buf;

        sb->cap = sb->grow_by * (capacity / sb->grow_by + 1);
        new_buf = alloc1char( sb->cap );

        if (sb->buf != NULL)
            memcpy( new_buf, sb->buf, (sb->len + 1) * sizeof(char) );
        else
            new_buf[0] = '\0';

        if (sb->buf)
            free1char(sb->buf);
        sb->buf = new_buf;
    }
}


char* backslash_encode( const char* str, const char* encode_chars, int newlines )
{
    char *result, *q;
    const char *p;
    ASSERT(str);

    /* the encoded string is, at max, double the original string */
    result = alloc1char( strlen(str) * 2 + 1 );

    for (p = str, q = result; *p; ++p, ++q) {

        /* encode "c" to "\\c" */
        if( strchr(encode_chars, *p) || *p == '\\' ) {
            *q++ = '\\';
            *q = *p;
            continue;
        }

        /* encode "\n" to \\n" and "\r" to \\r" */
        if (newlines) {
            if (*p == '\n') {
                *q++ = '\\';
                *q = 'n';
                continue;
            }
            if (*p == '\r') {
                *q++ = '\\';
                *q = 'r';
                continue;
            }
        }

        /* just copy the current char */
        *q = *p;
    }
    *q = '\0';

    return result;
}

char* backslash_decode( const char* str, const char* encode_chars, int newlines )
{
    char *result, *q;
    const char *p, *n;
    ASSERT(str);

    result = alloc1char( strlen(str) + 1 );

    for (p = str, q = result; *p; ++p, ++q) {

        if (*p == '\\') {
            n = p + 1;
            if(! *n)
                break;

            /* decode "\\c" to "c" */
            if ( *n == '\\' || strchr(encode_chars, *n) ) {
                *q = *n;
                p++;
                continue;
            }

            /* decode "\\n" to "\n" and "\\r" to "\r" */
            if( newlines ) {
                if (*n == 'n') {
                    *q = '\n';
                    ++p;
                    continue;
                }
                if (*n == 'r') {
                    *q = '\r';
                    ++p;
                    continue;
                }
            }
        }

        /* just copy the current char */
        *q = *p;
    }
    *q = '\0';

    return result;
}

/**
 * Returns a pointer to the first non backslash-encoded occurrence of the 'ch' character
 * in the 'str' string. If no match if found the function returns NULL.
 * Eg:
 *   strchr_encoded("\\aa", 'a') returns a pointer to the second 'a' character.
 *   strchr_encoded("\\a ", 'a') returns NULL.
 */
char* strchr_encoded( const char* str, char ch )
{
    const char* p;

    /* search for the next non backslash-prefixed 'c' */
    for (p = str; *p; ++p) {
        if (*p == ch && p != str && *(p - 1) != '\\')
            return (char*) p;
    }

    return NULL;
}


#define _HTTOSTR_ENCODE_CHARS "{}= "

/**
 * Creates a string representing the content of a named string se_hash
 * and append it to the string pointed by 'str'.  eg: "se_hash-name
 * {key1=val1 key2=val2}"
 */
int ht_to_str( se_hash ht, const char* name, char** str, int* len )
{
    strbuf sb;
    htiter hti;
    const char* key;
    void* val;
    char *ht_name_encoded, *ht_str, *new_str;
    int ht_len, i, str_len;

    ASSERT(ht);
    ASSERT(name);
    ASSERT(str);
    ASSERT(len);

    /* create a string like "hash_name{key1=val1 key2=val2 ...}"*/
    sb = create_strbuf();
    if (*len > 0)
        sb_append(sb, " ");
    ht_name_encoded = backslash_encode( name, _HTTOSTR_ENCODE_CHARS, 0 );
    sb_append(sb, ht_name_encoded);
    free1char(ht_name_encoded);
    sb_append(sb, "{");
    hti = create_htiter( ht );
    i = 0;
    while( hti_hasnext(hti) ) {
        char *key_encoded, *val_encoded;
        
        hti_next( hti, &key, &val );

        key_encoded = backslash_encode( (const char*)key,
                                       _HTTOSTR_ENCODE_CHARS, 1 );
        val_encoded = backslash_encode( (const char*)val,
                                       _HTTOSTR_ENCODE_CHARS, 1 );

        sb_append( sb, key_encoded );
        sb_append( sb, "=" );
        sb_append( sb, val_encoded );
        if( ++i < ht_size(ht) )
            sb_append( sb, " " );

        free1char( key_encoded );
        free1char( val_encoded );
    }
    destroy_htiter(hti);
    sb_append (sb, "}");
    ht_len = sb_length(sb);
    ht_str = sb_release(sb);
    destroy_strbuf(sb);

    /* append the string to str and update len */
    str_len = *len > 0 ? strlen(*str) : 0;
    if( *len - str_len >= ht_len ) {
        strcpy( *str + str_len, ht_str );
        *len += ht_len;
    } else {
        *len = ht_len + str_len + 1;
        new_str = alloc1char( *len );
        memcpy( new_str, *str, str_len * sizeof(char) );
        memcpy( new_str + str_len, ht_str, (ht_len + 1) * sizeof(char) );
        if (*str)
            free1char(*str);
        *str = new_str;
    }
    free1char(ht_str);

    return 0;
}

/**
 * Parses a string obtained by calling the ht_to_str() function and
 * returns a se_hash with the coresponding values.
 * On error the function returns NULL.
 * eg: "se_hash-name-1 {key1=val1 key2=val2} se_hash-name-2 {key1=val1 key2=val2}"
 */
se_hash ht_from_str( const char* str )
{
    const char *s1, *s2, *ob, *cb, *sp;
    char *ht_name, *ht_name_decoded, *ht_str, *keyval;

    se_hash ht, ht_inner;
    ASSERT(str);

    ht = create_hash();

    /* split by '}' */
    s1 = str;
    while ((cb = strchr_encoded(s1, '}'))) {
        ob = strchr_encoded(s1, '{');
        if(cb < ob) ERROR(("Corrupted hash string:[%s]", str));
        if(! ob) {
            htiter hti = create_htiter(ht);
            while(hti_hasnext(hti)) {
                const char* key;
                hti_next( hti, &key, (void**)&ht_inner);
                destroy_hash_and_keyval(ht_inner, 1, 1);
            }
            destroy_htiter( hti );
            destroy_hash_and_keyval(ht, 1, 0);
            ERROR(("ht_from_str: Missing '{' from [%s]", str));
        }

        ht_inner = create_hash();

        /* parse and decode se_hash name */
        ht_name = strndup( s1, ob - s1 );
        ht_name_decoded = backslash_decode( ht_name, _HTTOSTR_ENCODE_CHARS, 0 );
        /*TRACE(("'%s':", ht_name_decoded ));*/
        free1char( ht_name );

        /* parse se_hash value */
        ht_str = strndup( ob + 1, cb - ob - 1 );

        /* split by ' ' */
        s2 = ht_str;
        while ((sp = strchr_encoded(s2, ' '))) {
            /* parse key & value */
            keyval = strndup( s2, sp - s2 );
            _util_parse_keyval( keyval, ht_inner );

            s2 = sp + 1;
        }
        /* parse the last key & value */
        keyval = strdup( s2 );
        _util_parse_keyval( keyval, ht_inner );

        /* add to se_hash */
        free1char(ht_str);
        ht_put( ht, ht_name_decoded, ht_inner );
        s1 = cb + 1;
        while (*s1 == ' ')
            ++s1;
    }

    return ht;
}

/* ht_from_str() helper function */
static void _util_parse_keyval( char* keyval, se_hash ht_inner )
{
    char *keyval_orig, *keyval_decoded, *equal, *key, *val, *old_val;

    keyval_orig = keyval;
    while (*keyval == ' ')
        ++keyval;

    /* decode 'key=val' */
    keyval_decoded = backslash_decode( keyval, _HTTOSTR_ENCODE_CHARS, 1 );
    free1char( keyval_orig );

    if(! strlen(keyval_decoded)) {
        free1char( keyval_decoded );
        return;
    }

    /* split 'key=val' in 'key' and 'val' */
    equal = strchr( keyval_decoded, '=' );
    if( equal ) {
        key = strndup( keyval_decoded, equal - keyval_decoded );
        val = strdup( equal + 1 );
    } else {
        key = strdup( keyval_decoded );
        val = strdup("");
    }
    free1char( keyval_decoded );

    /* add 'key' & 'val' to the se_hash */
    /*TRACE(("'%s'='%s'", key, val ));*/
    old_val = ht_put( ht_inner, key, val );
    if (old_val) {
        free1char( old_val );
        free1char( key );
    }
}

/**
 * Destroys an se_hash obtained by calling the ht_from_str() function.
 */
void destroy_htht( se_hash ht )
{
    se_hash ht_inner;
    htiter hti;
    const char *key;
    ASSERT(ht);

    hti = create_htiter( ht );
    while( hti_hasnext(hti) ) {
        hti_next( hti, &key, (void**)&ht_inner );

        destroy_hash_and_entries( ht_inner, 1 );
    }
    destroy_htiter( hti );

    destroy_hash_and_keyval( ht, 1, 0 );
}


/**
 * Returns the name of the current UNIX user.
 * The user name is obtained by querying the "USER" environment variable.
 * If the "USER" variable is not defined the "<unkn-user>" string is returned.
 * The caller should call free on the returned string.
 */
char* get_user(void)
{
    char *user = getenv("USER");
    if(! user) {
        WARN(( "Cannot obtain username: getenv(\"USER\") failed" ));
        return strdup("<unkn_user>");
    }
    return strdup(user);
}

static char _forced_host_name[256];
static int _have_forced_host_name = 0;
/** force a hostname - just for this app - useful for debugging and testing */
void set_hostname(const char* hn)
{
    if(hn) {
        strncpy(_forced_host_name, hn, sizeof(_forced_host_name)-1);
        _forced_host_name[sizeof(_forced_host_name)-1] = 0;
        _have_forced_host_name = 1;
    } else {
        _have_forced_host_name = 0;
    }
}

/**
 * Returns the name of the current host.
 * On error the "localhost" string is returned.
 * The caller should call free on the returned string.
 */
char* get_hostname(void)
{
    char host[256];
    int rc;
    if(_have_forced_host_name) {
        return strdup(_forced_host_name);
    }

    rc = gethostname(host, 256);

    if (rc) {
        int err = errno;
        WARN(("Cannot obtain hostname: gethostname() failed: errno=%d - %s",
              err, strerror(err)));
        strcpy(host, "localhost");
    }

    return strdup(host);
}


/**
 * Returns a string representing the current local datetime.
 * The caller should call free on the returned string.
 */
char* get_datetime(void)
{
    return get_adatetime(time(NULL));
}

char* get_adatetime(time_t t)
{
    char* datetime = alloc1char(1024);
    ctime_r( &t, datetime );
    str_trim( datetime );

    return datetime;
}

/**
 * Returns the current working directory.
 * The caller should call free on the returned string.
 */
char* get_cwd(void)
{
    size_t size = 256;
    char* buff = alloc1char(size * sizeof(char));
    for(;;) {
        char* cwd = getcwd(buff, size);
        if(cwd != NULL) {
            return cwd;
        } else {
            if(errno == ERANGE) {
                free(buff);
                size += 64;
                buff = alloc1char(size * sizeof(char));
            } else {
                int err = errno;
                WARN(("Cannot obtain current directory: errno=%d - %s",
                      err, strerror(err)));
                return strdup("/"); /* return root - portability issue */
            }
        }
    }

    /* should never be reached */
    /* VERIFY(0); */
    /* return strdup("/"); */
}

char* get_tmp_dir()
{
    char* tmp = getenv("TMPDIR");
    if( ! tmp ) tmp = getenv("TEMPDIR");

    if( ! tmp ) {
#ifdef P_tmpdir
        tmp = const_cast<char*>(P_tmpdir);
#else
        tmp = "/tmp";
#endif
    }

    return strdup(tmp);
}

/**
 * Generates a random (type 4) Universal Unique IDentifier.
 * Use this function only after init() was called.
 * The caller should call free on the returned string.
 */
char* generate_uuid(void)
{
    static const int  UUID_LENGTH = 16;
    static const int  UUID_STR_LENGTH = 16 * 2 + 4 + 1;
    static const char HEX_CHARS[] = "0123456789abcdef";
    static const byte INDEX_TYPE = 6;
    static const byte INDEX_VARIATION = 8;
    static const byte TYPE_RANDOM_BASED = 4;

    byte uuid[16];
    char *uuid_str, *p;
    int hex, i;
    int num_got = 0, rand = 0, done = 0;

    /* generate 128-bit random number */
    while (! done) {
        for(i = 0; i < 4; ++i) {
            if (num_got == UUID_LENGTH) {
                done = 1;
                break;
            }
            rand = (i == 0 ? (int) random() : rand >> 8);
            uuid[num_got++] = (byte) rand;
        }
    }

    /* set various bits such as type */
    uuid[INDEX_TYPE] &= (byte) 0x0F;
    uuid[INDEX_TYPE] |= (byte) (TYPE_RANDOM_BASED << 4);
    uuid[INDEX_VARIATION] &= (byte) 0x3F;
    uuid[INDEX_VARIATION] |= (byte) 0x80;

    /* convert the byte array into a UUID formatted string */
    uuid_str = alloc1char( UUID_STR_LENGTH );
    for(i = 0, p = uuid_str; i < UUID_LENGTH; ++i) {

        if (i == 4 || i == 6 || i == 8 || i == 10)
            *p++ = '-';

        hex = uuid[i] & 0xFF;
        *p++ = HEX_CHARS[hex >> 4];
        *p++ = HEX_CHARS[hex & 0x0F];
    }
    *p = '\0';

    /* should be freed by the caller */
    return uuid_str;
}

/**
 * Determines if the process with the specified pid exists.
 */
int process_exists(int pid, int* exists)
{
    int rc;
    ASSERT(pid > 0);
    ASSERT(exists);
    *exists = 0;

    /* this only performs error checking, no actual signal is sent to the process */
    rc = kill( pid, 0 );

    if( rc == 0 ) {
        /* process exists */
        *exists = 1;
        return CODE_SUCCESS;
    }
    if( rc == -1 && errno == ESRCH ) {
        /* no such process */
        *exists = 0;
        return CODE_SUCCESS;
    }

    /* error (no permissions ?) */
    return CODE_ERROR;
}


/**
 * Convenient wrapper for the "system" function.
 */
int se_system(const char* command, ...)
{
    va_list args;
    char cmd[2048];
    int ret;

    va_start( args, command );
    vsnprintf( cmd, sizeof(cmd), command, args );
    va_end( args );
    INFOV((10, "About to execute [%s]", cmd));
    ret = system(cmd);
#ifdef _XOPEN_SOURCE
    if( WIFSIGNALED(ret) &&
        ( WTERMSIG(ret) == SIGINT ||
          WTERMSIG(ret) == SIGQUIT )) {
        WARN(("We received SIGINT or SIGQUIT while executing [%s]", cmd));
    }
#endif

    if( ret == -1 ) {
        WARN(("Call to system failed for command [%s]", cmd));
    } else {
        INFOV((20, "Finished [%s] with exit code %d", cmd, ret));
    }
    return ret;
}


/**
 * Conveniently formats a number of bytes (eg 2*1024 => "2KB" )
 * The caller should call free1char on the returned string.
 */
char* format_nbytes( off_t bytes )
{
    if (bytes < (off_t)1024)
        return asprintf( "%d B", (int)bytes );
    if (bytes < (off_t)1024*(off_t)1024)
        return asprintf( "%.3g KB", (double)bytes/1024.0 );
    if (bytes < (off_t)1024*(off_t)1024*(off_t)1024)
        return asprintf( "%.3g MB", (double)bytes/1024.0/1024.0 );
    if (bytes < (off_t)1024*(off_t)1024*(off_t)1024*(off_t)1024)
        return asprintf( "%.3g GB", (double)bytes/1024.0/1024.0/1024.0 );

    return asprintf( "%.3g TB", (double)bytes/1024.0/1024.0/1024.0/1024.0 );
}

/**
 * Conveniently formats a number of seconds (eg 70.2 => "1m10.2s" )
 * The caller should call free1char on the returned string.
 */
char* format_secs( double total_sec )
{
    int m = (int)(total_sec/60) % 60;
    int h = (int)(total_sec/3600) % 24;
    int d = (int)(total_sec/86400) % 7;
    int w = (int)(total_sec/604800);
    double s = total_sec - w*604800 - d*86400 - h*3600 - m*60;

    if (total_sec < 60)
        return asprintf( "%.3fs", s );
    if (total_sec < 3600)
        return asprintf( "%dm%.3fs", m, s );
    if (total_sec < 86400)
        return asprintf( "%dh%dm%.3fs", h, m, s );
    if (total_sec < 604800)
        return asprintf( "%dd%dh%dm%.3fs", d, h, m, s );

    return asprintf( "%dw%dd%dh%dm%.3fs", w, d, h, m, s );
}


char* center_string(const char* str, char fill, int n, int truncate)
{
    char* ret = NULL;

    if(n <= 0) {
        if(truncate) {
            ret = strdup("");
        } else {
            if(str == NULL || *str == 0) {
                ret = strdup("");
            } else {
                ret = strdup(str);
            }
        }
    }
    else if(str == NULL || *str == 0) {
        ret = alloc1char((size_t)n+1);
        memset(ret, fill, (size_t)n);
        ret[(size_t)n] = 0;
    }
    else {
        size_t nend_fill, nbeg_fill;
        size_t len = strlen(str);
        size_t nout = (size_t)n;

        if(len > nout) {
            if(truncate) {
                len = nout;
            } else {
                nout = len;
            }
        }

        ret = alloc1char(nout+1);

        nbeg_fill = (nout - len)/2;
        nend_fill = (nout - len - nbeg_fill);

        memset(ret              , fill, nbeg_fill);
        memcpy(ret+nbeg_fill    , str, len);
        memset(ret+nbeg_fill+len, fill, nend_fill);
        ret[nout] = 0;
    }

    return ret;
}

/**
 * Some systems have stricmp().  Others have strcasecmp().  Because
 * there is no consistency, we will define our own.
 */
int stricmp(const char *a, const char *b)
{
    while( *a!=0 && *b!=0 && tolower(*a) == tolower(*b) ) {
        ++a;
        ++b;
    }

    return tolower(*a) - tolower(*b);
}

int strnicmp(const char *a, const char *b, size_t _n)
{
    ssize_t n = (ssize_t)_n;
    while( n-- > 0 && *a!=0 && *b!=0 && tolower(*a) == tolower(*b) ) {
        a++; 
        b++;
    }

    return n<0 ? 0 : tolower(*a) - tolower(*b);
}

char* human_readable_size(int64_t size)
{
    if( size < 1024 ) {
        return asprintf("%Ld B", (long long int)size);
    } else if(size < (int64_t)1024*(int64_t)1024) {
        return asprintf("%.2f KB", (double)size/1024.0);
    } else if(size < (int64_t)1024*(int64_t)1024*(int64_t)1024) {
        return asprintf("%.2f MB", (double)size/1024.0/1024.0);
    } else if(size < (int64_t)1024*(int64_t)1024*(int64_t)1024*(int64_t)1024) {
        return asprintf("%.2f GB", (double)size/1024.0/1024.0/1024.0);
    } else {
        return asprintf("%.2f TB", (double)size/1024.0/1024.0/1024.0/1024.0);
    }
}

void print_human_readable_size(FILE* f, int64_t size)
{
    char* s = human_readable_size(size);
    fprintf(f, "%s", s);
    free(s);
}


static unsigned char _etoa[] = {
    /*        0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F */

    /* 0 */   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    /* 1 */   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    /* 2 */   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    /* 3 */   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    /* 4 */  32,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0, 46, 60, 40, 43,124,
    /* 5 */  38,  0,  0,  0,  0,  0,  0,  0,  0,  0, 33, 36, 42, 41, 59, 94,
    /* 6 */  45, 47,  0,  0,  0,  0,  0,  0,  0,  0,  0, 44, 37, 95, 62, 63,
    /* 7 */   0,  0,  0,  0,  0,  0,  0,  0,  0, 96, 58, 35, 64, 39, 61, 34,
    /* 8 */   0, 97, 98, 99,100,101,102,103,104,105,  0,  0,  0,  0,  0,  0,
    /* 9 */   0,106,107,108,109,110,111,112,113,114,  0,  0,  0,  0,  0,  0,
    /* A */   0,126,115,116,117,118,119,120,121,122,  0,  0,  0,  0,  0,  0,
    /* B */   0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,  0,
    /* C */ 123, 65, 66, 67, 68, 69, 70, 71, 72, 73,  0,  0,  0,  0,  0,  0,
    /* D */ 125, 74, 75, 76, 77, 78, 79, 80, 81, 82,  0,  0,  0,  0,  0,  0,
    /* E */  92,  0, 83, 84, 85, 86, 87, 88, 89, 90,  0,  0,  0,  0,  0,  0,
    /* F */  48, 49, 50, 51, 52, 53, 54, 55, 56, 57,  0,  0,  0,  0,  0,  0
};


static unsigned char _atoe[] = {
    /*        0   1   2   3   4   5   6   7   8   9   A   B   C   D   E   F */

    /* 0 */   0,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    /* 1 */ 255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,255,
    /* 2 */  64, 90,127,123, 91,108, 80,125, 77, 93, 92, 78,107, 96, 75, 97,
    /* 3 */ 240,241,242,243,244,245,246,247,248,249,122, 94, 76,126,110,111,
    /* 4 */ 124,193,194,195,196,197,198,199,200,201,209,210,211,212,213,214,
    /* 5 */ 215,216,217,226,227,228,229,230,231,232,233,255,224,255, 95,109,
    /* 6 */ 121,129,130,131,132,133,134,135,136,137,145,146,147,148,149,150,
    /* 7 */ 151,152,153,162,163,164,165,166,167,168,169,192, 79,208,161,255
};


void ebcdic2ascii( char *s, size_t len )
{
    while( len-- ) { s[len] = (char)_etoa[ (unsigned char)s[len] ]; }
}


void ascii2ebcdic( char *s, size_t len )
{
    while( len-- ) { s[len] = (char)_atoe[ (unsigned char) s[len] % 0x80 ]; }
}

size_t remove_nans(float* a, size_t n, int copy_last, float replacement)
{
    size_t i, cnt;
    for(cnt = i = 0; i < n; ++i) {
        float v = a[i];
        if( ! std::isfinite(v) || v > 1.0E+20 || v < -1.0E20 ) {
            a[i] = replacement;
            ++cnt;
        } else if( copy_last ) {
            replacement = v;
        }
    }

    return cnt;
}