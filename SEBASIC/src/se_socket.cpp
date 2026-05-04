#include "../include/se_socket.h"

static int  _socket_get_intopt( int socket_fd, int level, int opt,
                                const char* optname );
static void _socket_set_intopt( int socket_fd, int level, int opt,
                                const char* optname, int optval );


/* Socket */

typedef struct
{
    int     fd;
    char*   host;
    int32_t host_ip, local_ip;
    int     port, local_port;
} _socket_t;

se_socket se_socket_connect( const char* host, int port )
{
    _socket_t* s;
    struct hostent* host_ent;
    struct sockaddr_in host_addr, local_host_addr;
    socklen_t sock_len;
    int rc;

    ASSERT(host);
    ASSERT(port > 0);

    /* init the socket_t struct */
    s = (_socket_t*)malloc(sizeof(_socket_t) );
    s->fd = -1;
    s->host = NULL;
    s->host_ip = s->local_ip = 0;
    s->port = s->local_port = -1;

    /* create the socket */
    s->fd = socket( AF_INET, SOCK_STREAM, IPPROTO_TCP);
    if( s->fd < 0 ) {
        int err = errno;
        WARN(("Cannot create TCP socket to '%s:%d': socket() failed: "
              "errno=%d - %s",
              host, port, err, strerror(err)));
        se_socket_close(s);
        return NULL;
    }

    /* resolve the hostname */
    host_ent = gethostbyname( host );
    if(! host_ent) {
        int err = h_errno;
        WARN(("Cannot resolve host '%s': gethostbyname() failed: "
              "h_errno=%d - %s",
              host, err, hstrerror(err)));
        
        se_socket_close(s);
        return NULL;
    }

    /* connect to host */
    memcpy( &s->host_ip, host_ent->h_addr_list[0], host_ent->h_length );

    memset( &host_addr, 0, sizeof(struct sockaddr_in) );
    host_addr.sin_addr.s_addr = s->host_ip;
    host_addr.sin_port = htons( (uint16_t)port );
    host_addr.sin_family = AF_INET;

    rc = connect( s->fd, (struct sockaddr*) &host_addr,
                  sizeof(struct sockaddr_in) );
    if( rc < 0 ) {
        int err = errno;
        WARN(("Cannot connect socket to '%s:%d': connect() failed: "
              "errno=%d - %s",
              host, port, err, strerror(err)));
        se_socket_close(s);
        return NULL;
    }

    /* keep host/local_host info in the socket_t structure */
    s->host = strdup(host);
    s->port = port;
    sock_len = sizeof(struct sockaddr_in);
    rc = getsockname( s->fd, (struct sockaddr*) &local_host_addr, &sock_len );
    if( rc < 0 ) {
        int err = errno;
        WARN(("Cannot connect socket: getsockname() failed: errno=%d - %s",
              err, strerror(err)));
        se_socket_close(s);
        return NULL;
    }
    s->local_port = ntohs( ((struct sockaddr_in*) &local_host_addr)->sin_port );
    s->local_ip   = ((struct sockaddr_in*) &local_host_addr)->sin_addr.s_addr;

    TRACE(("Socket connected to %s:%d.", host, port));
    return (se_socket) s;
}

int se_socket_detach( se_socket _s )
{
    _socket_t* s = (_socket_t*) _s;
    int fd;
    ASSERT(s);

    fd = s->fd;

    if (s->host) free(s->host);
    free(s);

    return fd;
}

void se_socket_close( se_socket _s )
{
    _socket_t* s = (_socket_t*) _s;
    int rc;
    ASSERT(s);

    if( s->fd >= 0 ) {
        do {
            rc = close( s->fd );
        } while (rc < 0 && errno == EINTR);

        if( rc < 0 ) {
            int err = errno;
            WARN(( "Problems closing socket: errno=%d - %s",
                    err, strerror(err) ));
        }
    }

    if (s->host) free(s->host);
    free(s);
}

int se_socket_read( const se_socket _s, void* buffer, size_t count )
{
    _socket_t* s = (_socket_t*) _s;
    ssize_t rc;
    size_t pos;
    byte* bytes = (byte*) buffer;

    ASSERT(s);
    ASSERT(buffer || !count);

    pos = 0;
    while( count > 0 ) {
        do {
            rc = read( s->fd, bytes + pos, count );
        } while( rc < 0 && errno == EINTR );

        if( rc < 0 ) {
            WARN(("Problems reading from socket: errno=%d - %s",
                  errno, strerror(errno)));
            return CODE_ERROR;
        }
        if( rc == 0 ) {
            WARN(("Unexpected end-of-file reading from socket"));
            return CODE_ERROR;
        }

        count -= rc;
        pos   += rc;
    }

    return CODE_SUCCESS;
}

ssize_t se_socket_read_available( const se_socket _s,
                                   void* buffer, size_t count )
{
    _socket_t* s = (_socket_t*) _s;
    ssize_t rc;

    ASSERT(s);
    ASSERT(buffer || !count);

    do {
        rc = read( s->fd, buffer, count );
    } while( rc < 0 && errno == EINTR );

    if( rc < 0 ) {
        WARN(("Problems reading from socket: errno=%d - %s",
              errno, strerror(errno)));
        rc = CODE_ERROR;
    }

    return rc;
}

int se_socket_write( const se_socket _s, const void* buffer, size_t count )
{
    _socket_t* s = (_socket_t*) _s;
    ssize_t rc;
    size_t pos;
    int retry;
    const byte* bytes = (byte*) buffer;

    ASSERT(s);
    ASSERT(buffer || !count);

    pos = 0;
    retry = 0;
    while( count > 0 ) {
        do {
            rc = write( s->fd, bytes + pos, count );
        } while( rc < 0 && errno == EINTR );

        if( rc < 0 ) {
            WARN(("Problems writing to socket: errno=%d - %s",
                  errno, strerror(errno)));
            return CODE_ERROR;
        }

        if( rc == 0 ) {
            if( retry > 10 ) {
                ERROR(("Writing zero bytes to socket after %d trials",
                       retry));
            } else {
                WARN(("Writing zero bytes to socket. Trying again."));
                if( ++retry % 4 == 0 )
                    sleep(1);
            }
        }

        count -= rc;
        pos   += rc;
    }

    return CODE_SUCCESS;
}

int se_socket_get_keepalive( const se_socket _s )
{
    _socket_t* s = (_socket_t*) _s;
    ASSERT(s);

    return _socket_get_intopt( s->fd, SOL_SOCKET, 
                               SO_KEEPALIVE, "SO_KEEPALIVE" );
}

int se_socket_get_recvbuffsize( const se_socket _s )
{
    _socket_t* s = (_socket_t*) _s;
    ASSERT(s);

    return _socket_get_intopt( s->fd, SOL_SOCKET, SO_RCVBUF, "SO_RCVBUF" );
}

int se_socket_get_sendbuffsize( const se_socket _s )
{
    _socket_t* s = (_socket_t*) _s;
    ASSERT(s);

    return _socket_get_intopt( s->fd, SOL_SOCKET, SO_SNDBUF, "SO_SNDBUF" );
}

void se_socket_set_keepalive( se_socket _s, int enabled )
{
    _socket_t* s = (_socket_t*) _s;
    ASSERT(s);

    _socket_set_intopt( s->fd, SOL_SOCKET, SO_KEEPALIVE, "SO_KEEPALIVE",
                        enabled );
}

void se_socket_set_recvbuffsize( se_socket _s, int size )
{
    _socket_t* s = (_socket_t*) _s;
    ASSERT(s);
    ASSERT(size >= 0);

    _socket_set_intopt( s->fd, SOL_SOCKET, SO_RCVBUF, "SO_RCVBUF", size );
}

void se_socket_set_sendbuffsize( se_socket _s, int size )
{
    _socket_t* s = (_socket_t*) _s;
    ASSERT(s);
    ASSERT(size >= 0);

    _socket_set_intopt( s->fd, SOL_SOCKET, SO_SNDBUF, "SO_SNDBUF", size );
}


/* Server socket */

typedef struct
{
    int   fd;
    char* local_address;
    int   local_port;
} _ssocket_t;

static _ssocket_t* _se_ssocket_open( const char* address, int port,
                                      int backlog, int next_freeport )
{
    _ssocket_t* ss;
    struct sockaddr_in host_addr;
    socklen_t sock_len;
    int rc;

    ASSERT(address);
    ASSERT(port > 0);
    ASSERT(backlog > 0);

    /* init the socket_t struct */
    ss = (_ssocket_t*) malloc( sizeof(_ssocket_t) );
    ss->fd = -1;
    ss->local_address = se_strdup(address);
    ss->local_port = -1;

    /* create the server socket */
    ss->fd = socket( AF_INET, SOCK_STREAM, IPPROTO_TCP);
    if( ss->fd < 0 ) {
        WARN(("Cannot create TCP socket: socket() failed: errno=%d - %s",
              errno, strerror(errno)));
        se_ssocket_close(ss);
        return NULL;
    }

    /* bind the socket to port */
    memset( &host_addr, 0, sizeof(struct sockaddr_in) );
    host_addr.sin_addr.s_addr = INADDR_ANY;
    host_addr.sin_port = htons( (uint16_t)port );
    host_addr.sin_family = AF_INET;
    rc = bind( ss->fd, (struct sockaddr*) &host_addr,
               sizeof(struct sockaddr_in) );
    if( rc < 0 ) {
        if(! next_freeport && errno == EADDRINUSE) {
             WARN(("Cannot bind server socket to port %d: bind() failed: "
                  "errno=%d - %s",
                  port, errno, strerror(errno)));
        }
        se_ssocket_close(ss);
        return NULL;
    }

    /* get local port no */
    sock_len = sizeof(struct sockaddr_in);
    rc = getsockname( ss->fd, (struct sockaddr*) &host_addr, &sock_len );
    if( rc < 0 ) {
        WARN(("Cannot create socket: getsockname() failed: errno=%d - %s",
              errno, strerror(errno)));
        se_ssocket_close(ss);
        return NULL;
    }
    ss->local_port = ntohs( ((struct sockaddr_in*) &host_addr)->sin_port );

    /* listen */
    rc = listen( ss->fd, backlog );
    if( rc < 0 ) {
        WARN(("Cannot open server socket on port %d: listen() failed: "
              "errno=%d - %s",
              port, errno, strerror(errno)));
        se_ssocket_close(ss);
        return NULL;
    }

    TRACE(("Server Socket listening at port %d.", ss->local_port));
    return ss;
}

se_ssocket se_ssocket_open( const char* address, int port, int backlog )
{
    return (se_ssocket) _se_ssocket_open( address, port, backlog, 0 );
}

se_ssocket se_ssocket_open_next_freeport( const char* address,
                                            int* port, int backlog,
                                            int max_tries )
{
    _ssocket_t* ss;
    int _port;
    (void)address;
    
    ASSERT(port);
    ASSERT(*port > 0);
    ASSERT(backlog > 0);
    ASSERT(max_tries > 0);

    _port = *port;
    while (1) {
        ss = _se_ssocket_open( "*", _port, backlog, 1 );

        if( ss == NULL && errno == EADDRINUSE && 
            _port < 65535 && max_tries > 0 ) {
            --max_tries;
            TRACE(("Address %s:%d already in use, trying next port",
                   address, _port));
            ++_port;
            continue;
        }

        if (ss)
            *port = _port;
        break;
    }

    return (se_ssocket) ss;
}

void se_ssocket_close( se_ssocket _ss )
{
    _ssocket_t* ss = (_ssocket_t*) _ss;
    int rc;
    ASSERT(ss);

    if( ss->fd >= 0 ) {
        do {
            rc = close( ss->fd );
        } while (rc < 0 && errno == EINTR);
        if (rc < 0) {
            WARN(("Problems closing server socket: errno=%d - %s",
                  errno, strerror(errno)));
        }
    }

    if (ss->local_address) free(ss->local_address);
    free(ss);
}

se_socket se_ssocket_accept( const se_ssocket _ss )
{
    const _ssocket_t* ss = (_ssocket_t*) _ss;
    _socket_t* s;
    struct sockaddr_in host_addr;
    socklen_t sock_len;

    ASSERT(ss);

    /* init the socket_t struct */
    s = (_socket_t*) malloc( sizeof(_socket_t) );
    s->fd = -1;
    s->host = NULL;
    s->host_ip = s->local_ip = 0;
    s->port = s->local_port = -1;

    /* accept a connection */
    sock_len = sizeof(struct sockaddr_in);
    memset( &host_addr, 0, sizeof(struct sockaddr_in) );
    host_addr.sin_addr.s_addr = INADDR_ANY;
    host_addr.sin_port = htons( (uint16_t)(ss->local_port) );
    host_addr.sin_family = AF_INET;
    do {
        s->fd = accept( ss->fd, (struct sockaddr*) &host_addr,
                        (socklen_t*) &sock_len );
    } while (s->fd < 0 && errno == EINTR);
    if( s->fd < 0 ) {
        WARN(("Server socket error: accept() failed: errno=%d - %s",
              errno, strerror(errno) ));
        se_socket_close(s);
        return NULL;
    }

    s->local_port = ss->local_port;
    s->port = ntohs( host_addr.sin_port );
    s->host_ip = host_addr.sin_addr.s_addr;
    TRACE(("Server Socket: client connected."));
    return (se_socket) s;
}


/* internal stuff */
static int _socket_get_intopt( int socket_fd, int level, int opt,
                               const char* optname )
{
    int optval, rc;
    socklen_t opt_len;
    ASSERT(socket_fd);

    optval = 0;
    opt_len = sizeof(int);

    rc = getsockopt( socket_fd, level, opt, (void*) &optval, &opt_len );
    if( rc < 0 ) {
        ERROR(("Error getting socket option '%s': errno=%d - %s",
               optname, errno, strerror(errno)));
    }

    return optval;
}

static void _socket_set_intopt( int socket_fd, int level, int opt,
                                const char* optname, int optval )
{
    int rc;
    ASSERT(socket_fd);
    ASSERT(optname);

    rc = setsockopt( socket_fd, level, opt, (void*) &optval, sizeof(int) );
    if( rc < 0 ) {
        ERROR(("Error setting socket option '%s': errno=%d - %s",
               optname, errno, strerror(errno)));
    }
}
