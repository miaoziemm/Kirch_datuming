#ifndef SE_SOCKET_H
#define SE_SOCKET_H
#include <stdlib.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netdb.h>
#include <errno.h>
#include "se_assert.h"
#include "se_log.h"
#include "se_alloc.h"
#include "se_util.h"

/**
 * Public socket handle.
 **/
typedef void* se_socket;

se_socket se_socket_connect( const char* host, int port );
int  se_socket_detach( se_socket s );
void se_socket_close ( se_socket s );

int se_socket_read( const se_socket s, void* buffer, size_t count );
ssize_t se_socket_read_available( const se_socket s,
                                   void* buffer, size_t count );
int se_socket_write( const se_socket s,
                      const void* buffer, size_t count );

int  se_socket_get_keepalive   ( const se_socket s );
int  se_socket_get_recvbuffsize( const se_socket s );
int  se_socket_get_sendbuffsize( const se_socket s );

void se_socket_set_keepalive   ( se_socket s, int enabled );
void se_socket_set_recvbuffsize( se_socket s, int size );
void se_socket_set_sendbuffsize( se_socket s, int size );


/**
 * Server socket
 **/
typedef void* se_ssocket;

se_ssocket se_ssocket_open( const char* address, int port, int backlog );
se_ssocket se_ssocket_open_next_freeport( const char* address,
                                            int* port, int backlog,
                                            int max_tries );
void se_ssocket_close( se_ssocket ss );

se_socket se_ssocket_accept( const se_ssocket ss );

#endif