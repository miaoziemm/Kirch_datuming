#ifndef SE_SU_PAR_H
#define SE_SU_PAR_H


#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <stdbool.h> 
#include <stdarg.h>
#include <errno.h>
#include <time.h>
#include <ctype.h>
#include <math.h>
#include <float.h>
#ifndef __cplusplus
#include <complex.h>
#endif
#include <limits.h>
#include <stdint.h>



/* GLOBAL DECLARATIONS */
extern int xargc; 
extern char **xargv;

/* FUNCTION PROTOTYPES */

#ifdef __cplusplus  /* if C++, specify external C linkage */
extern "C" {
#endif

/* getpar parameter parsing */
void initargs(int argc, char **argv);
int getparint(char *name, int *p);
int getparint64(char *name, int64_t *p);
int getparuint(char *name, unsigned int *p);
int getparshort(char *name, short *p);
int getparushort(char *name, unsigned short *p);
int getparlong(char *name, long *p);
int getparulong(char *name, unsigned long *p);
int getparfloat(char *name, float *p);
int getpardouble(char *name, double *p);
int getparstring(char *name, char **p);
int getparstringarray(char *name, char **p);
int getnparint(int n, char *name, int *p);

int getnparuint(int n, char *name, unsigned int *p);
int getnparshort(int n, char *name, short *p);
int getnparushort(int n, char *name, unsigned short *p);
int getnparlong(int n, char *name, long *p);
int getnparulong(int n, char *name, unsigned long *p);
int getnparfloat(int n, char *name, float *p);
int getnpardouble(int n, char *name, double *p);
int getnparstring(int n, char *name, char **p);
int getnparstringarray(int n, char *name, char **p);
int getnpar(int n, char *name, char *type, void *ptr);
int countparname(char *name);
int countparval(char *name);
int countnparval(int n, char *name);
void checkpars( void );

/* string to numeric conversion with error checking */
short eatoh(char *s);
unsigned short eatou(char *s);
int eatoi(char *s);
unsigned int eatop(char *s);
long eatol(char *s);
unsigned long eatov(char *s);
float eatof(char *s);
double eatod(char *s);

#ifdef __cplusplus  /* if C++ (external C linkage is being specified) */
}
#endif


void err(char *fmt, ...);
void warn(char *fmt, ...);


#endif