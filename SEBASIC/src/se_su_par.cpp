#include "../include/se_su_par.h"

/* allocate a 1-d array */
static void *alloc1 (size_t n1, size_t size)
{
    void *p;

    if ((p=malloc(n1*size))==NULL)
	return NULL;
    return p;
}




/*getpar.c  Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.			*/
/*********************** self documentation **********************/
/*****************************************************************************
GETPARS - Functions to GET PARameterS from the command line. Numeric
	parameters may be single values or arrays of int, uint,
	short, ushort, long, ulong, float, or double.  Single character
	strings (type string or char *) may also be gotten.
	Arrays of strings, delimited by, but not containing
	commas are permitted.

The functions are:

initargs 	Makes command line args available to subroutines (re-entrant).
		Every par program starts with this call!
getparint		get integers
getparuint		get unsigned integers
getparshort		get short integers
getparushort		get unsigned short integers
getparlong		get long integers
getparulong		get unsigned long integers
getparfloat		get float
getpardouble		get double
getparstring		get a single string
getparstringarray	get string array (fields delimited by commas)
getpar			get parameter by type
getnparint		get n'th occurrence of integer
getnparuint		get n'th occurrence of unsigned int
getnparshort		get n'th occurrence of short integer
getnparushort		get n'th occurrence of unsigned short int
getnparlong		get n'th occurrence of long integer
getnparulong		get n'th occurrence of unsigned long int
getnparfloat		get n'th occurrence of float
getnpardouble		get n'th occurrence of double
getnparstring		get n'th occurrence of string
getnparstringarray	get n'th occurrence of string array
getnpar			get n'th occurrence by type
countparname		return the number of times a parameter names is used
countparval		return the number of values in the last occurrence
				of a parameter
countnparval		return the number of values in the n'th occurrence
				of a parameter
getPar			Promax compatible version of getpar
checkpars()		check the argument list for typos
******************************************************************************
Function Prototypes:
void initargs (int argc, char **argv);
int getparint (char *name, int *p);
int getparuint (char *name, unsigned int *p);
int getparshort (char *name, short *p);
int getparushort (char *name, unsigned short *p);
int getparlong (char *name, long *p);
int getparulong (char *name, unsigned long *p);
int getparfloat (char *name, float *p);
int getpardouble (char *name, double *p);
int getparstring (char *name, char **p);
int getparstringarray (char *name, char **p);
int getnparint (int n, char *name, int *p);
int getnparuint (int n, char *name, unsigned int *p);
int getnparshort (int n, char *name, short *p);
int getnparushort (int n, char *name, unsigned short *p);
int getnparlong (int n, char *name, long *p);
int getnparulong (int n, char *name, unsigned long *p);
int getnparfloat (int n, char *name, float *p);
int getnpardouble (int n, char *name, double *p);
int getnparstring (int n, char *name, char **p);
int getnparstringarray (int n, char *name, char **p);
int getnpar (int n, char *name, char *type, void *ptr);
int countparname (char *name);
int countparval (char *name);
int countnparval (int n, char *name);
void getPar(char *name, char *type, void *ptr);
void checkpars( void );

******************************************************************************
Notes:
Here are some usage examples:

	... if integer n not specified, then default to zero.
	if (!getparint("n", &n)) n = 0;

	... if array of floats vx is specified, then
	if (nx=countparval("vx")) {
		... allocate space for array
		vx = alloc1float(nx);
		... and get the floats
		getparfloat("vx",vx);
	}

The command line for the above examples might look like:
	progname n=35 vx=3.21,4,9.5
	Every par program starts with this call!

More examples are provided in the DTEST code at the end of this file.

The functions: eatoh, eatou, eatol, eatov, eatoi, eatop used
below are versions of atoi that check for overflow.  The source
file for these functions is atopkge.c.

******************************************************************************
Authors:
Rob Clayton & Jon Claerbout, Stanford University, 1979-1985
Shuki Ronen & Jack Cohen, Colorado School of Mines, 1985-1990
Dave Hale, Colorado School of Mines, 05/29/90
Credit to John E. Anderson for re-entrant initargs 03/03/94
*****************************************************************************/
/**************** end self doc ********************************/

/* parameter table */
typedef struct {
    char *name;		/* external name of parameter	*/
    char *asciival;		/* ascii value of parameter	*/
} pointer_table;

/* global variables declared and used internally */
static pointer_table *argtbl;	/* parameter table		*/
static int nargs;		/* number of args that parse	*/
static int tabled = false;	/* true when parameters tabled 	*/
static size_t targc;		/* total number of args		*/
static char **targv;		/* pointer to arg strings	*/
static char *argstr;		/* storage for command line	*/

/* functions declared and used internally */
static int getparindex (int n, char *name);
static void getparinit(void);
static void tabulate (size_t argc, char **argv);
static char *getpfname (void);
static int ccount (char c, char *s);

/*--------------------------------------------------------------------*\
  These variables are used by checkpars() to warn of parameter typos.
  par= does not use getpar() so we need to store that as the first
  parameter name.  lheaders= is buried in fgettr() so we intialize
  that also
  \*--------------------------------------------------------------------*/

#define PAR_NAMES_MAX 512

static char* par_names[PAR_NAMES_MAX];
static int   par_count=0;
static int parcheck = 0;

int xargc;
char **xargv;

/* make command line args available to subroutines -- re-entrant version */
void initargs(int argc, char **argv)
{
    memset( par_names ,0 ,sizeof(par_names) );
    par_names[0] = (char *)"par";
    par_names[1] = (char *)"lheader";
    par_count=2;

    xargc = argc; xargv = argv;
    if(tabled==true){
	free(argstr);
	free(targv);
	free(argtbl);
    }
    tabled =  false;
    return;
}

void strchop(char *s, char *t)
/***********************************************************************
strchop - chop off the tail end of a string "s" after a "," returning
	  the front part of "s" as "t".
	  ************************************************************************
Notes:
Based on strcpy in Kernighan and Ritchie's C [ANSI C] book, p. 106.
************************************************************************
Author: CWP: John Stockwell and Jack K. Cohen, July 1995
***********************************************************************/
{

    while ( (*s != ',') && (*s != '\0') ) {
	*t++ = *s++;
    }
    *t='\0';
}

/* functions to get values for the last occurrence of a parameter name */
int getparint (char *name, int *ptr)
{
    return getnpar(0,name,(char *)"i",ptr);
}

int getparint64 (char *name, int64_t *ptr)
{
    return getnpar(0,name,(char *)"i",ptr);
}

int getparuint (char *name, unsigned int *ptr)
{
    return getnpar(0,name,(char *)"p",ptr);
}
int getparshort (char *name, short *ptr)
{
    return getnpar(0,name,(char *)"h",ptr);
}
int getparushort (char *name, unsigned short *ptr)
{
    return getnpar(0,name,(char *)"u",ptr);
}
int getparlong (char *name, long *ptr)
{
    return getnpar(0,name,(char *)"l",ptr);
}
int getparulong (char *name, unsigned long *ptr)
{
    return getnpar(0,name,(char *)"v",ptr);
}
int getparfloat (char *name, float *ptr)
{
    return getnpar(0,name,(char *)"f",ptr);
}
int getpardouble (char *name, double *ptr)
{
    return getnpar(0,name,(char *)"d",ptr);
}
int getparstring (char *name, char **ptr)
{
    return getnpar(0,name,(char *)"s",ptr);
}
int getparstringarray (char *name, char **ptr)
{
    return getnpar(0,name,(char *)"a",ptr);
}
int getpar (char *name, char *type, void *ptr)
{
    return getnpar(0,name,type,ptr);
}

/* functions to get values for the n'th occurrence of a parameter name */
int getnparint (int n, char *name, int *ptr)
{
    return getnpar(n,name,(char *)"i",ptr);
}
int getnparuint (int n, char *name, unsigned int *ptr)
{
    return getnpar(n,name,(char *)"p",ptr);
}
int getnparshort (int n, char *name, short *ptr)
{
    return getnpar(n,name,(char *)"h",ptr);
}
int getnparushort (int n, char *name, unsigned short *ptr)
{
    return getnpar(n,name,(char *)"u",ptr);
}
int getnparlong (int n, char *name, long *ptr)
{
    return getnpar(n,name,(char *)"l",ptr);
}
int getnparulong (int n, char *name, unsigned long *ptr)
{
    return getnpar(n,name,(char *)"v",ptr);
}
int getnparfloat (int n, char *name, float *ptr)
{
    return getnpar(n,name,(char *)"f",ptr);
}
int getnpardouble (int n, char *name, double *ptr)
{
    return getnpar(n,name,(char *)"d",ptr);
}
int getnparstring (int n, char *name, char **ptr)
{
    return getnpar(n,name,(char *)"s",ptr);
}
int getnparstringarray (int n, char *name, char **ptr)
{
    return getnpar(n,name,(char *)"a",ptr);
}
int getnpar(int n, char *name, char *type, void *ptr)
{
    int i;			/* index of name in symbol table	*/
    int j;		  /* index for par_names[]		*/
    int nval;		/* number of parameter values found	*/
    char *aval;		/* ascii field of symbol		*/

/*--------------------------------------------------------------------*\
  getpar gets called in loops reading traces in some programs.  So
  check for having seen this name before. Also make sure we don't
  walk off the end of the table.
  \*--------------------------------------------------------------------*/

    if( parcheck && strcmp( "lheader" ,name ) ){
	fprintf( stderr ,"getpar() call after checkpars(): %s\n" ,name );
    }

    for( j=0; j<par_count; j++ ){
	if( !strcmp( par_names[j] ,name ) ){
	    break;
	}
    }

    if( j >= par_count && par_count < PAR_NAMES_MAX ){
	par_names[par_count++] = name;
    }

    if(  par_count == PAR_NAMES_MAX ){
	fprintf( stderr, " %s exceeded PAR_NAMES_MAX %d \n" ,xargv[0] ,PAR_NAMES_MAX );
    }

    if (xargc == 1) return 0;
    if (!tabled) getparinit();/* Tabulate command line and parfile */
    i = getparindex(n,name);/* Get parameter index */
    if (i < 0) return 0;	/* Not there */

    if (0 == ptr) {
	fprintf(stderr,"%s: getnpar called with 0 pointer, type = %s\n", __FILE__,type);
    }
	  

    /*
     * handle string type as a special case, since a string
     * may contain commas.
     */
    if (type[0]=='s') {
	*((char**)ptr) = argtbl[i].asciival;
	return 1;
    }

    /* convert vector of ascii values to numeric values */
    for (nval=0,aval=argtbl[i].asciival; *aval; nval++) {
	switch (type[0]) {
	case 'i':
	    *(int*)ptr = eatoi(aval);
	    ptr = (int*)ptr+1;
	    break;
	case 'p':
	    *(unsigned int*)ptr = eatop(aval);
	    ptr = (unsigned int*)ptr+1;
	    break;
	case 'h':
	    *(short*)ptr = eatoh(aval);
	    ptr = (short*)ptr+1;
	    break;
	case 'u':
	    *(unsigned short*)ptr = eatou(aval);
	    ptr = (unsigned short*)ptr+1;
	    break;
	case 'l':
	    *(long*)ptr = eatol(aval);
	    ptr = (long*)ptr+1;
	    break;
	case 'v':
	    *(unsigned long*)ptr = eatov(aval);
	    ptr = (unsigned long*)ptr+1;
	    break;
	case 'f':
	    *(float*)ptr = eatof(aval);
	    ptr = (float*)ptr+1;
	    break;
	case 'd':
	    *(double*)ptr = eatod(aval);
	    ptr = (double*)ptr+1;
	    break;
	case 'a':
	{ char *tmpstr=(char *)"";
		tmpstr = (char *)alloc1(strlen(aval)+1,1);

		strchop(aval,tmpstr);
		*(char**)ptr = tmpstr;
		ptr=(char **)ptr + 1;
	}
	break;
	default:
	    fprintf(stderr,"%s: invalid parameter type = %s",
		    __FILE__,type);
	}
	while (*aval++ != ',') {
	    if (!*aval) break;
	}
    }
    return nval;
}

void checkpars( void ){

    int i;
    int j;
    char buf[256];

    if( getparint( (char *)"verbose" ,&i ) && i == 1 ){

#ifdef SUXDR
	fprintf( stderr ,"Using Big Endian SU data format w/ XDR.\n" );
#else
	fprintf( stderr ,"Using native byte order SU data format w/o XDR.\n" );
#endif

    }

    for( j=1; j<xargc; j++){

	for( i=0; i<par_count; i++ ){
	    sprintf( buf ,"%s=" ,par_names[i] );

	    if( !strncmp( buf ,xargv[j] ,strlen(buf) ) ){
		break;
	    }
	}
	if( i == par_count && strchr( xargv[j] ,'=' ) ){
	    fprintf( stderr ,"Unknown %s argument %s\n" ,xargv[0] ,xargv[j] );
	}

    }

    parcheck = 1;
}

/* return number of occurrences of parameter name */
int countparname (char *name)
{
    int i,nname;

    if (xargc == 1) return 0;
    if (!tabled) getparinit();
    for (i=0,nname=0; i<nargs; ++i)
	if (!strcmp(name,argtbl[i].name)) ++nname;
    return nname;
}

/* return number of values in n'th occurrence of parameter name */
int countnparval (int n, char *name)
{
    int i;

    if (xargc == 1) return 0;
    if (!tabled) getparinit();
    i = getparindex(n,name);
    if (i>=0)
	return ccount(',',argtbl[i].asciival) + 1;
    else
	return 0;
}

/* return number of values in last occurrence of parameter name */
int countparval (char *name)
{
    return countnparval(0,name);
}


/*
 * Return the index of the n'th occurrence of a parameter name,
 * except if n==0, return the index of the last occurrence.
 * Return -1 if the specified occurrence does not exist.
 */
static int getparindex (int n, char *name)
{
    int i;
    if (n==0) {
	for (i=nargs-1; i>=0; --i)
	    if (!strcmp(name,argtbl[i].name)) break;
	return i;
    } else {
	for (i=0; i<nargs; ++i)
	    if (!strcmp(name,argtbl[i].name))
		if (--n==0) break;
	if (i<nargs)
	    return i;
	else
	    return -1;
    }
}

/* Initialize getpar */
static void getparinit (void)
{
    static char *pfname;	/* name of parameter file		*/
    FILE *pffd=NULL;	/* file id of parameter file		*/
    int pflen;		/* length of parameter file in bytes	*/
    int parfile;		/* parfile existence flag		*/
    int argstrlen=0;
    char *pargstr;		/* storage for parameter file args	*/
    int nread=0;		/* bytes fread				*/
    int i, j;		/* counters				*/
    int start = true;
    int debug = false;
    int valencete = false;

    tabled = true;		/* remember table is built		*/


    /* Check if xargc was initiated */

    if(!xargc)
	fprintf(stderr, "%s: xargc=%d -- not initiated in main\n", __FILE__, xargc);

    /* Space needed for command lines */

    for (i = 1, argstrlen = 0; i < xargc; i++) {
	argstrlen += strlen(xargv[i]) + 1;
    }

    /* Get parfile name if there is one */

    if ((pfname = getpfname())) {
	parfile = true;
    } else {
	parfile = false;
    }

    if (parfile) {
	pffd = fopen(pfname, "r");

	/* Get the length */
	fseek(pffd, 0, SEEK_END);

	pflen = ftell(pffd);

	rewind(pffd);
	argstrlen += pflen;
    } else {
	pflen = 0;
    }

/*--------------------------------------------------------------------*\
  Allocate space for command line and parameter file. The pointer
  table could be as large as the string buffer, but no larger.

  The parser logic has been completely rewritten to prevent bad
  input from crashing the program.

  Reginald H. Beardsley			    rhb@acm.org
  \*--------------------------------------------------------------------*/

    argstr = (char *) alloc1(argstrlen+1, 1);
    targv = (char **) alloc1((argstrlen+1)/4,sizeof(char*));

    if (parfile) {
	/* Read the parfile */

	nread = fread(argstr, 1, pflen, pffd);
	if (nread != pflen) {
	    fprintf(stderr,"%s: fread only %d bytes out of %d from %s\n",
		    __FILE__,  nread, pflen, pfname);
	}
	fclose(pffd);


    } 

    /* force input to valid 7 bit ASCII */

    for( i=0; i<nread; i++ ){
	argstr[i] &= 0x7F;
    }

    /* tokenize the input */

    j = 0;

    for( i=0; i<nread; i++ ){

	/* look for start of token */

	if( start ){

	    /* getpars.c:475: warning: subscript has type `char' */
	    if( isgraph( (int)argstr[i] ) ){
		targv[j] = &(argstr[i]);
		start = !start;
		j++;

	    }else{
		argstr[i] = 0;

	    }

	    /* terminate token */

/* getpars.c:487: warning: subscript has type `char' */
	}else if( !valencete && isspace( (int)argstr[i] ) ){
	    argstr[i] = 0;
	    start = !start;
	}

	/* toggle valencete semaphore */

	if( argstr[i] == '\'' || argstr[i] == '\"' ){
	    valencete = !valencete;

	}

    }

    /* display all tokens */

    if( debug ){

	i=0;
	while( i < j && targv[i] != 0 ){
	    if( strlen( targv[i] ) ){
		fprintf( stderr ,"%d -> %s\n" ,i ,targv[i] );
	    }
	    i++;

	}
    }

    /* discard non-parameter tokens */

    i=0;
    targc=0;
    while( i < j && targv[i] != 0 ){
	if( strchr( targv[i] ,'=' ) ){
	    targv[targc] = targv[i];
	    targc++;
	}
	i++;
    }

    /* Copy command line arguments */

    for (j = 1, pargstr = argstr + pflen + 1; j < xargc; j++) {
	strcpy(pargstr,xargv[j]);
	targv[targc++] = pargstr;
	pargstr += strlen(xargv[j]) + 1;
    }

    /* Allocate space for the pointer table */

    argtbl = (pointer_table*) alloc1(targc, sizeof(pointer_table));

    /* Tabulate targv */

    tabulate(targc, targv);

    return;
}
#define PFNAME "par="
/* Get name of parameter file */
static char *getpfname (void)
{
    int i;
    size_t pfnamelen;

    pfnamelen = strlen(PFNAME);
    for (i = xargc-1 ; i > 0 ; i--) {
	if(!strncmp(PFNAME, xargv[i], pfnamelen)
	   && strlen(xargv[i]) != pfnamelen) {
	    return xargv[i] + pfnamelen;
	}
    }
    return NULL;
}

#define iswhite(c)	((c) == ' ' || (c) == '\t' || (c) == '\n')

/* Install symbol table */
static void tabulate (size_t argc, char **argv)
{
    int i;
    char *eqptr;
    int debug=false;

    for (i = 0, nargs = 0 ; i < (int)argc; i++) {
	eqptr = strchr(argv[i], '=');
	if (eqptr) {
	    argtbl[nargs].name = argv[i];
	    argtbl[nargs].asciival = eqptr + 1;
	    *eqptr = (char)0;

	    /* Debugging dump */
	    if( debug ){
		fprintf(stderr,
			"argtbl[%d]: name=%s asciival=%s\n",
			nargs,argtbl[nargs].name,argtbl[nargs].asciival);

	    }
	    nargs++;
	}
    }
    return;
}

/* Count characters in a string */
static int ccount (char c, char *s)
{
    int i, count;
    for (i = 0, count = 0; s[i] != 0; i++)
	if(s[i] == c) count++;
    return count;
}



/* eatoh - convert string s to short integer {SHRT_MIN:SHRT_MAX} */
short eatoh(char *s)
{
    long n = strtol(s, NULL, 10);	
    return (short) n;
}
/* eatou - convert string s to unsigned short integer {0:USHRT_MAX} */
unsigned short eatou(char *s)
{
    unsigned long n = strtoul(s, NULL, 10);
    return (unsigned short) n;
}
/* eatoi - convert string s to integer {INT_MIN:INT_MAX} */
int eatoi(char *s)
{
    long n = strtol(s, NULL, 10);
    return (int) n;
}
/* eatop - convert string s to unsigned integer {0:UINT_MAX} */
unsigned int eatop(char *s)
{
    unsigned long n = strtoul(s, NULL, 10);
    return (unsigned int) n;
}
/* eatol - convert string s to long integer {LONG_MIN:LONG_MAX} */
long eatol(char *s)
{
    long n = strtol(s, NULL, 10);
    return n;
}
/* eatov - convert string s to unsigned long {0:ULONG_MAX} */
unsigned long eatov(char *s)
{
    unsigned long n = strtoul(s, NULL, 10);
    return n;
}
/* eatof - convert string s to float {-FLT_MAX:FLT_MAX} */
float eatof(char *s)
{
    float x = strtod(s, NULL);
    return (float) x;
}
/* eatod - convert string s to double {-DBL_MAX:DBL_MAX} */
double eatod(char *s)
{
    double x = strtod(s, NULL);
    return x;
}

/* Copyright (c) Colorado School of Mines, 2011.*/
/* All rights reserved.                       */
/* errpkge.c
   err	 print warning on application program SE_ERROR and die
   warn print warning on application program SE_ERROR 
   Examples:
   err("Cannot divide %f by %f", x, y);
   warn("fmax = %f exceeds half nyquist= %f", fmax, 0.25/dt);
 
   if (NULL == (fp = fopen(xargv[1], "r")))
   err("can't open %s", xargv[1]);
   ...
   if (-1 == close(fd))
   err("close failed");
*/
void err(char *fmt, ...)
{
    va_list args;

    if (EOF == fflush(stdout)) {
	fprintf(stderr, "\nerr: fflush failed on stdout");
    }
    fprintf(stderr, "\n%s: ", xargv[0]);
    va_start(args,fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
    exit(EXIT_FAILURE);
}


void warn(char *fmt, ...)
{
    va_list args;

    if (EOF == fflush(stdout)) {
	fprintf(stderr, "\nwarn: fflush failed on stdout");
    }
    fprintf(stderr, "\n%s: ", xargv[0]);
    va_start(args,fmt);
    vfprintf(stderr, fmt, args);
    va_end(args);
    fprintf(stderr, "\n");
    return;
}
