#ifndef SE_MACRO_H
#define SE_MACRO_H
/** Safe definition for no operation */
#define NOOP do{;}while(0);

#define CONCAT(prefix, name) prefix ## name

#define CONCAT3(a,b,c) a ## b ## c

#endif