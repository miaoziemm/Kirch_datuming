#ifndef SE_RETURN_CODE_H
#define SE_RETURN_CODE_H

/** Defines returns codes commonly used among the library */
typedef enum {
    CODE_SUCCESS           =  0, /**< Evreything went fine - return OK */
    CODE_ERROR        = -1, /**< Generic error code - check the secific function for details */
    CODE_TIMEDOUT     = -2, /**< Used to indicated a time-out condition */
    CODE_TOOMANY      = -3  /**< Indicate an overflow condition,溢出条件 such as a full fixed size queue */
}return_code_t;

#endif