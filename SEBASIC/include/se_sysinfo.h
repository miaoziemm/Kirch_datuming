#ifndef SE_SYSINFO_H
#define SE_SYSINFO_H

#include <cstdint>
#include <cstddef>
#include "se_type.h"
#include "se_log.h"

/**
 * Returns the number of CPUs as seen by the OS, or -1 on failure.
 * The implementation is OS dependent (there's no POSIX function for this).
 */
int se_get_ncpu( void );

/**
 * Returns the load values, as reported by the kernel. The
 * implementation is system dependent and at this time it works only
 * for Linux.
 *
 * \param[out] l1 will contain the load averaged over the last minute.
 * \param[out] l5 will contain the load averaged over the last 5 minutes.
 * \param[out] l15 will contain the load averaged over the last 15 minutes.
 */
int se_get_load( double* l1, double* l5, double* l15 );

/**
 * Returns memory size in bytes installed on this machine (not
 * including the swap), as reported by the operating system. Although
 * the returned value is in bytes, the operating system usually rounds
 * it to KB.
 *
 * This implementation is OS dependent.
 *
 * \return installed memory in bytes.
 */
int64_t se_sys_get_mem_total(void);

/**
 * Returns total physical memory in bytes. Cross-platform helper used by other
 * modules. Returns 0 on failure.
 */
size_t get_total_memory(void);

int se_get_physical_number_of_cores( void );

#endif