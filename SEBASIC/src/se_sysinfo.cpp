#include "../include/se_sysinfo.h"

#ifdef __linux__
    #include <unistd.h>
    #include <fstream>
    #include <sstream>
    #include <string>
    #include <sys/sysinfo.h>
#elif __APPLE__
    #include <unistd.h>
    #include <sys/types.h>
    #include <sys/sysctl.h>
    #include <mach/mach.h>
    #include <mach/vm_statistics.h>
    #include <mach/mach_host.h>
    #include <stdlib.h>
#endif

/**
 * Returns the number of CPUs as seen by the OS, or -1 on failure.
 * The implementation is OS dependent (there's no POSIX function for this).
 */
int se_get_ncpu(void) {
#ifdef __linux__
    return sysconf(_SC_NPROCESSORS_ONLN);
#elif __APPLE__
    int ncpu;
    size_t len = sizeof(ncpu);
    if (sysctlbyname("hw.ncpu", &ncpu, &len, NULL, 0) == 0) {
        return ncpu;
    }
    return -1;
#else
    return -1;  // Unsupported platform
#endif
}

/**
 * Returns the load values, as reported by the kernel. The
 * implementation is system dependent and at this time it works only
 * for Linux.
 *
 * \param[out] l1 will contain the load averaged over the last minute.
 * \param[out] l5 will contain the load averaged over the last 5 minutes.
 * \param[out] l15 will contain the load averaged over the last 15 minutes.
 */
int se_get_load(double* l1, double* l5, double* l15) {
    if (l1 == NULL || l5 == NULL || l15 == NULL) {
        return -1;
    }

#ifdef __linux__
    std::ifstream loadavg_file("/proc/loadavg");
    if (!loadavg_file.is_open()) {
        return -1;
    }
    
    std::string line;
    if (std::getline(loadavg_file, line)) {
        std::istringstream iss(line);
        if (iss >> *l1 >> *l5 >> *l15) {
            return 0;
        }
    }
    return -1;
#elif __APPLE__
    double load[3];
    if (getloadavg(load, 3) == 3) {
        *l1 = load[0];
        *l5 = load[1];
        *l15 = load[2];
        return 0;
    }
    return -1;
#else
    return -1;  // Unsupported platform
#endif
}

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
int64_t se_sys_get_mem_total(void) {
#ifdef __linux__
    struct sysinfo info;
    if (sysinfo(&info) == 0) {
        return (int64_t)info.totalram * info.mem_unit;
    }
    return -1;
#elif __APPLE__
    int64_t mem_size;
    size_t len = sizeof(mem_size);
    if (sysctlbyname("hw.memsize", &mem_size, &len, NULL, 0) == 0) {
        return mem_size;
    }
    return -1;
#else
    return -1;  // Unsupported platform
#endif
}

/**
 * Returns the physical number of CPU cores.
 */
int se_get_physical_number_of_cores(void) {
#ifdef __linux__
    // Try to read from /proc/cpuinfo
    std::ifstream cpuinfo("/proc/cpuinfo");
    if (!cpuinfo.is_open()) {
        // Fallback to logical CPUs
        return se_get_ncpu();
    }
    
    int physical_cores = 0;
    int current_physical_id = -1;
    int current_core_id = -1;
    std::string line;
    
    while (std::getline(cpuinfo, line)) {
        if (line.find("physical id") == 0) {
            std::size_t colon_pos = line.find(':');
            if (colon_pos != std::string::npos) {
                current_physical_id = std::stoi(line.substr(colon_pos + 1));
            }
        } else if (line.find("core id") == 0) {
            std::size_t colon_pos = line.find(':');
            if (colon_pos != std::string::npos) {
                current_core_id = std::stoi(line.substr(colon_pos + 1));
            }
        } else if (line.find("processor") == 0 && current_physical_id >= 0 && current_core_id >= 0) {
            // Count unique physical_id + core_id combinations
            physical_cores++;
            current_physical_id = -1;
            current_core_id = -1;
        }
    }
    
    // If we couldn't parse properly, fallback to logical CPUs
    return physical_cores > 0 ? physical_cores : se_get_ncpu();
    
#elif __APPLE__
    int physical_cores;
    size_t len = sizeof(physical_cores);
    if (sysctlbyname("hw.physicalcpu", &physical_cores, &len, NULL, 0) == 0) {
        return physical_cores;
    }
    // Fallback to logical CPUs
    return se_get_ncpu();
#else
    return -1;  // Unsupported platform
#endif
}

// Cross-platform helper to return total physical memory in bytes.
// Returns 0 on failure.
size_t get_total_memory() {
#if defined(_WIN32) || defined(_WIN64)
    MEMORYSTATUSEX status;
    status.dwLength = sizeof(status);
    if (GlobalMemoryStatusEx(&status)) {
        return (size_t)status.ullTotalPhys;
    } else {
        fprintf(stderr, "GlobalMemoryStatusEx failed\n");
        return 0;
    }
#elif defined(__APPLE__) || defined(__MACH__)
    int mib[2];
    uint64_t mem;
    size_t len = sizeof(mem);

    mib[0] = CTL_HW;
    mib[1] = HW_MEMSIZE;

    if (sysctl(mib, 2, &mem, &len, NULL, 0) == -1) {
        perror("sysctl");
        return 0;
    }

    return (size_t)mem;
#elif defined(__linux__)
    FILE *file = fopen("/proc/meminfo", "r");
    if (file == NULL) {
        perror("fopen");
        return 0;
    }

    char line[256];
    size_t mem_total_kb = 0;

    while (fgets(line, sizeof(line), file)) {
        if (sscanf(line, "MemTotal: %zu kB", &mem_total_kb) == 1) {
            mem_total_kb *= 1024; // convert to bytes
            break;
        }
    }

    fclose(file);
    return mem_total_kb;
#else
    fprintf(stderr, "Unsupported platform\n");
    return 0;
#endif
}