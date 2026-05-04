#ifndef SE_BASIC_H
#define SE_BASIC_H

#include <iostream>
#include <string>
#include <vector>
#include <map>
#include <memory>
#include <algorithm>
#include <functional>
#include <stdexcept>
#include <cmath>
#include <cstdint>
#include <limits>
#include <cassert>
#include <ccomplex>
#include <cstddef>
#include <cstdlib>
#include <unistd.h>
#include <iostream>
#include <cerrno>
#include <cstring>
#include <errno.h>

#include <sys/time.h>
#include <dirent.h>
#include <pwd.h>


#include "json.hpp"
#include "sqlite3.h"
#include "se_log.h"
#include "se_par.h"
#include "se_version.h"
#include "se_alloc.h"
#include "se_basic_math.h"
#include "se_struct_inter.h"
#include "se_type.h"
#include "se_minmax.h"
#include "se_return_code.h"
#include "se_macro.h"
#include "se_assert.h"
#include "se_coordinate_system.h"
#include "se_grid.h"
#include "se_byte.h"
#include "se_hash.h"
#include "se_array.h"
#include "se_rsf.h"
#include "se_util.h"
#include "se_socket.h"
#include "se_stack.h"
#include "se_exe.h"
#include "se_interpolation.h"
#include "se_timer.h"
#include "se_expression.h"
#include "se_vconv.h"
#include "se_agc.h"
#include "se_bandpass.h"
#include "se_max_factor.h"
#include "se_sort.h"
#include "se_ray.h"
#include "se_gb.h"
#include "se_vel.h"
#include "se_vel_cfg.h"
#include "se_module_par_desc.h"
#include "se_ode_cfg.h"
#include "se_trace_ray.h"
#include "se_trace_gb.h"
#include "se_allpass.h"
#include "se_apfilt.h"
#include "se_quantiler.h"
#include "se_sysinfo.h"
#include "se_geometry.h"
#include "se_array_bit.h"
#include "se_sparse_interpolation.h"
#include "se_horizon_model.h"
#include "se_pick_file.h"
#include "se_math_physics.h"
#include "se_su_fft.h"
#include "se_green_func.h"
#include "se_wavelet.h"
#include "se_thread.h"
#include "se_ray_trace_rk.h"

// Coherence / semblance style utilities
#include "se_coh_func.h"

#include "se_su_par.h"

// 移除冲突的包含，使用统一的复数类型定义
#ifdef __cplusplus
#include <complex>
// 定义 complex 为 std::complex<float>，与项目中的 real 类型保持一致
typedef std::complex<real> complex;
#else
#include <complex.h>
// 对于 C 代码，使用 C99 的复数类型
typedef float complex complex;
#endif

typedef nlohmann::json sejson;

#ifdef SE_USE_OMP
#include <omp.h>
#endif


// Global OpenMP thread auto-config: by default, set threads to maximum cores
// This runs for every C++ translation unit including this header; calling it
// multiple times is harmless. If user explicitly sets OMP_NUM_THREADS, we respect it.
#if defined(__cplusplus) && defined(SE_USE_OMP) && defined(_OPENMP)
namespace se_internal {
static inline void se_init_omp_threads_max() {
	const char* env = std::getenv("OMP_NUM_THREADS");
	if (!env || !*env) {
		omp_set_dynamic(0);
		int n = omp_get_num_procs();
		if (n > 0) omp_set_num_threads(n);
	}
}
struct OmpInitOnce {
	OmpInitOnce() { se_init_omp_threads_max(); }
};
static OmpInitOnce __se_omp_init_once;
} // namespace se_internal
#endif


#endif