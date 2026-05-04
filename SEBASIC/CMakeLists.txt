# =============================================================================
# SEBASIC Library
# =============================================================================
cmake_minimum_required(VERSION 3.16)

if(NOT DEFINED SE_INSTALL_EXPORT_NAME)
    set(SE_INSTALL_EXPORT_NAME SeismicEngineTargets)
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Header files
set(SEBASIC_HEADERS
    include/se_alloc.h
    include/se_basic.h
    include/se_log.h
    include/se_par.h
    include/se_version.h
    include/json.hpp
    include/se_type.h
    include/se_minmax.h
    include/se_basic_math.h
    include/se_math_physics.h
    include/se_struct_inter.h
    include/se_return_code.h
    include/se_macro.h
    include/se_assert.h
    include/se_coordinate_system.h
    include/se_grid.h
    include/se_byte.h
    include/se_hash.h
    include/se_rsf.h
    include/se_array.h
    include/se_util.h
    include/se_socket.h
    include/sqlite3.h
    include/se_stack.h
    include/se_exe.h
    include/se_interpolation.h
    include/se_timer.h
    include/se_expression.h
    include/se_vconv.h
    include/se_agc.h
    include/se_bandpass.h
    include/se_max_factor.h
    include/se_sort.h
    include/se_vel.h
    include/se_ray.h
    include/se_ray_trace_rk.h
    include/se_wavelet.h
    include/se_green_func.h
    include/se_su_fft.h
    include/se_module_par_desc.h
    include/se_ode_cfg.h
    include/se_trace_ray.h
    include/se_trace_gb.h
    include/se_gb.h
    include/se_allpass.h
    include/se_apfilt.h
    include/se_quantiler.h
    include/se_sysinfo.h
    include/se_geometry.h
    include/se_vel_cfg.h
    include/se_array_bit.h
    include/se_sparse_interpolation.h
    include/se_horizon_model.h
    include/se_pick_file.h
    include/se_coh_func.h
)

# Source files
set(SEBASIC_SOURCES
    src/se_alloc.cpp
    src/se_log.cpp
    src/se_basic_math.cpp
    src/se_math_physics.cpp
    src/se_assert.cpp
    src/se_byte.cpp
    src/se_hash.cpp
    src/se_array.cpp
    src/se_util.cpp
    src/se_socket.cpp
    src/se_stack.cpp
    src/se_exe.cpp
    src/se_interpolation.cpp
    src/se_timer.cpp
    src/se_expression.cpp
    src/se_vconv.cpp
    src/se_agc.cpp
    src/se_bandpass.cpp
    src/se_max_factor.cpp
    src/se_sort.cpp
    src/se_vel.cpp
    src/se_ray.cpp
    src/se_ray_trace2d_rk4.cpp
    src/se_ray_trace3d_rk4.cpp
    src/se_wavelet.cpp
    src/se_green_func.cpp
    src/se_su_fft.cpp
    src/se_module_par_desc.cpp
    src/se_ode_cfg.cpp
    src/se_trace_ray.cpp
    src/se_trace_gb.cpp
    src/se_allpass.cpp
    src/se_apfilt.cpp
    src/se_quantiler.cpp
    src/se_sysinfo.cpp
    src/se_geometry.cpp
    src/se_vel_cfg.cpp
    src/se_array_bit.cpp
    src/se_sparse_interpolation.cpp
    src/se_horizon_model.cpp
    src/se_pick_file.cpp
    src/se_coh_func.cpp
)

# SQLite source file (handled separately to suppress warnings)
set(SQLITE_SOURCES
    src/sqlite3.c
)

# Create static library
add_library(sebasic STATIC ${SEBASIC_SOURCES} ${SQLITE_SOURCES} ${SEBASIC_HEADERS})

# Link required system libraries for SQLite3
find_library(DL_LIBRARY dl)
if(DL_LIBRARY)
    target_link_libraries(sebasic PUBLIC ${DL_LIBRARY})
endif()

# Also ensure pthread support for SQLite3
find_package(Threads)
if(Threads_FOUND)
    target_link_libraries(sebasic PUBLIC Threads::Threads)
endif()

# Suppress warnings for SQLite source file - 兼容老编译器版本的警告选项检查
include(CheckCCompilerFlag)
include(CheckCXXCompilerFlag)
set(SQLITE_WARNING_FLAGS "")

# 基本的警告抑制选项，大多数编译器都支持
set(SQLITE_WARNING_FLAGS "${SQLITE_WARNING_FLAGS} -w")

# 检查 C 编译器支持的警告选项（SQLite3 是 C 代码）
# Only add -Wno-cast-function-type when supported (GCC >= 8 or Clang)
check_c_compiler_flag("-Wno-cast-function-type" HAS_C_WNO_CAST_FUNCTION_TYPE)
if(HAS_C_WNO_CAST_FUNCTION_TYPE)
    if(CMAKE_C_COMPILER_ID STREQUAL "GNU")
        if(CMAKE_C_COMPILER_VERSION VERSION_GREATER_EQUAL "8.0")
            set(SQLITE_WARNING_FLAGS "${SQLITE_WARNING_FLAGS} -Wno-cast-function-type")
        endif()
    elseif(CMAKE_C_COMPILER_ID MATCHES "Clang")
        set(SQLITE_WARNING_FLAGS "${SQLITE_WARNING_FLAGS} -Wno-cast-function-type")
    endif()
endif()

check_c_compiler_flag("-Wno-implicit-fallthrough" HAS_C_WNO_IMPLICIT_FALLTHROUGH)
if(HAS_C_WNO_IMPLICIT_FALLTHROUGH)
    set(SQLITE_WARNING_FLAGS "${SQLITE_WARNING_FLAGS} -Wno-implicit-fallthrough")
endif()

check_c_compiler_flag("-Wno-unused-const-variable" HAS_C_WNO_UNUSED_CONST_VARIABLE)
if(HAS_C_WNO_UNUSED_CONST_VARIABLE)
    set(SQLITE_WARNING_FLAGS "${SQLITE_WARNING_FLAGS} -Wno-unused-const-variable")
endif()

check_c_compiler_flag("-Wno-maybe-uninitialized" HAS_C_WNO_MAYBE_UNINITIALIZED)
if(HAS_C_WNO_MAYBE_UNINITIALIZED)
    set(SQLITE_WARNING_FLAGS "${SQLITE_WARNING_FLAGS} -Wno-maybe-uninitialized")
endif()

check_c_compiler_flag("-Wno-return-local-addr" HAS_C_WNO_RETURN_LOCAL_ADDR)
if(HAS_C_WNO_RETURN_LOCAL_ADDR)
    set(SQLITE_WARNING_FLAGS "${SQLITE_WARNING_FLAGS} -Wno-return-local-addr")
endif()

# 添加针对数组边界警告的抑制（GCC 4.6+）
check_c_compiler_flag("-Wno-array-bounds" HAS_C_WNO_ARRAY_BOUNDS)
if(HAS_C_WNO_ARRAY_BOUNDS)
    set(SQLITE_WARNING_FLAGS "${SQLITE_WARNING_FLAGS} -Wno-array-bounds")
endif()

# 老版本 GCC 的额外兼容性选项
if(CMAKE_COMPILER_IS_GNUCC AND CMAKE_C_COMPILER_VERSION VERSION_LESS "5.0")
    # 对于 GCC 4.x 版本，添加更保守的警告抑制
    set(SQLITE_WARNING_FLAGS "${SQLITE_WARNING_FLAGS} -Wno-unused-variable -Wno-unused-function")
endif()

set_source_files_properties(${SQLITE_SOURCES} PROPERTIES
    COMPILE_FLAGS "${SQLITE_WARNING_FLAGS}"
)

# Set library properties
set_target_properties(sebasic PROPERTIES
    OUTPUT_NAME "sebasic"
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    POSITION_INDEPENDENT_CODE ON
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Public include directories for the library
target_include_directories(sebasic PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS sebasic
    EXPORT ${SE_INSTALL_EXPORT_NAME}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

message(STATUS "SEBASIC library configured")