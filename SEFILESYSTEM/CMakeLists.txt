# =============================================================================
# SEFILESYSTEM Library
# =============================================================================
cmake_minimum_required(VERSION 3.16)

if(NOT DEFINED SE_INSTALL_EXPORT_NAME)
    set(SE_INSTALL_EXPORT_NAME SeismicEngineTargets)
endif()

# Include directories
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Header files
set(SEFS_HEADERS
    include/se_coordinate_system_par.h
    include/se_fs.h
    include/se_fs_io.h
    include/se_fs_segy.h
    include/se_fs_segy_bin_types.h
    include/se_fs_segy_header.h
    include/se_fs_segy_index.h
    include/se_fs_segy_trchead_eval.h
    include/se_fs_sep.h
    include/se_fs_sep_fake.h
    include/se_fs_struct.h
    include/se_io.h
    include/se_jcs.h
    include/se_module_par_desc_sep.h
    include/se_par_sep.h
    include/se_uri.h
    include/sesi_map.h
    include/se_par_getset.h
    include/se_vel_cfg_sep.h
    include/se_ode_cfg_sep.h
)

# Source files
set(SEFS_SOURCES
    src/se_coordinate_system_par.cpp
    src/se_fs_io.cpp
    src/se_fs_segy_bin_types.cpp
    src/se_fs_segy_data.cpp
    src/se_fs_segy_header.cpp
    src/se_fs_segy_index.cpp
    src/se_fs_segy_io.cpp
    src/se_fs_segy_trchead_eval.cpp
    src/se_fs_segy_util.cpp
    src/se_fs_sep.cpp
    src/se_io.cpp
    src/se_jcs.cpp
    src/se_module_par_desc_sep.cpp
    src/se_par_sep.cpp
    src/se_uri.cpp
    src/sesi_map.cpp
    src/se_par_getset.cpp
    src/se_vel_cfg_sep.cpp
    src/se_ode_cfg_sep.cpp
)

# Create static library
add_library(sefs STATIC ${SEFS_SOURCES} ${SEFS_HEADERS})

# Set library properties
set_target_properties(sefs PROPERTIES
    OUTPUT_NAME "sefs"
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    VERSION ${PROJECT_VERSION}
    SOVERSION ${PROJECT_VERSION_MAJOR}
)

# Public include directories for the library
target_include_directories(sefs PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)

install(TARGETS sefs
    EXPORT ${SE_INSTALL_EXPORT_NAME}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

message(STATUS "SEFILESYSTEM library configured")