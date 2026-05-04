# =============================================================================
# Seismic Engine Project - Unified CMake Configuration
# =============================================================================
cmake_minimum_required(VERSION 3.10)

# Set policies for modern CMake behavior
if(POLICY CMP0074)
    cmake_policy(SET CMP0074 NEW)  # Use <PackageName>_ROOT variables
endif()
if(POLICY CMP0077)
    cmake_policy(SET CMP0077 NEW)  # option() honors normal variables
endif()
if(POLICY CMP0079)
    cmake_policy(SET CMP0079 NEW)  # Allow linking to targets in other directories
endif()

# Project information
project(SeismicEngine 
    VERSION 1.0.0 
    DESCRIPTION "Seismic data processing and analysis engine"
    LANGUAGES C CXX
)

# Standard install directory layout
include(GNUInstallDirs)
set(SE_INSTALL_EXPORT_NAME SeismicEngineTargets)

# Configure RPATH so installed binaries find bundled libraries
set(CMAKE_SKIP_BUILD_RPATH FALSE)
set(CMAKE_BUILD_WITH_INSTALL_RPATH FALSE)
set(CMAKE_INSTALL_RPATH_USE_LINK_PATH TRUE)
if(APPLE)
    set(CMAKE_INSTALL_RPATH "@loader_path/../${CMAKE_INSTALL_LIBDIR}")
else()
    set(CMAKE_INSTALL_RPATH "$ORIGIN/../${CMAKE_INSTALL_LIBDIR}")
endif()

# =============================================================================
# Global Build Configuration
# =============================================================================

# Check compiler support for C++17
include(CheckCXXCompilerFlag)
check_cxx_compiler_flag("-std=c++17" COMPILER_SUPPORTS_CXX17)
check_cxx_compiler_flag("-std=c++1z" COMPILER_SUPPORTS_CXX1Z)

# Set C and C++ standards with fallback options
set(CMAKE_C_STANDARD 99)
set(CMAKE_C_STANDARD_REQUIRED ON)

if(COMPILER_SUPPORTS_CXX17)
    set(CMAKE_CXX_STANDARD 17)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    message(STATUS "C++17 support detected")
elseif(COMPILER_SUPPORTS_CXX1Z)
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++1z")
    message(STATUS "Using -std=c++1z for C++17 support")
else()
    # Fallback to C++14 if C++17 is not available
    set(CMAKE_CXX_STANDARD 14)
    set(CMAKE_CXX_STANDARD_REQUIRED ON)
    message(WARNING "C++17 not supported, falling back to C++14")
endif()

# Disable compiler extensions for better portability
set(CMAKE_CXX_EXTENSIONS OFF)

# Set default build type
if(NOT CMAKE_BUILD_TYPE)
    set(CMAKE_BUILD_TYPE Release)
endif()

# =============================================================================
# Build Options
# =============================================================================
option(SE_USE_OMP "Enable OpenMP support" ON)
option(SE_USE_MPI "Enable MPI support and build MPI programs" ON)
option(BUILD_EXAMPLES "Build example programs" ON)
option(SE_CUSTOM_MPI_PATH "Use custom MPI installation path" OFF)
option(SE_CUSTOM_OMP_PATH "Use custom OpenMP installation path" OFF)
option(SE_BUILD_GUI "Build GUI components (if any)" ON)
option(SE_BUNDLE_DEPENDENCIES "Bundle shared library dependencies with installed binaries" ON)
option(SE_USE_EIGEN "Enable Eigen-based coherence kernels" OFF)

# Custom paths
set(SE_MPI_ROOT "" CACHE PATH "Custom MPI installation root directory")
set(SE_OMP_ROOT "" CACHE PATH "Custom OpenMP installation root directory")

# Register runtime executables for post-install dependency bundling
function(se_register_runtime_target target_name)
    if(NOT TARGET ${target_name})
        message(FATAL_ERROR "Runtime target ${target_name} must be defined before registration")
    endif()
    install(TARGETS ${target_name}
        EXPORT ${SE_INSTALL_EXPORT_NAME}
        RUNTIME DESTINATION ${CMAKE_INSTALL_BINDIR}
        BUNDLE DESTINATION ${CMAKE_INSTALL_BINDIR}
    )
    set_property(GLOBAL APPEND PROPERTY SE_RUNTIME_TARGETS ${target_name})
endfunction()

# Discover Qt plugin directory for GUI deployments when available
set(SE_QT_PLUGIN_DIR "" CACHE PATH "Qt plugin directory used to bundle platform plugins")
find_program(QTPATHS_EXECUTABLE NAMES qtpaths6 qtpaths-qt6 qtpaths5 qtpaths-qt5 qtpaths
    DOC "Helper executable to query Qt installation paths")
if(NOT SE_QT_PLUGIN_DIR AND QTPATHS_EXECUTABLE)
    execute_process(
        COMMAND ${QTPATHS_EXECUTABLE} --plugin-dir
        OUTPUT_VARIABLE _autodetected_qt_plugin_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE
        ERROR_QUIET
        RESULT_VARIABLE _qtpaths_result
    )
    if(_qtpaths_result EQUAL 0 AND _autodetected_qt_plugin_dir)
        set(SE_QT_PLUGIN_DIR "${_autodetected_qt_plugin_dir}" CACHE PATH
            "Qt plugin directory used to bundle platform plugins" FORCE)
    endif()
endif()

# =============================================================================
# Global Compiler Settings
# =============================================================================

# Base compiler flags
set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -Wall -Wextra")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -Wextra")

# Debug flags
set(CMAKE_C_FLAGS_DEBUG "-g -O0 -DDEBUG")
set(CMAKE_CXX_FLAGS_DEBUG "-g -O0 -DDEBUG")

# Release flags
set(CMAKE_C_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")

# =============================================================================
# Output Directory Configuration
# =============================================================================
set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/bin)
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)
set(CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# =============================================================================
# Global Include Directories
# =============================================================================
include_directories(${CMAKE_CURRENT_SOURCE_DIR})

# =============================================================================
# OpenMP Configuration
# =============================================================================
if(SE_USE_OMP)
    # Set custom OpenMP search path if provided
    if(SE_CUSTOM_OMP_PATH AND SE_OMP_ROOT)
        message(STATUS "Using custom OpenMP path: ${SE_OMP_ROOT}")
        set(OpenMP_ROOT ${SE_OMP_ROOT})
        
        # Try to find specific OpenMP library files
        find_library(OpenMP_C_LIBRARY 
            NAMES omp gomp iomp5
            PATHS ${SE_OMP_ROOT}/lib ${SE_OMP_ROOT}/lib64
            NO_DEFAULT_PATH
        )
        find_library(OpenMP_CXX_LIBRARY 
            NAMES omp gomp iomp5
            PATHS ${SE_OMP_ROOT}/lib ${SE_OMP_ROOT}/lib64
            NO_DEFAULT_PATH
        )
        
        # Set include directories
        if(EXISTS ${SE_OMP_ROOT}/include)
            set(OpenMP_C_INCLUDE_DIR ${SE_OMP_ROOT}/include)
            set(OpenMP_CXX_INCLUDE_DIR ${SE_OMP_ROOT}/include)
        endif()
    endif()
    
    # Try to find OpenMP with better error handling
    if(CMAKE_VERSION VERSION_LESS "3.9")
        # For older CMake versions, use the old FindOpenMP
        find_package(OpenMP QUIET)
    else()
        # For newer CMake versions, use the improved FindOpenMP
        find_package(OpenMP QUIET)
    endif()
    
    if(OpenMP_C_FOUND OR OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found, enabling parallel processing")
        add_compile_definitions(SE_USE_OMP)
        
        # Apply OpenMP flags globally
        if(OpenMP_C_FOUND)
            set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
            if(OpenMP_C_LIBRARIES)
                link_libraries(${OpenMP_C_LIBRARIES})
            endif()
        endif()
        if(OpenMP_CXX_FOUND)
            set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
            if(OpenMP_CXX_LIBRARIES)
                link_libraries(${OpenMP_CXX_LIBRARIES})
            endif()
        endif()
        
        # Add include directories
        if(OpenMP_C_INCLUDE_DIR)
            include_directories(${OpenMP_C_INCLUDE_DIR})
        endif()
        if(OpenMP_CXX_INCLUDE_DIR)
            include_directories(${OpenMP_CXX_INCLUDE_DIR})
        endif()
    else()
        message(WARNING "OpenMP not found. Building without parallel processing support.")
        message(STATUS "You can try:")
        message(STATUS "  - Install OpenMP development packages")
        message(STATUS "  - Use -DSE_CUSTOM_OMP_PATH=ON -DSE_OMP_ROOT=/path/to/openmp")
        message(STATUS "  - Disable OpenMP with -DSE_USE_OMP=OFF")
        set(SE_USE_OMP OFF)
    endif()
endif()

# =============================================================================
# MPI Configuration
# =============================================================================
if(SE_USE_MPI)
    # Set custom MPI search path if provided
    if(SE_CUSTOM_MPI_PATH AND SE_MPI_ROOT)
        message(STATUS "Using custom MPI path: ${SE_MPI_ROOT}")
        set(MPI_ROOT ${SE_MPI_ROOT})
        set(MPI_HOME ${SE_MPI_ROOT})
        
        # Try to find MPI with custom paths first
        find_program(MPI_C_COMPILER 
            NAMES mpicc 
            PATHS ${SE_MPI_ROOT}/bin 
            NO_DEFAULT_PATH
        )
        find_program(MPI_CXX_COMPILER 
            NAMES mpicxx mpic++ 
            PATHS ${SE_MPI_ROOT}/bin 
            NO_DEFAULT_PATH
        )
        find_program(MPIEXEC_EXECUTABLE 
            NAMES mpirun mpiexec 
            PATHS ${SE_MPI_ROOT}/bin 
            NO_DEFAULT_PATH
        )
        
        # If not found in custom path, fall back to system search
        if(NOT MPI_C_COMPILER)
            find_program(MPI_C_COMPILER NAMES mpicc)
        endif()
        if(NOT MPI_CXX_COMPILER)
            find_program(MPI_CXX_COMPILER NAMES mpicxx mpic++)
        endif()
        if(NOT MPIEXEC_EXECUTABLE)
            find_program(MPIEXEC_EXECUTABLE NAMES mpirun mpiexec)
        endif()
    endif()
    
    # Find MPI package (only C and C++ components, no Fortran)
    find_package(MPI COMPONENTS C CXX QUIET)
    
    if(MPI_FOUND OR (MPI_C_FOUND AND MPI_CXX_FOUND))
        message(STATUS "MPI found")
        if(MPI_C_COMPILER)
            message(STATUS "MPI C compiler: ${MPI_C_COMPILER}")
        endif()
        if(MPI_CXX_COMPILER)
            message(STATUS "MPI C++ compiler: ${MPI_CXX_COMPILER}")
        endif()
        if(MPIEXEC_EXECUTABLE)
            message(STATUS "MPI executable: ${MPIEXEC_EXECUTABLE}")
        endif()
        
        # Disable MPI C++ binding warnings (system library issue) when supported
        include(CheckCXXCompilerFlag)
        check_cxx_compiler_flag("-Wno-cast-function-type" COMPILER_SUPPORTS_WNO_CAST_FUNCTION_TYPE)
        # GCC supports -Wno-cast-function-type since ~8.0; older compilers warn about unknown option
        if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU")
            if(CMAKE_CXX_COMPILER_VERSION VERSION_GREATER_EQUAL "8.0")
                if(COMPILER_SUPPORTS_WNO_CAST_FUNCTION_TYPE)
                    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-cast-function-type")
                endif()
            endif()
        elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
            if(COMPILER_SUPPORTS_WNO_CAST_FUNCTION_TYPE)
                set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-cast-function-type")
            endif()
        endif()
    else()
        message(WARNING "MPI not found but SE_USE_MPI is enabled.")
        message(STATUS "You can try:")
        message(STATUS "  - Install MPI development packages:")
        message(STATUS "    Ubuntu/Debian: sudo apt-get install libopenmpi-dev")
        message(STATUS "    CentOS/RHEL: sudo yum install openmpi-devel")
        message(STATUS "  - Use custom MPI path: -DSE_CUSTOM_MPI_PATH=ON -DSE_MPI_ROOT=/path/to/mpi")
        message(STATUS "  - Disable MPI support: -DSE_USE_MPI=OFF")
        message(STATUS "  Continuing without MPI support...")
        set(SE_USE_MPI OFF)
    endif()
endif()

# =============================================================================
# Eigen Configuration
# =============================================================================
if(SE_USE_EIGEN)
    find_package(Eigen3 3.3 QUIET)
    if(Eigen3_FOUND)
        message(STATUS "Eigen3 found")
        add_compile_definitions(SE_USE_EIGEN)
        include_directories(${EIGEN3_INCLUDE_DIRS})
        set(SE_EIGEN_FOUND TRUE)
    else()
        message(WARNING "Eigen3 not found while SE_USE_EIGEN=ON; disabling Eigen support")
        set(SE_USE_EIGEN OFF CACHE BOOL "Enable Eigen-based coherence kernels" FORCE)
        set(SE_EIGEN_FOUND FALSE)
    endif()
endif()

# =============================================================================
# Function: Create Library with Standard Settings
# =============================================================================
function(se_add_library target_name)
    cmake_parse_arguments(SE_LIB "" "VERSION" "SOURCES;HEADERS;DEPENDENCIES" ${ARGN})
    
    if(NOT SE_LIB_VERSION)
        set(SE_LIB_VERSION ${PROJECT_VERSION})
    endif()
    
    # Create static library
    add_library(${target_name} STATIC ${SE_LIB_SOURCES} ${SE_LIB_HEADERS})
    
    # Set standard properties
    set_target_properties(${target_name} PROPERTIES
        OUTPUT_NAME ${target_name}
        ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
        VERSION ${SE_LIB_VERSION}
        SOVERSION ${PROJECT_VERSION_MAJOR}
    )
    
    # Add dependencies if specified
    if(SE_LIB_DEPENDENCIES)
        target_link_libraries(${target_name} ${SE_LIB_DEPENDENCIES})
    endif()
endfunction()

# =============================================================================
# Subdirectory Configuration
# =============================================================================

# Set flag to indicate building as part of main project
set(SE_BUILDING_MAIN_PROJECT TRUE)

# Add subdirectories for each module
add_subdirectory(SEBASIC)
add_subdirectory(SEFILESYSTEM)
add_subdirectory(SERECKIRCH)
add_subdirectory(src)

message(STATUS "=============================================================================")
message(STATUS "Seismic Engine Build Configuration Summary")
message(STATUS "=============================================================================")
message(STATUS "Version: ${PROJECT_VERSION}")
message(STATUS "Build Type: ${CMAKE_BUILD_TYPE}")
message(STATUS "C Compiler: ${CMAKE_C_COMPILER}")
message(STATUS "CXX Compiler: ${CMAKE_CXX_COMPILER}")
message(STATUS "Install Prefix: ${CMAKE_INSTALL_PREFIX}")
message(STATUS "OpenMP Support: ${SE_USE_OMP}")
message(STATUS "MPI Support: ${SE_USE_MPI}")
message(STATUS "Build Examples: ${BUILD_EXAMPLES}")
message(STATUS "=============================================================================")

get_property(SE_ALL_RUNTIME GLOBAL PROPERTY SE_RUNTIME_TARGETS)
if(SE_ALL_RUNTIME)
    list(REMOVE_DUPLICATES SE_ALL_RUNTIME)
endif()

if(SE_BUNDLE_DEPENDENCIES AND SE_ALL_RUNTIME)
    file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/cmake)
    set(SE_DEPLOY_SCRIPT ${CMAKE_CURRENT_BINARY_DIR}/cmake/BundleDependencies.cmake)
    set(SE_DEPLOY_RUNTIME_TARGETS "${SE_ALL_RUNTIME}")
    set(SE_DEPLOY_QT_PLUGIN_DIR "${SE_QT_PLUGIN_DIR}")
    set(SE_DEPLOY_QT_LIBRARY_DIR "${SE_QT_LIBRARY_DIR}")
    set(SE_DEPLOY_ADDITIONAL_DIRS "")
    if(SE_DEPLOY_QT_LIBRARY_DIR)
        list(APPEND SE_DEPLOY_ADDITIONAL_DIRS "${SE_DEPLOY_QT_LIBRARY_DIR}")
    endif()
    set(SE_BUNDLE_FLAT_LAYOUT OFF)
    if(APPLE)
        set(_se_has_app_bundle OFF)
        foreach(_se_runtime ${SE_ALL_RUNTIME})
            get_target_property(_se_is_bundle ${_se_runtime} MACOSX_BUNDLE)
            if(_se_is_bundle)
                set(_se_has_app_bundle ON)
                break()
            endif()
        endforeach()
        if(NOT _se_has_app_bundle)
            set(SE_BUNDLE_FLAT_LAYOUT ON)
        endif()
    endif()
    configure_file(
        ${CMAKE_CURRENT_SOURCE_DIR}/cmake/BundleDependencies.cmake.in
        ${SE_DEPLOY_SCRIPT}
        @ONLY
    )
    install(SCRIPT ${SE_DEPLOY_SCRIPT})
endif()

# Export CMake targets for package consumers
install(EXPORT SeismicEngineTargets
    FILE SeismicEngineTargets.cmake
    NAMESPACE SeismicEngine::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SeismicEngine
)

# Basic CPack configuration for source tarball style releases
set(CPACK_GENERATOR "TGZ")
set(CPACK_PACKAGE_NAME "SeismicEngine")
set(CPACK_PACKAGE_VERSION ${PROJECT_VERSION})
set(CPACK_PACKAGE_CONTACT "changzhimiao")
set(CPACK_PACKAGE_FILE_NAME "${CPACK_PACKAGE_NAME}-${CPACK_PACKAGE_VERSION}-${CMAKE_SYSTEM_NAME}")
include(CPack)