# =============================================================================
# RECKIRCH Library
# =============================================================================
cmake_minimum_required(VERSION 3.16)

project(reckirch LANGUAGES CXX)

include(GNUInstallDirs)

if(NOT DEFINED SE_INSTALL_EXPORT_NAME)
    set(SE_INSTALL_EXPORT_NAME SeismicEngineTargets)
endif()

# =============================================================================
# Header files
# =============================================================================
set(RECKIRCH_HEADERS
    include/decart.h
    include/eikods.h
    include/mutter.h
    include/neighbors.h
    include/pqueue.h
    include/se_reckirch.h
    include/tinterpsr.h
    include/tpow.h
    include/kirmig.h
    include/kirdat.h
)

# =============================================================================
# Source files
# =============================================================================
set(RECKIRCH_SOURCES
    src/decart.cpp
    src/eikods.cpp
    src/mutter.cpp
    src/neighbors.cpp
    src/pqueue.cpp
    src/tinterpsr.cpp
    src/tpow.cpp
    src/kirmig.cpp
    src/kirdat.cpp
)

# =============================================================================
# Create static library
# =============================================================================
add_library(reckirch STATIC
    ${RECKIRCH_SOURCES}
    ${RECKIRCH_HEADERS}
)

# =============================================================================
# C++ standard
# =============================================================================
target_compile_features(reckirch PUBLIC cxx_std_11)

# =============================================================================
# Include directories
# =============================================================================
target_include_directories(reckirch PUBLIC
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
)

# =============================================================================
# Library properties
# =============================================================================
set_target_properties(reckirch PROPERTIES
    OUTPUT_NAME "reckirch"
    ARCHIVE_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib
    POSITION_INDEPENDENT_CODE ON
)

# =============================================================================
# Install rules
# =============================================================================
install(TARGETS reckirch
    EXPORT ${SE_INSTALL_EXPORT_NAME}
    ARCHIVE DESTINATION ${CMAKE_INSTALL_LIBDIR}
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/include/
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

message(STATUS "RECKIRCH library configured")