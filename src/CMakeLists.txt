# CMakeLists.txt for src directory
cmake_minimum_required(VERSION 3.10)

if(NOT DEFINED SE_INSTALL_EXPORT_NAME)
    set(SE_INSTALL_EXPORT_NAME SeismicEngineTargets)
endif()

# 包含父目录的头文件路径
include_directories(${CMAKE_SOURCE_DIR})

# 查找 sebasic 库
find_library(SEBASIC_LIB
    NAMES sebasic
    PATHS ${CMAKE_BINARY_DIR}/lib
    NO_DEFAULT_PATH
)

# 如果没有找到库，添加依赖到 sebasic 目标
if(NOT SEBASIC_LIB)
    set(SEBASIC_LIB sebasic)
endif()

# 查找 SEFILESYSTEM 库
find_library(SEFILESYSTEM_LIB
    NAMES sefs
    PATHS ${CMAKE_BINARY_DIR}/lib
    NO_DEFAULT_PATH
)

# 如果没有找到库，添加依赖到 sefilesystem 目标
if(NOT SERECKIRCH_LIB)
    set(SERECKIRCH_LIB reckirch)
endif()

# 查找 SEFILESYSTEM 库
find_library(SERECKIRCH_LIB
    NAMES reckirch
    PATHS ${CMAKE_BINARY_DIR}/lib
    NO_DEFAULT_PATH
)

# 如果没有找到库，添加依赖到 sefilesystem 目标
if(NOT SERECKIRCH_LIB)
    set(SERECKIRCH_LIB reckirch)
endif()


# 编译可执行文件
add_executable(test test.cpp)
add_executable(kirchmig2d kirchmig2d.cpp)
add_executable(eikods2d eikods2d.cpp)
add_executable(kirchdat2d kirchdat2d.cpp)
add_executable(shotgen shotgen.cpp)

# 查找系统库（SQLite3 依赖）
find_library(DL_LIBRARY dl)
find_package(Threads REQUIRED)

# 链接库
target_link_libraries(test ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchmig2d ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(eikods2d ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchdat2d ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(shotgen ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)

# 安装所有可执行文件到 bin 目录
install(TARGETS
    test
    kirchmig2d
    eikods2d
    kirchdat2d
    shotgen
)

# 确保依赖的库先编译完成
add_dependencies(test sebasic sefs reckirch)
add_dependencies(kirchmig2d sebasic sefs reckirch)
add_dependencies(eikods2d sebasic sefs reckirch)
add_dependencies(kirchdat2d sebasic sefs reckirch)
add_dependencies(shotgen sebasic sefs reckirch)
