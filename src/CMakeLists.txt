# CMakeLists.txt for src directory
cmake_minimum_required(VERSION 3.10)

if(NOT DEFINED SE_INSTALL_EXPORT_NAME)
    set(SE_INSTALL_EXPORT_NAME SeismicEngineTargets)
endif()

# 包含父目录的头文件路径
include_directories(${CMAKE_SOURCE_DIR})

# 查找 sebasic 库
if(TARGET sebasic)
    set(SEBASIC_LIB sebasic)
else()
    find_library(SEBASIC_LIB
        NAMES sebasic
        PATHS ${CMAKE_BINARY_DIR}/lib
        NO_DEFAULT_PATH
    )
    if(NOT SEBASIC_LIB)
        set(SEBASIC_LIB sebasic)
    endif()
endif()

# 查找 SEFILESYSTEM 库
if(TARGET sefs)
    set(SEFILESYSTEM_LIB sefs)
else()
    find_library(SEFILESYSTEM_LIB
        NAMES sefs
        PATHS ${CMAKE_BINARY_DIR}/lib
        NO_DEFAULT_PATH
    )
    if(NOT SEFILESYSTEM_LIB)
        set(SEFILESYSTEM_LIB sefs)
    endif()
endif()

# 查找 SERECKIRCH 库
if(TARGET reckirch)
    set(SERECKIRCH_LIB reckirch)
else()
    find_library(SERECKIRCH_LIB
        NAMES reckirch
        PATHS ${CMAKE_BINARY_DIR}/lib
        NO_DEFAULT_PATH
    )
    if(NOT SERECKIRCH_LIB)
        set(SERECKIRCH_LIB reckirch)
    endif()
endif()


# 编译可执行文件
# add_executable(test test.cpp)
add_executable(kirchmig2d kirchmig2d.cpp)
add_executable(eikods2d eikods2d.cpp)
add_executable(kirchdat2d kirchdat2d.cpp)
add_executable(kirchdat2d_auto kirchdat2d_auto.cpp)
add_executable(kirchdat2d_freq kirchdat2d_frequency_cs.cpp)
add_executable(shotgen shotgen.cpp)
add_executable(greensol greensol.cpp)
add_executable(greensol_auto greensol_auto.cpp)
add_executable(layercom layercom.cpp)
add_executable(window window.cpp)
add_executable(layervel layervel.cpp)
add_executable(kirchdat2d_auto_cs kirchdat2d_auto_cs.cpp)
add_executable(rsf_synth_gen rsf_synth_gen.cpp)
add_executable(rsf_l2cmp rsf_l2cmp.cpp)

# 查找系统库（SQLite3 依赖）
find_library(DL_LIBRARY dl)
find_package(Threads REQUIRED)

# 链接库
# target_link_libraries(test ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchmig2d ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(eikods2d ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchdat2d ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchdat2d_auto ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchdat2d_freq ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} fftw3f ${DL_LIBRARY} Threads::Threads)
target_link_libraries(shotgen ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(greensol ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(greensol_auto ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(layercom ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(window ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(layervel ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchdat2d_auto_cs ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(rsf_synth_gen ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(rsf_l2cmp ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)

target_include_directories(kirchdat2d_freq PRIVATE ${CMAKE_SOURCE_DIR}/fftw/api)

# 安装所有可执行文件到 bin 目录
install(TARGETS
    # test
    kirchmig2d
    eikods2d
    kirchdat2d
    kirchdat2d_auto
    kirchdat2d_freq
    shotgen
    greensol
    greensol_auto
    layercom
    window
    layervel
    kirchdat2d_auto_cs
    rsf_synth_gen
    rsf_l2cmp
)

# 确保依赖的库先编译完成
# add_dependencies(test sebasic sefs reckirch)
add_dependencies(kirchmig2d sebasic sefs reckirch)
add_dependencies(eikods2d sebasic sefs reckirch)
add_dependencies(kirchdat2d sebasic sefs reckirch)
add_dependencies(kirchdat2d_auto sebasic sefs reckirch)
add_dependencies(kirchdat2d_freq sebasic sefs reckirch fftw3f)
add_dependencies(shotgen sebasic sefs reckirch)
add_dependencies(greensol sebasic sefs reckirch)
add_dependencies(greensol_auto sebasic sefs reckirch)
add_dependencies(layercom sebasic sefs reckirch)
add_dependencies(window sebasic sefs reckirch)
add_dependencies(layervel sebasic sefs reckirch)
add_dependencies(kirchdat2d_auto_cs sebasic sefs reckirch)
add_dependencies(rsf_synth_gen sebasic sefs reckirch)
add_dependencies(rsf_l2cmp sebasic sefs reckirch)
