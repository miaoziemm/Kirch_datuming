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
add_executable(layer_kirchhoff_test layer_kirchhoff_test.cpp)
add_executable(layer_wavefield_fill layer_wavefield_fill.cpp)
add_executable(layered_frequency_wavefield layered_frequency_wavefield.cpp)
add_executable(layered_model_gen layered_model_gen.cpp)
add_executable(shotgen shotgen.cpp)
add_executable(greensol greensol.cpp)
add_executable(greensol_auto greensol_auto.cpp)
add_executable(layercom layercom.cpp)
add_executable(window window.cpp)
add_executable(layervel layervel.cpp)
add_executable(kirchdat2d_auto_cs kirchdat2d_auto_cs.cpp)
add_executable(rsf_synth_gen rsf_synth_gen.cpp)
add_executable(rsf_l2cmp rsf_l2cmp.cpp)
add_executable(kirchmig2d_cig kirchmig2d_cig.cpp)
add_executable(sewave2d sewave2d_app.cpp ../SEWAVE/sewave2d.cpp)
add_executable(sewave2d_synth sewave2d_synth.cpp ../SEWAVE/sewave2d.cpp)

# 查找系统库（SQLite3 依赖）
find_library(DL_LIBRARY dl)
find_package(Threads REQUIRED)
find_package(OpenMP QUIET)

# 链接库
# target_link_libraries(test ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchmig2d ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(eikods2d ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchdat2d ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchdat2d_auto ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchdat2d_freq ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} fftw3f ${DL_LIBRARY} Threads::Threads)
target_link_libraries(layer_kirchhoff_test ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} fftw3f ${DL_LIBRARY} Threads::Threads)
target_link_libraries(layer_wavefield_fill ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(layered_frequency_wavefield ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} fftw3f ${DL_LIBRARY} Threads::Threads)
target_link_libraries(layered_model_gen ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(shotgen ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(greensol ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(greensol_auto ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(layercom ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(window ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(layervel ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchdat2d_auto_cs ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(rsf_synth_gen ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(rsf_l2cmp ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(kirchmig2d_cig ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} ${SERECKIRCH_LIB} ${DL_LIBRARY} Threads::Threads)
target_link_libraries(sewave2d ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} fftw3f ${DL_LIBRARY} Threads::Threads)
target_link_libraries(sewave2d_synth ${SEFILESYSTEM_LIB} ${SEBASIC_LIB} fftw3f ${DL_LIBRARY} Threads::Threads)
if(OpenMP_CXX_FOUND)
    target_link_libraries(sewave2d OpenMP::OpenMP_CXX)
    target_link_libraries(sewave2d_synth OpenMP::OpenMP_CXX)
endif()
target_include_directories(kirchdat2d_freq PRIVATE ${CMAKE_SOURCE_DIR}/fftw/api)
target_include_directories(layer_kirchhoff_test PRIVATE ${CMAKE_SOURCE_DIR}/fftw/api)
target_include_directories(layered_frequency_wavefield PRIVATE ${CMAKE_SOURCE_DIR}/fftw/api)
target_include_directories(sewave2d PRIVATE ${CMAKE_SOURCE_DIR}/fftw/api)
target_include_directories(sewave2d_synth PRIVATE ${CMAKE_SOURCE_DIR}/fftw/api)

# 安装所有可执行文件到 bin 目录
install(TARGETS
    # test
    kirchmig2d
    eikods2d
    kirchdat2d
    kirchdat2d_auto
    kirchdat2d_freq
    layer_kirchhoff_test
    layer_wavefield_fill
    layered_frequency_wavefield
    layered_model_gen
    shotgen
    greensol
    greensol_auto
    layercom
    window
    layervel
    kirchdat2d_auto_cs
    rsf_synth_gen
    rsf_l2cmp
    kirchmig2d_cig
    sewave2d
    sewave2d_synth
)

# 确保依赖的库先编译完成
# add_dependencies(test sebasic sefs reckirch)
add_dependencies(kirchmig2d sebasic sefs reckirch)
add_dependencies(eikods2d sebasic sefs reckirch)
add_dependencies(kirchdat2d sebasic sefs reckirch)
add_dependencies(kirchdat2d_auto sebasic sefs reckirch)
add_dependencies(kirchdat2d_freq sebasic sefs reckirch fftw3f)
add_dependencies(layer_kirchhoff_test sebasic sefs reckirch fftw3f)
add_dependencies(layer_wavefield_fill sebasic sefs)
add_dependencies(layered_frequency_wavefield sebasic sefs reckirch fftw3f)
add_dependencies(layered_model_gen sebasic sefs)
add_dependencies(shotgen sebasic sefs reckirch)
add_dependencies(greensol sebasic sefs reckirch)
add_dependencies(greensol_auto sebasic sefs reckirch)
add_dependencies(layercom sebasic sefs reckirch)
add_dependencies(window sebasic sefs reckirch)
add_dependencies(layervel sebasic sefs reckirch)
add_dependencies(kirchdat2d_auto_cs sebasic sefs reckirch)
add_dependencies(rsf_synth_gen sebasic sefs reckirch)
add_dependencies(rsf_l2cmp sebasic sefs reckirch)
add_dependencies(kirchmig2d_cig sebasic sefs reckirch)
add_dependencies(sewave2d sebasic sefs fftw3f)
add_dependencies(sewave2d_synth sebasic sefs fftw3f)
