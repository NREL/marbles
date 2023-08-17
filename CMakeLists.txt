cmake_minimum_required (VERSION 3.18 FATAL_ERROR)
project(Marbles CXX)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON CACHE INTERNAL "")
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/CMake")
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)
include(marbles-utils)

option(MARBLES_ENABLE_ALL_WARNINGS "Show most warnings for most compilers" ON)
option(MARBLES_ENABLE_CLANG_TIDY "Compile with clang-tidy static analysis" OFF)
option(MARBLES_ENABLE_CPPCHECK "Enable cppcheck static analysis target" OFF)
option(MARBLES_ENABLE_FCOMPARE "Enable building fcompare when not testing" OFF)
option(MARBLES_TEST_WITH_FCOMPARE "Check test plots against gold files" OFF)
option(MARBLES_SAVE_GOLDS "Provide a directory in which to save golds during testing" OFF)
option(MARBLES_ENABLE_FPE_TRAP_FOR_TESTS "Enable FPE trapping in tests" ON)
option(MARBLES_ENABLE_MPI "Enable MPI" OFF)
option(MARBLES_ENABLE_OPENMP "Enable OpenMP" OFF)
option(MARBLES_ENABLE_CUDA "Enable CUDA" OFF)
option(MARBLES_ENABLE_HIP "Enable HIP" OFF)
option(MARBLES_ENABLE_SYCL "Enable Intel OneAPI SyCL" OFF)
option(MARBLES_ENABLE_TINY_PROFILE "Enable AMReX TinyProfile support" ON)
set(MARBLES_PRECISION "DOUBLE" CACHE STRING "Floating point precision SINGLE or DOUBLE")

# Third party libraries
option(MARBLES_ENABLE_HDF5 "Enable HDF5 library" OFF)
option(MARBLES_ENABLE_HDF5_ZFP "Enable ZFP compression in HDF5 library" OFF)

#Options for C++
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

if(MARBLES_ENABLE_CUDA)
  enable_language(CUDA)
  if(CMAKE_CUDA_COMPILER_VERSION VERSION_LESS "10.0")
    message(FATAL_ERROR "Your nvcc version is ${CMAKE_CUDA_COMPILER_VERSION} which is unsupported."
      "Please use CUDA toolkit version 10.0 or newer.")
  endif()
  if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    set(CMAKE_CUDA_ARCHITECTURES 70)
  endif()
endif()

if(MARBLES_ENABLE_HIP)
  find_package(HIP REQUIRED)
  if(NOT DEFINED AMD_ARCH)
    set(AMD_ARCH "gfx90a")
  endif()
endif()

if(MARBLES_TEST_WITH_FCOMPARE)
  set(MARBLES_ENABLE_FCOMPARE ON)
endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "AppleClang")
  set(MARBLES_ENABLE_FPE_TRAP_FOR_TESTS OFF)
  message(WARNING "Disabling FPE trapping for tests when using AppleClang.")
endif()

init_amrex()

if(MARBLES_ENABLE_MPI)
  find_package(MPI REQUIRED)
endif()

message(STATUS "Marbles Information:")
message(STATUS "CMAKE_SYSTEM_NAME = ${CMAKE_SYSTEM_NAME}")
message(STATUS "CMAKE_CXX_COMPILER_ID = ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "CMAKE_CXX_COMPILER_VERSION = ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}")

include(set_rpath)

set(marbles_lib_name "marbles_obj")
set(marbles_exe_name "marbles")

add_library(${marbles_lib_name} OBJECT)
add_executable(${marbles_exe_name})

init_code_checks()
if(CLANG_TIDY_EXE)
  set_target_properties(${marbles_lib_name} ${marbles_api_lib} ${marbles_exe_name}
                        PROPERTIES CXX_CLANG_TIDY "${CLANG_TIDY_EXE};--config-file=${CMAKE_SOURCE_DIR}/.clang-tidy")
endif()

include(set_compile_flags)

if(MARBLES_ENABLE_HDF5)
  set(CMAKE_PREFIX_PATH ${HDF5_DIR} ${CMAKE_PREFIX_PATH})
  target_compile_definitions(${marbles_lib_name} PUBLIC MARBLES_USE_HDF5)
endif()

if(MARBLES_ENABLE_HDF5_ZFP)
  if(NOT MARBLES_ENABLE_HDF5)
    message(FATAL_ERROR "HDF5 must be enabled to use ZFP compression in HDF5 plot files")
  endif()
  set(CMAKE_PREFIX_PATH ${H5Z_ZFP_DIR} ${CMAKE_PREFIX_PATH})
  target_compile_definitions(${marbles_lib_name} PUBLIC MARBLES_USE_HDF5_ZFP)
endif()

generate_version_info()
add_subdirectory(Source)

if(MARBLES_ENABLE_CUDA)
  set(mtargets "${marbles_lib_name};${marbles_exe_name}")
  foreach(tgt IN LISTS mtargets)
    set_cuda_build_properties(${tgt})
  endforeach()
endif()

enable_testing()
include(CTest)
add_subdirectory(Tests)

#Define what we want to be installed during a make install
install(
  TARGETS ${marbles_exe_name}
  ${marbles_api_lib} ${marbles_lib_name} buildInfo${marbles_lib_name}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME DESTINATION bin
  ARCHIVE DESTINATION lib
  LIBRARY DESTINATION lib)

install(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION lib/cmake/${PROJECT_NAME})

configure_package_config_file(
  CMake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION lib/cmake/${PROJECT_NAME})
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  DESTINATION lib/cmake/${PROJECT_NAME})
install(FILES
  DESTINATION ${CMAKE_INSTALL_DATADIR}/cmake/Modules)
