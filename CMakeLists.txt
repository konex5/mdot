cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

# set(PROJECT_VERSION "0.0.1")

project(
  "mdot"
  LANGUAGES C CXX
  VERSION "${PROJECT_VERSION}"
  DESCRIPTION "MDOT - Matrix Dot")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_VERBOSE_MAKEFILE ON)
set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})


list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)

# ##############################################################################
# OPTIONS
# ##############################################################################

option(PROJECT_BENCHMARK "Compile the performance in ctest" ON)
option(PROJECT_INTEGRATION "Compile the integration tests in ctest" OFF)
option(PROJECT_UNIT_TESTS "Compile the unit tests in ctest" ON)
option(PROJECT_SANDBOX "Compile the sandbox folder" OFF)
option(PROJECT_STATIC "Build Static Libraries" OFF)
option(PROJECT_MEMCHECK "Enable memcheck test with valgrind" OFF)

# ##############################################################################
# FLAGS
# ##############################################################################

set(PROJECT_CXX_FLAGS -Wall -Wextra)
set(PROJECT_CXX_FLAGS_DEBUG -g -fno-inline -fno-omit-frame-pointer) # -pedantic
# set(PROJECT_CXX_FLAGS_RELWITHDEBUG -g -O2) set(PROJECT_CXX_FLAGS_RELEASE -O3)
# set(PROJECT_CXX_FLAGS_MINSIZERELEASE -O2) -fno-inline)
# set(PROJECT_CXX_FLAGS_RELEASE -mfma -mavx) # -q openmp-simd -mftree-vectorize

add_compile_options(
  ${PROJECT_CXX_FLAGS}
  "$<$<CONFIG:Debug>:${PROJECT_CXX_FLAGS_DEBUG}>"
  "$<$<CONFIG:Release>:${PROJECT_CXX_FLAGS_RELEASE}>"
  "$<$<CONFIG:RelWithDebug>:${PROJECT_CXX_FLAGS_RELWITHDEBUG}>"
  "$<$<CONFIG:MinSizeRelease>:${PROJECT_CXX_FLAGS_MINSIZERELEASE}>")

# ##############################################################################
# FIND PACKAGES
# ##############################################################################
# * Boost
find_package(
  Boost REQUIRED
  COMPONENTS filesystem system program_options unit_test_framework
  REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

mark_as_advanced(
  Boost_SYSTEM_LIBRARY_RELEASE Boost_FILESYSTEM_LIBRARY_RELEASE
  Boost_PROGRAM_OPTIONS_LIBRARY_RELEASE
  Boost_UNIT_TEST_FRAMEWORK_LIBRARY_RELEASE)

# * TBB
find_package(TBB REQUIRED tbb tbbmalloc_proxy)

# ##############################################################################
# INCLUDE
# ##############################################################################
configure_file("${PROJECT_SOURCE_DIR}/src/utils/static.hpp.in"
               "${PROJECT_BINARY_DIR}/utils/static.hpp")
include_directories("${PROJECT_SOURCE_DIR}/src" ${PROJECT_BINARY_DIR})

# ##############################################################################
# PROJECT
# ##############################################################################
add_subdirectory(src/mdot)

# ##############################################################################
# TESTS and miscelaneous
# ##############################################################################
include(CTest)
# enable_testing()
if(${PROJECT_UNIT_TESTS})
  add_subdirectory(tests/unit)
endif()
if(${PROJECT_SANDBOX})
  add_subdirectory(tests/sandbox)
endif()
if(${PROJECT_BENCHMARK})
  add_subdirectory(tests/benchmark)
endif()
if(${PROJECT_INTEGRATION})
  add_subdirectory(tests/integration)
endif()
