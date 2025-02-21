cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

set(PROJECT_VERSION "0.0.1")

project(
  "golden_cpp"
  LANGUAGES CXX
  VERSION "${PROJECT_VERSION}"
  DESCRIPTION "CPP Golden Project")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})

# ##############################################################################
# OPTIONS
# ##############################################################################

option(PROJECT_TESTS "Build Cpp Tests" ON)
option(PROJECT_SANDBOX "Build Sandbox" OFF)
option(PROJECT_BENCHMARK "Build Benchmark" OFF)
option(PROJECT_SANITIZER_ASAN "Sanitizer ASAN" OFF)
option(PROJECT_SANITIZER_MSAN "Sanitizer MSAN" OFF)
option(PROJECT_SANITIZER_TSAN "Sanitizer TSAN" OFF)

# if(CMAKE_CXX_COMPILER_ID STREQUAL "Clang") elseif(CMAKE_CXX_COMPILER_ID
# STREQUAL "GNU") elseif(CMAKE_CXX_COMPILER_ID STREQUAL "Intel")
# elseif(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC") endif()

set(PROJECT_CXX_FLAGS_DEBUG -fno-omit-frame-pointer -fno-inline)
set(PROJECT_CXX_FLAGS_RELEASE -mfma -mavx) # SIMD? -mftree-vectorize
set(PROJECT_CXX_FLAGS_RELWITHDEBUG ${PROJECT_CXX_FLAGS_RELEASE} -g)
set(PROJECT_CXX_FLAGS_MINSIZERELEASE ${PROJECT_CXX_FLAGS_RELEASE} -g)

add_compile_options(
  ${CMAKE_CXX_FLAGS}
  "$<$<CONFIG:Debug>:${PROJECT_CXX_FLAGS_DEBUG}>"
  "$<$<CONFIG:Release>:${PROJECT_CXX_FLAGS_RELEASE}>"
  "$<$<CONFIG:RelWithDebug>:${PROJECT_CXX_FLAGS_RELWITHDEBUG}>"
  "$<$<CONFIG:MinSizeRelease>:${PROJECT_CXX_FLAGS_MINSIZERELEASE}>")

if(${PROJECT_ASAN})
  set(CMAKE_CXX_FLAGS "-O -g -fsanitize=address")
endif()

if(${PROJECT_TSAN})
  set(CMAKE_CXX_FLAGS "-O -g -fsanitize=thread -fopenmp")
endif()

if(${PROJECT_MSAN})
  set(CMAKE_CXX_FLAGS
      "-fsanitize=memory -fPIE -pie -fno-omit-frame-pointer -g -O2")
endif()

# ##############################################################################
# FIND PACKAGES
# ##############################################################################
# * Boost
find_package(
  Boost REQUIRED
  COMPONENTS filesystem system program_options unit_test_framework
  REQUIRED)
include_directories(SYSTEM ${Boost_INCLUDE_DIR})

# ##############################################################################
# INCLUDE
# ##############################################################################
# configure_file("${PROJECT_SOURCE_DIR}/src/linalg/utils/path.hpp.in"
# "${PROJECT_BINARY_DIR}/path.hpp")
include_directories("${PROJECT_SOURCE_DIR}/src" ${PROJECT_BINARY_DIR})

# ##############################################################################
# PROJECT
# ##############################################################################
add_subdirectory(src)

# ##############################################################################
# TESTS and miscelaneous
# ##############################################################################
if(${PROJECT_TESTS})
  include(CTest)
  enable_testing()
  add_subdirectory(tests/unit)
endif()
if(${PROJECT_SANDBOX})
  add_subdirectory(tests/sandbox)
endif()
if(${PROJECT_BENCHMARK})
  add_subdirectory(tests/benchmark)
endif()
if(${PROJECT_VALIDATION})
  add_subdirectory(tests/validation)
endif()
