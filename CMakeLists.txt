cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

set(PROJECT_VERSION "0.0.1")

project(
  "golden-cpp"
  LANGUAGES CXX C
  VERSION "${PROJECT_VERSION}"
  DESCRIPTION "CPP Golden Project")

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)

# ##############################################################################
# OPTIONS
# ##############################################################################

option(PROJECT_BENCHMARK "Compile the performance in ctest" OFF)
option(PROJECT_INTEGRATION "Compile the integration tests in ctest" ON)
option(PROJECT_UNIT_TESTS "Compile the unit tests in ctest" ON)
option(PROJECT_SANDBOX "Compile the sandbox folder" OFF)
option(PROJECT_STATIC "Build Static Libraries" OFF)

set(PROJECT_CXX_FLAGS_DEBUG -fno-omit-frame-pointer -fno-inline)
set(PROJECT_CXX_FLAGS_RELEASE -mfma -mavx) # SIMD? -mftree-vectorize
set(PROJECT_CXX_FLAGS_RELWITHDEBUG ${PROJECT_CXX_FLAGS_RELEASE} -g)
set(PROJECT_CXX_FLAGS_MINSIZERELEASE ${PROJECT_CXX_FLAGS_RELEASE})

add_compile_options(
  ${CMAKE_CXX_FLAGS}
  "$<$<CONFIG:Debug>:${PROJECT_CXX_FLAGS_DEBUG}>"
  "$<$<CONFIG:Release>:${PROJECT_CXX_FLAGS_RELEASE}>"
  "$<$<CONFIG:RelWithDebug>:${PROJECT_CXX_FLAGS_RELWITHDEBUG}>"
  "$<$<CONFIG:MinSizeRelease>:${PROJECT_CXX_FLAGS_MINSIZERELEASE}>")

# option(PROJECT_SANITIZER_ASAN "Sanitizer ASAN" OFF)
# option(PROJECT_SANITIZER_MSAN "Sanitizer MSAN" OFF)
# option(PROJECT_SANITIZER_TSAN "Sanitizer TSAN" OFF)

# if(${PROJECT_ASAN}) set(CMAKE_CXX_FLAGS "-O -g -fsanitize=address") endif()

# if(${PROJECT_TSAN}) set(CMAKE_CXX_FLAGS "-O -g -fsanitize=thread -fopenmp")
# endif()

# if(${PROJECT_MSAN}) set(CMAKE_CXX_FLAGS "-fsanitize=memory -fPIE -pie
# -fno-omit-frame-pointer -g -O2") endif()

if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND APPLE)
  set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wa,-q")
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
add_subdirectory(src/golden)
add_subdirectory(src/silver)

# ##############################################################################
# TESTS and miscelaneous
# ##############################################################################
include(CTest)
enable_testing()
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
