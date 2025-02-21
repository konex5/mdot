cmake_minimum_required(VERSION 3.1.0 FATAL_ERROR)

set(PROJECT_VERSION "0.0.1")

project(
  "golden-cpp"
  LANGUAGES CXX C
  VERSION "${PROJECT_VERSION}"
  DESCRIPTION "CPP Golden Project")

set(CMAKE_C_COMPILER $ENV{CC})
set(CMAKE_CXX_COMPILER $ENV{CXX})

list(APPEND CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/CMake)

# ##############################################################################
# OPTIONS
# ##############################################################################

option(PROJECT_BENCHMARK "Compile the performance in ctest" OFF)
option(PROJECT_INTEGRATION "Compile the integration tests in ctest" OFF)
option(PROJECT_UNIT_TESTS "Compile the unit tests in ctest" ON)
option(PROJECT_SANDBOX "Compile the sandbox folder" OFF)
option(PROJECT_STATIC "Build Static Libraries" OFF)
option(PROJECT_MEMCHECK "Enable memcheck test with valgrind" OFF)

# ##############################################################################
# FLAGS
# ##############################################################################

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

add_compile_options(-Wall -Wextra)
# #
# https://stackoverflow.com/questions/45955272/modern-way-to-set-compiler-flags-in-cross-platform-cmake-project
if(CMAKE_CXX_COMPILER_ID STREQUAL "GNU" AND APPLE)
  add_compile_options("-Wa" "-q")
else()
  add_compile_options(
    $<$<CONFIG:DEBUG>:-Werror> $<$<CONFIG:DEBUG>:-pedantic>
    $<$<CONFIG:DEBUG>:-Wunused-macros> $<$<CONFIG:DEBUG>:-Wshadow>)
  add_compile_options($<$<CONFIG:RELEASE>:-mfma> $<$<CONFIG:RELEASE>:-mavx>
  )# ?=> -mftree-vectorize
  if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
    add_compile_options(
      $<$<CONFIG:ASAN>:-fsanitize=address> $<$<CONFIG:ASAN>:-g>
      $<$<CONFIG:ASAN>:-O1> $<$<CONFIG:ASAN>:-fno-omit-frame-pointer>)
    add_compile_options(
      $<$<CONFIG:MSAN>:-g>
      $<$<CONFIG:MSAN>:-O2>
      $<$<CONFIG:MSAN>:-fsanitize=memory>
      $<$<CONFIG:MSAN>:-fsanitize-memory-track-origins=2>
      $<$<CONFIG:MSAN>:-fno-omit-frame-pointer>
      $<$<CONFIG:MSAN>:-fno-optimize-sibling-calls>)
    add_compile_options($<$<CONFIG:TSAN>:-g> $<$<CONFIG:TSAN>:-O1>
                        $<$<CONFIG:TSAN>:-fsanitize=thread>)
    add_compile_options(
      $<$<CONFIG:UBSAN>:-g> $<$<CONFIG:UBSAN>:-fno-omit-frame-pointer>
      $<$<CONFIG:UBSAN>:-fsanitize=undefined>)
    add_compile_options(
      $<$<CONFIG:CODECOV>:-g> $<$<CONFIG:CODECOV>:-fprofile-instr-generate>
      $<$<CONFIG:CODECOV>:-fcoverage-mapping> $<$<CONFIG:CODECOV>:-fno-inline>)
  else()
    add_compile_options($<$<CONFIG:ASAN>:-g> $<$<CONFIG:ASAN>:-Wshadow>)
    add_compile_options($<$<CONFIG:TSAN>:-g> $<$<CONFIG:TSAN>:-Wshadow>)
    add_compile_options($<$<CONFIG:CODECOV>:-g> $<$<CONFIG:CODECOV>:-Wshadow>)
    add_compile_options($<$<CONFIG:CODECOV>:-g> $<$<CONFIG:CODECOV>:-Wshadow>)
  endif()
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
add_subdirectory(src/golden)
add_subdirectory(src/silver)

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
