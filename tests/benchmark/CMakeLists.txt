
set(BENCHMARKS_NAMES matrix_size)

foreach(BENCHMARK_NAME ${BENCHMARKS_NAMES})
  add_executable(bench_${BENCHMARK_NAME} bench_${BENCHMARK_NAME}.cpp)
  target_compile_definitions(
    bench_${BENCHMARK_NAME} PUBLIC BOOST_TEST_MODULE=${BENCHMARK_NAME}
                                  BOOST_TEST_DYN_LINK=TRUE)
  target_link_libraries(
    bench_${BENCHMARK_NAME} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY}
    ${Boost_BENCHMARK_FRAMEWORK_LIBRARY} Boost::system Boost::program_options tbb)
  add_test(unit_${BENCHMARK_NAME} bench_${BENCHMARK_NAME})
endforeach()
