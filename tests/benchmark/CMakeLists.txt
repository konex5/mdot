
set(BENCHMARKS_NAMES bloc_norm matrix_size)

foreach(BENCHMARK_NAME ${BENCHMARKS_NAMES})
  add_executable(bench_${BENCHMARK_NAME} bench_${BENCHMARK_NAME}.cpp)
  target_link_libraries(
    bench_${BENCHMARK_NAME} ${Boost_FILESYSTEM_LIBRARY} ${Boost_SYSTEM_LIBRARY}
    ${Boost_BENCHMARK_FRAMEWORK_LIBRARY} Boost::system Boost::program_options tbb)
endforeach()
