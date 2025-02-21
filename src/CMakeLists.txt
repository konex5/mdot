add_executable(cli_golden_cpp cli_golden_cpp.cpp)
# add_library(golden SHARED src/golden.cpp src/silver/silver.hpp)
# target_link_libraries(golden ${BOOST_LIBRARIES})

install(
  TARGETS cli_golden_cpp
  CONFIGURATIONS Debug
  DESTINATION debug/bin)

install(
  TARGETS cli_golden_cpp
  CONFIGURATIONS Release
  DESTINATION bin)
