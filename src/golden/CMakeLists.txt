set(libraries gold)

foreach(library ${libraries})

  if(${PROJECT_STATIC})
    add_library(${library} STATIC core/${library}.cpp)
  else()
    add_library(${library} SHARED core/${library}.cpp)
  endif()
  target_link_libraries(${library}) # ${BOOST_LIBRARIES} Boost::system
                                    # Boost::filesystem

endforeach()

add_executable(cli_golden_cpp cli/golden.cpp)
# target_link_libraries(golden ${BOOST_LIBRARIES})

install(
  TARGETS cli_golden_cpp
  CONFIGURATIONS Release
  DESTINATION bin)
