set(libraries mat)

foreach(library ${libraries})

  if(${PROJECT_STATIC})
    add_library(${library} STATIC core/${library}.cpp)
  else()
    add_library(${library} SHARED core/${library}.cpp)
  endif()
  target_link_libraries(${library} ${TBB_TBB_LIBRARY_RELEASE} Boost::system
  				   Boost::filesystem)

endforeach()

set(clicommands climdot cliprintoperator create_maximal_entangled_state)

foreach(clicommand ${clicommands})

  add_executable(cli_${clicommand} cli/${clicommand}.cpp)
  target_link_libraries(cli_${clicommand} ${BOOST_LIBRARIES} Boost::filesystem Boost::program_options)

  install(
    TARGETS cli_${clicommand}
    CONFIGURATIONS Release
    DESTINATION bin)

endforeach()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/mdot/include/ DESTINATION include/mdot/include)
