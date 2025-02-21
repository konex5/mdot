set(libraries mat)

foreach(library ${libraries})

  if(${PROJECT_STATIC})
    add_library(${library} STATIC core/${library}.cpp)
  else()
    add_library(${library} SHARED core/${library}.cpp)
  endif()
  target_link_libraries(${library}) # ${BOOST_LIBRARIES} Boost::system
                                    # Boost::filesystem

endforeach()

set(clicommands climdot cliprintoperator)

foreach(clicommand ${clicommands})

  add_executable(cli_${clicommand} cli/${clicommand}.cpp)
  target_link_libraries(cli_${clicommand} ${BOOST_LIBRARIES}
                        Boost::program_options)

  install(
    TARGETS cli_${clicommand}
    CONFIGURATIONS Release
    DESTINATION bin)

endforeach()

install(DIRECTORY ${PROJECT_SOURCE_DIR}/src/mdot/include/ DESTINATION include/mdot/include)
