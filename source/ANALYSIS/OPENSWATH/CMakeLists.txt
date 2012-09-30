cmake_minimum_required (VERSION 2.6)

project(OpenSwath)

set( BinOut ${CMAKE_SOURCE_DIR}/${CMAKE_PROJECT_NAME}out/${} )
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${BinOut} )
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${BinOut} )
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${BinOut} )

find_package(Doxygen)
#if(DOXYGEN_FOUND)
#  configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile ${CMAKE_CURRENT_BINARY_DIR}/Doxyfile @ONLY)
#  add_custom_target(doc
#    ${DOXYGEN_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/Doxyfile
#    WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
#    COMMENT "Generating API documentation with Doxygen" VERBATIM)
#endif(DOXYGEN_FOUND)

find_package( Boost )
if(Boost_FOUND)
  SET(inc ${inc} ${Boost_INCLUDE_DIRS})
  SET(lib ${lib} ${Boost_LIBRARIES})
else()
  message(STATUS ">>>>> Couldn't find Boost <<<<< ")
endif()

SET(inc ${inc} ${CMAKE_SOURCE_DIR}/include)


message(STATUS "+++++++++++++include             : ${inc}")

ENABLE_TESTING()
add_subdirectory(TEST)
add_subdirectory(OPENSWATHALGO)
