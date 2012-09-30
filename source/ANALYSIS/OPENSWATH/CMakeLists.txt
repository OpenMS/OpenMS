cmake_minimum_required (VERSION 2.6)

project(OpenSwath)

# set( BinOut ${CMAKE_SOURCE_DIR}/${CMAKE_PROJECT_NAME}out/${} )
set( TestBinOut ${CMAKE_SOURCE_DIR}/source/TEST/bin/${} )
set( LibOut ${CMAKE_SOURCE_DIR}/lib/${} )
SET (CMAKE_LIBRARY_OUTPUT_DIRECTORY ${LibOut} )
SET (CMAKE_RUNTIME_OUTPUT_DIRECTORY ${TestBinOut} )
SET (CMAKE_ARCHIVE_OUTPUT_DIRECTORY ${LibOut} )
SET(inc ${inc} ${CMAKE_SOURCE_DIR}/include)

add_subdirectory(OPENSWATHALGO)

