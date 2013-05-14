cmake_minimum_required(VERSION 2.8)
project( PYOPENMS )

find_package(PythonInterp REQUIRED)

find_program( CYTHON_EXECUTABLE NAMES cython )

include( FindPackageHandleStandardArgs )
FIND_PACKAGE_HANDLE_STANDARD_ARGS( Cython REQUIRED_VARS CYTHON_EXECUTABLE )

mark_as_advanced( CYTHON_EXECUTABLE )

# copy files

FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyopenms)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pxds)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/addons)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/converters)

FILE(GLOB _python_files "pyopenms/*.py")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyopenms)

FILE(GLOB _python_files "pxds/*.pxd")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pxds)

FILE(GLOB _python_files "addons/*.pyx")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/addons)

FILE(GLOB _python_files "converters/*.py")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/converters)

FILE(COPY setup.py DESTINATION ${CMAKE_BINARY_DIR})

MESSAGE("PYOK " ${PYOK})

add_custom_target(pyopenms COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_egg)
