#IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
    #IF (WIN32)
        #MESSAGE(STATUS "bulding debug version on Windows not supported yet")
        #RETURN()
    #ENDIF()
#ENDIF()

find_package(PythonInterp REQUIRED)

# find out python version info
execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import sys; print '%s.%s' % sys.version_info[:2]"
     OUTPUT_VARIABLE PY_VERSION
     OUTPUT_STRIP_TRAILING_WHITESPACE
)

MESSAGE(STATUS "found python ${PY_VERSION}")

# windows support restrecited to pyhton2.7 at the moment !
IF (WIN32)
    IF (NOT PY_VERSION STREQUAL "2.7")
        MESSAGE(STATUS "need python 2.7 on windows")
        RETURN()
    ENDIF()

    IF (NOT MSVC90)
        MESSAGE(STATUS "need visual c++ 2008 compiler for building python 2.7 extensions")
        RETURN()
    ENDIF()


    include(InstallRequiredSystemLibraries)
    SET(MSVCR90DLL ${MSVC90_CRT_DIR}/msvcr90.dll)
    IF (NOT EXISTS ${MSVCR90DLL})
        MESSAGE(STATUS "missed msvcr90.dll - need visual c++ 2008 runtime (called vcredist)")
        RETURN()
    ENDIF()
    SET(MSVCP90DLL ${MSVC90_CRT_DIR}/msvcp90.dll)
    IF (NOT EXISTS ${MSVCP90DLL})
        MESSAGE(STATUS "missed msvcp90.dll - need visual c++ 2008 runtime (called vcredist)")
        RETURN()
    ENDIF()

ENDIF(WIN32)


MESSAGE(STATUS "Looking for cython")
find_program( CYTHON_EXECUTABLE NAMES cython )

SET(CYTHON-MISSING FALSE)
IF (DEFINED CYTHON_EXECUTABLE-NOTFOUND)
	SET(CYTHON-MISSING TRUE)
ENDIF()

IF (CYTHON-MISSING)
	MESSAGE(STATUS "Looking for cython - not found")
ELSE()
	MESSAGE(STATUS "Looking for cython - found")
ENDIF()

###### autwowrap check ########

MESSAGE(STATUS "Looking for autowrap")
execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import autowrap"
     RESULT_VARIABLE AUTOWRAP_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)

SET(AUTOWRAP-VERSION-OK FALSE)

IF(AUTOWRAP_MISSING)
	MESSAGE(STATUS "Looking for autowrap - not found")
ELSE()
	MESSAGE(STATUS "Looking for autowrap - found")
    execute_process(
        COMMAND
        ${PYTHON_EXECUTABLE} -c "import autowrap; exit(autowrap.version >= (0,2,7))"
        RESULT_VARIABLE AUTOWRAP_VERSION_OK
        ERROR_QUIET
        OUTPUT_QUIET
    )
    IF(AUTOWRAP_VERSION_OK)
        MESSAGE(STATUS "Looking for autowrap - version ok")
        SET(AUTOWRAP-VERSION-OK TRUE)
    ELSE()
        MESSAGE(STATUS "Looking for autowrap - version before 0.2.7, please upgrade")
    ENDIF()
ENDIF()

MESSAGE(STATUS "Looking for nose testing framework")
execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import nose"
     RESULT_VARIABLE NOSE_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)


SET(NOSE-MISSING TRUE)
IF( NOSE_MISSING EQUAL 0)
    SET(NOSE-MISSING FALSE)
ENDIF()
IF(NOSE_MISSING)
	MESSAGE(STATUS "Looking for nose testing framework - not found")
ELSE()
	MESSAGE(STATUS "Looking for nose testing framework - found")
ENDIF()


MESSAGE(STATUS "Looking for numpy")
execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import numpy"
     RESULT_VARIABLE NUMPY_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)

SET(NUMPY-MISSING TRUE)
IF( NUMPY_MISSING EQUAL 0)
  SET(NUMPY-MISSING FALSE)
ENDIF()
IF(NUMPY_MISSING)
	MESSAGE(STATUS "Looking for numpy - not found")
ELSE()
	MESSAGE(STATUS "Looking for numpy - found")
ENDIF()


IF (NUMPY-MISSING OR CYTHON-MISSING OR NOT AUTOWRAP-VERSION-OK OR NOSE-MISSING)
    MESSAGE(FATAL_ERROR "needed Python modules not found or out of date")
    RETURN()
ENDIF()


# copy files
# MESSAGE(STATUS ${CMAKE_BINARY_DIR}/pyOpenMS/tests/unittests)

FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/unittests)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/memoryleaktests)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/integration_tests)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyTOPP)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pxds)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/addons)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/converters)

FILE(GLOB _python_files "pyOpenMS/pyopenms/*.py")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)

FILE(GLOB _python_files "pyOpenMS/pyopenms/*.sh")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)

FILE(GLOB _python_files "pyOpenMS/pyTOPP/*.py")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyTOPP)

FILE(GLOB _python_files "pyOpenMS/tests/unittests/*")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/tests/unittests)

FILE(GLOB _python_files "pyOpenMS/tests/*.mzXML")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/tests)

FILE(GLOB _python_files "pyOpenMS/tests/memoryleaktests/*")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/tests/memoryleaktests)

FILE(GLOB _python_files "pyOpenMS/tests/integration_tests/*")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/tests/integration_tests)

FILE(GLOB _files "pyOpenMS/pxds/*.pxd")
FILE(COPY ${_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pxds)

FILE(GLOB _python_files "pyOpenMS/addons/*.pyx")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/addons)

FILE(GLOB _python_files "pyOpenMS/converters/*.py")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/converters)

FILE(COPY pyOpenMS/License.txt DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
FILE(COPY pyOpenMS/MANIFEST.in DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/)
FILE(COPY pyOpenMS/setup.py DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS)
FILE(COPY pyOpenMS/distribute_setup.py DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS)
FILE(COPY pyOpenMS/version.py DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS)
FILE(COPY pyOpenMS/version.py DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
FILE(COPY pyOpenMS/run_nose.py DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS)
FILE(COPY pyOpenMS/run_memleaks.py DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS)

IF (WIN32)
    FILE(COPY ${MSVCR90DLL} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
    FILE(COPY ${MSVCP90DLL} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
    SET(FOUND_XERCES FALSE)
    FOREACH(CONTRIB_PATH ${CONTRIB_DIR})
        IF (EXISTS ${CONTRIB_PATH}/lib/xerces-c_3_0.dll)
            FILE(COPY ${CONTRIB_PATH}/lib/xerces-c_3_0.dll DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
            SET(FOUND_XERCES TRUE)
        ENDIF()
    ENDFOREACH()
    IF (NOT FOUND_XERCES)
        MESSAGE(STATUS "cound not find xerces dll in contrib dir")
        RETURN()
    ENDIF()
ENDIF()


# write variables for setup.py as python script

set(ENVPATH ${CMAKE_BINARY_DIR}/pyOpenMS/env.py)

FILE(WRITE ${ENVPATH} OPEN_MS_SRC="${CMAKE_SOURCE_DIR}" "\n" )
FILE(APPEND ${ENVPATH} OPEN_MS_BUILD_DIR="${CMAKE_BINARY_DIR}" "\n")

FILE(APPEND ${ENVPATH} OPEN_MS_CONTRIB_BUILD_DIRS=\")
FOREACH(CONTRIB_PATH ${CONTRIB_DIR})
	FILE(APPEND ${ENVPATH} ${CONTRIB_PATH} ";")
ENDFOREACH()
FILE(APPEND ${ENVPATH} "\"\n")

FILE(APPEND ${ENVPATH} QT_HEADERS_DIR="${QT_HEADERS_DIR}" "\n")
FILE(APPEND ${ENVPATH} QT_LIBRARY_DIR="${QT_LIBRARY_DIR}" "\n")
FILE(APPEND ${ENVPATH} QT_QTCORE_INCLUDE_DIR="${QT_QTCORE_INCLUDE_DIR}" "\n")
FILE(APPEND ${ENVPATH} MSVCR90DLL="${MSVCR90DLL}" "\n")
FILE(APPEND ${ENVPATH} MSVCP90DLL="${MSVCP90DLL}" "\n")
FILE(APPEND ${ENVPATH} OPEN_MS_BUILD_TYPE="${CMAKE_BUILD_TYPE}" "\n")
IF (WIN32)
    IF (CMAKE_BUILD_TYPE STREQUAL "Debug")
        FILE(APPEND ${ENVPATH} OPEN_MS_LIB="${OpenMS_BINARY_DIR}/bin/OpenMSd.dll" "\n")
        FILE(APPEND ${ENVPATH} OPEN_SWATH_ALGO_LIB="${OpenMS_BINARY_DIR}/bin/OpenSwathAlgod.dll" "\n")
    ELSE()
        FILE(APPEND ${ENVPATH} OPEN_MS_LIB="${OpenMS_BINARY_DIR}/bin/OpenMS.dll" "\n")
        FILE(APPEND ${ENVPATH} OPEN_SWATH_ALGO_LIB="${OpenMS_BINARY_DIR}/bin/OpenSwathAlgo.dll" "\n")
    ENDIF()
ENDIF()

# create targets in makefile 
add_custom_target(pyopenms
	COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext --inplace
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
add_dependencies(pyopenms OpenMS)

add_custom_target(pyopenms_bdist_egg 
	COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_egg
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
add_dependencies(pyopenms_bdist_egg OpenMS)

add_custom_target(pyopenms_bdist 
	COMMAND ${PYTHON_EXECUTABLE} setup.py bdist  --format=zip
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
add_dependencies(pyopenms_bdist OpenMS)

add_custom_target(pyopenms_rpm 
	COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_rpm  
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
add_dependencies(pyopenms_rpm OpenMS)


enable_testing()
add_test(NAME test_pyopenms 
         COMMAND ${PYTHON_EXECUTABLE} run_nose.py
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS 
        )
IF(NOT WIN32)
    set_tests_properties(test_pyopenms PROPERTIES ENVIRONMENT 
                         "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
                     #set_target_properties(pyopenms PROPERTIES ENVIRONMENT 
                     #"LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
ENDIF()
