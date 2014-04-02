# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2012.
#
# This software is released under a three-clause BSD license:
#  * Redistributions of source code must retain the above copyright
#    notice, this list of conditions and the following disclaimer.
#  * Redistributions in binary form must reproduce the above copyright
#    notice, this list of conditions and the following disclaimer in the
#    documentation and/or other materials provided with the distribution.
#  * Neither the name of any author or any participating institution
#    may be used to endorse or promote products derived from this software
#    without specific prior written permission.
# For a full list of authors, refer to the file AUTHORS.
# --------------------------------------------------------------------------
# THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
# AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
# INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
# EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
# PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
# OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
# WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
# OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
# ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
#
# --------------------------------------------------------------------------
# $Maintainer: Hannes Röst $
# $Authors: Hannes Röst, Uwe Schmitt $
# --------------------------------------------------------------------------

#------------------------------------------------------------------------------
# find and handle python
find_package(PythonInterp REQUIRED)

# find out python version info
execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import sys; print '%s.%s' % sys.version_info[:2]"
     OUTPUT_VARIABLE PY_VERSION
     OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Python found at ${PYTHON_EXECUTABLE} with version ${PY_VERSION} (if this is wrong, configure with -DPYTHON_EXECUTABLE:FILEPATH=/path/to/python)")

#------------------------------------------------------------------------------
# Windows support restricted to Python 2.7 at the moment!
#  * we also require VC 2008 and the corresponding vcredist to be installed
if(WIN32)
    if(NOT MSVC90)
        message(STATUS "Need Visual C++ 2008 compiler for building Python 2.[67] extensions")
        message(FATAL_ERROR "Either reconfigure with Visual Studio 9 2008 generator or disable pyOpenMS.")
    endif()

    include(InstallRequiredSystemLibraries)
    set(MSVCR90DLL ${MSVC90_CRT_DIR}/msvcr90.dll)
    if(NOT EXISTS ${MSVCR90DLL})
        message(STATUS "Missing msvcr90.dll - Visual C++ 2008 Runtime (called vcredist)")
        message(FATAL_ERROR "Please install VC 2008 runtime package or disable pyOpenMS.")
    endif()
    set(MSVCP90DLL ${MSVC90_CRT_DIR}/msvcp90.dll)
    if(NOT EXISTS ${MSVCP90DLL})
        message(STATUS "Missing msvcp90.dll - Visual C++ 2008 Runtime (called vcredist)")
        message(FATAL_ERROR "Please install VC 2008 runtime package or disable pyOpenMS.")
    endif()
endif(WIN32)

#------------------------------------------------------------------------------
# Find Cython
find_program( CYTHON_EXECUTABLE NAMES cython )

set(CYTHON_MISSING FALSE)
if(DEFINED CYTHON_EXECUTABLE-NOTFOUND)
	set(CYTHON_MISSING TRUE)
endif()

if(CYTHON_MISSING)
	message(FATAL_ERROR "Looking for cython - not found")
else()
  # find out cython version info
  execute_process(
       COMMAND
       ${PYTHON_EXECUTABLE} -c "import Cython; print Cython.__version__"
       OUTPUT_VARIABLE CYTHON_VERSION
       OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  message(STATUS "Looking for cython - found version ${CYTHON_VERSION}")
endif()

#------------------------------------------------------------------------------
# Check for autowrap Cython
execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import autowrap"
     RESULT_VARIABLE AUTOWRAP_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)

set(AUTOWRAP_VERSION_OK FALSE)

if(AUTOWRAP_MISSING)
	message(STATUS "Looking for autowrap - not found")
else()
    execute_process(
        COMMAND
        ${PYTHON_EXECUTABLE} -c "import autowrap; exit(autowrap.version >= (0, 3, 2))"
        RESULT_VARIABLE _AUTOWRAP_VERSION_OK
        ERROR_QUIET
        OUTPUT_QUIET
    )
    execute_process(
        COMMAND
        ${PYTHON_EXECUTABLE} -c "import autowrap; print '%d.%d.%d' % (autowrap.version)"
        OUTPUT_VARIABLE AUTOWRAP_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(_AUTOWRAP_VERSION_OK)
        message(STATUS "Looking for autowrap - found autowrap ${AUTOWRAP_VERSION}, version ok")
        set(AUTOWRAP_VERSION_OK TRUE)
    else()
        message(STATUS "Found autowrap version ${AUTOWRAP_VERSION}. The version is to old (>= 0.3.2 is required)")
        message(FATAL_ERROR "Please upgrade autowrap or disable pyOpenMS.")
    endif()
endif()

#------------------------------------------------------------------------------
# get the qt version
execute_process(
    COMMAND ${QT_QMAKE_EXECUTABLE} -v
    OUTPUT_VARIABLE QT_QMAKE_VERSION_INFO
    )


#------------------------------------------------------------------------------
# Check for Nose Test Framework
execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import nose"
     RESULT_VARIABLE _NOSE_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)

set(NOSE_MISSING TRUE)
if( _NOSE_MISSING EQUAL 0)
    set(NOSE_MISSING FALSE)
endif()
if(NOSE_MISSING)
	message(FATAL_ERROR "Looking for nose testing framework - not found")
else()
	message(STATUS "Looking for nose testing framework - found")
endif()

#------------------------------------------------------------------------------
# Check for Numpy
execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import numpy"
     RESULT_VARIABLE _NUMPY_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)

set(NUMPY_MISSING TRUE)
if( _NUMPY_MISSING EQUAL 0)
  set(NUMPY_MISSING FALSE)
endif()
if(NUMPY_MISSING)
	message(FATAL_ERROR "Looking for numpy - not found")
else()
	message(STATUS "Looking for numpy - found")
endif()

#------------------------------------------------------------------------------
# Handle missing libraries (this should never be reached, as the individual
#  parts should fire FATAL_ERRORs if something is missing)
if(NUMPY_MISSING OR CYTHON_MISSING OR NOT AUTOWRAP_VERSION_OK OR NOSE_MISSING)
  message(FATAL_ERROR "Required Python modules not found or out of date")
endif()

#------------------------------------------------------------------------------
# clean python build directory from former cmake run (if exists)
# this can contain older versions of openms shared lib and might confuse
# the linker when working on pyopenms

FILE(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/pyOpenMS/build")
FILE(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/pyOpenMS/dist")
# OpenMS
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/OpenMSd.dll")
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/OpenMS.dll")
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libOpenMS.so")
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libOpenMS.dylib")
# OpenSwathAlgo
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/OpenSwathAlgod.dll")
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/OpenSwathAlgo.dll")
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libOpenSwathAlgo.so")
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libOpenSwathAlgo.dylib")
# SuperHirn
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/SuperHirnd.dll")
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/SuperHirn.dll")
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libSuperHirn.so")
FILE(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libSuperHirn.dylib")

#------------------------------------------------------------------------------
# copy files
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/unittests)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/memoryleaktests)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/integration_tests)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyTOPP)
FILE(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/extra_includes)

FILE(GLOB _python_files "${OPENMS_HOST_DIRECTORY}/pyOpenMS/pyopenms/*.py")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)

FILE(GLOB _python_files "${OPENMS_HOST_DIRECTORY}/pyOpenMS/pyopenms/*.sh")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)

FILE(GLOB _python_files "${OPENMS_HOST_DIRECTORY}/pyOpenMS/pyTOPP/*.py")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyTOPP)

FILE(GLOB _python_files "${OPENMS_HOST_DIRECTORY}/pyOpenMS/tests/unittests/*")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/tests/unittests)

FILE(GLOB _python_files "${OPENMS_HOST_DIRECTORY}/pyOpenMS/tests/*.mzXML")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/tests)

FILE(GLOB _python_files "${OPENMS_HOST_DIRECTORY}/pyOpenMS/tests/memoryleaktests/*")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/tests/memoryleaktests)

FILE(GLOB _python_files "${OPENMS_HOST_DIRECTORY}/pyOpenMS/tests/integration_tests/*")
FILE(COPY ${_python_files} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/tests/integration_tests)

# list of files required for the pyOpenMS build system
set(_pyopenms_files
  MANIFEST.in
  README.rst
  setup.py
  create_cpp_extension.py
  distribute_setup.py
  version.py
  run_nose.py
  run_memleaks.py
  doCythonCompileOnly.py
)

foreach(pyfile ${_pyopenms_files})
  configure_file(${OPENMS_HOST_DIRECTORY}/pyOpenMS/${pyfile} ${CMAKE_BINARY_DIR}/pyOpenMS/${pyfile} COPYONLY)
endforeach()

set(_sub_pyopenms_files
  License.txt
  version.py
)

# list of files located in pyOpenMS/pyopenms
foreach(pyfile ${_sub_pyopenms_files})
  configure_file(${OPENMS_HOST_DIRECTORY}/pyOpenMS/${pyfile} ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/${pyfile} COPYONLY)
endforeach()

if(WIN32)
    FILE(COPY ${MSVCR90DLL} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
    FILE(COPY ${MSVCP90DLL} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
    set(FOUND_XERCES FALSE)
    foreach(CONTRIB_PATH ${CONTRIB_DIR})
        if(EXISTS ${CONTRIB_PATH}/lib/xerces-c_3_1.dll)
            FILE(COPY ${CONTRIB_PATH}/lib/xerces-c_3_1.dll DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
            set(FOUND_XERCES TRUE)
        endif()
    endforeach()
    if(NOT FOUND_XERCES)
        message(STATUS "could not find xerces dll in contrib dir")
        RETURN()
    endif()
endif()

#------------------------------------------------------------------------------
# write variables for setup.py as Python script into pyOpenMS/env.py

set(_env_py_in ${OPENMS_HOST_DIRECTORY}/pyOpenMS/env.py.in)
set(_env_py ${CMAKE_BINARY_DIR}/pyOpenMS/env.py)

#------------------------------------------------------------------------------
# If there are other, external libraries that we would like to link, we can
# specify them here:
set(INCLUDE_DIRS_EXTEND "")
set(LIBRARIES_EXTEND "")
set(LIBRARY_DIRS_EXTEND "")

if (WITH_CRAWDAD)
  set(INCLUDE_DIRS_EXTEND "\"${CRAWDAD_INCLUDE_DIRS}\", \"${CRAWDAD_INCLUDE_DIRS}/msmat\", ${INCLUDE_DIRS_EXTEND}")
  set(LIBRARIES_EXTEND "\"Crawdad\", ${LIBRARIES_EXTEND} ${LIBRARIES_EXTEND}")
  set(LIBRARY_DIRS_EXTEND "\"${CRAWDAD_INCLUDE_DIRS}\", ${LIBRARY_DIRS_EXTEND}")
endif()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# collection of include dirs that is necessary to compile pyOpenMS
set(PYOPENMS_INCLUDE_DIRS
  ${OpenSwathAlgo_INCLUDE_DIRECTORIES}
  ${OpenMS_INCLUDE_DIRECTORIES}
  ${SuperHirn_INCLUDE_DIRECTORIES}
)

set(OPEN_MS_BUILD_TYPE ${CMAKE_BUILD_TYPE})
set(OPEN_MS_LIB "")
set(OPEN_SWATH_ALGO_LIB "")
set(SUPERHIRN_LIB "")

if(WIN32)
    if(CMAKE_BUILD_TYPE STREQUAL "Debug")
      set(OPEN_MS_LIB "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug/OpenMSd.dll")
      set(OPEN_SWATH_ALGO_LIB "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug/OpenSwathAlgod.dll")
      set(SUPERHIRN_LIB "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Debug/SuperHirnd.dll")
    else()
      set(OPEN_MS_LIB "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release/OpenMS.dll")
      set(OPEN_SWATH_ALGO_LIB "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release/OpenSwathAlgo.dll")
      set(SUPERHIRN_LIB "${CMAKE_RUNTIME_OUTPUT_DIRECTORY}/Release/SuperHirn.dll")
    endif()
endif()

configure_file(${_env_py_in} ${_env_py} @ONLY)

#------------------------------------------------------------------------------
# create targets in makefile

add_custom_target(pyopenms_create_cpp
	COMMAND ${PYTHON_EXECUTABLE} create_cpp_extension.py
	DEPENDS OpenMS SuperHirn
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )

add_custom_target(pyopenms
	COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext
	COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_egg
	COMMAND ${PYTHON_EXECUTABLE} setup.py bdist --format=zip
	COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext --inplace
	DEPENDS pyopenms_create_cpp
	WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )

add_dependencies(pyopenms OpenMS)

###########################################################################
#####                      Testing pyOpenMS                           #####
###########################################################################

# Original test using the "run_nose.py" script, testing all unittests at once
# => this is suboptimal for ctest and cdash because we dont see which tests
# actually have gone wrong. Thus we add additional tests below ...
enable_testing()
add_test(NAME test_pyopenms_unittests
         COMMAND ${PYTHON_EXECUTABLE} run_nose.py
         WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS
        )
if(NOT WIN32)
    set_tests_properties(test_pyopenms_unittests PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
endif()

# Please add your test here when you decide to write a new testfile in the tests/unittests folder
set(pyopenms_unittest_testfiles
test000.py
test_BaselineFiltering.py
test_ChromatogramExtractor.py
test_Convexhull.py
testCVTermList.py
test_DIAScoring.py
test_FileIO.py
test_Isobaric_Quantitation.py
testLightTargetedExperiment.py
test_MRMFeatureFinderScoring.py
test_MSSpectrumAndRichSpectrum.py
test_OpenSwathDataStructures.py
test_Smoothing.py
testSpecialCases.py
test_SpectraFilter.py
test_SpectrumAccessOpenMS.py
test_TraML.py
test_MzMLConsumer.py
test_MzXMLConsumer.py
)

# Please add your test here when you decide to write a new testfile in the tests/unittests folder
set(pyopenms_integrationtest_testfiles
test_MRMRTNormalizer.py
)

# Loop through all the test files
foreach (t ${pyopenms_unittest_testfiles})
  add_test(NAME "pyopenms_unittest_${t}"
    COMMAND ${PYTHON_EXECUTABLE} -c  "import nose; nose.run_exit()" ${CMAKE_BINARY_DIR}/pyOpenMS/tests/unittests/${t} -s -v)
  if(NOT WIN32)
    set_tests_properties("pyopenms_unittest_${t}" PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib"
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
  endif()
endforeach(t)

foreach (t ${pyopenms_integrationtest_testfiles})
  add_test(NAME "pyopenms_integrationtest_${t}"
    COMMAND ${PYTHON_EXECUTABLE} -c  "import nose; nose.run_exit()" ${CMAKE_BINARY_DIR}/pyOpenMS/tests/integration_tests/${t} -s -v)
  if(NOT WIN32)
    set_tests_properties("pyopenms_integrationtest_${t}" PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib"
      WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
  endif()
endforeach(t)

# Finally add the memory leaks test (in folder tests/memoryleaktests/)
add_test(NAME pyopenms_test_memoryleaktests
  COMMAND ${PYTHON_EXECUTABLE} -c  "import nose; nose.run_exit()" ${CMAKE_BINARY_DIR}/pyOpenMS/tests/memoryleaktests/ -s -v
  WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
if(NOT WIN32)
    set_tests_properties(pyopenms_test_memoryleaktests PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
endif()
