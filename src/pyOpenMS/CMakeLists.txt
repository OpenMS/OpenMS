# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
# $Authors: Hannes Röst, Uwe Schmitt, Stephan Aiche $
# --------------------------------------------------------------------------

project("pyOpenMS")
cmake_minimum_required(VERSION 2.8.3 FATAL_ERROR)

#------------------------------------------------------------------------------
# helper to copy/configure files from source to build
function(_copy_assets _asset_dir _regex _target_dir)
  file(GLOB _assets "${_asset_dir}/${_regex}")
  foreach(_asset_file ${_assets})
    # make path relative
    file(RELATIVE_PATH _relative_path ${_asset_dir} ${_asset_file})
    configure_file(${_asset_file} ${_target_dir}/${_relative_path} COPYONLY)
  endforeach()
endfunction()

#------------------------------------------------------------------------------
# find and handle python
find_package(PythonInterp REQUIRED)

# find out python version info
execute_process(
     COMMAND
     #${PYTHON_EXECUTABLE} -c "import sys; print('%s'% sys.version[:2])"
     ${PYTHON_EXECUTABLE} -c "from __future__ import print_function; import sys; print('%s.%s' % sys.version_info[:2])"
     OUTPUT_VARIABLE PY_VERSION
     OUTPUT_STRIP_TRAILING_WHITESPACE
)

message(STATUS "Python found at ${PYTHON_EXECUTABLE} with version ${PY_VERSION} (if this is wrong, configure with -DPYTHON_EXECUTABLE:FILEPATH=/path/to/python)")

#------------------------------------------------------------------------------
# Windows support requires that the correct Python version is matched to the correct MSVS version:
# * Python 2.[67] needs to be compiled with MSVS 2008
# * Python 3.[34] needs to be compiled with MSVS 2010
# * Python 3.[56] needs to be compiled with MSVS 2015
if(WIN32)
    if(MSVC90 AND (${PY_VERSION} STREQUAL "2.6" OR ${PY_VERSION} STREQUAL "2.7" ) )
        message(STATUS "Need Visual C++ 2008 compiler for building Python 2.[67] extensions -- ok")
    elseif(MSVC10 AND (${PY_VERSION} STREQUAL "3.3" OR ${PY_VERSION} STREQUAL "3.4" ))
        message(STATUS "Need Visual C++ 2010 compiler for building Python 3.[34] extensions -- ok")
    elseif(MSVC14 AND (${PY_VERSION} STREQUAL "3.5" OR ${PY_VERSION} STREQUAL "3.6" ))
        message(STATUS "Need Visual C++ 2015 compiler for building Python 3.[56] extensions -- ok")
    else()
        message(STATUS "Need Visual C++ 2008 compiler for building Python 2.[67] extensions")
        message(STATUS "Need Visual C++ 2010 compiler for building Python 3.[34] extensions")
        message(STATUS "Need Visual C++ 2015 compiler for building Python 3.[56] extensions")
        message(FATAL_ERROR "Either reconfigure with Visual Studio 2008, 2010 or 2015 generator, select a different Python version with -DPYTHON_EXECUTABLE:FILEPATH=/path/to/python or disable pyOpenMS.")
    endif()

    # Try to find MSVS runtime libs
    include(InstallRequiredSystemLibraries)

    if("${MSVS_RTLIBS}" STREQUAL "")
        set(MSVS_RTLIBS "${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}")
    endif()
     
    if("${MSVS_RTLIBS}" STREQUAL "")
        message(FATAL_ERROR "Did not find MSVS runtime, either provide with -DMSVS_RTLIBS='lib1;lib2' or disable pyOpenMS.")
    else()
    message(STATUS "Found MSVS runtime: ${MSVS_RTLIBS}")
    endif()

endif(WIN32)

#------------------------------------------------------------------------------
# Find Cython
find_program( CYTHON_EXECUTABLE NAMES cython )

set(CYTHON_MISSING FALSE)
if(DEFINED CYTHON_EXECUTABLE-NOTFOUND)
	set(CYTHON_MISSING TRUE)
endif()

# check minor version of Cython (0.xx.yy)
# Working versions (tested manually) are 0.25.2
set(CYTHON_VERSION_OK FALSE)
if(CYTHON_MISSING)
	message(FATAL_ERROR "Looking for cython - not found")
else()
  execute_process(
      COMMAND
      ${PYTHON_EXECUTABLE} -c "import Cython; exit(int(Cython.__version__.split('.')[1]) > 25 or (int(Cython.__version__.split('.')[1]) == 25 and int(Cython.__version__.split('.')[2]) >= 2))"
      RESULT_VARIABLE _CYTHON_VERSION_OK
      ERROR_QUIET
      OUTPUT_QUIET
  )
  execute_process(
      COMMAND
      ${PYTHON_EXECUTABLE} -c "from __future__ import print_function; import Cython; print (Cython.__version__)"
      OUTPUT_VARIABLE CYTHON_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
  )
  if(_CYTHON_VERSION_OK)
      message(STATUS "Looking for Cython - found Cython ${CYTHON_VERSION}, version ok")
      set(CYTHON_VERSION_OK TRUE)
  else()
      message(STATUS "Found Cython version ${CYTHON_VERSION}. The version is too old (>= 0.25.2 is required)")
      message(FATAL_ERROR "Please upgrade Cython or disable pyOpenMS.")
  endif()
endif()

#------------------------------------------------------------------------------
# Check for autowrap
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
        ${PYTHON_EXECUTABLE} -c "import autowrap; exit(autowrap.version >= (0, 13, 1))"
        RESULT_VARIABLE _AUTOWRAP_VERSION_OK
        ERROR_QUIET
        OUTPUT_QUIET
    )
    execute_process(
        COMMAND
        ${PYTHON_EXECUTABLE} -c "from __future__ import print_function; import autowrap; print ('%d.%d.%d' % (autowrap.version))"
        OUTPUT_VARIABLE AUTOWRAP_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(_AUTOWRAP_VERSION_OK)
        message(STATUS "Looking for autowrap - found autowrap ${AUTOWRAP_VERSION}, version ok")
        set(AUTOWRAP_VERSION_OK TRUE)
    else()
        message(STATUS "Found autowrap version ${AUTOWRAP_VERSION}. The version is too old (>= 0.13.1 is required)")
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
# Check for setuptools

execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import setuptools"
     RESULT_VARIABLE SETUPTOOLS_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)

set(SETUPTOOLS_VERSION_OK FALSE)

if(SETUPTOOLS_MISSING)
	message(STATUS "Looking for setuptools - not found")
else()
    execute_process(
        COMMAND
        ${PYTHON_EXECUTABLE} -c "import setuptools; exit(int(setuptools.__version__.split('.')[0]) >= 12)"
        RESULT_VARIABLE _SETUPTOOLS_VERSION_OK
        ERROR_QUIET
        OUTPUT_QUIET
    )
    execute_process(
        COMMAND
        ${PYTHON_EXECUTABLE} -c "from __future__ import print_function;import setuptools; print(setuptools.__version__)"
        OUTPUT_VARIABLE SETUPTOOLS_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(_SETUPTOOLS_VERSION_OK)
        message(STATUS "Looking for setuptools - found setuptools ${SETUPTOOLS_VERSION}, version ok")
        set(SETUPTOOLS_VERSION_OK TRUE)
    else()
        message(STATUS "Found setuptools version ${SETUPTOOLS_VERSION}. The version is too old (>= 12.0 is required)")
        message(FATAL_ERROR "Please upgrade setuptools or disable pyOpenMS.")
    endif()
endif()

#------------------------------------------------------------------------------
# Check for python wheel
execute_process(
     COMMAND
     ${PYTHON_EXECUTABLE} -c "import wheel"
     RESULT_VARIABLE WHEEL_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)

if(WHEEL_MISSING)
  message(FATAL_ERROR "Looking for python wheel - not found")
endif()

#------------------------------------------------------------------------------
# Handle missing libraries (this should never be reached, as the individual
#  parts should fire FATAL_ERRORs if something is missing)
if(NUMPY_MISSING OR CYTHON_MISSING OR NOT CYTHON_VERSION_OK OR NOT AUTOWRAP_VERSION_OK OR NOSE_MISSING
   OR SETUPTOOLS_MISSING OR WHEEL_MISSING)
  message(FATAL_ERROR "Required Python modules not found or out of date")
endif()

#------------------------------------------------------------------------------
# clean python build directory from former cmake run (if exists)
# this can contain older versions of openms shared lib and might confuse
# the linker when working on pyopenms

file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/pyOpenMS/build")
file(REMOVE_RECURSE "${CMAKE_BINARY_DIR}/pyOpenMS/dist")
# OpenMS
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/OpenMSd.dll")
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/OpenMS.dll")
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libOpenMS.so")
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libOpenMS.dylib")
# OpenSwathAlgo
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/OpenSwathAlgod.dll")
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/OpenSwathAlgo.dll")
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libOpenSwathAlgo.so")
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libOpenSwathAlgo.dylib")
# SuperHirn
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/SuperHirnd.dll")
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/SuperHirn.dll")
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libSuperHirn.so")
file(REMOVE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/libSuperHirn.dylib")

#------------------------------------------------------------------------------
# copy/configure files
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/unittests)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/memoryleaktests)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/integration_tests)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/share)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyTOPP)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/extra_includes)

_copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "*.py" ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/)
_copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "*.sh" ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/)
_copy_assets("${PROJECT_SOURCE_DIR}/pyTOPP/" "*.py" ${CMAKE_BINARY_DIR}/pyOpenMS/pyTOPP/)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/unittests/" "*" ${CMAKE_BINARY_DIR}/pyOpenMS/tests/unittests)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/" "*.mzXML" ${CMAKE_BINARY_DIR}/pyOpenMS/tests)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/memoryleaktests/" "*" ${CMAKE_BINARY_DIR}/pyOpenMS/tests/memoryleaktests)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/integration_tests/" "*" ${CMAKE_BINARY_DIR}/pyOpenMS/tests/integration_tests)
file(COPY ${PROJECT_SOURCE_DIR}/../../share DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/)

# list of files required for the pyOpenMS build system
set(_pyopenms_files
  MANIFEST.in
  README.rst
  setup.py
  create_cpp_extension.py
  version.py
  run_nose.py
  run_memleaks.py
  doCythonCompileOnly.py
  fix_pyopenms_dependencies_on_mac.sh
)

foreach(pyfile ${_pyopenms_files})
  configure_file(${PROJECT_SOURCE_DIR}/${pyfile} ${CMAKE_BINARY_DIR}/pyOpenMS/${pyfile} COPYONLY)
endforeach()

set(_sub_pyopenms_files
  License.txt
  version.py
)

# list of files located in pyOpenMS/pyopenms
foreach(pyfile ${_sub_pyopenms_files})
  configure_file(${PROJECT_SOURCE_DIR}/${pyfile} ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/${pyfile} COPYONLY)
endforeach()

if(WIN32)
    # copy all runtime files (dont preserve permissions, due to issues when deleting them from tmp folders)
    foreach (_runtime ${MSVS_RTLIBS})
      file(COPY ${_runtime} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms NO_SOURCE_PERMISSIONS)
    endforeach()
    
    set(FOUND_XERCES FALSE)
    foreach(CONTRIB_PATH ${CMAKE_PREFIX_PATH})
        if(EXISTS ${CONTRIB_PATH}/lib/xerces-c_3_1.dll)
            file(COPY ${CONTRIB_PATH}/lib/xerces-c_3_1.dll DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
            set(FOUND_XERCES TRUE)
            message(STATUS "found xerces-c-3_1_dll")
        endif()
    endforeach()
    if(NOT FOUND_XERCES)
        message(FATAL_ERROR "could not find xerces-c_3_1.dll in contrib dir")
    endif()
endif()

find_boost()
#------------------------------------------------------------------------------
# write variables for setup.py as Python script into pyOpenMS/env.py
#  1 thread for compilation and using 8 modules seems like a reasonable number
#  for now -- note that you need change __init__.py if you change the number of
#  modules...
if(NOT PY_NUM_THREADS)
  set(PY_NUM_THREADS 1)
endif()
if(NOT PY_NUM_MODULES)
  set(PY_NUM_MODULES 8)
endif()
set(_env_py_in ${PROJECT_SOURCE_DIR}/env.py.in)
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

# TODO: this is a workaround; we need to clarify why exactly pyOpenMS needs
#       access to the directories and not the library locations itself
set(CONTRIB_DIR ${CMAKE_PREFIX_PATH})
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

set(PY_EXTRA_ARGS "")
if(PY_SINGLE_THREADED)
  message(STATUS "Turning off multi-threaded python module compilation")
  set(PY_EXTRA_ARGS "${PY_EXTRA_ARGS}" "--single-threaded")
endif()
if(PY_NO_OPTIMIZATION)
  message(STATUS "Turning off optimization for faster python module compile time")
  set(PY_EXTRA_ARGS "${PY_EXTRA_ARGS}" "--no-optimization")
endif()
message(STATUS "Py extra args ${PY_EXTRA_ARGS}")

IF(${CMAKE_SYSTEM_NAME} MATCHES "Darwin")
    add_custom_target(pyopenms
        COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext ${PY_EXTRA_ARGS}
        COMMAND ${CMAKE_BINARY_DIR}/pyOpenMS/fix_pyopenms_dependencies_on_mac.sh
        COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_egg
        COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_wheel
        COMMAND ${PYTHON_EXECUTABLE} setup.py bdist --format=zip
        COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext --inplace
        DEPENDS pyopenms_create_cpp
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
ELSE()
    add_custom_target(pyopenms 
        COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext ${PY_EXTRA_ARGS}
        COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_egg
        COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_wheel
        COMMAND ${PYTHON_EXECUTABLE} setup.py bdist --format=zip
        COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext --inplace
        DEPENDS pyopenms_create_cpp
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
ENDIF()

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
  test_ChromatogramExtractorAlgorithm.py
  test_Convexhull.py
  testCVTermList.py
  test_DIAScoring.py
  test_FileIO.py
  test_Isobaric_Quantitation.py
  testLightTargetedExperiment.py
  test_MRMFeatureFinderScoring.py
  test_MSNumpressCoder.py
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

# Please add your test here when you decide to write a new testfile in the tests/integration_tests folder
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
if(NOT PY_MEMLEAK_DISABLE)
  add_test(NAME pyopenms_test_memoryleaktests
    COMMAND ${PYTHON_EXECUTABLE} -c  "import nose; nose.run_exit()" ${CMAKE_BINARY_DIR}/pyOpenMS/tests/memoryleaktests/ -s -v
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
  if(NOT WIN32)
      set_tests_properties(pyopenms_test_memoryleaktests PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
  endif()
endif()

