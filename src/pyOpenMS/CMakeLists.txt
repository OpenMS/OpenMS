# --------------------------------------------------------------------------
#                   OpenMS -- Open-Source Mass Spectrometry
# --------------------------------------------------------------------------
# Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
# ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
cmake_minimum_required(VERSION 3.1 FATAL_ERROR)

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
# See https://wiki.python.org/moin/WindowsCompilers
# Windows support requires that the correct Python version is matched to the
# correct MSVS version:
# * Python 2.[67] needs to be compiled with MSVS 2008
# * Python 3.[34] needs to be compiled with MSVS 2010
# * Python 3.[5678] needs to be compiled with MSVS 2015
if(WIN32)
    if(SKIP_WIN_COMPILERCHECK)
        message(STATUS "Skipping MS Visual C++ compiler check (not recommended)")
    elseif(MSVC90 AND (${PY_VERSION} STREQUAL "2.6" OR ${PY_VERSION} STREQUAL "2.7" ) )
        message(STATUS "Need Visual C++ 2008 compiler for building Python 2.[67] extensions -- ok")
    elseif(MSVC10 AND (${PY_VERSION} STREQUAL "3.3" OR ${PY_VERSION} STREQUAL "3.4" ))
        message(STATUS "Need Visual C++ 2010 compiler for building Python 3.[34] extensions -- ok")
    elseif(MSVC14 AND (${PY_VERSION} STREQUAL "3.5" OR ${PY_VERSION} STREQUAL "3.6" OR ${PY_VERSION} STREQUAL "3.7" OR ${PY_VERSION} STREQUAL "3.8" ))
        message(STATUS "Need Visual C++ 2015 compiler for building Python 3.[5678] extensions -- ok")
    else()
        message(STATUS "Need Visual C++ 2008 compiler for building Python 2.[67] extensions")
        message(STATUS "Need Visual C++ 2010 compiler for building Python 3.[34] extensions")
        message(STATUS "Need Visual C++ 2015 compiler for building Python 3.[5678] extensions")
        message(STATUS "To skip this check, please configure with -DSKIP_WIN_COMPILERCHECK=On")
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
        ${PYTHON_EXECUTABLE} -c "import autowrap; exit(autowrap.version >= (0, 18, 1))"
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
        message(STATUS "Found autowrap version ${AUTOWRAP_VERSION}. The version is too old (>= 0.18.1 is required)")
        message(FATAL_ERROR "Please upgrade autowrap or disable pyOpenMS.")
    endif()
endif()


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

#-----------------------------------------------------------------------------
# since boost 1.69 there seem to be symbols visible/imported after linking it
# statically into OpenMS. Therefore we need to link to it for pyOpenMS as well.
# Just using the OpenMS dependencies is hard since the CMake variables are
# a) mangled together with generator expressions (debug and release)
# b) do not specify if it was a static or dynamic library (on Unix this is easy to
# test from the extension, but not on Windows)
# c) may include recursively imported targets like Qt which list Qt::Core instead
# of the actual path
# If we only could get rid of distutils in general and do everything with CMake..
if(NOT WIN32)
  find_boost(regex) ## variables are not global -> find again
  
  # the following does not work since the FindBoost does not annotate all target_properties,
  # when run without BOOST_USE_STATIC On or Off.
  #set(OBJECTS_EXTEND "\$<GENEX_EVAL:\$<\$<STREQUAL:\$<TARGET_LINKER_FILE_SUFFIX:Boost::regex>,\".a\">:\$<TARGET_LINKER_FILE:Boost::regex>>>")
  
  # the logic for static vs dynamic is in setup.py
  set(LIBRARIES_TO_BE_PARSED_EXTEND "\$<TARGET_LINKER_FILE:Boost::regex>")
  
  # the following gets all dependent libraries of OpenMS
  # for the current config but filters out imported targets unfortunately.
  # We would need something recursive
  #set(ALL_OPENMS_DEPENDENCIES "\$<FILTER:\$<TARGET_GENEX_EVAL:OpenMS,\$<TARGET_PROPERTY:OpenMS,LINK_LIBRARIES>>,INCLUDE,/>")

  file(GENERATE
       OUTPUT ${CMAKE_BINARY_DIR}/pyOpenMS/env.py
       INPUT ${CMAKE_BINARY_DIR}/pyOpenMS/env.py)
endif()

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

# set data variable nightly builds use repository last change date
# yyyy-mm-dd 08:04:17 +0200 -> yyyymmdd
string(REPLACE "-" "" OPENMS_GIT_LC_DATE_REPLACED ${OPENMS_GIT_LC_DATE})
set(OPENMS_GIT_LC_DATE_REPLACED_LIST "${OPENMS_GIT_LC_DATE_REPLACED}")
separate_arguments(OPENMS_GIT_LC_DATE_REPLACED_LIST)
list(GET OPENMS_GIT_LC_DATE_REPLACED_LIST 0 OPENMS_GIT_LC_DATE_FORMAT)

###debug
set(OPENMS_GIT_LC_DATE ${OPENMS_GIT_LC_DATE})
###debug 

set(OPENMS_GIT_LC_DATE_FORMAT ${OPENMS_GIT_LC_DATE_FORMAT})

set(_env_py_in ${PROJECT_SOURCE_DIR}/env.py.in)
set(_env_py ${CMAKE_BINARY_DIR}/pyOpenMS/env.py)

#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# add additional sysroot information (osx)
if (APPLE)
  set(SYSROOT_OSX_PATH ${CMAKE_OSX_SYSROOT})
endif()
#------------------------------------------------------------------------------

#------------------------------------------------------------------------------
# collection of include dirs that is necessary to compile pyOpenMS

# Find Qt5 includes for pyOpenMS
find_package(Qt5 COMPONENTS Core Network REQUIRED)
get_target_property(QT_QMAKE_EXECUTABLE Qt5::qmake IMPORTED_LOCATION)
exec_program(${QT_QMAKE_EXECUTABLE} ARGS "-query QT_INSTALL_HEADERS" OUTPUT_VARIABLE QTHEADERS)
set(QT_INCLUDE_DIR ${QTHEADERS} CACHE INTERNAL "" FORCE)
exec_program(${QT_QMAKE_EXECUTABLE} ARGS "-v" OUTPUT_VARIABLE QT_QMAKE_VERSION_INFO)
exec_program(${QT_QMAKE_EXECUTABLE} ARGS "-query QT_INSTALL_LIBS" OUTPUT_VARIABLE QT_INSTALL_LIBS)
exec_program(${QT_QMAKE_EXECUTABLE} ARGS "-query QT_INSTALL_BINS" OUTPUT_VARIABLE QT_INSTALL_BINS)

set(PYOPENMS_INCLUDE_DIRS
  ${OpenSwathAlgo_INCLUDE_DIRECTORIES}
  ${OpenMS_INCLUDE_DIRECTORIES}
  ${SuperHirn_INCLUDE_DIRECTORIES}
  ${QT_INCLUDE_DIR}
)

list(REMOVE_DUPLICATES PYOPENMS_INCLUDE_DIRS)

set(CONTRIB_DIR ${CMAKE_PREFIX_PATH})
set(OPEN_MS_BUILD_TYPE ${CMAKE_BUILD_TYPE})

add_custom_target(
  prepare_pyopenms_libs
  DEPENDS OpenMS SuperHirn
)

# assemble the libraries

# copy the direct dependencies
set(PYOPENMS_DEPENDENCIES OpenMS OpenSwathAlgo SuperHirn)

foreach (PYOPENMS_DEPENDENCY ${PYOPENMS_DEPENDENCIES})
	add_custom_command(
		TARGET prepare_pyopenms_libs POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${PYOPENMS_DEPENDENCY}> ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms
	)
endforeach()
	
# fixup these dependencies (i.e. get other dynamic dependencies recursively)
if (APPLE) ## On APPLE use our script because the executables need to be relinked
## this is done before setup.py so that we fixup the resulting pyopenms.so's too
else()
    ## Assemble common required non-system libraries
    ## Note that we do not need the QT plugins or QTGui libraries since we do not include GUI tools here.
    foreach (PYOPENMS_DEPENDENCY ${PYOPENMS_DEPENDENCIES})
        add_custom_command(
            TARGET prepare_pyopenms_libs POST_BUILD
            COMMAND ${CMAKE_COMMAND} -DDEPS="$<TARGET_FILE:${PYOPENMS_DEPENDENCY}>" -DTARGET="${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/" -DLOOKUP_DIRS="${OPENMS_CONTRIB_LIBS}/lib\;${QT_INSTALL_BINS}\;${QT_INSTALL_LIBS}" -P ${PROJECT_SOURCE_DIR}/pyopenms_copy_deps.cmake
        )
    endforeach()

    if(WIN32)
        # copy all runtime files (do not preserve permissions, due to issues when deleting them from tmp folders)
        foreach (_runtime ${MSVS_RTLIBS})
            file(COPY ${_runtime} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms NO_SOURCE_PERMISSIONS)
        endforeach()
    endif()
endif()

# write configured variables into env.py
configure_file(${_env_py_in} ${_env_py} @ONLY)

#------------------------------------------------------------------------------
# create targets in makefile

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND PY_NO_OUTPUT)
  add_custom_target(pyopenms_create_cpp
    COMMAND ${PYTHON_EXECUTABLE} create_cpp_extension.py 2> /dev/null
    DEPENDS prepare_pyopenms_libs
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
ELSE()
  add_custom_target(pyopenms_create_cpp
    COMMAND ${PYTHON_EXECUTABLE} create_cpp_extension.py
    DEPENDS prepare_pyopenms_libs
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
ENDIF()


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


IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND PY_NO_OUTPUT)
  add_custom_target(pyopenms_build
            COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext ${PY_EXTRA_ARGS} 2> /dev/null
            DEPENDS pyopenms_create_cpp
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
ELSE()
  add_custom_target(pyopenms_build
            COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext ${PY_EXTRA_ARGS}
            DEPENDS pyopenms_create_cpp
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
ENDIF()

add_dependencies(pyopenms_build OpenMS)

if(APPLE)
add_custom_command(
  TARGET pyopenms_create_cpp POST_BUILD
  COMMAND ${CMAKE_SOURCE_DIR}/cmake/MacOSX/fix_dependencies.rb -l ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms -e "@rpath/" -v
  )
endif()

IF(${CMAKE_SYSTEM_NAME} MATCHES "Linux" AND PY_NO_OUTPUT)
  add_custom_target(pyopenms
            COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_egg 2> /dev/null
            COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_wheel 2> /dev/null
            COMMAND ${PYTHON_EXECUTABLE} setup.py bdist --format=zip 2> /dev/null
            COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext --inplace 2> /dev/null
            DEPENDS pyopenms_build
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
ELSE()
  add_custom_target(pyopenms
            COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_egg
            COMMAND ${PYTHON_EXECUTABLE} setup.py bdist_wheel
            COMMAND ${PYTHON_EXECUTABLE} setup.py bdist --format=zip
            COMMAND ${PYTHON_EXECUTABLE} setup.py build_ext --inplace
            DEPENDS pyopenms_build
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
ENDIF()


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
  test_tutorial.py
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

