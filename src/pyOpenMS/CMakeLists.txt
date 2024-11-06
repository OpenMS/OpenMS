# Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Hannes Röst $
# $Authors: Hannes Röst, Uwe Schmitt, Stephan Aiche $
# --------------------------------------------------------------------------
## 3.15 to enable new find policies for FindPython (note the discussion in https://gitlab.kitware.com/cmake/cmake/-/issues/21186)
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("pyOpenMS" LANGUAGES C CXX)

if(NOT TARGET OpenMS)
  find_package(OpenMS)
  #TODO when separating pyopenms from OpenMS, specify a pyopenms version and the compatible OpenMS version to be searched above.
  # For now we just hope that the right OpenMS version from this same git repo is installed and found and
  # therefore give pyopenms here the same version. CF_OPENMS_PACKAGE_VERSION is used in env.py for setup.py
  set(CF_OPENMS_PACKAGE_VERSION ${OpenMS_VERSION})
  #TODO create own module copy under pyopenms or somehow point to installed copy based on the found OpenMS variables
  include(../../cmake/Modules/GetGitRevisionDescription.cmake)
  git_short_info(OPENMS_GIT_SHORT_REFSPEC OPENMS_GIT_SHORT_SHA1 OPENMS_GIT_LC_DATE)
endif()

option(NO_DEPENDENCIES "Do not bundle dependencies into pyOpenMS build folder and wheel" OFF)
option(NO_SHARE "Do not bundle share into pyOpenMS build folder and wheel" OFF)

#------------------------------------------------------------------------------
# helper to copy/configure files from source to build
function(_copy_assets _asset_dir _regex _target_dir)
  file(GLOB _assets "${_asset_dir}/${_regex}")
  foreach(_asset_file ${_assets})
    # make path relative
    file(RELATIVE_PATH _relative_path ${_asset_dir} ${_asset_file})
    # check if the asset is a directory
    if(NOT IS_DIRECTORY ${_asset_file})
      configure_file(${_asset_file} ${_target_dir}/${_relative_path} COPYONLY)
    endif()
  endforeach()
endfunction()

# for backwards compatibility
if (PYTHON_EXECUTABLE)
  set(Python_EXECUTABLE ${PYTHON_EXECUTABLE})
endif()

#------------------------------------------------------------------------------
# find and handle python (Development for macros to build python modules)
find_package(Python COMPONENTS Interpreter Development.Module)

message(STATUS "Python found at ${Python_EXECUTABLE} with version ${Python_VERSION} (if this is wrong, see https://cmake.org/cmake/help/latest/module/FindPython.html and configure with e.g. -DPython_EXECUTABLE:=/path/to/python)")

if(UNIX AND NOT APPLE AND NOT CYGWIN)
    set(LINUX TRUE)
endif()

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
    elseif(MSVC90 AND (${Python_VERSION} STREQUAL "2.6" OR ${Python_VERSION} STREQUAL "2.7" ))
        message(STATUS "Need Visual C++ 2008 compiler for building Python 2.[67] extensions -- ok")
    elseif(MSVC10 AND (${Python_VERSION} STREQUAL "3.3" OR ${Python_VERSION} STREQUAL "3.4" ))
        message(STATUS "Need Visual C++ 2010 compiler for building Python 3.[34] extensions -- ok")
    elseif(MSVC14 AND (${Python_VERSION} VERSION_GREATER_EQUAL "3.5" ))
        message(STATUS "Need Visual C++ 2015/2017/2019 compiler for building Python 3.5+ extensions -- ok")
    else()
        message(STATUS "Need Visual C++ 2008 compiler for building Python 2.[67] extensions")
        message(STATUS "Need Visual C++ 2010 compiler for building Python 3.[34] extensions")
        message(STATUS "Need Visual C++ 2015/2017/2019 compiler for building Python 3.5+ extensions")
        message(STATUS "To skip this check, please configure with -DSKIP_WIN_COMPILERCHECK=On")
        message(FATAL_ERROR "Either reconfigure with the correct Visual Studio generator, select a different Python version with -DPython_EXECUTABLE:FILEPATH=/path/to/python or disable pyOpenMS.")
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
# Find Cython as module
execute_process(
  COMMAND
  ${Python_EXECUTABLE} -c "import Cython"
  RESULT_VARIABLE CYTHON_MODULE_NOTFOUND
  ERROR_QUIET
  OUTPUT_QUIET
)

if (CYTHON_MODULE_NOTFOUND)
    message(FATAL_ERROR "Cython import failed for the given Python executable.")
    set(CYTHON_MODULE_VERSION_CHECK_OK 0)
else()
    # Working versions (tested manually) are 0.25.2+
    execute_process(
      COMMAND
      ${Python_EXECUTABLE} -c "import Cython; exit(int(Cython.__version__.split('.')[0]) >= 3 or int(Cython.__version__.split('.')[1]) > 25 or (int(Cython.__version__.split('.')[1]) == 25 and int(Cython.__version__.split('.')[2]) >= 2))"
      RESULT_VARIABLE CYTHON_MODULE_VERSION_CHECK_OK
      ERROR_QUIET
      OUTPUT_QUIET
    )

    if(CYTHON_MODULE_VERSION_CHECK_OK)
        execute_process(
          COMMAND
          ${Python_EXECUTABLE} -c "from __future__ import print_function; import Cython; print (Cython.__version__)"
          OUTPUT_VARIABLE CYTHON_VERSION
          OUTPUT_STRIP_TRAILING_WHITESPACE
        )
        message(STATUS "Looking for Cython - found Cython version ${CYTHON_VERSION}. Version ok")
    else()
        message(STATUS "Looking for Cython - found Cython version ${CYTHON_VERSION}. The version is too old (>= 0.25.2 is required)")
        message(FATAL_ERROR "Please upgrade Cython or disable pyOpenMS.")
    endif()
endif()

#------------------------------------------------------------------------------
# Check for autowrap
execute_process(
     COMMAND
     ${Python_EXECUTABLE} -c "import autowrap"
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
        ${Python_EXECUTABLE} -c "import autowrap; exit(autowrap.version >= (0, 22, 11))"
        RESULT_VARIABLE _AUTOWRAP_VERSION_OK
        ERROR_QUIET
        OUTPUT_QUIET
    )
    execute_process(
        COMMAND
        ${Python_EXECUTABLE} -c "from __future__ import print_function; import autowrap; print ('%d.%d.%d' % (autowrap.version))"
        OUTPUT_VARIABLE AUTOWRAP_VERSION
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    execute_process(
        COMMAND
        ${Python_EXECUTABLE} -c "from __future__ import print_function; import autowrap; print (''.join([autowrap.__path__[0],'/data_files/autowrap/']))"
        OUTPUT_VARIABLE AUTOWRAP_INCLUDES
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    if(_AUTOWRAP_VERSION_OK)
        message(STATUS "Looking for autowrap - found autowrap ${AUTOWRAP_VERSION}, version ok")
        set(AUTOWRAP_VERSION_OK TRUE)
    else()
        message(STATUS "Found autowrap version ${AUTOWRAP_VERSION}. The version is too old (>= 0.22.11 is required)")
        message(FATAL_ERROR "Please upgrade autowrap or disable pyOpenMS.")
    endif()
endif()


#------------------------------------------------------------------------------
# Check for pytest Test Framework
execute_process(
     COMMAND
     ${Python_EXECUTABLE} -c "import pytest"
     RESULT_VARIABLE _PYTEST_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)

set(PYTEST_MISSING TRUE)
if( _PYTEST_MISSING EQUAL 0)
    set(PYTEST_MISSING FALSE)
endif()
if(PYTEST_MISSING)
	message(WARNING "Looking for pytest testing framework - not found. Testing pyopenms will not work.")
else()
	message(STATUS "Looking for pytest testing framework - found")
endif()

#------------------------------------------------------------------------------
# Check for Numpy (TODO since CMake 3.17 this can be done in the same call with the FindPython module)
execute_process(
     COMMAND
     ${Python_EXECUTABLE} -c "import numpy"
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
     ${Python_EXECUTABLE} -c "import setuptools"
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
        ${Python_EXECUTABLE} -c "import setuptools; exit(int(setuptools.__version__.split('.')[0]) >= 12)"
        RESULT_VARIABLE _SETUPTOOLS_VERSION_OK
        ERROR_QUIET
        OUTPUT_QUIET
    )
    execute_process(
        COMMAND
        ${Python_EXECUTABLE} -c "from __future__ import print_function;import setuptools; print(setuptools.__version__)"
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
     ${Python_EXECUTABLE} -c "import wheel"
     RESULT_VARIABLE WHEEL_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)

if(WHEEL_MISSING)
  message(FATAL_ERROR "Looking for python wheel - not found")
endif()

#------------------------------------------------------------------------------
# Check for pandas which is required to run pyopenms tests
execute_process(
     COMMAND
     ${Python_EXECUTABLE} -c "import pandas"
     RESULT_VARIABLE PANDAS_MISSING
     ERROR_QUIET
     OUTPUT_QUIET
)

if(PANDAS_MISSING)
  message(FATAL_ERROR "Looking for pandas - not found")
endif()

#------------------------------------------------------------------------------
# Handle missing libraries (this should never be reached, as the individual
#  parts should fire FATAL_ERRORs if something is missing)
if(NUMPY_MISSING OR CYTHON_MODULE_NOTFOUND OR NOT CYTHON_MODULE_VERSION_CHECK_OK OR NOT AUTOWRAP_VERSION_OK
   OR SETUPTOOLS_MISSING OR WHEEL_MISSING OR PANDAS_MISSING)
  message(FATAL_ERROR "Required Python modules not found or out of date")
endif()
    
# Find NumPy headers
exec_program(${Python_EXECUTABLE}
    ARGS "-c \"import numpy; print(numpy.get_include())\""
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    RETURN_VALUE NUMPY_NOT_FOUND
    )
if(NUMPY_NOT_FOUND)
    message(FATAL_ERROR "NumPy headers not found")
endif()

# Find Python headers
exec_program(${Python_EXECUTABLE}
    ARGS "-c \"import sysconfig; print(sysconfig.get_paths()['include'])\""
    OUTPUT_VARIABLE PYTHON_INCLUDE_DIRS
    RETURN_VALUE PYTHON_INCLUDE_DIRS_NOT_FOUND
    )
if(PYTHON_INCLUDE_DIRS_NOT_FOUND)
    message(FATAL_ERROR "Python headers not found")
endif()


# Get suffix for python C(++) extensions
exec_program(${Python_EXECUTABLE}
    ARGS "-c \"import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))\""
    OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
    RETURN_VALUE PYTHON_EXT_SUFFIX_NOT_FOUND
    )
if(PYTHON_EXT_SUFFIX_NOT_FOUND)
    message(FATAL_ERROR "Could not determine suffix for Python C extensions.")
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

#------------------------------------------------------------------------------
# copy/configure files
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/unittests)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/memoryleaktests)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/tests/integration_tests)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyTOPP)
file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/extra_includes)

_copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "*.py" ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/)
_copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "*.sh" ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/)
_copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "*.pyi" ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/)
_copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "py.typed" ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/)
_copy_assets("${PROJECT_SOURCE_DIR}/pyTOPP/" "*.py" ${CMAKE_BINARY_DIR}/pyOpenMS/pyTOPP/)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/unittests/" "*" ${CMAKE_BINARY_DIR}/pyOpenMS/tests/unittests)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/" "*.mzXML" ${CMAKE_BINARY_DIR}/pyOpenMS/tests)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/memoryleaktests/" "*" ${CMAKE_BINARY_DIR}/pyOpenMS/tests/memoryleaktests)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/integration_tests/" "*" ${CMAKE_BINARY_DIR}/pyOpenMS/tests/integration_tests)


if(NOT NO_SHARE)
  file(MAKE_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/share)
  file(COPY ${OPENMS_SHARE_DIR} DESTINATION ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/share)
endif()

# list of files required for the pyOpenMS build system
set(_pyopenms_files
  MANIFEST.in
  README.rst
  setup.py
  create_cpp_extension.py
  mac_fix_dependencies.rb
  doCythonCompileOnly.py
)

foreach(pyfile ${_pyopenms_files})
  configure_file(${PROJECT_SOURCE_DIR}/${pyfile} ${CMAKE_BINARY_DIR}/pyOpenMS/${pyfile} COPYONLY)
endforeach()


#------------------------------------------------------------------------------
# If there are other, external libraries that we would like to link, we can
# specify them here:
set(INCLUDE_DIRS_EXTEND "")
set(LIBRARIES_EXTEND "")
set(LIBRARY_DIRS_EXTEND "")

if (WITH_CRAWDAD)
  set(INCLUDE_DIRS_EXTEND ${CRAWDAD_INCLUDE_DIRS} ${CRAWDAD_INCLUDE_DIRS}/msmat ${INCLUDE_DIRS_EXTEND})
  set(LIBRARIES_EXTEND "Crawdad" ${LIBRARIES_EXTEND})
  set(LIBRARY_DIRS_EXTEND ${CRAWDAD_INCLUDE_DIRS} ${LIBRARY_DIRS_EXTEND})
endif()

if (MT_ENABLE_OPENMP)
  find_package(OpenMP COMPONENTS CXX)
endif()
if (OpenMP_FOUND)
  list(APPEND INCLUDE_DIRS_EXTEND ${OpenMP_CXX_INCLUDE_DIR})
# Unfortunately this does not work, since CMake somehow calls the library on macOS "libomp".
# Although the library name should clearly be just "omp". We do it manually in setup.py now.
# Might need regular updates..
# foreach(LIB ${OpenMP_CXX_LIB_NAMES}) # such as gomp, pthreads, omp, ...
#   list(APPEND LIBRARIES_EXTEND ${LIB})
#   #TODO maybe we need to add library and include paths as well. Let's hope they are in already registered folders
# endforeach()
endif()

# TODO the next section and most of env.py should not be needed anymore if CMake based build works reliably.
#  It is commented out for now.
##-----------------------------------------------------------------------------
## since boost 1.69 there seem to be symbols visible/imported after linking it
## statically into OpenMS. Therefore we need to link to it for pyOpenMS as well.
## Just using the OpenMS dependencies is hard since the CMake variables are
## a) mangled together with generator expressions (debug and release)
## b) do not specify if it was a static or dynamic library (on Unix this is easy to
## test from the extension, but not on Windows)
## c) may include recursively imported targets like Qt which list Qt::Core instead
## of the actual path

#if(NOT WIN32)

  ## the following does not work since the FindBoost does not annotate all target_properties,
  ## when run without BOOST_USE_STATIC On or Off.
  ##set(OBJECTS_EXTEND "\$<GENEX_EVAL:\$<\$<STREQUAL:\$<TARGET_LINKER_FILE_SUFFIX:Boost::regex>,\".a\">:\$<TARGET_LINKER_FILE:Boost::regex>>>")
  
  ## the logic for static vs dynamic is in setup.py
  #set(LIBRARIES_TO_BE_PARSED_EXTEND "\$<TARGET_LINKER_FILE:Boost::regex>" "\$<TARGET_LINKER_FILE:XercesC::XercesC>" )
  
  ## the following gets all dependent libraries of OpenMS
  ## for the current config but filters out imported targets unfortunately.
  ## We would need something recursive
  ##set(ALL_OPENMS_DEPENDENCIES "\$<FILTER:\$<TARGET_GENEX_EVAL:OpenMS,\$<TARGET_PROPERTY:OpenMS,LINK_LIBRARIES>>,INCLUDE,/>")
#endif()

file(GENERATE
     OUTPUT ${CMAKE_BINARY_DIR}/pyOpenMS/env.py
     INPUT ${CMAKE_BINARY_DIR}/pyOpenMS/env.py)
       
#------------------------------------------------------------------------------
# write variables for setup.py as Python script into pyOpenMS/env.py
#  1 thread for compilation and using 8 modules seems like a reasonable number
#  for now
if(NOT PY_NUM_THREADS)
  set(PY_NUM_THREADS 1)
endif()
if(NOT PY_NUM_MODULES)
  set(PY_NUM_MODULES 8)
endif()

if ("${PY_NUM_MODULES}" GREATER "1")
    set(PYCPPFILES "")
    set(PYPYXFILES "")
    foreach(NUM RANGE 1 ${PY_NUM_MODULES})
      list(APPEND PYCPPFILES "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/_pyopenms_${NUM}.cpp")
      list(APPEND PYPYXFILES "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/_pyopenms_${NUM}.pyx")
    endforeach()
else()
    set(PYCPPFILES "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/_pyopenms.cpp")
    set(PYPYXFILES "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/_pyopenms.pyx")
endif()


# set data variable nightly builds use repository last change date
# yyyy-mm-dd 08:04:17 +0200 -> yyyymmdd
string(REPLACE "-" "" OPENMS_GIT_LC_DATE_REPLACED ${OPENMS_GIT_LC_DATE})
set(OPENMS_GIT_LC_DATE_REPLACED_LIST "${OPENMS_GIT_LC_DATE_REPLACED}")
separate_arguments(OPENMS_GIT_LC_DATE_REPLACED_LIST)
list(GET OPENMS_GIT_LC_DATE_REPLACED_LIST 0 OPENMS_GIT_LC_DATE_FORMAT)

set(_env_py_in ${PROJECT_SOURCE_DIR}/env.py.in)
set(_env_py ${CMAKE_BINARY_DIR}/pyOpenMS/env.py)

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
  ${QT_INCLUDE_DIR}
  ${PROJECT_SOURCE_DIR}/extra_includes
)

list(REMOVE_DUPLICATES PYOPENMS_INCLUDE_DIRS)

## use / instead of \ because a path might end in / and thus might generate invalid python code in env.py: r"C:\dev\contrib_build;c:\dev\Qt5.6.2_\5.6\msvc2015_64\"
file(TO_CMAKE_PATH "${CMAKE_PREFIX_PATH}" CONTRIB_DIR)
set(OPEN_MS_BUILD_TYPE ${CMAKE_BUILD_TYPE})

add_custom_target(
  prepare_pyopenms_libs
  DEPENDS OpenMS
)

if(NOT NO_DEPENDENCIES)
    # assemble the libraries

    # copy the direct dependencies
    set(PYOPENMS_DEPENDENCIES OpenMS OpenSwathAlgo)

    foreach (PYOPENMS_DEPENDENCY ${PYOPENMS_DEPENDENCIES})
    	add_custom_command(
    		TARGET prepare_pyopenms_libs POST_BUILD
    		COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${PYOPENMS_DEPENDENCY}> ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms
    	)
    endforeach()
    	
    # fixup these dependencies (i.e. get other dynamic dependencies recursively)
    if (APPLE) ## On APPLE use our script because in the libraries the "install_names" and "rpaths" need to be renamed
    ## this is done after "setup.py build_ext" but before setup.py bdist_wheel so that we fixup the resulting pyopenms.so's too
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
endif()

# write configured variables into env.py
configure_file(${_env_py_in} ${_env_py} @ONLY)


#------------------------------------------------------------------------------
# create targets in makefile

IF(LINUX AND PY_NO_OUTPUT)
  add_custom_command(
    OUTPUT ${PYCPPFILES} ${PYPYXFILES}
    DEPENDS OpenMS # fake dependency to see if a header file has changed (TODO should probably depend on pxds)
    COMMAND ${Python_EXECUTABLE} create_cpp_extension.py 2> /dev/null
    COMMENT "Creating C++ extension with autowrap and Cython"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
ELSE()
  add_custom_command(
    OUTPUT ${PYCPPFILES} ${PYPYXFILES}
    DEPENDS OpenMS # fake dependency to see if a header file has changed (TODO should probably depend on pxds)
    COMMAND ${Python_EXECUTABLE} create_cpp_extension.py
    COMMENT "Creating C++ extension with autowrap and Cython"
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
ENDIF()

# this acts as a mutex in make, to not have race conditions on the above command in the next parallel targets
add_custom_target(compile_pxds DEPENDS ${PYCPPFILES} ${PYPYXFILES})

set(PY_EXTRA_ARGS "")

if(PY_NO_OPTIMIZATION)
  message(STATUS "Turning off optimization for faster python module compile time")
  set(PY_EXTRA_ARGS "${PY_EXTRA_ARGS}" "--no-optimization")
endif()

message(STATUS "Py extra args ${PY_EXTRA_ARGS}")

set(PY_RPATH ${CMAKE_INSTALL_RPATH})
if (APPLE)
 # TODO maybe we need one for linux too
 list(APPEND PY_RPATH "@loader_path/")
endif()

# actually on my Mac I didnt need that but on GH actions it fails
if (LINUX OR APPLE) # see https://github.com/cython/cython/issues/3380
  set(VIS_PRESET "default")
  set(INLINES_HIDDEN 0)
else()
  set(VIS_PRESET "hidden")
  set(INLINES_HIDDEN 1)
endif()

if(LINUX OR APPLE)
  set(EXTRA_WARNINGS "-Wno-unused-member-function" "-Wno-builtin-macro-redefined" "-Wno-unused" "-Wno-zero-as-null-pointer-constant")
  if(${Python_VERSION} VERSION_GREATER_EQUAL "3.8" AND ${Python_VERSION} VERSION_LESS "3.9")
    list(APPEND EXTRA_WARNINGS "-Wno-deprecated-declarations") # https://github.com/cython/cython/issues/3474
  endif()
else()
  set(EXTRA_WARNINGS "")
endif()

if ("${PY_NUM_MODULES}" GREATER "1")
    add_custom_target(pyopenms_compile)
    foreach(nr RANGE 1 ${PY_NUM_MODULES})
      #list(APPEND PYOPENMS_SRC "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/_pyopenms_${nr}.cpp")
      Python_add_library(pyopenms_${nr} MODULE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/_pyopenms_${nr}.cpp" )
      target_link_libraries(pyopenms_${nr} PRIVATE OpenMS)
      target_compile_definitions(pyopenms_${nr} PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
      target_include_directories(pyopenms_${nr} SYSTEM PUBLIC ${PYOPENMS_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIR} ${AUTOWRAP_INCLUDES})
      target_compile_options(pyopenms_${nr} PRIVATE ${EXTRA_WARNINGS})
      set_target_properties(pyopenms_${nr} PROPERTIES
        OUTPUT_NAME _pyopenms_${nr}
        PREFIX "" # if you build a shared module it adds lib as prefix
        SUFFIX ${PYTHON_EXT_SUFFIX} # since CMake 3.17 you can add WITH_SOABI to Python_add_library
        CXX_VISIBILITY_PRESET ${VIS_PRESET}
        VISIBILITY_INLINES_HIDDEN ${INLINES_HIDDEN}
        INSTALL_RPATH "${PY_RPATH}"
        BUILD_WITH_INSTALL_RPATH 1
        BUILD_WITH_INSTALL_NAME_DIR 1
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
      )
      add_dependencies(pyopenms_${nr} compile_pxds) # this is necessary for make (not ninja), to avoid race conditions
      add_dependencies(pyopenms_compile pyopenms_${nr})
      add_custom_command(TARGET pyopenms_${nr}
                   POST_BUILD
                   COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/fix_pyi_imports.py ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/_pyopenms_${nr}.pyi
                   COMMENT "Fixing typing information for _pyopenms_${nr}")
    endforeach()
else()
    Python_add_library(pyopenms_compile MODULE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/_pyopenms.cpp" )
    target_link_libraries(pyopenms_compile PRIVATE OpenMS)
    target_compile_definitions(pyopenms_compile PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
    target_include_directories(pyopenms_compile SYSTEM PUBLIC ${PYOPENMS_INCLUDE_DIRS} ${PYTHON_INCLUDE_DIRS} ${NUMPY_INCLUDE_DIR} ${AUTOWRAP_INCLUDES})
    target_compile_options(pyopenms_compile PRIVATE ${EXTRA_WARNINGS})
    set_target_properties(pyopenms_compile PROPERTIES
        OUTPUT_NAME _pyopenms
        PREFIX "" # if you build a shared module it adds lib as prefix
        SUFFIX ${PYTHON_EXT_SUFFIX} # since CMake 3.17 you can add WITH_SOABI to Python_add_library
        CXX_VISIBILITY_PRESET ${VIS_PRESET}
        VISIBILITY_INLINES_HIDDEN ${INLINES_HIDDEN}
        INSTALL_RPATH "${PY_RPATH}"
        BUILD_WITH_INSTALL_RPATH 1
        BUILD_WITH_INSTALL_NAME_DIR 1
        ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        LIBRARY_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        RUNTIME_OUTPUT_DIRECTORY_RELEASE "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        LIBRARY_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
        RUNTIME_OUTPUT_DIRECTORY_DEBUG "${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms"
    )
    add_custom_command(TARGET pyopenms_compile
    POST_BUILD
    COMMAND ${Python_EXECUTABLE} ${PROJECT_SOURCE_DIR}/fix_pyi_imports.py ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/pyopenms.pyi #TODO Check if really necessary
    COMMENT "Fixing typing information for pyopenms")
endif()

## use no-binary=pyopenms since we built the binaries already with cmake
IF(LINUX AND PY_NO_OUTPUT)
    add_custom_target(pyopenms
            #COMMAND ${Python_EXECUTABLE} setup.py bdist_egg 2> /dev/null
            COMMAND ${Python_EXECUTABLE} -m pip wheel -w dist --no-deps --no-binary=pyopenms .
            #COMMAND ${Python_EXECUTABLE} setup.py bdist --format=zip 2> /dev/null
            #COMMAND ${Python_EXECUTABLE} setup.py build_ext --inplace 2> /dev/null
            COMMENT "Packaging pyopenms wheel."
            DEPENDS prepare_pyopenms_libs pyopenms_compile
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
ELSE()
    IF(APPLE AND NOT NO_DEPENDENCIES)
        add_custom_target(pyopenms
            COMMAND ./mac_fix_dependencies.rb -l ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms -e "@rpath/" -f -v
            # ad hoc signature is represented by a dash
            # https://developer.apple.com/documentation/security/seccodesignatureflags/kseccodesignatureadhoc
            COMMAND codesign --force --sign "-" ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/Qt*
            COMMAND codesign --force --sign "-" ${CMAKE_BINARY_DIR}/pyOpenMS/pyopenms/*.dylib
            #COMMAND ${Python_EXECUTABLE} setup.py bdist_egg
            COMMAND ${Python_EXECUTABLE} -m pip wheel -w dist --no-deps --no-binary=pyopenms .
            #COMMAND ${Python_EXECUTABLE} setup.py bdist --format=zip
            #COMMAND ${Python_EXECUTABLE} setup.py build_ext --inplace
            COMMENT "Packaging pyopenms wheel."
            DEPENDS prepare_pyopenms_libs pyopenms_compile
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
    ELSE()
        add_custom_target(pyopenms
            #COMMAND ${Python_EXECUTABLE} setup.py bdist_egg
            COMMAND ${Python_EXECUTABLE} -m pip wheel -w dist --no-deps --no-binary=pyopenms .
            #COMMAND ${Python_EXECUTABLE} setup.py bdist --format=zip
            #COMMAND ${Python_EXECUTABLE} setup.py build_ext --inplace
            COMMENT "Packaging pyopenms wheel."
            DEPENDS prepare_pyopenms_libs pyopenms_compile
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS )
    ENDIF()
ENDIF()

###########################################################################
#####                      Testing pyOpenMS                           #####
###########################################################################

# Pytest, testing all unittests at once
# => this is suboptimal for ctest and cdash because we don't see which tests
# actually have gone wrong. Thus we add additional tests below ...
enable_testing()
add_test(NAME test_pyopenms_unittests
         COMMAND ${Python_EXECUTABLE} -m pytest tests/unittests
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
  test_AcquisitionInfo.py
)

# Please add your test here when you decide to write a new testfile in the tests/integration_tests folder
set(pyopenms_integrationtest_testfiles
test_MRMRTNormalizer.py
)

# Loop through all the test files
foreach (t ${pyopenms_unittest_testfiles})
  add_test(NAME "pyopenms_unittest_${t}"
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_BINARY_DIR}/pyOpenMS/tests/unittests/${t}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
  if(NOT WIN32)
    set_tests_properties("pyopenms_unittest_${t}" PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
  endif()
endforeach(t)

foreach (t ${pyopenms_integrationtest_testfiles})
  add_test(NAME "pyopenms_integrationtest_${t}"
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_BINARY_DIR}/pyOpenMS/tests/integration_tests/${t}
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
  if(NOT WIN32)
    set_tests_properties("pyopenms_integrationtest_${t}" PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
  endif()
endforeach(t)

# Finally add the memory leaks test (in folder tests/memoryleaktests/)
if(NOT PY_MEMLEAK_DISABLE)
  add_test(NAME pyopenms_test_memoryleaktests
    COMMAND ${Python_EXECUTABLE} -m pytest ${CMAKE_BINARY_DIR}/pyOpenMS/tests/memoryleaktests/testAll.py
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}/pyOpenMS)
  if(NOT WIN32)
      set_tests_properties(pyopenms_test_memoryleaktests PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
  endif()
endif()
