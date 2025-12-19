# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Hannes Röst $
# $Authors: Hannes Röst, Uwe Schmitt, Stephan Aiche $
# --------------------------------------------------------------------------
#
# pyOpenMS CMake build configuration
#
# This CMakeLists.txt can be used in two modes:
# 1. As part of the OpenMS build (when PYOPENMS=ON)
# 2. Standalone with py-build-cmake (pip install .)
#
# For standalone builds, set OpenMS_ROOT environment variable to point to
# an OpenMS build or install directory.
#
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("pyOpenMS" LANGUAGES C CXX)

#------------------------------------------------------------------------------
# Find OpenMS
#------------------------------------------------------------------------------
if(NOT TARGET OpenMS)
  # Standalone build - find installed/built OpenMS
  find_package(OpenMS REQUIRED)
  set(CF_OPENMS_PACKAGE_VERSION ${OpenMS_VERSION})

  # Try to get git info if available
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/Modules/GetGitRevisionDescription.cmake")
    include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/Modules/GetGitRevisionDescription.cmake)
    git_short_info(OPENMS_GIT_SHORT_REFSPEC OPENMS_GIT_SHORT_SHA1 OPENMS_GIT_LC_DATE)
  else()
    set(OPENMS_GIT_LC_DATE "unknown")
  endif()
endif()

#------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------
option(NO_DEPENDENCIES "Do not bundle dependencies (use wheel repair tools instead)" ON)
option(NO_SHARE "Do not bundle share data into wheel" OFF)
option(PY_NO_OPTIMIZATION "Disable optimization for faster compile time" OFF)
option(PY_NO_OUTPUT "Suppress output on Linux" OFF)
option(PY_MEMLEAK_DISABLE "Disable memory leak tests" OFF)

if(NOT PY_NUM_THREADS)
  set(PY_NUM_THREADS 4)
endif()
if(NOT PY_NUM_MODULES)
  set(PY_NUM_MODULES 8)
endif()

#------------------------------------------------------------------------------
# Platform detection
#------------------------------------------------------------------------------
if(UNIX AND NOT APPLE AND NOT CYGWIN)
    set(LINUX TRUE)
endif()

#------------------------------------------------------------------------------
# Find Python
#------------------------------------------------------------------------------
# py-build-cmake sets Python_EXECUTABLE and hints
find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)

message(STATUS "Python found at ${Python_EXECUTABLE} with version ${Python_VERSION}")
message(STATUS "Building pyOpenMS against OpenMS ${CF_OPENMS_PACKAGE_VERSION}")
message(STATUS "Number of submodules: ${PY_NUM_MODULES}")
message(STATUS "Number of concurrent threads: ${PY_NUM_THREADS}")

#------------------------------------------------------------------------------
# Windows compiler check
#------------------------------------------------------------------------------
if(WIN32)
    if(SKIP_WIN_COMPILERCHECK)
        message(STATUS "Skipping MS Visual C++ compiler check")
    elseif(MSVC AND (${Python_VERSION} VERSION_GREATER_EQUAL "3.5"))
        message(STATUS "Visual C++ 2015 or later for Python 3.5+ -- ok")
    else()
        message(STATUS "Need Visual C++ 2015 or later for building Python 3.5+ extensions")
        message(STATUS "To skip this check, configure with -DSKIP_WIN_COMPILERCHECK=On")
        message(FATAL_ERROR "Incompatible compiler version")
    endif()
endif()

#------------------------------------------------------------------------------
# Find Python package include directories
#------------------------------------------------------------------------------
# NumPy
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import numpy; print(numpy.get_include())"
    OUTPUT_VARIABLE NUMPY_INCLUDE_DIR
    ERROR_VARIABLE NUMPY_ERROR
    RESULT_VARIABLE NUMPY_RESULT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(NUMPY_RESULT)
    message(FATAL_ERROR "Failed to find numpy: ${NUMPY_ERROR}")
endif()
message(STATUS "NumPy include dir: ${NUMPY_INCLUDE_DIR}")

# Autowrap
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import autowrap; print(autowrap.__path__[0] + '/data_files/autowrap/')"
    OUTPUT_VARIABLE AUTOWRAP_INCLUDE_DIR
    ERROR_VARIABLE AUTOWRAP_ERROR
    RESULT_VARIABLE AUTOWRAP_RESULT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(AUTOWRAP_RESULT)
    message(FATAL_ERROR "Failed to find autowrap: ${AUTOWRAP_ERROR}")
endif()
message(STATUS "Autowrap include dir: ${AUTOWRAP_INCLUDE_DIR}")

# Python extension suffix
execute_process(
    COMMAND ${Python_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_config_var('EXT_SUFFIX'))"
    OUTPUT_VARIABLE PYTHON_EXT_SUFFIX
    RESULT_VARIABLE PYTHON_EXT_SUFFIX_RESULT
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
if(PYTHON_EXT_SUFFIX_RESULT)
    message(FATAL_ERROR "Could not determine Python extension suffix")
endif()

#------------------------------------------------------------------------------
# Find Qt6 (required for OpenMS)
#------------------------------------------------------------------------------
find_package(Qt6 COMPONENTS Core Network REQUIRED)
get_target_property(QT_QMAKE_EXECUTABLE Qt6::qmake IMPORTED_LOCATION)
execute_process(COMMAND ${QT_QMAKE_EXECUTABLE} -query QT_INSTALL_HEADERS OUTPUT_VARIABLE QT_INCLUDE_DIR OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${QT_QMAKE_EXECUTABLE} -query QT_INSTALL_LIBS OUTPUT_VARIABLE QT_INSTALL_LIBS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${QT_QMAKE_EXECUTABLE} -query QT_INSTALL_BINS OUTPUT_VARIABLE QT_INSTALL_BINS OUTPUT_STRIP_TRAILING_WHITESPACE)
execute_process(COMMAND ${QT_QMAKE_EXECUTABLE} -v OUTPUT_VARIABLE QT_QMAKE_VERSION_INFO OUTPUT_STRIP_TRAILING_WHITESPACE)

#------------------------------------------------------------------------------
# OpenMP support (optional)
#------------------------------------------------------------------------------
if(MT_ENABLE_OPENMP)
  find_package(OpenMP COMPONENTS CXX)
endif()

#------------------------------------------------------------------------------
# Collect include directories
#------------------------------------------------------------------------------
set(PYOPENMS_INCLUDE_DIRS
  ${OpenSwathAlgo_INCLUDE_DIRECTORIES}
  ${OpenMS_INCLUDE_DIRECTORIES}
  ${QT_INCLUDE_DIR}
  ${PROJECT_SOURCE_DIR}/extra_includes
)
list(REMOVE_DUPLICATES PYOPENMS_INCLUDE_DIRS)

if(OpenMP_FOUND)
  list(APPEND PYOPENMS_INCLUDE_DIRS ${OpenMP_CXX_INCLUDE_DIR})
endif()

#------------------------------------------------------------------------------
# Setup build directories
#------------------------------------------------------------------------------
set(PYOPENMS_BUILD_DIR "${CMAKE_BINARY_DIR}/pyOpenMS")
file(MAKE_DIRECTORY ${PYOPENMS_BUILD_DIR})
file(MAKE_DIRECTORY ${PYOPENMS_BUILD_DIR}/pyopenms)
file(MAKE_DIRECTORY ${PYOPENMS_BUILD_DIR}/extra_includes)

# Helper function to copy assets
function(_copy_assets _asset_dir _regex _target_dir)
  file(GLOB _assets "${_asset_dir}/${_regex}")
  foreach(_asset_file ${_assets})
    file(RELATIVE_PATH _relative_path ${_asset_dir} ${_asset_file})
    if(NOT IS_DIRECTORY ${_asset_file})
      configure_file(${_asset_file} ${_target_dir}/${_relative_path} COPYONLY)
    endif()
  endforeach()
endfunction()

# Copy Python source files
_copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "*.py" ${PYOPENMS_BUILD_DIR}/pyopenms/)
_copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "*.pyi" ${PYOPENMS_BUILD_DIR}/pyopenms/)
_copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "py.typed" ${PYOPENMS_BUILD_DIR}/pyopenms/)
_copy_assets("${PROJECT_SOURCE_DIR}/extra_includes/" "*" ${PYOPENMS_BUILD_DIR}/extra_includes/)

# Copy share data if requested and directory exists
# Note: When building standalone with py-build-cmake, the share dir may not exist
# In that case, share data should be included separately or installed system-wide
if(NOT NO_SHARE AND DEFINED OPENMS_SHARE_DIR AND EXISTS "${OPENMS_SHARE_DIR}")
  file(MAKE_DIRECTORY ${PYOPENMS_BUILD_DIR}/pyopenms/share)
  file(COPY ${OPENMS_SHARE_DIR} DESTINATION ${PYOPENMS_BUILD_DIR}/pyopenms/share)
endif()

#------------------------------------------------------------------------------
# Generate env.py for create_cpp_extension.py
# This is a simplified version - just what's needed for code generation
#------------------------------------------------------------------------------
set(OPEN_MS_VERSION "${CF_OPENMS_PACKAGE_VERSION}")
if(NOT OPENMS_GIT_LC_DATE)
  set(OPENMS_GIT_LC_DATE "unknown")
endif()
string(REPLACE "-" "" OPENMS_GIT_LC_DATE_REPLACED ${OPENMS_GIT_LC_DATE})
set(OPENMS_GIT_LC_DATE_REPLACED_LIST "${OPENMS_GIT_LC_DATE_REPLACED}")
separate_arguments(OPENMS_GIT_LC_DATE_REPLACED_LIST)
list(GET OPENMS_GIT_LC_DATE_REPLACED_LIST 0 OPENMS_GIT_LC_DATE_FORMAT)

file(WRITE ${PYOPENMS_BUILD_DIR}/env.py
"# Auto-generated by CMake for create_cpp_extension.py
PYOPENMS_SRC_DIR = '${PROJECT_SOURCE_DIR}'
OPEN_MS_BUILD_TYPE = '${CMAKE_BUILD_TYPE}'
OPEN_MS_VERSION = '${OPEN_MS_VERSION}'
PY_NUM_THREADS = '${PY_NUM_THREADS}'
PY_NUM_MODULES = '${PY_NUM_MODULES}'
QT_QMAKE_VERSION_INFO = '''${QT_QMAKE_VERSION_INFO}'''
# Legacy variables (not used but kept for compatibility)
OPEN_MS_CONTRIB_BUILD_DIRS = ''
OPEN_MS_LIB = 'OpenMS'
OPEN_SWATH_ALGO_LIB = 'OpenSwathAlgo'
OPEN_MS_BUILD_DIR = '${CMAKE_BINARY_DIR}'
Boost_MAJOR_VERSION = ''
Boost_MINOR_VERSION = ''
MSVS_RTLIBS = ''
")

# Copy create_cpp_extension.py
configure_file(${PROJECT_SOURCE_DIR}/create_cpp_extension.py ${PYOPENMS_BUILD_DIR}/create_cpp_extension.py COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/fix_pyi_imports.py ${PYOPENMS_BUILD_DIR}/fix_pyi_imports.py COPYONLY)

#------------------------------------------------------------------------------
# Define generated files
#------------------------------------------------------------------------------
if("${PY_NUM_MODULES}" GREATER "1")
    set(PYCPPFILES "")
    set(PYPYXFILES "")
    foreach(NUM RANGE 1 ${PY_NUM_MODULES})
      list(APPEND PYCPPFILES "${PYOPENMS_BUILD_DIR}/pyopenms/_pyopenms_${NUM}.cpp")
      list(APPEND PYPYXFILES "${PYOPENMS_BUILD_DIR}/pyopenms/_pyopenms_${NUM}.pyx")
    endforeach()
else()
    set(PYCPPFILES "${PYOPENMS_BUILD_DIR}/pyopenms/_pyopenms.cpp")
    set(PYPYXFILES "${PYOPENMS_BUILD_DIR}/pyopenms/_pyopenms.pyx")
endif()

#------------------------------------------------------------------------------
# Code generation with autowrap and Cython
#------------------------------------------------------------------------------
set(PYGEN_SENTINEL "${PYOPENMS_BUILD_DIR}/.cpp_extension_generated")

if(LINUX AND PY_NO_OUTPUT)
  add_custom_command(
    OUTPUT ${PYGEN_SENTINEL}
    DEPENDS OpenMS
    COMMAND ${Python_EXECUTABLE} create_cpp_extension.py 2> /dev/null
    COMMAND ${CMAKE_COMMAND} -E touch ${PYGEN_SENTINEL}
    COMMENT "Creating C++ extension with autowrap and Cython"
    WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
  )
else()
  add_custom_command(
    OUTPUT ${PYGEN_SENTINEL}
    DEPENDS OpenMS
    COMMAND ${Python_EXECUTABLE} create_cpp_extension.py
    COMMAND ${CMAKE_COMMAND} -E touch ${PYGEN_SENTINEL}
    COMMENT "Creating C++ extension with autowrap and Cython"
    WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
  )
endif()

add_custom_command(
  OUTPUT ${PYCPPFILES} ${PYPYXFILES}
  DEPENDS ${PYGEN_SENTINEL}
  COMMENT "Generated files depend on successful code generation"
)

add_custom_target(compile_pxds DEPENDS ${PYCPPFILES} ${PYPYXFILES})

#------------------------------------------------------------------------------
# Compiler settings
#------------------------------------------------------------------------------
set(PY_RPATH ${CMAKE_INSTALL_RPATH})
if(APPLE)
  list(APPEND PY_RPATH "@loader_path/")
endif()

# Visibility settings
if(LINUX OR APPLE)
  set(VIS_PRESET "default")
  set(INLINES_HIDDEN 0)
else()
  set(VIS_PRESET "hidden")
  set(INLINES_HIDDEN 1)
endif()

# Warning suppressions
if(LINUX OR APPLE)
  set(EXTRA_WARNINGS "-Wno-builtin-macro-redefined" "-Wno-unused" "-Wno-zero-as-null-pointer-constant")
  if("${CMAKE_CXX_COMPILER_ID}" MATCHES "Clang")
    list(APPEND EXTRA_WARNINGS "-Wno-unused-member-function")
  endif()
else()
  set(EXTRA_WARNINGS "")
endif()

#------------------------------------------------------------------------------
# Build extension modules
#------------------------------------------------------------------------------
if("${PY_NUM_MODULES}" GREATER "1")
    add_custom_target(pyopenms_compile ALL)
    foreach(nr RANGE 1 ${PY_NUM_MODULES})
      Python_add_library(pyopenms_${nr} MODULE "${PYOPENMS_BUILD_DIR}/pyopenms/_pyopenms_${nr}.cpp")
      target_link_libraries(pyopenms_${nr} PRIVATE OpenMS)
      target_compile_definitions(pyopenms_${nr} PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
      target_include_directories(pyopenms_${nr} SYSTEM PUBLIC
        ${PYOPENMS_INCLUDE_DIRS}
        ${NUMPY_INCLUDE_DIR}
        ${AUTOWRAP_INCLUDE_DIR}
      )
      target_compile_options(pyopenms_${nr} PRIVATE ${EXTRA_WARNINGS})
      set_target_properties(pyopenms_${nr} PROPERTIES
        OUTPUT_NAME _pyopenms_${nr}
        PREFIX ""
        SUFFIX ${PYTHON_EXT_SUFFIX}
        CXX_VISIBILITY_PRESET ${VIS_PRESET}
        VISIBILITY_INLINES_HIDDEN ${INLINES_HIDDEN}
        INSTALL_RPATH "${PY_RPATH}"
        BUILD_WITH_INSTALL_RPATH 1
        BUILD_WITH_INSTALL_NAME_DIR 1
        LIBRARY_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms"
        RUNTIME_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms"
      )
      add_dependencies(pyopenms_${nr} compile_pxds)
      add_dependencies(pyopenms_compile pyopenms_${nr})

      # Fix .pyi imports
      add_custom_command(TARGET pyopenms_${nr} POST_BUILD
        COMMAND ${Python_EXECUTABLE} ${PYOPENMS_BUILD_DIR}/fix_pyi_imports.py ${PYOPENMS_BUILD_DIR}/pyopenms/_pyopenms_${nr}.pyi
        COMMENT "Fixing typing information for _pyopenms_${nr}"
      )

      # Install for py-build-cmake
      install(TARGETS pyopenms_${nr}
        COMPONENT python_modules
        LIBRARY DESTINATION pyopenms
        RUNTIME DESTINATION pyopenms
      )
    endforeach()
else()
    Python_add_library(pyopenms_compile MODULE "${PYOPENMS_BUILD_DIR}/pyopenms/_pyopenms.cpp")
    target_link_libraries(pyopenms_compile PRIVATE OpenMS)
    target_compile_definitions(pyopenms_compile PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
    target_include_directories(pyopenms_compile SYSTEM PUBLIC
      ${PYOPENMS_INCLUDE_DIRS}
      ${NUMPY_INCLUDE_DIR}
      ${AUTOWRAP_INCLUDE_DIR}
    )
    target_compile_options(pyopenms_compile PRIVATE ${EXTRA_WARNINGS})
    set_target_properties(pyopenms_compile PROPERTIES
      OUTPUT_NAME _pyopenms
      PREFIX ""
      SUFFIX ${PYTHON_EXT_SUFFIX}
      CXX_VISIBILITY_PRESET ${VIS_PRESET}
      VISIBILITY_INLINES_HIDDEN ${INLINES_HIDDEN}
      INSTALL_RPATH "${PY_RPATH}"
      BUILD_WITH_INSTALL_RPATH 1
      BUILD_WITH_INSTALL_NAME_DIR 1
      LIBRARY_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms"
      RUNTIME_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms"
    )
    add_dependencies(pyopenms_compile compile_pxds)

    add_custom_command(TARGET pyopenms_compile POST_BUILD
      COMMAND ${Python_EXECUTABLE} ${PYOPENMS_BUILD_DIR}/fix_pyi_imports.py ${PYOPENMS_BUILD_DIR}/pyopenms/pyopenms.pyi
      COMMENT "Fixing typing information for pyopenms"
    )

    install(TARGETS pyopenms_compile
      COMPONENT python_modules
      LIBRARY DESTINATION pyopenms
      RUNTIME DESTINATION pyopenms
    )
endif()

# Create pyopenms alias target for CI compatibility (cibuild.cmake expects this target)
add_custom_target(pyopenms DEPENDS pyopenms_compile)

#------------------------------------------------------------------------------
# Install Python files for py-build-cmake
#------------------------------------------------------------------------------
install(DIRECTORY ${PYOPENMS_BUILD_DIR}/pyopenms/
  DESTINATION pyopenms
  COMPONENT python_modules
  FILES_MATCHING
  PATTERN "*.py"
  PATTERN "*.pyi"
  PATTERN "py.typed"
  PATTERN "__pycache__" EXCLUDE
)

# Install share data if included
if(NOT NO_SHARE AND EXISTS "${PYOPENMS_BUILD_DIR}/pyopenms/share")
  install(DIRECTORY ${PYOPENMS_BUILD_DIR}/pyopenms/share
    DESTINATION pyopenms
    COMPONENT python_modules
  )
endif()

#------------------------------------------------------------------------------
# Legacy pyopenms target (for backwards compatibility with old build system)
# This target is NOT used by py-build-cmake
#------------------------------------------------------------------------------
if(NOT NO_DEPENDENCIES)
    add_custom_target(prepare_pyopenms_libs DEPENDS OpenMS)
    set(PYOPENMS_DEPENDENCIES OpenMS OpenSwathAlgo)

    foreach(PYOPENMS_DEPENDENCY ${PYOPENMS_DEPENDENCIES})
      add_custom_command(TARGET prepare_pyopenms_libs POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E copy $<TARGET_FILE:${PYOPENMS_DEPENDENCY}> ${PYOPENMS_BUILD_DIR}/pyopenms
      )
    endforeach()

    if(NOT APPLE AND NOT WIN32)
      foreach(PYOPENMS_DEPENDENCY ${PYOPENMS_DEPENDENCIES})
        add_custom_command(TARGET prepare_pyopenms_libs POST_BUILD
          COMMAND ${CMAKE_COMMAND}
            -DDEPS="$<TARGET_FILE:${PYOPENMS_DEPENDENCY}>"
            -DTARGET="${PYOPENMS_BUILD_DIR}/pyopenms/"
            -DLOOKUP_DIRS="${OPENMS_CONTRIB_LIBS}/lib\;${OPENMS_CONTRIB_LIBS}/lib64\;${QT_INSTALL_BINS}\;${QT_INSTALL_LIBS}"
            -P ${PROJECT_SOURCE_DIR}/pyopenms_copy_deps.cmake
        )
      endforeach()
    endif()
endif()

# Legacy wheel packaging target (not used by py-build-cmake)
add_custom_target(pyopenms_wheel
  COMMAND ${Python_EXECUTABLE} -m pip wheel -w dist --no-deps --no-build-isolation .
  COMMENT "Packaging pyopenms wheel (legacy target)"
  DEPENDS pyopenms_compile
  WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
)

#------------------------------------------------------------------------------
# Testing
#------------------------------------------------------------------------------
enable_testing()

# Copy test files
file(MAKE_DIRECTORY ${PYOPENMS_BUILD_DIR}/tests/unittests)
file(MAKE_DIRECTORY ${PYOPENMS_BUILD_DIR}/tests/memoryleaktests)
file(MAKE_DIRECTORY ${PYOPENMS_BUILD_DIR}/tests/integration_tests)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/unittests/" "*" ${PYOPENMS_BUILD_DIR}/tests/unittests)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/" "*.mzXML" ${PYOPENMS_BUILD_DIR}/tests)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/memoryleaktests/" "*" ${PYOPENMS_BUILD_DIR}/tests/memoryleaktests)
_copy_assets("${PROJECT_SOURCE_DIR}/tests/integration_tests/" "*" ${PYOPENMS_BUILD_DIR}/tests/integration_tests)

# All unittests
add_test(NAME test_pyopenms_unittests
  COMMAND ${Python_EXECUTABLE} -m pytest tests/unittests
  WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
)
if(NOT WIN32)
  set_tests_properties(test_pyopenms_unittests PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
endif()

# Individual unittest files
set(pyopenms_unittest_testfiles
  test_ProteaseDigestion.py
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
  testLocale.py
  test_MRMFeatureFinderScoring.py
  test_MSNumpressCoder.py
  test_MSSpectrum.py
  test_OpenSwathDataStructures.py
  test_Smoothing.py
  testSpecialCases.py
  test_SpectraFilter.py
  test_SpectrumAccessOpenMS.py
  test_TraML.py
  test_MzMLConsumer.py
  test_MzXMLConsumer.py
  test_AcquisitionInfo.py
  test_FLASHDeconv.py
)

set(pyopenms_integrationtest_testfiles
  test_MRMRTNormalizer.py
)

foreach(t ${pyopenms_unittest_testfiles})
  add_test(NAME "pyopenms_unittest_${t}"
    COMMAND ${Python_EXECUTABLE} -m pytest ${PYOPENMS_BUILD_DIR}/tests/unittests/${t}
    WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
  )
  if(NOT WIN32)
    set_tests_properties("pyopenms_unittest_${t}" PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
  endif()
endforeach()

foreach(t ${pyopenms_integrationtest_testfiles})
  add_test(NAME "pyopenms_integrationtest_${t}"
    COMMAND ${Python_EXECUTABLE} -m pytest ${PYOPENMS_BUILD_DIR}/tests/integration_tests/${t}
    WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
  )
  if(NOT WIN32)
    set_tests_properties("pyopenms_integrationtest_${t}" PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
  endif()
endforeach()

if(NOT PY_MEMLEAK_DISABLE)
  add_test(NAME pyopenms_test_memoryleaktests
    COMMAND ${Python_EXECUTABLE} -m pytest ${PYOPENMS_BUILD_DIR}/tests/memoryleaktests/testAll.py
    WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
  )
  if(NOT WIN32)
    set_tests_properties(pyopenms_test_memoryleaktests PROPERTIES ENVIRONMENT "LD_LIBRARY_PATH=${CMAKE_BINARY_DIR}/lib")
  endif()
endif()
