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

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("pyOpenMS" LANGUAGES CXX)

#------------------------------------------------------------------------------
# Find OpenMS
#------------------------------------------------------------------------------
if(NOT TARGET OpenMS)
  # Standalone build - find installed/built OpenMS
  find_package(OpenMS REQUIRED)
  # TODO the following is just a fallback.
  # CF_OPENMS_PACKAGE_VERSION is what is written into pyopenms/_version.py
  # which is the old way of getting the version.
  # We will now switch to reading the version from the pyopenms pyproject.toml
  # directly.
  # The problem is that we might not always want the exact same version as OpenMS and
  # previously we did some logic in CMakeLists.txt or setup.py to add some kind of nightly
  # dev specifier to the end of it. We need to think about a better versioning strategy here.
  set(CF_OPENMS_PACKAGE_VERSION ${OpenMS_VERSION})

  # Try to get git info if available
  if(EXISTS "${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/Modules/GetGitRevisionDescription.cmake")
    include(${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/Modules/GetGitRevisionDescription.cmake)
    git_short_info(OPENMS_GIT_SHORT_REFSPEC OPENMS_GIT_SHORT_SHA1 OPENMS_GIT_LC_DATE)
  else()
    set(OPENMS_GIT_LC_DATE "unknown")
  endif()
else()
  find_package(Qt6 COMPONENTS Core Network REQUIRED)
endif()

#------------------------------------------------------------------------------
# Options
#------------------------------------------------------------------------------
option(NO_DEPENDENCIES "Do not bundle dependencies (use wheel repair tools instead)" ON)
option(NO_SHARE "Do not bundle share data into wheel" OFF)
option(PY_NO_OUTPUT "Suppress output on Linux" OFF)

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
find_package(Python COMPONENTS Interpreter Development.Module NumPy REQUIRED)

message(STATUS "Python found at ${Python_EXECUTABLE} with version ${Python_VERSION}")
message(STATUS "Building pyOpenMS against OpenMS ${CF_OPENMS_PACKAGE_VERSION}")
message(STATUS "Number of submodules: ${PY_NUM_MODULES}")
message(STATUS "Number of concurrent threads: ${PY_NUM_THREADS}")

#------------------------------------------------------------------------------
# Check share directory
# OPENMS_SHARE_DIR should be set by find_package(OpenMS) from OpenMSConfig.cmake
# If not set, fall back to relative path from source tree (for cibuildwheel builds)
#------------------------------------------------------------------------------
if(DEFINED OPENMS_SHARE_DIR AND EXISTS "${OPENMS_SHARE_DIR}")
  message(STATUS "OpenMS share directory from config: ${OPENMS_SHARE_DIR}")
else()
  # Try relative path from pyOpenMS source directory (../../share/OpenMS)
  get_filename_component(_POTENTIAL_SHARE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../../share/OpenMS" ABSOLUTE)
  if(EXISTS "${_POTENTIAL_SHARE_DIR}")
    set(OPENMS_SHARE_DIR "${_POTENTIAL_SHARE_DIR}")
    message(STATUS "Using relative OpenMS share directory: ${OPENMS_SHARE_DIR}")
  else()
    unset(OPENMS_SHARE_DIR)
    message(WARNING "OPENMS_SHARE_DIR not found - share data will not be included in wheel")
  endif()
endif()

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

# Use Qt6 version from the find module instead of querying qmake
set(QT_QMAKE_VERSION_INFO "${Qt6_VERSION}")

# Get Autowrap include dirs
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
  file(COPY ${OPENMS_SHARE_DIR} DESTINATION ${PYOPENMS_BUILD_DIR}/pyopenms/share
       PATTERN "examples" EXCLUDE)
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
")

# Copy create_cpp_extension.py and pyproject.toml
configure_file(${PROJECT_SOURCE_DIR}/create_cpp_extension.py ${PYOPENMS_BUILD_DIR}/create_cpp_extension.py COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/fix_pyi_imports.py ${PYOPENMS_BUILD_DIR}/fix_pyi_imports.py COPYONLY)
configure_file(${PROJECT_SOURCE_DIR}/pyproject.toml ${PYOPENMS_BUILD_DIR}/pyproject.toml COPYONLY)

#------------------------------------------------------------------------------
# Define generated files
#------------------------------------------------------------------------------
# Files generated by create_cpp_extension.py
set(PYGEN_EXTRA_FILES
    "${PYOPENMS_BUILD_DIR}/pyopenms/_version.py"
    "${PYOPENMS_BUILD_DIR}/pyopenms/_qt_version_info.py"
    "${PYOPENMS_BUILD_DIR}/pyopenms/_all_modules.py"
    "${PYOPENMS_BUILD_DIR}/pyopenms/_all_modules.pyi"
)

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
    # Fake dependency on OpenMS to trigger rebuild when headers change
    # TODO: Should depend on .pxd files instead for more accurate rebuilds
    DEPENDS OpenMS
    COMMAND ${Python_EXECUTABLE} create_cpp_extension.py 2> /dev/null
    COMMAND ${CMAKE_COMMAND} -E touch ${PYGEN_SENTINEL}
    COMMENT "Creating C++ extension with autowrap and Cython"
    WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
  )
else()
  add_custom_command(
    OUTPUT ${PYGEN_SENTINEL}
    # Fake dependency on OpenMS to trigger rebuild when headers change
    # TODO: Should depend on .pxd files instead for more accurate rebuilds
    DEPENDS OpenMS
    COMMAND ${Python_EXECUTABLE} create_cpp_extension.py
    COMMAND ${CMAKE_COMMAND} -E touch ${PYGEN_SENTINEL}
    COMMENT "Creating C++ extension with autowrap and Cython"
    WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
  )
endif()

add_custom_command(
  OUTPUT ${PYCPPFILES} ${PYPYXFILES} ${PYGEN_EXTRA_FILES}
  DEPENDS ${PYGEN_SENTINEL}
  COMMENT "Generated files depend on successful code generation"
)

add_custom_target(compile_pxds DEPENDS ${PYCPPFILES} ${PYPYXFILES})

#------------------------------------------------------------------------------
# Compiler settings
#------------------------------------------------------------------------------

# Option to prepare modules for wheel repair tools (delocate, auditwheel, delvewheel)
# Uses INSTALL_RPATH_USE_LINK_PATH to automatically embed library paths that repair tools can follow
option(PYOPENMS_PREPARE_WHEEL_REPAIR "Prepare modules for wheel repair tool compatibility" OFF)
# TODO remove once we disable the old way of building wheels and just use CMAKE_INSTALL_RPATH_USE_LINK_PATH

if(NOT PYOPENMS_PREPARE_WHEEL_REPAIR)
  set(PY_RPATH ${CMAKE_INSTALL_RPATH})
  if(APPLE)
    # TODO: Maybe we need a similar rpath for Linux too (e.g., $ORIGIN)
    list(APPEND PY_RPATH "@loader_path/")
  endif()
  message(STATUS "pyOpenMS: PY_RPATH = '${PY_RPATH}'")
  message(STATUS "pyOpenMS: CMAKE_INSTALL_RPATH was = '${CMAKE_INSTALL_RPATH}'")
else()
  # When building wheels, use INSTALL_RPATH_USE_LINK_PATH to automatically embed
  # absolute paths to linked libraries (OpenMS), which repair tools will follow
  set(PY_RPATH "")
  message(STATUS "Preparing modules for wheel repair: INSTALL_RPATH_USE_LINK_PATH enabled")
endif()

# Visibility settings
if(LINUX OR APPLE)
  # if the following vars are not set from outside, set them:
  if(NOT DEFINED CMAKE_CXX_VISIBILITY_PRESET)
    set(CMAKE_CXX_VISIBILITY_PRESET "default")
  endif()
  if(NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN)
    set(CMAKE_VISIBILITY_INLINES_HIDDEN 0)
  endif()
else()
  if(NOT DEFINED CMAKE_CXX_VISIBILITY_PRESET)
    set(CMAKE_CXX_VISIBILITY_PRESET "hidden")
  endif()
  if(NOT DEFINED CMAKE_VISIBILITY_INLINES_HIDDEN)
    set(CMAKE_VISIBILITY_INLINES_HIDDEN 1)
  endif()
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
      Python_add_library(pyopenms_${nr} MODULE "${PYOPENMS_BUILD_DIR}/pyopenms/_pyopenms_${nr}.cpp" WITH_SOABI)
      target_link_libraries(pyopenms_${nr} PRIVATE OpenMS Python::NumPy)
      target_compile_definitions(pyopenms_${nr} PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
      target_include_directories(pyopenms_${nr} SYSTEM PRIVATE
        ${AUTOWRAP_INCLUDE_DIR}
        ${PYOPENMS_BUILD_DIR}/extra_includes
      )
      target_compile_options(pyopenms_${nr} PRIVATE ${EXTRA_WARNINGS})
      set_target_properties(pyopenms_${nr} PROPERTIES
        OUTPUT_NAME _pyopenms_${nr}
        LIBRARY_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms"
        RUNTIME_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms"
      )
      if(NOT PYOPENMS_PREPARE_WHEEL_REPAIR)
        set_target_properties(pyopenms_${nr} PROPERTIES
          INSTALL_RPATH "${PY_RPATH}"
          BUILD_WITH_INSTALL_RPATH 1
          BUILD_WITH_INSTALL_NAME_DIR 1
        )
      else()
        # For wheel repair: use INSTALL_RPATH_USE_LINK_PATH to auto-embed library paths
        # Repair tools will follow these paths to bundle libs, then rewrite to relative paths
        set_target_properties(pyopenms_${nr} PROPERTIES
          INSTALL_RPATH_USE_LINK_PATH TRUE
          INSTALL_NAME_DIR ""
        )
      endif()
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
    target_link_libraries(pyopenms_compile PRIVATE OpenMS Python::NumPy)
    target_compile_definitions(pyopenms_compile PRIVATE NPY_NO_DEPRECATED_API=NPY_1_7_API_VERSION)
    target_include_directories(pyopenms_compile SYSTEM PUBLIC
      ${AUTOWRAP_INCLUDE_DIR}
      ${PYOPENMS_BUILD_DIR}/extra_includes
    )
    target_compile_options(pyopenms_compile PRIVATE ${EXTRA_WARNINGS})
    set_target_properties(pyopenms_compile PROPERTIES
      OUTPUT_NAME _pyopenms
      LIBRARY_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms"
      RUNTIME_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms"
    )
    if(NOT PYOPENMS_PREPARE_WHEEL_REPAIR)
      set_target_properties(pyopenms_compile PROPERTIES
        INSTALL_RPATH "${PY_RPATH}"
        BUILD_WITH_INSTALL_RPATH 1
        BUILD_WITH_INSTALL_NAME_DIR 1
      )
    else()
      # For wheel repair: use INSTALL_RPATH_USE_LINK_PATH to auto-embed library paths
      # Repair tools will follow these paths to bundle libs, then rewrite to relative paths
      set_target_properties(pyopenms_compile PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH TRUE
        INSTALL_NAME_DIR ""
      )
    endif()
    add_dependencies(pyopenms_compile compile_pxds)

    # Fix .pyi imports (single-module case uses _pyopenms.pyi)
    add_custom_command(TARGET pyopenms_compile POST_BUILD
      COMMAND ${Python_EXECUTABLE} ${PYOPENMS_BUILD_DIR}/fix_pyi_imports.py ${PYOPENMS_BUILD_DIR}/pyopenms/_pyopenms.pyi
      COMMENT "Fixing typing information for _pyopenms"
    )

    install(TARGETS pyopenms_compile
      COMPONENT python_modules
      LIBRARY DESTINATION pyopenms
      RUNTIME DESTINATION pyopenms
    )
endif()

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
# Define dependencies list (used by legacy and py-build-cmake paths)
#------------------------------------------------------------------------------
if(NOT NO_DEPENDENCIES)
  set(PYOPENMS_DEPENDENCIES OpenMS OpenSwathAlgo)
endif()

# Install OpenMS libraries to staging directory for py-build-cmake (before fix_dependencies.rb runs)
if(APPLE AND NOT NO_DEPENDENCIES)
  foreach(PYOPENMS_DEPENDENCY ${PYOPENMS_DEPENDENCIES})
    install(FILES $<TARGET_FILE:${PYOPENMS_DEPENDENCY}>
      DESTINATION pyopenms
      COMPONENT python_modules
    )
  endforeach()
endif()

# On macOS with NO_DEPENDENCIES=OFF, fix install_names in installed files
if(APPLE AND NOT NO_DEPENDENCIES)
  install(CODE "
    execute_process(
      COMMAND ${CMAKE_CURRENT_SOURCE_DIR}/../../cmake/MacOSX/fix_dependencies.rb -l \"\${CMAKE_INSTALL_PREFIX}/pyopenms\" -e \"@rpath/\" -f -v
      RESULT_VARIABLE _fix_result
    )
    if(_fix_result)
      message(WARNING \"fix_dependencies.rb failed with code \${_fix_result}\")
    endif()
    # Codesign the libraries
    execute_process(
      COMMAND codesign --force --sign \"-\" \"\${CMAKE_INSTALL_PREFIX}/pyopenms/Qt*\"
      ERROR_QUIET
    )
    execute_process(
      COMMAND /bin/sh -c \"codesign --force --sign '-' \\\"\${CMAKE_INSTALL_PREFIX}/pyopenms\\\"/*.dylib\"
      RESULT_VARIABLE _codesign_result
    )
    if(_codesign_result)
      message(STATUS \"Codesigning completed (some files may have been skipped)\")
    endif()
  "
  COMPONENT python_modules
  )
endif()