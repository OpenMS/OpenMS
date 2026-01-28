# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------


cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("OpenMS_src")

#------------------------------------------------------------------------------
# OpenSWATHAlgo lib
#------------------------------------------------------------------------------
add_subdirectory(openswathalgo)

#------------------------------------------------------------------------------
# OpenMS core lib
#------------------------------------------------------------------------------
add_subdirectory(openms)

#------------------------------------------------------------------------------
# pyOpenMS
#------------------------------------------------------------------------------
# Note: pyOpenMS is optional and can be disabled with -DPYOPENMS=OFF
option(PYOPENMS "Build pyOpenMS Python bindings" OFF)
option(WITH_UV "Install Python dependencies with uv" ON)

if(PYOPENMS)
  if(WITH_UV)
    # Discover uv executable
    find_program(UV_EXECUTABLE uv)
    if(NOT UV_EXECUTABLE)
      message(FATAL_ERROR "uv executable not found. Install uv (https://docs.astral.sh/uv/) and ensure it's on PATH.")
    endif()
    message(STATUS "Found uv: ${UV_EXECUTABLE}")

    set(Python_FIND_VIRTUALENV FIRST)
    # Install only build dependencies (not pyOpenMS itself)
    message(STATUS "Installing pyOpenMS build dependencies with uv...")
    execute_process(
      COMMAND ${UV_EXECUTABLE} sync --only-group dev
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/pyOpenMS
      RESULT_VARIABLE _uv_install_result
      ERROR_VARIABLE _uv_install_error
    )
    if(NOT _uv_install_result EQUAL 0)
      message(WARNING "Failed to install pyOpenMS build dependencies with uv: ${_uv_install_error}")
    endif()
    if(WIN32)
      set(Python_EXECUTABLE "${CMAKE_CURRENT_LIST_DIR}/pyOpenMS/.venv/Scripts/python.exe")
    else()
      set(Python_EXECUTABLE "${CMAKE_CURRENT_LIST_DIR}/pyOpenMS/.venv/bin/python")
    endif()
  elseif(NOT DEFINED Python_EXECUTABLE)
    set(Python_EXECUTABLE "python")
  endif()

  message(STATUS "pyOpenMS Python environment: ${Python_EXECUTABLE}")

  # Build pyOpenMS Python extension modules directly as part of the main build
  # This ensures all CMake variables are naturally inherited
  add_subdirectory(pyOpenMS)

  # Optional: Package pyOpenMS wheel from already-built modules
  # Usage: cmake --build . --target pyopenms_wheel
  if(WITH_UV)
    add_custom_target(pyopenms_wheel
      COMMAND ${UV_EXECUTABLE} build --wheel
        --out-dir ${CMAKE_BINARY_DIR}/pyopenms_wheels
        --config-setting=override=cmake.0.build_path=${CMAKE_BINARY_DIR}/py-build-cmake_cache/pyopenms
        --config-setting=override=cmake.0.options.PYOPENMS_USE_PREBUILT=ON
        --config-setting=override=cmake.0.options.PYOPENMS_PREBUILT_DIR=${CMAKE_BINARY_DIR}/pyOpenMS
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/pyOpenMS
      DEPENDS pyopenms
      COMMENT "Packaging pyOpenMS wheel from existing build"
    )
  else()
    add_custom_target(pyopenms_wheel
      COMMAND ${Python_EXECUTABLE} -m build --wheel
        --skip-dependency-check --no-isolation
        --outdir ${CMAKE_BINARY_DIR}/pyopenms_wheels
        --config-setting=override=cmake.0.build_path=${CMAKE_BINARY_DIR}/py-build-cmake_cache/pyopenms
        --config-setting=override=cmake.0.options.PYOPENMS_USE_PREBUILT=ON
        --config-setting=override=cmake.0.options.PYOPENMS_PREBUILT_DIR=${CMAKE_BINARY_DIR}/pyOpenMS
      WORKING_DIRECTORY ${CMAKE_SOURCE_DIR}/src/pyOpenMS
      DEPENDS pyopenms
      COMMENT "Packaging pyOpenMS wheel from existing build"
    )
  endif()
endif()

if(WITH_GUI)
  #------------------------------------------------------------------------------
  # OpenMS gui lib
  #------------------------------------------------------------------------------
  add_subdirectory(openms_gui)
endif()

#------------------------------------------------------------------------------
# TOPP tools
#------------------------------------------------------------------------------
add_subdirectory(topp)

#------------------------------------------------------------------------------
# Tests
#------------------------------------------------------------------------------
add_subdirectory(tests)
