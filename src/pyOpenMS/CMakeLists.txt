# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# pyOpenMS CMake build configuration (nanobind)
#
# This CMakeLists.txt can be used in two modes:
# 1. As part of the OpenMS build (when PYOPENMS=ON)
# 2. Standalone with py-build-cmake (pip install .)
#
# Bindings are compiled from pre-committed generated sources in
# bindings/generated/*.cpp — no code generation at build time.
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.20 FATAL_ERROR)
project("pyOpenMS" LANGUAGES CXX)

# C++17 is required for nanobind
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# --------------------------------------------------------------------------
# Options
# --------------------------------------------------------------------------
# Fast-path: reuse prebuilt pyOpenMS artifacts and skip rebuild/configure work
option(PYOPENMS_USE_PREBUILT "Reuse existing pyOpenMS build outputs" OFF)
set(PYOPENMS_PREBUILT_DIR "" CACHE PATH "Path to prebuilt pyOpenMS build directory (containing pyopenms folder)")

option(NO_DEPENDENCIES "Do not bundle dependencies (use wheel repair tools instead)" ON)
option(NO_SHARE "Do not bundle share data into wheel" OFF)
option(PYOPENMS_GENERATE_STUBS "Generate .pyi stub files" ON)
option(PYOPENMS_PREPARE_WHEEL_REPAIR "Set RPATH for wheel repair tools (auditwheel/delocate)" OFF)

set(PYOPENMS_DOMAINS
    kernel metadata chemistry analysis featurefinder
    format processing datastructures ml misc
)

# Set default build directory (may be overridden by PYOPENMS_USE_PREBUILT)
if(NOT DEFINED PYOPENMS_BUILD_DIR OR "${PYOPENMS_BUILD_DIR}" STREQUAL "")
  set(PYOPENMS_BUILD_DIR "${CMAKE_BINARY_DIR}/pyOpenMS")
endif()

if(PYOPENMS_USE_PREBUILT)
  if(NOT PYOPENMS_PREBUILT_DIR)
    message(FATAL_ERROR "PYOPENMS_PREBUILT_DIR must be set when PYOPENMS_USE_PREBUILT=ON")
  endif()
  # Override build dir to use prebuilt location
  set(PYOPENMS_BUILD_DIR "${PYOPENMS_PREBUILT_DIR}")

  # Provide a stub target so dependent targets remain satisfied
  add_custom_target(pyopenms_compile ALL)
else()

  # --------------------------------------------------------------------------
  # Platform detection
  # --------------------------------------------------------------------------
  if(UNIX AND NOT APPLE AND NOT CYGWIN)
    set(LINUX TRUE)
  endif()

  # --------------------------------------------------------------------------
  # Find dependencies
  # --------------------------------------------------------------------------
  if(NOT TARGET OpenMS)
    # Standalone build - find installed/built OpenMS
    find_package(OpenMS REQUIRED)
    set(CF_OPENMS_PACKAGE_VERSION ${OpenMS_VERSION})
  endif()

  # Handle both namespaced and non-namespaced target names
  if(TARGET OpenMS::OpenMS)
    set(OPENMS_TARGET OpenMS::OpenMS)
    set(OPENSWATHALGO_TARGET OpenMS::OpenSwathAlgo)
  elseif(TARGET OpenMS)
    set(OPENMS_TARGET OpenMS)
    set(OPENSWATHALGO_TARGET OpenSwathAlgo)
  else()
    message(FATAL_ERROR "OpenMS target not found")
  endif()

  # --------------------------------------------------------------------------
  # Find Python
  # --------------------------------------------------------------------------
  # py-build-cmake sets Python_EXECUTABLE and hints
  # If requested, run `uv sync` first in the build dir so that uv creates
  # a virtualenv there. This must happen before calling find_package(Python)
  # so that CMake can pick up the created interpreter.
  if(WITH_UV AND DEFINED UV_EXECUTABLE AND NOT DEFINED PY_BUILD_CMAKE_VERSION)
    message(STATUS "Found uv executable: ${UV_EXECUTABLE} — ensuring .venv in ${PYOPENMS_BUILD_DIR}")
    file(MAKE_DIRECTORY "${PYOPENMS_BUILD_DIR}")
    # Copy pyproject.toml before uv sync so that uv can read the project configuration
    configure_file(${PROJECT_SOURCE_DIR}/pyproject.toml ${PYOPENMS_BUILD_DIR}/pyproject.toml COPYONLY)
    if(DEFINED PYOPENMS_UV_PYTHON_VERSION AND NOT "${PYOPENMS_UV_PYTHON_VERSION}" STREQUAL "")
      message(STATUS "Invoking uv sync with Python version: ${PYOPENMS_UV_PYTHON_VERSION}")
      execute_process(
        COMMAND ${UV_EXECUTABLE} sync --only-group dev --python ${PYOPENMS_UV_PYTHON_VERSION}
        WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
        RESULT_VARIABLE _uv_install_result
        ERROR_VARIABLE _uv_install_error
      )
    else()
      execute_process(
        COMMAND ${UV_EXECUTABLE} sync --only-group dev
        WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
        RESULT_VARIABLE _uv_install_result
        ERROR_VARIABLE _uv_install_error
      )
    endif()
    if(NOT _uv_install_result EQUAL 0)
      message(WARNING "uv sync failed: ${_uv_install_error}")
    endif()

    if(WIN32)
      set(_py_venv_python "${PYOPENMS_BUILD_DIR}/.venv/Scripts/python.exe")
    else()
      set(_py_venv_python "${PYOPENMS_BUILD_DIR}/.venv/bin/python")
    endif()
    if(EXISTS "${_py_venv_python}")
      set(Python_EXECUTABLE "${_py_venv_python}")
      message(STATUS "Using venv Python at ${Python_EXECUTABLE}")
    else()
      message(STATUS "No venv Python found at ${_py_venv_python}; falling back to system Python discovery")
    endif()
  endif()

  find_package(Python 3.9 REQUIRED COMPONENTS Interpreter Development.Module NumPy)

  # Find or fetch nanobind
  if(NOT TARGET nanobind::nanobind)
    # Try find_package first (installed via pip or system package)
    find_package(nanobind CONFIG QUIET)
    if(NOT nanobind_FOUND)
      include(FetchContent)
      FetchContent_Declare(
          nanobind
          GIT_REPOSITORY https://github.com/wjakob/nanobind.git
          GIT_TAG        v2.4.0
      )
      FetchContent_MakeAvailable(nanobind)
    endif()
  endif()

  message(STATUS "Python found at ${Python_EXECUTABLE} with version ${Python_VERSION}")
  message(STATUS "Building pyOpenMS (nanobind) against OpenMS ${CF_OPENMS_PACKAGE_VERSION}")
  if(NO_DEPENDENCIES)
    message(STATUS "pyOpenMS will be built WITHOUT bundling dependencies (NO_DEPENDENCIES=ON)")
  else()
    message(STATUS "pyOpenMS will be built WITH dependencies bundled (NO_DEPENDENCIES=OFF)")
  endif()

  # --------------------------------------------------------------------------
  # Check share directory
  # --------------------------------------------------------------------------
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

  # --------------------------------------------------------------------------
  # RPATH settings
  # --------------------------------------------------------------------------
  if(NOT PYOPENMS_PREPARE_WHEEL_REPAIR)
    set(PY_RPATH ${CMAKE_INSTALL_RPATH})
    if(APPLE)
      list(APPEND PY_RPATH "@loader_path/")
    elseif(LINUX)
      list(APPEND PY_RPATH "$ORIGIN/")
    endif()
    message(STATUS "pyOpenMS: PY_RPATH = '${PY_RPATH}'")
  else()
    set(PY_RPATH "")
    message(STATUS "Preparing modules for wheel repair: INSTALL_RPATH_USE_LINK_PATH enabled")
  endif()

  # --------------------------------------------------------------------------
  # Setup build directories
  # --------------------------------------------------------------------------
  file(MAKE_DIRECTORY ${PYOPENMS_BUILD_DIR})
  file(MAKE_DIRECTORY ${PYOPENMS_BUILD_DIR}/pyopenms)
  file(MAKE_DIRECTORY ${PYOPENMS_BUILD_DIR}/pyopenms/addons)

  # --------------------------------------------------------------------------
  # Pre-committed generated sources (no generator execution at build time)
  # --------------------------------------------------------------------------
  set(PYOPENMS_BINDINGS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bindings/generated")

  set(PYOPENMS_GENERATED_SOURCES "${PYOPENMS_BINDINGS_DIR}/main_module.cpp")
  foreach(domain ${PYOPENMS_DOMAINS})
    list(APPEND PYOPENMS_GENERATED_SOURCES "${PYOPENMS_BINDINGS_DIR}/bind_${domain}.cpp")
  endforeach()

  # Verify generated sources exist
  foreach(_src ${PYOPENMS_GENERATED_SOURCES})
    if(NOT EXISTS "${_src}")
      message(FATAL_ERROR "Generated binding source not found: ${_src}\n"
        "Run the generator first: python -m generator --pxd-dir pxds --addons-dir addons "
        "--output-dir bindings/generated --openms-include-dir <include-dirs>")
    endif()
  endforeach()

  # --------------------------------------------------------------------------
  # Build split nanobind modules for parallel compilation
  # --------------------------------------------------------------------------

  # Helper function to set common module properties
  function(_pyopenms_module_props target_name)
    target_include_directories(${target_name} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/bindings/type_casters
    )
    target_link_libraries(${target_name} PRIVATE
      ${OPENMS_TARGET}
      ${OPENSWATHALGO_TARGET}
    )
    set_target_properties(${target_name} PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms"
      RUNTIME_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms"
    )
    if(NOT PYOPENMS_PREPARE_WHEEL_REPAIR)
      set_target_properties(${target_name} PROPERTIES
        INSTALL_RPATH "${PY_RPATH}"
      )
      if(APPLE)
        set_target_properties(${target_name} PROPERTIES
          BUILD_WITH_INSTALL_NAME_DIR 1
        )
      endif()
    else()
      set_target_properties(${target_name} PROPERTIES
        INSTALL_RPATH_USE_LINK_PATH TRUE
        INSTALL_NAME_DIR ""
      )
    endif()
    install(TARGETS ${target_name}
      LIBRARY DESTINATION pyopenms
      RUNTIME DESTINATION pyopenms
      COMPONENT python_modules
    )
  endfunction()

  # Build each domain submodule
  foreach(domain ${PYOPENMS_DOMAINS})
    set(MODULE_NAME "_pyopenms_${domain}")

    nanobind_add_module(${MODULE_NAME}
      NB_DOMAIN pyopenms  # Share types across modules
      ${PYOPENMS_BINDINGS_DIR}/bind_${domain}.cpp
    )
    _pyopenms_module_props(${MODULE_NAME})
  endforeach()

  # Main entry module
  nanobind_add_module(_pyopenms
    NB_DOMAIN pyopenms
    ${PYOPENMS_BINDINGS_DIR}/main_module.cpp
  )
  _pyopenms_module_props(_pyopenms)

  # --------------------------------------------------------------------------
  # Arrow zero-copy export (only when OpenMS was built with Parquet/Arrow)
  # --------------------------------------------------------------------------
  get_target_property(_OPENMS_DEFS ${OPENMS_TARGET} INTERFACE_COMPILE_DEFINITIONS)
  set(_HAS_PARQUET FALSE)
  if(_OPENMS_DEFS)
    foreach(_def ${_OPENMS_DEFS})
      if(_def MATCHES "^WITH_PARQUET")
        set(_HAS_PARQUET TRUE)
      endif()
    endforeach()
  endif()
  # Fall back: check if WITH_PARQUET was set in the OpenMS build cache
  if(NOT _HAS_PARQUET AND WITH_PARQUET)
    set(_HAS_PARQUET TRUE)
  endif()
  if(_HAS_PARQUET)
    message(STATUS "pyOpenMS: Building _arrow_zerocopy module (WITH_PARQUET detected)")
    nanobind_add_module(_arrow_zerocopy
      NB_DOMAIN pyopenms
      ${CMAKE_CURRENT_SOURCE_DIR}/bindings/arrow_zerocopy.cpp
    )
    _pyopenms_module_props(_arrow_zerocopy)
  else()
    message(STATUS "pyOpenMS: Skipping _arrow_zerocopy module (WITH_PARQUET not found)")
  endif()

  # --------------------------------------------------------------------------
  # Aggregate compile target
  # --------------------------------------------------------------------------
  add_custom_target(pyopenms_compile ALL)
  add_dependencies(pyopenms_compile _pyopenms)
  foreach(domain ${PYOPENMS_DOMAINS})
    add_dependencies(pyopenms_compile "_pyopenms_${domain}")
  endforeach()
  if(_HAS_PARQUET)
    add_dependencies(pyopenms_compile _arrow_zerocopy)
  endif()

  # --------------------------------------------------------------------------
  # Copy Python package files to build directory
  # --------------------------------------------------------------------------

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

  # Core Python files
  _copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "*.py" "${PYOPENMS_BUILD_DIR}/pyopenms/")
  _copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "py.typed" "${PYOPENMS_BUILD_DIR}/pyopenms/")

  # Addon files
  _copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/addons/" "*.py" "${PYOPENMS_BUILD_DIR}/pyopenms/addons/")

  # Copy share data if requested and directory exists
  if(NOT NO_SHARE AND DEFINED OPENMS_SHARE_DIR AND EXISTS "${OPENMS_SHARE_DIR}")
    file(MAKE_DIRECTORY ${PYOPENMS_BUILD_DIR}/pyopenms/share)
    file(COPY ${OPENMS_SHARE_DIR} DESTINATION ${PYOPENMS_BUILD_DIR}/pyopenms/share
        PATTERN "examples" EXCLUDE)
  endif()

  # --------------------------------------------------------------------------
  # Stub generation
  # --------------------------------------------------------------------------
  if(PYOPENMS_GENERATE_STUBS)
    set(_STUB_DEPS _pyopenms)
    foreach(domain ${PYOPENMS_DOMAINS})
      list(APPEND _STUB_DEPS "_pyopenms_${domain}")
    endforeach()

    add_custom_target(pyopenms_stubs
      COMMAND ${Python_EXECUTABLE} -m nanobind.stubgen
              -m pyopenms
              -r                          # recurse into submodules
              -P                          # include private members
              -O ${PYOPENMS_BUILD_DIR}/pyopenms
              -M ${PYOPENMS_BUILD_DIR}/pyopenms/py.typed
              -p ${CMAKE_CURRENT_SOURCE_DIR}/stubs.patch
      COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fix_stubs.py
              ${PYOPENMS_BUILD_DIR}/pyopenms
      WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
      DEPENDS ${_STUB_DEPS}
      COMMENT "Generating .pyi stub files for pyopenms"
    )
  endif()

  # --------------------------------------------------------------------------
  # Platform-specific dependency copying and fixing
  # --------------------------------------------------------------------------
  if(NOT NO_DEPENDENCIES)
    set(PYOPENMS_DEPENDENCIES ${OPENMS_TARGET} ${OPENSWATHALGO_TARGET})

    # Step 1: Copy main OpenMS libraries (all platforms)
    add_custom_target(pyopenms_copy_deps ALL
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${OPENMS_TARGET}> "${PYOPENMS_BUILD_DIR}/pyopenms/"
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${OPENSWATHALGO_TARGET}> "${PYOPENMS_BUILD_DIR}/pyopenms/"
      COMMENT "Copying OpenMS dependencies to build dir"
    )
    add_dependencies(pyopenms_copy_deps pyopenms_compile)

    # Step 2: Platform-specific transitive dependency handling
    if(APPLE)
      add_custom_target(pyopenms_fix_deps ALL
        COMMAND ${PROJECT_SOURCE_DIR}/mac_fix_dependencies.rb -l "${PYOPENMS_BUILD_DIR}/pyopenms" -e "@rpath/" -r "@loader_path/" -f -v
        COMMAND codesign --force --sign - ${PYOPENMS_BUILD_DIR}/pyopenms/Qt* || true
        COMMAND codesign --force --sign - ${PYOPENMS_BUILD_DIR}/pyopenms/*.dylib || true
        COMMENT "Fixing deps and codesigning pyopenms binaries in build dir"
      )
      add_dependencies(pyopenms_fix_deps pyopenms_copy_deps)
    else()
      # Linux/Windows: Copy transitive dependencies (Qt, contrib libs, etc.)
      foreach(PYOPENMS_DEPENDENCY ${PYOPENMS_DEPENDENCIES})
        add_custom_command(
          TARGET pyopenms_copy_deps POST_BUILD
          COMMAND ${CMAKE_COMMAND} -DDEPS="$<TARGET_FILE:${PYOPENMS_DEPENDENCY}>" -DTARGET="${PYOPENMS_BUILD_DIR}/pyopenms/" -DLOOKUP_DIRS="${OPENMS_CONTRIB_LIBS}/lib\;${OPENMS_CONTRIB_LIBS}/lib64\;${QT_INSTALL_BINS}\;${QT_INSTALL_LIBS}" -P ${PROJECT_SOURCE_DIR}/pyopenms_copy_deps.cmake
          COMMENT "Copying ${PYOPENMS_DEPENDENCY} transitive dependencies"
        )
      endforeach()

      if(WIN32)
        # Try to find MSVS runtime libs
        include(InstallRequiredSystemLibraries)
        if("${MSVS_RTLIBS}" STREQUAL "")
          set(MSVS_RTLIBS "${CMAKE_INSTALL_SYSTEM_RUNTIME_LIBS}")
        endif()
        foreach(_runtime ${MSVS_RTLIBS})
          add_custom_command(
            TARGET pyopenms_copy_deps POST_BUILD
            COMMAND ${CMAKE_COMMAND} -E copy_if_different ${_runtime} "${PYOPENMS_BUILD_DIR}/pyopenms/"
            COMMENT "Copying ${_runtime}"
          )
        endforeach()
      endif()
    endif()
  endif()

  install(DIRECTORY ${PYOPENMS_BUILD_DIR}/pyopenms/
    DESTINATION pyopenms
    COMPONENT python_modules
    FILES_MATCHING
    PATTERN "*.py"
    PATTERN "*.pyi"
    PATTERN "py.typed"
    PATTERN "Qt*"
    PATTERN "__pycache__" EXCLUDE
  )
endif() # PYOPENMS_USE_PREBUILT

# --------------------------------------------------------------------------
# Install prebuilt artifacts (when PYOPENMS_USE_PREBUILT=ON)
# --------------------------------------------------------------------------
if(PYOPENMS_USE_PREBUILT)
  install(DIRECTORY ${PYOPENMS_BUILD_DIR}/pyopenms/
    DESTINATION pyopenms
    COMPONENT python_modules
    FILES_MATCHING
    PATTERN "*.py"
    PATTERN "*.pyi"
    PATTERN "py.typed"
    PATTERN "*.so"
    PATTERN "*.dylib"
    PATTERN "*.pyd"
    PATTERN "Qt*"
    PATTERN "__pycache__" EXCLUDE
  )
endif()

# Install share data if included
if(NOT NO_SHARE AND EXISTS "${PYOPENMS_BUILD_DIR}/pyopenms/share")
  install(DIRECTORY ${PYOPENMS_BUILD_DIR}/pyopenms/share
    DESTINATION pyopenms
    COMPONENT python_modules
  )
endif()

# --------------------------------------------------------------------------
# Unified pyopenms target
# --------------------------------------------------------------------------
if(NOT PYOPENMS_USE_PREBUILT)
  add_custom_target(pyopenms ALL)
  if(NOT NO_DEPENDENCIES)
    if(APPLE)
      add_dependencies(pyopenms pyopenms_fix_deps)
    else()
      add_dependencies(pyopenms pyopenms_copy_deps)
    endif()
  else()
    add_dependencies(pyopenms pyopenms_compile)
  endif()
else()
  add_custom_target(pyopenms ALL)
  add_dependencies(pyopenms pyopenms_compile)
endif()

# --------------------------------------------------------------------------
# Tests (only when building as part of OpenMS, not standalone with py-build-cmake)
# --------------------------------------------------------------------------
if(TARGET ${OPENMS_TARGET} OR TARGET OpenMS)
  # Copy tests to build directory so pytest doesn't import from source tree
  file(COPY ${PROJECT_SOURCE_DIR}/tests DESTINATION ${PYOPENMS_BUILD_DIR})

  enable_testing()
  add_test(
    NAME pyopenms_tests
    COMMAND ${Python_EXECUTABLE} -m pytest tests
    WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
  )
endif()

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
message(STATUS "")
message(STATUS "pyOpenMS Configuration Summary:")
message(STATUS "  Python:            ${Python_EXECUTABLE} (${Python_VERSION})")
message(STATUS "  NumPy:             ${Python_NumPy_INCLUDE_DIRS}")
if(nanobind_SOURCE_DIR)
  message(STATUS "  nanobind:          ${nanobind_SOURCE_DIR}")
else()
  message(STATUS "  nanobind:          (found via find_package)")
endif()
message(STATUS "  OpenMS:            ${OpenMS_DIR}")
message(STATUS "  Domain modules:    ${PYOPENMS_DOMAINS}")
message(STATUS "  Prebuilt:          ${PYOPENMS_USE_PREBUILT}")
message(STATUS "  No dependencies:   ${NO_DEPENDENCIES}")
message(STATUS "  Wheel repair:      ${PYOPENMS_PREPARE_WHEEL_REPAIR}")
message(STATUS "  Generate stubs:    ${PYOPENMS_GENERATE_STUBS}")
message(STATUS "  Parquet/Arrow:     ${_HAS_PARQUET}")
message(STATUS "")
