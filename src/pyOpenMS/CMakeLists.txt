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
# Bindings are hand-maintained nanobind C++ sources in bindings/*.cpp.
# --------------------------------------------------------------------------

# Required when this file is the top-level project (py-build-cmake / pip install .).
# Keep in sync with the minimum in the root CMakeLists.txt.
cmake_minimum_required(VERSION 3.24 FATAL_ERROR)
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
    spectrum chromatogram experiment
)

# Set default build directory (may be overridden by PYOPENMS_USE_PREBUILT)
if(NOT DEFINED PYOPENMS_BUILD_DIR OR "${PYOPENMS_BUILD_DIR}" STREQUAL "")
  set(PYOPENMS_BUILD_DIR "${CMAKE_BINARY_DIR}/pyOpenMS")
endif()

# --------------------------------------------------------------------------
# Windows DLL search path
# --------------------------------------------------------------------------
# Windows has no RPATH, and since Python 3.8 an extension module's dependent
# DLLs are no longer looked up on PATH. pyopenms/__init__.py therefore registers
# a list of directories with os.add_dll_directory(), read from the
# PYOPENMS_DLL_PATH environment variable. Both stub generation and the test
# suite have to import the freshly built package, so both need that list;
# assembling it in one place keeps the two from drifting apart.
#
# The directories are the ones OpenMS' own packaging searches for runtime
# dependencies (see RUNTIME_DEP_SEARCH_DIRS in cmake/package_general.cmake):
# every prefix with both /bin and /lib appended. Listing both halves matters
# here -- on Windows CMake installs shared libraries to <prefix>/bin (RUNTIME)
# and import libraries to <prefix>/lib (ARCHIVE), so a dependency installed the
# CMake way puts its DLL in contrib/bin, while the prebuilt contrib keeps its
# DLLs next to the import libraries in contrib/lib. Naming only one of the two
# silently loses whichever dependency lives in the other; libcurl, the one
# contrib library the wheel workflow builds from source, lands in contrib/bin.
#
# Qt is deliberately absent: cmake_findExternalLibs.cmake guards the entire Qt
# block with if(WITH_GUI) and the core library links no Qt target, so QT_INSTALL_BINS
# is not even defined in the builds that feed pyOpenMS and no extension module
# has a Qt dependency to resolve.
function(_pyopenms_windows_dll_dirs _out_var)
  if(NOT WIN32)
    set(${_out_var} "" PARENT_SCOPE)
    return()
  endif()

  # Build-tree locations of the OpenMS libraries. Generator expressions, so the
  # per-configuration directory of a multi-config generator is resolved.
  set(_dirs "")
  foreach(_target IN ITEMS ${OPENMS_TARGET} ${OPENSWATHALGO_TARGET})
    if(TARGET ${_target})
      list(APPEND _dirs "$<TARGET_FILE_DIR:${_target}>")
    endif()
  endforeach()

  # Third-party DLLs. A ";"-separated environment variable expands into several
  # list elements here, which is exactly what CMAKE_PREFIX_PATH needs: CMake
  # reads that one from the environment for find_package(), but does not make it
  # visible as a CMake variable, and a standalone wheel build only ever sets the
  # environment variable.
  foreach(_prefix IN ITEMS ${OPENMS_CONTRIB_LIBS} ${CMAKE_PREFIX_PATH} $ENV{CMAKE_PREFIX_PATH})
    if(_prefix)
      # Normalise separators by hand. file(TO_CMAKE_PATH) is not usable here:
      # it splits its argument on the *host's* path-list separator, which turns
      # a native "C:\\dir" into the two entries "C" and "/dir" whenever the
      # host separator is ":".
      string(REPLACE "\\" "/" _prefix "${_prefix}")
      list(APPEND _dirs "${_prefix}/bin" "${_prefix}/lib")
    endif()
  endforeach()

  # Anything the caller supplied is kept rather than overwritten.
  if(NOT "$ENV{PYOPENMS_DLL_PATH}" STREQUAL "")
    list(APPEND _dirs $ENV{PYOPENMS_DLL_PATH})
  endif()

  list(REMOVE_DUPLICATES _dirs)
  # The result is one os.pathsep-separated string, not a CMake list: the ";" has
  # to survive into the command line, which is what $<SEMICOLON> is for.
  string(REPLACE ";" "$<SEMICOLON>" _joined "${_dirs}")
  set(${_out_var} "${_joined}" PARENT_SCOPE)
endfunction()

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

  # Keep in sync with requires-python in pyproject.toml (pandas 3 also needs >= 3.11)
  find_package(Python 3.11 REQUIRED COMPONENTS Interpreter Development.Module NumPy)
  find_package(Eigen3 REQUIRED)

  # Propagate Python_EXECUTABLE to parent scope so the pyopenms_wheel target can use it
  set(Python_EXECUTABLE "${Python_EXECUTABLE}" PARENT_SCOPE)

  # --------------------------------------------------------------------------
  # Find or fetch nanobind
  # --------------------------------------------------------------------------
  # pyOpenMS requires nanobind 3: the custom type casters in
  # bindings/type_casters/ implement the nanobind 3 caster protocol
  # (from_python(handle, uint32_t, cleanup_list *)). Keep the minimum in sync
  # with the nanobind requirement in pyproject.toml.
  set(PYOPENMS_NANOBIND_MIN_VERSION "3.0.1")
  set(PYOPENMS_NANOBIND_FETCH_TAG "v3.0.1")

  if(COMMAND nanobind_add_module)
    # An enclosing project already set nanobind up. Check it instead of silently
    # bypassing the version requirement. nanobind_VERSION is only defined when
    # nanobind was located through find_package; otherwise read the version
    # macros nanobind.h defines, so the check compares a complete version rather
    # than merely establishing that this is some nanobind 3.
    set(_pyopenms_nb_version "${nanobind_VERSION}")
    if(NOT _pyopenms_nb_version)
      set(_pyopenms_nb_header "${NB_DIR}/include/nanobind/nanobind.h")
      if(EXISTS "${_pyopenms_nb_header}")
        file(READ "${_pyopenms_nb_header}" _pyopenms_nb_header_content)
        set(_pyopenms_nb_parts "")
        foreach(_pyopenms_nb_component MAJOR MINOR PATCH)
          # A component can legitimately be 0, so collect matches and count them
          # instead of testing each captured value for truth.
          if(_pyopenms_nb_header_content MATCHES
             "#[ \t]*define[ \t]+NB_VERSION_${_pyopenms_nb_component}[ \t]+([0-9]+)")
            list(APPEND _pyopenms_nb_parts "${CMAKE_MATCH_1}")
          endif()
        endforeach()
        list(LENGTH _pyopenms_nb_parts _pyopenms_nb_part_count)
        if(_pyopenms_nb_part_count EQUAL 3)
          list(JOIN _pyopenms_nb_parts "." _pyopenms_nb_version)
        endif()
      endif()
    endif()

    if(NOT _pyopenms_nb_version)
      message(FATAL_ERROR
        "pyOpenMS requires nanobind >= ${PYOPENMS_NANOBIND_MIN_VERSION}, but the version of the "
        "nanobind supplied by the enclosing project could not be determined "
        "(NB_DIR='${NB_DIR}'). Supply it via "
        "find_package(nanobind ${PYOPENMS_NANOBIND_MIN_VERSION} CONFIG) so the version is known.")
    elseif(_pyopenms_nb_version VERSION_LESS PYOPENMS_NANOBIND_MIN_VERSION)
      message(FATAL_ERROR
        "pyOpenMS requires nanobind >= ${PYOPENMS_NANOBIND_MIN_VERSION}, but the enclosing "
        "project supplied nanobind ${_pyopenms_nb_version}.")
    endif()
    set(PYOPENMS_NANOBIND_VERSION "${_pyopenms_nb_version}")
  else()
    # Ask the *selected* interpreter where its nanobind lives, so the frontend we
    # compile against is the one belonging to the interpreter that later runs
    # nanobind.stubgen. A caller-supplied nanobind_ROOT/nanobind_DIR wins, but it
    # still has to satisfy the version requirement below.
    if(NOT nanobind_ROOT AND NOT nanobind_DIR)
      execute_process(
        COMMAND "${Python_EXECUTABLE}" -m nanobind --cmake_dir
        OUTPUT_STRIP_TRAILING_WHITESPACE
        OUTPUT_VARIABLE nanobind_ROOT
        RESULT_VARIABLE _nanobind_cmake_dir_result
        ERROR_QUIET)
      if(NOT _nanobind_cmake_dir_result EQUAL 0)
        message(STATUS "No nanobind importable from ${Python_EXECUTABLE}")
        set(nanobind_ROOT "")
      endif()
    endif()

    # Request the minimum explicitly: nanobind ships a SameMajorVersion config,
    # so this accepts >= 3.0.1 and < 4 and rejects an installed 2.x, which the
    # previous unversioned find_package would have selected.
    find_package(nanobind ${PYOPENMS_NANOBIND_MIN_VERSION} CONFIG QUIET)
    if(nanobind_FOUND)
      set(PYOPENMS_NANOBIND_VERSION "${nanobind_VERSION}")
    else()
      message(STATUS
        "No nanobind >= ${PYOPENMS_NANOBIND_MIN_VERSION} found for ${Python_EXECUTABLE}; "
        "fetching nanobind ${PYOPENMS_NANOBIND_FETCH_TAG}")
      include(FetchContent)
      FetchContent_Declare(
          nanobind
          GIT_REPOSITORY https://github.com/wjakob/nanobind.git
          GIT_TAG        ${PYOPENMS_NANOBIND_FETCH_TAG}
      )
      FetchContent_MakeAvailable(nanobind)
      set(PYOPENMS_NANOBIND_VERSION "${PYOPENMS_NANOBIND_FETCH_TAG} (fetched)")
    endif()
  endif()
  message(STATUS "pyOpenMS: using nanobind ${PYOPENMS_NANOBIND_VERSION} from ${NB_DIR}")

  # Patch nanobind's __forceinline → inline on MSVC to avoid O(n²) optimizer
  # blowup that causes multi-hour compile times on large binding files.
  # NB_DIR is set by nanobind-config.cmake regardless of install method
  # (FetchContent or pip). See https://github.com/wjakob/nanobind/discussions/791
  if(MSVC AND DEFINED NB_DIR)
    set(_nb_defs_h "${NB_DIR}/include/nanobind/nb_defs.h")
    if(EXISTS "${_nb_defs_h}")
      file(READ "${_nb_defs_h}" _nb_defs)
      string(FIND "${_nb_defs}" "__forceinline" _fi_pos)
      if(NOT _fi_pos EQUAL -1)
        string(REPLACE "__forceinline" "inline" _nb_defs "${_nb_defs}")
        file(WRITE "${_nb_defs_h}" "${_nb_defs}")
        message(STATUS "Patched nanobind nb_defs.h: __forceinline → inline (MSVC compile-time fix)")
      endif()
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
  # Hand-maintained nanobind binding sources
  # --------------------------------------------------------------------------
  set(PYOPENMS_BINDINGS_DIR "${CMAKE_CURRENT_SOURCE_DIR}/bindings")

  set(PYOPENMS_BINDING_SOURCES "${PYOPENMS_BINDINGS_DIR}/main_module.cpp")
  foreach(domain ${PYOPENMS_DOMAINS})
    list(APPEND PYOPENMS_BINDING_SOURCES "${PYOPENMS_BINDINGS_DIR}/bind_${domain}.cpp")
  endforeach()

  # Verify binding sources exist
  foreach(_src ${PYOPENMS_BINDING_SOURCES})
    if(NOT EXISTS "${_src}")
      message(FATAL_ERROR "Binding source not found: ${_src}")
    endif()
  endforeach()

  # --------------------------------------------------------------------------
  # Build split nanobind modules for parallel compilation
  # --------------------------------------------------------------------------

  # Helper function to set common module properties
  function(_pyopenms_module_props target_name)
    target_include_directories(${target_name} PRIVATE
      ${CMAKE_CURRENT_SOURCE_DIR}/bindings/type_casters
      ${EIGEN3_INCLUDE_DIR}
    )
    target_link_libraries(${target_name} PRIVATE
      ${OPENMS_TARGET}
      ${OPENSWATHALGO_TARGET}
      Eigen3::Eigen
    )
    # The trailing "$<0:>" is an (empty) generator expression: it makes CMake use
    # this path verbatim instead of appending a per-configuration subdirectory on
    # multi-config generators (Visual Studio, Xcode). Without it the modules land
    # in pyopenms/Release/, where pyopenms/__init__.py (pkgutil.iter_modules, which
    # does not recurse) cannot find them, and `import pyopenms` fails with
    # "Failed to import any pyOpenMS bindings".
    set_target_properties(${target_name} PROPERTIES
      LIBRARY_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms/$<0:>"
      RUNTIME_OUTPUT_DIRECTORY "${PYOPENMS_BUILD_DIR}/pyopenms/$<0:>"
    )
    if(NOT PYOPENMS_PREPARE_WHEEL_REPAIR)
      set_target_properties(${target_name} PROPERTIES
        INSTALL_RPATH "${PY_RPATH}"
      )
      if(NOT NO_DEPENDENCIES)
        # When dependencies are bundled alongside the modules, use the install
        # RPATH at build time so $ORIGIN/ finds the co-located libraries.
        set_target_properties(${target_name} PROPERTIES
          BUILD_WITH_INSTALL_RPATH 1
        )
      endif()
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
      LIBRARY DESTINATION pyopenms COMPONENT python_modules
      RUNTIME DESTINATION pyopenms COMPONENT python_modules
    )
  endfunction()

  # Build each domain submodule
  foreach(domain ${PYOPENMS_DOMAINS})
    set(MODULE_NAME "_pyopenms_${domain}")

    nanobind_add_module(${MODULE_NAME}
      NB_DOMAIN pyopenms  # Share types across modules
      NOMINSIZE            # Use CMake's default -O3 instead of nanobind's -Os
      ${PYOPENMS_BINDINGS_DIR}/bind_${domain}.cpp
    )
    _pyopenms_module_props(${MODULE_NAME})
  endforeach()

  # Main entry module
  nanobind_add_module(_pyopenms
    NB_DOMAIN pyopenms
    NOMINSIZE
    ${PYOPENMS_BINDINGS_DIR}/main_module.cpp
  )
  _pyopenms_module_props(_pyopenms)

  # --------------------------------------------------------------------------
  # Arrow zero-copy export
  # --------------------------------------------------------------------------
  # In standalone builds, OPENMS_ARROW_TARGET is not exported by OpenMSConfig.cmake.
  # Find Arrow ourselves so we can link the C++ bridge functions (ExportSchema, etc.).
  if(NOT OPENMS_ARROW_TARGET)
    # Don't pin the version in find_package — Arrow's ConfigVersion uses
    # SameMajorVersion compatibility, so requesting 23 rejects 24.x. Check the
    # version manually after the call. Mirrors the main cmake path in
    # cmake/cmake_findExternalLibs.cmake.
    find_package(Arrow CONFIG REQUIRED)
    if(Arrow_VERSION VERSION_LESS 23)
      message(FATAL_ERROR "Apache Arrow >= 23 required (found ${Arrow_VERSION}).")
    endif()
    if(ARROW_USE_STATIC AND TARGET Arrow::arrow_static)
      set(OPENMS_ARROW_TARGET Arrow::arrow_static)
    elseif(NOT ARROW_USE_STATIC AND TARGET Arrow::arrow_shared)
      set(OPENMS_ARROW_TARGET Arrow::arrow_shared)
    elseif(TARGET Arrow::arrow_static)
      set(OPENMS_ARROW_TARGET Arrow::arrow_static)
    elseif(TARGET Arrow::arrow_shared)
      set(OPENMS_ARROW_TARGET Arrow::arrow_shared)
    else()
      message(FATAL_ERROR "No suitable Arrow target found for _arrow_zerocopy")
    endif()
  endif()

  message(STATUS "pyOpenMS: Building _arrow_zerocopy module")
  nanobind_add_module(_arrow_zerocopy
    NB_DOMAIN pyopenms
    NOMINSIZE
    ${CMAKE_CURRENT_SOURCE_DIR}/bindings/arrow_zerocopy.cpp
  )
  _pyopenms_module_props(_arrow_zerocopy)
  # arrow_zerocopy.cpp uses Arrow C++ API (ExportSchema, ImportRecordBatch, etc.)
  # for zero-copy Table conversion between C++ and PyArrow.
  target_link_libraries(_arrow_zerocopy PRIVATE ${OPENMS_ARROW_TARGET})

  # --------------------------------------------------------------------------
  # Aggregate compile target
  # --------------------------------------------------------------------------
  add_custom_target(pyopenms_compile ALL)
  add_dependencies(pyopenms_compile _pyopenms)
  foreach(domain ${PYOPENMS_DOMAINS})
    add_dependencies(pyopenms_compile "_pyopenms_${domain}")
  endforeach()
  add_dependencies(pyopenms_compile _arrow_zerocopy)

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
  # py.typed is deliberately NOT copied here. It advertises that the package
  # ships type information, so it must only appear when the stubs were actually
  # produced -- nanobind.stubgen writes the marker itself (-M below). Clear a
  # marker and stubs left behind by an earlier build in the same directory that
  # did generate them, so disabling stub generation cannot leave a wheel
  # claiming types it no longer contains.
  if(NOT PYOPENMS_GENERATE_STUBS)
    file(REMOVE "${PYOPENMS_BUILD_DIR}/pyopenms/py.typed")
    file(GLOB _pyopenms_stale_stubs "${PYOPENMS_BUILD_DIR}/pyopenms/*.pyi")
    if(_pyopenms_stale_stubs)
      file(REMOVE ${_pyopenms_stale_stubs})
    endif()
  endif()

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
    # nanobind.stubgen is imported from ${Python_EXECUTABLE}, the same interpreter
    # used for nanobind discovery above. If that interpreter cannot import a
    # nanobind matching the frontend the modules were compiled against (which
    # happens when the FetchContent fallback kicked in), the generated stubs would
    # come from a different nanobind release than the modules.
    execute_process(
      COMMAND "${Python_EXECUTABLE}" -c "import nanobind; print(nanobind.__version__)"
      OUTPUT_STRIP_TRAILING_WHITESPACE
      OUTPUT_VARIABLE _pyopenms_stubgen_nanobind_version
      RESULT_VARIABLE _pyopenms_stubgen_nanobind_result
      ERROR_QUIET)
    if(NOT _pyopenms_stubgen_nanobind_result EQUAL 0)
      message(WARNING
        "PYOPENMS_GENERATE_STUBS=ON, but ${Python_EXECUTABLE} cannot import nanobind; "
        "stub generation will fail. Install nanobind >= ${PYOPENMS_NANOBIND_MIN_VERSION} "
        "for that interpreter or configure with -DPYOPENMS_GENERATE_STUBS=OFF.")
    elseif(_pyopenms_stubgen_nanobind_version VERSION_LESS PYOPENMS_NANOBIND_MIN_VERSION)
      message(WARNING
        "Stub generation uses nanobind ${_pyopenms_stubgen_nanobind_version} from "
        "${Python_EXECUTABLE}, but the modules are built against nanobind "
        "${PYOPENMS_NANOBIND_VERSION}. Install nanobind >= "
        "${PYOPENMS_NANOBIND_MIN_VERSION} for that interpreter.")
    endif()

    # nanobind.stubgen works by importing the freshly built package, so on
    # Windows every extension module has to be loadable at this point -- Linux
    # and macOS resolve OpenMS through RPATH instead.
    _pyopenms_windows_dll_dirs(_pyopenms_stub_dll_dirs)
    set(_pyopenms_stub_cmd_prefix "")
    if(_pyopenms_stub_dll_dirs)
      set(_pyopenms_stub_cmd_prefix ${CMAKE_COMMAND} -E env "PYOPENMS_DLL_PATH=${_pyopenms_stub_dll_dirs}")
      # Report what the environment contributed. When stubgen cannot import the
      # package on Windows the message is only ever "DLL load failed ... The
      # specified module could not be found", and whether a caller-supplied
      # value arrived is the one thing it can never tell you. Deliberately a
      # configure-time message and not a build-time `cmake -E echo` of the
      # resolved list: that value holds ";"-separated directories and generator
      # expressions containing "(", both of which the build tool's shell
      # interprets, which fails with a syntax error and takes the stub target
      # down with it.
      if("$ENV{PYOPENMS_DLL_PATH}" STREQUAL "")
        message(STATUS "pyOpenMS stubs: no PYOPENMS_DLL_PATH in the environment; "
                       "using the OpenMS target, contrib and CMAKE_PREFIX_PATH directories")
      else()
        message(STATUS "pyOpenMS stubs: PYOPENMS_DLL_PATH from the environment: $ENV{PYOPENMS_DLL_PATH}")
      endif()
    endif()

    # Every module the package is expected to provide. pyopenms/__init__.py
    # reports a failed domain import as a warning and only raises when *every*
    # module failed, so stubgen can otherwise succeed on a partial import and
    # write stubs with whole domains silently missing.
    set(_pyopenms_expected_modules _pyopenms)
    foreach(_domain IN LISTS PYOPENMS_DOMAINS)
      list(APPEND _pyopenms_expected_modules "_pyopenms_${_domain}")
    endforeach()
    if(TARGET _arrow_zerocopy)
      list(APPEND _pyopenms_expected_modules _arrow_zerocopy)
    endif()

    add_custom_target(pyopenms_stubs
      COMMAND ${_pyopenms_stub_cmd_prefix} ${Python_EXECUTABLE}
              ${CMAKE_CURRENT_SOURCE_DIR}/check_module_imports.py
              ${PYOPENMS_BUILD_DIR}/pyopenms ${_pyopenms_expected_modules}
      COMMAND ${_pyopenms_stub_cmd_prefix} ${Python_EXECUTABLE} -m nanobind.stubgen
              -m pyopenms
              -r                          # recurse into submodules
              -P                          # include private members
              -O ${PYOPENMS_BUILD_DIR}/pyopenms
              -M ${PYOPENMS_BUILD_DIR}/pyopenms/py.typed
              -p ${CMAKE_CURRENT_SOURCE_DIR}/stubs.patch
      COMMAND ${Python_EXECUTABLE} ${CMAKE_CURRENT_SOURCE_DIR}/fix_stubs.py
              ${PYOPENMS_BUILD_DIR}/pyopenms
      WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
      DEPENDS pyopenms_compile
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
          COMMAND ${CMAKE_COMMAND} -DDEPS="$<TARGET_FILE:${PYOPENMS_DEPENDENCY}>" -DTARGET="${PYOPENMS_BUILD_DIR}/pyopenms/" -DLOOKUP_DIRS="${OPENMS_CONTRIB_LIBS}/bin\;${OPENMS_CONTRIB_LIBS}/lib\;${OPENMS_CONTRIB_LIBS}/lib64" -P ${PROJECT_SOURCE_DIR}/pyopenms_copy_deps.cmake
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

    if(PYOPENMS_GENERATE_STUBS)
      # Stub generation imports pyopenms, which is only possible once the
      # bundled libraries sit next to the modules.
      if(APPLE)
        add_dependencies(pyopenms_stubs pyopenms_fix_deps)
      else()
        add_dependencies(pyopenms_stubs pyopenms_copy_deps)
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
  if(PYOPENMS_GENERATE_STUBS)
    # Without this the target is reachable from neither ALL nor the aggregate
    # target, and every wheel build skips stub generation silently.
    add_dependencies(pyopenms pyopenms_stubs)
  endif()
else()
  add_custom_target(pyopenms ALL)
  add_dependencies(pyopenms pyopenms_compile)
endif()

# --------------------------------------------------------------------------
# Tests (only when building as part of OpenMS, not standalone with py-build-cmake)
# --------------------------------------------------------------------------
if((DEFINED OPENMS_TARGET AND TARGET ${OPENMS_TARGET}) OR TARGET OpenMS)
  # Copy tests to build directory so pytest doesn't import from source tree
  file(COPY ${PROJECT_SOURCE_DIR}/tests DESTINATION ${PYOPENMS_BUILD_DIR})

  enable_testing()

  # Split pyopenms tests into separate CTest entries for better CDash reporting.
  # Use -q (quiet) and --tb=short to avoid cluttering CDash with verbose passing output.
  # Override pyproject.toml's addopts which uses -v (verbose).
  set(_PYTEST_ARGS -m pytest)

  # On Windows, when NO_DEPENDENCIES=ON the OpenMS runtime libraries are not
  # copied next to the compiled extension modules, so the test process needs the
  # same PYOPENMS_DLL_PATH the stub step uses (see _pyopenms_windows_dll_dirs
  # above). We set it with `cmake -E env` so the TARGET_FILE_DIR generator
  # expression -- the per-configuration OpenMS bin folder -- is resolved; being
  # an environment variable it also reaches the child processes some tests spawn
  # to `import pyopenms`. Empty, and therefore a no-op, off Windows.
  set(_pyopenms_test_cmd_prefix "")
  _pyopenms_windows_dll_dirs(_pyopenms_dll_dirs)
  if(_pyopenms_dll_dirs)
    set(_pyopenms_test_cmd_prefix ${CMAKE_COMMAND} -E env "PYOPENMS_DLL_PATH=${_pyopenms_dll_dirs}")
  endif()

  add_test(
    NAME pyopenms_unittests
    COMMAND ${_pyopenms_test_cmd_prefix} ${Python_EXECUTABLE} ${_PYTEST_ARGS}
      --ignore=tests/unittests/test_docstrings.py
      --ignore=tests/unittests/test_docstring_format.py
      tests/unittests
    WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
  )

  # Mirror the BrukerTimsFile_test pattern: forward optional CMake-configured
  # paths to real Bruker .d directories as environment variables, so
  # test_BrukerTimsFile.py can skip when they are unset.
  if(OPENTIMS_DDA_TEST_DATA OR OPENTIMS_DIA_TEST_DATA)
    set(_pyopenms_test_env "")
    if(OPENTIMS_DDA_TEST_DATA)
      list(APPEND _pyopenms_test_env "OPENTIMS_DDA_TEST_DATA=${OPENTIMS_DDA_TEST_DATA}")
    endif()
    if(OPENTIMS_DIA_TEST_DATA)
      list(APPEND _pyopenms_test_env "OPENTIMS_DIA_TEST_DATA=${OPENTIMS_DIA_TEST_DATA}")
    endif()
    set_tests_properties(pyopenms_unittests PROPERTIES
      ENVIRONMENT "${_pyopenms_test_env}"
    )
  endif()

  add_test(
    NAME pyopenms_docstring_tests
    COMMAND ${_pyopenms_test_cmd_prefix} ${Python_EXECUTABLE} ${_PYTEST_ARGS}
      tests/unittests/test_docstrings.py
      tests/unittests/test_docstring_format.py
    WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
  )

  add_test(
    NAME pyopenms_integration_tests
    COMMAND ${_pyopenms_test_cmd_prefix} ${Python_EXECUTABLE} ${_PYTEST_ARGS} tests/integration_tests
    WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
  )

  # Top-level test files (ProForma, modification discovery, matrix zerocopy, etc.)
  file(GLOB _pyopenms_toplevel_tests "${PYOPENMS_BUILD_DIR}/tests/test_*.py")
  if(_pyopenms_toplevel_tests)
    add_test(
      NAME pyopenms_feature_tests
      COMMAND ${_pyopenms_test_cmd_prefix} ${Python_EXECUTABLE} ${_PYTEST_ARGS} ${_pyopenms_toplevel_tests}
      WORKING_DIRECTORY ${PYOPENMS_BUILD_DIR}
    )
  endif()
endif()

# --------------------------------------------------------------------------
# Summary
# --------------------------------------------------------------------------
message(STATUS "")
message(STATUS "pyOpenMS Configuration Summary:")
message(STATUS "  Python:            ${Python_EXECUTABLE} (${Python_VERSION})")
message(STATUS "  NumPy:             ${Python_NumPy_INCLUDE_DIRS}")
if(nanobind_SOURCE_DIR)
  message(STATUS "  nanobind:          ${PYOPENMS_NANOBIND_VERSION} (${nanobind_SOURCE_DIR})")
else()
  message(STATUS "  nanobind:          ${PYOPENMS_NANOBIND_VERSION} (${NB_DIR})")
endif()
message(STATUS "  OpenMS:            ${OpenMS_DIR}")
message(STATUS "  Domain modules:    ${PYOPENMS_DOMAINS}")
message(STATUS "  Prebuilt:          ${PYOPENMS_USE_PREBUILT}")
message(STATUS "  No dependencies:   ${NO_DEPENDENCIES}")
message(STATUS "  Wheel repair:      ${PYOPENMS_PREPARE_WHEEL_REPAIR}")
message(STATUS "  Generate stubs:    ${PYOPENMS_GENERATE_STUBS}")
message(STATUS "  Parquet/Arrow:     ON")
message(STATUS "")
