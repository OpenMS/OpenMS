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
# Keep in sync with the minimum in the root CMakeLists.txt. Note that this is the
# floor for a plain pyOpenMS build; split mode (PYOPENMS_SPLIT_MODE, the default)
# additionally requires CMake 3.26 for FindPython's Development.SABIModule and
# checks for it below. That requirement is deliberately scoped to pyOpenMS -- the
# C++-only OpenMS build still needs no more than 3.24.
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

# Split mode compiles every extension against Python's stable ABI and moves
# nanobind's version-specific runtime out into the separate nanobind-backend
# package, so a single cp311-abi3 wheel per platform serves every supported
# CPython. Turning it off builds interpreter-specific modules instead, which is
# useful locally but must not be packaged as an abi3 wheel (see the guard below).
option(PYOPENMS_SPLIT_MODE "Build stable-ABI (abi3) extensions in nanobind split mode" ON)

# The stable-ABI floor has to be one number rather than two that happen to agree:
# CMake turns STABLE_ABI_VERSION into Py_LIMITED_API, py-build-cmake turns
# abi3_minimum_cpython_version into the wheel's cp<N>-abi3 tag, and nothing
# connects the two. Read the value out of pyproject.toml so the default here is
# that same number instead of a second copy that can drift away from it.
set(_pyopenms_abi3_floor "")
set(_pyopenms_pyproject "${CMAKE_CURRENT_SOURCE_DIR}/pyproject.toml")
if(EXISTS "${_pyopenms_pyproject}")
  file(STRINGS "${_pyopenms_pyproject}" _pyopenms_abi3_lines
       REGEX "^[ \t]*abi3_minimum_cpython_version[ \t]*=[ \t]*[0-9]+")
  # Exactly one key, or we cannot tell which one py-build-cmake will honour.
  list(LENGTH _pyopenms_abi3_lines _pyopenms_abi3_count)
  if(_pyopenms_abi3_count EQUAL 1 AND
     "${_pyopenms_abi3_lines}" MATCHES "=[ \t]*([0-9])([0-9]+)")
    set(_pyopenms_abi3_floor "${CMAKE_MATCH_1}.${CMAKE_MATCH_2}")
  endif()
endif()

set(PYOPENMS_STABLE_ABI_VERSION "${_pyopenms_abi3_floor}" CACHE STRING
    "CPython stable-ABI floor for split-mode builds (default: abi3_minimum_cpython_version from pyproject.toml)")

if(PYOPENMS_SPLIT_MODE)
  # Development.SABIModule arrived in CMake 3.26 and nanobind's split mode
  # requires that target on Windows. Only pyOpenMS needs it, so this is checked
  # here rather than raising the C++-only OpenMS minimum in the root project.
  if(CMAKE_VERSION VERSION_LESS 3.26)
    message(FATAL_ERROR
      "pyOpenMS split mode needs CMake 3.26 or newer for the Development.SABIModule "
      "component of find_package(Python); this is CMake ${CMAKE_VERSION}. Configure "
      "with -DPYOPENMS_SPLIT_MODE=OFF to build interpreter-specific modules instead.")
  endif()

  if(NOT PYOPENMS_STABLE_ABI_VERSION)
    message(FATAL_ERROR
      "Could not read a single abi3_minimum_cpython_version key from "
      "${_pyopenms_pyproject}, so the stable-ABI floor for split mode is unknown. Fix "
      "that key, or pass -DPYOPENMS_STABLE_ABI_VERSION=<major.minor> explicitly.")
  endif()

  # Overriding the floor only moves Py_LIMITED_API; the wheel keeps the tag
  # py-build-cmake derives from pyproject.toml. Packaging the two out of step
  # yields a cp311-abi3 wheel whose modules need a newer interpreter and fail to
  # import on 3.11. The knob is harmless outside packaging, so only refuse it
  # when py-build-cmake is the one driving this configure.
  if(DEFINED PY_BUILD_CMAKE_VERSION AND _pyopenms_abi3_floor AND
     NOT PYOPENMS_STABLE_ABI_VERSION VERSION_EQUAL "${_pyopenms_abi3_floor}")
    message(FATAL_ERROR
      "PYOPENMS_STABLE_ABI_VERSION is ${PYOPENMS_STABLE_ABI_VERSION}, but pyproject.toml "
      "tags the wheel for Python ${_pyopenms_abi3_floor} via abi3_minimum_cpython_version. "
      "Packaging that combination would produce a wheel advertising an ABI it was not built "
      "against. Change abi3_minimum_cpython_version as well if you really want a different "
      "stable-ABI floor.")
  endif()
else()
  # py-build-cmake tags the wheel cp311-abi3 from pyproject.toml, which only
  # describes a split-mode build. Refuse to package interpreter-specific modules
  # under that tag rather than shipping a wheel that lies about its ABI.
  if(DEFINED PY_BUILD_CMAKE_VERSION)
    message(FATAL_ERROR
      "PYOPENMS_SPLIT_MODE=OFF builds interpreter-specific extensions, but "
      "pyproject.toml tags the wheel as cp311-abi3. Building a wheel this way would "
      "mislabel it. Override [tool.py-build-cmake.wheel].python_abi too if you really "
      "want a linked-mode wheel.")
  endif()
endif()

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
# Windows DLL search path for in-tree imports
# --------------------------------------------------------------------------
# Windows has no RPATH, and since Python 3.8 an extension module's dependent
# DLLs are resolved only from the module's own directory, System32 and the
# directories registered with os.add_dll_directory() -- never from PATH.
# pyopenms/__init__.py registers everything listed in the PYOPENMS_DLL_PATH
# environment variable, so every in-tree `import pyopenms` (the ctest suite and
# the stub generator) has to be handed the directories holding OpenMS.dll and
# its transitive dependencies.
#
# Collected once here rather than at each use site, because the two hand-rolled
# lists that existed were each missing a different piece. Both were built from
# OPENMS_CONTRIB_LIBS and QT_INSTALL_BINS -- neither of which is defined when
# src/pyOpenMS is the top-level project under py-build-cmake -- and neither
# looked in a contrib *bin* directory at all. On Windows CMake installs DLLs to
# <prefix>/bin (RUNTIME) and import libraries to <prefix>/lib (ARCHIVE), so a
# contrib built with BUILD_SHARED_LIBS=ON puts its DLLs where a lib-only list
# can never find them. That is not hypothetical: the wheel workflow builds
# libcurl shared into contrib/bin, and libcurl.dll is a load-time import of
# OpenMS.dll, so every extension module failed to load with the loader's
# unhelpful "The specified module could not be found" (#10141).
function(_pyopenms_windows_dll_dirs out_var)
  # Caller-supplied directories come first: the environment is the only way an
  # enclosing build or a CI job can point at something CMake cannot know about.
  set(_dirs "")
  if(NOT "$ENV{PYOPENMS_DLL_PATH}" STREQUAL "")
    # A Windows PATH-style value, so ';' already makes it a CMake list.
    list(APPEND _dirs "$ENV{PYOPENMS_DLL_PATH}")
  endif()

  # Exact locations of the two libraries the modules link against. For an
  # imported OpenMS::OpenMS (the standalone wheel build) this resolves into the
  # install tree; in-tree it is the build output directory. Guarded because a
  # PYOPENMS_USE_PREBUILT configure never resolves the target names.
  foreach(_target IN ITEMS "${OPENMS_TARGET}" "${OPENSWATHALGO_TARGET}")
    if(NOT _target STREQUAL "" AND TARGET ${_target})
      list(APPEND _dirs "$<TARGET_FILE_DIR:${_target}>")
    endif()
  endforeach()

  # Every prefix this build was pointed at contributes both of its runtime
  # directories. Same transform cmake/package_general.cmake uses to feed
  # install(RUNTIME_DEPENDENCY_SET), and it is what picks up contrib/bin
  # without anyone having to name it. CMAKE_PREFIX_PATH is read as both a CMake
  # variable and an environment variable: find_package() consults both, but
  # only the variable is visible to ordinary script code, and the wheel
  # workflow supplies the environment one.
  set(_prefixes ${CMAKE_PREFIX_PATH})
  if(NOT "$ENV{CMAKE_PREFIX_PATH}" STREQUAL "")
    list(APPEND _prefixes "$ENV{CMAKE_PREFIX_PATH}")
  endif()
  foreach(_prefix IN LISTS _prefixes)
    if(NOT _prefix STREQUAL "")
      list(APPEND _dirs "${_prefix}/bin" "${_prefix}/lib")
    endif()
  endforeach()

  # An in-tree build knows the contrib prefix directly; a standalone one does not.
  if(OPENMS_CONTRIB_LIBS)
    list(APPEND _dirs
      "${OPENMS_CONTRIB_LIBS}/bin"
      "${OPENMS_CONTRIB_LIBS}/lib"
      "${OPENMS_CONTRIB_LIBS}/lib64")
  endif()
  # Only defined by a WITH_GUI build of the enclosing project (Qt is not linked
  # by the core library), so this is inert for the wheel build rather than the
  # dead weight it looked like when it was the only contrib-adjacent entry.
  if(QT_INSTALL_BINS)
    list(APPEND _dirs "${QT_INSTALL_BINS}")
  endif()
  if(QT_INSTALL_LIBS)
    list(APPEND _dirs "${QT_INSTALL_LIBS}")
  endif()

  # Only removes textual duplicates: two generator expressions that expand to the
  # same directory both survive. Registering a directory twice is harmless -- the
  # loop in pyopenms/__init__.py just calls os.add_dll_directory() again -- and
  # deduplicating after expansion is not worth it, because $<REMOVE_DUPLICATES>
  # would have to parse paths that can legitimately contain '(' and ')'.
  list(REMOVE_DUPLICATES _dirs)
  # '$<SEMICOLON>' rather than a literal ';': the whole list travels as ONE
  # command argument, and a bare ';' would split it into separate arguments.
  # Directories that do not exist are skipped by pyopenms/__init__.py, so the
  # list may name prefixes that have no bin/ or lib/ subdirectory.
  list(JOIN _dirs "$<SEMICOLON>" _joined)
  set(${out_var} "${_joined}" PARENT_SCOPE)
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
  set(PYOPENMS_PYTHON_COMPONENTS Interpreter Development.Module NumPy)
  if(PYOPENMS_SPLIT_MODE)
    # nanobind links split-mode extensions against Python::SABIModule on Windows.
    list(APPEND PYOPENMS_PYTHON_COMPONENTS Development.SABIModule)
  endif()
  find_package(Python 3.11 REQUIRED COMPONENTS ${PYOPENMS_PYTHON_COMPONENTS})

  # Targeting a stable ABI newer than the interpreter building us cannot work:
  # the Limited API its headers expose stops at their own version, so
  # Py_LIMITED_API would promise symbols that are not there. This is the one
  # reader of the ABI floor that cannot derive it -- cibuildwheel's TOML has no
  # interpolation, so [tool.cibuildwheel].build in pyproject.toml names the
  # frontend literally -- and therefore the one place it can drift. Fail here
  # rather than compile against a promise the headers do not keep.
  if(PYOPENMS_SPLIT_MODE AND Python_VERSION VERSION_LESS PYOPENMS_STABLE_ABI_VERSION)
    message(FATAL_ERROR
      "Split mode targets the CPython ${PYOPENMS_STABLE_ABI_VERSION} stable ABI, but this "
      "build runs Python ${Python_VERSION}, whose headers do not have it. Build with "
      "${PYOPENMS_STABLE_ABI_VERSION} or newer. When packaging, the interpreter is chosen by "
      "[tool.cibuildwheel].build in pyproject.toml (or CIBW_BUILD), which has to move with "
      "abi3_minimum_cpython_version.")
  endif()

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

  # Options shared by every pyOpenMS extension. Keeping them in one variable is
  # what guarantees all 15 modules land in the same nanobind domain and, in split
  # mode, target the same stable ABI -- a module added without them would silently
  # get its own type registry or an interpreter-specific tag.
  set(PYOPENMS_NB_MODULE_OPTIONS
    NB_DOMAIN pyopenms  # Share types across modules
    NOMINSIZE           # Use CMake's default -O3 instead of nanobind's -Os
  )
  if(PYOPENMS_SPLIT_MODE)
    # STABLE_ABI is what actually defines Py_LIMITED_API and gives the module the
    # .abi3 suffix; STABLE_ABI_VERSION only selects which CPython ABI to target,
    # and BACKEND_MODULE is what makes this a split build at all. All three are
    # needed -- STABLE_ABI_VERSION alone would leave a full-API module behind an
    # abi3 wheel tag.
    list(APPEND PYOPENMS_NB_MODULE_OPTIONS
      BACKEND_MODULE nanobind_backend
      STABLE_ABI
      STABLE_ABI_VERSION ${PYOPENMS_STABLE_ABI_VERSION}
    )
  endif()

  # Build each domain submodule
  foreach(domain ${PYOPENMS_DOMAINS})
    set(MODULE_NAME "_pyopenms_${domain}")

    nanobind_add_module(${MODULE_NAME}
      ${PYOPENMS_NB_MODULE_OPTIONS}
      ${PYOPENMS_BINDINGS_DIR}/bind_${domain}.cpp
    )
    _pyopenms_module_props(${MODULE_NAME})
  endforeach()

  # Main entry module
  nanobind_add_module(_pyopenms
    ${PYOPENMS_NB_MODULE_OPTIONS}
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
    ${PYOPENMS_NB_MODULE_OPTIONS}
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

  # Core Python files. py.typed is deliberately NOT copied here: it is the marker
  # that tells type checkers this package ships inline types, and nanobind.stubgen
  # writes it itself (its -M flag, below) only after the stubs were produced.
  # Copying it unconditionally made every wheel advertise types even when stub
  # generation was off or skipped -- which is what shipped, since the stub target
  # was never actually built (#10141).
  _copy_assets("${PROJECT_SOURCE_DIR}/pyopenms/" "*.py" "${PYOPENMS_BUILD_DIR}/pyopenms/")

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
    set(_pyopenms_stubgen_problem "")
    if(NOT _pyopenms_stubgen_nanobind_result EQUAL 0)
      set(_pyopenms_stubgen_problem
          "${Python_EXECUTABLE} cannot import nanobind")
    elseif(_pyopenms_stubgen_nanobind_version VERSION_LESS PYOPENMS_NANOBIND_MIN_VERSION)
      # One quoted argument: two would make this a CMake list and print as 'a;b'.
      set(_pyopenms_stubgen_problem
          "${Python_EXECUTABLE} has nanobind ${_pyopenms_stubgen_nanobind_version}, but the modules are built against nanobind ${PYOPENMS_NANOBIND_VERSION}")
    endif()

    # Now that the stub target is actually built (it is a dependency of the
    # aggregate pyopenms target), this is the difference between a build that
    # fails and one that quietly ships no stubs -- so it has to be decided here
    # rather than left as a warning about something that "will fail" later.
    #
    # Packaging is the case where a missing stub set is a defect in the artifact,
    # and nanobind is a declared build requirement there, so refuse to configure.
    # Any other build (a developer's in-tree OpenMS build, where nanobind may have
    # come from the FetchContent fallback and be unimportable from this
    # interpreter) keeps working, without stubs. Same split as the stable-ABI
    # guards above: fatal only when py-build-cmake is driving the configure.
    if(_pyopenms_stubgen_problem)
      if(DEFINED PY_BUILD_CMAKE_VERSION)
        message(FATAL_ERROR
          "Cannot generate .pyi stubs: ${_pyopenms_stubgen_problem}. A wheel built now "
          "would ship no type information. Install nanobind >= "
          "${PYOPENMS_NANOBIND_MIN_VERSION} for that interpreter, or configure with "
          "-DPYOPENMS_GENERATE_STUBS=OFF to build a deliberately untyped wheel.")
      endif()
      message(WARNING
        "Disabling stub generation: ${_pyopenms_stubgen_problem}. Install nanobind >= "
        "${PYOPENMS_NANOBIND_MIN_VERSION} for that interpreter to get .pyi files.")
      # A normal variable shadowing the cache entry, so the rest of this file --
      # including the configuration summary -- reports what will really happen.
      set(PYOPENMS_GENERATE_STUBS OFF)
    endif()
  endif()

  if(PYOPENMS_GENERATE_STUBS)
    # nanobind.stubgen generates stubs by IMPORTING the freshly built package, so
    # the modules and every native library behind them must be loadable at this
    # point in the build. On Windows that needs an explicit search path (see the
    # helper above); elsewhere RPATH already covers it and the variable is inert.
    set(_pyopenms_stub_cmd_prefix "")
    if(WIN32)
      _pyopenms_windows_dll_dirs(_pyopenms_stub_dll_dirs)
      set(_pyopenms_stub_cmd_prefix ${CMAKE_COMMAND} -E env
          "PYOPENMS_DLL_PATH=${_pyopenms_stub_dll_dirs}")
    endif()

    add_custom_target(pyopenms_stubs
      # Preflight. stubgen imports the package as a whole, and pyopenms/__init__.py
      # only raises when NOT ONE module imports -- a package missing half its
      # domains still produces stubs, silently describing an incomplete API. Check
      # every expected module first, and report which native library is missing
      # instead of leaving the loader's "The specified module could not be found".
      COMMAND ${_pyopenms_stub_cmd_prefix} ${Python_EXECUTABLE}
              ${CMAKE_CURRENT_SOURCE_DIR}/check_bindings.py
              --package-dir ${PYOPENMS_BUILD_DIR}/pyopenms
              --domains ${PYOPENMS_DOMAINS}
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
      # Required: PYOPENMS_DLL_PATH is one argument holding a ';'-separated list,
      # which the build tool would otherwise split (or hand to a shell that chokes
      # on it).
      VERBATIM
    )
  else()
    # Stubs are off. Make sure a build directory reused from an earlier run with
    # them on does not still carry .pyi files and the py.typed marker into the
    # wheel, where they would describe a package that no longer matches.
    file(GLOB _pyopenms_stale_stubs
         "${PYOPENMS_BUILD_DIR}/pyopenms/*.pyi"
         "${PYOPENMS_BUILD_DIR}/pyopenms/py.typed")
    if(_pyopenms_stale_stubs)
      message(STATUS "pyOpenMS: removing stale stub files (PYOPENMS_GENERATE_STUBS=OFF)")
      file(REMOVE ${_pyopenms_stale_stubs})
    endif()
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
          # contrib/bin is where a contrib built with BUILD_SHARED_LIBS=ON puts its
          # DLLs on Windows; without it the transitive copy silently misses them
          # (#10141). QT_INSTALL_* are only defined by a WITH_GUI build and are
          # empty entries otherwise, which pyopenms_copy_deps.cmake skips.
          COMMAND ${CMAKE_COMMAND} -DDEPS="$<TARGET_FILE:${PYOPENMS_DEPENDENCY}>" -DTARGET="${PYOPENMS_BUILD_DIR}/pyopenms/" -DLOOKUP_DIRS="${OPENMS_CONTRIB_LIBS}/bin\;${OPENMS_CONTRIB_LIBS}/lib\;${OPENMS_CONTRIB_LIBS}/lib64\;${QT_INSTALL_BINS}\;${QT_INSTALL_LIBS}" -P ${PROJECT_SOURCE_DIR}/pyopenms_copy_deps.cmake
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
  # The last step that puts a file into ${PYOPENMS_BUILD_DIR}/pyopenms. Stub
  # generation imports the package from there, so it has to run after this one:
  # when dependencies are bundled, the libraries the modules need are only
  # alongside them once the copy/fix step has finished.
  if(NOT NO_DEPENDENCIES)
    if(APPLE)
      set(_pyopenms_package_ready pyopenms_fix_deps)
    else()
      set(_pyopenms_package_ready pyopenms_copy_deps)
    endif()
  else()
    set(_pyopenms_package_ready pyopenms_compile)
  endif()
  add_dependencies(pyopenms ${_pyopenms_package_ready})

  # Stub generation is part of building pyOpenMS. It used to be neither in ALL nor
  # reachable from this target, so `.pyi` files were never produced during a wheel
  # build even with PYOPENMS_GENERATE_STUBS=ON -- the option only defined the
  # target (#10141).
  if(PYOPENMS_GENERATE_STUBS)
    add_dependencies(pyopenms_stubs ${_pyopenms_package_ready})
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
  # copied next to the compiled extension modules. Python >= 3.8 does not search
  # PATH for an extension module's dependent DLLs, and Windows has no RPATH, so
  # pyopenms/__init__.py registers the directories holding OpenMS.dll and its
  # transitive dependencies via os.add_dll_directory(), reading them from the
  # PYOPENMS_DLL_PATH environment variable. We set that variable for the test
  # process with `cmake -E env` so the TARGET_FILE_DIR generator expression (the
  # per-configuration OpenMS bin folder) is resolved; being an environment
  # variable it also reaches child processes that some tests spawn to
  # `import pyopenms`. Harmless no-op when the dependencies are already bundled
  # (NO_DEPENDENCIES=OFF) and on non-Windows.
  #
  # Same directory list as the stub generator uses -- in particular it now covers
  # the contrib *bin* directory, where a contrib built with shared libraries puts
  # its DLLs (#10141); this list previously named only contrib/lib.
  set(_pyopenms_test_cmd_prefix "")
  if(WIN32)
    _pyopenms_windows_dll_dirs(_pyopenms_dll_dirs)
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
