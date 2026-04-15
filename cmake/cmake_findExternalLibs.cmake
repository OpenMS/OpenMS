# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------


#------------------------------------------------------------------------------
# This cmake file handles finding external libs for OpenMS
#------------------------------------------------------------------------------


#------------------------------------------------------------------------------
# find libs (for linking)
# On Windows:
#   * on windows we need the *.lib versions (dlls alone won't do for linking)
#   * never mix Release/Debug versions of libraries. Leads to strange segfaults,
#     stack corruption etc, due to different runtime libs ...
# compiler-wise: use the same compiler for contrib and OpenMS!
find_package(XercesC REQUIRED)

#------------------------------------------------------------------------------
# BOOST
set(OpenMS_BOOST_COMPONENTS date_time regex CACHE INTERNAL "Boost components for core lib")
find_boost(iostreams ${OpenMS_BOOST_COMPONENTS})

if(Boost_FOUND)
  message(STATUS "Found Boost version ${Boost_MAJOR_VERSION}.${Boost_MINOR_VERSION}.${Boost_SUBMINOR_VERSION}" )
  set(CF_OPENMS_BOOST_VERSION_MAJOR ${Boost_MAJOR_VERSION})
  set(CF_OPENMS_BOOST_VERSION_MINOR ${Boost_MINOR_VERSION})
  set(CF_OPENMS_BOOST_VERSION_SUBMINOR ${Boost_SUBMINOR_VERSION})
  set(CF_OPENMS_BOOST_VERSION ${Boost_VERSION})

  get_target_property(location Boost::iostreams LOCATION)
  get_target_property(target_type Boost::iostreams TYPE)
  if (target_type STREQUAL "STATIC_LIBRARY" AND location MATCHES "^/usr/local/")
    message(WARNING "Statically linked Boost from system installations like brew, are not fully supported yet.
Either use '-DBOOST_USE_STATIC=OFF' to use the shared library or build boost with our contrib. Nonetheless,
we are going to try to continue building.")
    get_target_property(libs Boost::iostreams INTERFACE_LINK_LIBRARIES)
    # If boost from brew, replace simple "link flags" like "-lzstd" with
    # find_package calls and their resulting imported targets
    # since boost CMake does not expose this transitive dependency as targets!
    # see https://github.com/boostorg/boost_install/issues/64
    foreach (lib ${libs})
      if (lib MATCHES "zstd")
        find_package(zstd)
      elseif (lib MATCHES "lzma")
        find_package(LibLZMA)
      endif()
    endforeach ()
    ##
    set_target_properties(Boost::iostreams
          PROPERTIES INTERFACE_LINK_LIBRARIES "BZip2::BZip2;ZLIB::ZLIB;zstd::libzstd_shared;LibLZMA::LibLZMA")
  endif()

  get_target_property(location Boost::regex LOCATION)
  get_target_property(target_type Boost::regex TYPE)
  if (target_type STREQUAL "STATIC_LIBRARY" AND location MATCHES "^/usr/local/")
    get_target_property(libs Boost::regex INTERFACE_LINK_LIBRARIES)
    # If boost from brew, replace simple "link flags" like "-lzstd" with
    # find_package calls and their resulting imported targets
    # since boost CMake does not expose this transitive dependency as targets!
    # see https://github.com/boostorg/boost_install/issues/64
    foreach (lib ${libs})
      if (lib MATCHES "icui18n")
        find_package(ICU COMPONENTS "data" "uc" "i18n")
      endif()
    endforeach ()
    ##
    set_target_properties(Boost::regex
            PROPERTIES INTERFACE_LINK_LIBRARIES "ICU::data;ICU:uc;ICU::i18n")
  endif()
else()
  message(FATAL_ERROR "Boost or one of its components not found!")
endif()

#------------------------------------------------------------------------------
# libsvm
# creates LibSVM target
find_package(LIBSVM 2.91 REQUIRED) ## will not overwrite LIBSVM_LIBRARY if defined already
if (LIBSVM_FOUND)
	set(CF_OPENMS_LIBSVM_VERSION_MAJOR ${LIBSVM_MAJOR_VERSION})
	set(CF_OPENMS_LIBSVM_VERSION_MINOR ${LIBSVM_MINOR_VERSION})
	set(CF_OPENMS_LIBSVM_VERSION ${LIBSVM_VERSION})
endif()

#------------------------------------------------------------------------------
# COIN-OR
# Our find module creates an imported CoinOR::CoinOR target
find_package(COIN)
if (COIN_FOUND)
  set(OPENMS_HAS_COINOR 1)
  set(LPTARGET "CoinOR::CoinOR")
else()
  #------------------------------------------------------------------------------
  # GLPK
  # creates GLPK::GLPK target
  find_package(GLPK)
  if (GLPK_FOUND)
    set(CF_OPENMS_GLPK_VERSION_MAJOR ${GLPK_VERSION_MAJOR})
    set(CF_OPENMS_GLPK_VERSION_MINOR ${GLPK_VERSION_MINOR})
    set(CF_OPENMS_GLPK_VERSION ${GLPK_VERSION_STRING})
    set(LPTARGET "GLPK::GLPK")
  else()
    message(FATAL_ERROR "Either COIN-OR or GLPK has to be available (COIN-OR takes precedence).")
  endif()
endif()

#------------------------------------------------------------------------------
# zlib
# creates ZLIB::ZLIB target
find_package(ZLIB REQUIRED)
# Detect zlib-ng (drop-in compat mode). zlib-ng produces byte-different but
# functionally equivalent compressed output, which causes FuzzyDiff-based
# tests to fail when comparing mzML binary arrays.
include(CheckSymbolExists)
set(CMAKE_REQUIRED_INCLUDES ${ZLIB_INCLUDE_DIRS})
check_symbol_exists(ZLIBNG_VERSION "zlib.h" _HAVE_ZLIB_NG)
if(_HAVE_ZLIB_NG)
  set(HAVE_ZLIB_NG TRUE CACHE BOOL "zlib-ng detected as system zlib" FORCE)
  message(STATUS "Detected zlib-ng -- will relax compressed-data test comparisons")
else()
  set(HAVE_ZLIB_NG FALSE CACHE BOOL "zlib-ng detected as system zlib" FORCE)
endif()

#------------------------------------------------------------------------------
# bzip2
find_package(BZip2 REQUIRED)

#------------------------------------------------------------------------------
# libzip (ZIP64 archive support)
# Uses our FindLibzip.cmake module which does a manual header+library search.
# We intentionally avoid CONFIG mode because libzip <= 1.10 ships a CMake
# targets file that references CLI tool binaries (zipcmp, zipmerge, ziptool)
# and raises FATAL_ERROR if those aren't installed (e.g. Ubuntu without
# libzip-tools).
find_package(Libzip REQUIRED)

#------------------------------------------------------------------------------
# Find Eigen
# creates Eigen3::Eigen3 package
# CMake is garbage https://gitlab.kitware.com/cmake/cmake/-/issues/24581
# We also have to check for 3.4.0 because the configversion file in 3.4.0 tells CMake that 5 is incompatible
# Try to find any Eigen3 in the range [3.4.0, 6), quietly
# The package is still called Eigen3.. don't ask!
find_package(Eigen3 3.4.0...<6 QUIET)

if(TARGET Eigen3::Eigen)
    message(STATUS "Found Eigen3 version in range 3.4.0...<6: ${Eigen3_VERSION}")
else()
    # Fall back to the usual version compatibility check (for old/broken configversion files)
    find_package(Eigen3 3.4.0 REQUIRED) # fail if not found now
    if(TARGET Eigen3::Eigen)
        message(STATUS "Found Eigen3 version compatible to 3.4.0: ${Eigen3_VERSION}")
    endif()
endif()

#------------------------------------------------------------------------------
# Find Crawdad libraries if requested
# cmake args: -DCrawdad_DIR=/path/to/Crawdad/ -DWITH_CRAWDAD=TRUE
## TODO check if necessary
if (WITH_CRAWDAD)
  message(STATUS "Will compile with Crawdad support: ${Crawdad_DIR}" )
  find_package(Crawdad REQUIRED)
  # find archive (static) version and add it to the OpenMS library
  find_library(Crawdad_LIBRARY
   NAMES Crawdad.a Crawdad
    HINTS ${Crawdad_DIR})
endif()

#------------------------------------------------------------------------------
# HDF5
if (WITH_HDF5)
  # For MSVC use static linking to the HDF5 libraries
  if(MSVC)
    set(HDF5_USE_STATIC_LIBRARIES ON)
  endif()
  find_package(HDF5 MODULE REQUIRED COMPONENTS C CXX)
endif()


#------------------------------------------------------------------------------
# libcurl (used for HTTP in OpenMS core; also needed by Apache Arrow 23+)
find_package(CURL REQUIRED)

#------------------------------------------------------------------------------
# Apache Arrow and Parquet (required dependency)
# Arrow 23+ required for parquet file format compatibility
find_package(Arrow 23 CONFIG REQUIRED)
find_package(Parquet 23 CONFIG REQUIRED)

# Arrow's CMake config may import nlohmann_json as a transitive dependency.
# If so, force the vendored extern to use the already-imported target instead
# of trying to add_library() a second target with the same name.
if(TARGET nlohmann_json::nlohmann_json)
  set(USE_EXTERNAL_JSON ON CACHE BOOL "Use an external nlohmann::json library" FORCE)
endif()

# Determine Arrow target based on ARROW_USE_STATIC preference
if(ARROW_USE_STATIC AND TARGET Arrow::arrow_static)
  set(OPENMS_ARROW_TARGET Arrow::arrow_static)
elseif(NOT ARROW_USE_STATIC AND TARGET Arrow::arrow_shared)
  set(OPENMS_ARROW_TARGET Arrow::arrow_shared)
elseif(TARGET Arrow::arrow_static)
  set(OPENMS_ARROW_TARGET Arrow::arrow_static)
elseif(TARGET Arrow::arrow_shared)
  set(OPENMS_ARROW_TARGET Arrow::arrow_shared)
else()
  message(FATAL_ERROR "No suitable Arrow target found")
endif()

# Determine Arrow Compute target (may be bundled into Arrow::arrow_* on some distros).
# We need compute kernels (e.g., "equal") for filter expression binding.
if(ARROW_USE_STATIC AND TARGET Arrow::arrow_compute_static)
  set(OPENMS_ARROW_COMPUTE_TARGET Arrow::arrow_compute_static)
elseif(NOT ARROW_USE_STATIC AND TARGET Arrow::arrow_compute_shared)
  set(OPENMS_ARROW_COMPUTE_TARGET Arrow::arrow_compute_shared)
elseif(TARGET Arrow::arrow_compute_static)
  set(OPENMS_ARROW_COMPUTE_TARGET Arrow::arrow_compute_static)
elseif(TARGET Arrow::arrow_compute_shared)
  set(OPENMS_ARROW_COMPUTE_TARGET Arrow::arrow_compute_shared)
else()
  # Fallback: compute might be a plain library without a CMake target.
  # Use platform-neutral search with optional hints.
  find_library(OPENMS_ARROW_COMPUTE_LIB
    NAMES arrow_compute
    HINTS ${CMAKE_PREFIX_PATH} ${ARROW_HOME}
    PATH_SUFFIXES lib lib64
  )
  if(OPENMS_ARROW_COMPUTE_LIB)
    set(OPENMS_ARROW_COMPUTE_TARGET ${OPENMS_ARROW_COMPUTE_LIB})
    message(STATUS "Using Arrow Compute library from find_library: ${OPENMS_ARROW_COMPUTE_LIB}")
  else()
    # Last resort: compute is part of the main Arrow target.
    set(OPENMS_ARROW_COMPUTE_TARGET ${OPENMS_ARROW_TARGET})
    message(STATUS "Using Arrow Compute from Arrow target: ${OPENMS_ARROW_TARGET} (no separate compute target/library found)")
  endif()
endif()

# Determine Parquet target based on ARROW_USE_STATIC preference
if(ARROW_USE_STATIC AND TARGET Parquet::parquet_static)
  set(OPENMS_PARQUET_TARGET Parquet::parquet_static)
elseif(NOT ARROW_USE_STATIC AND TARGET Parquet::parquet_shared)
  set(OPENMS_PARQUET_TARGET Parquet::parquet_shared)
elseif(TARGET Parquet::parquet_static)
  set(OPENMS_PARQUET_TARGET Parquet::parquet_static)
elseif(TARGET Parquet::parquet_shared)
  set(OPENMS_PARQUET_TARGET Parquet::parquet_shared)
else()
  message(FATAL_ERROR "No suitable Parquet target found")
endif()

message(STATUS "Using Arrow target: ${OPENMS_ARROW_TARGET}")
message(STATUS "Using Arrow Compute target: ${OPENMS_ARROW_COMPUTE_TARGET}")
message(STATUS "Using Parquet target: ${OPENMS_PARQUET_TARGET}")

# Optional Arrow Dataset (for predicate pushdown)
find_package(ArrowDataset CONFIG QUIET)
if(ArrowDataset_FOUND)
  if(ARROW_USE_STATIC AND TARGET ArrowDataset::arrow_dataset_static)
    set(OPENMS_ARROW_DATASET_TARGET ArrowDataset::arrow_dataset_static)
  elseif(NOT ARROW_USE_STATIC AND TARGET ArrowDataset::arrow_dataset_shared)
    set(OPENMS_ARROW_DATASET_TARGET ArrowDataset::arrow_dataset_shared)
  elseif(TARGET ArrowDataset::arrow_dataset_static)
    set(OPENMS_ARROW_DATASET_TARGET ArrowDataset::arrow_dataset_static)
  elseif(TARGET ArrowDataset::arrow_dataset_shared)
    set(OPENMS_ARROW_DATASET_TARGET ArrowDataset::arrow_dataset_shared)
  endif()

  if(OPENMS_ARROW_DATASET_TARGET)
    message(STATUS "Using Arrow Dataset target: ${OPENMS_ARROW_DATASET_TARGET}")
    # Arrow Dataset (static) may pull in libxml2 symbols; link explicitly.
    # This avoids missing xmlBufferFree at runtime when dataset pushdown is enabled.
    find_package(LibXml2 REQUIRED)
  endif()
endif()

#------------------------------------------------------------------------------
# Done finding contrib libraries
#------------------------------------------------------------------------------

#except for the contrib libs, prefer shared libraries
if(NOT MSVC AND NOT APPLE)
	set(CMAKE_FIND_LIBRARY_SUFFIXES ".so;.a")
endif()

#------------------------------------------------------------------------------
# PTHREAD
#------------------------------------------------------------------------------
# Prefer the -pthread compiler flag to be consistent with SQLiteCpp and avoid
# rebuilds
# TODO Do we even need this, when OpenMP is not active?
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package (Threads REQUIRED)


#------------------------------------------------------------------------------
# QT (only needed for GUI)
#------------------------------------------------------------------------------
SET(QT_MIN_VERSION "6.1.0")

if (WITH_GUI)
  find_package(Qt6 ${QT_MIN_VERSION} COMPONENTS Core QUIET)

  IF (Qt6Core_FOUND)
    message(STATUS "Found Qt ${Qt6Core_VERSION}")
  ELSE()
    message(FATAL_ERROR "Qt6Core not found — required when WITH_GUI=ON. Use -DWITH_GUI=OFF to build without GUI.")
  ENDIF()

  # --------------------------------------------------------------------------
  # Find additional Qt libs
  #---------------------------------------------------------------------------
  set (TEMP_OpenMS_GUI_QT_COMPONENTS Gui Widgets Svg OpenGLWidgets)

  # On macOS the platform plugin of QT requires PrintSupport. We link
  # so it's packaged via the bundling/dependency tools/scripts
  if (APPLE)
    set (TEMP_OpenMS_GUI_QT_COMPONENTS ${TEMP_OpenMS_GUI_QT_COMPONENTS} PrintSupport)
  endif()

  set(OpenMS_GUI_QT_COMPONENTS ${TEMP_OpenMS_GUI_QT_COMPONENTS} CACHE INTERNAL "QT components for GUI lib")

  if(NOT NO_WEBENGINE_WIDGETS)
    set(OpenMS_GUI_QT_COMPONENTS_OPT WebEngineWidgets)
  endif()

  find_package(Qt6 REQUIRED COMPONENTS ${OpenMS_GUI_QT_COMPONENTS})

  IF (NOT Qt6Widgets_FOUND OR NOT Qt6Gui_FOUND OR NOT Qt6Svg_FOUND)
    message(STATUS "Qt6Widgets not found!")
    message(FATAL_ERROR "To find a custom Qt installation use: cmake <..more options..> -DCMAKE_PREFIX_PATH='<path_to_parent_folder_of_lib_folder_withAllQt6Libs>' <src-dir>")
  ENDIF()
  find_package(Qt6 QUIET COMPONENTS ${OpenMS_GUI_QT_COMPONENTS_OPT})

  # TODO only works if WebEngineWidgets is the only optional component
  set(OpenMS_GUI_QT_FOUND_COMPONENTS_OPT)
  if(Qt6WebEngineWidgets_FOUND)
    list(APPEND OpenMS_GUI_QT_FOUND_COMPONENTS_OPT "WebEngineWidgets")
  else()
    message(WARNING "Qt6WebEngineWidgets not found or disabled, disabling JS Views in TOPPView!")
  endif()

  set(OpenMS_GUI_DEP_LIBRARIES "OpenMS")

  foreach(COMP IN LISTS OpenMS_GUI_QT_COMPONENTS)
    list(APPEND OpenMS_GUI_DEP_LIBRARIES "Qt6::${COMP}")
  endforeach()

  foreach(COMP IN LISTS OpenMS_GUI_QT_FOUND_COMPONENTS_OPT)
    list(APPEND OpenMS_GUI_DEP_LIBRARIES "Qt6::${COMP}")
  endforeach()

endif()

#------------------------------------------------------------------------------
# opentims (Bruker TimsTOF .d file reading)
if (WITH_OPENTIMS)
  # Enable C language for bundled ZSTD fallback (zstddeclib.c)
  enable_language(C)

  find_package(Opentims QUIET)

  if(Opentims_FOUND)
    message(STATUS "opentims: using system installation (${Opentims_LIBRARIES})")

    # If the installed library was built with OPENTIMS_LINK_SQLITE_STATICALLY,
    # that compile definition is propagated via the imported target's
    # INTERFACE_COMPILE_DEFINITIONS. OpenMS must then supply sqlite3 headers
    # and the library, just as it does in the FetchContent path.
    # When the library was found via the manual search path (no CMake config
    # package), INTERFACE_COMPILE_DEFINITIONS may not be set; in that case we
    # conservatively inject sqlite3 because we cannot know how it was built.
    get_target_property(_opentims_defs opentims::opentims_cpp INTERFACE_COMPILE_DEFINITIONS)
    if(_opentims_defs MATCHES "OPENTIMS_LINK_SQLITE_STATICALLY"
       OR ((_opentims_defs MATCHES "NOTFOUND" OR _opentims_defs STREQUAL "") AND NOT opentims_FOUND))
      target_include_directories(opentims::opentims_cpp INTERFACE
        "${CMAKE_SOURCE_DIR}/src/openms/extern/SQLiteCpp/sqlite3")
      target_link_libraries(opentims::opentims_cpp INTERFACE sqlite3)
      message(STATUS "opentims: injecting OpenMS sqlite3 (library was built with static sqlite)")
    endif()
  else()
    # No system install found — fetch and build from source.
    message(STATUS "opentims: system installation not found, fetching from git")
    include(FetchContent)

    FetchContent_Declare(
      opentims
      GIT_REPOSITORY https://github.com/michalsta/opentims.git
      GIT_TAG        v1.2.0b4
    )

    # Build opentims as a C++ static library, not a Python module.
    # Use OpenMS's own sqlite3 instead of opentims's runtime dlopen.
    set(OPENTIMS_BUILD_LIB              ON  CACHE BOOL "" FORCE)
    set(OPENTIMS_BUILD_PYTHON           OFF CACHE BOOL "" FORCE)
    set(OPENTIMS_LINK_SQLITE_STATICALLY ON  CACHE BOOL "" FORCE)

    # OpenMS sets BUILD_SHARED_LIBS=true (to build libOpenMS.so). Without this
    # override, opentims's CMakeLists.txt would build libopentims_cpp.so instead
    # of libopentims_cpp.a. A shared opentims would land in _deps/opentims-build/
    # rather than the system lib path, so the build-time linker would fail with
    # "undefined reference to TimsDataHandle" when linking test binaries against
    # libOpenMS.so (because libopentims_cpp.so's directory is not in the test's
    # -L search path). Force static so all opentims symbols get absorbed into
    # libOpenMS.so at link time.
    set(_openms_saved_build_shared_libs ${BUILD_SHARED_LIBS})
    set(BUILD_SHARED_LIBS OFF)
    FetchContent_MakeAvailable(opentims)
    set(BUILD_SHARED_LIBS ${_openms_saved_build_shared_libs})

    # Provide OpenMS's sqlite3 headers to opentims so it compiles with
    # OPENTIMS_LINK_SQLITE_STATICALLY (direct sqlite3 calls vs. dlopen).
    # We intentionally do NOT call target_link_libraries(opentims_cpp PRIVATE sqlite3)
    # here: adding the CMake sqlite3 target as a PRIVATE dep of a STATIC library
    # would require sqlite3 to appear in the CMake export sets, which creates
    # conflicts with SQLiteCpp's own sqlite3 export. The sqlite3 symbols are
    # resolved at final link time through OpenMS's own SQLiteCpp dependency
    # (SQLiteCpp PUBLIC-links sqlite3, so libsqlite3.a already appears in
    # libOpenMS.so's link command after libopentims_cpp.a).
    target_include_directories(opentims_cpp PRIVATE
      "${CMAKE_SOURCE_DIR}/src/openms/extern/SQLiteCpp/sqlite3")

    # ZSTD: prefer system; fall back to opentims's bundled decoder.
    set(_OPENTIMS_SRC "${opentims_SOURCE_DIR}/src/opentims++")
    find_package(zstd QUIET)
    if(TARGET zstd::libzstd_shared)
      target_link_libraries(opentims_cpp PRIVATE zstd::libzstd_shared)
      message(STATUS "opentims: using system zstd (shared)")
    elseif(TARGET zstd::libzstd_static)
      target_link_libraries(opentims_cpp PRIVATE zstd::libzstd_static)
      message(STATUS "opentims: using system zstd (static)")
    else()
      target_sources(opentims_cpp PRIVATE "${_OPENTIMS_SRC}/zstd/zstddeclib.c")
      target_include_directories(opentims_cpp PRIVATE "${_OPENTIMS_SRC}/zstd")
      message(STATUS "opentims: using bundled zstd decoder (system zstd not found)")
    endif()

    # Suppress warnings from third-party code
    target_compile_options(opentims_cpp PRIVATE $<IF:$<CXX_COMPILER_ID:MSVC>,/w,-w>)

    # Expose a namespaced alias so downstream CMakeLists always use
    # opentims::opentims_cpp regardless of how the library was obtained.
    add_library(opentims::opentims_cpp ALIAS opentims_cpp)

    # Add opentims_cpp to both the install-tree and build-tree export sets
    # (same pattern as Evergreen, IsoSpec, and other bundled deps).
    # install_library() → OpenMSTargets install export
    # openms_register_export_target() → _OPENMS_EXPORT_TARGETS build-tree export()
    install_library(opentims_cpp)
    openms_register_export_target(opentims_cpp)

    message(STATUS "opentims: built from source (${opentims_SOURCE_DIR})")
  endif()
endif()
#------------------------------------------------------------------------------
