# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Oliver Kohlbacher $
# $Authors: Oliver Kohlbacher $
# --------------------------------------------------------------------------

# Vendored mzPeak C++ library (https://github.com/okohlbacher/mzpeak-openms, a fork of
# https://github.com/OpenMS/mzpeak; see README.OpenMS.md for the commit and the local changes).
#
# Built as a standalone STATIC library (OpenMS_mzpeak) with its OWN compile flags and linked PRIVATE into
# libOpenMS (see ../../CMakeLists.txt), like the vendored Percolator: the OpenMS target's warning-as-error
# flags do not reach the upstream code, and its headers never appear in an installed OpenMS header
# (MzPeakFile wraps the library behind a pimpl). Only MzPeakFile.cpp includes the library.

project("OpenMS_mzpeak")

include(${CMAKE_CURRENT_LIST_DIR}/sources.cmake)

add_library(OpenMS_mzpeak STATIC ${OpenMS_mzpeak_sources})

# The library is C++23 (std::unreachable, std::format, std::ranges); OpenMS itself
# compiles with cxx_std_23.
target_compile_features(OpenMS_mzpeak PRIVATE cxx_std_23)
# Linked into the shared libOpenMS, with its symbol visibility: nothing of the library or of the Boost.JSON compiled
# into it is exported from libOpenMS (a consumer with its own Boost.JSON or mzPeak must not bind to these)
set_target_properties(OpenMS_mzpeak PROPERTIES POSITION_INDEPENDENT_CODE ON
                                               CXX_VISIBILITY_PRESET hidden
                                               VISIBILITY_INLINES_HIDDEN ON)

# The library's headers (#include "mzpeak/..."), PUBLIC so that the one OpenMS translation unit that includes
# them (MzPeakFile.cpp) finds them through the link. BUILD_INTERFACE only: never installed.
target_include_directories(OpenMS_mzpeak PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>")
# src/util/numpress.cpp uses the MSNumpress copy that libOpenMS compiles (OpenMS/FORMAT/MSNUMPRESS), instead of
# the library's meson subproject, so that libOpenMS holds one copy of the codec.
target_include_directories(OpenMS_mzpeak PRIVATE
                           "${OPENMS_HOST_DIRECTORY}/src/openms/include")

# Dependencies: Arrow + Parquet (as libOpenMS links them), libzip, Boost (header-only here: Boost.JSON is
# compiled into this library from boost/json/src.hpp, see sources.cmake; Boost.TypeIndex, Boost.Range,
# Boost.Iostreams templates, Boost.Math); no newer Boost than OpenMS requires (1.81).
target_link_libraries(OpenMS_mzpeak PRIVATE
                      ${OPENMS_ARROW_TARGET}
                      ${OPENMS_PARQUET_TARGET}
                      libzip::zip
                      Boost::boost
                      Threads::Threads)
if (WIN32)
  # <windows.h> reaches in through Arrow; keep its min/max macros out (as the library's meson.build does)
  target_compile_definitions(OpenMS_mzpeak PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN)
endif()
if (MSVC)
  # Boost.JSON is compiled as one translation unit (src/boost_json_src.cpp), as OpenMS' own targets get it
  target_compile_options(OpenMS_mzpeak PRIVATE /bigobj)
endif()

# The library's own warning settings (meson.build): warnings are not errors here.
if (NOT MSVC)
  target_compile_options(OpenMS_mzpeak PRIVATE -Wall -Wextra -Wno-unused-parameter)
endif()
