# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow $
# $Authors: Andreas Bertsch, Chris Bielow, Stephan Aiche $
# --------------------------------------------------------------------------

# Vendored, namespace-isolated Percolator subtree (OpenMS::Internal::Percolator).
#
# This is built as a standalone STATIC library with its OWN compile flags and
# linked PRIVATE into libOpenMS (see ../../CMakeLists.txt). Compiling it as a
# separate target — instead of folding the sources into libOpenMS — is what
# isolates this upstream code from the OpenMS target's "treat warning as error"
# flags (MSVC /we4100 + /we4189, GCC -Werror).
#
# Previously the sources were added to OpenMS_sources and the per-file
# suppressions (/wd4100 ...) had to "out-order" the target-wide /we4100 on the
# same command line. CMake does not guarantee that ordering across generators /
# versions, so a clean build was effectively luck: in some setups /we4100 won
# and the vendored code failed to compile. A dedicated target never receives
# /we4100 in the first place, so the suppressions are no longer order-sensitive.
#
# Side benefit: the upstream sources show up under this target in IDEs instead
# of cluttering the OpenMS target's file list at the project root.

project("OpenMS_Percolator")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# Source list (also consumed by the upstream-sync tooling). It PREPENDs this
# directory to each entry, so it resolves correctly when included from here.
include(${CMAKE_CURRENT_LIST_DIR}/sources.cmake)

add_library(OpenMS_Percolator STATIC ${OpenMS_percolator_sources})

# Compile at the same C++ standard as libOpenMS so the vendored headers — which
# are also pulled into OpenMS's Percolator pimpl translation unit — behave
# identically on both sides of the link.
target_compile_features(OpenMS_Percolator PRIVATE cxx_std_23)

# The vendored headers include each other as #include "X.h" (no subdirectory).
# Expose this directory PUBLIC so the dependency propagates to libOpenMS's pimpl
# TU (PercolatorImpl.h) through the link below. BUILD_INTERFACE only: these
# headers are OpenMS-internal and never installed.
target_include_directories(OpenMS_Percolator PUBLIC
                           "$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>")

# External dependencies used by the vendored code: Eigen + Boost (header-only
# here) and OpenMP for the #pragma omp regions. PRIVATE: these are build-time
# usage requirements of this lib, not part of anything it exposes.
target_link_libraries(OpenMS_Percolator PRIVATE Eigen3::Eigen Boost::boost)
if (OPENMP_FOUND)
  target_link_libraries(OpenMS_Percolator PRIVATE OpenMP::OpenMP_CXX)
endif()

# Linked into the shared libOpenMS, so the objects must be position-independent.
set_target_properties(OpenMS_Percolator PROPERTIES POSITION_INDEPENDENT_CODE ON)

# Keep the vendored symbols out of libOpenMS's dynamic symbol table, matching
# the previous in-tree -fvisibility=hidden treatment of these sources.
set_target_properties(OpenMS_Percolator PROPERTIES
                      CXX_VISIBILITY_PRESET hidden
                      VISIBILITY_INLINES_HIDDEN 1)

# Suppress upstream-originated warnings on this target ONLY. GCC/Clang use
# -Wno-*; MSVC needs /wd<code> and chokes on the -Wno-* spellings (D8021).
# Because this is a separate target, it never inherits OpenMS's /we* (errors),
# so these are plain suppressions with no ordering dependency.
if (MSVC)
  target_compile_options(OpenMS_Percolator PRIVATE
    /wd4018 /wd4100 /wd4189 /wd4244 /wd4267 /wd4996 /wd4458)
else()
  target_compile_options(OpenMS_Percolator PRIVATE
    -Wno-sign-compare -Wno-unused-parameter -Wno-deprecated-declarations
    -Wno-unused-variable -Wno-overloaded-virtual -Wno-reorder -Wno-reorder-ctor
    -Wno-uninitialized -Wno-ignored-qualifiers -Wno-unused-but-set-variable)
endif()

# Note: intentionally NOT installed/exported. libOpenMS is always shared
# (BUILD_SHARED_LIBS is forced ON) and links this PRIVATE, so the archive is
# fully absorbed into libOpenMS.so and never appears in the OpenMS export set.
