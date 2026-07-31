# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Divyansh Lalotra $
# $Authors: Divyansh Lalotra $
# --------------------------------------------------------------------------

project("OpenMS_benchmark_tests")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# Benchmark tests for large datasets (Issue #8788).
#
# These tests download large external datasets (multi-GB) and are meant to run
# on a schedule (e.g. nightly), not on every PR. They are enabled explicitly
# at configure time with -DENABLE_BENCHMARK_TESTING=ON (default OFF).
#
# Naming scheme (test names must be distinct):
#   BENCH_<tool>_<dataset>[_<num>]
#
# The BENCH_ prefix distinguishes benchmark tests from TOPP_/unit tests in
# nightly builds and CDash, mirroring the rationale documented for the TOPP_
# prefix in src/tests/topp/CMakeLists.txt (see issue #8788).
#
# This is the PR 1 skeleton: one registered test that proves the CTest wiring.
# Dataset bootstrap (FetchContent + URL_HASH), metric extraction, baseline
# comparison and CI scheduling land in later PRs.

#------------------------------------------------------------------------------
# General settings
#------------------------------------------------------------------------------
set(TESTS_TEMP_DIR "${PROJECT_BINARY_DIR}/tmp_path")
file(MAKE_DIRECTORY "${TESTS_TEMP_DIR}")

#------------------------------------------------------------------------------
# Skeleton test (placeholder)
#
# Proves the benchmark test harness is registered with CTest and selectable
# via `ctest -L benchmark`. Will be replaced by a real benchmark
# (download -> run -> extract -> compare -> report) in later PRs.
#------------------------------------------------------------------------------
add_test(NAME "BENCH_SKELETON"
         COMMAND ${CMAKE_COMMAND} -E echo "benchmark test skeleton registered")
set_tests_properties("BENCH_SKELETON" PROPERTIES LABELS "benchmark")
