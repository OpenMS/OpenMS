# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)
project("OpenMS_class_tests_support")

#------------------------------------------------------------------------------
# libOpenMS-specific hooks for the library-agnostic class-test framework:
# deterministic UniqueIdGenerator seeding at START_TEST and XML schema
# validation of temporary files at END_TEST (see OpenMSTestSupport.cpp).
#
# Built as an OBJECT library: its static registrar must end up in every test
# executable, which a plain static library would not guarantee (linkers drop
# unreferenced object files from archives).
add_library(OpenMSTestSupport OBJECT OpenMSTestSupport.cpp)

# PUBLIC: consuming test executables inherit the framework and, transitively, libOpenMS
target_link_libraries(OpenMSTestSupport PUBLIC OpenMSTestFramework)

openms_add_library_compiler_flags(OpenMSTestSupport)
