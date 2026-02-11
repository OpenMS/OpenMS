# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Chris Bielow $
# $Authors: Chris Bielow $
# --------------------------------------------------------------------------


project("OpenMS-simde")
cmake_minimum_required(VERSION 3.10.0 FATAL_ERROR)

add_library(SIMDe INTERFACE)
target_include_directories(SIMDe SYSTEM INTERFACE ${PROJECT_SOURCE_DIR})
