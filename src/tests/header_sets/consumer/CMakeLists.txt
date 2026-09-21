# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.24)
project(OpenMSHeaderSetsConsumer CXX)
include("${TARGETS_FILE}")
add_executable(consumer main.cpp)
target_link_libraries(consumer PRIVATE OpenMS::OpenMSHeaderProbe)
