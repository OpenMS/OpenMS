# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche, Chris Bielow, Hannes Röst $
# $Authors: Hannes Röst, Stephan Aiche $
# --------------------------------------------------------------------------


project("OpenSWATHAlgo")
cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

#------------------------------------------------------------------------------
# Find Boost lib
#------------------------------------------------------------------------------
find_boost()

if(NOT Boost_FOUND)
  message(FATAL_ERROR "Boost was not found!")
endif()

#------------------------------------------------------------------------------
# get information regarding the openswath files
include(source/OPENSWATHALGO/OpenSwathAlgoFiles.cmake)

#------------------------------------------------------------------------------
# add the library
openms_add_library(TARGET_NAME OpenSwathAlgo
                   SOURCE_FILES ${OpenSwathAlgoFiles}
                   HEADER_FILES ${OpenSwathAlgoHeaders}
                   INTERNAL_INCLUDES ${PROJECT_SOURCE_DIR}/include ${PROJECT_BINARY_DIR}/include
                   EXTERNAL_INCLUDES ${Boost_INCLUDE_DIRS}
                   PRIVATE_LINK_LIBRARIES Eigen3::Eigen
                   DLL_EXPORT_PATH "OpenMS/OPENSWATHALGO/")

openms_doc_path("${PROJECT_SOURCE_DIR}/include")
