# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

cmake_minimum_required(VERSION 3.15 FATAL_ERROR)

# include helper functions 
include ( ${SCRIPT_DIR}common.cmake )

set(required_variables "PLATFORM;TARGET_DIR;TEMPLATE_FOLDER")
check_variables(required_variables)

configure_file(${TEMPLATE_FOLDER}/binaries_${PLATFORM}.ini ${TARGET_DIR}/binaries.ini COPYONLY)
