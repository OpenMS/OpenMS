# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------


add_library(GTE INTERFACE)
target_include_directories(GTE INTERFACE .)

get_filename_component(GTE_FOLDER "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(APPEND ${TPLICENSEFILE} "\nGeometric Tools (Mathematics):\n")
cat(${GTE_FOLDER}/LICENSE ${TPLICENSEFILE})
