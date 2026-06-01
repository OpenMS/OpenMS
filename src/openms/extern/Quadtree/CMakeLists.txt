# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# 
# --------------------------------------------------------------------------
# $Maintainer: Julianus Pfeuffer $
# $Authors: Julianus Pfeuffer $
# --------------------------------------------------------------------------


add_library(Quadtree INTERFACE)
target_include_directories(Quadtree INTERFACE "include")

get_filename_component(QUADTREE_FOLDER "${CMAKE_CURRENT_LIST_FILE}" PATH)
file(APPEND ${TPLICENSEFILE} "\nQuadtree:\n")
cat(${QUADTREE_FOLDER}/LICENSE ${TPLICENSEFILE})
