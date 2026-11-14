# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Timo Sachsenberg $
# $Authors: Timo Sachsenberg $
# --------------------------------------------------------------------------

# Copies the parts of share/OpenMS the WebAssembly module packages (see CMakeLists.txt).
# Run with: cmake -DSOURCE=<share/OpenMS> -DDESTINATION=<dir> -DDIRECTORIES=<a;b> -DEXCLUDE=<a/x;b/y> -P stage_share.cmake

file(REMOVE_RECURSE "${DESTINATION}")
foreach(dir IN LISTS DIRECTORIES)
  file(COPY "${SOURCE}/${dir}" DESTINATION "${DESTINATION}"
       PATTERN ".*" EXCLUDE)
endforeach()
foreach(file IN LISTS EXCLUDE)
  file(REMOVE "${DESTINATION}/${file}")
endforeach()
