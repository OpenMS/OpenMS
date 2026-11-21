# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
# $Maintainer: Timo Sachsenberg $

# Keep the loader identity and Mach-O current version consistent.
function(openms_set_library_soversion target)
  set_target_properties(${target} PROPERTIES
    VERSION "${OPENMS_PACKAGE_VERSION}"
    SOVERSION "${OPENMS_PACKAGE_VERSION}")
endfunction()
