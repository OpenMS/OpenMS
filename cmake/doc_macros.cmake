# Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

# clear variable to avoid accumulation
set(OPENMS_DOCUMENTATION_DIRECTORIES ""
  CACHE
  INTERNAL
  "List of paths to be searched when the API documentation is generated")


#------------------------------------------------------------------------------
# Registers the given path in the documentation system
#
# @param path_to_document Path containing header files that should be documented
macro(openms_doc_path path_to_document)
  set(OPENMS_DOCUMENTATION_DIRECTORIES
    ${path_to_document}
    ${OPENMS_DOCUMENTATION_DIRECTORIES}
    CACHE
    INTERNAL
    "List of paths to be searched when the API documentation is generated")
endmacro()
