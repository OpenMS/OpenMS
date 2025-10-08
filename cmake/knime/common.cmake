# Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
# SPDX-License-Identifier: BSD-3-Clause
#
# --------------------------------------------------------------------------
# $Maintainer: Stephan Aiche $
# $Authors: Stephan Aiche $
# --------------------------------------------------------------------------

function(check_variables variable_list)
  foreach(req IN LISTS ${variable_list})
    if(NOT DEFINED ${req})
      message(FATAL_ERROR "The variable ${req} was not set.")
    endif()
  endforeach()
endfunction()
