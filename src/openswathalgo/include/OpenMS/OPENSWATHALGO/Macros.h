// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <cassert>

// Simple implementation of PRE and POST conditions using assert (should be on
// during debug mode and off during release mode) with a informative message
// which is printed alongside the dump.
// see http://stackoverflow.com/questions/3692954/add-custom-messages-in-assert
// "Since a pointer "is true" if it's non-null, you can use the &&-operator to
// chain and display the message".
#define OPENSWATH_PRECONDITION(condition, message)\
  assert( (condition) && (message));

#define OPENSWATH_POSTCONDITION(condition, message)\
  assert( (condition) && (message));

