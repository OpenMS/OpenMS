// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#pragma once

namespace OpenDIA
{
  /**
    @brief Run the tool-private OpenDIA command-line workflow.

    @param[in] argc Number of command-line arguments.
    @param[in] argv Command-line arguments, including the executable name.
    @return The TOPP exit code.
  */
  int run(int argc, const char** argv);
}
