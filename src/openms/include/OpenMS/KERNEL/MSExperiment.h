// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm, Tom Waschischeck $
// --------------------------------------------------------------------------

#pragma once
/*!
    @note just to resolve issue 
          for the files where msexperiment is used
*/
#include <OpenMS/KERNEL/MSRun.h>

namespace OpenMS{
    /// @brief MSExperiment is deprecated, use MSRun instead
    using MSExperiment = MSRun;
}
