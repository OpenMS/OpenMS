// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once
/*!
    @note just to resolve issue 
          for the files where OnDiscMSExperiment is used
*/
#include <OpenMS/KERNEL/OnDiscMSExperiment.h>

namespace OpenMS{
    /// @brief OnDiscMSExperiment is deprecated, use OnDiscMSRun instead
    using OnDiscMSExperiment = OnDiscMSRun;
}