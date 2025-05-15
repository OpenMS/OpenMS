// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Witold Wolski $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>
#include <map>
#include <vector>
#include <string>

namespace OpenSwath
{

  struct LightCompound;
  struct LightTargetedExperiment;
  struct LightTransition;

  struct OPENSWATHALGO_DLLAPI TransitionHelper
  {

    static void convert(LightTargetedExperiment& lte,
                        std::map<std::string,
                        std::vector<LightTransition> >& transmap);


    // TODO : remove and explain German comments
    // spiegel
    static bool findPeptide(const LightTargetedExperiment& lte,
                            const std::string& peptideRef,
                            LightCompound& pep);
  };

} //end namespace

