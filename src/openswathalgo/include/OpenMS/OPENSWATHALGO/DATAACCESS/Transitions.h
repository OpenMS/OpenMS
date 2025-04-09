// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/OPENSWATHALGO/OpenSwathAlgoConfig.h>
#include <string>
#include <vector>
#include <OpenMS/OPENSWATHALGO/DATAACCESS/TransitionExperiment.h>

namespace OpenSwath
{

  struct OPENSWATHALGO_DLLAPI Peptide
  {
    double rt;
    int charge;
    std::string sequence;
    std::string id;
    int getChargeState() const
    {
      return charge;
    }

    std::vector<LightModification> modifications;
    std::vector<LightTransition> transitions;
  };

  struct OPENSWATHALGO_DLLAPI Protein
  {
    std::string id;
    std::string sequence;
    std::vector<Peptide> peptides;
  };

  struct OPENSWATHALGO_DLLAPI TargetedExperiment
  {
    std::vector<Protein> proteins;
  };

}

