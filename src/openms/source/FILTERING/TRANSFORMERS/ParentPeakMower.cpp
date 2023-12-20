// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//

#include <OpenMS/FILTERING/TRANSFORMERS/ParentPeakMower.h>

using namespace std;

namespace OpenMS
{
  ParentPeakMower::ParentPeakMower() :
    DefaultParamHandler("ParentPeakMower")
  {
    defaults_.setValue("window_size", 2.0, "The size of the m/z window where the peaks are removed, +/- window_size.");
    defaults_.setValue("default_charge", 2, "If the precursor has no charge set, the default charge is assumed.");
    defaults_.setValue("clean_all_charge_states", 1, "Set to 1 if precursor ions of all possible charge states should be removed.", {"advanced"});
    defaults_.setValue("consider_NH3_loss", 1, "Whether NH3 loss peaks from the precursor should be removed.");
    defaults_.setValue("consider_H2O_loss", 1, "Whether H2O loss peaks from the precursor should be removed.");
    defaults_.setValue("reduce_by_factor", 0, "Reduce the intensities of the precursor and related ions by a given factor (set 'set_to_zero' to 0).", {"advanced"});
    defaults_.setValue("factor", 1000.0, "Factor which is used to reduce the intensities if 'reduce_by_factor' is selected.", {"advanced"});
    defaults_.setValue("set_to_zero", 1, "Reduce the intensities of the precursor and related ions to zero.", {"advanced"});
    defaultsToParam_();
  }

  ParentPeakMower::~ParentPeakMower() = default;

  ParentPeakMower::ParentPeakMower(const ParentPeakMower & source) :
    DefaultParamHandler(source)
  {
  }

  ParentPeakMower & ParentPeakMower::operator=(const ParentPeakMower & source)
  {
    if (this != &source)
    {
      DefaultParamHandler::operator=(source);
    }
    return *this;
  }

  void ParentPeakMower::filterPeakSpectrum(PeakSpectrum & spectrum)
  {
    filterSpectrum(spectrum);
  }

  void ParentPeakMower::filterPeakMap(PeakMap & exp)
  {
    for (PeakMap::Iterator it = exp.begin(); it != exp.end(); ++it)
    {
      filterSpectrum(*it);
    }
  }

}
