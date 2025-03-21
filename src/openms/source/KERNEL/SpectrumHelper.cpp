// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest g$
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/KERNEL/SpectrumHelper.h>

namespace OpenMS
{

  void copySpectrumMeta(const MSSpectrum & input, MSSpectrum & output, bool clear_spectrum)
  {
    // clear old spectrum first before copying
    if (clear_spectrum) output.clear(true);

    // copy the spectrum meta data
    output.SpectrumSettings::operator=(input);
    output.setRT(input.getRT());
    output.setDriftTime(input.getDriftTime());
    output.setDriftTimeUnit(input.getDriftTimeUnit());
    output.setMSLevel(input.getMSLevel());
    output.setName(input.getName());
  }
}

