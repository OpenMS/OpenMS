// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/OPENSWATHALGO/DATAACCESS/ISpectrumAccess.h>

#include <OpenMS/KERNEL/MSSpectrum.h>

namespace OpenMS
{
  /**
  @brief The SpectrumAddition is used to add up individual spectra 
  
  It uses the given sampling rate to resample the spectra in m/z domain and
  then add them up. This may lead to a certain inaccuracy, especially if a
  inappropriate resampling rate is chosen.

  */
  class OPENMS_DLLAPI SpectrumAddition
  {

public:

    /// adds up a list of Spectra by resampling them and then addition of intensities
    static OpenSwath::SpectrumPtr addUpSpectra(const std::vector<OpenSwath::SpectrumPtr>& all_spectra,
                                               double sampling_rate,
                                               bool filter_zeros);

    /// adds up a list of Spectra by resampling them and then addition of intensities
    static OpenMS::MSSpectrum addUpSpectra(const std::vector< OpenMS::MSSpectrum>& all_spectra,
                                           double sampling_rate,
                                           bool filter_zeros);

  };
}


