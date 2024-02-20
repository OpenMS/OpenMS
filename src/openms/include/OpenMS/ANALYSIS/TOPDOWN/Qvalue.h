// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHDeconvHelperStructs.h>
#include <OpenMS/KERNEL/Peak1D.h>

namespace OpenMS
{
  class PeakGroup;

  /**
@brief   Qvalue : contains functions to calculate Qvalues from deconvolution quality score
@ingroup Topdown
*/

  class OPENMS_DLLAPI Qvalue
  {
  public:
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// Calculate and perform a batch update of peak group qvalues using Qscores of target and dummy peak groups in deconvolved spectra, when FDR report is necessary.
    /// @param deconvolved_spectra target and decoy deconvolved spectra
    /// @return the noise decoy weight for decoy output
    static double updatePeakGroupQvalues(std::vector<DeconvolvedSpectrum>& deconvolved_spectra);

  private:
  };
} // namespace OpenMS