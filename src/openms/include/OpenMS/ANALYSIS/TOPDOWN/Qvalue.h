// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
    void static updatePeakGroupQvalues(std::vector<DeconvolvedSpectrum>& deconvolved_spectra
                                       //std::vector<DeconvolvedSpectrum>& deconvolved_decoy_spectra,
                                       //bool in_features
                                       );

  private:
    /// get a bin number given qvalue. qvalue is calculated per bin (bin from 0 to 1).
    static uint getBinNumber(double qscore, uint total_bin_number);
    /// get the qvalue corresponding to a bin number
    static double getBinValue(uint bin_number, uint total_bin_number);
    /// get the setQscore distribution
    static std::vector<double> getDistribution(const std::vector<double>& qscores, uint bin_number);
    /// get the weights of different dummy types.
    static std::vector<double> getDistributionWeights(const std::vector<double>& mixed_dist, const std::vector<std::vector<double>>& comp_dists, uint num_iterations = 100);
  };
} // namespace OpenMS