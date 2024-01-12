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
#include <OpenMS/METADATA/Precursor.h>


namespace OpenMS
{
  class PeakGroup;

  /**
@brief   setQscore : quality score for PeakGroup. This class is being updated.
   For now, simply it calculate the setQscore using a fixed weight vector.
   The weight vector has been determined by logistic regression.
   But afterwards, the training part for the setQscore should be added in here.
   Or other technique such as deep learning would be used.
   This class also contains tsv output function. The tsv file contains features of PeakGroups which are used for training.
@ingroup Topdown
*/

  class OPENMS_DLLAPI Qscore
  {
  public:
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// get QScore for a peak group of specific abs_charge
    static double getQscore(const PeakGroup* pg, bool is_profile = false, double cv = -1);

    static void writeAttCsvFromDummy(const DeconvolvedSpectrum& deconvolved_spectrum, std::fstream& f);

    static void writeAttCsvFromDummyHeader(std::fstream& f);

  private:
    /// convert a peak group to a feature vector for setQscore calculation


    static std::vector<double> toFeatureVector_(const PeakGroup* pg);

    static std::vector<double> weight_centroid_;
    static std::vector<double> weight_profile_;
    static std::vector<double> weight_CV_;
    // the weights for per cosine, SNR, PPM error, charge score, and intercept.
    // ======================================
    // IsotopeCosine                  40.7425
    // ChargeCosine                      0.205
    // MassSNR1                        -0.1984
    // ChargeSNR1                        0.213
    // Intercept                      -40.3701
    //
    // IsotopeCosine    -55.8387
    // ChargeCosine        0.0253
    // MassSNR1            0.2473
    // ChargeSNR2         -0.6765
    // Intercept          55.8594

  };
} // namespace OpenMS
