// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Kyowon Jeong $
// $Authors: Kyowon Jeong $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/ANALYSIS/TOPDOWN/DeconvolvedSpectrum.h>
#include <OpenMS/ANALYSIS/TOPDOWN/FLASHHelperClasses.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/METADATA/Precursor.h>


namespace OpenMS
{
  class PeakGroup;

  /**
@brief scoring functions for PeakGroup.
   For now, only Qscore has been implemented. For Qscore,
   the weight vector has been determined by logistic regression.
   In the future, other technique such as deep learning would be used.
@ingroup Topdown
*/

  class OPENMS_DLLAPI PeakGroupScoring
  {
  public:
    typedef FLASHHelperClasses::LogMzPeak LogMzPeak;

    /// get the number of features used for Qsocre
    static int getQscoreFeatureCount();

    /// get QScore for a peak group of specific abs_charge
    static double getQscore(const PeakGroup* pg);

    /// Write Csv file for Qscore training.
    static void writeAttCsvForQscoreTraining(const DeconvolvedSpectrum& deconvolved_spectrum, std::fstream& f);

    /// Write Csv file for Qscore training header.
    static void writeAttCsvForQscoreTrainingHeader(std::fstream& f);

    /// get Deep learning based peak group score. Not implemented yet.
    static double getDLscore(PeakGroup* pg, const MSSpectrum& spec, const FLASHHelperClasses::PrecalculatedAveragine& avg, double tol);

    /// convert a peak group to a feature vector for setQscore calculation
    static std::vector<double> toFeatureVector(const PeakGroup* pg);

  private:
    /// the weights for Qscore calculation
    static std::vector<double> weight_;

    /// charge and isotope counts for DL scoring.
    static const int charge_count_for_DL_scoring_ = 11;
    static const int iso_count_for_DL_scoring_ = 13;
  };
} // namespace OpenMS