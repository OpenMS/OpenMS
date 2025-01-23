// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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
@brief   Qscore : quality score for PeakGroup. This class is being updated.
   For now, simply it calculate the Qscore using a fixed weight vector.
   The weight vector has been determined by logistic regression.
   But afterwards, the training part for the Qscore should be added in here.
   Or other technique such as deep learning would be used.
   This class also contains tsv output function. The tsv file contains features of PeakGroups which are used for training.
@ingroup Topdown
*/

  class OPENMS_DLLAPI Qscore
  {
  public:
    typedef FLASHDeconvHelperStructs::LogMzPeak LogMzPeak;

    /// get QScore for a peak group of specific abs_charge
    static float getQscore(const PeakGroup* pg);

  private:
    /// convert a peak group to a feature vector for Qscore calculation
    static std::vector<double> toFeatureVector_(const PeakGroup* pg);
  };
} // namespace OpenMS