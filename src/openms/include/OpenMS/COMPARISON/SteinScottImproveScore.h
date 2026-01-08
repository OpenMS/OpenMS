// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Vipul Patel $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/COMPARISON/PeakSpectrumCompareFunctor.h>
#include <cmath>

namespace OpenMS
{
  /**
      @brief Similarity score based of Stein & Scott

      This is a pairwise based score function. The spectrum contains peaks, and
      each peak can be defined by two values (mz and the intensity).  The score
      function takes the sum of the product of the peak intensities from
      Spectrum 1 and Spectrum 2, only if the mz-ratio distance between the two
      spectra is smaller than a given window size. In the default status, the
      window size is (accuracy of the mass spectrometer). This sum is
      normalised by dividing it with a distance function. sqrt(sum of the
      Intensity of square Spectrum1 sum of the Intensity of square Spectrum2).
      This is all based on SteinScott score.  To distinguish the close from the
      distant spectra an additional term is calculated. It denotes the expected
      value of both Spectra under the random placement of all peaks, within
      the given mass-to-charge range. The probability that two peaks with
      randomized intensity values lie within two epsilon of each other is a
      constant. This constant is proportional to epsilon. So the additional
      term is the sum over all peaks of Spectrum 1 and Spectrum 2 of the
      products of their intensities multiplied with the constant.


      The details of the score can be found in:
      Signal Maps for Mass Spectrometry-based
      Comparative Proteomics
      Amol Prakash, Parag Mallick , Jeffrey Whiteaker, Heidi Zhang,
      Amanda Paulovich, Mark Flory, Hookeun Lee, Ruedi Aebersold,
      and Benno Schwikowski

      @htmlinclude OpenMS_SteinScottImproveScore.parameters

      @ingroup SpectraComparison
  */
  class OPENMS_DLLAPI SteinScottImproveScore :
    public PeakSpectrumCompareFunctor
  {

public:
    /// default constructor
    SteinScottImproveScore();
    /// copy constructor
    SteinScottImproveScore(const SteinScottImproveScore & source);
    /// destructor
    ~SteinScottImproveScore() override;
    /// assignment operator
    SteinScottImproveScore & operator=(const SteinScottImproveScore & source);
    /**
        @brief Similarity pairwise score

        This function return the similarity score of two Spectra based on SteinScott.
    */
    double operator()(const PeakSpectrum & spec1, const PeakSpectrum & spec2) const override;
    /**
        @brief Similarity pairwise score itself

        This function return the similarity score of itself based on SteinScott.
    */
    double operator()(const PeakSpectrum & spec) const override;
  };
}


