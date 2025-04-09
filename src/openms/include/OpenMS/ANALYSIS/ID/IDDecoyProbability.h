// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Andreas Bertsch$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/MATH/STATISTICS/GammaDistributionFitter.h>
#include <OpenMS/MATH/STATISTICS/GaussFitter.h>

#include <vector>

namespace OpenMS
{
  /**
    @brief IDDecoyProbability calculates probabilities using decoy approach

    This class calculates the probabilities using a target decoy approach. Like
    in peptide prophet the forward distribution is modeled using a gaussian
    distribution the reverse scores are modeled using a gamma distribution.

    @htmlinclude OpenMS_IDDecoyProbability.parameters

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI IDDecoyProbability :
    public DefaultParamHandler
  {
public:

    /// Default constructor
    IDDecoyProbability();

    /// Copy constructor
    IDDecoyProbability(const IDDecoyProbability & rhs);

    /// Destructor
    ~IDDecoyProbability() override;

    /// assignment operator
    IDDecoyProbability & operator=(const IDDecoyProbability & rhs);

    /**	@brief Converts the forward and reverse identification into probabilities

              @param prob_ids Output of the algorithm which includes identifications with probability based scores
              @param fwd_ids Input parameter which represents the identifications of the forward search
              @param rev_ids Input parameter which represents the identifications of the reversed search
      */
    void apply(std::vector<PeptideIdentification> & prob_ids,
               const std::vector<PeptideIdentification> & fwd_ids,
               const std::vector<PeptideIdentification> & rev_ids);

    void apply(std::vector<PeptideIdentification> & ids);

protected:

    /** @brief struct to be used to store a transformation (used for fitting)


    */
    struct Transformation_
    {
      Transformation_() :
        max_intensity(0),
        diff_score(0),
        min_score(0),
        max_score(0),
        max_intensity_bin(0)
      {
      }

      Transformation_(const Transformation_ & rhs) :
        max_intensity(rhs.max_intensity),
        diff_score(rhs.diff_score),
        min_score(rhs.min_score),
        max_score(rhs.max_score),
        max_intensity_bin(rhs.max_intensity_bin)
      {
      }

      Transformation_ & operator=(const Transformation_ & rhs)
      {
        if (this != &rhs)
        {
          max_intensity = rhs.max_intensity;
          diff_score = rhs.diff_score;
          min_score = rhs.min_score;
          max_score = rhs.max_score;
          max_intensity_bin = rhs.max_intensity_bin;
        }
        return *this;
      }

      double max_intensity;
      double diff_score;
      double min_score;
      double max_score;
      Size max_intensity_bin;
    };

    // normalizes histograms
    void normalizeBins_(const std::vector<double> & scores, std::vector<double> & binned, Transformation_ & trafo);

    // returns the probability of given score with the transformations of reverse and forward searches and the results of the fits
    double getProbability_(const Math::GammaDistributionFitter::GammaDistributionFitResult & result_gamma,
                               const Transformation_ & gamma_trafo,
                               const Math::GaussFitter::GaussFitResult & result_gauss,
                               const Transformation_ & gauss_trafo,
                               double score);


    void generateDistributionImage_(const std::vector<double> & ids, const String & formula, const String & filename);

    void generateDistributionImage_(const std::vector<double> & all_ids, const Transformation_ & all_trans, const String & fwd_formula, const String & rev_formula, const String & filename);


    void apply_(std::vector<PeptideIdentification> & ids, const std::vector<double> & rev_scores, const std::vector<double> & fwd_scores, const std::vector<double> & all_scores);

  };

} // namespace OpenMS

