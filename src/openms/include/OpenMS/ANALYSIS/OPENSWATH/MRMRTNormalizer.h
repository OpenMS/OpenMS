// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

#include <cstddef> // for size_t & ptrdiff_t
#include <vector>
#include <string>

namespace OpenMS
{

  /**
    @brief The MRMRTNormalizer will find retention time peptides in data.

    This tool will take a description of RT peptides and their normalized
    retention time to write out a transformation file on how to transform the
    RT space into the normalized space.

    The principle is adapted from the following publication:
    Escher, C. et al. (2012), Using iRT, a normalized retention time for more 
    targeted measurement of peptides. Proteomics, 12: 1111-1121.
  */
  class OPENMS_DLLAPI MRMRTNormalizer
  {

protected:
  
    /**
      @brief This function computes a candidate outlier peptide by iteratively
       leaving one peptide out to find the one which results in the maximum R^2
       of a first order linear regression of the remaining ones. The data points
       are submitted as two vectors of doubles (x- and y-coordinates).

      @return The position of the candidate outlier peptide as supplied by the
       vector is returned.

      @exception Exception::UnableToFit is thrown if fitting cannot be performed
    */
    static int jackknifeOutlierCandidate_(const std::vector<double>& x, const std::vector<double>& y);

    /**
      @brief This function computes a candidate outlier peptide by computing
       the residuals of all points to the linear fit and selecting the one with
       the largest deviation. The data points are submitted as two vectors of
       doubles (x- and y-coordinates).

      @return The position of the candidate outlier peptide as supplied by the
       vector is returned.

      @exception Exception::UnableToFit is thrown if fitting cannot be performed
    */
    static int residualOutlierCandidate_(const std::vector<double>& x, const std::vector<double>& y);

public:
 
    /**
      @brief This function removes potential outliers in a linear regression dataset.

       Two thresholds need to be defined, first a lower R^2 limit to accept the
       regression for the RT normalization and second, the lower limit of peptide
       coverage. The algorithms then selects candidate outlier peptides using the
       RANSAC outlier detection algorithm and returns the corrected set of peptides
       if the two thresholds are satisfied.

      @param pairs Input data (paired data of type <experimental_rt, theoretical_rt>)
      @param rsq_limit Minimal R^2 required
      @param coverage_limit Minimal coverage required (if the number of points
       falls below this fraction, the algorithm aborts)
      @param max_iterations Maximum iterations for the RANSAC algorithm
      @param max_rt_threshold Maximum deviation from fit for the retention time.
       This must be in the unit of the second dimension (e.g. theoretical_rt).
      @param sampling_size The number of data points to sample for the RANSAC algorithm.

      @return A vector of pairs is returned if the R^2 limit was reached without
       reaching the coverage limit. If the limits are reached, an exception is
       thrown.

      @exception Exception::UnableToFit is thrown if fitting cannot be
      performed (rsq_limit and coverage_limit cannot be fulfilled)
    */ 
    static std::vector<std::pair<double, double> > removeOutliersRANSAC(const std::vector<std::pair<double, double> >& pairs,
                                                                        double rsq_limit,
                                                                        double coverage_limit,
                                                                        size_t max_iterations,
                                                                        double max_rt_threshold,
                                                                        size_t sampling_size);

    /**
      @brief This function removes potential outliers in a linear regression dataset.

       Two thresholds need to be defined, first a lower R^2 limit to accept the
       regression for the RT normalization and second, the lower limit of peptide
       coverage. The algorithms then selects candidate outlier peptides and applies
       the Chauvenet's criterion on the assumption that the residuals are normal
       distributed to determine whether the peptides can be removed. This is done
       iteratively until both limits are reached.

      @param pairs Input data (paired data of type <experimental_rt, theoretical_rt>)
      @param rsq_limit Minimal R^2 required
      @param coverage_limit Minimal coverage required (the number of points
      falls below this fraction, the algorithm aborts)
      @param use_chauvenet Whether to only remove outliers that fulfill
      Chauvenet's criterion for outliers (otherwise it will remove any outlier
      candidate regardless of the criterion)
      @param method Outlier detection method ("iter_jackknife" or "iter_residual")

      @return A vector of pairs is returned if the R^2 limit was reached without
       reaching the coverage limit. If the limits are reached, an exception is
       thrown.

      @exception Exception::UnableToFit is thrown if fitting cannot be
      performed (rsq_limit and coverage_limit cannot be fulfilled)
    */
    static std::vector<std::pair<double, double> > removeOutliersIterative(const std::vector<std::pair<double, double> >& pairs,
                                                                           double rsq_limit, 
                                                                           double coverage_limit, 
                                                                           bool use_chauvenet,
                                                                           const std::string& method);

    /**
      @brief This function computes Chauvenet's criterion probability for a vector
       and a value whose position is submitted.

      @return Chauvenet's criterion probability
    */
    static double chauvenet_probability(const std::vector<double>& residuals, int pos);

    /**
      @brief This function computes Chauvenet's criterion for a vector and a value
       whose position is submitted.

      @return TRUE, if Chauvenet's criterion is fulfilled and the outlier can be removed.
    */
    static bool chauvenet(const std::vector<double>& residuals, int pos);

    /**
      * @brief Computes coverage of the RT normalization peptides over the whole RT range, ensuring that each bin has enough peptides
      *
      * @param rtRange The (estimated) full RT range in iRT space (theoretical RT)
      * @param pairs The RT normalization peptide pairs (pair = experimental RT / theoretical RT)
      * @param nrBins The number of bins to be used
      * @param minPeptidesPerBin The minimal number of peptides per bin to be used to be considered full
      * @param minBinsFilled The minimal number of bins needed to be full
      *
      * @return Whether more than the minimal number of bins are covered
      *
    */
    static bool computeBinnedCoverage(const std::pair<double,double> & rtRange, 
                                      const std::vector<std::pair<double, double> > & pairs,
                                      int nrBins, 
                                      int minPeptidesPerBin,
                                      int minBinsFilled);

  };

}

