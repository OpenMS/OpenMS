// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/DATASTRUCTURES/StringUtils.h>
#include <vector>


namespace OpenMS
{
  namespace Math
  {
    /**
      @brief Implements a fitter for the Gumbel distribution.

      This class fits a Gumbel distribution to a number of data points.
      The results as well as the initial guess are specified using the struct
      GumbelDistributionFitResult.

      The formula with the fitted parameters can be transformed into a
      gnuplot formula using getGnuplotFormula() after fitting.

      @ingroup Math
      */
    class OPENMS_DLLAPI GumbelMaxLikelihoodFitter
    {

public:

      /// struct to represent the parameters of a gumbel distribution
      struct OPENMS_DLLAPI GumbelDistributionFitResult
      {
        GumbelDistributionFitResult(double local_a, double local_b) :
          a(local_a),
          b(local_b)
        {
        }

        /// location parameter a
        double a;
        /// scale parameter b
        double b;

        /// log density at @p x. Note: despite its name, the density is fully normalized.
        double log_eval_no_normalize(double x) const;

        /// (normalized) density at @p x
        double eval(double x) const;

        /// mode of the distribution (equals the location parameter a)
        double mode() const;
      };

      /// Default constructor
      GumbelMaxLikelihoodFitter();
      /// Default constructor
      GumbelMaxLikelihoodFitter(GumbelDistributionFitResult init);
      /// Destructor
      virtual ~GumbelMaxLikelihoodFitter();

      /// sets the gumbel distribution start parameters a and b for the fitting
      void setInitialParameters(const GumbelDistributionFitResult & result);

      /**
          @brief Fits a Gumbel (maximum) distribution to the data points @p x with weights @p w by weighted maximum likelihood.

          The scale parameter is the unique root of the (monotone) profile score equation, which is solved by
          safeguarded Newton iterations; the location parameter then follows in closed form.
          Data points with non-finite values or non-positive weights are ignored.
          If no data point remains, the initial parameters are returned.
          The fitted parameters are used as initial parameters of the next call.

          @param[in] x Input x values
          @param[in] w Input weights (same length as @p x)

          @exception Exception::InvalidParameter is thrown if @p x and @p w differ in size
          @exception Exception::UnableToFit is thrown if fitting cannot be performed (e.g. the data have no spread)
      */
      GumbelDistributionFitResult fitWeighted(const std::vector<double> & x, const std::vector<double> & w);

protected:

      GumbelDistributionFitResult init_param_;

private:
      /// Copy constructor (not implemented)
      GumbelMaxLikelihoodFitter(const GumbelMaxLikelihoodFitter & rhs);
      /// assignment operator (not implemented)
      GumbelMaxLikelihoodFitter & operator=(const GumbelMaxLikelihoodFitter & rhs);
    };
  }
}
