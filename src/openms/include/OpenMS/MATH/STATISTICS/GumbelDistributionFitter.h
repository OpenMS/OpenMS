// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------
//
#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

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
    class OPENMS_DLLAPI GumbelDistributionFitter
    {
public:

      /// struct to represent the parameters of a gumbel distribution
      struct GumbelDistributionFitResult
      {
        GumbelDistributionFitResult(double local_a = 1.0, double local_b = 2.0) :
          a(local_a),
          b(local_b)
        {
        }

        /// location parameter a
        double a;
        /// scale parameter b
        double b;

        double eval(double x) const;
        double log_eval_no_normalize(double x) const ;
      };

      /// Default constructor
      GumbelDistributionFitter();
      /// Destructor
      virtual ~GumbelDistributionFitter();

      /// sets the gumbel distribution start parameters a and b for the fitting
      void setInitialParameters(const GumbelDistributionFitResult & result);

      /**
          @brief Fits a gumbel distribution to the given data points

          @param points Input parameter which represents the point used for the fitting

          @exception Exception::UnableToFit is thrown if fitting cannot be performed
      */
      GumbelDistributionFitResult fit(std::vector<DPosition<2> > & points) const;

      /**
          @brief Fits a gumbel distribution to the given data x values. Fills a
          weighted histogram first and generates y values.

          @param x Input x values
          @param w Input weights

          @exception Exception::UnableToFit is thrown if fitting cannot be performed
      */
      GumbelDistributionFitResult fitWeighted(const std::vector<double> & x, const std::vector<double> & w);

protected:

      GumbelDistributionFitResult init_param_;

private:
      /// Copy constructor (not implemented)
      GumbelDistributionFitter(const GumbelDistributionFitter & rhs);
      /// assignment operator (not implemented)
      GumbelDistributionFitter & operator=(const GumbelDistributionFitter & rhs);
    };
  }
}

