// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
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
      @brief Implements a fitter for the Gamma distribution.

      This class fits a Gamma distribution to a number of data points.
      The results as well as the initial guess are specified using the struct
      GammaDistributionFitResult.
     
      @note We actually fit a slightly customized version of the gamma distribution
      that is 0.0 if the parameters b or p are <= 0.0. With this modification we 
      can still use an unconstrained optimization algorithm.

      @ingroup Math
    */
    class OPENMS_DLLAPI GammaDistributionFitter
    {
public:

      /// struct to represent the parameters of a gamma distribution
      struct GammaDistributionFitResult
      {
public:

        GammaDistributionFitResult(double bIn, double pIn) :
          b(bIn),
          p(pIn)
        {
        }

        /// parameter b of the gamma distribution
        double b;

        /// parameter p of the gamma distribution
        double p;
      };

      /// Default constructor
      GammaDistributionFitter();
      /// Destructor
      virtual ~GammaDistributionFitter();

      /// sets the gamma distribution start parameters b and p for the fitting
      void setInitialParameters(const GammaDistributionFitResult& result);

      /**
          @brief Fits a gamma distribution to the given data points

          @param points Input parameter which represents the point used for the fitting

          @exception Exception::UnableToFit is thrown if fitting cannot be performed
      */
      GammaDistributionFitResult fit(const std::vector<DPosition<2> >& points) const;

protected:

      GammaDistributionFitResult init_param_;

private:
      /// Copy constructor (not implemented to prevent usage)
      GammaDistributionFitter(const GammaDistributionFitter& rhs);
      /// assignment operator (not implemented to prevent usage)
      GammaDistributionFitter& operator=(const GammaDistributionFitter& rhs);
    };
  }
}

