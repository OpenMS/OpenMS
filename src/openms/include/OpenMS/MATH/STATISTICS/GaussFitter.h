// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

#pragma once

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <cmath>
#include <vector>

namespace OpenMS
{
  namespace Math
  {
    /**
        @brief Implements a fitter for Gaussian functions

        This class fits a Gaussian distribution to a number of data points.
        The results as well as the initial guess are specified using the struct GaussFitResult.

        The complete Gaussian formula with the fitted parameters can be transformed into a
        gnuplot formula using getGnuplotFormula after fitting.

        @ingroup Math
    */
    class OPENMS_DLLAPI GaussFitter
    {
public:

      /// struct of parameters of a Gaussian distribution
      struct OPENMS_DLLAPI GaussFitResult
      {
public:
        GaussFitResult()
        : A(-1.0), x0(-1.0), sigma(-1.0) {}
        GaussFitResult(double a, double x, double s)
        : A(a), x0(x), sigma(s) {}

        /// parameter A of Gaussian distribution (amplitude)
        double A;

        /// parameter x0 of Gaussian distribution (center position)
        double x0;

        /// parameter sigma of Gaussian distribution (width)
        double sigma;


        /**
          @brief Evaluate the current density Gaussian model at the specified point.

          Returns the intensities (i.e. probabilities scaled by the factor 'A') of the PDF at the given positions.
          This function can be called with any set of parameters, e.g. the initial parameters (to get a 'before-fit' status),
          or after fitting.
        */
        double eval(double x) const;

        /**
          @brief Evaluate the current log density of the Gaussian model at the specified point.

          Returns the intensities (i.e. probabilities scaled by the factor 'A') of the PDF at the given positions.
          This function can be called with any set of parameters, e.g. the initial parameters (to get a 'before-fit' status),
          or after fitting.
        */
        double log_eval_no_normalize(double x) const;

      private:
        double halflogtwopi = 0.5*log(2.0*Constants::PI);
      };

      /// Constructor
      GaussFitter();

      /// Destructor
      virtual ~GaussFitter();

      /// sets the initial parameters used by the fit method as initial guess for the Gaussian
      void setInitialParameters(const GaussFitResult& result);

      /**
          @brief Fits a Gaussian distribution to the given data points

          @param points the data points used for the Gaussian fitting

          @exception Exception::UnableToFit is thrown if fitting cannot be performed
      */
      GaussFitResult fit(std::vector<DPosition<2> > & points) const;

      /**
        @brief Evaluate the current Gaussian model at the specified points.

        Returns the intensities (i.e. probabilities scaled by the factor 'A') of the PDF at the given positions.
        This function can be called with any set of parameters, e.g. the initial parameters (to get a 'before-fit' status),
        or after fitting.

      */
      static std::vector<double> eval(const std::vector<double>& evaluation_points, const GaussFitResult& model);

protected:

      GaussFitResult init_param_;

private:
     /// Copy constructor (not implemented)
      GaussFitter(const GaussFitter & rhs);

      /// Assignment operator (not implemented)
      GaussFitter & operator=(const GaussFitter & rhs);
    };
  }
}

