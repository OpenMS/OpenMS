// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Andreas Bertsch $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_MATH_STATISTICS_GAMMADISTRIBUTIONFITTER_H
#define OPENMS_MATH_STATISTICS_GAMMADISTRIBUTIONFITTER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

#include <vector>

// gsl includes
#include <gsl/gsl_rng.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_multifit_nlin.h>


namespace OpenMS
{
  namespace Math
  {
    /**
      @brief Implements a fitter for the Gamma distribution.

      This class fits a Gamma distribution to a number of data points.
      The results as well as the initial guess are specified using the struct
          GammaDistributionFitResult.

      The formula with the fitted parameters can be transformed into a
      gnuplot formula using getGnuplotFormulai() after fitting.

          The implementation is done using GSL fitting algorithms.

          @ingroup Math
      */
    class OPENMS_DLLAPI GammaDistributionFitter
    {
public:

      /// struct to represent the parameters of a gamma distribution
      struct GammaDistributionFitResult
      {
public:

        GammaDistributionFitResult() :
          b(1.0),
          p(5.0)
        {
        }

        GammaDistributionFitResult(const GammaDistributionFitResult & rhs) :
          b(rhs.b),
          p(rhs.p)
        {
        }

        GammaDistributionFitResult & operator=(const GammaDistributionFitResult & rhs)
        {
          if (this != &rhs)
          {
            b = rhs.b;
            p = rhs.p;
          }
          return *this;
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
      void setInitialParameters(const GammaDistributionFitResult & result);

      /**
          @brief Fits a gamma distribution to the given data points

          @param points Input parameter which represents the point used for the fitting

          @exception Exception::UnableToFit is thrown if fitting cannot be performed
      */
      GammaDistributionFitResult fit(std::vector<DPosition<2> > & points);

      /// returns the gnuplot formula of the fitted gamma distribution
      const String & getGnuplotFormula() const;

protected:

      static int gammaDistributionFitterf_(const gsl_vector * x, void * params, gsl_vector * f);

      static int gammaDistributionFitterdf_(const gsl_vector * x, void * params, gsl_matrix * J);

      static int gammaDistributionFitterfdf_(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J);

      void printState_(size_t iter, gsl_multifit_fdfsolver * s);

      GammaDistributionFitResult init_param_;

      String gnuplot_formula_;

private:
      /// Copy constructor (not implemented)
      GammaDistributionFitter(const GammaDistributionFitter & rhs);
      /// assignment operator (not implemented)
      GammaDistributionFitter & operator=(const GammaDistributionFitter & rhs);
    };
  }
}

#endif
