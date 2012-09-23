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

#ifndef OPENMS_MATH_STATISTICS_GAUSSFITTER_H
#define OPENMS_MATH_STATISTICS_GAUSSFITTER_H

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
        @brief Implements a fitter for gaussian functions

        This class fits a gaussian distribution to a number of data points.
        The results as well as the initial guess are specified using the struct GaussFitResult.

        The complete gaussian formula with the fitted parameters can be transformed into a
        gnuplot formula using getGnuplotFormula after fitting.

        The fitting is implemented using GSL fitting algorithms.

        @ingroup Math
    */
    class OPENMS_DLLAPI GaussFitter
    {
public:

      /// struct of parameters of a gaussian distribution
      struct GaussFitResult
      {
public:

        /// parameter A of gaussian distribution (amplitude)
        double A;

        /// parameter x0 of gaussian distribution (left/right shift)
        double x0;

        /// parameter sigma of gaussian distribution (width)
        double sigma;
      };

      /// Default constructor
      GaussFitter();

      /// Destructor
      virtual ~GaussFitter();

      /// sets the initial parameters used by the fit method as inital guess for the gaussian
      void setInitialParameters(const GaussFitResult & result);

      /**
          @brief Fits a gaussian distribution to the given data points

          @param points the data points used for the gaussian fitting

          @exception Exception::UnableToFit is thrown if fitting cannot be performed
      */
      GaussFitResult fit(std::vector<DPosition<2> > & points);

      /// return the gnuplot formula of the gaussian
      const String & getGnuplotFormula() const;

protected:

      static int gaussFitterf_(const gsl_vector * x, void * params, gsl_vector * f);

      static int gaussFitterdf_(const gsl_vector * x, void * params, gsl_matrix * J);

      static int gaussFitterfdf_(const gsl_vector * x, void * params, gsl_vector * f, gsl_matrix * J);

      void printState_(size_t iter, gsl_multifit_fdfsolver * s);

      GaussFitResult init_param_;

      String gnuplot_formula_;

private:

      /// Copy constructor (not implemented)
      GaussFitter(const GaussFitter & rhs);

      /// Assignment operator (not implemented)
      GaussFitter & operator=(const GaussFitter & rhs);
    };
  }
}

#endif
