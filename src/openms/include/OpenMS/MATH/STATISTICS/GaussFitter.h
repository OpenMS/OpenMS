// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
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
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------
//

#ifndef OPENMS_MATH_STATISTICS_GAUSSFITTER_H
#define OPENMS_MATH_STATISTICS_GAUSSFITTER_H

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/DATASTRUCTURES/DPosition.h>

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
      struct GaussFitResult
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

#endif
