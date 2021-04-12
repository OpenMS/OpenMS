// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
      GumbelDistributionFitResult fit(std::vector<DPosition<2> > & points);

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

