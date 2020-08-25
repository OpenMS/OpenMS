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
      GammaDistributionFitResult fit(const std::vector<DPosition<2> >& points);

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

