// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/MATH/STATISTICS/CumulativeBinomial.h>
#include <boost/math/special_functions/binomial.hpp>

#include <iostream>
#include <numeric>

namespace OpenMS
{

  namespace Math
  {

    double CumulativeBinomial::compute(Size n, Size k, double p)
    {
      double p_cumul = 0.0;
      if (p < 1e-99) return static_cast<double>(k == 0);
      if (1 - p < 1e-99) return static_cast<double>(k != n);
      if (k > n)  return 1.0;

      for (Size j = 0; j < k; ++j)
      {
        double coeff = 0;

        try
        {
          coeff = boost::math::binomial_coefficient<double>(static_cast<unsigned int>(n), static_cast<unsigned int>(j));
        }
        catch (std::overflow_error const& /*e*/)
        {
          // not sure if a warning is appropriate here, since if it happens, it will happen very often for the same spectrum and flood the stdout
//          std::cout << "Warning: Binomial coefficient for match-odds score has overflowed! Setting value to the maximal double value." << std::endl;
//          std::cout << "binomial_coefficient was called with N = " << n << " and k = " << j << std::endl;
          coeff = std::numeric_limits<double>::max();
        }

        p_cumul += coeff * pow(p,  static_cast<int>(j)) * pow((1-p), static_cast<int>((n-j)));
      }

      // A result of p_cumul >= 1.0 does not make sense theoretically, but it might reach 1.0 because of insufficient precision,
      // solved by using largest value smaller than 1.0
      if (p_cumul >= 1.0)
      {
        // TODO: C11 change to nexttoward
//        p_cumul = nexttoward(1.0, 0.0);
        p_cumul = 1.0 - std::numeric_limits<double>::epsilon();
      }

      return p_cumul;
    }
  } // namespace Math

} // namespace OpenMS
