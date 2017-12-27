// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
#ifndef OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSDIFFFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSDIFFFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <cmath>

namespace OpenMS
{
  /**
    @brief NeutralLossDiffFilter returns the total intensity ob peak pairs whose m/z difference can be explained by a neutral loss

        @htmlinclude OpenMS_NeutralLossDiffFilter.parameters

        @ingroup SpectraFilter
  */
  class OPENMS_DLLAPI NeutralLossDiffFilter :
    public FilterFunctor
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    NeutralLossDiffFilter();

    /// copy constructor
    NeutralLossDiffFilter(const NeutralLossDiffFilter & source);

    /// destructor
    ~NeutralLossDiffFilter() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    NeutralLossDiffFilter & operator=(const NeutralLossDiffFilter & source);
    // @}

    // @name Accessors
    // @{
    ///
    static FilterFunctor * create() { return new NeutralLossDiffFilter(); }

    ///
    template <typename SpectrumType>
    double apply(SpectrumType & spectrum)
    {
      double tolerance = (double)param_.getValue("tolerance");
      double isodiff = 0;
      //iterate over all peaks
      for (int i = 0; i < (int)spectrum.size(); ++i)
      {
        for (int j = 1; i - j >= 0; ++j)
        {
          double pos_diff = std::fabs(spectrum[i - j].getPosition()[0] - spectrum[i].getPosition()[0]);
          if (std::fabs(pos_diff - 18) < tolerance || std::fabs(pos_diff - 17) < tolerance)   // water and ammonium
          {
            isodiff += spectrum[i - j].getIntensity() + spectrum[i].getIntensity();
          }
          else
          {
            if (pos_diff > 18 + tolerance)
            {
              break;
            }
          }
        }
      }

      return isodiff;
    }

    ///
    static const String getProductName()
    {
      return "NeutralLossDiffFilter";
    }

    // @}

  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSDIFFFILTER_H
