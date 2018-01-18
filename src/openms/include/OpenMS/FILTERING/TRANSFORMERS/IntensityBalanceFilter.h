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
#ifndef OPENMS_FILTERING_TRANSFORMERS_INTENSITYBALANCEFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_INTENSITYBALANCEFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <utility>

namespace OpenMS
{
  /**
      @brief IntensityBalanceFilter divides the m/z-range into ten regions and sums the
             intensity in these regions.

      The result is the intensity of the two bins with the highest intensity minus the intensity of the seven bins with lowest intensity.

        @ingroup SpectraFilter
  */
  class OPENMS_DLLAPI IntensityBalanceFilter :
    public FilterFunctor
  {

public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    IntensityBalanceFilter();

    /// copy constructor
    IntensityBalanceFilter(const IntensityBalanceFilter & source);

    /// destructor
    ~IntensityBalanceFilter() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    IntensityBalanceFilter & operator=(const IntensityBalanceFilter & source);
    // @}

    // @name Accessors
    // @{
    ///
    static FilterFunctor * create() { return new IntensityBalanceFilter(); }

    ///
    template <typename SpectrumType>
    double apply(SpectrumType & spectrum)
    {
      double bands = 10;
      std::multimap<double, Size> band_intensity;
      double parentmass = 0.0;
      if (!spectrum.getPrecursors().empty()) parentmass = spectrum.getPrecursors()[0].getMZ();
      Size j = 0;
      for (Size i = 0; i < bands; ++i)
      {
        double intensity = 0;

        //bern 2004 says to only check between 300 and size
        //but that seems inappropriate for small peptides (smallest is ca 450)
        while (j < spectrum.size() && spectrum[j].getPosition()[0] < (parentmass - 300) / bands * (i + 1) + 300)
        {
          intensity += spectrum[j++].getIntensity();
        }
        band_intensity.insert(std::make_pair(intensity, i));
      }
      j = 0;
      double total_intensity = 0;
      double twobiggest = 0;
      double sevensmallest = 0;
      for (std::multimap<double, Size>::reverse_iterator mmrit = band_intensity.rbegin(); mmrit != band_intensity.rend(); ++mmrit, ++j)
      {
        total_intensity += mmrit->first;
        //take the two biggest
        if (j < 2)
        {
          twobiggest += mmrit->first;
        }
        //take the seven smallest
        if (j > 2)
        {
          sevensmallest += mmrit->first;
        }
      }

      return (twobiggest - sevensmallest) / total_intensity;
    }

    ///
    static const String getProductName()
    {
      return "IntensityBalanceFilter";
    }

    // @}

  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_INTENSITYBALANCEFILTER_H
