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
#ifndef OPENMS_FILTERING_TRANSFORMERS_ISOTOPEDIFFFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_ISOTOPEDIFFFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>
#include <cmath>

namespace OpenMS
{
  /**
    @brief IsotopeDiffFilter returns total intensity of peak pairs that could result from isotope peaks

        @htmlinclude OpenMS_IsotopeDiffFilter.parameters

        @ingroup SpectraFilter
  */
  class OPENMS_DLLAPI IsotopeDiffFilter :
    public FilterFunctor
  {

public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    IsotopeDiffFilter();

    /// copy constructor
    IsotopeDiffFilter(const IsotopeDiffFilter & source);

    /// destructor
    ~IsotopeDiffFilter() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    IsotopeDiffFilter & operator=(const IsotopeDiffFilter & source);
    // @}

    // @name Accessors
    // @{
    ///
    static FilterFunctor * create() { return new IsotopeDiffFilter(); }

    ///
    template <typename SpectrumType>
    double apply(SpectrumType & spectrum)
    {
      double tolerance = (double)param_.getValue("tolerance");
      double isodiff = 0;

      //iterate over all peaks
      for (Size i = 0; i < spectrum.size(); ++i)
      {
        for (Size j = 1; i + j < spectrum.size(); ++j)
        {
          double pos_ij = spectrum[i + j].getPosition()[0];
          double pos_i = spectrum[i].getPosition()[0];
          if (std::fabs(pos_ij - pos_i + 1) < tolerance)
          {
            isodiff += spectrum[i].getIntensity() + spectrum[i + j].getIntensity();
          }
          else
          {
            if (std::fabs(spectrum[i + j].getPosition()[0] - spectrum[i].getPosition()[0]) > 1 + tolerance)
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
      return "IsotopeDiffFilter";
    }

    // @}

private:
  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_ISOTOPEDIFFFILTER_H
