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

#ifndef OPENMS_FILTERING_TRANSFORMERS_GOODDIFFFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_GOODDIFFFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <string>
#include <cmath>

namespace OpenMS
{
  /**
    @brief GoodDiffFilter counts the number ob peak pairs whose m/z difference can be explained by a amino acid loss

        @htmlinclude OpenMS_GoodDiffFilter.parameters

        @ingroup SpectraFilter
  */
  class OPENMS_DLLAPI GoodDiffFilter :
    public FilterFunctor
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    GoodDiffFilter();

    /// copy constructor
    GoodDiffFilter(const GoodDiffFilter & source);

    /// destructor
    ~GoodDiffFilter() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    GoodDiffFilter & operator=(const GoodDiffFilter & source);
    // @}

    // @name Accessors
    // @{
    ///
    static FilterFunctor * create() { return new GoodDiffFilter(); }

    ///
    template <typename SpectrumType>
    double apply(SpectrumType & spectrum)
    {
      double tolerance = (double)param_.getValue("tolerance");
      double gooddiff = 0;
      //iterate over all peaks
      double totaldiff = 0;
      for (Size i = 0; i < spectrum.size(); ++i)
      {
        //look for each peakdifference that is in range of aa residuemasses (56/187), if it could be a aa (aamass)
        for (Size j = i; i + j < spectrum.size(); ++j)
        {
          double diff =  spectrum[i + j].getPosition()[0] - spectrum[i].getPosition()[0];
          if (diff < 56)
          {
            continue;
          }

          if (diff > 187)
          {
            j = spectrum.size();
          }
          else
          {
            totaldiff += spectrum[i + j].getIntensity() + spectrum[i].getIntensity();
            std::map<double, char>::const_iterator aait = aamass_.lower_bound(diff);
            if (aait == aamass_.end())
            {
              continue;
            }
            //look for aamasses that fit diff
            if (fabs(aait->first - diff) <= tolerance)
            {
              gooddiff += spectrum[i + j].getIntensity()  + spectrum[i].getIntensity();
            }
            else
            {
              ++aait;
              if ((aait) != aamass_.end() && fabs((aait)->first - diff) <= tolerance)
              {
                gooddiff += spectrum[i + j].getIntensity() + spectrum[i].getIntensity();
              }
            }
          }
        }
      }

      return gooddiff / totaldiff;
    }

    ///
    static const String getProductName()
    {
      return "GoodDiffFilter";
    }

    // @}


private:

    /// list of unique amino acid masses
    std::map<double, char> aamass_;
  };
}

#endif // OPENMS_FILTERING_TRANSFORMERS_GOODDIFFFILTER_H
