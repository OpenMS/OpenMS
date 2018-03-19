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
#ifndef OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTFILTER_H
#define OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTFILTER_H

#include <OpenMS/FILTERING/TRANSFORMERS/FilterFunctor.h>

#include <map>
#include <cmath>

namespace OpenMS
{
  /**
    @brief total intensity of peak pairs that could result from complementing fragments of charge state 1

        @htmlinclude OpenMS_ComplementFilter.parameters

        @ingroup SpectraFilter
  */
  class OPENMS_DLLAPI ComplementFilter :
    public FilterFunctor
  {
public:

    // @name Constructors and Destructors
    //@{
    /// standard constructor
    ComplementFilter();

    /// copy constructor
    ComplementFilter(const ComplementFilter & source);

    /// destructor
    ~ComplementFilter() override;
    //@}

    // @name Operators
    //@{
    /// assignment operator
    ComplementFilter & operator=(const ComplementFilter & source);
    //@}

    // @name Accessors
    //@{
    static FilterFunctor * create() { return new ComplementFilter(); }

    /// returns the total intensity of peak pairs which could result from complementing fragments
    template <typename SpectrumType>
    double apply(SpectrumType & spectrum)
    {
      if (spectrum.size() < 2)
      {
        return 0;
      }
      double tolerance = (double)param_.getValue("tolerance");
      double parentmass = 0.0;
      if (!spectrum.getPrecursors().empty()) parentmass = spectrum.getPrecursors()[0].getMZ();
      double result(0);

      spectrum.sortByPosition();

      /// @improvement think about an correct fast algorithm, not just an heuristic (Andreas)
      Size j = spectrum.size() - 1;
      for (Size i = 0; i < spectrum.size() && i <= j; /*++i*/)
      {
        double sum = spectrum[i].getPosition()[0] + spectrum[j].getPosition()[0];

        if (std::fabs(sum - parentmass) < tolerance)
        {
          result += spectrum[i].getIntensity() + spectrum[j].getIntensity();
        }

        if (sum < parentmass)
        {
          ++i;
        }
        else
        {
          if (sum > parentmass)
          {
            --j;
          }
        }
      }

      return result;
    }

    /// returns the name for registration at the factory
    static const String getProductName()
    {
      return "ComplementFilter";
    }

    //@}

  };
}
#endif // OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTFILTER_H
