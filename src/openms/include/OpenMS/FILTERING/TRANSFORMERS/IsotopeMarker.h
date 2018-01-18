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
#ifndef OPENMS_FILTERING_TRANSFORMERS_ISOTOPEMARKER_H
#define OPENMS_FILTERING_TRANSFORMERS_ISOTOPEMARKER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopeDistribution.h>

#include <map>
#include <vector>
#include <cmath>
#include <utility>

namespace OpenMS
{

  /**
      @brief IsotopeMarker marks peak pairs which could represent an ion and its isotope

        @todo implement a real isotope marking here with isotopedistributions and fitting (Andreas)

        @htmlinclude OpenMS_IsotopeMarker.parameters

        @ingroup PeakMarker
  */
  class OPENMS_DLLAPI IsotopeMarker :
    public PeakMarker
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    IsotopeMarker();

    /// copy constructor
    IsotopeMarker(const IsotopeMarker & source);

    /// destructor
    ~IsotopeMarker() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    IsotopeMarker & operator=(const IsotopeMarker & source);
    // @}

    // @name Accessors
    // @{
    ///
    static PeakMarker * create() { return new IsotopeMarker(); }

    ///
    template <typename SpectrumType>
    void apply(std::map<double, bool> & marked, SpectrumType & spectrum)
    {
      double mzvariation = (double)param_.getValue("mz_variation");
      double invariation = (double)param_.getValue("in_variation");
      Size marks = param_.getValue("marks");

      spectrum.sortByPosition();

      std::map<double, Size> isotopemarks;        // possible isotopes

      for (Size i = 0; i < spectrum.size(); ++i)
      {
        double mz = spectrum[i].getPosition()[0];
        double intensity = spectrum[i].getIntensity();
        Size j = i + 1;

        //std::vector<std::pair<double, double> > isotopes = SpectrumGenerator::instance()->isotopepeaks(mz, intensity);
        CoarseIsotopeDistribution id;
        id.estimateFromPeptideWeight(mz);

        while (j < spectrum.size() && spectrum[j].getPosition()[0] <= mz + 3 + mzvariation)
        {
          double curmz = spectrum[j].getPosition()[0];
          double curIntensity = spectrum[j].getIntensity();
          UInt iso = (UInt)(curmz - mz + 0.499999);
          if (iso > 0 && curmz - mz - iso > mzvariation)
          {
            ++j;
            continue;
          }
          if (std::fabs(id.begin()->getIntensity() * intensity - curIntensity) < invariation * id.begin()->getIntensity() * intensity)
          {
            isotopemarks[mz]++;
            isotopemarks[curmz]++;
          }
          ++j;
        }
      }

      for (std::map<double, Size>::const_iterator cmit = isotopemarks.begin(); cmit != isotopemarks.end(); ++cmit)
      {
        if (cmit->second >= marks)
        {
          marked.insert(std::pair<double, bool>(cmit->first, true));
        }
      }
      return;
    }

    ///
    static const String getProductName()
    {
      return "IsotopeMarker";
    }

    // @}

  };

}

#endif //OPENMS_FILTERING_TRANSFORMERS_ISOTOPEMARKER_H
