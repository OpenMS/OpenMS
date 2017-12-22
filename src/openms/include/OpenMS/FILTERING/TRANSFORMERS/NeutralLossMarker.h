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
#ifndef OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSMARKER_H
#define OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSMARKER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>

#include <map>
#include <cmath>

namespace OpenMS
{
  /**
    @brief NeutralLossMarker marks peak pairs which could represent an ion an its neutral loss (water, ammonia)

        @htmlinclude OpenMS_NeutralLossMarker.parameters

        @ingroup PeakMarker
  */
  class OPENMS_DLLAPI NeutralLossMarker :
    public PeakMarker
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    NeutralLossMarker();

    /// copy constructor
    NeutralLossMarker(const NeutralLossMarker & source);

    /// destructor
    ~NeutralLossMarker() override;
    // @}

    // @name Operators
    // @{
    /// assignment operator
    NeutralLossMarker & operator=(const NeutralLossMarker & source);
    // @}

    // @name Accessors
    // @{
    ///
    static PeakMarker * create() { return new NeutralLossMarker(); }

    ///
    template <typename SpectrumType>
    void apply(std::map<double, bool> & marked, SpectrumType & spectrum)
    {
      // how often a peak needs to be marked to be returned
      double marks = (double)param_.getValue("marks");
      double tolerance = (double)param_.getValue("tolerance");
      std::map<double, SignedSize> ions_w_neutrallosses;
      spectrum.sortByPosition();
      for (Size i = 0; i < spectrum.size(); ++i)
      {
        double mz = spectrum[i].getPosition()[0];
        double intensity = spectrum[i].getIntensity();
        SignedSize j = i - 1;
        while (j >= 0)
        {
          double curmz = spectrum[j].getPosition()[0];
          double curIntensity = spectrum[j].getIntensity();

          // check for peak that's a water or ammonia away
          if (std::fabs(mz - curmz - 17) < tolerance || std::fabs(mz - curmz - 18) < tolerance)
          {
            // neutral loss peak should be smaller
            if (curIntensity < intensity)
            {
              ions_w_neutrallosses[mz]++;
              // neutral loss peak not marked
              //ions_w_neutrallosses[curmz]++;
            }
          }
          else
          {
            if (mz - curmz > 18.3)
            {
              break;
            }
          }
          --j;
        }
      }

      for (std::map<double, SignedSize>::const_iterator cmit = ions_w_neutrallosses.begin(); cmit != ions_w_neutrallosses.end(); ++cmit)
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
      return "NeutralLossMarker";
    }

    // @}

  };

}
#endif //OPENMS_FILTERING_TRANSFORMERS_NEUTRALLOSSMARKER_H
