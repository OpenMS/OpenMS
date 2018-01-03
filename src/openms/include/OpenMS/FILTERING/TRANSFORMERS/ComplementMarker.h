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
#ifndef OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTMARKER_H
#define OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTMARKER_H

#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>

#include <map>
#include <cmath>

namespace OpenMS
{
  /**
    @brief ComplementMarker marks peak pairs which could represent y - b ion pairs

        @htmlinclude OpenMS_ComplementMarker.parameters

        @ingroup PeakMarker
  */
  class OPENMS_DLLAPI ComplementMarker :
    public PeakMarker
  {
public:

    // @name Constructors and Destructors
    //@{
    /// standard constructor
    ComplementMarker();

    /// copy constructor
    ComplementMarker(const ComplementMarker & source);

    /// destructor
    ~ComplementMarker() override;
    //@}

    // @name Operators
    //@{
    /// assignment operator
    ComplementMarker & operator=(const ComplementMarker & source);
    //@}

    // @name Accessors
    //@{
    ///
    static PeakMarker * create() { return new ComplementMarker(); }

    ///
    template <typename SpectrumType>
    void apply(std::map<double, bool> marked, SpectrumType & spectrum)
    {
      if (spectrum.size() < 2)
      {
        return;
      }

      // how often a peak needs to be marked to be returned
      double marks = (double)param_.getValue("marks");
      double parentmass = 0.0;
      if (!spectrum.getPrecursors().empty()) parentmass = spectrum.getPrecursors()[0].getMZ();
      double tolerance = (double)param_.getValue("tolerance");
      std::map<double, int> matching_b_y_ions;

      spectrum.sortByPosition();

      SignedSize j = spectrum.size() - 1;
      for (Size i = 0; i < spectrum.size(); ++i)
      {
        while (j >= 0 && spectrum[j].getPosition()[0] > (parentmass - spectrum[i].getPosition()[0]) + tolerance)
        {
          j--;
        }

        // just takes the first matching ion; todo take all
        if (j >= 0 && std::fabs(spectrum[i].getPosition()[0] + spectrum[j].getPosition()[0] - parentmass) < tolerance)
        {
          matching_b_y_ions[spectrum[i].getPosition()[0]]++;
          matching_b_y_ions[spectrum[j].getPosition()[0]]++;
          j--;
        }
      }

      for (std::map<double, int>::const_iterator cmit = matching_b_y_ions.begin(); cmit != matching_b_y_ions.end(); ++cmit)
      {
        if (cmit->second >= marks)
        {
          marked.insert(std::pair<double, bool>(cmit->first, true));
        }
      }
    }

    /// returns the name to register at the factory
    static const String getProductName()
    {
      return "ComplementMarker";
    }

    //@}

  };

}
#endif //OPENMS_FILTERING_TRANSFORMERS_COMPLEMENTMARKER_H
