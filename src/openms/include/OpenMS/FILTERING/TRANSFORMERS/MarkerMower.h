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
// $Maintainer: Mathias Walzer $
// $Authors: $
// --------------------------------------------------------------------------
//
#ifndef OPENMS_FILTERING_TRANSFORMERS_MARKERMOWER_H
#define OPENMS_FILTERING_TRANSFORMERS_MARKERMOWER_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/FILTERING/TRANSFORMERS/PeakMarker.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <vector>
#include <map>

namespace OpenMS
{
  /**
    @brief MarkerMower uses PeakMarker to find peaks, those that are not marked get removed

    @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI MarkerMower :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{
    /// default constructor
    MarkerMower();
    /// destructor
    ~MarkerMower() override;

    /// copy constructor
    MarkerMower(const MarkerMower & source);
    /// assignment operator
    MarkerMower & operator=(const MarkerMower & source);
    // @}

    // @name Accessors
    // @{
    ///
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      typedef typename SpectrumType::Iterator Iterator;

      std::map<double, int> marks;
      for (std::vector<PeakMarker *>::const_iterator cvit = markers_.begin(); cvit != markers_.end(); ++cvit)
      {
        std::map<double, bool> marked;
        (*cvit)->apply(marked, spectrum);
        for (std::map<double, bool>::const_iterator cmit = marked.begin(); cmit != marked.end(); ++cmit)
        {
          if (cmit->second)
          {
            marks[cmit->first]++;
          }
        }
      }

      for (Iterator it = spectrum.begin(); it != spectrum.end(); )
      {
        if (marks[it->getMZ()] > 0)
        {
          ++it;
        }
        else
        {
          it = spectrum.erase(it);
        }
      }
    }

    void filterPeakSpectrum(PeakSpectrum & spectrum);

    void filterPeakMap(PeakMap & exp);

    static const String getProductName()
    {
      return "MarkerMower";
    }

    /// insert new Marker (violates the DefaultParamHandler interface)
    void insertmarker(PeakMarker * peak_marker);

    //TODO reimplement DefaultParamHandler::updateMembers_()

    // @}

private:
    /// used peak markers
    std::vector<PeakMarker *> markers_;

  };
}
#endif // OPENMS_COMPARISON_CLUSTERING_MARKERMOWER_H
