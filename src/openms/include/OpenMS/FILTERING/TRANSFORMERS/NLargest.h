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
#ifndef OPENMS_FILTERING_TRANSFORMERS_NLARGEST_H
#define OPENMS_FILTERING_TRANSFORMERS_NLARGEST_H

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

namespace OpenMS
{

  /**
    @brief NLargest removes all but the n largest peaks

    @htmlinclude OpenMS_NLargest.parameters

    @ingroup SpectraPreprocessers
  */
  class OPENMS_DLLAPI NLargest :
    public DefaultParamHandler
  {
public:

    // @name Constructors and Destructors
    // @{

    /// default constructor
    NLargest();
    /// detailed constructor
    NLargest(UInt n);
    /// destructor
    ~NLargest() override;

    /// copy constructor
    NLargest(const NLargest & source);
    /// assignment operator
    NLargest & operator=(const NLargest & source);

    // @}

    ///
    template <typename SpectrumType>
    void filterSpectrum(SpectrumType & spectrum)
    {
      if (spectrum.size() <= peakcount_) return;

      // sort by reverse intensity
      spectrum.sortByIntensity(true);

      // keep the n largest peaks if more than n are present
      std::vector<Size> indices;
      for (Size i = 0; i != peakcount_; ++i)
      {
        indices.push_back(i);
      }
      spectrum.select(indices);
    }

    void filterPeakSpectrum(PeakSpectrum & spectrum);

    void filterPeakMap(PeakMap & exp);

    //TODO reimplement DefaultParamHandler::updateMembers_()

    // @}

protected:
    void updateMembers_() override;
    UInt peakcount_;

    /// handles the initialization of the default parameters for the 2 constructors
    void init_();

  };

}
#endif //OPENMS_FILTERING_TRANSFORMERS_NLARGEST_H
