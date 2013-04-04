// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#ifndef OPENMS_VISUAL_TOPPVIEWSPECTRAVIEWBEHAVIOR_H
#define OPENMS_VISUAL_TOPPVIEWSPECTRAVIEWBEHAVIOR_H

#include <OpenMS/METADATA/SpectrumSettings.h>
#include <OpenMS/VISUAL/LayerData.h>
#include <vector>
#include <OpenMS/VISUAL/TOPPViewBehaviorInterface.h>

namespace OpenMS
{
  class TOPPViewBase;

  /**
  @brief Behavior of TOPPView in spectra view mode.
  */
  class TOPPViewSpectraViewBehavior :
    public TOPPViewBehaviorInterface
  {
    Q_OBJECT
    ///@name Type definitions
    //@{
    //Feature map type
    typedef LayerData::FeatureMapType FeatureMapType;
    //Feature map managed type
    typedef LayerData::FeatureMapSharedPtrType FeatureMapSharedPtrType;

    //Consensus feature map type
    typedef LayerData::ConsensusMapType ConsensusMapType;
    //Consensus  map managed type
    typedef LayerData::ConsensusMapSharedPtrType ConsensusMapSharedPtrType;

    //Peak map type
    typedef LayerData::ExperimentType ExperimentType;
    //Main managed data type (experiment)
    typedef LayerData::ExperimentSharedPtrType ExperimentSharedPtrType;
    ///Peak spectrum type
    typedef ExperimentType::SpectrumType SpectrumType;
    //@}

public:
    /// Construct the behaviour with its parent
    TOPPViewSpectraViewBehavior(TOPPViewBase * parent);

public slots:
    /// Behavior for showSpectrumAs1D
    virtual void showSpectrumAs1D(int index);

    /// Behavior for showSpectrumAs1D
    virtual void showSpectrumAs1D(std::vector<int, std::allocator<int> > indices);

    /// Behavior for activate1DSpectrum
    virtual void activate1DSpectrum(int index);

    /// Behavior for activate1DSpectrum
    virtual void activate1DSpectrum(std::vector<int, std::allocator<int> > indices);

    /// Behavior for deactivate1DSpectrum
    virtual void deactivate1DSpectrum(int index);

    /// Slot for behavior activation
    virtual void activateBehavior();

    /// Slot for behavior deactivation
    virtual void deactivateBehavior();
private:
    TOPPViewBase * tv_;
  };
}
#endif // OPENMS_VISUAL_TOPPVIEWSPECTRAVIEWBEHAVIOR_H
