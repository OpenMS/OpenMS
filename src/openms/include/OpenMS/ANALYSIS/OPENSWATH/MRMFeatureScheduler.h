// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h> // OPENMS_DLLAPI
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <map>
#include <vector>

namespace OpenMS
{
  /**
    Class used to schedule multiple calls to `MRMFeatureSelector`. It helps with
    settings the parameters for each call to the selector, through the
    `setSchedulerParameters()` method. The class offers a generic scheduler
    (where the user is supposed to pass a `MRMFeatureSelector` derived object) and
    two specialized versions (Score and QMIP).
  */
  class OPENMS_DLLAPI MRMFeatureScheduler
  {
public:
    MRMFeatureScheduler() = default;
    ~MRMFeatureScheduler() = default;

    /**
      Calls `feature_selector.select_MRMFeature()` feeding it the parameters found in `parameters_`.
      It calls said method `parameters_.size()` times, using the result of each cycle as input
      for the next cycle.

      @param[in] feature_selector Base class for the feature selector to use
      @param[in] features Input features
      @param[out] selected_features Selected features
    */
    void scheduleMRMFeatures(
      const MRMFeatureSelector& feature_selector,
      const FeatureMap& features,
      FeatureMap& selected_features
    ) const;

    /// Calls `scheduleMRMFeatures()` using a `MRMFeatureSelectorScore` selector
    void scheduleMRMFeaturesScore(const FeatureMap& features, FeatureMap& selected_features) const;

    /// Calls `scheduleMRMFeatures()` using a `MRMFeatureSelectorQMIP` selector
    void scheduleMRMFeaturesQMIP(const FeatureMap& features, FeatureMap& selected_features) const;

    /// Setter for the scheduler's parameters
    void setSchedulerParameters(const std::vector<MRMFeatureSelector::SelectorParameters>& parameters);

    /// Getter for the scheduler's parameters
    std::vector<MRMFeatureSelector::SelectorParameters>& getSchedulerParameters(void);

private:
    /**
      Parameters for a single call to the scheduler.

      The scheduler goes through each element of this vector. Each of these elements
      contains the parameters' values for a single run of the chosen selector.
    */
    std::vector<MRMFeatureSelector::SelectorParameters> parameters_;
  };
}
