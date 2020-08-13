// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

#include <OpenMS/ANALYSIS/OPENSWATH/MRMBatchFeatureSelector.h>
#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>
#include <OpenMS/KERNEL/FeatureMap.h>

namespace OpenMS
{
  void MRMBatchFeatureSelector::batchMRMFeatures(
    const MRMFeatureSelector& feature_selector,
    const FeatureMap& features,
    FeatureMap& selected_features,
    const std::vector<MRMFeatureSelector::SelectorParameters>& parameters
  )
  {
    FeatureMap input_features = features;
    selected_features.clear();
    for (const MRMFeatureSelector::SelectorParameters& params : parameters)
    {
      feature_selector.selectMRMFeature(input_features, selected_features, params);
      input_features = selected_features;
    }
  }

  void MRMBatchFeatureSelector::batchMRMFeaturesQMIP(
    const FeatureMap& features,
    FeatureMap& selected_features,
    const std::vector<MRMFeatureSelector::SelectorParameters>& parameters
  )
  {
    MRMFeatureSelectorQMIP feature_selector;
    batchMRMFeatures(feature_selector, features, selected_features, parameters);
  }

  void MRMBatchFeatureSelector::batchMRMFeaturesScore(
    const FeatureMap& features,
    FeatureMap& selected_features,
    const std::vector<MRMFeatureSelector::SelectorParameters>& parameters
  )
  {
    MRMFeatureSelectorScore feature_selector;
    batchMRMFeatures(feature_selector, features, selected_features, parameters);
  }
}
