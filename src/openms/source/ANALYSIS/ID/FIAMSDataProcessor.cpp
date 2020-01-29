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
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>
#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>
#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>
using namespace OpenMS;
using namespace std;


/// default constructor
FIAMSDataProcessor::FIAMSDataProcessor(float resolution, float min_mz, float max_mz, float bin_step)
    : 
    resolution_(resolution),
    min_mz_(min_mz),
    max_mz_(max_mz),
    bin_step_(bin_step),
    mzs_(),
    bin_sizes_()
  {
  size_t n_bins = static_cast<int> (max_mz_ - min_mz_) / bin_step_;
  for (size_t i = 0; i < n_bins; i++) {
      mzs_.push_back(i*bin_step_);
      bin_sizes_.push_back(mzs_[i] / (resolution_*4.0));
  }
}

/// default destructor
FIAMSDataProcessor::~FIAMSDataProcessor() {
}

/// copy constructor
FIAMSDataProcessor::FIAMSDataProcessor(const FIAMSDataProcessor& source) :
  resolution_(source.resolution_),
  min_mz_(source.min_mz_),
  max_mz_(source.max_mz_),
  bin_step_(source.bin_step_),
  mzs_(source.mzs_),
  bin_sizes_(source.bin_sizes_)
  {
  }

/// assignment operator
FIAMSDataProcessor& FIAMSDataProcessor::operator=(const FIAMSDataProcessor& rhs) {
  if (this == &rhs) return *this;
  resolution_ = rhs.resolution_;
  min_mz_ = rhs.min_mz_;
  max_mz_ = rhs.max_mz_;
  bin_step_ = rhs.bin_step_;
  mzs_ = rhs.mzs_;
  bin_sizes_ = rhs.bin_sizes_;
  return *this;
}

void FIAMSDataProcessor::cutForTime(const MSExperiment & experiment, vector<MSSpectrum> & output, float n_seconds) {
    for (auto s : experiment.getSpectra()) {
        if (s.getRT() < n_seconds) output.push_back(s);
    }
}

MSSpectrum FIAMSDataProcessor::mergeAlongTime(
  const std::vector<MSSpectrum> & input
  ) {
    MSSpectrum output;
    for (size_t i = 1; i < mzs_.size() - 1; i++) {
        OpenMS::MSSpectrum full_spectrum = OpenMS::SpectrumAddition::addUpSpectra(
            input, bin_sizes_[i], false
        );
        for (auto it = full_spectrum.begin(); it != full_spectrum.end(); ++it) {
            if (it->getMZ() > mzs_[i+1]) break;
            if (it->getMZ() >= mzs_[i]) output.push_back(*it);
        }
    }
    output.sortByPosition();
    return output;
}
