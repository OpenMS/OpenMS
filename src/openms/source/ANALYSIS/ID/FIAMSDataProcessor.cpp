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
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace OpenMS;
using namespace std;


/// default constructor
FIAMSDataProcessor::FIAMSDataProcessor(
  String filename, 
  String dir_input, 
  String dir_output, 
  float resolution, 
  String polarity, 
  String db_mapping, 
  String db_struct, 
  String positive_adducts, 
  String negative_adducts, 
  bool store_progress,
  float min_mz, 
  float max_mz, 
  float bin_step
  )
    : 
    filename_(filename),
    dir_input_(dir_input),
    dir_output_(dir_output),
    resolution_(resolution),
    polarity_(polarity),
    min_mz_(min_mz),
    max_mz_(max_mz),
    bin_step_(bin_step),
    db_mapping_(db_mapping),
    db_struct_(db_struct),
    positive_adducts_(positive_adducts),
    negative_adducts_(negative_adducts),
    store_progress_(store_progress),
    mzs_(),
    bin_sizes_()
  {
  size_t n_bins = static_cast<int> (max_mz_ - min_mz_) / bin_step_;
  for (size_t i = 0; i < n_bins; i++) {
      mzs_.push_back(i*bin_step_);
      bin_sizes_.push_back(mzs_[i] / (resolution_*4.0));
  }
  loadExperiment_();
}

/// default destructor
FIAMSDataProcessor::~FIAMSDataProcessor() {
}

/// copy constructor
FIAMSDataProcessor::FIAMSDataProcessor(const FIAMSDataProcessor& source) :
  filename_(source.filename_),
  dir_input_(source.dir_input_),
  dir_output_(source.dir_output_),
  resolution_(source.resolution_),
  polarity_(source.polarity_),
  min_mz_(source.min_mz_),
  max_mz_(source.max_mz_),
  bin_step_(source.bin_step_),
  db_mapping_(source.db_mapping_),
  db_struct_(source.db_struct_),
  positive_adducts_(source.positive_adducts_),
  negative_adducts_(source.negative_adducts_),
  mzs_(source.mzs_),
  store_progress_(source.store_progress_),
  bin_sizes_(source.bin_sizes_)
  {
  }

/// assignment operator
FIAMSDataProcessor& FIAMSDataProcessor::operator=(const FIAMSDataProcessor& rhs) {
  if (this == &rhs) return *this;
  filename_ = rhs.filename_;
  dir_input_ = rhs.dir_input_;
  dir_output_ = rhs.dir_output_;
  resolution_ = rhs.resolution_;
  polarity_ = rhs.polarity_;
  min_mz_ = rhs.min_mz_;
  max_mz_ = rhs.max_mz_;
  bin_step_ = rhs.bin_step_;
  db_mapping_ = rhs.db_mapping_;
  db_struct_ = rhs.db_struct_;
  positive_adducts_ = rhs.positive_adducts_;
  negative_adducts_ = rhs.negative_adducts_;
  store_progress_ = rhs.store_progress_;
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
    if (store_progress_) {
        MzMLFile mzml;
        MSExperiment exp;
        exp.addSpectrum(output);
        mzml.store(dir_output_ + filename_ + "_merged.mzML", exp);
    }
    return output;
}

FeatureMap FIAMSDataProcessor::extractPeaks(const MSSpectrum & input) {
  MSSpectrum spectrum(input);

  SavitzkyGolayFilter filter;
  filter.filter(spectrum);

  MorphologicalFilter morph_filter;
  morph_filter.filter(spectrum);

  PeakPickerHiRes picker;

  MSSpectrum picked;
  picker.pick(spectrum, picked);

  if (store_progress_) {
      MzMLFile mzml;
      MSExperiment exp;
      exp.addSpectrum(picked);
      mzml.store(dir_output_ + filename_ + "_picked.mzML", exp);
  }

  FeatureMap output;
  for (auto it = picked.begin(); it != picked.end(); ++it) {
      Feature f;
      f.setIntensity(it->getIntensity());
      f.setMZ(it->getMZ());
      output.push_back(f);
      output[0].setMetaValue("scan_polarity", polarity_);
  }
  return output;
}

void FIAMSDataProcessor::runAccurateMassSearch(FeatureMap & input, OpenMS::MzTab & output) {
  Param ams_param;
  ams_param.setValue("ionization_mode", "auto");
  ams_param.setValue("mass_error_value", 1e+06 / (resolution_*2));
  ams_param.setValue("db:mapping", ListUtils::create<String>(db_mapping_));
  ams_param.setValue("db:struct", ListUtils::create<String>(db_struct_));
  ams_param.setValue("positive_adducts", positive_adducts_);
  ams_param.setValue("negative_adducts", negative_adducts_);

  AccurateMassSearchEngine ams;
  ams.setParameters(ams_param);
  ams.init();

  ams.run(input, output);
}

void FIAMSDataProcessor::run(
  float n_seconds,
  OpenMS::MzTab & output
  ) {
  vector<MSSpectrum> output_cut;
  cutForTime(experiment_, output_cut, n_seconds);
  MSSpectrum merged_spectrum = mergeAlongTime(output_cut);
  FeatureMap picked_features = extractPeaks(merged_spectrum);
  runAccurateMassSearch(picked_features, output);
  OpenMS::MzTabFile mztab_outfile;
  mztab_outfile.store(dir_output_ + filename_ + "_" + String(n_seconds) + ".mzTab", output);
}

void FIAMSDataProcessor::loadExperiment_() {
  MzMLFile mzml;
  mzml.load(dir_input_ + filename_ + ".mzML", experiment_);
}

/// Get filename
const String FIAMSDataProcessor::getFilename() {
  return filename_;
}

/// Get input directory
const String FIAMSDataProcessor::getInputDir() {
  return dir_input_;
}

/// Get output directory
const String FIAMSDataProcessor::getOutputDir() {
  return dir_output_;
}

/// Get resolution
const float FIAMSDataProcessor::getResolution() {
  return resolution_;
}

/// Get polarity
const String FIAMSDataProcessor::getPolarity() {
  return polarity_;
}

/// Get minimum mass-to-charge
const float FIAMSDataProcessor::getMinMZ(){
  return min_mz_;
}

/// Get maximum mass-to-charge
const float FIAMSDataProcessor::getMaxMZ(){
  return max_mz_;
}

/// Get the sliding bin step
const float FIAMSDataProcessor::getBinStep(){
  return bin_step_;
}

/// Get the path to the db:mapping for passing to AccurateMassSearch
const String FIAMSDataProcessor::getDBMapping(){
  return db_mapping_;
}

/// Get the path to the db:struct for passing to AccurateMassSearch
const String FIAMSDataProcessor::getDBStruct(){
  return db_struct_;
}

/// Get the path to the positive adducts for passing to AccurateMassSearch
const String FIAMSDataProcessor::getPositiveAdducts(){
  return positive_adducts_;
}

/// Get the path to the negative adducts for passing to AccurateMassSearch
const String FIAMSDataProcessor::getNegativeAdducts(){
  return negative_adducts_;
}

/// Get mass-to-charge ratios to base the sliding window upon
const std::vector<float> FIAMSDataProcessor::getMZs(){
  return mzs_;
}

/// Get the sliding bin sizes
const std::vector<float> FIAMSDataProcessor::getBinSizes(){
  return bin_sizes_;
}
