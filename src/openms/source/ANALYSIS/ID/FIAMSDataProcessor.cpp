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
// $Authors: Svetlana Kutuzova, Douglas McCloskey $
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
#include <OpenMS/FORMAT/CsvFile.h>

#include <boost/algorithm/string.hpp>
#ifdef _OPENMP
#include <omp.h>
#endif

namespace OpenMS {
  /// default constructor
  FIAMSDataProcessor::FIAMSDataProcessor() :
      DefaultParamHandler("FIAMSDataProcessor"),
      mzs_(),
      bin_sizes_()
    {
    defaults_.setValue("filename", "fiams", "The filename to use for naming the output files");
    defaults_.setValue("dir_output", "", "The path to the directory where the output files will be placed");

    defaults_.setValue("resolution", 120000.0, "The instrument settings: resolution");

    defaults_.setValue("polarity", "positive", "The instrument settings: polarity");
    defaults_.setValidStrings("polarity", ListUtils::create<String>("positive,negative"));

    defaults_.setValue("max_mz", 1500, "Maximum mz");

    defaults_.setValue("bin_step", 20, "The size of the step to recalculated the bin size used for adding up spectra along the time axis");

    defaults_.setValue("db:mapping", ListUtils::create<String>("CHEMISTRY/HMDBMappingFile.tsv"), "For the accurate mass search. Database input file(s), containing three tab-separated columns of mass, formula, identifier. "
                                                                      "If 'mass' is 0, it is re-computed from the molecular sum formula. "
                                                                      "By default CHEMISTRY/HMDBMappingFile.tsv in OpenMS/share is used! If empty, the default will be used.");
    defaults_.setValue("db:struct", ListUtils::create<String>("CHEMISTRY/HMDB2StructMapping.tsv"), "For the accurate mass search. Database input file(s), containing four tab-separated columns of identifier, name, SMILES, INCHI."
                                                                        "The identifier should match with mapping file. SMILES and INCHI are reported in the output, but not used otherwise. "
                                                                        "By default CHEMISTRY/HMDB2StructMapping.tsv in OpenMS/share is used! If empty, the default will be used.");
    defaults_.setValue("positive_adducts", "CHEMISTRY/PositiveAdducts.tsv", "For the accurate mass search. This file contains the list of potential positive adducts that will be looked for in the database. "
                                                                                  "Edit the list if you wish to exclude/include adducts. "
                                                                                  "By default CHEMISTRY/PositiveAdducts.tsv in OpenMS/share is used! If empty, the default will be used.", ListUtils::create<String>("advanced"));
    defaults_.setValue("negative_adducts", "CHEMISTRY/NegativeAdducts.tsv", "For the accurate mass search. This file contains the list of potential negative adducts that will be looked for in the database. "
                                                                                  "Edit the list if you wish to exclude/include adducts. "
                                                                                  "By default CHEMISTRY/NegativeAdducts.tsv in OpenMS/share is used! If empty, the default will be used.", ListUtils::create<String>("advanced"));

    defaults_.setValue("store_progress", "true", "If the intermediate files should be stored in the output directory");
    defaultsToParam_();
  }

  /// Copy constructor
  FIAMSDataProcessor::FIAMSDataProcessor(const FIAMSDataProcessor& source) :
    DefaultParamHandler(source),
    mzs_(source.mzs_),
    bin_sizes_(source.bin_sizes_)
  {
  }

  /// assignment operator
  FIAMSDataProcessor& FIAMSDataProcessor::operator=(const FIAMSDataProcessor& rhs) {
    if (this == &rhs) return *this;
    mzs_ = rhs.mzs_;
    bin_sizes_ = rhs.bin_sizes_;
    return *this;
  }

  void FIAMSDataProcessor::updateMembers_() {
    float max_mz_ = param_.getValue("max_mz");
    float bin_step_ = param_.getValue("bin_step");
    float resolution_ = static_cast<float>(param_.getValue("resolution"));
    size_t n_bins = static_cast<int> (max_mz_ / bin_step_);
    mzs_.clear();
    bin_sizes_.clear();
    for (size_t i = 0; i < n_bins; i++) {
        mzs_.push_back((i+1)*bin_step_);
        bin_sizes_.push_back(mzs_[i] / (resolution_*4.0));
    }
  }

  void FIAMSDataProcessor::cutForTime(const MSExperiment & experiment, std::vector<MSSpectrum> & output, const float & n_seconds) {
      for (const auto & s : experiment.getSpectra()) {
          if (s.getRT() < n_seconds) output.push_back(s);
      }
  }

  MSSpectrum FIAMSDataProcessor::mergeAlongTime(
    const std::vector<MSSpectrum> & input
    ) {
      MSSpectrum output;
      for (size_t i = 0; i < mzs_.size() - 1; i++) {
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

  MSSpectrum FIAMSDataProcessor::extractPeaks(const MSSpectrum & input) {
    MSSpectrum spectrum(input);

    SavitzkyGolayFilter filter;
    filter.filter(spectrum);

    MorphologicalFilter morph_filter;
    morph_filter.filter(spectrum);

    PeakPickerHiRes picker;

    MSSpectrum picked;
    picker.pick(spectrum, picked);

    return picked;
  }

  FeatureMap FIAMSDataProcessor::convertToFeatureMap(const MSSpectrum & input) {
    String polarity_ = param_.getValue("polarity");
    FeatureMap output;
    for (auto it = input.begin(); it != input.end(); ++it) {
        Feature f;
        f.setIntensity(it->getIntensity());
        f.setMZ(it->getMZ());
        f.setMetaValue("scan_polarity", polarity_);
        output.push_back(f);
    }
    return output;
  }

  void FIAMSDataProcessor::runAccurateMassSearch(FeatureMap & input, OpenMS::MzTab & output) {
    Param ams_param;
    ams_param.setValue("ionization_mode", "auto");
    ams_param.setValue("mass_error_value", 1e+06 / (static_cast<float>(param_.getValue("resolution"))*2));
    ams_param.setValue("db:mapping", param_.getValue("db:mapping"));
    ams_param.setValue("db:struct", param_.getValue("db:struct"));
    ams_param.setValue("positive_adducts", param_.getValue("positive_adducts"));
    ams_param.setValue("negative_adducts", param_.getValue("negative_adducts"));

    AccurateMassSearchEngine ams;
    ams.setParameters(ams_param);
    ams.init();

    ams.run(input, output);
  }

  MSSpectrum FIAMSDataProcessor::trackNoise(const MSSpectrum & input) {
    SignalToNoiseEstimatorMedianRapid sne(50);
    MSSpectrum output;
    if (input.size() == 0) {
      return output;
    }
    std::vector<double> mzs, intensities;
    for (auto it = input.begin(); it != input.end(); ++it)
    {
        mzs.push_back(it->getMZ());
        intensities.push_back(it->getIntensity());
    }
    SignalToNoiseEstimatorMedianRapid::NoiseEstimator e = sne.estimateNoise(mzs, intensities);
    
    for (auto it = input.begin(); it != input.end(); ++it)
    {
        Peak1D peak;
        peak.setMZ(it->getMZ());
        peak.setIntensity(e.get_noise_value(it->getMZ()));
        output.push_back(peak);
    }
    return output;
  }

  void FIAMSDataProcessor::run(const MSExperiment & experiment, const float & n_seconds, OpenMS::MzTab & output) {
    std::vector<MSSpectrum> output_cut;
    cutForTime(experiment, output_cut, n_seconds);
    MSSpectrum merged_spectrum = mergeAlongTime(output_cut);
    MSSpectrum picked_spectrum = extractPeaks(merged_spectrum);
    String postfix = String(static_cast<int>(n_seconds));
    String dir_output_ = param_.getValue("dir_output");
    String filename_ = param_.getValue("filename");
    if (param_.getValue("store_progress").toBool()) {
      storeSpectrum_(merged_spectrum, dir_output_ + "/" + filename_ + "_merged_" + postfix + ".mzML");
      storeSpectrum_(picked_spectrum, dir_output_ + "/" + filename_ + "_picked_" + postfix + ".mzML");
    }
    FeatureMap picked_features = convertToFeatureMap(picked_spectrum);
    MSSpectrum signal_to_noise = trackNoise(picked_spectrum);
    storeSpectrum_(signal_to_noise, dir_output_ + "/" + filename_ + "_signal_to_noise_" + postfix + ".mzML");
    runAccurateMassSearch(picked_features, output);
    OpenMS::MzTabFile mztab_outfile;
    mztab_outfile.store(dir_output_ + "/" + filename_ + "_" + postfix + ".mzTab", output);
  }

  void FIAMSDataProcessor::storeSpectrum_(const MSSpectrum & input, String filename) {
      MzMLFile mzml;
      MSExperiment exp;
      exp.addSpectrum(input);
      mzml.store(filename, exp);
  }

  /// Get mass-to-charge ratios to base the sliding window upon
  const std::vector<float> FIAMSDataProcessor::getMZs(){
    return mzs_;
  }

  /// Get the sliding bin sizes
  const std::vector<float> FIAMSDataProcessor::getBinSizes(){
    return bin_sizes_;
  }

  /// default constructor
  FIAMSScheduler::FIAMSScheduler(
    String filename,
    String base_dir
  )
    : 
    filename_(filename),
    base_dir_(base_dir),
    samples_()
  {
  loadSamples_();
  }

  /// assignment operator
  FIAMSScheduler& FIAMSScheduler::operator=(const FIAMSScheduler& rhs) {
    if (this == &rhs) return *this;
    filename_ = rhs.filename_;
    base_dir_ = rhs.base_dir_;
    samples_ = rhs.samples_;
    return *this;
  }

  void FIAMSScheduler::loadSamples_() {
    CsvFile csv_file(filename_, ',');
    StringList headers;
    csv_file.getRow(0, headers);
    StringList row;
    for (Size i = 1; i < csv_file.rowCount(); ++i) {
      csv_file.getRow(i, row);
      std::map<String, String> mapping;
      for (Size j = 0; j < headers.size(); ++j) {
        mapping[headers[j]] = row[j];
      }
      samples_.push_back(mapping);
    }
  }

  void FIAMSScheduler::loadExperiment_(const String & dir_input, const String & filename, MSExperiment output) {
    MzMLFile mzml;
    mzml.load(dir_input + "/" + filename + ".mzML", output);
  }

  void FIAMSScheduler::run() {
    #pragma omp parallel for
    for (size_t i = 0; i < samples_.size(); ++i) {
      MSExperiment exp;
      loadExperiment_(base_dir_ + samples_[i].at("dir_input"), samples_[i].at("filename"), exp);

      FIAMSDataProcessor fia_processor;
      Param p;
      p.setValue("filename", samples_[i].at("filename"));
      p.setValue("dir_output", base_dir_ + samples_[i].at("dir_output"));
      p.setValue("resolution", std::stof(samples_[i].at("resolution")));
      p.setValue("polarity", samples_[i].at("charge"));
      p.setValue("db:mapping", ListUtils::create<String>(base_dir_ + samples_[i].at("db_mapping")));
      p.setValue("db:struct", ListUtils::create<String>(base_dir_ + samples_[i].at("db_struct")));
      p.setValue("positive_adducts", base_dir_ + samples_[i].at("positive_adducts"));
      p.setValue("negative_adducts", base_dir_ + samples_[i].at("negative_adducts"));
      fia_processor.setParameters(p);

      String time = samples_[i].at("time");
      std::vector<String> times;
      boost::split(times, time, boost::is_any_of(";"));
      for (size_t j = 0; j < times.size(); ++j) {
        std::cout << "Started " << samples_[i].at("filename") << " for " << times[j] << " seconds" << std::endl;
        MzTab mztab_output;
        fia_processor.run(exp, std::stof(times[j]), mztab_output);
        std::cout << "Finished " << samples_[i].at("filename") << " for " << times[j] << " seconds" << std::endl;
      }
    }
  }

  const std::vector<std::map<String, String>> FIAMSScheduler::getSamples() {
    return samples_;
  }

  const String FIAMSScheduler::getBaseDir() {
    return base_dir_;
  }
}