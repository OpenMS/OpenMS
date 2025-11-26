// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/DATASTRUCTURES/DBoundingBox.h>
#include <OpenMS/IONMOBILITY/FAIMSHelper.h>
#include <OpenMS/IONMOBILITY/IMDataConverter.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/MATH/MathFunctions.h>
#include <OpenMS/SYSTEM/StopWatch.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
#include <set>
#include <limits>

using namespace std;

namespace OpenMS
{

Biosaur2Algorithm::Biosaur2Algorithm() :
  DefaultParamHandler("Biosaur2Algorithm")
{
  defaults_.setValue("mini", 1.0, "Minimum intensity threshold");
  defaults_.setMinFloat("mini", 0.0);

  defaults_.setValue("minmz", 350.0, "Minimum m/z value");
  defaults_.setMinFloat("minmz", 0.0);

  defaults_.setValue("maxmz", 1500.0, "Maximum m/z value");
  defaults_.setMinFloat("maxmz", 0.0);

  defaults_.setValue("htol", 8.0, "Mass accuracy in ppm for combining peaks into hills");
  defaults_.setMinFloat("htol", 0.0);

  defaults_.setValue("itol", 8.0, "Mass accuracy in ppm for isotopic patterns");
  defaults_.setMinFloat("itol", 0.0);

  defaults_.setValue("hvf", 1.3, "Hill valley factor for splitting hills");
  defaults_.setMinFloat("hvf", 1.0);

  defaults_.setValue("ivf", 5.0, "Isotope valley factor for splitting isotope patterns");
  defaults_.setMinFloat("ivf", 1.0);

  defaults_.setValue("minlh", 2, "Minimum number of scans for a hill");
  defaults_.setMinInt("minlh", 1);

  defaults_.setValue("pasefmini", 100.0, "Minimum combined intensity for PASEF/TIMS clusters after m/z–ion-mobility centroiding.");
  defaults_.setMinFloat("pasefmini", 0.0);

  defaults_.setValue("pasefminlh", 1, "Minimum number of raw points per PASEF/TIMS cluster during centroiding.");
  defaults_.setMinInt("pasefminlh", 1);

  defaults_.setValue("cmin", 1, "Minimum charge state");
  defaults_.setMinInt("cmin", 1);

  defaults_.setValue("cmax", 6, "Maximum charge state");
  defaults_.setMinInt("cmax", 1);

  defaults_.setValue("iuse", 0, "Number of isotopes for intensity calculation (0=mono only, -1=all, 1=mono+first, etc.)");
  defaults_.setMinInt("iuse", -1);

  defaults_.setValue("nm", "false", "Negative mode (affects neutral mass calculation)");
  defaults_.setValidStrings("nm", {"true", "false"});

  defaults_.setValue("tof", "false", "Enable TOF-specific intensity filtering");
  defaults_.setValidStrings("tof", {"true", "false"});

  defaults_.setValue("profile", "false", "Enable profile mode processing (centroid spectra using PeakPickerHiRes)");
  defaults_.setValidStrings("profile", {"true", "false"});

  defaults_.setValue("paseftol", 0.0, "Ion mobility accuracy (in the same units as the ion-mobility array) for linking peaks into hills and grouping isotopes (0 = disable IM-based gating).");
  defaults_.setMinFloat("paseftol", 0.0);

  defaults_.setValue("use_hill_calib", "false", "Enable automatic hill mass tolerance calibration");
  defaults_.setValidStrings("use_hill_calib", {"true", "false"});

  defaults_.setValue("ignore_iso_calib", "false", "Disable automatic isotope mass error calibration");
  defaults_.setValidStrings("ignore_iso_calib", {"true", "false"});

  defaults_.setValue("hrttol", 10.0, "Maximum allowed RT difference (in seconds) between monoisotopic hill apex and isotope hill apex when assembling isotope patterns (0 disables RT gating).");
  defaults_.setMinFloat("hrttol", 0.0);

  defaults_.setValue("convex_hulls", "bounding_box",
                     "Representation of feature convex hulls in the output FeatureMap. "
                     "'bounding_box' stores a single RT–m/z bounding box per feature "
                     "(smaller featureXML, no per-trace detail), "
                     "whereas 'mass_traces' stores one convex hull per contributing hill using all mass-trace points "
                     "(larger featureXML, preserves detailed trace shape).");
  defaults_.setValidStrings("convex_hulls", {"mass_traces", "bounding_box"});

  defaultsToParam_();
  updateMembers_();
}

void Biosaur2Algorithm::updateMembers_()
{
  mini_ = param_.getValue("mini");
  minmz_ = param_.getValue("minmz");
  maxmz_ = param_.getValue("maxmz");
  htol_ = param_.getValue("htol");
  itol_ = param_.getValue("itol");
  hvf_ = param_.getValue("hvf");
  ivf_ = param_.getValue("ivf");
  minlh_ = static_cast<Size>(param_.getValue("minlh"));
  cmin_ = param_.getValue("cmin");
  cmax_ = param_.getValue("cmax");
  pasefmini_ = param_.getValue("pasefmini");
  pasefminlh_ = static_cast<Size>(param_.getValue("pasefminlh"));
  iuse_ = param_.getValue("iuse");
  negative_mode_ = param_.getValue("nm").toBool();
  tof_mode_ = param_.getValue("tof").toBool();
  profile_mode_ = param_.getValue("profile").toBool();
  use_hill_calib_ = param_.getValue("use_hill_calib").toBool();
  ignore_iso_calib_ = param_.getValue("ignore_iso_calib").toBool();
  paseftol_ = param_.getValue("paseftol");
  hrttol_ = param_.getValue("hrttol");
  convex_hull_mode_ = param_.getValue("convex_hulls").toString();

  OPENMS_LOG_DEBUG << "Biosaur2Algorithm parameters after updateMembers_: "
                   << "mini=" << mini_
                   << ", minmz=" << minmz_
                   << ", maxmz=" << maxmz_
                   << ", htol=" << htol_
                   << ", itol=" << itol_
                   << ", hvf=" << hvf_
                   << ", ivf=" << ivf_
                   << ", minlh=" << minlh_
                   << ", cmin=" << cmin_
                   << ", cmax=" << cmax_
                   << ", pasefmini=" << pasefmini_
                   << ", pasefminlh=" << pasefminlh_
                   << ", iuse=" << iuse_
                   << ", negative_mode=" << negative_mode_
                   << ", tof_mode=" << tof_mode_
                   << ", profile_mode=" << profile_mode_
                   << ", use_hill_calib=" << use_hill_calib_
                   << ", ignore_iso_calib=" << ignore_iso_calib_
                   << ", paseftol=" << paseftol_
                   << ", hrttol=" << hrttol_
                   << ", convex_hulls=" << convex_hull_mode_
                   << endl;
}

void Biosaur2Algorithm::setMSData(const MSExperiment& ms_data)
{
  ms_data_ = ms_data;
}

void Biosaur2Algorithm::setMSData(MSExperiment&& ms_data)
{
  ms_data_ = std::move(ms_data);
}

MSExperiment& Biosaur2Algorithm::getMSData()
{
  return ms_data_;
}

const MSExperiment& Biosaur2Algorithm::getMSData() const
{
  return ms_data_;
}

void Biosaur2Algorithm::run(FeatureMap& feature_map)
{
  vector<Hill> tmp_hills;
  vector<PeptideFeature> tmp_features;
  run(feature_map, tmp_hills, tmp_features);
}

void Biosaur2Algorithm::run(FeatureMap& feature_map,
                            vector<Hill>& hills,
                            vector<PeptideFeature>& peptide_features)
{
  feature_map.clear(true);
  hills.clear();
  peptide_features.clear();

  // Filter to keep only MS1 spectra
  ms_data_.getSpectra().erase(
    remove_if(ms_data_.begin(), ms_data_.end(),
              [](const MSSpectrum& s) { return s.getMSLevel() != 1; }),
    ms_data_.end());
  
  if (profile_mode_)
  {
    centroidProfileSpectra_(ms_data_);
  }
  
  if (tof_mode_)
  {
    processTOF_(ms_data_);
  }

  OPENMS_LOG_INFO << "Loaded " << ms_data_.size() << " MS1 spectra" << endl;

  if (ms_data_.empty())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No MS1 spectra found in input experiment!", "");
  }

  // Build FAIMS-aware processing groups (one per CV or a single non-FAIMS group).
  std::vector<std::pair<double, MSExperiment>> groups =
    IMDataConverter::splitByFAIMSCV(std::move(ms_data_));

  const Size n_groups = groups.size();
  vector<vector<Hill>> hills_per_group(n_groups);
  vector<vector<PeptideFeature>> features_per_group(n_groups);
  vector<FeatureMap> fmap_per_group(n_groups);

  const double original_paseftol = paseftol_;

  // Parallelize processing across FAIMS groups. Each group is handled
  // independently using its own local hills/features containers.
  #ifdef _OPENMP
  #pragma omp parallel for schedule(dynamic)
  #endif
  for (int i = 0; i < static_cast<int>(n_groups); ++i)
  {
    auto& group_pair = groups[static_cast<Size>(i)];
    double cv = group_pair.first;
    MSExperiment& group_exp = group_pair.second;

    vector<Hill> local_hills;
    vector<PeptideFeature> local_features;

    processFAIMSGroup_(cv, group_exp, original_paseftol, local_hills, local_features);

    FeatureMap local_map;
    if (!local_hills.empty() && !local_features.empty())
    {
      local_map = convertToFeatureMap_(local_features, local_hills);
    }

    hills_per_group[static_cast<Size>(i)] = std::move(local_hills);
    features_per_group[static_cast<Size>(i)] = std::move(local_features);
    fmap_per_group[static_cast<Size>(i)] = std::move(local_map);
  }

  // Restore original paseftol_ value for subsequent calls.
  paseftol_ = original_paseftol;

  vector<Hill> all_hills;
  vector<PeptideFeature> all_features;
  FeatureMap combined_feature_map;

  for (Size i = 0; i < n_groups; ++i)
  {
    all_hills.insert(all_hills.end(),
                     hills_per_group[i].begin(), hills_per_group[i].end());
    all_features.insert(all_features.end(),
                        features_per_group[i].begin(), features_per_group[i].end());

    FeatureMap& gm = fmap_per_group[i];
    if (!gm.empty())
    {
      combined_feature_map.insert(combined_feature_map.end(), gm.begin(), gm.end());
    }
  }

  hills = std::move(all_hills);
  peptide_features = std::move(all_features);

  feature_map = std::move(combined_feature_map);
  feature_map.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
  feature_map.ensureUniqueId();
  feature_map.getProteinIdentifications().resize(1);
}

double Biosaur2Algorithm::calculatePPM_(double mz1, double mz2) const
{
  if (fabs(mz2) < 1e-10) return 0.0;
  return Math::getPPM(mz1, mz2);
}

double Biosaur2Algorithm::calculateMedian_(const vector<double>& values) const
{
  if (values.empty()) return 0.0;
  vector<double> sorted = values;
  return Math::median(sorted.begin(), sorted.end(), false);
}

bool Biosaur2Algorithm::shouldThrowForMissingIM_(const MSSpectrum& spectrum) const
{
  // The Python reference implementation degrades gracefully when ion-mobility
  // arrays are missing by simply disabling IM-based gating. To mirror this
  // behavior, we never treat missing per-peak IM arrays as a hard error here.
  (void)spectrum;
  return false;
}

double Biosaur2Algorithm::cosineCorrelation1D_(const vector<double>& v1,
                                               const vector<double>& v2) const
{
  Size n = min(v1.size(), v2.size());
  if (n == 0) return 0.0;

  double dot = 0.0;
  double n1 = 0.0;
  double n2 = 0.0;
  for (Size i = 0; i < n; ++i)
  {
    dot += v1[i] * v2[i];
    n1 += v1[i] * v1[i];
    n2 += v2[i] * v2[i];
  }
  if (n1 == 0.0 || n2 == 0.0) return 0.0;
  return dot / (sqrt(n1) * sqrt(n2));
}

pair<double, Size> Biosaur2Algorithm::checkingCosCorrelationForCarbon_(
  const vector<double>& theor_full,
  const vector<double>& exp_full,
  double thresh) const
{
  if (exp_full.size() <= 1 || theor_full.empty())
  {
    return make_pair(0.0, Size(1));
  }

  double theor_total_sum = accumulate(theor_full.begin(), theor_full.end(), 0.0);
  if (theor_total_sum <= 0.0)
  {
    return make_pair(0.0, Size(1));
  }

  double best_cor = 0.0;
  Size best_pos = 1;

  int pos = static_cast<int>(exp_full.size());
  while (pos != 1)
  {
    Size suit_len = min(static_cast<Size>(pos), theor_full.size());
    vector<double> theor(theor_full.begin(), theor_full.begin() + suit_len);
    vector<double> exp(exp_full.begin(), exp_full.begin() + pos);

    double theor_partial_sum = accumulate(theor.begin(), theor.end(), 0.0);
    double averagine_explained =
      (theor_total_sum > 0.0) ? (theor_partial_sum / theor_total_sum) : 0.0;

    double cor = cosineCorrelation1D_(theor, exp);

    if (averagine_explained >= 0.5 && cor >= thresh)
    {
      if (cor > best_cor)
      {
        best_cor = cor;
        best_pos = static_cast<Size>(pos);
      }
      break;
    }

    --pos;
  }

  return make_pair(best_cor, best_pos);
}

vector<double> Biosaur2Algorithm::meanFilter_(const vector<double>& data, Size window) const
{
  vector<double> result(data.size());
  if (data.empty()) return result;

  // Treat 'window' as the half-width of the kernel, i.e.
  // kernel length = 2 * window + 1, mirroring the Python meanfilt
  // implementation with NumPy's 'same' convolution (zero padding).
  Size half_window = window;
  Size kernel_len = 2 * half_window + 1;
  if (kernel_len == 0) return result;

  for (Size i = 0; i < data.size(); ++i)
  {
    double sum = 0.0;
    // Convolution with implicit zeros outside [0, data.size()).
    for (Size k = 0; k < kernel_len; ++k)
    {
      // Corresponding index in the input signal
      // centered at position i.
      long j = static_cast<long>(i) + static_cast<long>(k) - static_cast<long>(half_window);
      if (j >= 0 && static_cast<Size>(j) < data.size())
      {
        sum += data[static_cast<Size>(j)];
      }
    }
    result[i] = sum / static_cast<double>(kernel_len);
  }

  return result;
}

pair<double, double> Biosaur2Algorithm::calibrateMass_(const vector<double>& mass_errors, double bin_width) const
{
  if (mass_errors.empty())
  {
    return make_pair(0.0, 10.0);
  }

  double min_error = *min_element(mass_errors.begin(), mass_errors.end());
  double max_error = *max_element(mass_errors.begin(), mass_errors.end());
  double mass_left = -min_error;
  double mass_right = max_error;

  int n_bins = static_cast<int>((mass_left + mass_right) / bin_width);
  if (n_bins < 5)
  {
    return make_pair(0.0, 10.0);
  }

  vector<double> bin_centers;
  vector<int> bin_counts(n_bins, 0);

  for (int i = 0; i < n_bins; ++i)
  {
    bin_centers.push_back(-mass_left + (i + 0.5) * bin_width);
  }

  for (double error : mass_errors)
  {
    int bin = static_cast<int>((error + mass_left) / bin_width);
    if (bin >= 0 && bin < n_bins)
    {
      bin_counts[bin]++;
    }
  }

  double sum_x = 0.0, sum_x2 = 0.0, sum_w = 0.0;
  for (size_t i = 0; i < bin_centers.size(); ++i)
  {
    double x = bin_centers[i];
    double w = bin_counts[i];
    sum_x += w * x;
    sum_x2 += w * x * x;
    sum_w += w;
  }

  if (sum_w < 10)
  {
    return make_pair(0.0, 10.0);
  }

  double mean = sum_x / sum_w;
  double variance = (sum_x2 / sum_w) - (mean * mean);
  double sigma = sqrt(max(variance, 0.01));

  if (fabs(mean) >= max(mass_left, mass_right))
  {
    return calibrateMass_(mass_errors, 0.25);
  }

  if (isinf(sigma) || isnan(sigma))
  {
    return make_pair(0.0, 10.0);
  }

  return make_pair(mean, sigma);
}

pair<vector<double>, Size> Biosaur2Algorithm::computeAveragine_(double neutral_mass,
                                                                double apex_intensity) const
{
  // Averagine-based theoretical isotope intensities (C-only binomial model),
  // using the same neutral-mass binning (100 Da bins) and binomial
  // parameters as the reference Cython/SciPy implementation.
  constexpr int averagine_mass_bin_step = 100;
  constexpr int averagine_max_mass_bin_index = 199; // 0..199 => 0..19900 Da

  static vector<vector<double>> averagine_table;
  static vector<Size> averagine_max_pos;
  static bool averagine_initialized = false;

  if (!averagine_initialized)
  {
    const double averagine_mass = 111.1254;
    const double averagine_C = 4.9384;
    const double p = 0.0107;

    averagine_table.assign(averagine_max_mass_bin_index + 1, vector<double>(10, 0.0));
    averagine_max_pos.assign(averagine_max_mass_bin_index + 1, 0);

    for (int bin_idx = 0; bin_idx <= averagine_max_mass_bin_index; ++bin_idx)
    {
      double neutral_mass_bin = static_cast<double>(bin_idx * averagine_mass_bin_step);
      int n_C = static_cast<int>(round(neutral_mass_bin / averagine_mass * averagine_C));
      n_C = max(n_C, 1);

      boost::math::binomial_distribution<double> dist(n_C, p);

      double sum_prob = 0.0;
      for (Size k = 0; k < 10; ++k)
      {
        double prob = (static_cast<int>(k) <= n_C) ?
                      boost::math::pdf(dist, static_cast<unsigned>(k)) : 0.0;
        averagine_table[bin_idx][k] = prob;
        sum_prob += prob;
      }

      if (sum_prob <= 0.0) sum_prob = 1.0;
      for (Size k = 0; k < 10; ++k)
      {
        averagine_table[bin_idx][k] /= sum_prob;
      }

      Size max_pos = distance(averagine_table[bin_idx].begin(),
                              max_element(averagine_table[bin_idx].begin(),
                                          averagine_table[bin_idx].end()));
      if (max_pos < 4) max_pos = 4;
      averagine_max_pos[bin_idx] = max_pos;
    }

    averagine_initialized = true;
  }

  // Map neutral mass onto the same 100 Da grid as the Python implementation.
  int bin_idx = static_cast<int>(floor(neutral_mass / static_cast<double>(averagine_mass_bin_step)));
  if (bin_idx < 0)
  {
    bin_idx = 0;
  }
  if (bin_idx > averagine_max_mass_bin_index)
  {
    bin_idx = averagine_max_mass_bin_index;
  }

  vector<double> theor(10, 0.0);
  const vector<double>& probs = averagine_table[bin_idx];

  // Scale normalized probabilities to the mono apex intensity,
  // preserving the same relative shape as in the Python code.
  double mono_prob = probs[0];
  if (mono_prob <= 0.0)
  {
    mono_prob = 1.0;
  }
  for (Size k = 0; k < 10; ++k)
  {
    theor[k] = apex_intensity * probs[k] / mono_prob;
  }

  Size max_pos = averagine_max_pos[bin_idx];
  return make_pair(theor, max_pos);
}

double Biosaur2Algorithm::computeHillMzStep_(const MSExperiment& exp,
                                             double htol_ppm,
                                             double min_intensity,
                                             double min_mz,
                                             double max_mz) const
{
  double max_mz_value = 0.0;
  for (Size scan_idx = 0; scan_idx < exp.size(); ++scan_idx)
  {
    const MSSpectrum& spectrum = exp[scan_idx];
    for (Size peak_idx = 0; peak_idx < spectrum.size(); ++peak_idx)
    {
      const Peak1D& peak = spectrum[peak_idx];
      double mz = peak.getMZ();
      double intensity = peak.getIntensity();
      if (intensity < min_intensity || mz < min_mz || mz > max_mz)
      {
        continue;
      }
      if (mz > max_mz_value)
      {
        max_mz_value = mz;
      }
    }
  }

  if (max_mz_value <= 0.0 || htol_ppm <= 0.0)
  {
    return 0.0;
  }

  return htol_ppm * 1e-6 * max_mz_value;
}

void Biosaur2Algorithm::processTOF_(MSExperiment& exp) const
{
  OPENMS_LOG_INFO << "Applying TOF-specific intensity filtering..." << endl;

  const double mz_bin_size = 50.0;
  map<int, vector<double>> intensity_bins;

  Size sample_size = min(Size(25), exp.size());
  for (Size i = 0; i < sample_size; ++i)
  {
    for (Size j = 0; j < exp[i].size(); ++j)
    {
      double mz = exp[i][j].getMZ();
      double intensity = exp[i][j].getIntensity();
      if (mz >= minmz_ && mz <= maxmz_)
      {
        int bin = static_cast<int>(mz / mz_bin_size);
        if (intensity <= 0.0)
        {
          continue;
        }
        double log_intensity = log10(intensity);
        if (!std::isfinite(log_intensity))
        {
          continue;
        }
        intensity_bins[bin].push_back(log_intensity);
      }
    }
  }

  map<int, double> bin_thresholds;
  for (auto& bin_pair : intensity_bins)
  {
    vector<double> finite_intensities;
    finite_intensities.reserve(bin_pair.second.size());
    for (double value : bin_pair.second)
    {
      if (std::isfinite(value))
      {
        finite_intensities.push_back(value);
      }
    }

    if (finite_intensities.size() >= 150)
    {
      double sum = accumulate(finite_intensities.begin(), finite_intensities.end(), 0.0);
      double mean = sum / finite_intensities.size();

      double sq_sum = 0.0;
      for (double val : finite_intensities)
      {
        sq_sum += (val - mean) * (val - mean);
      }
      double std_dev = sqrt(sq_sum / finite_intensities.size());

      if (!std::isfinite(mean) || !std::isfinite(std_dev))
      {
        continue;
      }

      bin_thresholds[bin_pair.first] = pow(10.0, mean + 2.0 * std_dev);
    }
  }

  Size total_peaks_before = 0;
  Size total_peaks_after = 0;

  for (auto& spectrum : exp)
  {
    total_peaks_before += spectrum.size();

    MSSpectrum filtered_spectrum;
    filtered_spectrum.setRT(spectrum.getRT());
    filtered_spectrum.setMSLevel(spectrum.getMSLevel());

    for (Size i = 0; i < spectrum.size(); ++i)
    {
      double mz = spectrum[i].getMZ();
      double intensity = spectrum[i].getIntensity();
      int bin = static_cast<int>(mz / mz_bin_size);

      double threshold = 150.0;
      if (bin_thresholds.find(bin) != bin_thresholds.end())
      {
        threshold = bin_thresholds[bin];
      }

      if (intensity >= threshold)
      {
        filtered_spectrum.push_back(spectrum[i]);
      }
    }

    filtered_spectrum.getFloatDataArrays() = spectrum.getFloatDataArrays();
    filtered_spectrum.getIntegerDataArrays() = spectrum.getIntegerDataArrays();
    filtered_spectrum.getStringDataArrays() = spectrum.getStringDataArrays();
    spectrum = filtered_spectrum;
    total_peaks_after += spectrum.size();
  }

  OPENMS_LOG_INFO << "TOF filtering: " << total_peaks_before
                  << " peaks -> " << total_peaks_after << " peaks" << endl;
}

void Biosaur2Algorithm::centroidProfileSpectra_(MSExperiment& exp) const
{
  OPENMS_LOG_INFO << "Centroiding profile spectra using PeakPickerHiRes..." << endl;

  PeakPickerHiRes picker;
  Param picker_param = picker.getParameters();
  picker_param.setValue("signal_to_noise", 0.0);
  picker.setParameters(picker_param);

  MSExperiment centroided_exp;
  Size total_peaks_before = 0;
  Size total_peaks_after = 0;

  for (Size i = 0; i < exp.size(); ++i)
  {
    total_peaks_before += exp[i].size();
    MSSpectrum centroided_spectrum;

    if (exp[i].getType() == SpectrumSettings::CENTROID)
    {
      centroided_spectrum = exp[i];
    }
    else
    {
      picker.pick(exp[i], centroided_spectrum);
      centroided_spectrum.setRT(exp[i].getRT());
      centroided_spectrum.setMSLevel(exp[i].getMSLevel());
      centroided_spectrum.setType(SpectrumSettings::CENTROID);
      if (exp[i].getDriftTime() >= 0)
      {
        centroided_spectrum.setDriftTime(exp[i].getDriftTime());
      }
    }

    centroided_exp.addSpectrum(centroided_spectrum);
    total_peaks_after += centroided_spectrum.size();
  }

  exp = centroided_exp;
  OPENMS_LOG_INFO << "Centroiding: " << total_peaks_before
                  << " profile points -> " << total_peaks_after << " centroided peaks" << endl;
}

void Biosaur2Algorithm::centroidPASEFData_(MSExperiment& exp, double mz_step, double pasef_tolerance) const
{
  if (mz_step <= 0.0 || pasef_tolerance <= 0.0)
  {
    return;
  }

  const double hill_mz_accuracy = htol_;
  const double ion_mobility_accuracy = pasef_tolerance;

  Size total_peaks_before = 0;
  Size total_peaks_after = 0;

  auto centroid_one_spectrum = [&](MSSpectrum& spectrum)
  {
    total_peaks_before += spectrum.size();

    const auto& fda = spectrum.getFloatDataArrays();
    auto it_im = getDataArrayByName(fda, Constants::UserParam::ION_MOBILITY);
    if (it_im == fda.end() || it_im->size() != spectrum.size())
    {
      // No usable per-peak ion-mobility data for this spectrum; leave it as-is.
      return;
    }

    const Size n_peaks = spectrum.size();
    if (n_peaks == 0)
    {
      return;
    }

    vector<double> mz_ar;
    vector<double> intensity_ar;
    vector<double> im_ar;
    mz_ar.reserve(n_peaks);
    intensity_ar.reserve(n_peaks);
    im_ar.reserve(n_peaks);

    for (Size i = 0; i < n_peaks; ++i)
    {
      double mz = spectrum[i].getMZ();
      double intensity = spectrum[i].getIntensity();
      if (intensity < mini_ || mz < minmz_ || mz > maxmz_)
      {
        continue;
      }
      double im = (*it_im)[i];
      mz_ar.push_back(mz);
      intensity_ar.push_back(intensity);
      im_ar.push_back(im);
    }

    if (mz_ar.empty())
    {
      spectrum.clear(false);
      return;
    }

    auto it_max_im = max_element(im_ar.begin(), im_ar.end());
    if (it_max_im == im_ar.end() || *it_max_im <= 0.0)
    {
      return;
    }

    const double ion_mobility_step = (*it_max_im) * ion_mobility_accuracy;
    if (ion_mobility_step <= 0.0)
    {
      return;
    }

    const Size n = mz_ar.size();
    vector<int> mz_fast(n);
    vector<int> im_fast(n);
    for (Size i = 0; i < n; ++i)
    {
      mz_fast[i] = static_cast<int>(mz_ar[i] / mz_step);
      im_fast[i] = static_cast<int>(im_ar[i] / ion_mobility_step);
    }

    // Sort by coarse m/z bin index, mirroring centroid_pasef_scan.
    vector<Size> order(n);
    iota(order.begin(), order.end(), 0);
    sort(order.begin(), order.end(),
         [&](Size a, Size b) { return mz_fast[a] < mz_fast[b]; });

    vector<double> mz_sorted(n);
    vector<double> intensity_sorted(n);
    vector<double> im_sorted(n);
    vector<int> mz_fast_sorted(n);
    vector<int> im_fast_sorted(n);
    for (Size pos = 0; pos < n; ++pos)
    {
      Size idx = order[pos];
      mz_sorted[pos] = mz_ar[idx];
      intensity_sorted[pos] = intensity_ar[idx];
      im_sorted[pos] = im_ar[idx];
      mz_fast_sorted[pos] = mz_fast[idx];
      im_fast_sorted[pos] = im_fast[idx];
    }

    vector<bool> banned(n, false);
    vector<double> mz_new;
    vector<double> intensity_new;
    vector<double> im_new;

    Size peak_idx = 0;
    const Size max_peak_idx = n;

    while (peak_idx < max_peak_idx)
    {
      vector<Size> tmp;

      if (!banned[peak_idx])
      {
        const double mass_accuracy_cur = mz_sorted[peak_idx] * 1e-6 * hill_mz_accuracy;
        const int mz_val_int = mz_fast_sorted[peak_idx];
        const int im_val_int = im_fast_sorted[peak_idx];

        tmp.push_back(peak_idx);

        Size peak_idx_2 = peak_idx + 1;
        while (peak_idx_2 < max_peak_idx)
        {
          if (!banned[peak_idx_2])
          {
            const int mz_val_int_2 = mz_fast_sorted[peak_idx_2];
            if (mz_val_int_2 - mz_val_int > 1)
            {
              break;
            }

            if (fabs(mz_sorted[peak_idx] - mz_sorted[peak_idx_2]) <= mass_accuracy_cur)
            {
              const int im_val_int_2 = im_fast_sorted[peak_idx_2];
              if (abs(im_val_int - im_val_int_2) <= 1)
              {
                if (fabs(im_sorted[peak_idx] - im_sorted[peak_idx_2]) <= ion_mobility_accuracy)
                {
                  tmp.push_back(peak_idx_2);
                  peak_idx = peak_idx_2;
                }
              }
            }
          }
          ++peak_idx_2;
        }
      }

      const Size l_new = tmp.size();
      if (l_new >= pasefminlh_)
      {
        if (l_new == 1)
        {
          const double i_val_new = intensity_sorted[peak_idx];
          if (i_val_new >= pasefmini_)
          {
            mz_new.push_back(mz_sorted[peak_idx]);
            intensity_new.push_back(i_val_new);
            im_new.push_back(im_sorted[peak_idx]);
            banned[peak_idx] = true;
          }
        }
        else
        {
          double i_val_new = 0.0;
          vector<double> all_mz;
          vector<double> all_im;
          all_mz.reserve(l_new);
          all_im.reserve(l_new);

          for (Size idx : tmp)
          {
            const double intensity = intensity_sorted[idx];
            i_val_new += intensity;
            all_mz.push_back(mz_sorted[idx]);
            all_im.push_back(im_sorted[idx]);
          }

          if (i_val_new >= pasefmini_)
          {
            double mz_weighted_sum = 0.0;
            double im_weighted_sum = 0.0;
            for (Size k = 0; k < l_new; ++k)
            {
              mz_weighted_sum += all_mz[k] * intensity_sorted[tmp[k]];
              im_weighted_sum += all_im[k] * intensity_sorted[tmp[k]];
            }
            const double mz_val_new = mz_weighted_sum / i_val_new;
            const double im_val_new = im_weighted_sum / i_val_new;

            mz_new.push_back(mz_val_new);
            intensity_new.push_back(i_val_new);
            im_new.push_back(im_val_new);

            for (Size idx : tmp)
            {
              banned[idx] = true;
            }
          }
        }
      }

      ++peak_idx;
    }

    // Replace spectrum peaks and ion-mobility array with centroided values.
    spectrum.clear(false);
    for (Size i = 0; i < mz_new.size(); ++i)
    {
      Peak1D peak;
      peak.setMZ(mz_new[i]);
      peak.setIntensity(intensity_new[i]);
      spectrum.push_back(peak);
    }

    MSSpectrum::FloatDataArrays new_fda;
    MSSpectrum::FloatDataArray im_array;
    im_array.setName(Constants::UserParam::ION_MOBILITY);
    im_array.assign(im_new.begin(), im_new.end());
    new_fda.push_back(im_array);
    spectrum.setFloatDataArrays(new_fda);

    total_peaks_after += spectrum.size();
  };

  for (auto& spectrum : exp)
  {
    if (spectrum.getMSLevel() != 1)
    {
      continue;
    }
    centroid_one_spectrum(spectrum);
  }

  OPENMS_LOG_INFO << "PASEF centroiding: " << total_peaks_before
                  << " peaks -> " << total_peaks_after << " centroided clusters" << endl;
}

void Biosaur2Algorithm::processFAIMSGroup_(double faims_cv,
                                           MSExperiment& group_exp,
                                           double original_paseftol,
                                           vector<Hill>& hills_out,
                                           vector<PeptideFeature>& features_out)
{
  if (group_exp.empty())
  {
    return;
  }

  // Determine ion-mobility array availability for this group.
  bool any_im_array = false;
  bool any_missing_im = false;
  for (const auto& spec : group_exp)
  {
    const auto& fda = spec.getFloatDataArrays();
    auto it_im = getDataArrayByName(fda, Constants::UserParam::ION_MOBILITY);
    if (it_im != fda.end() && !it_im->empty())
    {
      any_im_array = true;
    }
    else
    {
      any_missing_im = true;
    }
  }

  // Compute mz_step for PASEF centroiding, mirroring Python's use of
  // hill mass accuracy and the maximum m/z in the group.
  double mz_step = 0.0;
  double max_mz_value = 0.0;
  for (const auto& spectrum : group_exp)
  {
    for (Size peak_idx = 0; peak_idx < spectrum.size(); ++peak_idx)
    {
      const Peak1D& peak = spectrum[peak_idx];
      double mz = peak.getMZ();
      double intensity = peak.getIntensity();
      if (intensity < mini_ || mz < minmz_ || mz > maxmz_)
      {
        continue;
      }
      if (mz > max_mz_value)
      {
        max_mz_value = mz;
      }
    }
  }
  if (max_mz_value > 0.0 && htol_ > 0.0)
  {
    mz_step = htol_ * 1e-6 * max_mz_value;
  }

  // Decide whether to use IM-based centroiding and gating for this group.
  const bool use_im_group = (any_im_array && !any_missing_im && original_paseftol > 0.0);
  if (use_im_group && mz_step > 0.0)
  {
    OPENMS_LOG_INFO << "Applying PASEF/TIMS centroiding for group (FAIMS CV="
                    << faims_cv << ") with paseftol=" << original_paseftol << endl;
    StopWatch stage_timer;
    stage_timer.start();
    centroidPASEFData_(group_exp, mz_step, original_paseftol);
    stage_timer.stop();
    OPENMS_LOG_INFO << "PASEF centroiding for group (FAIMS CV=" << faims_cv
                    << ") took " << stage_timer.toString() << endl;
  }
  else
  {
    // Mirror Python behavior: if ion mobility is not consistently available
    // (or absent entirely), disable IM-based gating.
    if (original_paseftol > 0.0 && any_im_array && any_missing_im)
    {
      OPENMS_LOG_WARN << "Disabling ion-mobility gating for group (FAIMS CV="
                      << faims_cv
                      << ") due to missing/partial IM arrays; proceeding in 1D m/z space."
                      << endl;
    }
  }

  // Hill mass calibration per group, analogous to Python's per-FAIMS pass.
  double calibrated_htol = htol_;
  if (use_hill_calib_)
  {
    OPENMS_LOG_INFO << "Performing hill mass tolerance calibration for group (FAIMS CV="
                    << faims_cv << ")..." << endl;
    StopWatch stage_timer;
    stage_timer.start();
    vector<double> mass_diffs;
    Size sample_size = min(group_exp.size(), Size(1000));
    Size start_idx = (group_exp.size() > 1000) ? (group_exp.size() / 2 - 500) : 0;

    MSExperiment calib_exp;
    for (Size i = start_idx; i < start_idx + sample_size && i < group_exp.size(); ++i)
    {
      calib_exp.addSpectrum(group_exp[i]);
    }

    vector<Hill> calib_hills = detectHills_(calib_exp, htol_, mini_, minmz_, maxmz_, /*use_im*/ false, &mass_diffs);
    (void)calib_hills;
    if (!mass_diffs.empty())
    {
      auto calib = calibrateMass_(mass_diffs);
      double calibrated_sigma = calib.second;
      calibrated_htol = min(htol_, 5.0 * calibrated_sigma);
      OPENMS_LOG_INFO << "Automatically optimized htol parameter for group (FAIMS CV="
                      << faims_cv << "): " << calibrated_htol
                      << " ppm (was " << htol_ << " ppm)" << endl;
    }
    stage_timer.stop();
    OPENMS_LOG_INFO << "Hill calibration for group (FAIMS CV=" << faims_cv
                    << ") took " << stage_timer.toString() << endl;
  }

  // Hill detection and processing for this group.
  StopWatch stage_timer;
  stage_timer.start();
  vector<Hill> group_hills = detectHills_(group_exp, calibrated_htol, mini_, minmz_, maxmz_, use_im_group);
  stage_timer.stop();
  OPENMS_LOG_INFO << "Hill detection for group (FAIMS CV=" << faims_cv
                  << ") found " << group_hills.size()
                  << " hills in " << stage_timer.toString() << endl;
  stage_timer.reset();
  stage_timer.start();
  group_hills = processHills_(group_hills, minlh_);
  stage_timer.stop();
  OPENMS_LOG_INFO << "Hill preprocessing for group (FAIMS CV=" << faims_cv
                  << ") kept " << group_hills.size()
                  << " hills in " << stage_timer.toString() << endl;
  stage_timer.reset();
  stage_timer.start();
  group_hills = splitHills_(group_hills, hvf_, minlh_);
  stage_timer.stop();
  OPENMS_LOG_INFO << "Hill splitting for group (FAIMS CV=" << faims_cv
                  << ") produced " << group_hills.size()
                  << " hills in " << stage_timer.toString() << endl;

  bool enable_isotope_calib = !ignore_iso_calib_;
  stage_timer.reset();
  stage_timer.start();
  vector<PeptideFeature> group_features =
    detectIsotopePatterns_(group_hills, itol_, cmin_, cmax_, negative_mode_, ivf_, iuse_, enable_isotope_calib, use_im_group);
  stage_timer.stop();
  OPENMS_LOG_INFO << "Isotope pattern detection for group (FAIMS CV=" << faims_cv
                  << ") produced " << group_features.size()
                  << " features in " << stage_timer.toString() << endl;

  hills_out = std::move(group_hills);
  features_out = std::move(group_features);
}

void Biosaur2Algorithm::linkScanToHills_(const MSSpectrum& spectrum,
                                         Size scan_idx,
                                         double htol_ppm,
                                         double min_intensity,
                                         double min_mz,
                                         double max_mz,
                                         double mz_step,
                                         bool use_im_global,
                                         vector<Hill>& hills,
                                         Size& hill_idx_counter,
                                         vector<Size>& prev_peak_to_hill,
                                         const MSSpectrum*& prev_spectrum_ptr,
                                         map<int, vector<int>>& prev_fast_dict,
                                         vector<int>& prev_im_bins,
                                         vector<double>* hill_mass_diffs) const
{
  double rt = spectrum.getRT();

  // Collect indices of peaks passing basic filters.
  vector<int> valid_indices;
  valid_indices.reserve(spectrum.size());
  for (Size peak_idx = 0; peak_idx < spectrum.size(); ++peak_idx)
  {
    const Peak1D& peak = spectrum[peak_idx];
    double mz = peak.getMZ();
    double intensity = peak.getIntensity();
    if (intensity < min_intensity || mz < min_mz || mz > max_mz)
    {
      continue;
    }
    valid_indices.push_back(static_cast<int>(peak_idx));
  }

  const Size len_mz = valid_indices.size();
  if (len_mz == 0)
  {
    // Reset linking if the current scan has no usable peaks.
    prev_fast_dict.clear();
    prev_peak_to_hill.clear();
    prev_spectrum_ptr = nullptr;
    return;
  }

  // Build intensity, m/z and optional ion-mobility-bin arrays for valid peaks.
  vector<double> intensities(len_mz);
  vector<double> mzs(len_mz);
  vector<int> im_bin_per_peak(spectrum.size(), 0);

  const auto& fda = spectrum.getFloatDataArrays();
  auto it_im = getDataArrayByName(fda, Constants::UserParam::ION_MOBILITY);
  const bool use_im_current = use_im_global && (it_im != fda.end());

  for (Size i = 0; i < len_mz; ++i)
  {
    const Size peak_idx = static_cast<Size>(valid_indices[i]);
    intensities[i] = spectrum[peak_idx].getIntensity();
    mzs[i] = spectrum[peak_idx].getMZ();
    if (use_im_current)
    {
      if (peak_idx >= it_im->size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Ion mobility array shorter than peak list.",
                                      String(peak_idx));
      }
      const double im = (*it_im)[peak_idx];
      int im_bin = (paseftol_ > 0.0) ? static_cast<int>(im / paseftol_) : 0;
      im_bin_per_peak[peak_idx] = im_bin;
    }
  }

  // Sort peaks in descending intensity (reference implementation behavior).
  vector<int> order(len_mz);
  iota(order.begin(), order.end(), 0);
  sort(order.begin(), order.end(),
       [&](int a, int b) { return intensities[a] > intensities[b]; });

  vector<int> basic_id_sorted(len_mz);
  vector<int> fast_array(len_mz);

  for (Size pos = 0; pos < len_mz; ++pos)
  {
    const int local_idx = order[pos];
    const int orig_idx = valid_indices[static_cast<Size>(local_idx)];
    basic_id_sorted[pos] = orig_idx;
    int fm = 0;
    if (mz_step > 0.0)
    {
      fm = static_cast<int>(mzs[static_cast<Size>(local_idx)] / mz_step);
    }
    fast_array[pos] = fm;
  }

  // Build fast lookup structure in m/z space for the current scan.
  map<int, vector<int>> fast_dict;
  for (Size pos = 0; pos < len_mz; ++pos)
  {
    fast_dict[fast_array[pos]].push_back(basic_id_sorted[pos]);
  }

  // Hill assignments for peaks in the current scan (indexed by peak index).
  vector<Size> curr_peak_to_hill(spectrum.size(),
                                 numeric_limits<Size>::max());
  set<int> banned_prev_idx_set;

  auto append_peak_to_hill = [&](Size hill_id, int peak_idx)
  {
    Hill& hill = hills[hill_id];
    const Size p_idx = static_cast<Size>(peak_idx);

    hill.scan_indices.push_back(scan_idx);
    hill.peak_indices.push_back(p_idx);

    const Peak1D& peak = spectrum[p_idx];
    const double mz = peak.getMZ();
    const double intensity = peak.getIntensity();

    hill.mz_values.push_back(mz);
    hill.intensities.push_back(intensity);
    hill.rt_values.push_back(rt);

    const double drift_time = spectrum.getDriftTime();
    hill.drift_times.push_back(drift_time);

    double ion_mobility = -1.0;
    const auto& fda_local = spectrum.getFloatDataArrays();
    auto it_im_local = getDataArrayByName(fda_local, Constants::UserParam::ION_MOBILITY);
    if (it_im_local != fda_local.end())
    {
      if (p_idx >= it_im_local->size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Ion mobility array shorter than peak list.",
                                      String(p_idx));
      }
      ion_mobility = (*it_im_local)[p_idx];
    }
    else if (drift_time >= 0)
    {
      if (shouldThrowForMissingIM_(spectrum))
      {
        throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            "Ion mobility array missing although drift times are present.");
      }
    }
    hill.ion_mobilities.push_back(ion_mobility);

    hill.length = hill.scan_indices.size();
    hill.rt_end = rt;
    if (hill.length == 1)
    {
      hill.rt_start = rt;
      hill.rt_apex = rt;
      hill.intensity_apex = intensity;
    }
    else if (intensity > hill.intensity_apex)
    {
      hill.intensity_apex = intensity;
      hill.rt_apex = rt;
    }
    hill.intensity_sum += intensity;

    curr_peak_to_hill[p_idx] = hill_id;
  };

  for (Size pos = 0; pos < len_mz; ++pos)
  {
    const int idx = basic_id_sorted[pos];
    const int fm = fast_array[pos];
    const int fi = use_im_current ? im_bin_per_peak[static_cast<Size>(idx)] : 0;

    // Collect candidate previous-scan peaks from neighboring m/z bins.
    bool flag1 = prev_fast_dict.find(fm) != prev_fast_dict.end();
    bool flag2 = prev_fast_dict.find(fm - 1) != prev_fast_dict.end();
    bool flag3 = prev_fast_dict.find(fm + 1) != prev_fast_dict.end();

    Size assigned_hill = numeric_limits<Size>::max();

    if (flag1 || flag2 || flag3)
    {
      vector<int> all_idx;
      if (flag1)
      {
        const auto& v = prev_fast_dict[fm];
        all_idx.insert(all_idx.end(), v.begin(), v.end());
        if (flag2)
        {
          const auto& v2 = prev_fast_dict[fm - 1];
          all_idx.insert(all_idx.end(), v2.begin(), v2.end());
        }
        if (flag3)
        {
          const auto& v3 = prev_fast_dict[fm + 1];
          all_idx.insert(all_idx.end(), v3.begin(), v3.end());
        }
      }
      else if (flag2)
      {
        const auto& v = prev_fast_dict[fm - 1];
        all_idx.insert(all_idx.end(), v.begin(), v.end());
        if (flag3)
        {
          const auto& v3 = prev_fast_dict[fm + 1];
          all_idx.insert(all_idx.end(), v3.begin(), v3.end());
        }
      }
      else if (flag3)
      {
        const auto& v = prev_fast_dict[fm + 1];
        all_idx.insert(all_idx.end(), v.begin(), v.end());
      }

      // Optional ion-mobility gating for previous-scan candidates.
      if (use_im_current && !prev_im_bins.empty())
      {
        vector<int> filtered;
        filtered.reserve(all_idx.size());
        for (int idx_prev : all_idx)
        {
          if (idx_prev < 0 || static_cast<Size>(idx_prev) >= prev_im_bins.size())
          {
            continue;
          }
          int prev_bin = prev_im_bins[static_cast<Size>(idx_prev)];
          if (abs(prev_bin - fi) <= 1)
          {
            filtered.push_back(idx_prev);
          }
        }
        all_idx.swap(filtered);
      }

      double best_intensity = 0.0;
      int best_idx_prev = -1;
      double best_mass_diff_with_sign = 0.0;

      const double mz_cur = spectrum[static_cast<Size>(idx)].getMZ();

      if (prev_spectrum_ptr != nullptr)
      {
        const MSSpectrum& prev_spectrum = *prev_spectrum_ptr;
        for (int idx_prev : all_idx)
        {
          if (idx_prev < 0 || static_cast<Size>(idx_prev) >= prev_spectrum.size())
          {
            continue;
          }
          if (banned_prev_idx_set.find(idx_prev) != banned_prev_idx_set.end())
          {
            continue;
          }

          const Peak1D& prev_peak = prev_spectrum[static_cast<Size>(idx_prev)];
          const double cur_intensity = prev_peak.getIntensity();
          const double mz_prev = prev_peak.getMZ();

          const double cur_mass_diff_with_sign =
            (mz_cur - mz_prev) / mz_cur * 1e6;
          const double cur_mass_diff = fabs(cur_mass_diff_with_sign);

          if (cur_mass_diff <= htol_ppm && cur_intensity >= best_intensity)
          {
            best_intensity = cur_intensity;
            best_idx_prev = idx_prev;
            best_mass_diff_with_sign = cur_mass_diff_with_sign;
          }
        }
      }

      if (best_idx_prev >= 0)
      {
        banned_prev_idx_set.insert(best_idx_prev);
        if (hill_mass_diffs != nullptr)
        {
          hill_mass_diffs->push_back(best_mass_diff_with_sign);
        }

        if (!prev_peak_to_hill.empty() &&
            static_cast<Size>(best_idx_prev) < prev_peak_to_hill.size())
        {
          assigned_hill = prev_peak_to_hill[static_cast<Size>(best_idx_prev)];
        }
      }
    }

    // Attach to existing hill if a suitable previous peak was found.
    if (assigned_hill != numeric_limits<Size>::max() &&
        assigned_hill < hills.size())
    {
      append_peak_to_hill(assigned_hill, idx);
    }
    else
    {
      // Start a new hill for this peak.
      Hill hill;
      hill.hill_idx = hill_idx_counter++;
      hills.push_back(hill);
      append_peak_to_hill(hill.hill_idx, idx);
    }
  }

  // Prepare state for linking from this scan to the next.
  prev_fast_dict = std::move(fast_dict);
  prev_peak_to_hill.assign(spectrum.size(), numeric_limits<Size>::max());
  prev_im_bins.assign(spectrum.size(), 0);
  for (Size i = 0; i < spectrum.size(); ++i)
  {
    prev_peak_to_hill[i] = curr_peak_to_hill[i];
    prev_im_bins[i] = im_bin_per_peak[i];
  }
  prev_spectrum_ptr = &spectrum;
}

vector<Biosaur2Algorithm::Hill> Biosaur2Algorithm::detectHills_(const MSExperiment& exp,
                                                                double htol_ppm,
                                                                double min_intensity,
                                                                double min_mz,
                                                                double max_mz,
                                                                bool use_im,
                                                                vector<double>* hill_mass_diffs) const
{
  vector<Hill> hills;
  Size hill_idx_counter = 0;

  // Bin width in m/z space, mirroring the Python reference implementation.
  const double mz_step = computeHillMzStep_(exp, htol_ppm, min_intensity, min_mz, max_mz);

  if (hill_mass_diffs != nullptr)
  {
    hill_mass_diffs->clear();
  }

  // Mapping from previous-scan peak index to hill index and ion-mobility bin.
  vector<Size> prev_peak_to_hill;
  const MSSpectrum* prev_spectrum_ptr = nullptr;
  map<int, vector<int>> prev_fast_dict;
  vector<int> prev_im_bins;
  const bool use_im_global = use_im;

  for (Size scan_idx = 0; scan_idx < exp.size(); ++scan_idx)
  {
    const MSSpectrum& spectrum = exp[scan_idx];
    linkScanToHills_(spectrum,
                     scan_idx,
                     htol_ppm,
                     min_intensity,
                     min_mz,
                     max_mz,
                     mz_step,
                     use_im_global,
                     hills,
                     hill_idx_counter,
                     prev_peak_to_hill,
                     prev_spectrum_ptr,
                     prev_fast_dict,
                     prev_im_bins,
                     hill_mass_diffs);
  }

  return hills;
}

vector<Biosaur2Algorithm::Hill> Biosaur2Algorithm::processHills_(const vector<Hill>& hills, Size min_length) const
{
  vector<Hill> processed;
  for (const auto& hill : hills)
  {
    if (hill.length >= min_length)
    {
      Hill processed_hill = hill;
      processed_hill.drift_time_median = calculateMedian_(hill.drift_times);

      // Recompute hill center m/z and ion mobility using intensity-weighted
      // averages to mirror the reference Biosaur2 implementation.
      if (!hill.mz_values.empty() && hill.mz_values.size() == hill.intensities.size())
      {
        double weighted_mz_sum = 0.0;
        double intensity_sum = 0.0;
        for (Size i = 0; i < hill.mz_values.size(); ++i)
        {
          const double intensity = hill.intensities[i];
          weighted_mz_sum += hill.mz_values[i] * intensity;
          intensity_sum += intensity;
        }
        if (intensity_sum > 0.0)
        {
          processed_hill.mz_weighted_mean = weighted_mz_sum / intensity_sum;
        }
      }

      if (!hill.ion_mobilities.empty() && hill.ion_mobilities.size() == hill.intensities.size())
      {
        double weighted_im_sum = 0.0;
        double intensity_sum_im = 0.0;
        for (Size i = 0; i < hill.ion_mobilities.size(); ++i)
        {
          const double intensity = hill.intensities[i];
          weighted_im_sum += hill.ion_mobilities[i] * intensity;
          intensity_sum_im += intensity;
        }
        if (intensity_sum_im > 0.0)
        {
          processed_hill.ion_mobility_median = weighted_im_sum / intensity_sum_im;
        }
        else
        {
          processed_hill.ion_mobility_median = calculateMedian_(hill.ion_mobilities);
        }
      }
      else
      {
        processed_hill.ion_mobility_median = calculateMedian_(hill.ion_mobilities);
      }
      processed.push_back(processed_hill);
    }
  }
  return processed;
}

vector<Biosaur2Algorithm::Hill> Biosaur2Algorithm::splitHills_(const vector<Hill>& hills, double hvf, Size min_length) const
{
  vector<Hill> split_hills;

  // Determine the next free hill index so that new split segments
  // receive unique IDs, similar to the Python implementation where
  // newly created segments get fresh hill indices.
  Size next_hill_idx = 0;
  for (const auto& h : hills)
  {
    next_hill_idx = max(next_hill_idx, h.hill_idx + 1);
  }

  for (const auto& hill : hills)
  {
    // Only attempt splitting for sufficiently long hills, mirroring the
    // reference implementation (length >= 2 * min_length).
    if (hill.length < 2 * min_length)
    {
      split_hills.push_back(hill);
      continue;
    }

    const Size hill_len = hill.length;
    vector<double> smoothed = meanFilter_(hill.intensities, 3);

    // First pass: identify candidate valley indices (min_idx_list) and
    // corresponding recheck positions, closely following split_peaks.
    vector<Size> min_idx_list;
    vector<Size> recheck_positions;

    const int min_len = static_cast<int>(min_length);
    const int c_len = static_cast<int>(hill_len) - min_len;
    int idx = min_len - 1;
    Size l_idx = 0;
    double min_val = 0.0;

    while (idx <= c_len)
    {
      if (!min_idx_list.empty() &&
          static_cast<Size>(idx) >= min_idx_list.back() + min_length)
      {
        l_idx = min_idx_list.back();
      }

      // Left and right maxima around the candidate valley.
      double valley_intensity = smoothed[static_cast<Size>(idx)];
      if (valley_intensity <= 0.0)
      {
        ++idx;
        continue;
      }

      double left_max = 0.0;
      for (Size j = l_idx; j < static_cast<Size>(idx); ++j)
      {
        left_max = max(left_max, smoothed[j]);
      }
      if (left_max == 0.0)
      {
        ++idx;
        continue;
      }
      double l_r = left_max / valley_intensity;
      if (l_r >= hvf)
      {
        double right_max = 0.0;
        for (Size j = static_cast<Size>(idx) + 1; j < hill_len; ++j)
        {
          right_max = max(right_max, smoothed[j]);
        }
        if (right_max > 0.0)
        {
          double r_r = right_max / valley_intensity;
          if (r_r >= hvf)
          {
            double mult_val = l_r * r_r;
            int include_factor = (l_r > r_r) ? 1 : 0;
            int candidate_pos_int = idx + include_factor;
            if (min_len <= candidate_pos_int && candidate_pos_int <= c_len)
            {
              Size candidate_pos = static_cast<Size>(candidate_pos_int);
              if (min_idx_list.empty() ||
                  candidate_pos >= min_idx_list.back() + min_length)
              {
                min_idx_list.push_back(candidate_pos);
                recheck_positions.push_back(static_cast<Size>(idx));
                min_val = mult_val;
              }
              else if (mult_val > min_val)
              {
                min_idx_list.back() = candidate_pos;
                recheck_positions.back() = static_cast<Size>(idx);
                min_val = mult_val;
              }
            }
          }
        }
      }
      ++idx;
    }

    // Second pass: recheck right-hand maxima for each candidate valley.
    vector<Size> final_splits;
    if (!min_idx_list.empty())
    {
      for (Size k = 0; k < min_idx_list.size(); ++k)
      {
        Size min_idx = min_idx_list[k];
        Size recheck_idx = recheck_positions[k];
        Size end_idx = (k + 1 < min_idx_list.size()) ? min_idx_list[k + 1]
                                                     : hill_len;

        if (recheck_idx + 1 >= end_idx || recheck_idx >= hill_len)
        {
          continue;
        }

        double recheck_val = smoothed[recheck_idx];
        if (recheck_val <= 0.0)
        {
          continue;
        }

        double right_max = 0.0;
        for (Size j = recheck_idx + 1; j < end_idx; ++j)
        {
          right_max = max(right_max, smoothed[j]);
        }
        if (right_max / recheck_val >= hvf)
        {
          final_splits.push_back(min_idx);
        }
      }
    }

    if (final_splits.empty())
    {
      // No valleys detected: keep original hill unchanged.
      split_hills.push_back(hill);
      continue;
    }

    // Construct segments defined by split positions.
    vector<Size> boundaries;
    boundaries.push_back(0);
    for (Size pos : final_splits)
    {
      boundaries.push_back(pos);
    }
    boundaries.push_back(hill_len);

    for (Size s = 0; s + 1 < boundaries.size(); ++s)
    {
      Size seg_start = boundaries[s];
      Size seg_end = boundaries[s + 1];
      if (seg_end <= seg_start) continue;

      Size seg_len = seg_end - seg_start;
      if (seg_len < min_length) continue;

      Hill new_hill = hill;
      new_hill.scan_indices.clear();
      new_hill.peak_indices.clear();
      new_hill.mz_values.clear();
      new_hill.intensities.clear();
      new_hill.rt_values.clear();
      new_hill.drift_times.clear();
      new_hill.ion_mobilities.clear();

      for (Size k = seg_start; k < seg_end; ++k)
      {
        new_hill.scan_indices.push_back(hill.scan_indices.at(k));
        new_hill.peak_indices.push_back(hill.peak_indices.at(k));
        new_hill.mz_values.push_back(hill.mz_values.at(k));
        new_hill.intensities.push_back(hill.intensities.at(k));
        new_hill.rt_values.push_back(hill.rt_values.at(k));
        if (k < hill.drift_times.size())
        {
          new_hill.drift_times.push_back(hill.drift_times.at(k));
        }
        if (k < hill.ion_mobilities.size())
        {
          new_hill.ion_mobilities.push_back(hill.ion_mobilities.at(k));
        }
      }

      if (new_hill.drift_times.size() != new_hill.scan_indices.size() ||
          new_hill.ion_mobilities.size() != new_hill.scan_indices.size())
      {
        throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                      "Split hill meta data arrays are inconsistent.",
                                      String(new_hill.scan_indices.size()));
      }

	      new_hill.length = new_hill.scan_indices.size();
	      new_hill.rt_start = new_hill.rt_values.front();
	      new_hill.rt_end = new_hill.rt_values.back();
	      // Recompute hill center m/z for the split segment using
	      // an intensity-weighted mean to mirror processHills_ and
	      // the reference Biosaur2 implementation.
	      double weighted_mz_sum = 0.0;
	      double intensity_sum = 0.0;
	      for (Size idx = 0; idx < new_hill.mz_values.size(); ++idx)
	      {
	        const double intensity = new_hill.intensities[idx];
	        weighted_mz_sum += new_hill.mz_values[idx] * intensity;
	        intensity_sum += intensity;
	      }
	      if (intensity_sum > 0.0)
	      {
	        new_hill.mz_weighted_mean = weighted_mz_sum / intensity_sum;
	      }
	      new_hill.intensity_sum = intensity_sum;
	      auto max_it = max_element(new_hill.intensities.begin(), new_hill.intensities.end());
	      Size apex_idx = distance(new_hill.intensities.begin(), max_it);
	      new_hill.intensity_apex = *max_it;
	      new_hill.rt_apex = new_hill.rt_values[apex_idx];

	      // Assign hill indices: keep the original ID for the first
	      // segment and assign fresh IDs to subsequent segments,
	      // mirroring the Python split_peaks behavior.
	      if (s == 0)
	      {
	        new_hill.hill_idx = hill.hill_idx;
	      }
	      else
	      {
	        new_hill.hill_idx = next_hill_idx++;
	      }

      split_hills.push_back(new_hill);
    }
  }

  return split_hills;
}

Size Biosaur2Algorithm::checkIsotopeValleySplit_(const vector<IsotopeCandidate>& isotopes,
                                                 const vector<Hill>& hills,
                                                 double ivf) const
{
  if (isotopes.size() <= 1)
  {
    return isotopes.size();
  }

  vector<double> isotope_intensities;
  isotope_intensities.reserve(isotopes.size());
  for (const auto& iso : isotopes)
  {
    const Hill& hill = *find_if(hills.begin(), hills.end(),
                                [&iso](const Hill& h) { return h.hill_idx == iso.hill_idx; });
    isotope_intensities.push_back(hill.intensity_apex);
  }

  vector<double> smoothed = meanFilter_(isotope_intensities, 3);
  Size max_pos = distance(smoothed.begin(), max_element(smoothed.begin(), smoothed.end()));
  Size min_check_pos = max(Size(4), max_pos + 1);

  for (Size i = min_check_pos; i < smoothed.size() - 1; ++i)
  {
    double local_min = smoothed[i];
    double right_max = *max_element(smoothed.begin() + i + 1, smoothed.end());
    if (local_min * ivf < right_max)
    {
      return i;
    }
  }

  return smoothed.size();
}

vector<Biosaur2Algorithm::PeptideFeature> Biosaur2Algorithm::detectIsotopePatterns_(vector<Hill>& hills,
                                                                                    double itol_ppm,
                                                                                    int min_charge,
                                                                                    int max_charge,
                                                                                    bool negative_mode,
                                                                                    double ivf,
                                                                                    int iuse,
                                                                                    bool enable_isotope_calib,
                                                                                    bool use_im) const
{
  vector<PeptideFeature> features;
  set<Size> used_hills;

  OPENMS_LOG_INFO << "Detecting isotope patterns..." << endl;

  sort(hills.begin(), hills.end(),
       [](const Hill& a, const Hill& b) { return a.mz_weighted_mean < b.mz_weighted_mean; });
  const double ISOTOPE_MASSDIFF = Constants::C13C12_MASSDIFF_U;

  map<int, pair<double, double>> isotope_calib_map;
  for (int ic = 1; ic <= 9; ++ic)
  {
    isotope_calib_map[ic] = make_pair(0.0, itol_ppm);
  }

  if (enable_isotope_calib)
  {
    OPENMS_LOG_INFO << "Performing isotope calibration..." << endl;
    map<int, vector<double>> isotope_errors;
    for (int ic = 1; ic <= 9; ++ic)
    {
      isotope_errors[ic] = vector<double>();
    }

    for (Size i = 0; i < hills.size(); ++i)
    {
      const Hill& mono_hill = hills[i];
      double mono_mz = mono_hill.mz_weighted_mean;

      for (int charge = max_charge; charge >= min_charge; --charge)
      {
        double mz_spacing = ISOTOPE_MASSDIFF / charge;
        bool found_first = false;

        for (int iso_num = 1; iso_num <= 9; ++iso_num)
        {
          double expected_mz = mono_mz + iso_num * mz_spacing;
          double mz_tolerance = expected_mz * itol_ppm * 1e-6;

          for (Size j = i + 1; j < hills.size(); ++j)
          {
            if (hills[j].mz_weighted_mean > expected_mz + mz_tolerance)
            {
              break;
            }

            double diff = fabs(hills[j].mz_weighted_mean - expected_mz);
            if (diff <= mz_tolerance)
            {
              if (mono_hill.length >= 3)
              {
                double mass_diff_ppm = calculatePPM_(hills[j].mz_weighted_mean, expected_mz);
                isotope_errors[iso_num].push_back(mass_diff_ppm);
                if (iso_num == 1)
                {
                  found_first = true;
                }
              }
              break;
            }
          }
        }

        if (found_first)
        {
          break;
        }
      }
    }

    for (int ic = 1; ic <= 3; ++ic)
    {
      if (isotope_errors[ic].size() >= 1000)
      {
        auto calib = calibrateMass_(isotope_errors[ic]);
        isotope_calib_map[ic] = calib;
        OPENMS_LOG_INFO << "Isotope " << ic << " calibration: shift="
                        << calib.first << " ppm, sigma=" << calib.second << " ppm" << endl;
      }
    }

    for (int ic = 4; ic <= 9; ++ic)
    {
      if (isotope_errors[ic].size() >= 1000)
      {
        isotope_calib_map[ic] = calibrateMass_(isotope_errors[ic]);
      }
      else if (ic > 1 && isotope_calib_map.find(ic - 1) != isotope_calib_map.end())
      {
        auto prev = isotope_calib_map[ic - 1];
        auto prev2 = isotope_calib_map.find(ic - 2) != isotope_calib_map.end() ?
                     isotope_calib_map[ic - 2] : make_pair(0.0, itol_ppm);

        double shift_delta = prev.first - prev2.first;
        double sigma_ratio = prev.second / max(prev2.second, 0.1);
        isotope_calib_map[ic] = make_pair(prev.first + shift_delta, prev.second * sigma_ratio);
      }
    }

    OPENMS_LOG_INFO << "Isotope 1 calibration: shift=" << isotope_calib_map[1].first
                    << " ppm, sigma=" << isotope_calib_map[1].second << " ppm" << endl;
  }

  // Build fast m/z lookup for hills (analogous to hills_mz_median_fast_dict).
  double max_mz_value = 0.0;
  for (const auto& h : hills)
  {
    max_mz_value = max(max_mz_value, h.mz_weighted_mean);
  }
  double mz_step = (max_mz_value > 0.0) ? (htol_ * 1e-6 * max_mz_value) : 0.0;

  struct FastHillEntry
  {
    Size hill_index;
    Size first_scan;
    Size last_scan;
  };

  map<int, vector<FastHillEntry>> hills_mz_fast;
  // Optional ion-mobility bins for hills (analogous to hills_im_median_fast_dict).
  vector<int> hill_im_bins(hills.size(), 0);

  for (Size idx = 0; idx < hills.size(); ++idx)
  {
    const Hill& h = hills[idx];
    if (h.scan_indices.empty()) continue;

    int mz_bin = (mz_step > 0.0) ? static_cast<int>(h.mz_weighted_mean / mz_step) : 0;
    Size first_scan = h.scan_indices.front();
    Size last_scan = h.scan_indices.back();

    // Register each hill in the central bin and its two neighbors
    // (mz_bin-1, mz_bin, mz_bin+1), analogous to the Python
    // hills_mz_median_fast_dict population.
    hills_mz_fast[mz_bin - 1].push_back({idx, first_scan, last_scan});
    hills_mz_fast[mz_bin].push_back({idx, first_scan, last_scan});
    hills_mz_fast[mz_bin + 1].push_back({idx, first_scan, last_scan});

    if (use_im && h.ion_mobility_median >= 0.0)
    {
      hill_im_bins[idx] = static_cast<int>(h.ion_mobility_median / paseftol_);
    }
  }

  // Helper to check if two scan index lists share at least one scan.
  auto hasScanOverlap = [](const vector<Size>& a, const vector<Size>& b) -> bool
  {
    Size i = 0, j = 0;
    while (i < a.size() && j < b.size())
    {
      if (a[i] == b[j]) return true;
      if (a[i] < b[j]) ++i;
      else ++j;
    }
    return false;
  };

  // Local pattern candidate structure mirroring the Python dictionaries.
  struct PatternCandidate
  {
    Size mono_index;                 // index into hills
    double mono_mz = 0.0;
    int charge = 0;
    double cos_cor_isotopes = 0.0;   // cosine in isotope-intensity space
    vector<IsotopeCandidate> isotopes;
    Size n_scans = 0;
  };

  vector<PatternCandidate> ready;

  // Charge ordering: high to low.
  vector<int> charges;
  for (int c = min_charge; c <= max_charge; ++c) charges.push_back(c);
  reverse(charges.begin(), charges.end());

  // Charge ban map as in the Python reference.
  map<int, vector<int>> charge_ban_map = {
    {8, {1, 2, 4}},
    {7, {1}},
    {6, {1, 2, 3}},
    {5, {1}},
    {4, {1, 2}},
    {3, {1}},
    {2, {1}},
    {1, {1}}
  };

  // Fast lookup from hill_idx to index into the hills vector
  // to avoid repeated linear scans when accessing apex intensities.
  map<Size, Size> hill_idx_to_index;
  for (Size idx = 0; idx < hills.size(); ++idx)
  {
    hill_idx_to_index[hills[idx].hill_idx] = idx;
  }

  // Initial isotope candidate generation (analogous to get_initial_isotopes).
  for (Size i = 0; i < hills.size(); ++i)
  {
    const Hill& mono_hill = hills[i];
    if (mono_hill.scan_indices.empty()) continue;

    double mono_mz = mono_hill.mz_weighted_mean;
    double mz_tol = itol_ppm * 1e-6 * mono_mz;

    Size hill_scans_1_number = mono_hill.length;
    const vector<Size>& scans1 = mono_hill.scan_indices;
    Size hill_scans_1_list_first = scans1.front();
    Size hill_scans_1_list_last = scans1.back();

    map<int, Size> banned_charges;

    for (int charge : charges)
    {
      vector<vector<IsotopeCandidate>> candidates;

      for (int iso_num = 1; iso_num <= 9; ++iso_num)
      {
        vector<IsotopeCandidate> tmp_candidates;
        double expected_mz = mono_mz + ISOTOPE_MASSDIFF * iso_num / static_cast<double>(charge);
        int m_to_check_fast = (mz_step > 0.0) ? static_cast<int>(expected_mz / mz_step) : 0;

	        auto it_bin = hills_mz_fast.find(m_to_check_fast);
	        if (it_bin != hills_mz_fast.end())
	        {
	          for (const auto& entry : it_bin->second)
	          {
	            Size idx2 = entry.hill_index;
	            if (idx2 == i) continue;
	            const Hill& hill2 = hills[idx2];
	            if (hill2.scan_indices.empty()) continue;

	            Size hill_scans_2_list_first = entry.first_scan;
	            Size hill_scans_2_list_last = entry.last_scan;

	            if (hill_scans_1_list_last < hill_scans_2_list_first ||
	                hill_scans_2_list_last < hill_scans_1_list_first)
	            {
	              continue;
	            }

	            // Optional ion mobility gating for isotope hills.
	            if (use_im && hill_im_bins[i] != 0 && hill_im_bins[idx2] != 0)
	            {
	              if (abs(hill_im_bins[i] - hill_im_bins[idx2]) > 1)
	              {
	                continue;
	              }
	            }

	            double mass_diff_abs = hill2.mz_weighted_mean - expected_mz;
	            if (fabs(mass_diff_abs) > mz_tol)
	            {
	              continue;
	            }

	            if (!hasScanOverlap(scans1, hill2.scan_indices))
	            {
	              continue;
	            }

	            double cos_cor_RT = cosineCorrelation_(mono_hill.intensities, mono_hill.scan_indices,
	                                                   hill2.intensities, hill2.scan_indices);
	            if (cos_cor_RT < 0.6)
	            {
	              continue;
	            }

	            IsotopeCandidate cand;
	            cand.hill_idx = hill2.hill_idx;
	            cand.isotope_number = iso_num;
	            cand.mass_diff_ppm = mass_diff_abs * 1e6 / expected_mz;
	            cand.cos_corr = cos_cor_RT;
	            tmp_candidates.push_back(cand);
	          }
	        }

        if (!tmp_candidates.empty())
        {
          candidates.push_back(tmp_candidates);
        }

        if (candidates.size() < static_cast<Size>(iso_num))
        {
          break;
        }
      }

      Size min_required = 1;
      auto it_ban = banned_charges.find(charge);
      if (it_ban != banned_charges.end())
      {
        min_required = it_ban->second;
      }

        if (candidates.size() >= min_required)
        {
          double hill_intensity_apex_1 = mono_hill.intensity_apex;
          auto averagine = computeAveragine_(mono_mz * charge, hill_intensity_apex_1);
          vector<double> all_theoretical_int = averagine.first;
          Size max_pos = averagine.second;

          // Enumerate all combinations of isotope candidates across orders.
          if (!candidates.empty())
          {
            vector<Size> indices(candidates.size(), 0);
            bool done = false;

            while (!done)
            {
              vector<IsotopeCandidate> picked_isotopes;
              picked_isotopes.reserve(candidates.size());
              for (Size k = 0; k < candidates.size(); ++k)
              {
                picked_isotopes.push_back(candidates[k][indices[k]]);
              }

              vector<double> all_exp_intensity;
              all_exp_intensity.reserve(picked_isotopes.size() + 1);
              all_exp_intensity.push_back(hill_intensity_apex_1);

              double local_minimum = 0.0;
              Size local_minimum_pos = 0;

              Size i_local_isotope = 1;
              for (const auto& iso_cand : picked_isotopes)
              {
                // Find apex intensity of isotope hill (fast lookup by hill_idx).
                double hill_intensity_apex_2 = 0.0;
                auto it_idx = hill_idx_to_index.find(iso_cand.hill_idx);
                if (it_idx != hill_idx_to_index.end())
                {
                  hill_intensity_apex_2 = hills[it_idx->second].intensity_apex;
                }

                if (i_local_isotope > max_pos)
                {
                  if (i_local_isotope == max_pos + 1 || hill_intensity_apex_2 < local_minimum)
                  {
                    local_minimum = hill_intensity_apex_2;
                    local_minimum_pos = i_local_isotope;
                  }
                  if (hill_intensity_apex_2 >= ivf * local_minimum)
                  {
                    if (local_minimum_pos + 1 < all_exp_intensity.size())
                    {
                      all_exp_intensity.resize(local_minimum_pos + 1);
                    }
                    break;
                  }
                }

                all_exp_intensity.push_back(hill_intensity_apex_2);
                ++i_local_isotope;
              }

              // Compute cosine correlation and optimal truncation in isotope-intensity space.
              auto cc = checkingCosCorrelationForCarbon_(all_theoretical_int, all_exp_intensity, 0.6);
              double cos_corr_iso = cc.first;
              Size best_pos = cc.second;

              if (cos_corr_iso > 0.0 && best_pos > 1)
              {
                Size n_iso_used = best_pos - 1;
                n_iso_used = min(n_iso_used, picked_isotopes.size());

                PatternCandidate pc;
                pc.mono_index = i;
                pc.mono_mz = mono_mz;
                pc.charge = charge;
                pc.cos_cor_isotopes = cos_corr_iso;
                pc.n_scans = hill_scans_1_number;

                pc.isotopes.assign(picked_isotopes.begin(),
                                   picked_isotopes.begin() + n_iso_used);

                if (!pc.isotopes.empty())
                {
                  ready.push_back(pc);

                  for (int ch_v : charge_ban_map[charge])
                  {
                    auto& ref = banned_charges[ch_v];
                    ref = max(ref, n_iso_used);
                  }
                }
              }

              // Next combination.
              Size pos_idx = 0;
              while (pos_idx < indices.size())
              {
                ++indices[pos_idx];
                if (indices[pos_idx] < candidates[pos_idx].size())
                {
                  break;
                }
                indices[pos_idx] = 0;
                ++pos_idx;
              }
              if (pos_idx == indices.size())
              {
                done = true;
              }
            }
          }
        }
    }
  }

  // Isotope mass error calibration based on the initial candidates.
  map<int, vector<double>> isotope_errors_ready;
  for (int ic = 1; ic <= 9; ++ic) isotope_errors_ready[ic] = vector<double>();

  // Optional RT-apex gating before isotope mass calibration:
  // discard isotope hills whose apex RT is farther than hrttol_ seconds
  // from the monoisotopic hill apex.
  if (hrttol_ > 0.0 && !ready.empty())
  {
    vector<PatternCandidate> rt_filtered;
    rt_filtered.reserve(ready.size());

    for (auto pc : ready)
    {
      const Hill& mono_hill = hills[pc.mono_index];
      double mono_rt_apex = mono_hill.rt_apex;

      vector<IsotopeCandidate> kept;
      kept.reserve(pc.isotopes.size());

      for (const auto& iso_cand : pc.isotopes)
      {
        auto it_h = hill_idx_to_index.find(iso_cand.hill_idx);
        if (it_h == hill_idx_to_index.end())
        {
          continue;
        }
        const Hill& iso_hill = hills[it_h->second];
        double iso_rt_apex = iso_hill.rt_apex;
        if (fabs(iso_rt_apex - mono_rt_apex) <= hrttol_)
        {
          kept.push_back(iso_cand);
        }
      }

      if (!kept.empty())
      {
        pc.isotopes.swap(kept);
        rt_filtered.push_back(pc);
      }
    }

    OPENMS_LOG_INFO << "After RT apex filter (hrttol=" << hrttol_
                    << " s), " << rt_filtered.size()
                    << " potential isotope patterns remain." << endl;

    ready.swap(rt_filtered);
  }

  for (const auto& pc : ready)
  {
    if (pc.n_scans < 3) continue;
    for (Size idx = 0; idx < pc.isotopes.size(); ++idx)
    {
      int iso_num = static_cast<int>(pc.isotopes[idx].isotope_number);
      if (iso_num >= 1 && iso_num <= 9)
      {
        isotope_errors_ready[iso_num].push_back(pc.isotopes[idx].mass_diff_ppm);
      }
    }
  }

  map<int, pair<double, double>> isotope_calib_map_ready;
  for (int ic = 1; ic <= 9; ++ic)
  {
    isotope_calib_map_ready[ic] = make_pair(0.0, itol_ppm);
  }

  if (enable_isotope_calib)
  {
    for (int ic = 1; ic <= 3; ++ic)
    {
      if (isotope_errors_ready[ic].size() >= 1000)
      {
        auto calib = calibrateMass_(isotope_errors_ready[ic]);
        isotope_calib_map_ready[ic] = calib;
      }
    }

    for (int ic = 4; ic <= 9; ++ic)
    {
      if (isotope_errors_ready[ic].size() >= 1000)
      {
        isotope_calib_map_ready[ic] = calibrateMass_(isotope_errors_ready[ic]);
      }
      else if (ic > 1)
      {
        auto prev = isotope_calib_map_ready[ic - 1];
        auto prev2 = (ic > 2) ? isotope_calib_map_ready[ic - 2] : make_pair(0.0, itol_ppm);
        double shift_delta = prev.first - prev2.first;
        double sigma_ratio = prev.second / max(prev2.second, 0.1);
        isotope_calib_map_ready[ic] = make_pair(prev.first + shift_delta, prev.second * sigma_ratio);
      }
    }
  }

  // Apply calibrated isotope mass filters and re-check isotope-intensity cosine.
  vector<PatternCandidate> filtered_ready;
  filtered_ready.reserve(ready.size());
  for (auto pc : ready)
  {
    if (enable_isotope_calib)
    {
      vector<IsotopeCandidate> tmp;
      for (const auto& cand : pc.isotopes)
      {
        int iso_num = static_cast<int>(cand.isotope_number);
        auto calib = isotope_calib_map_ready[iso_num];
        if (fabs(cand.mass_diff_ppm - calib.first) <= 5.0 * calib.second)
        {
          tmp.push_back(cand);
        }
        else
        {
          break;
        }
      }
      pc.isotopes.swap(tmp);
    }

    if (pc.isotopes.empty())
    {
      continue;
    }

    // Recompute cosine correlation in isotope-intensity space and
    // potentially truncate the series, analogous to the second pass
    // in process_features_iteration.
    const Hill& mono_hill = hills[pc.mono_index];
    double mono_mz = mono_hill.mz_weighted_mean;
    double hill_intensity_apex_1 = mono_hill.intensity_apex;

    auto averagine = computeAveragine_(mono_mz * pc.charge, hill_intensity_apex_1);
    vector<double> all_theoretical_int = averagine.first;

    vector<double> all_exp_intensity;
    all_exp_intensity.reserve(pc.isotopes.size() + 1);
    all_exp_intensity.push_back(hill_intensity_apex_1);

    for (const auto& iso_cand : pc.isotopes)
    {
      double hill_intensity_apex_2 = 0.0;
      auto it_idx = hill_idx_to_index.find(iso_cand.hill_idx);
      if (it_idx != hill_idx_to_index.end())
      {
        hill_intensity_apex_2 = hills[it_idx->second].intensity_apex;
      }
      all_exp_intensity.push_back(hill_intensity_apex_2);
    }

    auto cc = checkingCosCorrelationForCarbon_(all_theoretical_int, all_exp_intensity, 0.6);
    double cos_corr_iso = cc.first;
    Size best_pos = cc.second;

    if (cos_corr_iso <= 0.0 || best_pos <= 1)
    {
      continue;
    }

    Size n_iso_used = best_pos - 1;
    n_iso_used = min(n_iso_used, pc.isotopes.size());
    pc.isotopes.resize(n_iso_used);
    pc.cos_cor_isotopes = cos_corr_iso;

    if (!pc.isotopes.empty())
    {
      filtered_ready.push_back(pc);
    }
  }

  // Final greedy selection of non-overlapping patterns (analogous to ready_final).
  sort(filtered_ready.begin(), filtered_ready.end(),
       [](const PatternCandidate& a, const PatternCandidate& b)
       {
         if (a.isotopes.size() != b.isotopes.size())
         {
           return a.isotopes.size() > b.isotopes.size();
         }
         return a.cos_cor_isotopes > b.cos_cor_isotopes;
       });

  // Build a lookup from hill_idx to Hill pointer for fast access to
  // apex intensities during the truncation / re-correlation step.
  map<Size, const Hill*> hill_lookup_for_iso;
  for (const auto& h : hills)
  {
    hill_lookup_for_iso[h.hill_idx] = &h;
  }

  set<Size> occupied_hills;
  for (const auto& pc_in : filtered_ready)
  {
    PatternCandidate pc = pc_in; // local copy that we may truncate
    const Hill& mono_hill = hills[pc.mono_index];
    const double mono_mz_center = mono_hill.mz_weighted_mean;

    // Skip patterns whose monoisotopic hill is already used.
    if (occupied_hills.find(mono_hill.hill_idx) != occupied_hills.end())
    {
      continue;
    }

    bool iso_conflict = false;
    for (const auto& iso : pc.isotopes)
    {
      if (occupied_hills.find(iso.hill_idx) != occupied_hills.end())
      {
        iso_conflict = true;
        break;
      }
    }

    if (iso_conflict)
    {
      // Try to keep only leading isotopes whose hills are not yet used,
      // mirroring the Python ready_final truncation behavior.
      vector<IsotopeCandidate> tmp_iso;
      for (const auto& iso : pc.isotopes)
      {
        if (occupied_hills.find(iso.hill_idx) == occupied_hills.end())
        {
          tmp_iso.push_back(iso);
        }
        else
        {
          break;
        }
      }
      if (tmp_iso.empty()) continue;

      // Recompute cosine correlation in isotope-intensity space for the
      // truncated pattern and potentially truncate further based on the
      // averagine-explained / correlation criterion, analogous to
      // checking_cos_correlation_for_carbon in the Python code.
      double mono_mz = mono_hill.mz_weighted_mean;
      double hill_intensity_apex_1 = mono_hill.intensity_apex;

      auto averagine = computeAveragine_(mono_mz * pc.charge, hill_intensity_apex_1);
      const vector<double>& theor_full = averagine.first;

      Size len = min(static_cast<Size>(tmp_iso.size()) + 1, theor_full.size());
      if (len <= 1) continue;

      vector<double> theor(len);
      vector<double> exp(len);

      for (Size k = 0; k < len; ++k)
      {
        theor[k] = theor_full[k];
      }

      exp[0] = hill_intensity_apex_1;
      for (Size k = 1; k < len; ++k)
      {
        Size iso_index = k - 1;
        auto it_h = hill_lookup_for_iso.find(tmp_iso[iso_index].hill_idx);
        if (it_h != hill_lookup_for_iso.end())
        {
          exp[k] = it_h->second->intensity_apex;
        }
        else
        {
          exp[k] = 0.0;
        }
      }

      auto cc = checkingCosCorrelationForCarbon_(theor, exp, 0.6);
      double cos_corr_iso = cc.first;
      Size best_pos = cc.second;

      if (cos_corr_iso <= 0.0 || best_pos <= 1)
      {
        continue;
      }

      Size n_iso_used = min(best_pos - 1, tmp_iso.size());
      tmp_iso.resize(n_iso_used);
      pc.isotopes = tmp_iso;
      pc.cos_cor_isotopes = cos_corr_iso;
    }

    if (pc.isotopes.empty())
    {
      continue;
    }

    // Debug sanity checks for the final isotope pattern (optional).
    for (const auto& iso_cand : pc.isotopes)
    {
      auto it_h = hill_lookup_for_iso.find(iso_cand.hill_idx);
      if (it_h == hill_lookup_for_iso.end())
      {
        continue;
      }
      const Hill* iso_hill_ptr = it_h->second;
      debugCheckIsotopeConsistency_("detectIsotopePatterns_",
                                    mono_mz_center,
                                    mono_hill.rt_apex,
                                    mono_hill.hill_idx,
                                    pc.charge,
                                    itol_ppm,
                                    *iso_hill_ptr,
                                    iso_cand.isotope_number);
    }

    // At this point, either there was no conflict or we have a
    // truncated pattern that still passes the cosine check.
    PeptideFeature feature;
    feature.mz = pc.mono_mz;
    feature.rt_start = mono_hill.rt_start;
    feature.rt_end = mono_hill.rt_end;
    feature.rt_apex = mono_hill.rt_apex;
    feature.intensity_apex = mono_hill.intensity_apex;
    feature.intensity_sum = mono_hill.intensity_sum;

    if (iuse != 0)
    {
      int isotopes_to_add = (iuse == -1) ? static_cast<int>(pc.isotopes.size())
                                         : min(static_cast<int>(pc.isotopes.size()), iuse);
      for (int iso_idx = 0; iso_idx < isotopes_to_add; ++iso_idx)
      {
        auto it_h = hill_lookup_for_iso.find(pc.isotopes[static_cast<Size>(iso_idx)].hill_idx);
        if (it_h != hill_lookup_for_iso.end())
        {
          const Hill* h = it_h->second;
          feature.intensity_apex += h->intensity_apex;
          feature.intensity_sum += h->intensity_sum;
        }
      }
    }

    feature.charge = pc.charge;
    feature.n_isotopes = pc.isotopes.size() + 1;
    feature.n_scans = mono_hill.length;
    feature.isotopes = pc.isotopes;
    feature.mono_hill_idx = mono_hill.hill_idx;
    feature.drift_time = mono_hill.drift_time_median;
    feature.ion_mobility = mono_hill.ion_mobility_median;

    double proton_mass = Constants::PROTON_MASS_U;
    if (negative_mode)
    {
      feature.mass_calib = pc.mono_mz * pc.charge + proton_mass * pc.charge;
    }
    else
    {
      feature.mass_calib = pc.mono_mz * pc.charge - proton_mass * pc.charge;
    }

    features.push_back(feature);
    occupied_hills.insert(mono_hill.hill_idx);
    for (const auto& iso : pc.isotopes)
    {
      occupied_hills.insert(iso.hill_idx);
    }
  }

  OPENMS_LOG_INFO << "Detected " << features.size() << " features with isotope patterns" << endl;
  return features;
}

FeatureMap Biosaur2Algorithm::convertToFeatureMap_(const vector<PeptideFeature>& features,
                                                   const vector<Hill>& hills) const
{
  FeatureMap feature_map;

   const bool use_mass_trace_hulls = (convex_hull_mode_ == "mass_traces");

  // Build a lookup from hill index to Hill pointer so we can reconstruct
  // per-isotope convex hulls for each peptide feature.
  map<Size, const Hill*> hill_lookup;
  for (const auto& h : hills)
  {
    hill_lookup[h.hill_idx] = &h;
  }

  for (const auto& f : features)
  {
    const Hill* mono_hill_ptr = nullptr;
    auto mono_it = hill_lookup.find(f.mono_hill_idx);
    if (mono_it != hill_lookup.end())
    {
      mono_hill_ptr = mono_it->second;
    }
    if (mono_hill_ptr != nullptr && !f.isotopes.empty())
    {
      for (const auto& iso : f.isotopes)
      {
        auto iso_it = hill_lookup.find(iso.hill_idx);
        if (iso_it == hill_lookup.end())
        {
          continue;
        }
        const Hill* iso_hill_ptr = iso_it->second;
        debugCheckIsotopeConsistency_("convertToFeatureMap_",
                                      mono_hill_ptr->mz_weighted_mean,
                                      mono_hill_ptr->rt_apex,
                                      mono_hill_ptr->hill_idx,
                                      f.charge,
                                      itol_,
                                      *iso_hill_ptr,
                                      iso.isotope_number);
      }
    }

    Feature feature;
    feature.setMZ(f.mz);
    feature.setRT(f.rt_apex);
    feature.setIntensity(f.intensity_apex);
    feature.setCharge(f.charge);
    feature.setOverallQuality(f.n_isotopes);

    // Collect the monoisotopic hill and all isotope hills contributing to
    // this feature and build one convex hull per hill, analogous to how
    // FeatureFinderCentroided uses mass-trace hulls.
    set<Size> pattern_hill_ids;
    pattern_hill_ids.insert(f.mono_hill_idx);
    for (const auto& iso : f.isotopes)
    {
      pattern_hill_ids.insert(iso.hill_idx);
    }

    DBoundingBox<2> feature_box;
    bool have_box_points = false;

    for (Size hill_id : pattern_hill_ids)
    {
      auto it = hill_lookup.find(hill_id);
      if (it == hill_lookup.end())
      {
        continue;
      }
      const Hill& hill = *(it->second);
      if (hill.rt_values.empty() || hill.mz_values.empty())
      {
        continue;
      }

      const Size n_pts = min(hill.rt_values.size(), hill.mz_values.size());
      if (n_pts == 0)
      {
        continue;
      }

      if (use_mass_trace_hulls)
      {
        ConvexHull2D::PointArrayType hull_points(n_pts);
        for (Size i = 0; i < n_pts; ++i)
        {
          hull_points[i][0] = hill.rt_values[i];
          hull_points[i][1] = hill.mz_values[i];
        }

        ConvexHull2D hull;
        hull.addPoints(hull_points);
        feature.getConvexHulls().push_back(hull);
      }
      else
      {
        for (Size i = 0; i < n_pts; ++i)
        {
          const double rt = hill.rt_values[i];
          const double mz = hill.mz_values[i];
          if (!have_box_points)
          {
            DBoundingBox<2>::PositionType p(rt, mz);
            feature_box = DBoundingBox<2>(p, p);
            have_box_points = true;
          }
          else
          {
            feature_box.enlarge(rt, mz);
          }
        }
      }
    }

    if (!use_mass_trace_hulls && have_box_points)
    {
      ConvexHull2D::PointArrayType hull_points(4);
      hull_points[0][0] = feature_box.minX();
      hull_points[0][1] = feature_box.minY();
      hull_points[1][0] = feature_box.maxX();
      hull_points[1][1] = feature_box.minY();
      hull_points[2][0] = feature_box.minX();
      hull_points[2][1] = feature_box.maxY();
      hull_points[3][0] = feature_box.maxX();
      hull_points[3][1] = feature_box.maxY();

      ConvexHull2D hull;
      hull.addPoints(hull_points);
      feature.getConvexHulls().push_back(hull);
    }

    // Fallback: if something went wrong while resolving hills, keep the
    // previous simple two-point hull to avoid features without any hull.
    if (feature.getConvexHulls().empty())
    {
      ConvexHull2D::PointArrayType hull_points(2);
      hull_points[0][0] = f.rt_start;
      hull_points[0][1] = f.mz;
      hull_points[1][0] = f.rt_end;
      hull_points[1][1] = f.mz;

      ConvexHull2D hull;
      hull.addPoints(hull_points);
      feature.getConvexHulls().push_back(hull);
    }

    feature.setMetaValue("mass_calib", f.mass_calib);
    feature.setMetaValue("n_isotopes", f.n_isotopes);
    feature.setMetaValue("n_scans", f.n_scans);
    feature.setMetaValue("intensity_sum", f.intensity_sum);
    if (f.drift_time >= 0)
    {
      feature.setMetaValue("FAIMS_compensation_voltage", f.drift_time);
    }
    if (f.ion_mobility >= 0)
    {
    feature.setMetaValue("ion_mobility", f.ion_mobility);
    }

    feature.ensureUniqueId();
    feature_map.push_back(feature);
  }

  feature_map.applyMemberFunction(&UniqueIdInterface::ensureUniqueId);
  feature_map.ensureUniqueId();
  feature_map.getProteinIdentifications().resize(1);
  return feature_map;
}

void Biosaur2Algorithm::debugCheckIsotopeConsistency_(const char* stage_label,
                                                      double mono_mz_center,
                                                      double mono_rt_apex,
                                                      Size mono_hill_idx,
                                                      int charge,
                                                      double itol_ppm,
                                                      const Hill& iso_hill,
                                                      Size isotope_number) const
{
  // RT-apex sanity check.
  if (hrttol_ > 0.0)
  {
    const double rt_delta = fabs(iso_hill.rt_apex - mono_rt_apex);
    if (rt_delta > hrttol_ + 1e-6)
    {
      OPENMS_LOG_WARN << "Biosaur2 isotope debug (" << stage_label << "): "
                      << "mono m/z " << mono_mz_center
                      << " (charge " << charge
                      << ", mono hill_idx=" << mono_hill_idx
                      << ") uses isotope hill_idx=" << iso_hill.hill_idx
                      << " with RT apex delta " << rt_delta
                      << " s > hrttol=" << hrttol_
                      << ". This indicates an inconsistent isotope assignment."
                      << endl;
    }
  }

  // m/z sanity check.
  if (charge > 0)
  {
    const double ISOTOPE_MASSDIFF = Constants::C13C12_MASSDIFF_U;
    const double expected_mz = mono_mz_center +
                               static_cast<double>(isotope_number) *
                                 ISOTOPE_MASSDIFF /
                                 static_cast<double>(charge);
    const double observed_mz = iso_hill.mz_weighted_mean;
    const double diff_ppm = Math::getPPM(observed_mz, expected_mz);
    const double mz_ppm_threshold = std::max(80.0, 10.0 * itol_ppm);

    if (fabs(diff_ppm) > mz_ppm_threshold)
    {
      OPENMS_LOG_WARN << "Biosaur2 isotope debug (" << stage_label << "): "
                      << "mono m/z " << mono_mz_center
                      << " (charge " << charge
                      << ", mono hill_idx=" << mono_hill_idx
                      << ") uses isotope #" << isotope_number
                      << " (hill_idx=" << iso_hill.hill_idx
                      << ") at m/z=" << observed_mz
                      << " which is " << diff_ppm
                      << " ppm away from expected " << expected_mz
                      << " (itol=" << itol_ppm << " ppm)."
                      << endl;
    }
  }
}

double Biosaur2Algorithm::cosineCorrelation_(const vector<double>& intensities1,
                                             const vector<Size>& scans1,
                                             const vector<double>& intensities2,
                                             const vector<Size>& scans2) const
{
  map<Size, double> map1, map2;
  for (Size i = 0; i < scans1.size(); ++i)
  {
    map1[scans1[i]] = intensities1[i];
  }
  for (Size i = 0; i < scans2.size(); ++i)
  {
    map2[scans2[i]] = intensities2[i];
  }

  double dot_product = 0.0;
  double norm1 = 0.0;
  double norm2 = 0.0;

  for (const auto& p1 : map1)
  {
    Size scan = p1.first;
    double i1 = p1.second;
    auto it = map2.find(scan);
    if (it != map2.end())
    {
      dot_product += i1 * it->second;
    }
    norm1 += i1 * i1;
  }

  for (const auto& p2 : map2)
  {
    norm2 += p2.second * p2.second;
  }

  if (norm1 == 0.0 || norm2 == 0.0)
  {
    return 0.0;
  }
  return dot_product / (sqrt(norm1) * sqrt(norm2));
}

void Biosaur2Algorithm::writeTSV(const vector<PeptideFeature>& features, const String& filename) const
{
  ofstream out(filename);
  if (!out)
  {
    throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
  }
  out << "massCalib\trtApex\tintensityApex\tintensitySum\tcharge\t"
      << "nIsotopes\tnScans\tmz\trtStart\trtEnd\tFAIMS\tIM" << endl;

  for (const auto& f : features)
  {
    out << f.mass_calib << "\t"
        << f.rt_apex << "\t"
        << f.intensity_apex << "\t"
        << f.intensity_sum << "\t"
        << f.charge << "\t"
        << f.n_isotopes << "\t"
        << f.n_scans << "\t"
        << f.mz << "\t"
        << f.rt_start << "\t"
        << f.rt_end << "\t"
        << f.drift_time << "\t"
        << f.ion_mobility << endl;
  }

  OPENMS_LOG_INFO << "Wrote " << features.size() << " features to TSV file: " << filename << endl;
}

void Biosaur2Algorithm::writeHills(const vector<Hill>& hills, const String& filename) const
{
  ofstream out(filename);
  if (!out)
  {
    throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
  }
  out << "hill_idx\tmz\trtStart\trtEnd\trtApex\tintensityApex\tintensitySum\tnScans" << endl;

  for (const auto& hill : hills)
  {
    out << hill.hill_idx << "\t"
        << hill.mz_weighted_mean << "\t"
        << hill.rt_start << "\t"
        << hill.rt_end << "\t"
        << hill.rt_apex << "\t"
        << hill.intensity_apex << "\t"
        << hill.intensity_sum << "\t"
        << hill.length << endl;
  }

  OPENMS_LOG_INFO << "Wrote " << hills.size() << " hills to: " << filename << endl;
}

} // namespace OpenMS
