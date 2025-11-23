// Copyright (c) 2002-present, OpenMS Inc.
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/FEATUREFINDER/Biosaur2Algorithm.h>

#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/IONMOBILITY/FAIMSHelper.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <map>
#include <numeric>
#include <set>

#ifdef _OPENMP
#  include <omp.h>
#endif

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

  defaults_.setValue("cmin", 1, "Minimum charge state");
  defaults_.setMinInt("cmin", 1);

  defaults_.setValue("cmax", 6, "Maximum charge state");
  defaults_.setMinInt("cmax", 1);

  defaults_.setValue("threads", 1, "Number of threads for parallel processing (0 = auto-detect)");
  defaults_.setMinInt("threads", 0);

  defaults_.setValue("iuse", 0, "Number of isotopes for intensity calculation (0=mono only, -1=all, 1=mono+first, etc.)");
  defaults_.setMinInt("iuse", -1);

  defaults_.setValue("nm", "false", "Negative mode (affects neutral mass calculation)");
  defaults_.setValidStrings("nm", {"true", "false"});

  defaults_.setValue("tof", "false", "Enable TOF-specific intensity filtering");
  defaults_.setValidStrings("tof", {"true", "false"});

  defaults_.setValue("profile", "false", "Enable profile mode processing (centroid spectra using PeakPickerHiRes)");
  defaults_.setValidStrings("profile", {"true", "false"});

  defaults_.setValue("use_hill_calib", "false", "Enable automatic hill mass tolerance calibration");
  defaults_.setValidStrings("use_hill_calib", {"true", "false"});

  defaults_.setValue("ignore_iso_calib", "false", "Disable automatic isotope mass error calibration");
  defaults_.setValidStrings("ignore_iso_calib", {"true", "false"});

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
  threads_ = param_.getValue("threads");
  iuse_ = param_.getValue("iuse");
  negative_mode_ = param_.getValue("nm").toBool();
  tof_mode_ = param_.getValue("tof").toBool();
  profile_mode_ = param_.getValue("profile").toBool();
  use_hill_calib_ = param_.getValue("use_hill_calib").toBool();
  ignore_iso_calib_ = param_.getValue("ignore_iso_calib").toBool();
}

void Biosaur2Algorithm::run(const MSExperiment& input, FeatureMap& feature_map)
{
  vector<Hill> tmp_hills;
  vector<PeptideFeature> tmp_features;
  run(input, feature_map, tmp_hills, tmp_features);
}

void Biosaur2Algorithm::run(const MSExperiment& input,
                            FeatureMap& feature_map,
                            vector<Hill>& hills,
                            vector<PeptideFeature>& peptide_features)
{
  MSExperiment exp = input;
  feature_map.clear(true);
  hills.clear();
  peptide_features.clear();

#ifdef _OPENMP
  int threads_to_use = threads_;
  if (threads_to_use == 0)
  {
    threads_to_use = omp_get_max_threads();
  }
  omp_set_num_threads(threads_to_use);
  OPENMS_LOG_INFO << "Using " << threads_to_use << " threads for parallel processing" << endl;
#else
  (void)threads_;
#endif

  exp.getSpectra().erase(
    remove_if(exp.begin(), exp.end(),
              [](const MSSpectrum& s) { return s.getMSLevel() != 1; }),
    exp.end());

  OPENMS_LOG_INFO << "Loaded " << exp.size() << " MS1 spectra" << endl;

  if (exp.empty())
  {
    throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "No MS1 spectra found in input experiment!", "");
  }

  auto cvs = FAIMSHelper::getCompensationVoltages(exp);
  if (!cvs.empty())
  {
    OPENMS_LOG_INFO << "Detected FAIMS data with " << cvs.size() << " compensation voltage(s)" << endl;
    for (const auto& cv : cvs)
    {
      OPENMS_LOG_INFO << "  CV: " << cv << " V" << endl;
    }
  }

  bool has_ion_mobility = false;
  if (!exp.empty())
  {
    const auto& fda = exp[0].getFloatDataArrays();
    auto it_im = getDataArrayByName(fda, Constants::UserParam::ION_MOBILITY);
    if (it_im != fda.end())
    {
      has_ion_mobility = true;
      OPENMS_LOG_INFO << "Detected ion mobility (PASEF/IM) data" << endl;
    }
  }
  (void)has_ion_mobility;

  if (profile_mode_)
  {
    centroidProfileSpectra(exp);
  }

  if (tof_mode_)
  {
    processTOF(exp);
  }

  double calibrated_htol = htol_;
  if (use_hill_calib_)
  {
    OPENMS_LOG_INFO << "Performing hill mass tolerance calibration..." << endl;
    vector<double> mass_diffs;
    Size sample_size = min(exp.size(), Size(1000));
    Size start_idx = (exp.size() > 1000) ? (exp.size() / 2 - 500) : 0;

    MSExperiment calib_exp;
    for (Size i = start_idx; i < start_idx + sample_size && i < exp.size(); ++i)
    {
      calib_exp.addSpectrum(exp[i]);
    }

    vector<Hill> calib_hills = detectHills(calib_exp, htol_, mini_, minmz_, maxmz_, &mass_diffs);
    (void)calib_hills;
    if (!mass_diffs.empty())
    {
      auto calib = calibrateMass(mass_diffs);
      double calibrated_sigma = calib.second;
      calibrated_htol = min(htol_, 5.0 * calibrated_sigma);
      OPENMS_LOG_INFO << "Automatically optimized htol parameter: "
                      << calibrated_htol << " ppm (was " << htol_ << " ppm)" << endl;
    }
  }

  hills = detectHills(exp, calibrated_htol, mini_, minmz_, maxmz_);
  hills = processHills(hills, minlh_);
  hills = splitHills(hills, hvf_, minlh_);

  bool enable_isotope_calib = !ignore_iso_calib_;
  peptide_features = detectIsotopePatterns(hills, itol_, cmin_, cmax_, negative_mode_, ivf_, iuse_, enable_isotope_calib);

  feature_map = convertToFeatureMap(peptide_features);
}

double Biosaur2Algorithm::calculatePPM(double mz1, double mz2) const
{
  if (fabs(mz2) < 1e-10)
  {
    return 0.0;
  }
  return (mz1 - mz2) / mz2 * 1e6;
}

double Biosaur2Algorithm::calculateMedian(const vector<double>& values) const
{
  if (values.empty())
  {
    return 0.0;
  }

  vector<double> sorted = values;
  return Math::median(sorted.begin(), sorted.end(), false);
}

bool Biosaur2Algorithm::shouldThrowForMissingIM(const MSSpectrum& spectrum) const
{
  // Only throw if this is true ion mobility (TIMS/PASEF), not FAIMS
  // FAIMS uses spectrum-level CV, not per-peak ion mobility arrays
  DriftTimeUnit dt_unit = spectrum.getDriftTimeUnit();
  return (dt_unit != DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE && dt_unit != DriftTimeUnit::NONE);
}

vector<double> Biosaur2Algorithm::meanFilter(const vector<double>& data, Size window) const
{
  vector<double> result(data.size());
  Size half_window = window / 2;

  for (Size i = 0; i < data.size(); ++i)
  {
    Size start = (i >= half_window) ? i - half_window : 0;
    Size end = min(i + half_window + 1, data.size());

    double sum = 0.0;
    for (Size j = start; j < end; ++j)
    {
      sum += data[j];
    }
    result[i] = sum / static_cast<double>(end - start);
  }

  return result;
}

pair<double, double> Biosaur2Algorithm::calibrateMass(const vector<double>& mass_errors, double bin_width) const
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
    return calibrateMass(mass_errors, 0.25);
  }

  if (isinf(sigma) || isnan(sigma))
  {
    return make_pair(0.0, 10.0);
  }

  return make_pair(mean, sigma);
}

void Biosaur2Algorithm::processTOF(MSExperiment& exp) const
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

void Biosaur2Algorithm::centroidProfileSpectra(MSExperiment& exp) const
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

vector<Biosaur2Algorithm::Hill> Biosaur2Algorithm::detectHills(const MSExperiment& exp,
                                                               double htol_ppm,
                                                               double min_intensity,
                                                               double min_mz,
                                                               double max_mz,
                                                               vector<double>* hill_mass_diffs) const
{
  vector<Hill> hills;
  Size hill_idx_counter = 0;

  for (Size scan_idx = 0; scan_idx < exp.size(); ++scan_idx)
  {
    const MSSpectrum& spectrum = exp[scan_idx];
    double rt = spectrum.getRT();

    for (Size peak_idx = 0; peak_idx < spectrum.size(); ++peak_idx)
    {
      const Peak1D& peak = spectrum[peak_idx];
      double mz = peak.getMZ();
      double intensity = peak.getIntensity();

      if (intensity < min_intensity || mz < min_mz || mz > max_mz)
      {
        continue;
      }

      bool added_to_existing = false;
      for (auto& hill : hills)
      {
        if (!hill.scan_indices.empty() && hill.scan_indices.back() == scan_idx - 1)
        {
          double ppm = fabs(calculatePPM(mz, hill.mz_median));
          if (ppm <= htol_ppm)
          {
            hill.scan_indices.push_back(scan_idx);
            hill.peak_indices.push_back(peak_idx);
            hill.mz_values.push_back(mz);
            hill.intensities.push_back(intensity);
            hill.rt_values.push_back(rt);

            double drift_time = spectrum.getDriftTime();
            hill.drift_times.push_back(drift_time);

            double ion_mobility = -1.0;
            const auto& fda = spectrum.getFloatDataArrays();
            auto it_im = getDataArrayByName(fda, Constants::UserParam::ION_MOBILITY);
            if (it_im != fda.end())
            {
              if (peak_idx >= it_im->size())
              {
                throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                              "Ion mobility array shorter than peak list.",
                                              String(peak_idx));
              }
              ion_mobility = (*it_im)[peak_idx];
            }
            else if (drift_time >= 0)
            {
              if (shouldThrowForMissingIM(spectrum))
              {
                throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                    "Ion mobility array missing although drift times are present.");
              }
            }
            hill.ion_mobilities.push_back(ion_mobility);

            hill.mz_median = calculateMedian(hill.mz_values);
            hill.rt_end = rt;
            hill.length = hill.scan_indices.size();
            hill.intensity_sum += intensity;
            if (intensity > hill.intensity_apex)
            {
              hill.intensity_apex = intensity;
              hill.rt_apex = rt;
            }

            added_to_existing = true;
            break;
          }
        }
      }

      if (!added_to_existing)
      {
        Hill hill;
        hill.scan_indices.push_back(scan_idx);
        hill.peak_indices.push_back(peak_idx);
        hill.mz_values.push_back(mz);
        hill.intensities.push_back(intensity);
        hill.rt_values.push_back(rt);
        hill.mz_median = mz;
        hill.rt_start = rt;
        hill.rt_end = rt;
        hill.rt_apex = rt;
        hill.intensity_apex = intensity;
        hill.intensity_sum = intensity;
        hill.length = 1;
        hill.hill_idx = hill_idx_counter++;

        double drift_time = spectrum.getDriftTime();
        hill.drift_times.push_back(drift_time);

        double ion_mobility = -1.0;
        const auto& fda = spectrum.getFloatDataArrays();
        auto it_im = getDataArrayByName(fda, Constants::UserParam::ION_MOBILITY);
        if (it_im != fda.end())
        {
          if (peak_idx >= it_im->size())
          {
            throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Ion mobility array shorter than peak list.",
                                          String(peak_idx));
          }
          ion_mobility = (*it_im)[peak_idx];
        }
        else if (drift_time >= 0)
        {
          if (shouldThrowForMissingIM(spectrum))
          {
            throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                                "Ion mobility array missing although drift times are present.");
          }
        }
        hill.ion_mobilities.push_back(ion_mobility);

        hills.push_back(hill);
      }
    }
  }

  if (hill_mass_diffs != nullptr)
  {
    for (const auto& hill : hills)
    {
      if (hill.length >= 2)
      {
        for (Size i = 1; i < hill.length; ++i)
        {
          double mz_prev = hill.mz_values[i - 1];
          double mz_curr = hill.mz_values[i];
          hill_mass_diffs->push_back(calculatePPM(mz_curr, mz_prev));
        }
      }
    }
  }

  return hills;
}

vector<Biosaur2Algorithm::Hill> Biosaur2Algorithm::processHills(const vector<Hill>& hills, Size min_length) const
{
  vector<Hill> processed;
  for (const auto& hill : hills)
  {
    if (hill.length >= min_length)
    {
      Hill processed_hill = hill;
      processed_hill.drift_time_median = calculateMedian(hill.drift_times);
      processed_hill.ion_mobility_median = calculateMedian(hill.ion_mobilities);
      processed.push_back(processed_hill);
    }
  }
  return processed;
}

vector<Biosaur2Algorithm::Hill> Biosaur2Algorithm::splitHills(const vector<Hill>& hills, double hvf, Size min_length) const
{
  vector<Hill> split_hills;

  for (const auto& hill : hills)
  {
    if (hill.length <= 2)
    {
      split_hills.push_back(hill);
      continue;
    }

    vector<double> smoothed = meanFilter(hill.intensities, max<Size>(3, hill.length / 5));
    Size start_idx = 0;
    bool splitting = false;

    for (Size i = 1; i + 1 < smoothed.size(); ++i)
    {
      if (smoothed[i] * hvf < smoothed[i - 1] && smoothed[i] * hvf < smoothed[i + 1])
      {
        splitting = true;
        if (i - start_idx + 1 >= min_length)
        {
          Hill new_hill = hill;
          new_hill.scan_indices.clear();
          new_hill.peak_indices.clear();
          new_hill.mz_values.clear();
          new_hill.intensities.clear();
          new_hill.rt_values.clear();
          new_hill.drift_times.clear();
          new_hill.ion_mobilities.clear();

          for (Size k = start_idx; k <= i; ++k)
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
          new_hill.mz_median = calculateMedian(new_hill.mz_values);
          new_hill.intensity_sum = accumulate(new_hill.intensities.begin(), new_hill.intensities.end(), 0.0);
          auto max_it = max_element(new_hill.intensities.begin(), new_hill.intensities.end());
          Size apex_idx = distance(new_hill.intensities.begin(), max_it);
          new_hill.intensity_apex = *max_it;
          new_hill.rt_apex = new_hill.rt_values[apex_idx];
          split_hills.push_back(new_hill);
        }
        start_idx = i;
      }
    }

    if (!splitting)
    {
      split_hills.push_back(hill);
    }
  }

  return split_hills;
}

Size Biosaur2Algorithm::checkIsotopeValleySplit(const vector<IsotopeCandidate>& isotopes,
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

  vector<double> smoothed = meanFilter(isotope_intensities, 3);
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

vector<Biosaur2Algorithm::PeptideFeature> Biosaur2Algorithm::detectIsotopePatterns(vector<Hill>& hills,
                                                                                   double itol_ppm,
                                                                                   int min_charge,
                                                                                   int max_charge,
                                                                                   bool negative_mode,
                                                                                   double ivf,
                                                                                   int iuse,
                                                                                   bool enable_isotope_calib) const
{
  vector<PeptideFeature> features;
  set<Size> used_hills;

  OPENMS_LOG_INFO << "Detecting isotope patterns..." << endl;

  sort(hills.begin(), hills.end(),
       [](const Hill& a, const Hill& b) { return a.mz_median < b.mz_median; });
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
      double mono_mz = mono_hill.mz_median;

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
            if (hills[j].mz_median > expected_mz + mz_tolerance)
            {
              break;
            }

            double diff = fabs(hills[j].mz_median - expected_mz);
            if (diff <= mz_tolerance)
            {
              if (mono_hill.length >= 3)
              {
                double mass_diff_ppm = calculatePPM(hills[j].mz_median, expected_mz);
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
        auto calib = calibrateMass(isotope_errors[ic]);
        isotope_calib_map[ic] = calib;
        OPENMS_LOG_INFO << "Isotope " << ic << " calibration: shift="
                        << calib.first << " ppm, sigma=" << calib.second << " ppm" << endl;
      }
    }

    for (int ic = 4; ic <= 9; ++ic)
    {
      if (isotope_errors[ic].size() >= 1000)
      {
        isotope_calib_map[ic] = calibrateMass(isotope_errors[ic]);
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

  for (Size i = 0; i < hills.size(); ++i)
  {
    if (used_hills.find(hills[i].hill_idx) != used_hills.end())
    {
      continue;
    }

    const Hill& mono_hill = hills[i];
    double mono_mz = mono_hill.mz_median;

    for (int charge = max_charge; charge >= min_charge; --charge)
    {
      double mz_spacing = ISOTOPE_MASSDIFF / charge;
      vector<IsotopeCandidate> isotopes;
      bool pattern_valid = true;

      for (int iso_num = 1; iso_num <= 9; ++iso_num)
      {
        double expected_mz = mono_mz + iso_num * mz_spacing;
        double mz_tolerance = expected_mz * itol_ppm * 1e-6;
        bool found = false;

        for (Size j = i + 1; j < hills.size(); ++j)
        {
          if (hills[j].mz_median > expected_mz + mz_tolerance)
          {
            break;
          }
          if (used_hills.find(hills[j].hill_idx) != used_hills.end())
          {
            continue;
          }

          double diff = fabs(hills[j].mz_median - expected_mz);
          if (diff <= mz_tolerance)
          {
            double cos_corr = cosineCorrelation(mono_hill.intensities, mono_hill.scan_indices,
                                                hills[j].intensities, hills[j].scan_indices);
            if (cos_corr >= 0.6)
            {
              IsotopeCandidate cand;
              cand.hill_idx = hills[j].hill_idx;
              cand.isotope_number = iso_num;
              cand.mass_diff_ppm = calculatePPM(hills[j].mz_median, expected_mz);
              cand.cos_corr = cos_corr;
              isotopes.push_back(cand);
              found = true;
              break;
            }
          }
        }

        if (!found && iso_num == 1)
        {
          pattern_valid = false;
          break;
        }
        else if (!found)
        {
          break;
        }
      }

      if (pattern_valid && !isotopes.empty())
      {
        if (enable_isotope_calib)
        {
          vector<IsotopeCandidate> filtered;
          for (const auto& cand : isotopes)
          {
            auto calib = isotope_calib_map[cand.isotope_number];
            if (fabs(cand.mass_diff_ppm - calib.first) <= 5.0 * calib.second)
            {
              filtered.push_back(cand);
            }
            else
            {
              break;
            }
          }
          isotopes = filtered;
        }

        if (!isotopes.empty())
        {
          Size keep = checkIsotopeValleySplit(isotopes, hills, ivf);
          if (keep < isotopes.size())
          {
            isotopes.resize(keep);
          }
        }

        if (!isotopes.empty())
        {
          PeptideFeature feature;
          feature.mz = mono_mz;
          feature.rt_start = mono_hill.rt_start;
          feature.rt_end = mono_hill.rt_end;
          feature.rt_apex = mono_hill.rt_apex;
          feature.intensity_apex = mono_hill.intensity_apex;
          feature.intensity_sum = mono_hill.intensity_sum;

          if (iuse != 0)
          {
            int isotopes_to_add = (iuse == -1) ? static_cast<int>(isotopes.size())
                                               : min(static_cast<int>(isotopes.size()), iuse);
            for (int iso_idx = 0; iso_idx < isotopes_to_add; ++iso_idx)
            {
              for (const auto& hill : hills)
              {
                if (hill.hill_idx == isotopes[iso_idx].hill_idx)
                {
                  feature.intensity_apex += hill.intensity_apex;
                  feature.intensity_sum += hill.intensity_sum;
                  break;
                }
              }
            }
          }

          feature.charge = charge;
          feature.n_isotopes = isotopes.size() + 1;
          feature.n_scans = mono_hill.length;
          feature.isotopes = isotopes;
          feature.mono_hill_idx = mono_hill.hill_idx;
          feature.drift_time = mono_hill.drift_time_median;
          feature.ion_mobility = mono_hill.ion_mobility_median;

          double proton_mass = Constants::PROTON_MASS_U;
          if (negative_mode)
          {
            feature.mass_calib = mono_mz * charge + proton_mass * charge;
          }
          else
          {
            feature.mass_calib = mono_mz * charge - proton_mass * charge;
          }

          features.push_back(feature);
          used_hills.insert(mono_hill.hill_idx);
          for (const auto& iso : isotopes)
          {
            used_hills.insert(iso.hill_idx);
          }

          break;
        }
      }
    }
  }

  OPENMS_LOG_INFO << "Detected " << features.size() << " features with isotope patterns" << endl;
  return features;
}

FeatureMap Biosaur2Algorithm::convertToFeatureMap(const vector<PeptideFeature>& features) const
{
  FeatureMap feature_map;

  for (const auto& f : features)
  {
    Feature feature;
    feature.setMZ(f.mz);
    feature.setRT(f.rt_apex);
    feature.setIntensity(f.intensity_apex);
    feature.setCharge(f.charge);
    feature.setOverallQuality(f.n_isotopes);

    ConvexHull2D hull;
    vector<DPosition<2>> hull_points;
    hull_points.emplace_back(f.rt_start, f.mz);
    hull_points.emplace_back(f.rt_end, f.mz);
    hull.setHullPoints(hull_points);
    feature.getConvexHulls().push_back(hull);

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

double Biosaur2Algorithm::cosineCorrelation(const vector<double>& intensities1,
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
  out << "hill_idx\tmz\trtStart\trtEnd\trtApex\tintensityApex\tintensitySum\tnScans" << endl;

  for (const auto& hill : hills)
  {
    out << hill.hill_idx << "\t"
        << hill.mz_median << "\t"
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
