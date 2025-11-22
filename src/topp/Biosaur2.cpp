// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2023.
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
// $Authors: Mark Ivanov, Timo Sachsenberg $
// --------------------------------------------------------------------------

/**
 * @page TOPP_Biosaur2 Biosaur2
 *
 * @brief Feature detection for LC-MS1 data
 *
 * This TOPP tool is a C++ reimplementation of the Biosaur2 feature detection algorithm.
 * It detects peptide features in centroided LC-MS1 data by:
 * 1. Grouping peaks across scans into "hills"
 * 2. Splitting hills at valley points
 * 3. Detecting isotope patterns
 * 4. Calculating feature properties (m/z, RT, intensity, charge)
 *
 * Reference:
 * Abdrakhimov, et al. Biosaur: An open-source Python software for liquid 
 * chromatography-mass spectrometry peptide feature detection with ion mobility support.
 * https://doi.org/10.1002/rcm.9045
 *
 * <B>The command line parameters of this tool are:</B>
 * @verbinclude TOPP_Biosaur2.cli
 * <B>INI file documentation of this tool:</B>
 * @htmlinclude TOPP_Biosaur2.html
 */

#include <OpenMS/APPLICATIONS/TOPPBase.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/DATASTRUCTURES/ConvexHull2D.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/IsotopeDistribution.h>
#include <OpenMS/CHEMISTRY/ISOTOPEDISTRIBUTION/CoarseIsotopePatternGenerator.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/MATH/STATISTICS/BasicStatistics.h>
#include <OpenMS/PROCESSING/CENTROIDING/PeakPickerHiRes.h>
#include <OpenMS/IONMOBILITY/FAIMSHelper.h>

#include <algorithm>
#include <map>
#include <vector>
#include <cmath>
#include <numeric>
#include <set>

#ifdef _OPENMP
#include <omp.h>
#endif

using namespace OpenMS;
using namespace std;

//-------------------------------------------------------------
// Doxygen docu
//-------------------------------------------------------------

/**
  @page TOPP_Biosaur2 Biosaur2

  @brief Feature detection for LC-MS1 data using the Biosaur2 algorithm.
*/

// We do not want this class to show up in the docu:
/// @cond TOPPCLASSES

class TOPPBiosaur2 :
  public TOPPBase
{
public:
  TOPPBiosaur2() :
    TOPPBase("Biosaur2", "Feature detection for LC-MS1 data", false)
  {
  }

protected:
  // Structure to represent a peak across multiple scans (a "hill")
  struct Hill
  {
    vector<Size> scan_indices;     // Indices of scans containing this hill
    vector<Size> peak_indices;     // Peak indices within each scan
    vector<double> mz_values;      // m/z values for each peak
    vector<double> intensities;    // Intensity values for each peak
    vector<double> rt_values;      // RT values for each scan
    vector<double> drift_times;    // Drift times (FAIMS compensation voltage)
    vector<double> ion_mobilities; // Ion mobility values (PASEF/IM data)
    double mz_median;              // Median m/z value
    double rt_start;               // Start RT
    double rt_end;                 // End RT
    double rt_apex;                // RT at apex
    double intensity_apex;         // Intensity at apex
    double intensity_sum;          // Sum of intensities
    double drift_time_median;      // Median drift time (FAIMS)
    double ion_mobility_median;    // Median ion mobility (PASEF/IM)
    Size length;                   // Number of scans
    Size hill_idx;                 // Unique hill identifier
  };

  // Structure to represent an isotope candidate
  struct IsotopeCandidate
  {
    Size hill_idx;                 // Hill index
    Size isotope_number;           // Isotope number (0=mono, 1=first 13C, etc)
    double mass_diff_ppm;          // Mass difference in ppm
    double cos_corr;               // Cosine correlation of RT profiles
  };

  // Structure to represent a feature
  struct PeptideFeature
  {
    double mz;                     // Monoisotopic m/z
    double rt_start;               // Start RT
    double rt_end;                 // End RT
    double rt_apex;                // RT at apex
    double intensity_apex;         // Intensity at apex
    double intensity_sum;          // Sum of intensities
    int charge;                    // Charge state
    Size n_isotopes;               // Number of isotopes
    Size n_scans;                  // Number of scans
    double mass_calib;             // Calibrated neutral mass
    double drift_time;             // Drift time (FAIMS compensation voltage)
    double ion_mobility;           // Ion mobility (PASEF/IM data)
    vector<IsotopeCandidate> isotopes; // Isotope information
    Size mono_hill_idx;            // Monoisotopic hill index
  };

  void registerOptionsAndFlags_() override
  {
    registerInputFile_("in", "<file>", "", "Input mzML file (centroided data)");
    setValidFormats_("in", ListUtils::create<String>("mzML"));
    
    registerOutputFile_("out", "<file>", "", "Output featureXML file");
    setValidFormats_("out", ListUtils::create<String>("featureXML"));
    
    registerOutputFile_("out_tsv", "<file>", "", "Optional: output TSV file (Biosaur2 format)", false);
    setValidFormats_("out_tsv", ListUtils::create<String>("tsv"));

    // Parameters matching Python biosaur2
    registerDoubleOption_("mini", "<value>", 1.0, "Minimum intensity threshold", false);
    setMinFloat_("mini", 0.0);
    
    registerDoubleOption_("minmz", "<value>", 350.0, "Minimum m/z value", false);
    setMinFloat_("minmz", 0.0);
    
    registerDoubleOption_("maxmz", "<value>", 1500.0, "Maximum m/z value", false);
    setMinFloat_("maxmz", 0.0);
    
    registerDoubleOption_("htol", "<value>", 8.0, "Mass accuracy in ppm for combining peaks into hills", false);
    setMinFloat_("htol", 0.0);
    
    registerDoubleOption_("itol", "<value>", 8.0, "Mass accuracy in ppm for isotopic patterns", false);
    setMinFloat_("itol", 0.0);
    
    registerDoubleOption_("hvf", "<value>", 1.3, "Hill valley factor for splitting hills", false);
    setMinFloat_("hvf", 1.0);
    
    registerDoubleOption_("ivf", "<value>", 5.0, "Isotope valley factor for splitting isotope patterns", false);
    setMinFloat_("ivf", 1.0);
    
    registerIntOption_("minlh", "<value>", 2, "Minimum number of scans for a hill", false);
    setMinInt_("minlh", 1);
    
    registerIntOption_("cmin", "<value>", 1, "Minimum charge state", false);
    setMinInt_("cmin", 1);
    
    registerIntOption_("cmax", "<value>", 6, "Maximum charge state", false);
    setMinInt_("cmax", 1);
    
    registerIntOption_("threads", "<value>", 1, "Number of threads for parallel processing (0 = auto-detect)", false);
    setMinInt_("threads", 0);
    
    registerIntOption_("iuse", "<value>", 0, "Number of isotopes for intensity calculation (0=mono only, -1=all, 1=mono+first, etc.)", false);
    setMinInt_("iuse", -1);
    
    registerFlag_("nm", "Negative mode (affects neutral mass calculation)", false);
    
    registerFlag_("tof", "Enable TOF-specific intensity filtering", false);
    
    registerFlag_("profile", "Enable profile mode processing (centroid spectra using PeakPickerHiRes)", false);
    
    registerFlag_("use_hill_calib", "Enable automatic hill mass tolerance calibration", false);
    
    registerFlag_("ignore_iso_calib", "Disable automatic isotope mass error calibration", false);
    
    registerFlag_("write_hills", "Write intermediate hills to TSV file", false);
    
    registerOutputFile_("out_hills", "<file>", "", "Optional: output hills TSV file", false);
    setValidFormats_("out_hills", ListUtils::create<String>("tsv"));
  }

  // Calculate ppm difference between two m/z values
  double calculatePPM(double mz1, double mz2)
  {
    if (fabs(mz2) < 1e-10) return 0.0; // Guard against division by zero
    return (mz1 - mz2) / mz2 * 1e6;
  }

  // Calculate median of a vector
  double calculateMedian(const vector<double>& values)
  {
    if (values.empty()) return 0.0;
    
    vector<double> sorted = values;
    sort(sorted.begin(), sorted.end());
    size_t n = sorted.size();
    
    if (n % 2 == 0)
    {
      return (sorted[n/2 - 1] + sorted[n/2]) / 2.0;
    }
    else
    {
      return sorted[n/2];
    }
  }

  // Simple moving average filter
  vector<double> meanFilter(const vector<double>& data, Size window)
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
      result[i] = sum / (end - start);
    }
    
    return result;
  }

  // Calibrate mass errors using histogram and Gaussian fitting
  // Returns pair of (mass_shift, mass_sigma)
  pair<double, double> calibrateMass(const vector<double>& mass_errors, double bin_width = 0.05)
  {
    if (mass_errors.empty()) return make_pair(0.0, 10.0);
    
    // Find range
    double min_error = *min_element(mass_errors.begin(), mass_errors.end());
    double max_error = *max_element(mass_errors.begin(), mass_errors.end());
    double mass_left = -min_error;
    double mass_right = max_error;
    
    // Create histogram
    int n_bins = static_cast<int>((mass_left + mass_right) / bin_width);
    if (n_bins < 5) return make_pair(0.0, 10.0); // Not enough bins
    
    vector<double> bin_centers;
    vector<int> bin_counts(n_bins, 0);
    
    for (int i = 0; i < n_bins; ++i)
    {
      bin_centers.push_back(-mass_left + (i + 0.5) * bin_width);
    }
    
    // Fill histogram
    for (double error : mass_errors)
    {
      int bin = static_cast<int>((error + mass_left) / bin_width);
      if (bin >= 0 && bin < n_bins)
      {
        bin_counts[bin]++;
      }
    }
    
    // Simple Gaussian fitting using moments
    // Calculate weighted mean and std from histogram
    double sum_x = 0.0, sum_x2 = 0.0, sum_w = 0.0;
    for (size_t i = 0; i < bin_centers.size(); ++i)
    {
      double x = bin_centers[i];
      double w = bin_counts[i];
      sum_x += w * x;
      sum_x2 += w * x * x;
      sum_w += w;
    }
    
    if (sum_w < 10) return make_pair(0.0, 10.0); // Not enough data
    
    double mean = sum_x / sum_w;
    double variance = (sum_x2 / sum_w) - (mean * mean);
    double sigma = sqrt(max(variance, 0.01)); // Ensure positive variance
    
    // Sanity checks
    if (fabs(mean) >= max(mass_left, mass_right))
    {
      // Mean is too extreme, try with wider bins
      return calibrateMass(mass_errors, 0.25);
    }
    
    if (isinf(sigma) || isnan(sigma))
    {
      return make_pair(0.0, 10.0);
    }
    
    return make_pair(mean, sigma);
  }

  // TOF-specific intensity filtering
  // Estimates noise distribution and filters low-intensity peaks
  void processTOF(MSExperiment& exp, double min_mz, double max_mz)
  {
    OPENMS_LOG_INFO << "Applying TOF-specific intensity filtering..." << endl;
    
    // Group m/z values into bins for noise estimation
    const double mz_bin_size = 50.0;
    map<int, vector<double>> intensity_bins;
    
    // Sample first 25 spectra to estimate noise distribution
    Size sample_size = min(Size(25), exp.size());
    for (Size i = 0; i < sample_size; ++i)
    {
      for (Size j = 0; j < exp[i].size(); ++j)
      {
        double mz = exp[i][j].getMZ();
        if (mz >= min_mz && mz <= max_mz)
        {
          int bin = static_cast<int>(mz / mz_bin_size);
          double log_intensity = log10(exp[i][j].getIntensity());
          intensity_bins[bin].push_back(log_intensity);
        }
      }
    }
    
    // Calculate threshold for each bin (mean + 2*std)
    map<int, double> bin_thresholds;
    for (auto& bin_pair : intensity_bins)
    {
      if (bin_pair.second.size() >= 150)
      {
        vector<double>& intensities = bin_pair.second;
        double sum = accumulate(intensities.begin(), intensities.end(), 0.0);
        double mean = sum / intensities.size();
        
        double sq_sum = 0.0;
        for (double val : intensities)
        {
          sq_sum += (val - mean) * (val - mean);
        }
        double std_dev = sqrt(sq_sum / intensities.size());
        
        bin_thresholds[bin_pair.first] = pow(10.0, mean + 2.0 * std_dev);
      }
    }
    
    // Apply filtering to all spectra
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
        
        double threshold = 150.0; // Default threshold
        if (bin_thresholds.find(bin) != bin_thresholds.end())
        {
          threshold = bin_thresholds[bin];
        }
        
        if (intensity >= threshold)
        {
          filtered_spectrum.push_back(spectrum[i]);
        }
      }
      
      spectrum = filtered_spectrum;
      total_peaks_after += spectrum.size();
    }
    
    OPENMS_LOG_INFO << "TOF filtering: " << total_peaks_before 
                    << " peaks -> " << total_peaks_after << " peaks" << endl;
  }

  // Centroid profile spectra using PeakPickerHiRes
  void centroidProfileSpectra(MSExperiment& exp)
  {
    OPENMS_LOG_INFO << "Centroiding profile spectra using PeakPickerHiRes..." << endl;
    
    PeakPickerHiRes picker;
    Param picker_param = picker.getParameters();
    
    // Configure peak picker for high-resolution data
    picker_param.setValue("signal_to_noise", 0.0); // Disable S/N filtering (we filter by intensity)
    picker.setParameters(picker_param);
    
    MSExperiment centroided_exp;
    Size total_peaks_before = 0;
    Size total_peaks_after = 0;
    
    for (Size i = 0; i < exp.size(); ++i)
    {
      total_peaks_before += exp[i].size();
      
      MSSpectrum centroided_spectrum;
      
      // Check if spectrum is already centroided
      if (exp[i].getType() == SpectrumSettings::CENTROID)
      {
        centroided_spectrum = exp[i];
      }
      else
      {
        // Pick peaks from profile spectrum
        picker.pick(exp[i], centroided_spectrum);
        centroided_spectrum.setRT(exp[i].getRT());
        centroided_spectrum.setMSLevel(exp[i].getMSLevel());
        centroided_spectrum.setType(SpectrumSettings::CENTROID);
        
        // Copy drift time if present (for FAIMS)
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

  // Cosine correlation between two RT profiles
  double cosineCorrelation(const vector<double>& intensities1, 
                          const vector<Size>& scans1,
                          const vector<double>& intensities2,
                          const vector<Size>& scans2)
  {
    // Build maps for efficient lookup
    map<Size, double> map1, map2;
    for (Size i = 0; i < scans1.size(); ++i)
    {
      map1[scans1[i]] = intensities1[i];
    }
    for (Size i = 0; i < scans2.size(); ++i)
    {
      map2[scans2[i]] = intensities2[i];
    }
    
    // Calculate cosine correlation
    double dot_product = 0.0;
    double norm1 = 0.0;
    double norm2 = 0.0;
    
    for (const auto& p1 : map1)
    {
      Size scan = p1.first;
      double i1 = p1.second;
      
      if (map2.find(scan) != map2.end())
      {
        double i2 = map2[scan];
        dot_product += i1 * i2;
      }
      norm1 += i1 * i1;
    }
    
    for (const auto& p2 : map2)
    {
      norm2 += p2.second * p2.second;
    }
    
    if (norm1 == 0.0 || norm2 == 0.0) return 0.0;
    
    return dot_product / (sqrt(norm1) * sqrt(norm2));
  }

  // Detect hills (peaks grouped across scans)
  // If collect_mass_diffs is true, returns mass differences for calibration
  vector<Hill> detectHills(const MSExperiment& exp, 
                          double htol_ppm,
                          double min_intensity,
                          double min_mz,
                          double max_mz,
                          vector<double>* mass_diffs = nullptr)
  {
    vector<Hill> all_hills;
    map<Size, Hill> active_hills; // Map from hill_idx to hill
    Size next_hill_idx = 0;
    
    OPENMS_LOG_INFO << "Detecting hills across " << exp.size() << " MS1 spectra..." << endl;
    
    for (Size scan_idx = 0; scan_idx < exp.size(); ++scan_idx)
    {
      const MSSpectrum& spectrum = exp[scan_idx];
      
      // Only process MS1 spectra
      if (spectrum.getMSLevel() != 1) continue;
      
      double rt = spectrum.getRT();
      
      // Create new hills map for this scan
      map<Size, Hill> new_active_hills;
      set<Size> matched_peaks;
      
      // Try to extend existing hills
      for (auto& hill_pair : active_hills)
      {
        Hill& hill = hill_pair.second;
        double target_mz = hill.mz_median;
        double mz_tolerance = target_mz * htol_ppm * 1e-6;
        
        // Find closest peak within tolerance
        double best_diff = mz_tolerance + 1.0;
        Size best_peak_idx = 0;
        bool found = false;
        
        for (Size peak_idx = 0; peak_idx < spectrum.size(); ++peak_idx)
        {
          // Skip peaks outside filter range
          double peak_mz = spectrum[peak_idx].getMZ();
          if (peak_mz < min_mz || peak_mz > max_mz) continue;
          if (spectrum[peak_idx].getIntensity() < min_intensity) continue;
          
          // Skip already matched peaks
          if (matched_peaks.find(peak_idx) != matched_peaks.end()) continue;
          
          double diff = abs(peak_mz - target_mz);
          if (diff <= mz_tolerance && diff < best_diff)
          {
            best_diff = diff;
            best_peak_idx = peak_idx;
            found = true;
          }
        }
        
        if (found)
        {
          // Extend hill
          matched_peaks.insert(best_peak_idx);
          
          // Collect mass difference for calibration if requested
          if (mass_diffs != nullptr)
          {
            double mass_diff_ppm = calculatePPM(spectrum[best_peak_idx].getMZ(), target_mz);
            mass_diffs->push_back(mass_diff_ppm);
          }
          
          hill.scan_indices.push_back(scan_idx);
          hill.peak_indices.push_back(best_peak_idx);
          hill.mz_values.push_back(spectrum[best_peak_idx].getMZ());
          hill.intensities.push_back(spectrum[best_peak_idx].getIntensity());
          hill.rt_values.push_back(rt);
          
          // Collect drift time if available (FAIMS)
          if (spectrum.getDriftTime() >= 0)
          {
            hill.drift_times.push_back(spectrum.getDriftTime());
          }
          
          // Collect ion mobility if available (PASEF/IM data)
          const auto& fda = spectrum.getFloatDataArrays();
          auto it_im = getDataArrayByName(fda, Constants::UserParam::ION_MOBILITY);
          if (it_im != fda.end() && best_peak_idx < it_im->size())
          {
            hill.ion_mobilities.push_back((*it_im)[best_peak_idx]);
          }
          
          // Update median m/z
          vector<double> mz_copy = hill.mz_values;
          hill.mz_median = calculateMedian(mz_copy);
          
          new_active_hills[hill.hill_idx] = hill;
        }
        else
        {
          // Hill ended, save it
          all_hills.push_back(hill);
        }
      }
      
      // Start new hills from unmatched peaks
      for (Size peak_idx = 0; peak_idx < spectrum.size(); ++peak_idx)
      {
        double peak_mz = spectrum[peak_idx].getMZ();
        double peak_int = spectrum[peak_idx].getIntensity();
        
        // Apply filters
        if (peak_mz < min_mz || peak_mz > max_mz) continue;
        if (peak_int < min_intensity) continue;
        if (matched_peaks.find(peak_idx) != matched_peaks.end()) continue;
        
        // Create new hill
        Hill new_hill;
        new_hill.hill_idx = next_hill_idx++;
        new_hill.scan_indices.push_back(scan_idx);
        new_hill.peak_indices.push_back(peak_idx);
        new_hill.mz_values.push_back(peak_mz);
        new_hill.intensities.push_back(peak_int);
        new_hill.rt_values.push_back(rt);
        new_hill.mz_median = peak_mz;
        
        // Collect drift time if available (FAIMS)
        if (spectrum.getDriftTime() >= 0)
        {
          new_hill.drift_times.push_back(spectrum.getDriftTime());
        }
        
        // Collect ion mobility if available (PASEF/IM data)
        const auto& fda = spectrum.getFloatDataArrays();
        auto it_im = getDataArrayByName(fda, Constants::UserParam::ION_MOBILITY);
        if (it_im != fda.end() && peak_idx < it_im->size())
        {
          new_hill.ion_mobilities.push_back((*it_im)[peak_idx]);
        }
        
        new_active_hills[new_hill.hill_idx] = new_hill;
      }
      
      active_hills = new_active_hills;
    }
    
    // Save remaining active hills
    for (auto& hill_pair : active_hills)
    {
      all_hills.push_back(hill_pair.second);
    }
    
    OPENMS_LOG_INFO << "Detected " << all_hills.size() << " initial hills" << endl;
    
    return all_hills;
  }

  // Process hills: calculate properties and filter by minimum length
  vector<Hill> processHills(vector<Hill>& hills, Size min_length)
  {
    vector<Hill> processed_hills;
    
    for (auto& hill : hills)
    {
      hill.length = hill.scan_indices.size();
      
      // Filter by minimum length
      if (hill.length < min_length) continue;
      
      // Calculate RT properties
      hill.rt_start = hill.rt_values.front();
      hill.rt_end = hill.rt_values.back();
      
      // Find apex
      Size apex_idx = 0;
      double max_intensity = 0.0;
      for (Size i = 0; i < hill.intensities.size(); ++i)
      {
        if (hill.intensities[i] > max_intensity)
        {
          max_intensity = hill.intensities[i];
          apex_idx = i;
        }
      }
      hill.rt_apex = hill.rt_values[apex_idx];
      hill.intensity_apex = max_intensity;
      
      // Calculate sum intensity
      hill.intensity_sum = accumulate(hill.intensities.begin(), hill.intensities.end(), 0.0);
      
      // Calculate median drift time if available (FAIMS)
      if (!hill.drift_times.empty())
      {
        vector<double> drift_copy = hill.drift_times;
        hill.drift_time_median = calculateMedian(drift_copy);
      }
      else
      {
        hill.drift_time_median = -1.0; // No drift time
      }
      
      // Calculate median ion mobility if available (PASEF/IM)
      if (!hill.ion_mobilities.empty())
      {
        vector<double> im_copy = hill.ion_mobilities;
        hill.ion_mobility_median = calculateMedian(im_copy);
      }
      else
      {
        hill.ion_mobility_median = -1.0; // No ion mobility
      }
      
      processed_hills.push_back(hill);
    }
    
    OPENMS_LOG_INFO << "Processed " << processed_hills.size() << " hills (after filtering by min length)" << endl;
    
    return processed_hills;
  }

  // Split hills at valley points
  vector<Hill> splitHills(vector<Hill>& hills, double hvf, Size min_length)
  {
    vector<Hill> split_hills;
    Size next_hill_idx = 0;
    
    for (auto& hill : hills)
    {
      if (hill.length < min_length * 2)
      {
        // Too short to split, keep as is
        hill.hill_idx = next_hill_idx++;
        split_hills.push_back(hill);
        continue;
      }
      
      // Apply smoothing
      vector<double> smoothed = meanFilter(hill.intensities, 3);
      
      // Find local minima
      vector<Size> split_points;
      
      for (Size i = min_length - 1; i < hill.length - min_length; ++i)
      {
        // Check if this is a valley
        double valley_int = smoothed[i];
        if (valley_int == 0.0) continue;
        
        // Find max intensity on left and right
        double left_max = *max_element(smoothed.begin(), smoothed.begin() + i);
        double right_max = *max_element(smoothed.begin() + i + 1, smoothed.end());
        
        // Check valley criteria
        if (left_max / valley_int >= hvf && right_max / valley_int >= hvf)
        {
          // Check if not too close to previous split point
          if (split_points.empty() || i >= split_points.back() + min_length)
          {
            split_points.push_back(i);
          }
        }
      }
      
      if (split_points.empty())
      {
        // No split needed
        hill.hill_idx = next_hill_idx++;
        split_hills.push_back(hill);
      }
      else
      {
        // Split into multiple hills
        Size start_idx = 0;
        for (Size split_idx : split_points)
        {
          Hill new_hill;
          new_hill.hill_idx = next_hill_idx++;
          
          // Copy data from start to split point
          for (Size i = start_idx; i <= split_idx; ++i)
          {
            new_hill.scan_indices.push_back(hill.scan_indices[i]);
            new_hill.peak_indices.push_back(hill.peak_indices[i]);
            new_hill.mz_values.push_back(hill.mz_values[i]);
            new_hill.intensities.push_back(hill.intensities[i]);
            new_hill.rt_values.push_back(hill.rt_values[i]);
          }
          
          // Recalculate properties
          vector<double> mz_copy = new_hill.mz_values;
          new_hill.mz_median = calculateMedian(mz_copy);
          new_hill.length = new_hill.scan_indices.size();
          
          if (new_hill.length >= min_length)
          {
            split_hills.push_back(new_hill);
          }
          
          start_idx = split_idx + 1;
        }
        
        // Add remaining part
        if (start_idx < hill.length)
        {
          Hill new_hill;
          new_hill.hill_idx = next_hill_idx++;
          
          for (Size i = start_idx; i < hill.length; ++i)
          {
            new_hill.scan_indices.push_back(hill.scan_indices[i]);
            new_hill.peak_indices.push_back(hill.peak_indices[i]);
            new_hill.mz_values.push_back(hill.mz_values[i]);
            new_hill.intensities.push_back(hill.intensities[i]);
            new_hill.rt_values.push_back(hill.rt_values[i]);
          }
          
          vector<double> mz_copy = new_hill.mz_values;
          new_hill.mz_median = calculateMedian(mz_copy);
          new_hill.length = new_hill.scan_indices.size();
          
          if (new_hill.length >= min_length)
          {
            split_hills.push_back(new_hill);
          }
        }
      }
    }
    
    // Recalculate properties for all split hills
    vector<Hill> final_hills = processHills(split_hills, min_length);
    
    OPENMS_LOG_INFO << "After splitting: " << final_hills.size() << " hills" << endl;
    
    return final_hills;
  }

  // Check if isotope pattern should be split based on valley detection
  // Returns the number of isotopes to keep (0 means reject the pattern)
  Size checkIsotopeValleySplit(const vector<IsotopeCandidate>& isotopes,
                                const vector<Hill>& hills,
                                double ivf)
  {
    if (isotopes.empty() || ivf <= 1.0) return isotopes.size();
    
    // Build intensity profile across isotopes
    vector<double> isotope_intensities;
    isotope_intensities.push_back(1.0); // Monoisotopic (normalized to 1.0)
    
    for (const auto& iso : isotopes)
    {
      // Find the hill for this isotope
      for (const auto& hill : hills)
      {
        if (hill.hill_idx == iso.hill_idx)
        {
          isotope_intensities.push_back(hill.intensity_apex);
          break;
        }
      }
    }
    
    if (isotope_intensities.size() < 3) return isotope_intensities.size() - 1;
    
    // Normalize intensities
    double max_intensity = *max_element(isotope_intensities.begin(), isotope_intensities.end());
    for (auto& intensity : isotope_intensities)
    {
      intensity /= max_intensity;
    }
    
    // Look for valleys starting from position 4 (or max position + 1)
    Size max_pos = 0;
    for (Size i = 0; i < isotope_intensities.size(); ++i)
    {
      if (isotope_intensities[i] > isotope_intensities[max_pos])
      {
        max_pos = i;
      }
    }
    
    Size min_check_pos = max(Size(4), max_pos + 1);
    
    // Check for valleys
    for (Size i = min_check_pos; i < isotope_intensities.size() - 1; ++i)
    {
      double local_min = isotope_intensities[i];
      double right_max = *max_element(isotope_intensities.begin() + i + 1, 
                                       isotope_intensities.end());
      
      if (local_min * ivf < right_max)
      {
        // Found a valley, split here
        return i;
      }
    }
    
    return isotope_intensities.size() - 1;
  }

  // Detect isotope patterns and create features
  vector<PeptideFeature> detectIsotopePatterns(vector<Hill>& hills,
                                               double itol_ppm,
                                               int min_charge,
                                               int max_charge,
                                               bool negative_mode,
                                               double ivf,
                                               int iuse,
                                               bool enable_isotope_calib)
  {
    vector<PeptideFeature> features;
    set<Size> used_hills;
    
    OPENMS_LOG_INFO << "Detecting isotope patterns..." << endl;
    
    // Sort hills by m/z for efficient searching
    sort(hills.begin(), hills.end(), 
         [](const Hill& a, const Hill& b) { return a.mz_median < b.mz_median; });
    
    // Use OpenMS constant for C13-C12 mass difference (needed in both passes)
    const double ISOTOPE_MASSDIFF = Constants::C13C12_MASSDIFF_U;
    
    // Isotope mass error calibration map
    // isotope_calib_map[isotope_number] = pair(mass_shift, mass_sigma)
    map<int, pair<double, double>> isotope_calib_map;
    
    // Initialize with default values
    for (int ic = 1; ic <= 9; ++ic)
    {
      isotope_calib_map[ic] = make_pair(0.0, itol_ppm);
    }
    
    // First pass: collect isotope mass errors for calibration if enabled
    if (enable_isotope_calib)
    {
      OPENMS_LOG_INFO << "Performing isotope calibration..." << endl;
      
      // Collect mass errors for each isotope position
      map<int, vector<double>> isotope_errors;
      for (int ic = 1; ic <= 9; ++ic)
      {
        isotope_errors[ic] = vector<double>();
      }
      
      // Quick pass to collect isotopes
      for (Size i = 0; i < hills.size(); ++i)
      {
        const Hill& mono_hill = hills[i];
        double mono_mz = mono_hill.mz_median;
        
        // Try different charge states
        for (int charge = max_charge; charge >= min_charge; --charge)
        {
          double mz_spacing = ISOTOPE_MASSDIFF / charge;
          
          // Look for first 9 isotopes
          bool found_first = false;
          for (int iso_num = 1; iso_num <= 9; ++iso_num)
          {
            double expected_mz = mono_mz + iso_num * mz_spacing;
            double mz_tolerance = expected_mz * itol_ppm * 1e-6;
            
            // Search for matching hill
            for (Size j = i + 1; j < hills.size(); ++j)
            {
              if (hills[j].mz_median > expected_mz + mz_tolerance) break;
              
              double diff = abs(hills[j].mz_median - expected_mz);
              if (diff <= mz_tolerance)
              {
                // Check if this is a good match (sufficient scans)
                if (mono_hill.length >= 3)
                {
                  double mass_diff_ppm = calculatePPM(hills[j].mz_median, expected_mz);
                  isotope_errors[iso_num].push_back(mass_diff_ppm);
                  
                  if (iso_num == 1) found_first = true;
                }
                break;
              }
            }
          }
          
          if (found_first) break; // Found valid charge state
        }
      }
      
      // Calibrate for isotopes 1-3 (most reliable)
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
      
      // Extrapolate for higher isotopes
      for (int ic = 4; ic <= 9; ++ic)
      {
        if (isotope_errors[ic].size() >= 1000)
        {
          auto calib = calibrateMass(isotope_errors[ic]);
          isotope_calib_map[ic] = calib;
        }
        else if (ic > 1 && isotope_calib_map.find(ic-1) != isotope_calib_map.end())
        {
          // Extrapolate from previous isotope
          auto prev = isotope_calib_map[ic-1];
          auto prev2 = isotope_calib_map.find(ic-2) != isotope_calib_map.end() ? 
                       isotope_calib_map[ic-2] : make_pair(0.0, itol_ppm);
          
          double shift_delta = prev.first - prev2.first;
          double sigma_ratio = prev.second / max(prev2.second, 0.1);
          
          isotope_calib_map[ic] = make_pair(prev.first + shift_delta, prev.second * sigma_ratio);
        }
      }
      
      OPENMS_LOG_INFO << "Isotope 1 calibration: shift=" << isotope_calib_map[1].first 
                      << " ppm, sigma=" << isotope_calib_map[1].second << " ppm" << endl;
    }
    
    // Second pass (or only pass): detect features with calibrated tolerances
    sort(hills.begin(), hills.end(), 
         [](const Hill& a, const Hill& b) { return a.mz_median < b.mz_median; });
    
    for (Size i = 0; i < hills.size(); ++i)
    {
      // Skip if already used
      if (used_hills.find(hills[i].hill_idx) != used_hills.end()) continue;
      
      const Hill& mono_hill = hills[i];
      double mono_mz = mono_hill.mz_median;
      
      // Try different charge states
      for (int charge = max_charge; charge >= min_charge; --charge)
      {
        double mz_spacing = ISOTOPE_MASSDIFF / charge;
        
        vector<IsotopeCandidate> isotopes;
        bool pattern_valid = true;
        
        // Look for isotopes
        for (int iso_num = 1; iso_num <= 9; ++iso_num)
        {
          double expected_mz = mono_mz + iso_num * mz_spacing;
          double mz_tolerance = expected_mz * itol_ppm * 1e-6;
          
          // Search for matching hill
          bool found = false;
          for (Size j = i + 1; j < hills.size(); ++j)
          {
            // Stop if we're too far
            if (hills[j].mz_median > expected_mz + mz_tolerance) break;
            
            // Skip if already used
            if (used_hills.find(hills[j].hill_idx) != used_hills.end()) continue;
            
            double diff = abs(hills[j].mz_median - expected_mz);
            if (diff <= mz_tolerance)
            {
              // Check RT overlap
              double cos_corr = cosineCorrelation(
                mono_hill.intensities, mono_hill.scan_indices,
                hills[j].intensities, hills[j].scan_indices
              );
              
              if (cos_corr >= 0.6)
              {
                IsotopeCandidate candidate;
                candidate.hill_idx = hills[j].hill_idx;
                candidate.isotope_number = iso_num;
                candidate.mass_diff_ppm = calculatePPM(hills[j].mz_median, expected_mz);
                candidate.cos_corr = cos_corr;
                
                isotopes.push_back(candidate);
                found = true;
                break;
              }
            }
          }
          
          if (!found && iso_num == 1)
          {
            // Must have at least first isotope
            pattern_valid = false;
            break;
          }
          else if (!found)
          {
            // No more isotopes found
            break;
          }
        }
        
        // Create feature if we have a valid pattern
        if (pattern_valid && !isotopes.empty())
        {
          // Apply isotope calibration filtering if enabled
          if (enable_isotope_calib)
          {
            vector<IsotopeCandidate> filtered_isotopes;
            for (const auto& cand : isotopes)
            {
              auto calib = isotope_calib_map[cand.isotope_number];
              double mass_shift = calib.first;
              double mass_sigma = calib.second;
              
              // Keep isotope if within 5*sigma of calibrated shift
              if (fabs(cand.mass_diff_ppm - mass_shift) <= 5.0 * mass_sigma)
              {
                filtered_isotopes.push_back(cand);
              }
              else
              {
                break; // Stop at first rejected isotope
              }
            }
            isotopes = filtered_isotopes;
          }
          
          // Check if pattern should be split using ivf
          if (!isotopes.empty())
          {
            Size isotopes_to_keep = checkIsotopeValleySplit(isotopes, hills, ivf);
            
            // Trim isotopes if needed
            if (isotopes_to_keep < isotopes.size())
            {
              isotopes.resize(isotopes_to_keep);
            }
          }
          
          // Only create feature if we still have isotopes
          if (!isotopes.empty())
          {
            PeptideFeature feature;
            feature.mz = mono_mz;
            feature.rt_start = mono_hill.rt_start;
            feature.rt_end = mono_hill.rt_end;
            feature.rt_apex = mono_hill.rt_apex;
            feature.intensity_apex = mono_hill.intensity_apex;
            feature.intensity_sum = mono_hill.intensity_sum;
            
            // Add intensities from isotopes based on iuse parameter
            if (iuse != 0)
            {
              int isotopes_to_add = (iuse == -1) ? isotopes.size() : min(static_cast<int>(isotopes.size()), iuse);
              for (int iso_idx = 0; iso_idx < isotopes_to_add; ++iso_idx)
              {
                // Find the hill for this isotope
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
            feature.n_isotopes = isotopes.size() + 1; // +1 for monoisotope
            feature.n_scans = mono_hill.length;
            feature.isotopes = isotopes;
            feature.mono_hill_idx = mono_hill.hill_idx;
            feature.drift_time = mono_hill.drift_time_median; // FAIMS drift time
            feature.ion_mobility = mono_hill.ion_mobility_median; // IM/PASEF ion mobility

            // Calculate neutral mass
            double proton_mass = 1.007276;
            if (negative_mode)
            {
              feature.mass_calib = mono_mz * charge + proton_mass * charge;
            }
            else
            {
              feature.mass_calib = mono_mz * charge - proton_mass * charge;
            }
            
            features.push_back(feature);
            
            // Mark hills as used
            used_hills.insert(mono_hill.hill_idx);
            for (const auto& iso : isotopes)
            {
              used_hills.insert(iso.hill_idx);
            }
            
            break; // Found valid charge state, move to next monoisotopic peak
          }
        }
      }
    }
    
    OPENMS_LOG_INFO << "Detected " << features.size() << " features with isotope patterns" << endl;
    
    return features;
  }

  // Write features to TSV file (Biosaur2 format)
  void writeTSV(const vector<PeptideFeature>& features, const String& filename)
  {
    ofstream out(filename);
    
    // Write header
    out << "massCalib\trtApex\tintensityApex\tintensitySum\tcharge\t"
        << "nIsotopes\tnScans\tmz\trtStart\trtEnd\tFAIMS\tIM" << endl;
    
    // Write features
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
    
    out.close();
    OPENMS_LOG_INFO << "Wrote " << features.size() << " features to TSV file: " << filename << endl;
  }

  // Convert features to OpenMS FeatureMap
  FeatureMap convertToFeatureMap(const vector<PeptideFeature>& features)
  {
    FeatureMap feature_map;
    
    for (const auto& f : features)
    {
      Feature feature;
      feature.setMZ(f.mz);
      feature.setRT(f.rt_apex);
      feature.setIntensity(f.intensity_apex);
      feature.setCharge(f.charge);
      
      // Set quality (use number of isotopes as quality metric)
      feature.setOverallQuality(f.n_isotopes);
      
      // Create convex hull
      ConvexHull2D hull;
      vector<DPosition<2>> hull_points;
      hull_points.push_back(DPosition<2>(f.rt_start, f.mz));
      hull_points.push_back(DPosition<2>(f.rt_end, f.mz));
      hull.setHullPoints(hull_points);
      feature.getConvexHulls().push_back(hull);
      
      // Add metadata
      feature.setMetaValue("mass_calib", f.mass_calib);
      feature.setMetaValue("n_isotopes", f.n_isotopes);
      feature.setMetaValue("n_scans", f.n_scans);
      feature.setMetaValue("intensity_sum", f.intensity_sum);
      
      // Add FAIMS drift time if available
      if (f.drift_time >= 0)
      {
        feature.setMetaValue("FAIMS_compensation_voltage", f.drift_time);
      }
      
      // Add ion mobility if available (PASEF/IM)
      if (f.ion_mobility >= 0)
      {
        feature.setMetaValue("ion_mobility", f.ion_mobility);
      }
      
      feature_map.push_back(feature);
    }
    
    // Set protein identifications (empty)
    feature_map.getProteinIdentifications().resize(1);
    
    return feature_map;
  }

  // Write hills to TSV file
  void writeHills(const vector<Hill>& hills, const String& filename)
  {
    ofstream out(filename);
    
    // Write header
    out << "hill_idx\tmz\trtStart\trtEnd\trtApex\tintensityApex\tintensitySum\tnScans" << endl;
    
    // Write hills
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
    
    out.close();
    OPENMS_LOG_INFO << "Wrote " << hills.size() << " hills to: " << filename << endl;
  }

  ExitCodes main_(int, const char**) override
  {
    //-------------------------------------------------------------
    // Parsing parameters
    //-------------------------------------------------------------
    String in = getStringOption_("in");
    String out = getStringOption_("out");
    String out_tsv = getStringOption_("out_tsv");
    String out_hills = getStringOption_("out_hills");
    
    double mini = getDoubleOption_("mini");
    double minmz = getDoubleOption_("minmz");
    double maxmz = getDoubleOption_("maxmz");
    double htol = getDoubleOption_("htol");
    double itol = getDoubleOption_("itol");
    double hvf = getDoubleOption_("hvf");
    double ivf = getDoubleOption_("ivf");
    Size minlh = getIntOption_("minlh");
    int cmin = getIntOption_("cmin");
    int cmax = getIntOption_("cmax");
    int threads = getIntOption_("threads");
    int iuse = getIntOption_("iuse");
    bool negative_mode = getFlag_("nm");
    bool tof_mode = getFlag_("tof");
    bool profile_mode = getFlag_("profile");
    bool use_hill_calib = getFlag_("use_hill_calib");
    bool ignore_iso_calib = getFlag_("ignore_iso_calib");
    bool write_hills = getFlag_("write_hills");
    
    // Set number of threads for OpenMP
#ifdef _OPENMP
    if (threads == 0)
    {
      threads = omp_get_max_threads();
    }
    omp_set_num_threads(threads);
    OPENMS_LOG_INFO << "Using " << threads << " threads for parallel processing" << endl;
#else
    if (threads > 1)
    {
      OPENMS_LOG_WARN << "OpenMP not available, using single thread" << endl;
    }
#endif
    
    //-------------------------------------------------------------
    // Reading input
    //-------------------------------------------------------------
    OPENMS_LOG_INFO << "Loading input file: " << in << endl;
    
    MSExperiment exp;
    MzMLFile mzml_file;
    mzml_file.load(in, exp);
    
    // Filter to MS1 only
    exp.getSpectra().erase(
      remove_if(exp.begin(), exp.end(), 
                [](const MSSpectrum& s) { return s.getMSLevel() != 1; }),
      exp.end()
    );
    
    OPENMS_LOG_INFO << "Loaded " << exp.size() << " MS1 spectra" << endl;
    
    if (exp.empty())
    {
      OPENMS_LOG_ERROR << "No MS1 spectra found in input file!" << endl;
      return ILLEGAL_PARAMETERS;
    }
    
    // Check for FAIMS data
    auto cvs = FAIMSHelper::getCompensationVoltages(exp);
    if (!cvs.empty())
    {
      OPENMS_LOG_INFO << "Detected FAIMS data with " << cvs.size() << " compensation voltage(s)" << endl;
      for (const auto& cv : cvs)
      {
        OPENMS_LOG_INFO << "  CV: " << cv << " V" << endl;
      }
    }
    
    // Check for ion mobility (PASEF/IM) data
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
    
    //-------------------------------------------------------------
    // Pre-processing
    //-------------------------------------------------------------
    
    // Centroid profile spectra if requested
    if (profile_mode)
    {
      centroidProfileSpectra(exp);
    }
    
    // Apply TOF processing if requested
    if (tof_mode)
    {
      processTOF(exp, minmz, maxmz);
    }
    
    //-------------------------------------------------------------
    // Feature detection
    //-------------------------------------------------------------
    
    // Hill calibration if requested
    double calibrated_htol = htol;
    if (use_hill_calib)
    {
      OPENMS_LOG_INFO << "Performing hill mass tolerance calibration..." << endl;
      
      vector<double> mass_diffs;
      Size sample_size = min(exp.size(), Size(1000));
      Size start_idx = (exp.size() > 1000) ? (exp.size() / 2 - 500) : 0;
      
      // Create a subset of the experiment for calibration
      MSExperiment calib_exp;
      for (Size i = start_idx; i < start_idx + sample_size && i < exp.size(); ++i)
      {
        calib_exp.addSpectrum(exp[i]);
      }
      
      // Detect hills and collect mass differences
      vector<Hill> calib_hills = detectHills(calib_exp, htol, mini, minmz, maxmz, &mass_diffs);
      
      if (!mass_diffs.empty())
      {
        auto calib = calibrateMass(mass_diffs);
        double calibrated_sigma = calib.second;
        
        // Update htol to be the minimum of current value and 5*sigma
        calibrated_htol = min(htol, 5.0 * calibrated_sigma);
        
        OPENMS_LOG_INFO << "Automatically optimized htol parameter: " 
                        << calibrated_htol << " ppm (was " << htol << " ppm)" << endl;
      }
    }
    
    // Step 1: Detect hills with calibrated tolerance
    vector<Hill> hills = detectHills(exp, calibrated_htol, mini, minmz, maxmz);
    
    // Step 2: Process hills (calculate properties, filter by length)
    hills = processHills(hills, minlh);
    
    // Step 3: Split hills at valleys
    hills = splitHills(hills, hvf, minlh);
    
    // Write hills if requested
    if (write_hills || !out_hills.empty())
    {
      String hills_file = out_hills;
      if (hills_file.empty())
      {
        String base = out;
        Size dot_pos = base.find_last_of('.');
        if (dot_pos != String::npos)
        {
          base = base.substr(0, dot_pos);
        }
        hills_file = base + ".hills.tsv";
      }
      writeHills(hills, hills_file);
    }
    
    // Step 4: Detect isotope patterns with calibration
    bool enable_isotope_calib = !ignore_iso_calib;
    vector<PeptideFeature> features = detectIsotopePatterns(hills, itol, cmin, cmax, negative_mode, ivf, iuse, enable_isotope_calib);
    
    //-------------------------------------------------------------
    // Writing output
    //-------------------------------------------------------------
    
    // Write FeatureXML
    FeatureMap feature_map = convertToFeatureMap(features);
    FeatureXMLFile feature_file;
    feature_file.store(out, feature_map);
    OPENMS_LOG_INFO << "Wrote " << features.size() << " features to: " << out << endl;
    
    // Write TSV if requested
    if (!out_tsv.empty())
    {
      writeTSV(features, out_tsv);
    }
    
    return EXECUTION_OK;
  }
};

int main(int argc, const char** argv)
{
  TOPPBiosaur2 tool;
  return tool.main(argc, argv);
}

/// @endcond
