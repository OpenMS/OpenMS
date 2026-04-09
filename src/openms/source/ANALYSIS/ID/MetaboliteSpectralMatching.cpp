// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Erhan Kenar $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>

#include <OpenMS/CONCEPT/Constants.h>

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <numeric>
#include <cmath>
#include <cctype>
#include <limits>
#include <boost/math/special_functions/factorials.hpp>

#include <boost/dynamic_bitset.hpp>

#include <OpenMS/PROCESSING/SPECTRAMERGING/SpectraMerger.h>
#include <OpenMS/PROCESSING/FILTERING/WindowMower.h>

using namespace std;

namespace OpenMS
{
  namespace
  {
    /**
     * @brief Calculate bonus for consecutive ion series in fragment annotations
     * 
     * Detects consecutive fragment ions of the same type (e.g., w3, w4, w5)
     * which provide strong evidence for correct sequence assignment.
     * 
     * @param annotations Vector of peak annotations containing ion type and position
     * @return Bonus score based on longest consecutive series found
     */
    double calculateConsecutiveSeriesBonus(const vector<PeptideHit::PeakAnnotation>& annotations)
    {
      if (annotations.empty()) return 0.0;
      
      // Group annotations by ion type and extract positions
      map<char, set<int>> ion_positions;
      
      for (const auto& ann : annotations)
      {
        String annotation = ann.annotation;  // Copy, not reference - avoid potential dangling reference
        if (annotation.empty()) continue;
        
        // Skip special cases (like 'precursor' or empty annotations)
        if (!std::isalpha(static_cast<unsigned char>(annotation[0]))) continue;
        
        char ion_type = annotation[0];
        
        // Extract position number from annotation (e.g., "w3" -> 3)
        annotation.erase(0, 1);  // Remove first character
        annotation.trim();  // Remove any whitespace
        
        //Remove charge indicators
        size_t plus_pos = annotation.find('+');
        if (plus_pos != String::npos)
        {
          annotation = annotation.substr(0, plus_pos);
        }
        size_t minus_pos = annotation.find('-');
        if (minus_pos != String::npos)
        {
          annotation = annotation.substr(0, minus_pos);
        }
        
        if (annotation.empty()) continue;
        
        int position;
        try
        {
          position = annotation.toInt();
        }
        catch (const Exception::ConversionError&)
        {
          throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
            "Failed to parse fragment ion annotation '" + ann.annotation + 
            "' - expected ion type followed by position number. Internal error - annotations should be well-formed.",
            ann.annotation);
        }
        ion_positions[ion_type].insert(position);
      }
      
      // Find the longest consecutive sequence for each ion type
      int max_consecutive = 0;
      
      for (const auto& ion_pair : ion_positions)
      {
        const set<int>& positions = ion_pair.second;
        if (positions.size() < 2) continue;
        
        int current_consecutive = 1;
        int prev_position = -1; // Use -1 to indicate no previous position
        
        for (int pos : positions)
        {
          if (prev_position >= 0 && pos == prev_position + 1)
          {
            current_consecutive++;
          }
          else if (prev_position >= 0)
          {
            // Broke the consecutive streak, update max and restart
            max_consecutive = max(max_consecutive, current_consecutive);
            current_consecutive = 1;
          }
          prev_position = pos;
        }
        
        // Don't forget to check the final consecutive run
        max_consecutive = max(max_consecutive, current_consecutive);
      }
      
      // Calculate multiplicative bonus factor
      // Consecutive series provide strong evidence, so boost score multiplicatively
      // Returns: 1.0 (no bonus) to ~2.0 (very long series)
      if (max_consecutive >= 2)
      {
        // 10% boost per consecutive ion, capped to avoid over-weighting
        double bonus_factor = 1.0 + (0.1 * max_consecutive);
        return min(bonus_factor, 2.0); // Cap at 2x multiplier
      }
      
      return 1.0; // No bonus = 1.0 multiplier (identity)
    }
  }

  SpectralMatch::SpectralMatch() :
    observed_precursor_mass_(),
    observed_precursor_rt_(),
    found_precursor_mass_(),
    found_precursor_charge_(),
    matching_score_(),
    observed_spectrum_idx_(),
    matching_spectrum_idx_(),
    observed_spectrum_native_id_(),
    primary_id_(),
    secondary_id_(),
    common_name_(),
    sum_formula_(),
    inchi_string_(),
    smiles_string_(),
    precursor_adduct_()
  {
  }

  /// Default destructor
  SpectralMatch::~SpectralMatch() = default;

  /// Copy constructor
  SpectralMatch::SpectralMatch(const SpectralMatch& sm) = default;

  /// Assignment operator
  SpectralMatch& SpectralMatch::operator=(const SpectralMatch& rhs)
  {
    if (this == &rhs) return *this;

    observed_precursor_mass_ = rhs.observed_precursor_mass_;
    observed_precursor_rt_ = rhs.observed_precursor_rt_;
    found_precursor_mass_ = rhs.found_precursor_mass_;
    found_precursor_charge_ = rhs.found_precursor_charge_;
    matching_score_ = rhs.matching_score_;
    observed_spectrum_idx_ = rhs.observed_spectrum_idx_;
    matching_spectrum_idx_ = rhs.matching_spectrum_idx_;
    observed_spectrum_native_id_ = rhs.observed_spectrum_native_id_;
    primary_id_ = rhs.primary_id_;
    secondary_id_ = rhs.secondary_id_;
    common_name_ = rhs.common_name_;
    sum_formula_ = rhs.sum_formula_;
    inchi_string_ = rhs.inchi_string_;
    smiles_string_ = rhs.smiles_string_;
    precursor_adduct_ = rhs.precursor_adduct_;

    return *this;
  }


  double SpectralMatch::getObservedPrecursorMass() const
  {
    return observed_precursor_mass_;
  }


  void SpectralMatch::setObservedPrecursorMass(const double& qmass)
  {
    observed_precursor_mass_ = qmass;
  }


  double SpectralMatch::getObservedPrecursorRT() const
  {
    return observed_precursor_rt_;
  }


  void SpectralMatch::setObservedPrecursorRT(const double& prt)
  {
    observed_precursor_rt_ = prt;
  }


  double SpectralMatch::getFoundPrecursorMass() const
  {
    return found_precursor_mass_;
  }


  void SpectralMatch::setFoundPrecursorMass(const double& fmass)
  {
    found_precursor_mass_ = fmass;
  }


  Int SpectralMatch::getFoundPrecursorCharge() const
  {
    return found_precursor_charge_;
  }


  void SpectralMatch::setFoundPrecursorCharge(const Int& pch)
  {
    found_precursor_charge_ = pch;
  }


  double SpectralMatch::getMatchingScore() const
  {
    return matching_score_;
  }


  void SpectralMatch::setMatchingScore(const double& mscore)
  {
    matching_score_ = mscore;
  }


  Size SpectralMatch::getObservedSpectrumIndex() const
  {
    return observed_spectrum_idx_;
  }


  void SpectralMatch::setObservedSpectrumIndex(const Size& obs_spec_idx)
  {
    observed_spectrum_idx_ = obs_spec_idx;
  }


  Size SpectralMatch::getMatchingSpectrumIndex() const
  {
    return matching_spectrum_idx_;
  }


  void SpectralMatch::setMatchingSpectrumIndex(const Size& match_spec_idx)
  {
    matching_spectrum_idx_ = match_spec_idx;
  }


  String SpectralMatch::getObservedSpectrumNativeID() const
  {
    return observed_spectrum_native_id_;
  }


  void SpectralMatch::setObservedSpectrumNativeID(const String& obs_spec_native_id)
  {
    observed_spectrum_native_id_ = obs_spec_native_id;
  }

  String SpectralMatch::getPrimaryIdentifier() const
  {
    return primary_id_;
  }


  void SpectralMatch::setPrimaryIdentifier(const String& pid)
  {
    primary_id_ = pid;
  }


  String SpectralMatch::getSecondaryIdentifier() const
  {
    return secondary_id_;
  }


  void SpectralMatch::setSecondaryIdentifier(const String& sid)
  {
    secondary_id_ = sid;
  }


  String SpectralMatch::getCommonName() const
  {
    return common_name_;
  }


  void SpectralMatch::setCommonName(const String& cname)
  {
    common_name_ = cname;
  }


  String SpectralMatch::getSumFormula() const
  {
    return sum_formula_;
  }


  void SpectralMatch::setSumFormula(const String& sf)
  {
    sum_formula_ = sf;
  }


  String SpectralMatch::getInchiString() const
  {
    return inchi_string_;
  }


  void SpectralMatch::setInchiString(const String& istr)
  {
    inchi_string_ = istr;
  }


  String SpectralMatch::getSMILESString() const
  {
    return smiles_string_;
  }


  void SpectralMatch::setSMILESString(const String& sstr)
  {
    smiles_string_ = sstr;
  }


  String SpectralMatch::getPrecursorAdduct() const
  {
    return precursor_adduct_;
  }


  void SpectralMatch::setPrecursorAdduct(const String& padd)
  {
    precursor_adduct_ = padd;
  }


  MetaboliteSpectralMatching::MetaboliteSpectralMatching() :
    DefaultParamHandler("MetaboliteSpectralMatching"), ProgressLogger()
  {
    defaults_.setValue("prec_mass_error_value", 100.0, "Error allowed for precursor ion mass.");
    defaults_.setValue("frag_mass_error_value", 500.0, "Error allowed for product ions.");

    defaults_.setValue("mass_error_unit", "ppm", "Unit of mass error (ppm or Da)");
    defaults_.setValidStrings("mass_error_unit", {"ppm","Da"});

    defaults_.setValue("report_mode", "top3", "Which results shall be reported: the top-three scoring ones or the best scoring one?");
    defaults_.setValidStrings("report_mode", {"top3","best"});

    defaults_.setValue("ionization_mode", "positive", "Positive or negative ionization mode?");
    defaults_.setValidStrings("ionization_mode", {"positive","negative"});

    defaults_.setValue("merge_spectra", "true", "Merge MS2 spectra with the same precursor mass.");
    defaults_.setValidStrings("merge_spectra", {"true","false"});

    defaultsToParam_();

    this->setLogType(CMD);
  }


  MetaboliteSpectralMatching::~MetaboliteSpectralMatching() = default;


  double MetaboliteSpectralMatching::computeHyperScore(
    double fragment_mass_error,
    bool fragment_mass_tolerance_unit_ppm,
    const MSSpectrum& exp_spectrum,
    const MSSpectrum& db_spectrum,
    double mz_lower_bound)
  {
    return computeHyperScore_(fragment_mass_error,
                              fragment_mass_tolerance_unit_ppm, exp_spectrum,
                              db_spectrum, nullptr, mz_lower_bound);
  }


  double MetaboliteSpectralMatching::computeHyperScore(
    double fragment_mass_error,
    bool fragment_mass_tolerance_unit_ppm,
    const MSSpectrum& exp_spectrum,
    const MSSpectrum& db_spectrum,
    vector<PeptideHit::PeakAnnotation>& annotations,
    double mz_lower_bound)
  {
    return computeHyperScore_(fragment_mass_error,
                              fragment_mass_tolerance_unit_ppm, exp_spectrum,
                              db_spectrum, &annotations, mz_lower_bound);
  }


  double MetaboliteSpectralMatching::computeHyperScore_(
    double fragment_mass_error,
    bool fragment_mass_tolerance_unit_ppm,
    const MSSpectrum& exp_spectrum,
    const MSSpectrum& db_spectrum,
    vector<PeptideHit::PeakAnnotation>* annotations,
    double mz_lower_bound)
  {
    if (exp_spectrum.empty()) return 0;

    // define m/z range to consider:
    double min_exp_mz = exp_spectrum[0].getMZ(); // lowest experimental m/z
    double mz_offset = fragment_mass_error;
    if (fragment_mass_tolerance_unit_ppm)
    {
      mz_offset = min_exp_mz * fragment_mass_error * 1e-6;
    }
    mz_lower_bound = max(mz_lower_bound, min_exp_mz - mz_offset);
    double max_exp_mz = exp_spectrum.back().getMZ(); // highest experimental m/z
    if (fragment_mass_tolerance_unit_ppm)
    {
      mz_offset = max_exp_mz * fragment_mass_error * 1e-6;
    }
    double mz_upper_bound = max_exp_mz + mz_offset;

    // for every DB (theoretical) peak in the valid m/z range, find the closest
    // matching experimental (observed) peak within the allowed tolerance;
    // in principle, multiple DB peaks can match to the same exp. peak:
    map<Size, vector<MSSpectrum::ConstIterator>> peak_matches;
    for (auto db_it = db_spectrum.MZBegin(mz_lower_bound);
         db_it != db_spectrum.MZEnd(mz_upper_bound); ++db_it)
    {
      double db_mz = db_it->getMZ();

      if (fragment_mass_tolerance_unit_ppm)
      {
        mz_offset = db_mz * fragment_mass_error * 1e-6;
      }

      Int index = exp_spectrum.findNearest(db_mz, mz_offset);
      if (index >= 0) peak_matches[index].push_back(db_it);
    }

    double dot_product = 0.0;
    for (const auto& match : peak_matches)
    {
      double db_intensity = 0.0;
      double best_mass_error = std::numeric_limits<double>::max();
      for (const auto& db_it : match.second)
      {
        db_intensity = max(db_intensity, double(db_it->getIntensity()));
        // Track the smallest mass error for this experimental peak
        double mass_error = abs(exp_spectrum[match.first].getMZ() - db_it->getMZ());
        best_mass_error = min(best_mass_error, mass_error);
      }
      
      // Calculate mass accuracy weight (1.0 = perfect match, 0.0 = at tolerance edge)
      double tolerance = fragment_mass_tolerance_unit_ppm ? 
                        exp_spectrum[match.first].getMZ() * fragment_mass_error * 1e-6 :
                        fragment_mass_error;
      double mass_accuracy_weight = 1.0 - (best_mass_error / tolerance);
      mass_accuracy_weight = max(0.0, mass_accuracy_weight); // Ensure non-negative
      
      dot_product += mass_accuracy_weight * db_intensity * exp_spectrum[match.first].getIntensity();
    }

    // return annotations for matching peaks?
    if ((annotations != nullptr) &&
        !db_spectrum.getStringDataArrays().empty() &&
        !db_spectrum.getIntegerDataArrays().empty())
    {
      for (const auto& match : peak_matches)
      {
        const auto& exp_peak = exp_spectrum[match.first];
        // potentially add several annotations for the same peak if there are
        // multiple matches for that peak:
        for (const auto& db_it : match.second)
        {
          PeptideHit::PeakAnnotation ann;
          Size index = db_it - db_spectrum.begin();
          ann.annotation = db_spectrum.getStringDataArrays()[0].at(index);
          ann.charge = db_spectrum.getIntegerDataArrays()[0].at(index);
          ann.mz = exp_peak.getMZ();
          ann.intensity = exp_peak.getIntensity();
          annotations->push_back(ann);
        }
      }
    }

    Size matched_ions_count = peak_matches.size(); // count obs. peaks only once
    double matched_ions_term = 0.0;

    // return score 0 if too few matched ions
    if (matched_ions_count < 3)
    {
      return matched_ions_term;
    }

    if (matched_ions_count <= boost::math::max_factorial<double>::value)
    {
      matched_ions_term = log(boost::math::factorial<double>(matched_ions_count));
    }
    else
    {
      matched_ions_term = log(boost::math::factorial<double>(boost::math::max_factorial<double>::value));
    }

    double hyperscore = log(dot_product) + matched_ions_term;
    
    // Apply consecutive ion series multiplier if we have annotations
    if (annotations != nullptr && !annotations->empty())
    {
      double series_multiplier = calculateConsecutiveSeriesBonus(*annotations);
      hyperscore *= series_multiplier;
    }
    
    if (hyperscore < 0) hyperscore = 0;

    return hyperscore;
  }


  double MetaboliteSpectralMatching::computeMVHScore(
    double fragment_mass_error,
    bool fragment_mass_tolerance_unit_ppm,
    const MSSpectrum& exp_spectrum,
    const MSSpectrum& db_spectrum,
    Size num_intensity_classes,
    double tic_fraction,
    double mz_lower_bound)
  {
    return computeMVHScore_(fragment_mass_error,
                            fragment_mass_tolerance_unit_ppm, exp_spectrum,
                            db_spectrum, num_intensity_classes, tic_fraction,
                            nullptr, mz_lower_bound);
  }


  double MetaboliteSpectralMatching::computeMVHScore(
    double fragment_mass_error,
    bool fragment_mass_tolerance_unit_ppm,
    const MSSpectrum& exp_spectrum,
    const MSSpectrum& db_spectrum,
    vector<PeptideHit::PeakAnnotation>& annotations,
    Size num_intensity_classes,
    double tic_fraction,
    double mz_lower_bound)
  {
    return computeMVHScore_(fragment_mass_error,
                            fragment_mass_tolerance_unit_ppm, exp_spectrum,
                            db_spectrum, num_intensity_classes, tic_fraction,
                            &annotations, mz_lower_bound);
  }


  double MetaboliteSpectralMatching::computeMVHScore_(
    double fragment_mass_error,
    bool fragment_mass_tolerance_unit_ppm,
    const MSSpectrum& exp_spectrum,
    const MSSpectrum& db_spectrum,
    Size num_intensity_classes,
    double tic_fraction,
    vector<PeptideHit::PeakAnnotation>* annotations,
    double mz_lower_bound)
  {
    if (exp_spectrum.empty()) return 0;
    if (num_intensity_classes < 1) num_intensity_classes = 1;
    if (num_intensity_classes > 7) num_intensity_classes = 7;
    if (tic_fraction < 0.5) tic_fraction = 0.5;
    if (tic_fraction > 1.0) tic_fraction = 1.0;

    // Step 1: Filter peaks by TIC fraction
    // Sort peaks by intensity to find TIC threshold
    vector<pair<double, Size>> intensity_index_pairs;
    intensity_index_pairs.reserve(exp_spectrum.size());
    double total_intensity = 0.0;
    for (Size i = 0; i < exp_spectrum.size(); ++i)
    {
      double inten = exp_spectrum[i].getIntensity();
      intensity_index_pairs.emplace_back(inten, i);
      total_intensity += inten;
    }
    
    // Sort by intensity (descending)
    sort(intensity_index_pairs.begin(), intensity_index_pairs.end(),
         [](const pair<double, Size>& a, const pair<double, Size>& b) {
           return a.first > b.first;
         });
    
    // Find peaks that account for tic_fraction of total intensity
    double cumulative_intensity = 0.0;
    double target_intensity = total_intensity * tic_fraction;
    vector<Size> retained_peak_indices;
    
    for (const auto& pair : intensity_index_pairs)
    {
      retained_peak_indices.push_back(pair.second);
      cumulative_intensity += pair.first;
      if (cumulative_intensity >= target_intensity) break;
    }
    
    // Sort retained indices by m/z for efficient matching
    sort(retained_peak_indices.begin(), retained_peak_indices.end());
    
    // Step 2: Classify retained peaks into intensity classes
    // Classes have sizes in ratio 1:2:4:8... (geometric progression with ratio 2)
    Size total_peaks = retained_peak_indices.size();
    if (total_peaks == 0) return 0;
    
    // Calculate minimum peaks needed for intensity classes
    Size min_peaks_for_classes = (1 << num_intensity_classes) - 1;
    if (total_peaks < min_peaks_for_classes)
    {
      // Not enough peaks for requested number of classes
      return 0;
    }
    
    // Calculate class sizes (smallest class first)
    vector<Size> class_sizes(num_intensity_classes);
    Size peaks_allocated = 0;
    for (Size i = 0; i < num_intensity_classes; ++i)
    {
      class_sizes[i] = 1 << i; // Powers of 2: 1, 2, 4, 8...
      peaks_allocated += class_sizes[i];
    }
    
    // Distribute remaining peaks to largest class
    if (total_peaks > peaks_allocated)
    {
      class_sizes[num_intensity_classes - 1] += (total_peaks - peaks_allocated);
    }
    
    // Step 3: Assign peaks to classes based on intensity ranking
    // Class 0 = highest intensity (A), Class 1 = medium (B), etc.
    vector<Size> peak_to_class(exp_spectrum.size(), SIZE_MAX);
    
    Size class_idx = 0;
    Size peaks_in_current_class = 0;
    for (Size i = 0; i < retained_peak_indices.size(); ++i)
    {
      Size peak_idx = retained_peak_indices[i];
      
      // Move to next class if current class is full
      while (class_idx < num_intensity_classes && 
             peaks_in_current_class >= class_sizes[class_idx])
      {
        class_idx++;
        peaks_in_current_class = 0;
      }
      
      if (class_idx < num_intensity_classes)
      {
        peak_to_class[peak_idx] = class_idx;
        peaks_in_current_class++;
      }
    }
    
    // Step 4: Match theoretical peaks to experimental peaks
    // Count total m/z locations and matches per class
    double min_exp_mz = exp_spectrum[0].getMZ();
    double mz_offset = fragment_mass_error;
    if (fragment_mass_tolerance_unit_ppm)
    {
      mz_offset = min_exp_mz * fragment_mass_error * 1e-6;
    }
    mz_lower_bound = max(mz_lower_bound, min_exp_mz - mz_offset);
    
    double max_exp_mz = exp_spectrum.back().getMZ();
    if (fragment_mass_tolerance_unit_ppm)
    {
      mz_offset = max_exp_mz * fragment_mass_error * 1e-6;
    }
    double mz_upper_bound = max_exp_mz + mz_offset;
    
    // Count matches per intensity class
    vector<Size> matches_per_class(num_intensity_classes, 0);
    map<Size, vector<MSSpectrum::ConstIterator>> peak_matches;
    double mass_accuracy_sum = 0.0; // Track cumulative mass accuracy for bonus term
    
    for (auto db_it = db_spectrum.MZBegin(mz_lower_bound);
         db_it != db_spectrum.MZEnd(mz_upper_bound); ++db_it)
    {
      double db_mz = db_it->getMZ();
      
      if (fragment_mass_tolerance_unit_ppm)
      {
        mz_offset = db_mz * fragment_mass_error * 1e-6;
      }
      
      Int index = exp_spectrum.findNearest(db_mz, mz_offset);
      if (index >= 0)
      {
        Size exp_idx = static_cast<Size>(index);
        peak_matches[exp_idx].push_back(db_it);
        
        // Count this match in the appropriate intensity class
        Size class_id = peak_to_class[exp_idx];
        if (class_id < num_intensity_classes)
        {
          matches_per_class[class_id]++;
        }
        
        // Calculate mass accuracy weight for this match
        double mass_error = abs(exp_spectrum[exp_idx].getMZ() - db_mz);
        double tolerance = fragment_mass_tolerance_unit_ppm ? 
                          db_mz * fragment_mass_error * 1e-6 :
                          fragment_mass_error;
        double mass_accuracy_weight = 1.0 - (mass_error / tolerance);
        mass_accuracy_sum += max(0.0, mass_accuracy_weight);
      }
    }
    
    // Return annotations if requested
    if ((annotations != nullptr) &&
        !db_spectrum.getStringDataArrays().empty() &&
        !db_spectrum.getIntegerDataArrays().empty())
    {
      for (const auto& match : peak_matches)
      {
        const auto& exp_peak = exp_spectrum[match.first];
        for (const auto& db_it : match.second)
        {
          PeptideHit::PeakAnnotation ann;
          Size index = db_it - db_spectrum.begin();
          ann.annotation = db_spectrum.getStringDataArrays()[0].at(index);
          ann.charge = db_spectrum.getIntegerDataArrays()[0].at(index);
          ann.mz = exp_peak.getMZ();
          ann.intensity = exp_peak.getIntensity();
          annotations->push_back(ann);
        }
      }
    }
    
    // Step 5: Compute MVH score using multivariate hypergeometric distribution
    // Calculate total number of theoretical peaks (M)
    Size total_theoretical_peaks = 0;
    for (auto db_it = db_spectrum.MZBegin(mz_lower_bound);
         db_it != db_spectrum.MZEnd(mz_upper_bound); ++db_it)
    {
      total_theoretical_peaks++;
    }
    
    if (total_theoretical_peaks == 0) return 0;
    
    // Calculate total number of possible m/z locations (T)
    // This is approximate based on mass range and tolerance
    double mass_range = mz_upper_bound - mz_lower_bound;
    double location_resolution = fragment_mass_error;
    if (fragment_mass_tolerance_unit_ppm)
    {
      // Use average m/z for ppm-based resolution estimate
      double avg_mz = (mz_lower_bound + mz_upper_bound) / 2.0;
      location_resolution = avg_mz * fragment_mass_error * 1e-6;
    }
    Size total_locations = static_cast<Size>(mass_range / location_resolution);
    if (total_locations < total_peaks + total_theoretical_peaks)
    {
      total_locations = total_peaks + total_theoretical_peaks;
    }
    
    // Count peaks per class and unmatched theoretical peaks
    vector<Size> peaks_per_class(num_intensity_classes, 0);
    for (Size idx : retained_peak_indices)
    {
      Size class_id = peak_to_class[idx];
      if (class_id < num_intensity_classes)
      {
        peaks_per_class[class_id]++;
      }
    }
    
    Size total_matched = 0;
    for (Size m : matches_per_class) total_matched += m;
    
    Size unmatched_theoretical = total_theoretical_peaks - total_matched;
    Size empty_locations = total_locations - total_peaks;
    
    // Compute log of MVH probability using log-space arithmetic
    // log P = sum(log(C(ti, mi))) - log(C(T, M))
    // where C(n, k) is binomial coefficient "n choose k"
    
    auto log_binomial = [](Size n, Size k) -> double {
      if (k > n) return -std::numeric_limits<double>::infinity();
      if (k == 0 || k == n) return 0.0;
      
      // Use lgamma for numerical stability
      return lgamma(n + 1) - lgamma(k + 1) - lgamma(n - k + 1);
    };
    
    double log_prob = 0.0;
    
    // Add contribution from each intensity class
    for (Size i = 0; i < num_intensity_classes; ++i)
    {
      log_prob += log_binomial(peaks_per_class[i], matches_per_class[i]);
    }
    
    // Add contribution from empty locations
    log_prob += log_binomial(empty_locations, unmatched_theoretical);
    
    // Subtract normalization term
    log_prob -= log_binomial(total_locations, total_theoretical_peaks);
    
    // Return negative log probability (higher = better match)
    double mvh_score = -log_prob;
    
    // Add mass accuracy bonus term (scaled to be comparable to log probability)
    // Average mass accuracy per match, scaled by number of matches
    Size total_matches = 0;
    for (Size m : matches_per_class) total_matches += m;
    if (total_matches > 0)
    {
      double avg_mass_accuracy = mass_accuracy_sum / total_matches;
      mvh_score += avg_mass_accuracy * log(1.0 + total_matches);
    }
    
    // Apply consecutive ion series multiplier if we have annotations
    if (annotations != nullptr && !annotations->empty())
    {
      double series_multiplier = calculateConsecutiveSeriesBonus(*annotations);
      mvh_score *= series_multiplier;
    }
    
    // Ensure non-negative score
    if (mvh_score < 0) mvh_score = 0;
    
    return mvh_score;
  }


  void MetaboliteSpectralMatching::run(PeakMap& msexp, PeakMap& spec_db, MzTab& mztab_out, String& out_spectra)
  {
    sort(spec_db.begin(), spec_db.end(), PrecursorMZLess);

    vector<double> mz_keys;

    // copy precursor m/z values to vector for searching
    for (Size spec_idx = 0; spec_idx < spec_db.size(); ++spec_idx)
    {
      mz_keys.push_back(spec_db[spec_idx].getPrecursors()[0].getMZ());
    }

    // remove potential noise peaks by selecting the ten most intense peak per 100 Da window
    WindowMower wm;
    Param wm_param;

    wm_param.setValue("windowsize", 20.0);
    wm_param.setValue("movetype", "slide");
    wm_param.setValue("peakcount", 5);
    wm.setParameters(wm_param);

    wm.filterPeakMap(msexp);

    // merge MS2 spectra with same precursor mass
    if (merge_spectra_)
    {
      SpectraMerger spme;
      spme.mergeSpectraPrecursors(msexp);
      wm.filterPeakMap(msexp);
    }

    // store the spectra if an output file path is given
    if (!out_spectra.empty())
    {
      FileHandler().storeExperiment(out_spectra, msexp, {FileTypes::MZML});
    }


    // container storing results
    vector<SpectralMatch> matching_results;

    bool fragment_error_unit_ppm(true);
    if (mz_error_unit_ == "Da") { fragment_error_unit_ppm = false; }

    for (Size spec_idx = 0; spec_idx < msexp.size(); ++spec_idx)
    {
      OPENMS_LOG_DEBUG << "merged spectrum no. " << spec_idx << " with #fragment ions: " << msexp[spec_idx].size() << endl;

      // iterate over all precursor masses
      for (Size prec_idx = 0; prec_idx < msexp[spec_idx].getPrecursors().size(); ++prec_idx)
      {
        // get precursor m/z
        double precursor_mz(msexp[spec_idx].getPrecursors()[prec_idx].getMZ());

        OPENMS_LOG_DEBUG << "precursor no. " << prec_idx << ": mz " << precursor_mz << " ";

        double prec_mz_lowerbound, prec_mz_upperbound;

        if (!fragment_error_unit_ppm) // Da
        {
          prec_mz_lowerbound = precursor_mz - precursor_mz_error_;
          prec_mz_upperbound = precursor_mz + precursor_mz_error_;
        }
        else // ppm
        {
          double ppm_offset(precursor_mz * 1e-6 * precursor_mz_error_);
          prec_mz_lowerbound = precursor_mz - ppm_offset;
          prec_mz_upperbound = precursor_mz + ppm_offset;
        }

        OPENMS_LOG_DEBUG << "lower mz: " << prec_mz_lowerbound << " ";
        OPENMS_LOG_DEBUG << "upper mz: " << prec_mz_upperbound << endl;

        vector<double>::const_iterator lower_it = lower_bound(mz_keys.begin(), mz_keys.end(), prec_mz_lowerbound);
        vector<double>::const_iterator upper_it = upper_bound(mz_keys.begin(), mz_keys.end(), prec_mz_upperbound);

        Size start_idx(lower_it - mz_keys.begin());
        Size end_idx(upper_it - mz_keys.begin());

        {
          String id_to_log = msexp[spec_idx].metaValueExists("GNPS_Spectrum_ID")
            ? msexp[spec_idx].getMetaValue("GNPS_Spectrum_ID").toString()
            : msexp[spec_idx].getNativeID();
          OPENMS_LOG_DEBUG << "identifying " << id_to_log << endl;
        }

        vector<SpectralMatch> partial_results;

        for (Size search_idx = start_idx; search_idx < end_idx; ++search_idx)
        {
          // do spectral matching
          // Debug: list all available metadata keys
          OPENMS_LOG_DEBUG << "Available metadata keys for spectrum " << search_idx << ":";
          std::vector<String> keys;
          spec_db[search_idx].getKeys(keys);
          for (const auto& key : keys)
          {
            OPENMS_LOG_DEBUG << " " << key;
          }
          OPENMS_LOG_DEBUG << endl;
          
          String metabolite_name = "";
          if (spec_db[search_idx].metaValueExists(Constants::UserParam::MSM_METABOLITE_NAME)) {
            metabolite_name = spec_db[search_idx]
                                .getMetaValue(Constants::UserParam::MSM_METABOLITE_NAME)
                                .toString();
          } else if (spec_db[search_idx].metaValueExists("GNPS_Spectrum_ID")) {
            metabolite_name = spec_db[search_idx]
                                .getMetaValue("GNPS_Spectrum_ID")
                                .toString();
          }

          OPENMS_LOG_DEBUG << "scanning "
                           << spec_db[search_idx].getPrecursors()[0].getMZ()
                           << " " << metabolite_name << endl;
          // check for charge state of precursor ions: do they match?
          if ( (ion_mode_ == "positive" && spec_db[search_idx].getPrecursors()[0].getCharge() < 0) || (ion_mode_ == "negative" && spec_db[search_idx].getPrecursors()[0].getCharge() > 0))
          {
            continue;
          }

          double hyperscore(computeHyperScore(fragment_mz_error_, fragment_error_unit_ppm, msexp[spec_idx], spec_db[search_idx], 0.0));

          OPENMS_LOG_DEBUG << " scored with " << hyperscore << endl;
          if (hyperscore > 0)
          {
            String massbank_id = "";
            String metabolite_name = "";
            
            if (spec_db[search_idx].metaValueExists("GNPS_Spectrum_ID")) {
              massbank_id = spec_db[search_idx].getMetaValue("GNPS_Spectrum_ID").toString();
            }
            if (spec_db[search_idx].metaValueExists("Metabolite_Name")) {
              metabolite_name = spec_db[search_idx].getMetaValue("Metabolite_Name").toString();
            } else if (spec_db[search_idx].metaValueExists("GNPS_Spectrum_ID")) {
              metabolite_name = massbank_id; // Use GNPS_Spectrum_ID as name if no Metabolite_Name
            }
            
            OPENMS_LOG_DEBUG << "  ** detected " << massbank_id << " " << metabolite_name << " scored with " << hyperscore << endl;

            // score result temporarily
            SpectralMatch tmp_match;
            tmp_match.setObservedPrecursorMass(precursor_mz);
            tmp_match.setFoundPrecursorMass(spec_db[search_idx].getPrecursors()[0].getMZ());
            double obs_rt = floor(msexp[spec_idx].getRT() * 10)/10.0;
            tmp_match.setObservedPrecursorRT(obs_rt);
            tmp_match.setFoundPrecursorCharge(spec_db[search_idx].getPrecursors()[0].getCharge());
            tmp_match.setMatchingScore(hyperscore);
            tmp_match.setObservedSpectrumIndex(spec_idx);
            tmp_match.setMatchingSpectrumIndex(search_idx);
            tmp_match.setObservedSpectrumNativeID(msexp[spec_idx].getNativeID());

            String primary_id_value;
            if (spec_db[search_idx].metaValueExists("GNPS_Spectrum_ID"))
            {
              primary_id_value = spec_db[search_idx].getMetaValue("GNPS_Spectrum_ID").toString();
            }
            else if (spec_db[search_idx].metaValueExists("Massbank_Accession_ID"))
            {
              primary_id_value = spec_db[search_idx].getMetaValue("Massbank_Accession_ID").toString();
            }
            else if (spec_db[search_idx].metaValueExists(Constants::UserParam::MSM_METABOLITE_NAME))
            {
              primary_id_value = spec_db[search_idx].getMetaValue(Constants::UserParam::MSM_METABOLITE_NAME).toString();
            }
            else
            {
              primary_id_value = spec_db[search_idx].getNativeID();
            }
            tmp_match.setPrimaryIdentifier(primary_id_value);
            tmp_match.setSecondaryIdentifier(spec_db[search_idx].getMetaValue("HMDB_ID"));
            tmp_match.setSumFormula(spec_db[search_idx].getMetaValue(Constants::UserParam::MSM_SUM_FORMULA));
            tmp_match.setCommonName(metabolite_name);
            tmp_match.setInchiString(spec_db[search_idx].getMetaValue(Constants::UserParam::MSM_INCHI_STRING));
            tmp_match.setSMILESString(spec_db[search_idx].getMetaValue(Constants::UserParam::MSM_SMILES_STRING));
            tmp_match.setPrecursorAdduct(spec_db[search_idx].getMetaValue(Constants::UserParam::MSM_PRECURSOR_ADDUCT));

            partial_results.push_back(tmp_match);
          }
        }

        // sort results by decreasing store
        sort(partial_results.begin(), partial_results.end(), SpectralMatchScoreGreater);

        // report mode: top3 or best?
        if (report_mode_ == "top3")
        {
          Size num_results(partial_results.size());

          Size last_result_idx = (num_results >= 3) ? 3 : num_results;

          for (Size result_idx = 0; result_idx < last_result_idx; ++result_idx)
          {
            OPENMS_LOG_DEBUG << "score: " << partial_results[result_idx].getMatchingScore() << " " << partial_results[result_idx].getMatchingSpectrumIndex() << endl;
            matching_results.push_back(partial_results[result_idx]);
          }
        }

        if (report_mode_ == "best")
        {
          if (!partial_results.empty())
          {
            matching_results.push_back(partial_results[0]);
          }
        }

      } // end precursor loop
    } // end spectra loop

    // write final results to MzTab
    exportMzTab_(matching_results, mztab_out);
  }


  /// protected methods

  void MetaboliteSpectralMatching::updateMembers_()
  {
    precursor_mz_error_ = (double)param_.getValue("prec_mass_error_value");
    fragment_mz_error_ = (double)param_.getValue("frag_mass_error_value");
    ion_mode_ = param_.getValue("ionization_mode").toString();

    mz_error_unit_ = param_.getValue("mass_error_unit").toString();
    report_mode_ = param_.getValue("report_mode").toString();
    merge_spectra_ = (bool)param_.getValue("merge_spectra").toBool();
  }


  /// private methods

  void MetaboliteSpectralMatching::exportMzTab_(const vector<SpectralMatch>& overall_results, MzTab& mztab_out)
  {
    // iterate the overall results table
    MzTabSmallMoleculeSectionRows all_sm_rows;

    for (Size id_idx = 0; id_idx < overall_results.size(); ++id_idx)
    {
      const SpectralMatch& current_id(overall_results[id_idx]);

      MzTabSmallMoleculeSectionRow mztab_row_record;

      // set the identifier field
      String hid_temp = current_id.getPrimaryIdentifier();
      MzTabString prim_id;
      prim_id.set(hid_temp);
      vector<MzTabString> id_dummy;
      id_dummy.push_back(prim_id);
      MzTabStringList string_dummy_list;
      string_dummy_list.set(id_dummy);

      mztab_row_record.identifier = string_dummy_list;

      // set the chemical formula field
      MzTabString chem_form;
      String form_temp = current_id.getSumFormula();
      chem_form.set(form_temp);

      mztab_row_record.chemical_formula = chem_form;

      // set the smiles field
      String smi_temp = current_id.getSMILESString();     // extract SMILES from struct mapping file
      MzTabString smi_string;
      smi_string.set(smi_temp);

      mztab_row_record.smiles = smi_string;

      // set the inchi_key field
      String inchi_temp = current_id.getInchiString();    // extract INCHIKEY from struct mapping file
      MzTabString inchi_key;
      inchi_key.set(inchi_temp);

      mztab_row_record.inchi_key = inchi_key;

      // set description field (we use it for the common name of the compound)
      String name_temp = current_id.getCommonName();
      MzTabString common_name;
      common_name.set(name_temp);

      mztab_row_record.description = common_name;

      // set mass_to_charge field (precursor mass here)
      double mz_temp = current_id.getFoundPrecursorMass();
      MzTabDouble mass_to_charge;
      mass_to_charge.set(mz_temp);

      mztab_row_record.exp_mass_to_charge = mass_to_charge;  //TODO: distinguish the experimental precursor mass and spectral library precursor mass (later should probably go into cv_opt_ column)

      // set charge field
      Int ch_temp = current_id.getFoundPrecursorCharge();
      MzTabInteger mcharge;
      mcharge.set(ch_temp);

      mztab_row_record.charge = mcharge;

      // set RT field
      double rt_temp = current_id.getObservedPrecursorRT();
      MzTabDouble rt_temp2;
      rt_temp2.set(rt_temp);
      vector<MzTabDouble> rt_temp3;
      rt_temp3.push_back(rt_temp2);
      MzTabDoubleList observed_rt;
      observed_rt.set(rt_temp3);

      mztab_row_record.retention_time = observed_rt;

      // set database field
      String dbname_temp = "MassBank";
      MzTabString dbname;
      dbname.set(dbname_temp);

      mztab_row_record.database = dbname;

      // set database_version field
      String dbver_temp = "Sep 27, 2013";
      MzTabString dbversion;
      dbversion.set(dbver_temp);

      mztab_row_record.database_version = dbversion;

      // set smallmolecule_abundance_sub
      // check if we deal with a feature or consensus feature
      vector<MzTabDouble> int_temp3;

      double int_temp(0.0);
      MzTabDouble int_temp2;
      int_temp2.set(int_temp);
      int_temp3.push_back(int_temp2);

      for (Size i = 0; i != int_temp3.size(); ++i)
      {
        mztab_row_record.smallmolecule_abundance_study_variable[i + 1] = int_temp3[i];
      }

      // set smallmolecule_abundance_stdev_sub; not applicable for a single feature intensity, however must be filled. Otherwise, the mzTab export fails.
      double stdev_temp(0.0);
      MzTabDouble stdev_temp2;
      stdev_temp2.set(stdev_temp);
      vector<MzTabDouble> stdev_temp3;

      stdev_temp3.push_back(stdev_temp2);

      for (Size i = 0; i != stdev_temp3.size(); ++i)
      {
        mztab_row_record.smallmolecule_abundance_stdev_study_variable[i + 1] = stdev_temp3[i];
      }

      // set smallmolecule_abundance_std_error_sub; not applicable for a single feature intensity, however must be filled. Otherwise, the mzTab export fails.
      double stderr_temp(0.0);
      MzTabDouble stderr_temp2;
      stderr_temp2.set(stderr_temp);
      vector<MzTabDouble> stderr_temp3;

      stderr_temp3.push_back(stderr_temp2);

      for (Size i = 0; i != stderr_temp3.size(); ++i)
      {
        mztab_row_record.smallmolecule_abundance_std_error_study_variable[i + 1] = stderr_temp3[i];
      }

      // optional columns:
      vector<MzTabOptionalColumnEntry> optionals;

      // ppm error
      double error_ppm(((current_id.getFoundPrecursorMass() - current_id.getObservedPrecursorMass())/current_id.getFoundPrecursorMass())*1e6);
      error_ppm = floor(error_ppm*100)/100;

      MzTabString ppmerr;
      ppmerr.set(String(error_ppm));
      MzTabOptionalColumnEntry col0;
      col0.first = "opt_ppm_error";
      col0.second = ppmerr;
      optionals.push_back(col0);

      // set found adduct ion
      String addion_temp = current_id.getPrecursorAdduct();
      MzTabString addion;
      addion.set(addion_temp);
      MzTabOptionalColumnEntry col1;
      col1.first = "opt_adduct_ion";
      col1.second = addion;
      optionals.push_back(col1);

      // set isotope similarity score
      double sim_score_temp = current_id.getMatchingScore();
      stringstream read_in;
      read_in << sim_score_temp;
      String sim_score_temp2(read_in.str());
      MzTabString sim_score;
      sim_score.set(sim_score_temp2);
      MzTabOptionalColumnEntry col2;
      col2.first = "opt_match_score";
      col2.second = sim_score;
      optionals.push_back(col2);

      // set secondary ID (here HMDB id)
      String sec_id = current_id.getSecondaryIdentifier();
      MzTabString sec_id_str;
      sec_id_str.set(sec_id);
      MzTabOptionalColumnEntry col3;
      col3.first = "opt_sec_id";
      col3.second = sec_id_str;
      optionals.push_back(col3);

      // set source spectra index
      // TODO: this should use spectra_ref column
      String source_idx = String(current_id.getObservedSpectrumIndex());
      MzTabString source_idx_str;
      source_idx_str.set(source_idx);
      MzTabOptionalColumnEntry col4;
      col4.first = "opt_source_idx";
      col4.second = source_idx_str;
      optionals.push_back(col4);

      // set spectrum native ID
      String spec_native_id = current_id.getObservedSpectrumNativeID();
      MzTabString spec_native_id_str;
      spec_native_id_str.set(spec_native_id);
      MzTabOptionalColumnEntry col5;
      col5.first = "opt_spec_native_id";
      col5.second = spec_native_id_str;
      optionals.push_back(col5);

      mztab_row_record.opt_ = optionals;

      all_sm_rows.push_back(mztab_row_record);
    }
    mztab_out.setSmallMoleculeSectionRows(all_sm_rows);
  }

} // closing namespace OpenMS
