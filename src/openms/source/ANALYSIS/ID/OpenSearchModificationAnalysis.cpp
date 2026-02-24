// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/OpenSearchModificationAnalysis.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <algorithm>
#include <cmath>
#include <fstream>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include <boost/math/distributions/normal.hpp>

namespace OpenMS
{

  std::pair<OpenSearchModificationAnalysis::DeltaMassHistogram, OpenSearchModificationAnalysis::DeltaMassToChargeCount>
  OpenSearchModificationAnalysis::analyzeDeltaMassPatterns(const PeptideIdentificationList& peptide_ids, 
                                                          bool use_smoothing, 
                                                          bool /*debug*/) const
  {
    // Constants
    constexpr double deltamass_tolerance = 0.0005;

    // Lambda to round values to the specified tolerance
    auto roundToTolerance = [](double value) {
      return std::round(value / deltamass_tolerance) * deltamass_tolerance;
    };

    // Data structures to store histogram and charge states
    DeltaMassHistogram histogram(FuzzyDoubleComparator(1e-9));
    DeltaMassToChargeCount charge_counts(FuzzyDoubleComparator(1e-9));
    std::unordered_map<double, std::unordered_set<int>> charge_states;

    // Process each peptide identification
    for (const auto& peptide_id : peptide_ids)
    {
      const auto& hits = peptide_id.getHits();
      for (const auto& hit : hits)
      {
        // Retrieve delta mass and charge
        if (!hit.metaValueExists("DeltaMass"))
          continue;
          
        double delta_mass = hit.getMetaValue("DeltaMass");
        int charge = hit.getCharge();

        // Ignore delta masses close to zero
        if (std::abs(delta_mass) <= DELTA_MASS_ZERO_THRESHOLD_)
          continue;

        // Round delta mass to bin similar values
        double rounded_mass = roundToTolerance(delta_mass);

        // Update histogram count
        histogram[rounded_mass] += 1.0;

        // Update unique charge count
        if (charge_states[rounded_mass].insert(charge).second)
        {
          charge_counts[rounded_mass] += 1;
        }
      }
    }

    // Prepare results
    std::pair<DeltaMassHistogram, DeltaMassToChargeCount> results{histogram, charge_counts};

    // Apply smoothing if requested
    if (use_smoothing)
    {
      DeltaMassHistogram smoothed_hist = smoothDeltaMassHistogram_(histogram, 0.0001);
      DeltaMassHistogram hist_maxima = findPeaksInHistogram_(smoothed_hist, 0.0, 3.0);

      // Update charge counts for the smoothed maxima
      DeltaMassToChargeCount smoothed_charge_counts(FuzzyDoubleComparator(1e-9));
      for (const auto& [mass, _] : hist_maxima)
      {
        smoothed_charge_counts[mass] = charge_counts[mass];
      }

      // Update results with smoothed data
      results = {hist_maxima, smoothed_charge_counts};
    }

    return results;
  }

  std::vector<OpenSearchModificationAnalysis::ModificationSummary>
  OpenSearchModificationAnalysis::mapDeltaMassesToModifications(const DeltaMassHistogram& delta_mass_histogram,
                                                               const DeltaMassToChargeCount& charge_histogram,
                                                               PeptideIdentificationList& peptide_ids,
                                                               double precursor_mass_tolerance,
                                                               bool precursor_mass_tolerance_unit_ppm,
                                                               const String& output_file) const
  {
    std::map<double, String, FuzzyDoubleComparator> mass_to_modification(FuzzyDoubleComparator(1e-9));
    std::map<String, ModificationPattern> modifications;
    std::map<double, String> histogram_found;

    // Load modifications from the database
    std::vector<String> modification_names;
    ModificationsDB* mod_db = ModificationsDB::getInstance();
    mod_db->getAllSearchModifications(modification_names);
    
    for (const String& mod_name : modification_names)
    {
      const ResidueModification* residue = mod_db->getModification(mod_name);
      String full_name = residue->getFullName();
      double diff_mono_mass = residue->getDiffMonoMass();
      
      if (full_name.find("substitution") == std::string::npos)
      {
        mass_to_modification[diff_mono_mass] = full_name;
      }
    }

    // Generate combinations of modifications
    std::map<double, String, FuzzyDoubleComparator> combo_modifications(FuzzyDoubleComparator(1e-9));
    for (auto it1 = mass_to_modification.begin(); it1 != mass_to_modification.end(); ++it1)
    {
      for (auto it2 = it1; it2 != mass_to_modification.end(); ++it2)
      {
        combo_modifications[it1->first + it2->first] = it1->second + "++" + it2->second;
      }
    }

    // Helper function to add or update modifications
    auto addOrUpdateModification = [&](const String& mod_name, double mass, double count, int num_charges)
    {
      if (modifications.find(mod_name) == modifications.end())
      {
        ModificationPattern pattern;
        pattern.masses.push_back(mass);
        pattern.count = count;
        pattern.num_charge_states = num_charges;
        modifications[mod_name] = pattern;
      }
      else
      {
        modifications[mod_name].count += count;
        modifications[mod_name].num_charge_states = std::max(num_charges, modifications[mod_name].num_charge_states);
      }
    };

    // Compute effective tolerance for matching delta masses to known modifications.
    // Cap at MAX_MOD_MAPPING_TOL_ to prevent overly broad matching in open search
    // (where precursor tolerance can be e.g. 500 Da).
    // For ppm mode, use the cap directly (ppm→Da requires a per-mass reference
    // which is not meaningful for a single tolerance value).
    const double effective_tol = precursor_mass_tolerance_unit_ppm
      ? MAX_MOD_MAPPING_TOL_
      : std::min(precursor_mass_tolerance, MAX_MOD_MAPPING_TOL_);
    const double epsilon = 1e-8;

    // Map delta masses to modifications
    for (const auto& hist_entry : delta_mass_histogram)
    {
      double cluster_mass = hist_entry.first;
      double count = hist_entry.second;
      double lower_bound = cluster_mass - effective_tol;
      double upper_bound = cluster_mass + effective_tol;

      // Search for modifications within bounds
      bool mapping_found = false;
      String mod_name;
      double mod_mass = 0.0;

      // Search in single modifications
      auto it_lower = mass_to_modification.lower_bound(lower_bound - epsilon);
      bool found_lower = false;
      if (it_lower != mass_to_modification.end() &&
          std::abs(it_lower->first - cluster_mass) <= effective_tol)
      {
        found_lower = true;
      }

      auto it_upper = mass_to_modification.upper_bound(upper_bound + epsilon);
      bool found_upper = false;
      if (it_upper != mass_to_modification.begin())
      {
        --it_upper;
        if (std::abs(it_upper->first - cluster_mass) <= effective_tol)
        {
          found_upper = true;
        }
      }

      // Compare results from lower_bound and upper_bound
      if (found_lower && found_upper)
      {
        if (it_lower->first == it_upper->first && it_lower->second == it_upper->second)
        {
          mod_name = it_lower->second;
          mod_mass = it_lower->first;
          histogram_found[mod_mass] = mod_name;
          mapping_found = true;
        }
        else
        {
          mod_name = it_lower->second + "//" + it_upper->second;
          mod_mass = cluster_mass;
          histogram_found[it_lower->first] = it_lower->second;
          histogram_found[it_upper->first] = it_upper->second;
          mapping_found = true;
        }
      }
      else
      {
        // Check if modification can be explained by known modifications
        for (const auto& hit : histogram_found)
        {
          if (std::abs(hit.first - cluster_mass) < effective_tol)
          {
            addOrUpdateModification(hit.second, hit.first, count, charge_histogram.at(cluster_mass));
            mapping_found = true;
            break;
          }
          // Check if modification can be explained by a +1 isotope variant
          else if (std::abs((hit.first + 1.0) - cluster_mass) < effective_tol)
          {
            String temp_mod_name = hit.second + "+1Da";
            addOrUpdateModification(temp_mod_name, hit.first + 1.0, count, charge_histogram.at(cluster_mass));
            histogram_found[hit.first + 1.0] = temp_mod_name;
            mapping_found = true;
            break;
          }
        }

        // Search in combination modifications
        if (!mapping_found)
        {
          auto it = combo_modifications.lower_bound(cluster_mass - epsilon);
          if (it != combo_modifications.end() && 
              std::abs(it->first - cluster_mass) <= effective_tol)
          {
            mod_name = it->second;
            mod_mass = it->first;
            mapping_found = true;
          }
        }
      }

      if (mapping_found)
      {
        // Skip if the matched modification mass is near zero (unmodified)
        if (std::abs(mod_mass) < DELTA_MASS_ZERO_THRESHOLD_)
          continue;
        addOrUpdateModification(mod_name, mod_mass, count, charge_histogram.at(cluster_mass));
      }
      else
      {
        // Unknown modification (cluster_mass is already filtered for near-zero in analyzeDeltaMassPatterns)
        String unknown_mod_name = "Unknown" + String(std::round(cluster_mass));
        addOrUpdateModification(unknown_mod_name, cluster_mass, count, charge_histogram.at(cluster_mass));
      }
    }

    // Collect all modification data into a vector
    std::vector<ModificationSummary> modification_summaries;
    
    for (const auto& mod_pair : modifications)
    {
      ModificationSummary summary;
      summary.count = static_cast<int>(std::round(mod_pair.second.count));
      summary.name = mod_pair.first;
      summary.num_charge_states = mod_pair.second.num_charge_states;
      summary.masses = mod_pair.second.masses;
      
      modification_summaries.push_back(summary);
    }

    // Sort modifications by (num_charge_states + count) in descending order
    std::sort(modification_summaries.begin(), modification_summaries.end(),
              [](const ModificationSummary& a, const ModificationSummary& b)
              {
                return (a.num_charge_states + a.count) > (b.num_charge_states + b.count);
              });

    // Add modifications to peptide identifications
    for (auto& peptide_id : peptide_ids)
    {
      auto& hits = peptide_id.getHits();
      for (auto& hit : hits)
      {
        if (!hit.metaValueExists("DeltaMass"))
          continue;
          
        double delta_mass = hit.getMetaValue("DeltaMass");
        String ptm = "";

        // Check if too close to zero
        if (std::abs(delta_mass) < DELTA_MASS_ZERO_THRESHOLD_)
        {
          hit.setMetaValue("PTM", ptm);
          continue;
        }

        bool found = false;
        // Check with error tolerance if already present in histogram
        for (const auto& entry : histogram_found)
        {
          if (std::abs(delta_mass - entry.first) < effective_tol)
          {
            ptm = entry.second;
            found = true;
            break;
          }
        }
        
        // Otherwise assign unknown
        if (!found)
        {
          ptm = "Unknown" + String(delta_mass);
        }
        
        hit.setMetaValue("PTM", ptm);
      }
    }

    // Write modification summary table if output file is specified
    if (!output_file.empty())
    {
      writeModificationSummary_(modification_summaries, output_file);
    }

    return modification_summaries;
  }

  std::vector<OpenSearchModificationAnalysis::ModificationSummary>
  OpenSearchModificationAnalysis::analyzeModifications(PeptideIdentificationList& peptide_ids,
                                                      double precursor_mass_tolerance,
                                                      bool precursor_mass_tolerance_unit_ppm,
                                                      bool use_smoothing,
                                                      const String& output_file) const
  {
    // Analyze delta mass patterns
    auto [histogram, charge_counts] = analyzeDeltaMassPatterns(peptide_ids, use_smoothing, false);
    
    // Map to modifications and annotate peptides
    return mapDeltaMassesToModifications(histogram, charge_counts, peptide_ids,
                                       precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm,
                                       output_file);
  }

  // Private helper functions

  double OpenSearchModificationAnalysis::gaussian_(double x, double sigma)
  {
    boost::math::normal_distribution<> normal_dist(0.0, sigma);
    return boost::math::pdf(normal_dist, x);
  }

  OpenSearchModificationAnalysis::DeltaMassHistogram 
  OpenSearchModificationAnalysis::smoothDeltaMassHistogram_(const DeltaMassHistogram& histogram, double sigma)
  {
    if (histogram.size() < 3)
    {
      return histogram; // Not enough data points for smoothing
    }

    DeltaMassHistogram smoothed_histogram(FuzzyDoubleComparator(1e-9));

    // Extract delta masses and counts into vectors for efficient access
    std::vector<double> deltas;
    std::vector<double> counts;
    deltas.reserve(histogram.size());
    counts.reserve(histogram.size());

    for (const auto& [delta, count] : histogram)
    {
      deltas.push_back(delta);
      counts.push_back(count);
    }

    const size_t n = deltas.size();
    std::vector<double> smoothed_counts(n, 0.0);

    // Perform Gaussian smoothing
    for (size_t i = 0; i < n; ++i)
    {
      double weight_sum = 0.0;

      for (size_t j = 0; j < n; ++j)
      {
        double mz_diff = deltas[i] - deltas[j];

        // Ignore points beyond 3 standard deviations
        if (std::abs(mz_diff) > 3.0 * sigma)
          continue;

        double weight = gaussian_(mz_diff, sigma);
        smoothed_counts[i] += weight * counts[j];
        weight_sum += weight;
      }

      if (weight_sum != 0.0)
      {
        smoothed_counts[i] /= weight_sum;
      }
    }

    // Populate the smoothed histogram
    for (size_t i = 0; i < n; ++i)
    {
      smoothed_histogram[deltas[i]] = smoothed_counts[i];
    }

    return smoothed_histogram;
  }

  OpenSearchModificationAnalysis::DeltaMassHistogram 
  OpenSearchModificationAnalysis::findPeaksInHistogram_(const DeltaMassHistogram& histogram, 
                                                       double count_threshold, 
                                                       double snr)
  {
    if (histogram.size() < 3)
    {
      return histogram; // Not enough data points to find peaks
    }

    DeltaMassHistogram peaks(FuzzyDoubleComparator(1e-9));

    // Extract counts to compute noise level (median count)
    std::vector<double> counts;
    counts.reserve(histogram.size());

    for (const auto& [_, count] : histogram)
    {
      counts.push_back(count);
    }

    // Calculate median as noise level
    std::nth_element(counts.begin(), counts.begin() + counts.size() / 2, counts.end());
    double noise_level = counts[counts.size() / 2];

    // Convert histogram to vector for indexed access
    std::vector<std::pair<double, double>> hist_vector(histogram.begin(), histogram.end());

    // Check each point except the first and last for local maxima
    for (size_t i = 1; i < hist_vector.size() - 1; ++i)
    {
      double prev_count = hist_vector[i - 1].second;
      double curr_count = hist_vector[i].second;
      double next_count = hist_vector[i + 1].second;

      // Check if current point is a local maximum
      if (curr_count >= prev_count && curr_count >= next_count &&
          curr_count > count_threshold &&
          curr_count / noise_level > snr)
      {
        peaks[hist_vector[i].first] = curr_count;
      }
    }

    return peaks;
  }

  void OpenSearchModificationAnalysis::writeModificationSummary_(const std::vector<ModificationSummary>& modifications,
                                                                const String& output_file) const
  {
    // Remove 'idxml' extension and add '_OutputTable.tsv'
    String output_table = output_file;
    if (output_table.hasSuffix(".idXML"))
    {
      output_table = output_table.substr(0, output_table.size() - 6) + "_OutputTable.tsv";
    }
    else if (output_table.hasSuffix(".idxml"))
    {
      output_table = output_table.substr(0, output_table.size() - 6) + "_OutputTable.tsv";
    }
    else
    {
      output_table += "_OutputTable.tsv";
    }

    std::ofstream output_stream(output_table);
    if (!output_stream.is_open())
    {
      OPENMS_LOG_ERROR << "Error opening file: " << output_table << std::endl;
      return;
    }

    output_stream << "Name\tMass\tModified Peptides (incl. charge variants)\tModified Peptides\n";

    for (const auto& mod_data : modifications)
    {
      output_stream << mod_data.name << '\t';

      // Output mass or masses
      if (mod_data.masses.size() < 2)
      {
        output_stream << mod_data.masses.at(0) << '\t';
      }
      else
      {
        output_stream << mod_data.masses.at(0) << "/" << mod_data.masses.at(1) << '\t';
      }

      // Output counts
      output_stream << mod_data.num_charge_states + mod_data.count << '\t'
                   << mod_data.count << '\n';
    }

    output_stream.close();
  }

  // New methods for structured statistics tables

  OpenSearchModificationAnalysis::OpenSearchAnalysisResult
  OpenSearchModificationAnalysis::analyzeModificationsWithStatistics(PeptideIdentificationList& peptide_ids,
                                                                     double precursor_mass_tolerance,
                                                                     bool precursor_mass_tolerance_unit_ppm,
                                                                     bool use_smoothing,
                                                                     const String& output_file) const
  {
    OpenSearchAnalysisResult result;

    // Analyze delta mass patterns
    auto [histogram, charge_counts] = analyzeDeltaMassPatterns(peptide_ids, use_smoothing, false);

    // Map to modifications and annotate peptides (also fills legacy summaries)
    result.summaries = mapDeltaMassesToModifications(histogram, charge_counts, peptide_ids,
                                                     precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm,
                                                     output_file);

    // Generate structured statistics tables
    result.delta_mass_stats = generateDeltaMassStatistics(histogram, charge_counts, peptide_ids,
                                                          precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);

    result.ptm_stats = generatePTMStatistics(peptide_ids, precursor_mass_tolerance, precursor_mass_tolerance_unit_ppm);

    // Write statistics tables if output file is specified
    if (!output_file.empty())
    {
      String base_name = output_file;
      if (base_name.hasSuffix(".idXML") || base_name.hasSuffix(".idxml"))
      {
        base_name = base_name.substr(0, base_name.size() - 6);
      }
      writeDeltaMassStatistics(result.delta_mass_stats, base_name + "_DeltaMassStats.tsv");
      writePTMStatistics(result.ptm_stats, base_name + "_PTMStats.tsv");
    }

    return result;
  }

  OpenSearchModificationAnalysis::DeltaMassStatistics
  OpenSearchModificationAnalysis::generateDeltaMassStatistics(const DeltaMassHistogram& histogram,
                                                              const DeltaMassToChargeCount& charge_histogram,
                                                              const PeptideIdentificationList& peptide_ids,
                                                              double precursor_mass_tolerance,
                                                              bool precursor_mass_tolerance_unit_ppm) const
  {
    DeltaMassStatistics stats;
    constexpr double min_mass_for_ppm = 0.1; // Minimum mass to avoid division issues with ppm

    // Effective tolerance for modification mass matching (capped for open search)
    const double mod_tol = precursor_mass_tolerance_unit_ppm
      ? MAX_MOD_MAPPING_TOL_
      : std::min(precursor_mass_tolerance, MAX_MOD_MAPPING_TOL_);

    // Build modification lookup
    auto mod_lookup = buildModificationMassLookup_();

    // Count total PSMs
    for (const auto& peptide_id : peptide_ids)
    {
      for (const auto& hit : peptide_id.getHits())
      {
        stats.total_psms++;
        if (hit.metaValueExists("DeltaMass"))
        {
          double delta_mass = hit.getMetaValue("DeltaMass");
          if (std::abs(delta_mass) <= DELTA_MASS_ZERO_THRESHOLD_)
          {
            stats.unmodified_psms++;
          }
          else
          {
            stats.modified_psms++;
          }
        }
      }
    }

    // Collect delta masses for median calculation
    std::vector<double> all_delta_masses;
    double sum_delta_mass = 0.0;

    // Convert histogram to entries
    for (const auto& [delta_mass, count] : histogram)
    {
      DeltaMassEntry entry;
      entry.delta_mass = delta_mass;
      entry.count = static_cast<int>(count);

      // Get charge state count
      auto charge_it = charge_histogram.find(delta_mass);
      if (charge_it != charge_histogram.end())
      {
        entry.num_charge_states = charge_it->second;
      }

      // Calculate percentage
      if (stats.modified_psms > 0)
      {
        entry.percentage = (static_cast<double>(entry.count) / stats.modified_psms) * 100.0;
      }

      // Compute tolerance in Da based on unit
      double tolerance_da = precursor_mass_tolerance;
      if (precursor_mass_tolerance_unit_ppm)
      {
        // For ppm, tolerance depends on the mass being compared
        // Use the delta_mass as reference, with a minimum to avoid issues near zero
        double reference_mass = std::max(min_mass_for_ppm, std::abs(delta_mass));
        tolerance_da = reference_mass * precursor_mass_tolerance * 1e-6;
      }

      // Count unique peptides
      entry.unique_peptides = countUniquePeptides_(peptide_ids, delta_mass, tolerance_da);

      // Try to map to known modification
      for (const auto& [mod_mass, mod_name] : mod_lookup)
      {
        if (std::abs(mod_mass - delta_mass) <= mod_tol)
        {
          entry.mapped_modification = mod_name;
          entry.is_known_modification = true;
          break;
        }
      }

      // Collect for statistics
      for (int i = 0; i < entry.count; ++i)
      {
        all_delta_masses.push_back(delta_mass);
        sum_delta_mass += delta_mass;
      }

      stats.entries.push_back(entry);
    }

    // Sort entries by count (descending)
    std::sort(stats.entries.begin(), stats.entries.end(),
              [](const DeltaMassEntry& a, const DeltaMassEntry& b) {
                return a.count > b.count;
              });

    // Calculate mean and median
    if (!all_delta_masses.empty())
    {
      stats.mean_delta_mass = sum_delta_mass / all_delta_masses.size();

      std::nth_element(all_delta_masses.begin(),
                       all_delta_masses.begin() + all_delta_masses.size() / 2,
                       all_delta_masses.end());
      stats.median_delta_mass = all_delta_masses[all_delta_masses.size() / 2];
    }

    return stats;
  }

  OpenSearchModificationAnalysis::PTMStatistics
  OpenSearchModificationAnalysis::generatePTMStatistics(const PeptideIdentificationList& peptide_ids,
                                                        double /* precursor_mass_tolerance */,
                                                        bool /* precursor_mass_tolerance_unit_ppm */) const
  {
    // Note: tolerance parameters are part of the API for consistency but not used here.
    // PTM annotations are already assigned by mapDeltaMassesToModifications() which applies the tolerance.
    // This function only aggregates statistics from existing PTM meta values.

    PTMStatistics stats;

    // Map to collect PTM data
    std::map<String, PTMEntry> ptm_map;
    std::map<String, std::unordered_set<std::string>> ptm_unique_peptides;
    std::map<String, std::unordered_set<int>> ptm_unique_charges;

    // Build modification lookup for theoretical masses
    auto mod_lookup = buildModificationMassLookup_();
    std::map<String, double> name_to_mass;
    for (const auto& [mass, name] : mod_lookup)
    {
      name_to_mass[name] = mass;
    }

    // Process all peptide identifications
    for (const auto& peptide_id : peptide_ids)
    {
      for (const auto& hit : peptide_id.getHits())
      {
        if (!hit.metaValueExists("PTM"))
          continue;

        String ptm_name = hit.getMetaValue("PTM");
        if (ptm_name.empty())
        {
          continue;
        }

        // Check if unknown modification
        if (ptm_name.hasPrefix("Unknown"))
        {
          stats.unknown_modification_psms++;
          continue;
        }

        stats.total_modified_psms++;

        // Get or create PTM entry
        if (ptm_map.find(ptm_name) == ptm_map.end())
        {
          PTMEntry entry;
          entry.name = ptm_name;

          // Get theoretical mass
          if (name_to_mass.find(ptm_name) != name_to_mass.end())
          {
            entry.theoretical_mass = name_to_mass[ptm_name];
          }

          // Get target residues
          entry.target_residues = getTargetResidues_(ptm_name);

          ptm_map[ptm_name] = entry;
        }

        PTMEntry& entry = ptm_map[ptm_name];
        entry.count++;

        // Track unique peptides
        std::string seq_str = hit.getSequence().toString();
        if (ptm_unique_peptides[ptm_name].insert(seq_str).second)
        {
          entry.unique_peptides++;
        }

        // Track unique charge states per PTM
        ptm_unique_charges[ptm_name].insert(hit.getCharge());

        // Analyze observed mass and residue frequency if available
        if (hit.metaValueExists("DeltaMass"))
        {
          double obs_mass = hit.getMetaValue("DeltaMass");
          // Running average for observed mass
          entry.observed_mass = ((entry.observed_mass * (entry.count - 1)) + obs_mass) / entry.count;

          // Count all residues in this peptide for residue frequency analysis
          const AASequence& seq = hit.getSequence();
          for (Size i = 0; i < seq.size(); ++i)
          {
            char residue = seq[i].getOneLetterCode()[0];
            entry.residue_counts[residue]++;
          }
        }
      }
    }

    // Convert map to vector and calculate derived statistics
    for (auto& [name, entry] : ptm_map)
    {
      // Set charge state count from tracked unique charges
      entry.num_charge_states = static_cast<int>(ptm_unique_charges[name].size());

      // Calculate mass deviation
      if (entry.theoretical_mass != 0.0)
      {
        entry.mass_deviation = entry.observed_mass - entry.theoretical_mass;
      }

      // Calculate percentage
      if (stats.total_modified_psms > 0)
      {
        entry.percentage = (static_cast<double>(entry.count) / stats.total_modified_psms) * 100.0;
      }

      stats.entries.push_back(entry);
    }

    // Sort by count (descending)
    std::sort(stats.entries.begin(), stats.entries.end(),
              [](const PTMEntry& a, const PTMEntry& b) {
                return a.count > b.count;
              });

    stats.num_unique_modifications = static_cast<int>(stats.entries.size());

    return stats;
  }

  std::map<char, int>
  OpenSearchModificationAnalysis::analyzeResidueFrequency(const PeptideIdentificationList& peptide_ids,
                                                          double delta_mass,
                                                          double tolerance) const
  {
    std::map<char, int> residue_counts;

    for (const auto& peptide_id : peptide_ids)
    {
      for (const auto& hit : peptide_id.getHits())
      {
        if (!hit.metaValueExists("DeltaMass"))
          continue;

        double hit_delta_mass = hit.getMetaValue("DeltaMass");
        if (std::abs(hit_delta_mass - delta_mass) > tolerance)
          continue;

        const AASequence& seq = hit.getSequence();
        for (Size i = 0; i < seq.size(); ++i)
        {
          char residue = seq[i].getOneLetterCode()[0];
          residue_counts[residue]++;
        }
      }
    }

    return residue_counts;
  }

  void OpenSearchModificationAnalysis::writeDeltaMassStatistics(const DeltaMassStatistics& stats,
                                                                const String& output_file) const
  {
    std::ofstream output_stream(output_file);
    if (!output_stream.is_open())
    {
      OPENMS_LOG_ERROR << "Error opening file: " << output_file << std::endl;
      return;
    }

    // Write header
    output_stream << "# Delta Mass Statistics Table\n";
    output_stream << "# Total PSMs: " << stats.total_psms << "\n";
    output_stream << "# Modified PSMs: " << stats.modified_psms << "\n";
    output_stream << "# Unmodified PSMs: " << stats.unmodified_psms << "\n";
    output_stream << "# Mean Delta Mass: " << stats.mean_delta_mass << "\n";
    output_stream << "# Median Delta Mass: " << stats.median_delta_mass << "\n";
    output_stream << "#\n";
    output_stream << "DeltaMass\tCount\tUniquePeptides\tChargeStates\tPercentage\tMappedModification\tIsKnown\n";

    for (const auto& entry : stats.entries)
    {
      output_stream << entry.delta_mass << '\t'
                    << entry.count << '\t'
                    << entry.unique_peptides << '\t'
                    << entry.num_charge_states << '\t'
                    << entry.percentage << '\t'
                    << entry.mapped_modification << '\t'
                    << (entry.is_known_modification ? "true" : "false") << '\n';
    }

    output_stream.close();
    OPENMS_LOG_INFO << "Delta mass statistics written to: " << output_file << std::endl;
  }

  void OpenSearchModificationAnalysis::writePTMStatistics(const PTMStatistics& stats,
                                                          const String& output_file) const
  {
    std::ofstream output_stream(output_file);
    if (!output_stream.is_open())
    {
      OPENMS_LOG_ERROR << "Error opening file: " << output_file << std::endl;
      return;
    }

    // Write header
    output_stream << "# PTM Statistics Table\n";
    output_stream << "# Total Modified PSMs: " << stats.total_modified_psms << "\n";
    output_stream << "# Unknown Modification PSMs: " << stats.unknown_modification_psms << "\n";
    output_stream << "# Unique Modifications: " << stats.num_unique_modifications << "\n";
    output_stream << "#\n";
    output_stream << "Name\tTheoreticalMass\tObservedMass\tMassDeviation\tCount\tUniquePeptides\tPercentage\tTargetResidues\tResidueFrequency\n";

    for (const auto& entry : stats.entries)
    {
      output_stream << entry.name << '\t'
                    << entry.theoretical_mass << '\t'
                    << entry.observed_mass << '\t'
                    << entry.mass_deviation << '\t'
                    << entry.count << '\t'
                    << entry.unique_peptides << '\t'
                    << entry.percentage << '\t'
                    << entry.target_residues << '\t';

      // Format residue frequency as "A:10,C:5,..."
      bool first = true;
      for (const auto& [residue, count] : entry.residue_counts)
      {
        if (!first) output_stream << ',';
        output_stream << residue << ':' << count;
        first = false;
      }
      output_stream << '\n';
    }

    output_stream.close();
    OPENMS_LOG_INFO << "PTM statistics written to: " << output_file << std::endl;
  }

  std::map<double, String, OpenSearchModificationAnalysis::FuzzyDoubleComparator>
  OpenSearchModificationAnalysis::buildModificationMassLookup_() const
  {
    std::map<double, String, FuzzyDoubleComparator> mass_to_mod; // uses default epsilon (1e-9)

    std::vector<String> modification_names;
    ModificationsDB* mod_db = ModificationsDB::getInstance();
    mod_db->getAllSearchModifications(modification_names);

    for (const String& mod_name : modification_names)
    {
      const ResidueModification* residue = mod_db->getModification(mod_name);
      String full_name = residue->getFullName();
      double diff_mono_mass = residue->getDiffMonoMass();

      // Skip substitutions
      if (full_name.find("substitution") == std::string::npos)
      {
        mass_to_mod[diff_mono_mass] = full_name;
      }
    }

    return mass_to_mod;
  }

  String OpenSearchModificationAnalysis::getTargetResidues_(const String& mod_name) const
  {
    // Split compound names (e.g. "Oxidation//Deamidated" or "Oxidation+1Da") and resolve each part
    std::vector<String> parts;
    if (mod_name.find("//") != std::string::npos)
    {
      mod_name.split("//", parts);
    }
    else
    {
      // Strip isotope suffixes like "+1Da", "+2Da" etc.
      String base_name = mod_name;
      auto pos = mod_name.find("+");
      if (pos != std::string::npos && pos > 0)
      {
        String suffix = mod_name.substr(pos);
        if (suffix.hasSuffix("Da"))
        {
          base_name = mod_name.substr(0, pos);
        }
      }
      parts.push_back(base_name);
    }

    std::vector<String> residues;
    ModificationsDB* mod_db = ModificationsDB::getInstance();

    for (const auto& part : parts)
    {
      String trimmed = part;
      trimmed.trim();
      if (trimmed.empty()) continue;

      try
      {
        const ResidueModification* mod = mod_db->getModification(trimmed);
        if (mod != nullptr)
        {
          String result;
          String origin = mod->getOrigin();
          if (!origin.empty() && origin != "X")
          {
            result = origin;
          }

          auto term_spec = mod->getTermSpecificity();
          if (term_spec == ResidueModification::N_TERM)
          {
            result = "N-term" + (result.empty() ? "" : "(" + result + ")");
          }
          else if (term_spec == ResidueModification::C_TERM)
          {
            result = "C-term" + (result.empty() ? "" : "(" + result + ")");
          }
          else if (term_spec == ResidueModification::PROTEIN_N_TERM)
          {
            result = "Protein N-term" + (result.empty() ? "" : "(" + result + ")");
          }
          else if (term_spec == ResidueModification::PROTEIN_C_TERM)
          {
            result = "Protein C-term" + (result.empty() ? "" : "(" + result + ")");
          }

          if (!result.empty())
          {
            residues.push_back(result);
          }
        }
      }
      catch (const Exception::BaseException&)
      {
        OPENMS_LOG_WARN << "Could not resolve target residues for modification: " << trimmed
                        << " (from compound name: " << mod_name << ")" << std::endl;
      }
    }

    // Join results from multiple parts
    String result;
    for (size_t i = 0; i < residues.size(); ++i)
    {
      if (i > 0) result += ",";
      result += residues[i];
    }
    return result;
  }

  int OpenSearchModificationAnalysis::countUniquePeptides_(const PeptideIdentificationList& peptide_ids,
                                                           double delta_mass,
                                                           double tolerance) const
  {
    std::unordered_set<std::string> unique_sequences;

    for (const auto& peptide_id : peptide_ids)
    {
      for (const auto& hit : peptide_id.getHits())
      {
        if (!hit.metaValueExists("DeltaMass"))
          continue;

        double hit_delta_mass = hit.getMetaValue("DeltaMass");
        if (std::abs(hit_delta_mass - delta_mass) <= tolerance)
        {
          unique_sequences.insert(hit.getSequence().toString());
        }
      }
    }

    return static_cast<int>(unique_sequences.size());
  }

} // namespace OpenMS
