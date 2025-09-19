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

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

namespace OpenMS
{

  std::pair<OpenSearchModificationAnalysis::DeltaMassHistogram, OpenSearchModificationAnalysis::DeltaMassToChargeCount>
  OpenSearchModificationAnalysis::analyzeDeltaMassPatterns(const PeptideIdentificationList& peptide_ids, 
                                                          bool use_smoothing, 
                                                          bool /*debug*/) const
  {
    // Constants
    constexpr double deltamass_tolerance = 0.0005;
    constexpr double delta_mass_zero_threshold = 0.05;

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
        if (std::abs(delta_mass) <= delta_mass_zero_threshold)
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

    // Map delta masses to modifications
    for (const auto& hist_entry : delta_mass_histogram)
    {
      double cluster_mass = hist_entry.first;
      double count = hist_entry.second;
      
      double lower_bound, upper_bound;
      const double epsilon = 1e-8;

      if (precursor_mass_tolerance_unit_ppm)
      {
        double tolerance = cluster_mass * precursor_mass_tolerance * 1e-6;
        lower_bound = cluster_mass - tolerance;
        upper_bound = cluster_mass + tolerance;
      }
      else
      {
        lower_bound = cluster_mass - precursor_mass_tolerance;
        upper_bound = cluster_mass + precursor_mass_tolerance;
      }

      // Search for modifications within bounds
      bool mapping_found = false;
      String mod_name;
      double mod_mass = 0.0;

      // Search in single modifications
      auto it_lower = mass_to_modification.lower_bound(lower_bound - epsilon);
      bool found_lower = false;
      if (it_lower != mass_to_modification.end() && 
          std::abs(it_lower->first - cluster_mass) <= precursor_mass_tolerance)
      {
        found_lower = true;
      }

      auto it_upper = mass_to_modification.upper_bound(upper_bound + epsilon);
      bool found_upper = false;
      if (it_upper != mass_to_modification.begin())
      {
        --it_upper;
        if (std::abs(it_upper->first - cluster_mass) <= precursor_mass_tolerance)
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
          if (std::abs(hit.first - cluster_mass) < precursor_mass_tolerance)
          {
            addOrUpdateModification(hit.second, hit.first, count, charge_histogram.at(cluster_mass));
            mapping_found = true;
            break;
          }
          // Check if modification can be explained by a +1 isotope variant
          else if (std::abs((hit.first + 1.0) - cluster_mass) < precursor_mass_tolerance)
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
              std::abs(it->first - cluster_mass) <= precursor_mass_tolerance / 10.0)
          {
            mod_name = it->second;
            mod_mass = it->first;
            mapping_found = true;
          }
        }
      }

      if (std::abs(mod_mass) < precursor_mass_tolerance) 
        continue; // Skip if closest mod_mass is too close to 0

      if (mapping_found)
      {
        addOrUpdateModification(mod_name, mod_mass, count, charge_histogram.at(cluster_mass));
      }
      else
      {
        // Unknown modification
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
        if (std::abs(delta_mass) < 0.05)
        {
          hit.setMetaValue("PTM", ptm);
          continue;
        }

        bool found = false;
        // Check with error tolerance if already present in histogram
        for (const auto& entry : histogram_found)
        {
          if (std::abs(delta_mass - entry.first) < precursor_mass_tolerance)
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
    return std::exp(-(x * x) / (2 * sigma * sigma)) / (sigma * std::sqrt(2 * M_PI));
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

} // namespace OpenMS