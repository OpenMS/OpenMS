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
#include <iterator>
#include <numeric>
#include <unordered_map>
#include <unordered_set>

#include <boost/math/distributions/normal.hpp>

namespace OpenMS
{

  namespace
  {
    /// The entry of a table that is ordered by mass (the key) whose mass is nearest to @p mass,
    /// or table.end() for an empty table. Of two equally distant entries the lower one is returned.
    template <typename Table>
    typename Table::const_iterator nearestEntry_(const Table& table, double mass)
    {
      const auto upper = table.lower_bound(mass); // the lowest entry that is not below the mass
      if (upper == table.begin())
      {
        return upper; // nothing below (table.end() if the table is empty)
      }
      const auto lower = std::prev(upper);
      if (upper == table.end() || std::abs(lower->first - mass) <= std::abs(upper->first - mass))
      {
        return lower;
      }
      return upper;
    }

    /// The entries of a table that is ordered by mass (the key) that name a cluster of mass @p mass:
    /// the nearest entry if it is within @p tolerance, and all entries within @p tolerance whose distance to
    /// @p mass is at most @p ambiguity_tolerance larger than that of the nearest one (these cannot be told apart).
    /// In order of mass. Empty if no entry is within @p tolerance.
    template <typename Table>
    std::vector<typename Table::const_iterator> namingEntries_(const Table& table, double mass, double tolerance, double ambiguity_tolerance)
    {
      std::vector<typename Table::const_iterator> entries;
      const auto nearest = nearestEntry_(table, mass);
      if (nearest == table.end() || std::abs(nearest->first - mass) > tolerance)
      {
        return entries;
      }

      const double reach = std::min(tolerance, std::abs(nearest->first - mass) + ambiguity_tolerance);
      constexpr double epsilon = 1e-8; // guard band when seeding a lookup at a window border
      for (auto it = table.lower_bound(mass - reach - epsilon); it != table.end() && it->first <= mass + reach + epsilon; ++it)
      {
        if (std::abs(it->first - mass) <= reach)
        {
          entries.push_back(it);
        }
      }
      return entries;
    }

    /// The names of the entries joined with "//" ("A//B")
    template <typename Entries>
    std::string joinedNames_(const Entries& entries)
    {
      std::string joined;
      for (const auto& entry : entries)
      {
        if (!joined.empty())
        {
          joined += "//";
        }
        joined += entry->second;
      }
      return joined;
    }
  } // namespace

  std::pair<OpenSearchModificationAnalysis::DeltaMassHistogram, OpenSearchModificationAnalysis::DeltaMassToChargeCount>
  OpenSearchModificationAnalysis::analyzeDeltaMassPatterns(const PeptideIdentificationList& peptide_ids, 
                                                          bool use_smoothing, 
                                                          bool /*debug*/) const
  {

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
        // Retrieve delta mass and charge (of the target hits)
        if (!hit.metaValueExists("DeltaMass") || hit.isDecoy())
          continue;

        double delta_mass = (double)hit.getMetaValue("DeltaMass");
        int charge = hit.getCharge();

        // Ignore delta masses close to zero
        if (isUnmodified(delta_mass))
          continue;

        // Round delta mass to bin similar values
        double rounded_mass = static_cast<double>(deltaMassBin_(delta_mass)) * DELTA_MASS_BIN_WIDTH_;

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
                                                               const std::string& output_file) const
  {
    std::map<double, std::string, FuzzyDoubleComparator> mass_to_modification(FuzzyDoubleComparator(1e-9));
    std::map<std::string, ModificationPattern> modifications;
    std::map<double, std::string> histogram_found;

    // Load modifications from the database
    std::vector<std::string> modification_names;
    const ModificationsDB* mod_db = ModificationsDB::getInstance();
    mod_db->getAllSearchModifications(modification_names);
    
    for (const std::string& mod_name : modification_names)
    {
      const ResidueModification* residue = mod_db->getModification(mod_name);
      std::string full_name = residue->getFullName();
      double diff_mono_mass = residue->getDiffMonoMass();
      
      if (!full_name.contains("substitution"))
      {
        mass_to_modification[diff_mono_mass] = full_name;
      }
    }

    // Names of the single modifications, in the order in which they are combined below
    std::vector<std::string> single_mod_names;
    single_mod_names.reserve(mass_to_modification.size());
    for (const auto& mod_entry : mass_to_modification)
    {
      single_mod_names.push_back(mod_entry.second);
    }

    // Generate combinations of modifications. Only the pair of indices into
    // single_mod_names is stored; the combined name is assembled on demand for an actual
    // match, since almost none of the ~n^2/2 combinations is ever looked up.
    std::map<double, std::pair<UInt32, UInt32>, FuzzyDoubleComparator> combo_modifications(FuzzyDoubleComparator(1e-9));
    UInt32 index1 = 0;
    for (auto it1 = mass_to_modification.begin(); it1 != mass_to_modification.end(); ++it1, ++index1)
    {
      UInt32 index2 = index1;
      for (auto it2 = it1; it2 != mass_to_modification.end(); ++it2, ++index2)
      {
        combo_modifications[it1->first + it2->first] = {index1, index2};
      }
    }

    // Helper function to add or update modifications
    auto addOrUpdateModification = [&](const std::string& mod_name, double mass, double count, int num_charges)
    {
      if (!modifications.contains(mod_name))
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

    // Map delta masses to modifications
    for (const auto& hist_entry : delta_mass_histogram)
    {
      double cluster_mass = hist_entry.first;
      double count = hist_entry.second;

      bool mapping_found = false;
      std::string mod_name;
      double mod_mass = 0.0;

      // Search in single modifications. The cluster is named after the modification nearest to its mass; the other
      // modifications within the tolerance are not named, as they are less likely to be what was observed. Only if
      // several are about equally close (within MOD_MAPPING_AMBIGUITY_TOL_ of the distance of the nearest one),
      // the data cannot tell them apart, and the cluster is named after all of them ("A//B").
      const auto candidates = namingEntries_(mass_to_modification, cluster_mass, effective_tol, MOD_MAPPING_AMBIGUITY_TOL_);
      if (!candidates.empty())
      {
        // Remember what was found: it explains clusters at higher masses (the +1 isotope variant) and
        // is what the peptide hits are annotated with. Only the modifications that name the cluster qualify.
        for (const auto& candidate : candidates)
        {
          histogram_found[candidate->first] = candidate->second;
        }

        if (candidates.size() == 1)
        {
          mod_name = candidates.front()->second;
          mod_mass = candidates.front()->first;
        }
        else
        {
          mod_name = joinedNames_(candidates);
          mod_mass = cluster_mass;
        }
        mapping_found = true;
      }
      else
      {
        // Check if modification can be explained by known modifications. Only entries within the tolerance of
        // cluster_mass (direct match) or of cluster_mass - 1.0 (+1 isotope variant) qualify, and then the nearest one.
        // effective_tol is capped at MAX_MOD_MAPPING_TOL_, so that the two windows are disjoint and the order of
        // the checks does not matter.
        const auto isotope_hit = nearestEntry_(histogram_found, cluster_mass - 1.0);
        if (isotope_hit != histogram_found.end() &&
            std::abs((isotope_hit->first + 1.0) - cluster_mass) < effective_tol)
        {
          // Explained by a +1 isotope variant
          std::string temp_mod_name = isotope_hit->second + "+1Da";
          addOrUpdateModification(temp_mod_name, isotope_hit->first + 1.0, count, charge_histogram.at(cluster_mass));
          histogram_found[isotope_hit->first + 1.0] = temp_mod_name;
          mapping_found = true;
        }

        if (!mapping_found)
        {
          const auto found_hit = nearestEntry_(histogram_found, cluster_mass);
          if (found_hit != histogram_found.end() &&
              std::abs(found_hit->first - cluster_mass) < effective_tol)
          {
            addOrUpdateModification(found_hit->second, found_hit->first, count, charge_histogram.at(cluster_mass));
            mapping_found = true;
          }
        }

        // Search in combination modifications: the nearest one, whether above or below the cluster mass
        if (!mapping_found)
        {
          const auto combo = nearestEntry_(combo_modifications, cluster_mass);
          if (combo != combo_modifications.end() &&
              std::abs(combo->first - cluster_mass) <= effective_tol)
          {
            mod_name = single_mod_names[combo->second.first] + "++" + single_mod_names[combo->second.second];
            mod_mass = combo->first;
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
        std::string unknown_mod_name = "Unknown" + StringUtils::toStr(std::round(cluster_mass));
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
          
        double delta_mass = (double)hit.getMetaValue("DeltaMass");
        std::string ptm;

        // Check if too close to zero
        if (std::abs(delta_mass) < DELTA_MASS_ZERO_THRESHOLD_)
        {
          hit.setMetaValue("PTM", ptm);
          continue;
        }

        bool found = false;
        // Check with error tolerance if already present in histogram: the entry nearest to delta_mass
        // is the one, if it is within effective_tol.
        const auto entry = nearestEntry_(histogram_found, delta_mass);
        if (entry != histogram_found.end() && std::abs(delta_mass - entry->first) < effective_tol)
        {
          ptm = entry->second;
          found = true;
        }

        // Otherwise assign unknown
        if (!found)
        {
          ptm = "Unknown" + StringUtils::toStr(delta_mass);
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
                                                      const std::string& output_file) const
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

    // Perform Gaussian smoothing. deltas is ascending (extracted from the ordered
    // histogram), so the points within 3 standard deviations form a contiguous
    // range located by binary search — an open-search histogram has 1e4-1e5 bins
    // and a full O(n^2) scan would dominate the analysis.
    const double cutoff = 3.0 * sigma;
    for (size_t i = 0; i < n; ++i)
    {
      double weight_sum = 0.0;

      const size_t j_begin = std::lower_bound(deltas.begin(), deltas.end(), deltas[i] - cutoff) - deltas.begin();
      const size_t j_end = std::upper_bound(deltas.begin(), deltas.end(), deltas[i] + cutoff) - deltas.begin();
      for (size_t j = j_begin; j < j_end; ++j)
      {
        double mz_diff = deltas[i] - deltas[j];

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
                                                                const std::string& output_file) const
  {
    // Remove 'idxml' extension and add '_OutputTable.tsv'
    std::string output_table = output_file;
    if (StringUtils::hasSuffix(output_table, ".idXML"))
    {
      output_table = StringUtils::substr(output_table, 0, output_table.size() - 6) + "_OutputTable.tsv";
    }
    else if (StringUtils::hasSuffix(output_table, ".idxml"))
    {
      output_table = StringUtils::substr(output_table, 0, output_table.size() - 6) + "_OutputTable.tsv";
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
                                                                     const std::string& output_file) const
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
      std::string base_name = output_file;
      if (StringUtils::hasSuffix(base_name, ".idXML") || StringUtils::hasSuffix(base_name, ".idxml"))
      {
        base_name = StringUtils::substr(base_name, 0, base_name.size() - 6);
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

    // Effective tolerance for modification mass matching (capped for open search)
    const double mod_tol = precursor_mass_tolerance_unit_ppm
      ? MAX_MOD_MAPPING_TOL_
      : std::min(precursor_mass_tolerance, MAX_MOD_MAPPING_TOL_);

    // Build modification lookup
    auto mod_lookup = buildModificationMassLookup_();

    // Count the target PSMs (decoys only counted) and, in the same pass, collect the (histogram bin, sequence) of each
    // target hit with a delta mass, for the unique peptide counts of the bins below.
    std::vector<std::pair<Int64, std::string>> bin_sequences;
    for (const auto& peptide_id : peptide_ids)
    {
      for (const auto& hit : peptide_id.getHits())
      {
        if (hit.isDecoy())
        {
          stats.decoy_psms++;
          continue;
        }
        stats.total_psms++;
        if (hit.metaValueExists("DeltaMass"))
        {
          double delta_mass = (double)hit.getMetaValue("DeltaMass");
          if (isUnmodified(delta_mass))
          {
            stats.unmodified_psms++;
          }
          else
          {
            stats.modified_psms++;
          }
          bin_sequences.emplace_back(deltaMassBin_(delta_mass), hit.getSequence().toString());
        }
      }
    }
    // The distinct sequences per bin
    std::sort(bin_sequences.begin(), bin_sequences.end());
    bin_sequences.erase(std::unique(bin_sequences.begin(), bin_sequences.end()), bin_sequences.end());
    std::unordered_map<Int64, int> unique_peptides_of_bin;
    for (const auto& [bin, sequence] : bin_sequences) { ++unique_peptides_of_bin[bin]; }
    bin_sequences = {};

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

      // Count the unique peptides of the PSMs of the entry: those whose delta mass falls into its histogram bin, as
      // they are counted in 'count' (with a smoothed histogram, those of the bin of its maximum). Not those within the
      // precursor mass tolerance: in an open search that is the whole open window, i.e. nearly all PSMs.
      const auto unique = unique_peptides_of_bin.find(deltaMassBin_(delta_mass));
      if (unique != unique_peptides_of_bin.end()) { entry.unique_peptides = unique->second; }

      // Try to map to known modification, named as in mapDeltaMassesToModifications(): after the nearest
      // one within mod_tol, or after all of them ("A//B") if several are about equally close
      const auto mod_candidates = namingEntries_(mod_lookup, delta_mass, mod_tol, MOD_MAPPING_AMBIGUITY_TOL_);
      if (!mod_candidates.empty())
      {
        entry.mapped_modification = joinedNames_(mod_candidates);
        entry.is_known_modification = true;
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
    std::map<std::string, PTMEntry> ptm_map;
    std::map<std::string, std::unordered_set<std::string>> ptm_unique_peptides;
    std::map<std::string, std::unordered_set<int>> ptm_unique_charges;

    // Build modification lookup for theoretical masses
    auto mod_lookup = buildModificationMassLookup_();
    std::map<std::string, double> name_to_mass;
    for (const auto& [mass, name] : mod_lookup)
    {
      name_to_mass[name] = mass;
    }

    // Process all peptide identifications
    for (const auto& peptide_id : peptide_ids)
    {
      for (const auto& hit : peptide_id.getHits())
      {
        if (!hit.metaValueExists("PTM") || hit.isDecoy()) // (of the target hits)
          continue;

        std::string ptm_name = StringUtils::toStr(hit.getMetaValue("PTM"));
        if (ptm_name.empty())
        {
          continue;
        }

        // Check if unknown modification
        if (StringUtils::hasPrefix(ptm_name, "Unknown"))
        {
          stats.unknown_modification_psms++;
          continue;
        }

        stats.total_modified_psms++;

        // Get or create PTM entry
        if (!ptm_map.contains(ptm_name))
        {
          PTMEntry entry;
          entry.name = ptm_name;

          // Get theoretical mass
          if (name_to_mass.contains(ptm_name))
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
          double obs_mass = (double)hit.getMetaValue("DeltaMass");
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

        double hit_delta_mass = (double)hit.getMetaValue("DeltaMass");
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
                                                                const std::string& output_file) const
  {
    std::ofstream output_stream(output_file);
    if (!output_stream.is_open())
    {
      OPENMS_LOG_ERROR << "Error opening file: " << output_file << std::endl;
      return;
    }

    // Write header
    output_stream << "# Delta Mass Statistics Table\n";
    output_stream << "# Counts of the target PSMs (the decoy PSMs are left out)\n";
    output_stream << "# Total PSMs: " << stats.total_psms << "\n";
    output_stream << "# Modified PSMs: " << stats.modified_psms << "\n";
    output_stream << "# Unmodified PSMs: " << stats.unmodified_psms << "\n";
    output_stream << "# Decoy PSMs: " << stats.decoy_psms << "\n";
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
                                                          const std::string& output_file) const
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

  std::map<double, std::string, OpenSearchModificationAnalysis::FuzzyDoubleComparator>
  OpenSearchModificationAnalysis::buildModificationMassLookup_() const
  {
    std::map<double, std::string, FuzzyDoubleComparator> mass_to_mod; // uses default epsilon (1e-9)

    std::vector<std::string> modification_names;
    const ModificationsDB* mod_db = ModificationsDB::getInstance();
    mod_db->getAllSearchModifications(modification_names);

    for (const std::string& mod_name : modification_names)
    {
      const ResidueModification* residue = mod_db->getModification(mod_name);
      std::string full_name = residue->getFullName();
      double diff_mono_mass = residue->getDiffMonoMass();

      // Skip substitutions
      if (!full_name.contains("substitution"))
      {
        mass_to_mod[diff_mono_mass] = full_name;
      }
    }

    return mass_to_mod;
  }

  std::string OpenSearchModificationAnalysis::getTargetResidues_(const std::string& mod_name) const
  {
    // Split compound names (e.g. "Oxidation//Deamidated" or "Oxidation+1Da") and resolve each part
    std::vector<std::string> parts;
    if (mod_name.contains("//"))
    {
      StringUtils::split(mod_name, "//", parts);
    }
    else
    {
      // Strip isotope suffixes like "+1Da", "+2Da" etc.
      std::string base_name = mod_name;
      auto pos = mod_name.find("+");
      if (pos != std::string::npos && pos > 0)
      {
        std::string suffix = StringUtils::substr(mod_name, pos);
        if (StringUtils::hasSuffix(suffix, "Da"))
        {
          base_name = StringUtils::substr(mod_name, 0, pos);
        }
      }
      parts.push_back(base_name);
    }

    std::vector<std::string> residues;
    const ModificationsDB* mod_db = ModificationsDB::getInstance();

    for (const auto& part : parts)
    {
      std::string trimmed = part;
      StringUtils::trim(trimmed);
      if (trimmed.empty()) continue;

      try
      {
        const ResidueModification* mod = mod_db->getModification(trimmed);
        if (mod != nullptr)
        {
          std::string result;
          std::string origin(1, mod->getOrigin());
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
    std::string result;
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

        double hit_delta_mass = (double)hit.getMetaValue("DeltaMass");
        if (std::abs(hit_delta_mass - delta_mass) <= tolerance)
        {
          unique_sequences.insert(hit.getSequence().toString());
        }
      }
    }

    return static_cast<int>(unique_sequences.size());
  }

} // namespace OpenMS
