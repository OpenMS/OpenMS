// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/DATASTRUCTURES/String.h>

#include <map>
#include <vector>
#include <unordered_map>
#include <unordered_set>

namespace OpenMS
{
  /**
   * @brief Utility class for analyzing modification patterns in open search results
   * 
   * This class provides functionality to analyze delta mass patterns from open search
   * peptide identifications, identify common modifications, and map them to known
   * modifications from the ModificationsDB. Originally extracted from SageAdapter.
   */
  class OPENMS_DLLAPI OpenSearchModificationAnalysis
  {
  public:
    
    /// Stores details of a modification pattern found in the data
    struct ModificationPattern
    {
      double count = 0.0;           ///< Number of peptides with this modification
      std::vector<double> masses;   ///< Masses associated with the modification
      int num_charge_states = 0;    ///< Number of different charge states observed
    };

    /// Data structure for modification summary output
    struct ModificationSummary
    {
      int count;                    ///< Modification rate (number of occurrences)
      String name;                  ///< Modification name
      int num_charge_states;        ///< Number of charge states
      std::vector<double> masses;   ///< Masses associated with the modification
    };

    /// Comparator for approximate comparison of double values
    struct FuzzyDoubleComparator
    {
      double epsilon;
      FuzzyDoubleComparator(double eps = 1e-9) : epsilon(eps) {}
      bool operator()(const double& a, const double& b) const 
      {
        return std::fabs(a - b) >= epsilon && a < b;
      }
    };

    /// Type definitions for delta mass analysis
    using DeltaMassHistogram = std::map<double, double, FuzzyDoubleComparator>;
    using DeltaMassToChargeCount = std::map<double, int, FuzzyDoubleComparator>;

    /// Default constructor
    OpenSearchModificationAnalysis() = default;

    /// Destructor
    ~OpenSearchModificationAnalysis() = default;

    /**
     * @brief Analyze delta mass patterns from peptide identifications
     * 
     * @param peptide_ids List of peptide identifications containing delta mass information
     * @param use_smoothing Whether to apply smoothing to the delta mass histogram
     * @param debug Enable debug output
     * @return Pair containing delta mass histogram and charge state counts
     */
    std::pair<DeltaMassHistogram, DeltaMassToChargeCount> 
    analyzeDeltaMassPatterns(const PeptideIdentificationList& peptide_ids, 
                            bool use_smoothing = false, 
                            bool debug = false) const;

    /**
     * @brief Map delta masses to known modifications and annotate peptides
     * 
     * @param delta_mass_histogram Histogram of delta masses
     * @param charge_histogram Charge state counts for each delta mass
     * @param peptide_ids List of peptide identifications to annotate (modified in-place)
     * @param precursor_mass_tolerance Mass tolerance for mapping
     * @param precursor_mass_tolerance_unit_ppm Whether tolerance is in ppm (true) or Da (false)
     * @param output_file Optional file path for writing modification summary table
     * @return List of modification summaries found
     */
    std::vector<ModificationSummary>
    mapDeltaMassesToModifications(const DeltaMassHistogram& delta_mass_histogram,
                                 const DeltaMassToChargeCount& charge_histogram,
                                 PeptideIdentificationList& peptide_ids,
                                 double precursor_mass_tolerance = 5.0,
                                 bool precursor_mass_tolerance_unit_ppm = true,
                                 const String& output_file = "") const;

    /**
     * @brief Complete analysis workflow: analyze patterns and map to modifications
     * 
     * @param peptide_ids List of peptide identifications (modified in-place)
     * @param precursor_mass_tolerance Mass tolerance for mapping
     * @param precursor_mass_tolerance_unit_ppm Whether tolerance is in ppm (true) or Da (false)
     * @param use_smoothing Whether to apply smoothing to delta mass histogram
     * @param output_file Optional file path for writing modification summary table
     * @return List of modification summaries found
     */
    std::vector<ModificationSummary>
    analyzeModifications(PeptideIdentificationList& peptide_ids,
                        double precursor_mass_tolerance = 5.0,
                        bool precursor_mass_tolerance_unit_ppm = true,
                        bool use_smoothing = false,
                        const String& output_file = "") const;

  private:
    
    /// Gaussian function for smoothing
    static double gaussian_(double x, double sigma);

    /// Smooth delta mass histogram using Gaussian kernel density estimation
    static DeltaMassHistogram smoothDeltaMassHistogram_(const DeltaMassHistogram& histogram, 
                                                       double sigma = 0.001);

    /// Find peaks in delta mass histogram based on count threshold and signal-to-noise ratio
    static DeltaMassHistogram findPeaksInHistogram_(const DeltaMassHistogram& histogram, 
                                                   double count_threshold = 0.0, 
                                                   double snr = 2.0);

    /// Write modification summary table to file
    void writeModificationSummary_(const std::vector<ModificationSummary>& modifications,
                                  const String& output_file) const;
  };

} // namespace OpenMS