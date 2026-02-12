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
   *
   * The class can generate two types of statistics tables:
   * 1. PTM Statistics Table - Shows known modifications mapped to delta masses
   * 2. Delta Mass Statistics Table - Shows raw delta mass distribution analysis
   *
   * These tables can be used for modification discovery in open search workflows.
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

    /// Statistics for a single delta mass bin in the histogram
    struct DeltaMassEntry
    {
      double delta_mass = 0.0;            ///< Central delta mass value
      int count = 0;                      ///< Number of PSMs with this delta mass
      int unique_peptides = 0;            ///< Number of unique peptide sequences
      int num_charge_states = 0;          ///< Number of different charge states observed
      double percentage = 0.0;            ///< Percentage of total PSMs
      String mapped_modification = "";    ///< Name of mapped modification (if any)
      bool is_known_modification = false; ///< Whether this maps to a known modification
    };

    /// Statistics for a mapped PTM
    struct PTMEntry
    {
      String name;                        ///< Modification name (e.g., "Oxidation (M)")
      double theoretical_mass = 0.0;      ///< Theoretical mass shift from ModificationsDB
      double observed_mass = 0.0;         ///< Mean observed delta mass
      double mass_deviation = 0.0;        ///< Deviation between theoretical and observed
      int count = 0;                      ///< Number of PSMs with this modification
      int unique_peptides = 0;            ///< Number of unique peptide sequences
      int num_charge_states = 0;          ///< Number of different charge states observed
      double percentage = 0.0;            ///< Percentage of total modified PSMs
      std::map<char, int> residue_counts; ///< Count per amino acid residue
      String target_residues;             ///< Target residues for this modification
    };

    /// Container for delta mass statistics table
    struct DeltaMassStatistics
    {
      std::vector<DeltaMassEntry> entries;  ///< All delta mass entries
      int total_psms = 0;                   ///< Total number of PSMs analyzed
      int modified_psms = 0;                ///< PSMs with non-zero delta mass
      int unmodified_psms = 0;              ///< PSMs with ~zero delta mass
      double mean_delta_mass = 0.0;         ///< Mean delta mass (excluding unmodified)
      double median_delta_mass = 0.0;       ///< Median delta mass (excluding unmodified)
    };

    /// Container for PTM statistics table
    struct PTMStatistics
    {
      std::vector<PTMEntry> entries;        ///< All PTM entries
      int total_modified_psms = 0;          ///< Total PSMs with mapped modifications
      int unknown_modification_psms = 0;    ///< PSMs with unknown modifications
      int num_unique_modifications = 0;     ///< Number of distinct modifications found
    };

    /// Combined result of open search modification analysis
    struct OpenSearchAnalysisResult
    {
      DeltaMassStatistics delta_mass_stats; ///< Delta mass histogram statistics
      PTMStatistics ptm_stats;              ///< Mapped PTM statistics
      std::vector<ModificationSummary> summaries; ///< Legacy modification summaries
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
     * @param[in,out] peptide_ids List of peptide identifications containing delta mass information
     * @param[in] use_smoothing Whether to apply smoothing to the delta mass histogram
     * @param[out] debug Enable debug output
     * @return Pair containing delta mass histogram and charge state counts
     */
    std::pair<DeltaMassHistogram, DeltaMassToChargeCount> 
    analyzeDeltaMassPatterns(const PeptideIdentificationList& peptide_ids, 
                            bool use_smoothing = false, 
                            bool debug = false) const;

    /**
     * @brief Map delta masses to known modifications and annotate peptides
     * 
     * @param[in] delta_mass_histogram Histogram of delta masses
     * @param[in] charge_histogram Charge state counts for each delta mass
     * @param[in,out] peptide_ids List of peptide identifications to annotate (modified in-place)
     * @param[in] precursor_mass_tolerance Mass tolerance for mapping
     * @param[in] precursor_mass_tolerance_unit_ppm Whether tolerance is in ppm (true) or Da (false)
     * @param[in] output_file Optional file path for writing modification summary table
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
     * @param[in] peptide_ids List of peptide identifications (modified in-place)
     * @param[in] precursor_mass_tolerance Mass tolerance for mapping
     * @param[in] precursor_mass_tolerance_unit_ppm Whether tolerance is in ppm (true) or Da (false)
     * @param[in] use_smoothing Whether to apply smoothing to delta mass histogram
     * @param[in] output_file Optional file path for writing modification summary table
     * @return List of modification summaries found
     */
    std::vector<ModificationSummary>
    analyzeModifications(PeptideIdentificationList& peptide_ids,
                        double precursor_mass_tolerance = 5.0,
                        bool precursor_mass_tolerance_unit_ppm = true,
                        bool use_smoothing = false,
                        const String& output_file = "") const;

    /**
     * @brief Complete analysis returning structured statistics tables
     *
     * This is the main entry point for fragment index open search modification discovery.
     * It performs a complete analysis workflow and returns structured tables containing:
     * - Delta mass statistics (histogram of mass shifts)
     * - PTM statistics (mapped modifications with residue localization)
     *
     * @param peptide_ids List of peptide identifications (modified in-place with PTM annotations)
     * @param precursor_mass_tolerance Mass tolerance for mapping delta masses to known modifications
     * @param precursor_mass_tolerance_unit_ppm Whether tolerance is in ppm (true) or Da (false)
     * @param use_smoothing Whether to apply Gaussian smoothing to delta mass histogram
     * @param output_file Optional file path for writing CSV/TSV output tables
     * @return OpenSearchAnalysisResult containing delta mass and PTM statistics tables
     */
    OpenSearchAnalysisResult
    analyzeModificationsWithStatistics(PeptideIdentificationList& peptide_ids,
                                       double precursor_mass_tolerance = 5.0,
                                       bool precursor_mass_tolerance_unit_ppm = true,
                                       bool use_smoothing = false,
                                       const String& output_file = "") const;

    /**
     * @brief Generate delta mass statistics table from histogram data
     *
     * Converts the raw delta mass histogram into a structured statistics table
     * with additional computed metrics like percentages and unique peptide counts.
     *
     * @param histogram Delta mass histogram from analyzeDeltaMassPatterns()
     * @param charge_histogram Charge state counts per delta mass
     * @param peptide_ids Peptide identifications for computing unique peptide counts
     * @param precursor_mass_tolerance Mass tolerance for grouping
     * @param precursor_mass_tolerance_unit_ppm Whether tolerance is in ppm
     * @return DeltaMassStatistics table with all entries and summary statistics
     */
    DeltaMassStatistics
    generateDeltaMassStatistics(const DeltaMassHistogram& histogram,
                                const DeltaMassToChargeCount& charge_histogram,
                                const PeptideIdentificationList& peptide_ids,
                                double precursor_mass_tolerance = 5.0,
                                bool precursor_mass_tolerance_unit_ppm = true) const;

    /**
     * @brief Generate PTM statistics table with residue localization
     *
     * Analyzes peptide identifications to generate a table of mapped PTMs
     * including residue-specific localization analysis.
     *
     * @param peptide_ids Peptide identifications with PTM annotations
     * @param precursor_mass_tolerance Mass tolerance for mapping
     * @param precursor_mass_tolerance_unit_ppm Whether tolerance is in ppm
     * @return PTMStatistics table with all mapped modifications
     */
    PTMStatistics
    generatePTMStatistics(const PeptideIdentificationList& peptide_ids,
                          double precursor_mass_tolerance = 5.0,
                          bool precursor_mass_tolerance_unit_ppm = true) const;

    /**
     * @brief Analyze which amino acid residues are associated with a delta mass
     *
     * For each delta mass, examines the peptide sequences to determine which
     * amino acids are most frequently present, helping to localize modifications.
     *
     * @param peptide_ids Peptide identifications with DeltaMass meta values
     * @param delta_mass Target delta mass to analyze
     * @param tolerance Mass tolerance for matching
     * @return Map from amino acid character to occurrence count
     */
    std::map<char, int>
    analyzeResidueFrequency(const PeptideIdentificationList& peptide_ids,
                            double delta_mass,
                            double tolerance = 0.01) const;

    /**
     * @brief Write delta mass statistics to a TSV file
     *
     * @param stats Delta mass statistics to write
     * @param output_file Output file path
     */
    void writeDeltaMassStatistics(const DeltaMassStatistics& stats,
                                  const String& output_file) const;

    /**
     * @brief Write PTM statistics to a TSV file
     *
     * @param stats PTM statistics to write
     * @param output_file Output file path
     */
    void writePTMStatistics(const PTMStatistics& stats,
                            const String& output_file) const;

  private:

    /// Maximum tolerance (Da) for matching delta masses to known modification masses.
    /// Prevents overly broad matching when the precursor tolerance is large (e.g. open search).
    static constexpr double MAX_MOD_MAPPING_TOL_ = 0.02; // Da

    /// Delta masses within this threshold of zero are considered unmodified.
    static constexpr double DELTA_MASS_ZERO_THRESHOLD_ = 0.05; // Da

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

    /// Build lookup table mapping mass differences to known modifications
    std::map<double, String, FuzzyDoubleComparator>
    buildModificationMassLookup_() const;

    /// Get target residues for a modification name
    String getTargetResidues_(const String& mod_name) const;

    /// Count unique peptide sequences matching a delta mass
    int countUniquePeptides_(const PeptideIdentificationList& peptide_ids,
                            double delta_mass,
                            double tolerance) const;
  };

} // namespace OpenMS
