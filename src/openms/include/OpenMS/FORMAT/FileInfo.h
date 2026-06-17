// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>
#include <OpenMS/CONCEPT/Types.h>           // Int / UInt / UInt64
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/MATH/StatisticFunctions.h> // Math::SummaryStatistics

#include <cstdint>
#include <iosfwd>
#include <map>
#include <optional>
#include <string>
#include <tuple>
#include <utility>
#include <vector>

namespace OpenMS
{
  /**
    @brief Library-level equivalent of the @ref TOPP_FileInfo tool.

    Loads a mass-spectrometry-related file (type auto-detected unless forced via
    Options::forced_type), extracts all file-level information into a structured
    @ref FileInfo::Result, and can render that result as the exact human-readable
    text / TSV that the FileInfo command line tool emits.

    Computation (run()) is fully separated from formatting (toText()/toTSV()), so
    the result is consumable directly from C++ and pyOpenMS without any stream.

    @code
    FileInfo fi;
    FileInfo::Result r = fi.runAll("data.consensusXML");
    std::cout << r.feature->num_features << "\n";
    std::cout << FileInfo::toText(r);     // == the FileInfo CLI output
    @endcode

    @ingroup FileIO
  */
  class OPENMS_DLLAPI FileInfo
  {
  public:
    FileInfo() = default;
    ~FileInfo() = default;

    // ---------------------------------------------------------------------
    // nested value types (all pyOpenMS-bindable: plain structs / vector / map)
    // ---------------------------------------------------------------------

    /// one [min,max] interval for a single dimension; present==false => "<none>"
    struct OPENMS_DLLAPI Range
    {
      bool present = false;
      double min = 0.0;
      double max = 0.0;
    };

    /// the four range dimensions FileInfo reports for a map / spectrum group
    struct OPENMS_DLLAPI RangeSet
    {
      Range rt;        ///< retention time (seconds)
      Range mz;        ///< mass-to-charge
      Range mobility;  ///< ion mobility (only filled when the source carries it)
      Range intensity;
      bool has_mobility = false; ///< whether the mobility dimension applies at all
    };

    /// MSExperiment carries four range categories; non-MSExperiment maps fill only @p combined
    struct OPENMS_DLLAPI Ranges
    {
      RangeSet combined;                              ///< spectra + chromatograms (or the whole map)
      RangeSet spectra_overall;                       ///< MSExperiment only
      std::map<UInt, RangeSet> per_ms_level;          ///< MSExperiment only, keyed by MS level
      RangeSet chromatograms;                         ///< MSExperiment only
      bool is_experiment = false;                     ///< true => the per-level / chrom sections apply
    };

    /// general header block (always populated)
    struct OPENMS_DLLAPI FileMeta
    {
      std::string file_name;
      FileTypes::Type file_type = FileTypes::UNKNOWN;
      std::string file_type_name;                     ///< FileTypes::typeToName(file_type)
    };

    /// experiment / instrument / sample / contact metadata (the @c -m block, peak files)
    struct OPENMS_DLLAPI ExperimentMeta
    {
      bool present = false;
      std::string document_id;
      std::string date;
      // sample
      std::string sample_name, sample_organism, sample_comment;
      // instrument
      std::string instrument_name, instrument_model, instrument_vendor;
      std::vector<std::string> ion_sources, mass_analyzers, detectors;
      // contacts
      struct OPENMS_DLLAPI Contact { std::string first_name, last_name, email; };
      std::vector<Contact> contacts;
    };

    /// one data-processing step (the @c -p block)
    struct OPENMS_DLLAPI ProcessingStep
    {
      std::string software_name, software_version, completion_time;
      std::vector<std::string> actions;
    };

    /// one named SummaryStatistics block; @p title is the label used in both renderers
    struct OPENMS_DLLAPI NamedStats
    {
      std::string title;
      Math::SummaryStatistics<std::vector<double>> stats;
    };

    /// peak-file (MSExperiment) specifics
    struct OPENMS_DLLAPI PeakInfo
    {
      std::string instrument_name;
      std::vector<std::pair<std::string, double>> mass_analyzers;  ///< (type, resolution)
      std::vector<Int> ms_levels;
      UInt64 total_peaks = 0;                                  ///< incl. chromatographic peaks
      UInt64 num_spectra = 0;
      std::map<Int, UInt64> spectra_per_ms_level;
      std::map<Int, std::string> peak_type_per_ms_level;      ///< annotation/estimate per MS level
      /// activation methods: (ms_level, method_full_name) -> count
      std::map<std::pair<Int, std::string>, UInt64> activation_methods;
      std::map<Int, UInt64> precursor_charges;                ///< charge -> count
      std::map<std::string, UInt64> float_arrays, int_arrays, string_arrays; ///< name -> count
      std::vector<double> faims_cvs;
      UInt64 num_chromatograms = 0;
      UInt64 num_chrom_peaks = 0;
      std::map<std::string, UInt64> chromatogram_types;       ///< type name -> count

      /// flattened (ms_level, method, count) view for ergonomic binding (mirrors the TSV columns)
      std::vector<std::tuple<Int, std::string, UInt64>> activationMethodsFlat() const;
    };

    /// feature / consensus specifics
    struct OPENMS_DLLAPI FeatureInfo
    {
      bool is_consensus = false;
      UInt64 num_features = 0;                       ///< number of (consensus) features
      double tic = 0.0;                              ///< total ion current (featureXML)
      std::map<Int, UInt64> charges;                 ///< charge -> count
      std::map<UInt64, UInt64> ids_per_element;      ///< #IDs -> #elements
      UInt64 assigned_ids = 0, unassigned_ids = 0;
      // consensus only:
      std::map<UInt64, UInt64> size_distribution;    ///< consensus-feature size -> count
      struct OPENMS_DLLAPI MapColumn { std::string filename, identifier, label; UInt64 size = 0; };
      std::vector<MapColumn> map_columns;
    };

    /// identification specifics (idXML / mzIdentML)
    struct OPENMS_DLLAPI IdentInfo
    {
      std::string db_name, db_version, taxonomy;
      std::vector<std::string> search_engines;       ///< "engine (version)"
      UInt64 num_runs = 0, protein_hits = 0, non_redundant_protein_hits = 0;
      UInt64 matched_spectra = 0, peptide_hits = 0;
      double psms_per_spectrum = 0.0, avg_peptide_length = 0.0;
      UInt64 non_redundant_peptides = 0, modified_tophits = 0;
      std::map<std::string, UInt64> modification_counts;
    };

    /// FASTA specifics
    struct OPENMS_DLLAPI FastaInfo
    {
      UInt64 num_sequences = 0, total_residues = 0;
      bool is_nucleic_acid = false;
      Math::SummaryStatistics<std::vector<double>> length_stats;
      std::map<char, UInt64> residue_counts;
      UInt64 seq_with_ambiguous = 0, dup_headers = 0, dup_sequences = 0;
      std::map<std::string, UInt64> ambiguity_counts; ///< buckets depend on is_nucleic_acid
    };

    /// mzTab specifics
    struct OPENMS_DLLAPI MzTabInfo
    {
      std::string version, mode, type;
      UInt64 psms = 0, peptides = 0, proteins = 0, oligonucleotides = 0, osms = 0,
             small_molecules = 0, nucleic_acids = 0;
    };

    /// the @c -v (schema/semantic validation) and @c -i (indexed mzML) blocks
    struct OPENMS_DLLAPI ValidationInfo
    {
      bool performed = false, supported = true, valid = false;
      std::string schema_version;
      std::string detail;                            ///< captured validator output (re-emitted verbatim)
      std::vector<std::string> warnings, errors;     ///< semantic validation
      // indexed-mzML (-i)
      bool index_checked = false, index_valid = false;
      UInt64 indexed_spectra = 0, indexed_chromatograms = 0;
    };

    /// the @c -c (corrupt data) block
    struct OPENMS_DLLAPI CorruptionInfo
    {
      bool performed = false;
      std::vector<std::string> errors, warnings;     ///< already-formatted message lines
    };

    /// detailed per-spectrum listing (the @c -d block); kept as pre-rendered lines
    struct OPENMS_DLLAPI DetailInfo
    {
      bool performed = false;
      std::vector<std::string> lines;
    };

    /// the master result aggregate
    struct OPENMS_DLLAPI Result
    {
      FileMeta                      meta;
      Ranges                        ranges;          ///< filled for peak / feature / consensus
      std::optional<PeakInfo>       peak;
      std::optional<FeatureInfo>    feature;
      std::optional<IdentInfo>      ident;
      std::optional<FastaInfo>      fasta;
      std::optional<MzTabInfo>      mztab;
      std::optional<ExperimentMeta> experiment_meta;
      std::vector<ProcessingStep>   processing;
      std::vector<NamedStats>       statistics;
      ValidationInfo                validation;
      CorruptionInfo                corruption;
      DetailInfo                    detail;
      std::string                   transformation_summary; ///< trafoXML: model + printSummary
      std::string                   targeted_summary;       ///< PQP: getSummary()

      /// Human-readable rendering, identical to the FileInfo CLI @c -out output (filled by run()).
      std::string                   text;
      /// TSV rendering, identical to the FileInfo CLI @c -out_tsv output (filled by run()).
      std::string                   tsv;
    };

    /// what to compute (mirrors the CLI flags; lets callers opt into expensive work)
    struct OPENMS_DLLAPI Options
    {
      FileTypes::Type forced_type = FileTypes::UNKNOWN; ///< == CLI -in_type
      bool meta = false;          ///< -m
      bool processing = false;    ///< -p
      bool statistics = false;    ///< -s
      bool detailed = false;      ///< -d
      bool check_corrupt = false; ///< -c
      bool validate = false;      ///< -v
      bool check_index = false;   ///< -i
    };

    // ---------------------------------------------------------------------
    // public API (all pyOpenMS-bindable)
    // ---------------------------------------------------------------------

    /**
      @brief Load @p filename (type auto-detected unless Options::forced_type set) and compute everything requested.
      @throw Exception::FileNotFound if @p filename does not exist
      @throw Exception::ParseError / format-specific exceptions if the file cannot be read
      (the loaders are not caught here; the CLI wrapper relies on TOPPBase to convert these to exit codes).
    */
    Result run(const std::string& filename, const Options& options);

    /// Load @p filename with default options (no extra flags). @see run(filename, options) for the exception contract.
    Result run(const std::string& filename) { return run(filename, Options()); }

    /// Convenience: compute all content metrics (meta/processing/statistics on; no v/i/d/c). @see run for throws.
    Result runAll(const std::string& filename);

    /// Return the human-readable rendering produced by run() (== the FileInfo CLI @c -out output, cached in @p r).
    /// @note @p options is accepted for API symmetry only; the rendering reflects the Options passed to run().
    static std::string toText(const Result& r, const Options& options);

    /// Return the human-readable rendering produced by run() (== the FileInfo CLI @c -out output, cached in @p r).
    static std::string toText(const Result& r) { return toText(r, Options()); }

    /// Return the TSV rendering produced by run() (== the FileInfo CLI @c -out_tsv output, cached in @p r).
    /// @note @p options is accepted for API symmetry only; the rendering reflects the Options passed to run().
    static std::string toTSV(const Result& r, const Options& options);

    /// Return the TSV rendering produced by run() (== the FileInfo CLI @c -out_tsv output, cached in @p r).
    static std::string toTSV(const Result& r) { return toTSV(r, Options()); }

  private:
    /// Single computation+rendering core: loads @p in (already type-detected as @p in_type),
    /// writes the human-readable / TSV report to @p os / @p os_tsv, and populates @p r.
    /// Mirrors the historical TOPPFileInfo::outputTo_ so CLI output stays byte-identical.
    void report_(const std::string& in, FileTypes::Type in_type, const Options& o,
                 std::ostream& os, std::ostream& os_tsv, Result& r);
  };
} // namespace OpenMS
