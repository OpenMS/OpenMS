// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// The results of ProSEAlgorithm.searchFiles() (bind_prose.cpp), which JavaScript holds: shared by the binding files
// that read them (bind_prose_tables.cpp: the results as typed columns, and the peaks of the identified spectra).

#include "js_types.h"

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <cstddef>
#include <cstdint>
#include <optional>
#include <string>
#include <vector>

namespace OpenMS::WASM
{
  /// Values appended in blocks of a fixed size: unlike a std::vector, it grows without moving its values, so that it
  /// never holds a second buffer while it grows (for the peaks of many spectra, in the 4 GB of a wasm32 module)
  template <typename T>
  struct BlockVector
  {
    static constexpr std::size_t block_size = std::size_t(1) << 20;
    std::vector<std::vector<T>> blocks;
    std::size_t size = 0;

    void push_back(T value)
    {
      if (blocks.empty() || blocks.back().size() == block_size)
      {
        blocks.emplace_back();
        blocks.back().reserve(block_size);
      }
      blocks.back().push_back(value);
      ++size;
    }

    /// Frees the capacity of the last block that is not used
    void shrink_to_fit()
    {
      if (!blocks.empty()) blocks.back().shrink_to_fit();
    }
  };

  /// The peaks of the identified spectra of a ProSEResultSet (spectra()), read from its input files
  struct ProSEStoredSpectra
  {
    std::vector<int32_t> psm_spectrum;     ///< per PSM row of tables(): the index of its spectrum, -1: none
    std::vector<uint32_t> offsets{0};      ///< the peaks of spectrum s: offsets[s] .. offsets[s + 1] of mz and intensity
    BlockVector<float> mz, intensity;      ///< as read, sorted by m/z
    std::vector<double> precursor_mz;      ///< of the first precursor; NaN: none
    std::vector<float> collision_energy;   ///< the activation energy of the first precursor; NaN: unknown (none, or 0)
  };

  /// The search of one input file of ProSEAlgorithm.searchFiles()
  struct ProSEFileResult
  {
    std::string path;                                   ///< as given
    ProSEAlgorithm::ExitCodes exit_code = ProSEAlgorithm::ExitCodes::EXECUTION_OK;
    ProSEAlgorithm::RunStatistics stats;                ///< of the PSMs of the file (after its FDR:PSM)
    ProSEAlgorithm::PercolatorRescoring percolator;
    std::optional<double> peptdeep_nce;                 ///< the NCE of PeptDeep's predictions (PeptDeepRescoring::NCE_META_VALUE)
  };

  /// The results of ProSEAlgorithm.searchFiles(), which JavaScript holds until delete(): the search of each input
  /// file and the pooled evaluation of all of them (ProSEAlgorithm::mergeFileResults(), the ProSE tool's -out_merged).
  struct ProSEResultSet
  {
    std::vector<ProSEFileResult> files;      ///< in input order
    /// the pooled run and the PSMs of all files; the meta value 'id_merge_index' of a PSM is the position of its file
    ProSEAlgorithm::MergedResult merged;
    std::vector<FASTAFile::FASTAEntry> database; ///< the proteins of the FASTA file searched
    bool rescore = false;                    ///< whether the search rescored the PSMs with Percolator ('rescore')
    bool open_search = false;
    /// of the PSMs of all files (open search)
    OpenSearchModificationAnalysis::OpenSearchAnalysisResult modification_analysis;
    std::string decoy_string;                ///< the decoy marker of the accessions (empty: target-only)
    bool decoy_is_prefix = true;             ///< whether decoy_string is a prefix (else a suffix)
    double psm_fdr = 0.0;                    ///< FDR:PSM
    double protein_fdr = 0.0;                ///< FDR:protein
    bool keep_spectra = true;                ///< ProSEFilesOptions.keepSpectra

    /// How ProSE built the decoys it searched, to give the decoy proteins their sequences (tables()): whether it
    /// generated them from the targets of the database (otherwise they are in it), and the settings it used
    struct DecoyGeneration
    {
      bool generated = false;                ///< the search's decoy mode: "generated"
      std::string enzyme;                    ///< parameter 'enzyme'
      bool unspecific = false;               ///< parameter 'peptide:enzyme_specificity' "none": whole proteins reversed
      bool keep_initial_met = false;         ///< parameter 'peptide:clip_nterm_methionine': an initial M stays in place
    } decoys;

    /// The Percolator features of the search engine of the runs of the files (their 'extra_features', in order of
    /// appearance), for the feature columns of tables()
    std::vector<std::string> extra_features;
    /// The precursor charges of the search (precursor:min_charge, precursor:max_charge), for the charge features of
    /// Percolator in tables()
    int min_charge = 0, max_charge = 0;

    /// The peaks of the spectra of the PSMs (collectSpectra())
    ProSEStoredSpectra stored_spectra;

    ProSEResultSet() = default;
    // no copies of the results by mistake (Embind copies a class returned by value)
    ProSEResultSet(const ProSEResultSet&) = delete;
    ProSEResultSet& operator=(const ProSEResultSet&) = delete;

    /// The summary of the search (ProSEMultiSummary)
    JS::ProSEMultiSummary summary() const;

    /// The pooled results as typed columns (ResultTables; bind_prose_tables.cpp)
    JS::ResultTables tables() const;

    /// Reads the peaks of the spectra of the PSM rows of tables() from the input files into stored_spectra
    /// (bind_prose_tables.cpp), at the end of searchFiles(), while the files are there. With keep_spectra: the MS2
    /// spectra of the rows whose q-value passes FDR:PSM (all rows without q-values, or with FDR:PSM 0), one per file
    /// and native ID, streamed from one file at a time; otherwise none (every row -1).
    void collectSpectra();

    /// The peaks of the identified spectra (StoredSpectra; bind_prose_tables.cpp), copied into JavaScript arrays
    JS::StoredSpectra spectra() const;

    /// For tests: replaces the native ID of PSM row @p row of tables() and collects the spectra again
    /// (collectSpectra()), which re-reads the input files (bind_prose_tables.cpp)
    void debugSetNativeId(int row, const std::string& native_id);

    /// Stores the pooled run and PSMs as idXML (the ProSE tool's -out_merged)
    void storeIdXML(const std::string& path) const;
  };
} // namespace OpenMS::WASM
