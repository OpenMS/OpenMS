// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

// The results of ProSEAlgorithm.searchFiles() (bind_prose.cpp), which JavaScript holds: shared by the binding files
// that read them (bind_prose_tables.cpp: the results as typed columns).

#include "js_types.h"

#include <OpenMS/ANALYSIS/ID/ProSEAlgorithm.h>
#include <OpenMS/FORMAT/FASTAFile.h>

#include <optional>
#include <string>
#include <vector>

namespace OpenMS::WASM
{
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

    ProSEResultSet() = default;
    // no copies of the results by mistake (Embind copies a class returned by value)
    ProSEResultSet(const ProSEResultSet&) = delete;
    ProSEResultSet& operator=(const ProSEResultSet&) = delete;

    /// The summary of the search (ProSEMultiSummary)
    JS::ProSEMultiSummary summary() const;

    /// The pooled results as typed columns (ResultTables; bind_prose_tables.cpp)
    JS::ResultTables tables() const;

    /// Stores the pooled run and PSMs as idXML (the ProSE tool's -out_merged)
    void storeIdXML(const std::string& path) const;
  };
} // namespace OpenMS::WASM
