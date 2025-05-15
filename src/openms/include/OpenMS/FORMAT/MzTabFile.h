// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/MzTab.h>

#include <OpenMS/METADATA/PeptideHit.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <vector>

namespace OpenMS
{
  class String;
  class SVOutStream;
/**
    @brief File adapter for MzTab files

    @ingroup FileIO
  */
  class OPENMS_DLLAPI MzTabFile
  {
  public:
    ///Default constructor
    MzTabFile();
    ///Destructor
    ~MzTabFile();

    typedef std::map<std::pair<String, String>, std::vector<PeptideHit> > MapAccPepType;

    // store MzTab file
    void store(const String& filename, const MzTab& mz_tab) const;

    // stream IDs to file
    void store(
          const String& filename,
          const std::vector<ProteinIdentification>& protein_identifications,
          const std::vector<PeptideIdentification>& peptide_identifications,
          bool first_run_inference_only,
          bool export_empty_pep_ids = false,
          bool export_all_psms = false,
          const String& title = "ID export from OpenMS");

    // stream ConsensusMap to file
    void store(
      const String& filename, 
      const ConsensusMap& cmap,
      const bool first_run_inference_only,
      const bool export_unidentified_features,
      const bool export_unassigned_ids,
      const bool export_subfeatures,
      const bool export_empty_pep_ids = false,
      const bool export_all_psms = false) const;

    // Set store behaviour of optional "reliability" and "uri" columns (default=no)
    void storeProteinReliabilityColumn(bool store);
    void storePeptideReliabilityColumn(bool store);
    void storePSMReliabilityColumn(bool store);
    void storeSmallMoleculeReliabilityColumn(bool store);
    void storeProteinUriColumn(bool store);
    void storePeptideUriColumn(bool store);
    void storePSMUriColumn(bool store);
    void storeSmallMoleculeUriColumn(bool store);
    void storeProteinGoTerms(bool store);

    // load MzTab file
    void load(const String& filename, MzTab& mz_tab);

  protected:
    bool store_protein_reliability_;
    bool store_peptide_reliability_;
    bool store_psm_reliability_;
    bool store_smallmolecule_reliability_;
    bool store_protein_uri_;
    bool store_peptide_uri_;
    bool store_psm_uri_;
    bool store_smallmolecule_uri_;
    bool store_protein_goterms_;
    bool store_nucleic_acid_reliability_;
    bool store_oligonucleotide_reliability_;
    bool store_osm_reliability_;
    bool store_nucleic_acid_uri_;
    bool store_oligonucleotide_uri_;
    bool store_osm_uri_;
    bool store_nucleic_acid_goterms_;

    void generateMzTabMetaDataSection_(const MzTabMetaData& map, StringList& sl) const;

    /// Needs a reference row to get the collected optional columns from the MetaValues
    /// TODO refactor this behaviour by e.g. storing it in the MzTab object
    String generateMzTabProteinHeader_(const MzTabProteinSectionRow& reference_row,
        const Size n_best_search_engine_scores,
        const std::vector<String>& optional_columns,
        const MzTabMetaData& meta,
        size_t& n_columns) const;

    String generateMzTabSectionRow_(const MzTabProteinSectionRow& row, const std::vector<String>& optional_columns, const MzTabMetaData& meta, size_t& n_columns) const;

    String generateMzTabPeptideHeader_(Size search_ms_runs, Size n_best_search_engine_scores, Size n_search_engine_score, Size assays, Size study_variables, const std::vector<String>& optional_columns, size_t& n_columns) const;

    String generateMzTabSectionRow_(const MzTabPeptideSectionRow& row, const std::vector<String>& optional_columns, const MzTabMetaData& meta, size_t& n_columns) const;

    String generateMzTabPSMHeader_(Size n_search_engine_scores, const std::vector<String>& optional_columns, size_t& n_columns) const;

    String generateMzTabSectionRow_(const MzTabPSMSectionRow& row, const std::vector<String>& optional_columns, const MzTabMetaData& meta, size_t& n_columns) const;

    String generateMzTabSmallMoleculeHeader_(Size search_ms_runs, Size n_best_search_engine_scores, Size n_search_engine_score, Size assays, Size study_variables, const std::vector<String>& optional_columns, size_t& n_columns) const;

    String generateMzTabSectionRow_(const MzTabSmallMoleculeSectionRow& row, const std::vector<String>& optional_columns, const MzTabMetaData& meta, size_t& n_columns) const;

    String generateMzTabNucleicAcidHeader_(Size search_ms_runs, Size n_best_search_engine_scores, Size n_search_engine_scores, const std::vector<String>& optional_columns, size_t& n_columns) const;

    String generateMzTabSectionRow_(const MzTabNucleicAcidSectionRow& row, const std::vector<String>& optional_columns, const MzTabMetaData& meta, size_t& n_columns) const;

    String generateMzTabOligonucleotideHeader_(Size search_ms_runs, Size n_best_search_engine_scores, Size n_search_engine_score, const std::vector<String>& optional_columns, size_t& n_columns) const;

    String generateMzTabSectionRow_(const MzTabOligonucleotideSectionRow& row, const std::vector<String>& optional_columns, const MzTabMetaData& meta, size_t& n_columns) const;

    String generateMzTabOSMHeader_(Size n_search_engine_scores, const std::vector<String>& optional_columns, size_t& n_columns) const;

    String generateMzTabSectionRow_(const MzTabOSMSectionRow& row, const std::vector<String>& optional_columns, const MzTabMetaData& meta, size_t& n_columns) const;

    /// Generate an mzTab section comprising multiple rows of the same type and perform sanity check
    template <typename SectionRow> void generateMzTabSection_(const std::vector<SectionRow>& rows, const std::vector<String>& optional_columns, const MzTabMetaData& meta, StringList& output, size_t n_header_columns) const
    {
      output.reserve(output.size() + rows.size() + 1);
      for (const auto& row : rows)
      {
        size_t n_section_columns = 0;
        output.push_back(generateMzTabSectionRow_(row, optional_columns, meta, n_section_columns));
        if (n_header_columns != n_section_columns)  throw Exception::Postcondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Header and content differs in columns. Please report this bug to the OpenMS developers.");
      }
    }

    // auxiliary functions

    /// Helper function for "generateMzTabSectionRow_" functions
    static void addOptionalColumnsToSectionRow_(const std::vector<String>& column_names, const std::vector<MzTabOptionalColumnEntry>& column_entries, StringList& output);

    // extract two integers from string (e.g. search_engine_score[1]_ms_run[2] -> 1,2)
    static std::pair<int, int> extractIndexPairsFromBrackets_(const String& s);

    static void sortPSM_(std::vector<PeptideIdentification>::iterator begin, std::vector<PeptideIdentification>::iterator end);

    static void keepFirstPSM_(std::vector<PeptideIdentification>::iterator begin, std::vector<PeptideIdentification>::iterator end);

    /// Extract protein and peptide identifications for each run. maps are assumed empty.
    static void partitionIntoRuns_(const std::vector<PeptideIdentification>& pep_ids,
                                   const std::vector<ProteinIdentification>& pro_ids,
                                   std::map<String, std::vector<PeptideIdentification> >& map_run_to_pepids,
                                   std::map<String, std::vector<ProteinIdentification> >& map_run_to_proids
                                   );


    /// create links from protein to peptides
    static void createProteinToPeptideLinks_(const std::map<String, std::vector<PeptideIdentification> >& map_run_to_pepids, MapAccPepType& map_run_accession_to_pephits);

    /// Extracts, if possible a unique protein accession for a peptide hit in mzTab format. Otherwise NA is returned
    static String extractProteinAccession_(const PeptideHit& peptide_hit);

    /// Extracts, modifications and positions of a peptide hit in mzTab format
    static String extractPeptideModifications_(const PeptideHit& peptide_hit);

    /// Map search engine identifier to CV, param etc.
    static String mapSearchEngineToCvParam_(const String& openms_search_engine_name);

    static String mapSearchEngineScoreToCvParam_(const String& openms_search_engine_name, double score, String score_type);

    static String extractNumPeptides_(const String& common_identifier, const String& protein_accession,
                                      const MapAccPepType& map_run_accession_to_peptides);

    // mzTab definition of distinct
    static String extractNumPeptidesDistinct_(String common_identifier, String protein_accession,
                                              const MapAccPepType& map_run_accession_to_peptides);

    // same as distinct but additional constraint of uniqueness (=maps to exactly one Protein)
    static String extractNumPeptidesUnambiguous_(String common_identifier, String protein_accession,
                                                 const MapAccPepType& map_run_accession_to_peptides);

    static std::map<String, Size> extractNumberOfSubSamples_(const std::map<String, std::vector<ProteinIdentification> >& map_run_to_proids);

    static void writePeptideHeader_(SVOutStream& output, std::map<String, Size> n_sub_samples);

    static void writeProteinHeader_(SVOutStream& output, std::map<String, Size> n_sub_samples);

    static void writeProteinData_(SVOutStream& output,
                                  const ProteinIdentification& prot_id,
                                  Size run_count,
                                  String input_filename,
                                  bool has_coverage,
                                  const MapAccPepType& map_run_accession_to_peptides,
                                  const std::map<String, Size>& map_run_to_num_sub
                                  );

  private:
    friend class MzTabMFile;
  };

} // namespace OpenMS

