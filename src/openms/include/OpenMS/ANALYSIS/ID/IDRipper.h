// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Immanuel Luhn, Leon Kuchenbecker$
// --------------------------------------------------------------------------
#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <unordered_map>


namespace OpenMS
{
  /**
    @brief Ripping protein/peptide identification according their file origin.

    Helper class, which is used by @ref TOPP_ProteinQuantifier. See there for further documentation.

    @htmlinclude OpenMS_IDRipper.parameters

    @ingroup Analysis_ID
  */
  class OPENMS_DLLAPI IDRipper :
    public DefaultParamHandler
  {
public:
    /// Possible input file encodings for the origin as used by different versions of IDMerger
    enum OriginAnnotationFormat { FILE_ORIGIN = 0, MAP_INDEX = 1, ID_MERGE_INDEX = 2, UNKNOWN_OAF = 3, SIZE_OF_ORIGIN_ANNOTATION_FORMAT = 4 };

    /// String representations for the OriginAnnotationFormat enum
    static const std::array<std::string, SIZE_OF_ORIGIN_ANNOTATION_FORMAT> names_of_OriginAnnotationFormat;

    /// Represents a set of IdentificationRuns
    struct OPENMS_DLLAPI IdentificationRuns
    {
        /// Maps a unique index to every IdentificationRun string representation (getIdentifier()).
        std::map<String, UInt> index_map;
        /// Maps the list of spectra data elements to every IdentificationRun index.
        std::vector<StringList> spectra_data;

        /// Generates a new IdentificationRuns object from a vector of ProteinIdentification objects.
        IdentificationRuns(const std::vector<ProteinIdentification>& prot_ids);
    };

    /// Identifies an IDRipper output file
    struct OPENMS_DLLAPI RipFileIdentifier
    {
        /// The numerical index of the source IdentificationRun
        UInt ident_run_idx{};
        /// The numerical index of the source file_origin / spectra_data element
        UInt file_origin_idx{};
        /// The output basename derived from the file_origin / spectra_data element
        String out_basename;
        /// The full length origin read from the file_origin / spectra_data element
        String origin_fullname;

        /// Constructs a new RipFileIdentifier object
        RipFileIdentifier(const IDRipper::IdentificationRuns& id_runs, 
          const PeptideIdentification& pep_id, 
          const std::map<String, UInt>& file_origin_map, 
          const IDRipper::OriginAnnotationFormat origin_annotation_fmt, 
          bool split_ident_runs);

        /// Get identification run index
        UInt getIdentRunIdx() const;

        /// Get file origin index
        UInt getFileOriginIdx() const;

        /// Get origin full name
        const String & getOriginFullname() const;

        /// Get output base name
        const String & getOutputBasename() const;
    };

    /// Provides a 'less' operation for RipFileIdentifiers that ignores the out_basename and origin_fullname members
    struct RipFileIdentifierIdxComparator
    {
        bool operator()(const RipFileIdentifier& left, const RipFileIdentifier& right) const;
    };

    /// Represents the content of an IDRipper output file
    struct OPENMS_DLLAPI RipFileContent
    {
        /// Protein identifications
        std::vector<ProteinIdentification> prot_idents;
        /// Peptide identifications
        std::vector<PeptideIdentification> pep_idents;
        /// Constructs a new RipFileContent object
        RipFileContent(const std::vector<ProteinIdentification>& prot_idents, const std::vector<PeptideIdentification>& pep_idents)
            : prot_idents(prot_idents), pep_idents(pep_idents) {}
        /// Get protein identifications
        const std::vector<ProteinIdentification> & getProteinIdentifications();
        /// Get peptide identifications
        const std::vector<PeptideIdentification> & getPeptideIdentifications();
    };

    /// Represents the result of an IDRipper process, a map assigning file content to output file identifiers
    typedef std::map<RipFileIdentifier, RipFileContent, RipFileIdentifierIdxComparator> RipFileMap;

    /// Default constructor
    IDRipper();

    /// Destructor
    ~IDRipper() override;

    /**
      @brief Ripping protein/peptide identification according their file origin

      Iteration over all @p peptides. For each annotated file origin create a map entry and store the
      respective @p peptides and @p proteins.

      @param ripped Contains the protein identification and peptide identification for each file origin annotated in proteins and peptides
      @param proteins Protein identification
      @param peptides Peptide identification annotated with file origin
      @param numeric_filenames If false, deduce output files using basenames of origin annotations. Throws an exception if they are not unique. If true, assemble output files based on numerical IDs only.
      @param split_ident_runs Split identification runs into different files.
    */
    void rip(
            RipFileMap& ripped,
            std::vector<ProteinIdentification>& proteins,
            std::vector<PeptideIdentification>& peptides,
            bool numeric_filenames,
            bool split_ident_runs);

    /**
      @brief Ripping protein/peptide identification according their file origin

      Iteration over all @p peptides. For each annotated file origin create a map entry and store the
      respective @p peptides and @p proteins.

      @param[out] rfis File info relating to @p rfcs
      @param[out] rfcs Contains the protein identification and peptide identification for each file origin annotated in proteins and peptides
      @param proteins Protein identification
      @param peptides Peptide identification annotated with file origin
      @param numeric_filenames If false, deduce output files using basenames of origin annotations. Throws an exception if they are not unique. If true, assemble output files based on numerical IDs only.
      @param split_ident_runs Split identification runs into different files.
    */
    // Autowrap compatible wrapper for rip(RipFileMap,...)
    void rip(
            std::vector<RipFileIdentifier>& rfis,
            std::vector<RipFileContent>& rfcs,
            std::vector<ProteinIdentification>& proteins,
            std::vector<PeptideIdentification>& peptides,
            bool numeric_filenames,
            bool split_ident_runs);

private:
    // Not implemented
    /// Copy constructor
    IDRipper(const IDRipper & rhs);

    // Not implemented
    /// Assignment
    IDRipper & operator=(const IDRipper & rhs);

    /// helper function, detects file origin annotation standard from collections of protein and peptide hits
    OriginAnnotationFormat detectOriginAnnotationFormat_(std::map<String, UInt> & file_origin_map, const std::vector<PeptideIdentification> & peptide_idents);
    /// helper function, extracts all protein hits that match the protein accession
    void getProteinHits_(std::vector<ProteinHit> & result, const std::unordered_map<String, const ProteinHit*> & acc2protein_hits, const std::set<String> & protein_accessions);
    /// helper function, returns the string representation of the peptide hit accession
    std::set<String> getProteinAccessions_(const std::vector<PeptideHit> & peptide_hits);
    /// helper function, returns the index of the protein identification for the given peptide identification based on the same identifier using id_runs as lookup
    int getProteinIdentification_(const PeptideIdentification& pep_ident, const IdentificationRuns& id_runs);
    /// helper function, register a potential output file basename to detect duplicate output basenames
    bool registerBasename_(std::map<String, std::pair<UInt, UInt> >& basename_to_numeric, const IDRipper::RipFileIdentifier& rfi);
    /// helper function, sets the value of mode to new_value and returns true if the old value was identical or unset (-1)
    bool setOriginAnnotationMode_(short& mode, short const new_value);
  };

} // namespace OpenMS
