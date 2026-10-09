// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>

#include <functional>
#include <iosfwd>
#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    class FeatureXMLHandler;
    class ConsensusXMLHandler;
  }
  /**
    @brief Used to load and store idXML files

    This class is used to load and store documents that implement
    the schema of idXML files.

    A documented schema for this format can be found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

    One file can contain several ProteinIdentification runs. Each run consists of peptide hits stored in
    PeptideIdentification and (optional) protein hits stored in Identification. Peptide and protein
    hits are connected via a string identifier. We use the search engine and the date as identifier.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI IdXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile,
    public ProgressLogger
  {
public:
    // both ConsensusXMLFile and FeatureXMLFile use some protected IdXML helper functions to parse identifications without code duplication
    friend class Internal::ConsensusXMLHandler;
    friend class Internal::FeatureXMLHandler;

    /// Constructor
    IdXMLFile();

    /**
        @brief Loads the identifications of an idXML file without identifier

        The information is read in and the information is stored in the
        corresponding variables

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const std::string& filename, std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids);

    /**
        @brief Loads the identifications of an idXML file

        The information is read in and the information is stored in the
        corresponding variables

        @exception Exception::FileNotFound is thrown if the file could not be opened
        @exception Exception::ParseError is thrown if an error occurs during parsing
    */
    void load(const std::string& filename, std::vector<ProteinIdentification>& protein_ids, PeptideIdentificationList& peptide_ids, std::string& document_id);

    /**
        @brief Stores the data in an idXML file

        The data is read in and stored in the file 'filename'. PeptideHits are sorted by score.
        Note that ranks are not stored and need to be reassigned after loading.

        @exception Exception::UnableToCreateFile is thrown if the file could not be created, or if writing it failed
                   (e.g. disk full or an I/O error).

        If store() fails, the file is not removed: a partial file may remain. An exception of OpenMS (e.g.
        Exception::ConversionError for a meta value that cannot be written) keeps its type, and its message names the
        file and says so (unless memory runs out while the note is added). Any other exception (e.g. std::bad_alloc) is
        raised unchanged.
    */
    void store(const std::string& filename, const std::vector<ProteinIdentification>& protein_ids, const PeptideIdentificationList& peptide_ids, const std::string& document_id = "");


protected:
    // Docu in base class
    void onEndElement(const char16_t* qname) override;

    // Docu in base class
    void onStartElement(const char16_t* qname, const Internal::XMLAttributes& attributes) override;

    /// Add data from ProteinGroups to a MetaInfoInterface
    /// Since it can be used during load and store, it needs to take a param for the current mode (LOAD/STORE)
    /// to throw appropriate warnings/errors
    void addProteinGroups_(MetaInfoInterface& meta, const std::vector<ProteinIdentification::ProteinGroup>& groups,
                           const std::string& group_name, const std::unordered_map<std::string, UInt>& accession_to_id, XMLHandler::ActionMode mode);

    /// Read and store ProteinGroup data
    void getProteinGroups_(std::vector<ProteinIdentification::ProteinGroup>& groups, const std::string& group_name);

    /**
      * Helper function to create the XML string for the amino acids before and after the peptide position in a protein.
      * Can be reused by e.g. ConsensusXML, FeatureXML to write PeptideHit elements  
      */
    static std::ostream& createFlankingAAXMLString_(const std::vector<PeptideEvidence> & pes, std::ostream& os);

    /**
      * Helper function to create the XML string for the position of the peptide in a protein.
      * Can be reused by e.g. ConsensusXML, FeatureXML to write PeptideHit elements  
      */
    static std::ostream& createPositionXMLString_(const std::vector<PeptideEvidence> & pes, std::ostream& os);


    /**
      * Helper function to write out fragment annotations as user param fragment_annotation
      */  
    static void writeFragmentAnnotations_(const std::string & tag_name, std::ostream & os, 
                                          const std::vector<PeptideHit::PeakAnnotation>& annotations, UInt indent);

    /**
      * Helper function to parse fragment annotations from string
      */  
    static void parseFragmentAnnotation_(const std::string& s, std::vector<PeptideHit::PeakAnnotation> & annotations);

    /**
      @brief Formats one block of the parallel peptide identification writer of store(): calls @p write on @p block_os.

      Exceptions are enabled on @p block_os (badbit and failbit) first. A failure of its stream buffer, such as a
      std::bad_alloc while a std::stringbuf grows, then propagates (the original exception is rethrown) instead of only
      setting badbit and leaving a silently truncated block that would be written as if it were complete.

      @exception std::exception whatever @p write or the stream buffer throws; std::ios_base::failure if the stream
                 fails without an exception of its own
    */
    static void formatBlock_(std::ostream& block_os, const std::function<void(std::ostream&)>& write);


    /// @name members for loading data
    //@{
    /// Pointer to fill in protein identifications
    std::vector<ProteinIdentification>* prot_ids_;
    /// Pointer to fill in peptide identifications
    PeptideIdentificationList* pep_ids_;
    /// Pointer to last read object with MetaInfoInterface
    MetaInfoInterface* last_meta_;
    /// Search parameters map (key is the "id")
    std::map<std::string, ProteinIdentification::SearchParameters> parameters_;
    /// Temporary search parameters variable
    ProteinIdentification::SearchParameters param_;
    /// Temporary id
    std::string id_;
    /// Temporary protein ProteinIdentification
    ProteinIdentification prot_id_;
    /// Temporary peptide ProteinIdentification
    PeptideIdentification pep_id_;
    /// Temporary protein hit
    ProteinHit prot_hit_;
    /// Temporary peptide hit
    PeptideHit pep_hit_;
    /// Temporary analysis result instance
    PeptideHit::PepXMLAnalysisResult current_analysis_result_;
    /// Temporary peptide evidences
    std::vector<PeptideEvidence> peptide_evidences_;
    /// Map from protein id to accession
    std::unordered_map<std::string, std::string> proteinid_to_accession_;
    /// Document identifier
    std::string* document_id_;
    /// true if a prot id is contained in the current run
    bool prot_id_in_run_;
    //@}
  };

} // namespace OpenMS

