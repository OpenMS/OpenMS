// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/OPTIONS/PeakFileOptions.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <unordered_map>
#include <map>

namespace OpenMS
{
namespace Internal
{
  /**
    @brief This class provides Input functionality for ConsensusMaps and Output functionality for
    alignments and quantitation.

    This class can be used to load the content of a consensusXML file into a ConsensusMap
    or to save the content of an ConsensusMap object into an XML file.

    A documented schema for this format can be found at https://github.com/OpenMS/OpenMS/tree/develop/share/OpenMS/SCHEMAS

  @todo Take care that unique ids are assigned properly by TOPP tools before calling ConsensusXMLFile::store().  There will be a message on OPENMS_LOG_INFO but we will make no attempt to fix the problem in this class.  (all developers)

    @ingroup FileIO
  */
  class OPENMS_DLLAPI ConsensusXMLHandler :
    public Internal::XMLHandler,
    public ProgressLogger
  {
public:
    ///Constructor
    ConsensusXMLHandler(ConsensusMap& map, const String& filename);
    ConsensusXMLHandler(const ConsensusMap& map, const String& filename);
    ///Destructor
    ~ConsensusXMLHandler() override;

    /// Mutable access to the options for loading/storing
    PeakFileOptions& getOptions();

    void setOptions(const PeakFileOptions&);

    /// Non-mutable access to the options for loading/storing
    const PeakFileOptions& getOptions() const;

    /// Docu in base class XMLHandler::writeTo
    void writeTo(std::ostream& os) override;

protected:

    // Docu in base class
    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

    // Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

    // Docu in base class
    void characters(const XMLCh* const chars, const XMLSize_t length) override;

    /// Writes a peptide identification to a stream (for assigned/unassigned peptide identifications)
    void writePeptideIdentification_(const String& filename, std::ostream& os, const PeptideIdentification& id, const String& tag_name, UInt indentation_level);

    /// Add data from ProteinGroups to a MetaInfoInterface
    /// Since it can be used during load and store, it needs to take a param for the current mode (LOAD/STORE)
    /// to throw appropriate warnings/errors
    void addProteinGroups_(MetaInfoInterface& meta, const std::vector<ProteinIdentification::ProteinGroup>& groups,
                           const String& group_name, const std::unordered_map<std::string, UInt>& accession_to_id,
                           const String& runid, XMLHandler::ActionMode mode);

    /// Read and store ProteinGroup data
    void getProteinGroups_(std::vector<ProteinIdentification::ProteinGroup>& groups, const String& group_name);

    /// Options that can be set
    PeakFileOptions options_;

    ///@name Temporary variables for parsing
    //@{
    ConsensusMap* consensus_map_;
    const ConsensusMap* cconsensus_map_;
    ConsensusFeature act_cons_element_;
    DPosition<2> pos_;
    double it_;
    //@}

    /// Pointer to last read object as a MetaInfoInterface, or null.
    MetaInfoInterface* last_meta_;
    /// Temporary protein ProteinIdentification
    ProteinIdentification prot_id_;
    /// Temporary peptide ProteinIdentification
    PeptideIdentification pep_id_;
    /// Temporary protein hit
    ProteinHit prot_hit_;
    /// Temporary peptide hit
    PeptideHit pep_hit_;
    /// Temporary peptide evidences
    std::vector<PeptideEvidence> peptide_evidences_;
    /// Map from protein id to accession
    std::map<String, String> proteinid_to_accession_;
    /// Map from search identifier concatenated with protein accession to id
    std::unordered_map<std::string, UInt> accession_to_id_;
    /// Map from identification run identifier to file xs:id (for linking peptide identifications to the corresponding run)
    std::map<String, String> identifier_id_;
    /// Map from file xs:id to identification run identifier (for linking peptide identifications to the corresponding run)
    std::map<String, String> id_identifier_;
    /// Temporary search parameters file
    ProteinIdentification::SearchParameters search_param_;

    UInt progress_;
  };
} // namespace Internal
} // namespace OpenMS

