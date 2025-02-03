// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>
#include <stack>

namespace OpenMS
{
  class String;
  class ProteinIdentification;

  /**
    @brief Used to load XTandemXML files

    This class is used to load documents that implement
    the schema of XTandemXML files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI XTandemXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Default constructor
    XTandemXMLFile();

    /// Destructor
    ~XTandemXMLFile() override;
    /**
      @brief loads data from an X! Tandem XML file

      @param filename the file to be loaded
      @param protein_identification protein identifications belonging to the whole experiment
      @param id_data the identifications with m/z and RT
      @param mod_def_set Fixed and variable modifications defined for the search. May be extended with additional (X! Tandem default) modifications if those are found in the file.

      This class serves to read in an X! Tandem XML file. The information can be
      retrieved via the load function.

      @ingroup FileIO
    */
    void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data, ModificationDefinitionsSet& mod_def_set);


protected:

    // Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

    // Docu in base class
    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

    // Docu in base class
    void characters(const XMLCh* const chars, const XMLSize_t /*length*/) override;

    XTandemXMLFile(const XTandemXMLFile& rhs);

    XTandemXMLFile& operator=(const XTandemXMLFile& rhs);

private:

    ProteinIdentification* protein_identification_;

    // true during "note" element containing protein accession
    bool is_protein_note_;

    // true during "note" element containing spectrum ID
    bool is_spectrum_note_;

    // true after non-new protein entries, so that with the next "protein note" the
    // accession will not be updated again
    bool skip_protein_acc_update_;

    // peptide hits per spectrum
    std::map<UInt, std::vector<PeptideHit> > peptide_hits_;

    // protein hits
    std::vector<ProteinHit> protein_hits_;

    // protein unique IDs (assigned by X! Tandem), to keep track of which proteins were already seen
    std::set<UInt> protein_uids_;

    // accession of the current protein
    String current_protein_;

    // charge of current peptide
    Int current_charge_;

    // X! Tandem ID of current peptide
    UInt current_id_;

    // tag
    String tag_;

    // start position of current peptide in protein sequence
    UInt current_start_;

    // stop position of current peptide in protein sequence
    UInt current_stop_;

    // previous peptide sequence
    String previous_seq_;

    // mapping from X! Tandem ID to spectrum ID
    std::map<UInt, String> spectrum_ids_;

    // modification definitions
    ModificationDefinitionsSet mod_def_set_;

    // modifications used by X! Tandem by default
    ModificationDefinitionsSet default_nterm_mods_;

    // the possible type attributes of the group tag elements
    enum class GroupType
    {
      MODEL,
      PARAMETERS,
      SUPPORT
    };

    // stack of types of the group elements
    // they can be nested (e.g. a support group in a model group)
    // parsing of child elements sometimes depends on the group type
    std::stack<GroupType> group_type_stack_;
    
  };

} // namespace OpenMS

