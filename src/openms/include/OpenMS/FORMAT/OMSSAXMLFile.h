// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/XMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/ModificationDefinitionsSet.h>

#include <vector>

namespace OpenMS
{
  class String;
  class ModificationsDB;

  /**
    @brief Used to load OMSSAXML files

    This class is used to load documents that implement
    the schema of OMSSAXML files.

    @ingroup FileIO
  */
  class OPENMS_DLLAPI OMSSAXMLFile :
    protected Internal::XMLHandler,
    public Internal::XMLFile
  {
public:

    /// Default constructor
    OMSSAXMLFile();

    /// Destructor
    ~OMSSAXMLFile() override;
    /**
      @brief loads data from a OMSSAXML file

      @param filename The file to be loaded
      @param protein_identification Protein identifications belonging to the whole experiment
      @param id_data The identifications with m/z and RT
      @param load_proteins If this flag is set to false, the protein identifications are not loaded
      @param load_empty_hits Many spectra will not return a hit. Report empty peptide identifications?

      This class serves to read in a OMSSAXML file. The information can be
      retrieved via the load function.

          @exception FileNotFound
          @exception ParseError

      @ingroup FileIO
    */
    void load(const String& filename,
              ProteinIdentification& protein_identification,
              std::vector<PeptideIdentification>& id_data,
              bool load_proteins = true,
              bool load_empty_hits = true);

    /// sets the valid modifications
    void setModificationDefinitionsSet(const ModificationDefinitionsSet& rhs);

protected:
    // Docu in base class
    void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

    // Docu in base class
    void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

    // Docu in base class
    void characters(const XMLCh* const chars, const XMLSize_t /*length*/) override;

private:

    OMSSAXMLFile(const OMSSAXMLFile& rhs);

    OMSSAXMLFile& operator=(const OMSSAXMLFile& rhs);

    /// reads the mapping file needed for modifications
    void readMappingFile_();

    /// the identifications (storing the peptide hits)
    std::vector<PeptideIdentification>* peptide_identifications_;

    ProteinHit actual_protein_hit_;

    PeptideHit actual_peptide_hit_;

    PeptideEvidence actual_peptide_evidence_;

    std::vector<PeptideEvidence> actual_peptide_evidences_;

    PeptideIdentification actual_peptide_id_;

    ProteinIdentification actual_protein_id_;

    String tag_;

    /// site of the actual modification (simple position in the peptide)
    UInt actual_mod_site_;

    /// type of the modification
    String actual_mod_type_;

    /// modifications of the peptide defined by site and type
    std::vector<std::pair<UInt, String> > modifications_;

    /// should protein hits be read from the file?
    bool load_proteins_;

    /// should empty peptide identifications be loaded or skipped?
    bool load_empty_hits_;

    /// modifications mapping file from OMSSA mod num to UniMod accession
    std::map<UInt, std::vector<const ResidueModification*> > mods_map_;

    /// modification mapping reverse, from the modification to the mod_num
    std::map<String, UInt> mods_to_num_;

    /// modification definitions set of the search, needed to annotate fixed modifications
    ModificationDefinitionsSet mod_def_set_;
  };

} // namespace OpenMS

