// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/PeptideEvidence.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/SpectrumMetaDataLookup.h>

#include <vector>

namespace OpenMS
{
  namespace Internal
  {
    /**
      @brief Handler that is used for parsing MascotXML data
    */
    class OPENMS_DLLAPI MascotXMLHandler :
      public XMLHandler
    {
public:
      /// Constructor
      MascotXMLHandler(ProteinIdentification& protein_identification,
                       std::vector<PeptideIdentification>& identifications,
                       const String& filename,
                       std::map<String, std::vector<AASequence> >& peptides,
                       const SpectrumMetaDataLookup& lookup);

      /// Destructor
      ~MascotXMLHandler() override;

      // Docu in base class
      void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

      // Docu in base class
      void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

      // Docu in base class
      void characters(const XMLCh* const chars, const XMLSize_t /*length*/) override;
      
      /// Split modification search parameter if for more than one amino acid specified e.g. Phospho (ST)
      static std::vector<String> splitModificationBySpecifiedAA(const String& mod);

private:

      ProteinIdentification& protein_identification_; ///< the protein identifications
      std::vector<PeptideIdentification>& id_data_; ///< the identifications (storing the peptide hits)
      ProteinHit actual_protein_hit_;
      PeptideHit actual_peptide_hit_;
      PeptideEvidence actual_peptide_evidence_;
      UInt peptide_identification_index_;
      String tag_;
      DateTime date_;
      String date_time_string_;
      UInt actual_query_;
      ProteinIdentification::SearchParameters search_parameters_;
      String identifier_;
      String actual_title_;
      std::map<String, std::vector<AASequence> >& modified_peptides_;

      StringList tags_open_; ///< tracking the current XML tree
      String character_buffer_; ///< filled by MascotXMLHandler::characters
      String major_version_;
      String minor_version_;
      
      // list of modifications, which cannot be set as fixed and needs
      // to be removed, because added from mascot as variable modification
      std::vector<String> remove_fixed_mods_;

      /// Helper object for looking up RT information
      const SpectrumMetaDataLookup& lookup_;

      /// Error for missing RT information already reported?
      bool no_rt_error_;
    };

  } // namespace Internal
} // namespace OpenMS

