// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer, Andreas Bertsch $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/METADATA/IdentificationHit.h>
#include <OpenMS/METADATA/Identification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/DigestionEnzymeProtein.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>
#include <map>

namespace OpenMS
{
  class ProgressLogger;

  namespace Internal
  {

    /**
        @brief XML STREAM handler for MzIdentMLFile

        In read-mode, this class will parse an MzIdentML XML file and append the input
        identifications to the provided PeptideIdentifications and ProteinIdentifications.

        @note Do not use this class. It is only needed in MzIdentMLFile.
        @note DOM and STREAM handler for MzIdentML have the same interface for legacy id structures.
    */
    class OPENMS_DLLAPI MzIdentMLHandler :
      public XMLHandler
    {
public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler
      MzIdentMLHandler(const Identification& id, const String& filename, const String& version, const ProgressLogger& logger);
      /// Constructor for a write-only handler for internal identification structures
      MzIdentMLHandler(const std::vector<ProteinIdentification>& pro_id, const std::vector<PeptideIdentification>& pep_id, const String& filename, const String& version, const ProgressLogger& logger);

      /// Constructor for a read-only handler
      MzIdentMLHandler(Identification& id, const String& filename, const String& version, const ProgressLogger& logger);
      /// Constructor for a read-only handler for internal identification structures
      MzIdentMLHandler(std::vector<ProteinIdentification>& pro_id, std::vector<PeptideIdentification>& pep_id, const String& filename, const String& version, const ProgressLogger& logger);

      /// Destructor
      ~MzIdentMLHandler() override;
      //@}


      // Docu in base class
      void endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname) override;

      // Docu in base class
      void startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes) override;

      // Docu in base class
      void characters(const XMLCh* const chars, const XMLSize_t length) override;

      //Docu in base class
      void writeTo(std::ostream& os) override;

protected:
      /// Progress logger
      const ProgressLogger& logger_;

      ///Controlled vocabulary (psi-ms from OpenMS/share/OpenMS/CV/psi-ms.obo)
      ControlledVocabulary cv_;
      ///Controlled vocabulary for modifications (unimod from OpenMS/share/OpenMS/CV/unimod.obo)
      ControlledVocabulary unimod_;

      //~ PeakMap* ms_exp_;

      ///XML tag parse element
      String tag_;

      ///Identification Item
      Identification* id_;
      ///internal Identification Item for proteins
      std::vector<ProteinIdentification>* pro_id_;
      ///Identification Item for peptides
      std::vector<PeptideIdentification>* pep_id_;

      const Identification* cid_;
      const std::vector<ProteinIdentification>* cpro_id_;
      const std::vector<PeptideIdentification>* cpep_id_;

      ///SpectrumIdentification Item
      SpectrumIdentification current_spectrum_id_;

      ///IdentificationHit Item
      IdentificationHit current_id_hit_;

      /// Handles CV terms
      void handleCVParam_(const String& parent_parent_tag, const String& parent_tag, const String& accession, /* const String& name, */ /* const String& value, */ const xercesc::Attributes& attributes, const String& cv_ref /* ,  const String& unit_accession="" */);

      /// Handles user terms
      void handleUserParam_(const String& parent_parent_tag, const String& parent_tag, const String& name, const String& type, const String& value);

      /// Writes user terms
      void writeMetaInfos_(String& s, const MetaInfoInterface& meta, UInt indent) const;

      /// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
      ControlledVocabulary::CVTerm getChildWithName_(const String& parent_accession, const String& name) const;

      /// Helper method that writes a source file
      //void writeSourceFile_(std::ostream& os, const String& id, const SourceFile& software);

      /// Helper method that writes the Enzymes
      void writeEnzyme_(String& s, const DigestionEnzymeProtein& enzy, UInt miss, UInt indent) const;

      /// Helper method that writes the modification search params (fixed or variable)
      void writeModParam_(String& s, const std::vector<String>& mod_names, bool fixed, UInt indent) const;

      /// Helper method that writes the FragmentAnnotations section of a spectrum identification
      void writeFragmentAnnotations_(String& s, const std::vector<PeptideHit::PeakAnnotation>& annotations, UInt indent, bool is_ppxl) const;

      /// Convenience method to remove the [] from OpenMS internal file uri representation
      String trimOpenMSfileURI(const String& file) const;

      /// Abstraction of PeptideHit loop for most PeptideHits
      void writePeptideHit(const PeptideHit& hit,
                                std::vector<PeptideIdentification>::const_iterator& it,
                                std::map<String, String>& pep_ids,
                                const String& cv_ns, std::set<String>& sen_set,
                                std::map<String, String>& sen_ids,
                                std::map<String, std::vector<String> >& pep_evis,
                                std::map<String, double>& pp_identifier_2_thresh,
                                String& sidres);

      /// Abstraction of PeptideHit loop for XL-MS data from OpenPepXL
      void writeXLMSPeptideHit(const PeptideHit& hit,
                                std::vector<PeptideIdentification>::const_iterator& it,
                                const String& ppxl_linkid, std::map<String, String>& pep_ids,
                                const String& cv_ns, std::set<String>& sen_set,
                                std::map<String, String>& sen_ids,
                                std::map<String, std::vector<String> >& pep_evis,
                                std::map<String, double>& pp_identifier_2_thresh,
                                double ppxl_crosslink_mass,
                                std::map<String, String>& ppxl_specref_2_element,
                                String& sid, bool alpha_peptide);

private:
      MzIdentMLHandler();
      MzIdentMLHandler(const MzIdentMLHandler& rhs);
      MzIdentMLHandler& operator=(const MzIdentMLHandler& rhs);
      std::map<String, AASequence> pep_sequences_;
      std::map<String, String> pp_identifier_2_sil_; ///< mapping peptide/proteinidentification identifier_ to spectrumidentificationlist
      std::map<String, String> sil_2_sdb_; ///< mapping spectrumidentificationlist to the search data bases
      std::map<String, String> sil_2_sdat_; ///< mapping spectrumidentificationlist to the search input
      std::map<String, String> ph_2_sdat_; ///< mapping identification runs (mapping PeptideIdentifications and ProteinIdentifications via .getIdentifier()) to spectra data
      std::map<String, String> sil_2_sip_; ///< mapping spectrumidentificationlist to the search protocol (where the params are at)
      AASequence actual_peptide_;
      Int current_mod_location_;
      ProteinHit actual_protein_;

    };
  } // namespace Internal
} // namespace OpenMS
