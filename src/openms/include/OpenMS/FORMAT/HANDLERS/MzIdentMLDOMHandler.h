// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer$
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/METADATA/CVTermList.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMDocumentType.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMNodeIterator.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMText.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>
#include <xercesc/framework/psvi/XSValue.hpp>
#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/XMLUni.hpp>

#include <list>
#include <map>
#include <stdexcept>
#include <string>
#include <vector>

// Error codes
//enum {
//   ERROR_ARGS = 1,
//   ERROR_XERCES_INIT,
//   ERROR_PARSE,
//   ERROR_EMPTY_DOCUMENT
//};

namespace OpenMS
{
  class ProgressLogger;

  namespace Internal
  {
    /**
        @brief XML DOM handler for MzIdentMLFile

        In read-mode, this class will parse an MzIdentML XML file and append the input
        identifications to the provided PeptideIdentifications and ProteinIdentifications.

        @note Do not use this class. It is only needed in MzIdentMLFile.
        @note DOM and STREAM handler for MzIdentML have the same interface for legacy id structures.
        @note Only upon destruction of this class it can be guaranteed that all
        data has been appended to the appropriate containers. Do not try to access the data before that.
    */
    class OPENMS_DLLAPI MzIdentMLDOMHandler
    {
public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler for internal identification structures
      MzIdentMLDOMHandler(const std::vector<ProteinIdentification>& pro_id, const std::vector<PeptideIdentification>& pep_id, const String& version, const ProgressLogger& logger);

      /// Constructor for a read-only handler for internal identification structures
      MzIdentMLDOMHandler(std::vector<ProteinIdentification>& pro_id, std::vector<PeptideIdentification>& pep_id, const String& version, const ProgressLogger& logger);

      /// Destructor
      virtual ~MzIdentMLDOMHandler();
      //@}

      /// Provides the functionality of reading a mzid with a handler object
      void readMzIdentMLFile(const std::string& mzid_file);
      /// Provides the functionality to write a mzid with a handler object
      void writeMzIdentMLFile(const std::string& mzid_file);

protected:
      /// Progress logger
      const ProgressLogger& logger_;

      ///Controlled vocabulary (psi-ms from OpenMS/share/OpenMS/CV/psi-ms.obo)
      ControlledVocabulary cv_;
      ///Controlled vocabulary for modifications (unimod from OpenMS/share/OpenMS/CV/unimod.obo)
      ControlledVocabulary unimod_;

      ///Internal +w Identification Item for proteins
      std::vector<ProteinIdentification>* pro_id_ = nullptr;
      ///Internal +w Identification Item for peptides
      std::vector<PeptideIdentification>* pep_id_ = nullptr;

      ///Internal -w Identification Item for proteins
      const std::vector<ProteinIdentification>* cpro_id_ = nullptr;
      ///Internal -w Identification Item for peptides
      const std::vector<PeptideIdentification>* cpep_id_ = nullptr;

      ///Internal version keeping
      const String schema_version_;

      /// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
      ControlledVocabulary::CVTerm getChildWithName_(const String& parent_accession, const String& name) const;

      /**@name Helper functions to build the internal id structures from the DOM tree */
      //@{

      /// First: CVparams, Second: userParams (independent of each other)
      std::pair<CVTermList, std::map<String, DataValue> > parseParamGroup_(xercesc::DOMNodeList* paramGroup);
      CVTerm parseCvParam_(xercesc::DOMElement* param);
      std::pair<String, DataValue> parseUserParam_(xercesc::DOMElement* param);
      void parseAnalysisSoftwareList_(xercesc::DOMNodeList* analysisSoftwareElements);
      void parseDBSequenceElements_(xercesc::DOMNodeList* dbSequenceElements);
      void parsePeptideElements_(xercesc::DOMNodeList* peptideElements);
      //AASequence parsePeptideSiblings_(xercesc::DOMNodeList* peptideSiblings);
      AASequence parsePeptideSiblings_(xercesc::DOMElement* peptide);
      void parsePeptideEvidenceElements_(xercesc::DOMNodeList* peptideEvidenceElements);
      void parseSpectrumIdentificationElements_(xercesc::DOMNodeList* spectrumIdentificationElements);
      void parseSpectrumIdentificationProtocolElements_(xercesc::DOMNodeList* spectrumIdentificationProtocolElements);
      void parseInputElements_(xercesc::DOMNodeList* inputElements);
      void parseSpectrumIdentificationListElements_(xercesc::DOMNodeList* spectrumIdentificationListElements);
      void parseSpectrumIdentificationItemSetXLMS(std::set<String>::const_iterator set_it, std::multimap<String, int> xl_val_map, xercesc::DOMElement* element_res, const String& spectrumID);
      void parseSpectrumIdentificationItemElement_(xercesc::DOMElement* spectrumIdentificationItemElement, PeptideIdentification& spectrum_identification, String& spectrumIdentificationList_ref);
      void parseProteinDetectionHypothesisElement_(xercesc::DOMElement* proteinDetectionHypothesisElement, ProteinIdentification& protein_identification);
      void parseProteinAmbiguityGroupElement_(xercesc::DOMElement* proteinAmbiguityGroupElement, ProteinIdentification& protein_identification);
      void parseProteinDetectionListElements_(xercesc::DOMNodeList* proteinDetectionListElements);
      static ProteinIdentification::SearchParameters findSearchParameters_(std::pair<CVTermList, std::map<String, DataValue> > as_params);
      //@}

      /**@name Helper functions to build a DOM tree from the internal id structures*/
      void buildCvList_(xercesc::DOMElement* cvElements);
      void buildAnalysisSoftwareList_(xercesc::DOMElement* analysisSoftwareElements);
      void buildSequenceCollection_(xercesc::DOMElement* sequenceCollectionElements);
      void buildAnalysisCollection_(xercesc::DOMElement* analysisCollectionElements);
      void buildAnalysisProtocolCollection_(xercesc::DOMElement* protocolElements);
      void buildInputDataCollection_(xercesc::DOMElement* inputElements);
      void buildEnclosedCV_(xercesc::DOMElement* parentElement, const String& encel, const String& acc, const String& name, const String& cvref);
      void buildAnalysisDataCollection_(xercesc::DOMElement* analysisElements);
      //@}


private:
      MzIdentMLDOMHandler();
      MzIdentMLDOMHandler(const MzIdentMLDOMHandler& rhs);
      MzIdentMLDOMHandler& operator=(const MzIdentMLDOMHandler& rhs);

      ///Struct to hold the used analysis software for that file
      struct AnalysisSoftware
      {
        String name;
        String version;
      };
      ///Struct to hold the PeptideEvidence information
      struct PeptideEvidence
      {
        int start;
        int stop;
        char pre;
        char post;
        bool idec;
      };
      ///Struct to hold the information from the DBSequence xml tag
      struct DBSequence
      {
        String sequence;
        String database_ref;
        String accession;
        CVTermList cvs;
      };
      ///Struct to hold the information from the SpectrumIdentification xml tag
      struct SpectrumIdentification
      {
        String spectra_data_ref;
        String search_database_ref;
        String spectrum_identification_protocol_ref;
        String spectrum_identification_list_ref;
      };
      ///Struct to hold the information from the ModificationParam xml tag
      struct ModificationParam
      {
        String fixed_mod;
        long double mass_delta;
        String residues;
        CVTermList modification_param_cvs;
        CVTermList specificities;
      };
      ///Struct to hold the information from the SpectrumIdentificationProtocol xml tag
      struct SpectrumIdentificationProtocol
      {
        CVTerm searchtype;
        String enzyme;
        CVTermList parameter_cvs;
        std::map<String, DataValue> parameter_ups;
//        std::vector<ModificationParam> modification_parameter;
        CVTermList modification_parameter;
        long double precursor_tolerance;
        long double fragment_tolerance;
        CVTermList threshold_cvs;
        std::map<String, DataValue> threshold_ups;
      };
      ///Struct to hold the information from the DatabaseInput xml tag
      struct DatabaseInput
      {
        String name;
        String location;
        String version;
        DateTime date;
      };

      XMLCh* xml_root_tag_ptr_;
      XMLCh* xml_cvparam_tag_ptr_;
      XMLCh* xml_name_attr_ptr_;

      xercesc::XercesDOMParser mzid_parser_;

      std::unique_ptr<XMLHandler> xml_handler_ = nullptr;

      //from AnalysisSoftware
      String search_engine_;
      String search_engine_version_;
      //mapping from AnalysisSoftware
      std::map<String, AnalysisSoftware> as_map_; ///< mapping AnalysisSoftware id -> AnalysisSoftware

      //mapping from DataCollection Inputs
      std::map<String, String> sr_map_; ///< mapping sourcefile id -> sourcefile location
      std::map<String, String> sd_map_; ///< mapping spectradata id -> spectradata location
      std::map<String, DatabaseInput> db_map_; ///< mapping database id -> DatabaseInput

      //mapping from SpectrumIdentification - SpectrumIdentification will be the new IdentificationRuns
      std::map<String, SpectrumIdentification> si_map_; ///< mapping SpectrumIdentification id -> SpectrumIdentification (id refs)
      std::map<String, size_t> si_pro_map_; ///< mapping SpectrumIdentificationList id -> index to ProteinIdentification in pro_id_

      //mapping from SpectrumIdentificationProtocol
      std::map<String, SpectrumIdentificationProtocol> sp_map_; ///< mapping SpectrumIdentificationProtocol id -> SpectrumIdentificationProtocol

      //mapping from SequenceCollection
      std::map<String, AASequence> pep_map_; ///< mapping Peptide id -> Sequence
      std::map<String, PeptideEvidence> pe_ev_map_; ///< mapping PeptideEvidence id -> PeptideEvidence
      std::map<String, String> pv_db_map_; ///< mapping PeptideEvidence id -> DBSequence id
      std::multimap<String, String> p_pv_map_; ///< mapping Peptide id -> PeptideEvidence id, multiple PeptideEvidences can have equivalent Peptides.
      std::map<String, DBSequence> db_sq_map_; ///< mapping DBSequence id -> Sequence

      std::list<std::list<String> > hit_pev_; ///< writing help only

      bool xl_ms_search_; ///< is true when reading a file containing Cross-Linking MS search results
      std::map<String, String> xl_id_donor_map_; ///< mapping Peptide id -> crosslink donor value
      //std::map<String, String> xl_id_acceptor_map_; ///< mapping Peptide id -> crosslink acceptor value
      std::map<String, String> xl_id_acceptor_map_; ///< mapping  peptide id of acceptor peptide -> crosslink acceptor value
      std::map<String, SignedSize> xl_donor_pos_map_; ///< mapping donor value -> cross-link modification location
      std::map<String, SignedSize> xl_acceptor_pos_map_; ///< mapping acceptor value -> cross-link modification location
      std::map<String, double> xl_mass_map_; ///< mapping Peptide id -> cross-link mass
      std::map<String, String> xl_mod_map_; ///< mapping peptide id -> cross-linking reagent name

    };
  } // namespace Internal
} // namespace OpenMS

