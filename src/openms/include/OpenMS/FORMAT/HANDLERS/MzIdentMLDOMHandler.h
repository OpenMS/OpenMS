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
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/METADATA/ProteinIdentification.h>


#include <list>

struct _xmlNode;
typedef struct _xmlNode xmlNode;
#include <map>
#include <set>
#include <stdexcept>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <vector>

namespace OpenMS
{
  class ProgressLogger;

  namespace Internal
  {
    struct MzIdentMLNode;

    /**
        @brief XML handler for MzIdentMLFile

        In read-mode, this class will parse an MzIdentML XML file and append the input
        identifications to the provided PeptideIdentifications and ProteinIdentifications.
        It reads the file with the SAX2 interface of libxml2, one SpectrumIdentificationResult
        at a time, and writes a file through a libxml2 DOM tree.

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
      MzIdentMLDOMHandler(const std::vector<ProteinIdentification>& pro_id, const PeptideIdentificationList& pep_id, const std::string& version, const ProgressLogger& logger);

      /// Constructor for a read-only handler for internal identification structures
      MzIdentMLDOMHandler(std::vector<ProteinIdentification>& pro_id, PeptideIdentificationList& pep_id, const std::string& version, const ProgressLogger& logger);

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
      PeptideIdentificationList* pep_id_ = nullptr;

      ///Internal -w Identification Item for proteins
      const std::vector<ProteinIdentification>* cpro_id_ = nullptr;
      ///Internal -w Identification Item for peptides
      const PeptideIdentificationList* cpep_id_ = nullptr;

      ///Internal version keeping
      const std::string schema_version_;

      /// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
      ControlledVocabulary::CVTerm getChildWithName_(const std::string& parent_accession, const std::string& name) const;

      /// Precompute the CV child-term sets used per PSM (constant across a file); shared by both constructors.
      void initScoreTermCaches_();

      /**@name Helper functions to build the internal id structures from the elements of the file */
      //@{

      /// First: CVparams, Second: userParams (independent of each other)
      std::pair<CVTermList, std::map<std::string, DataValue> > parseParamGroup_(const std::vector<MzIdentMLNode*>& paramGroup);
      CVTerm parseCvParam_(MzIdentMLNode* param);
      std::pair<std::string, DataValue> parseUserParam_(MzIdentMLNode* param);
      void parseAnalysisSoftwareList_(const std::vector<MzIdentMLNode*>& analysisSoftwareElements);
      void parseDBSequenceElements_(const std::vector<MzIdentMLNode*>& dbSequenceElements);
      void parsePeptideElements_(const std::vector<MzIdentMLNode*>& peptideElements);
      //AASequence parsePeptideSiblings_(const std::vector<MzIdentMLNode*>& peptideSiblings);
      AASequence parsePeptideSiblings_(MzIdentMLNode* peptide);
      /// Best-effort inference of a modification's mzIdentML @c location (0 = N-terminus, @p peptide_length + 1 = C-terminus)
      /// for the case where the optional @c location attribute is absent. The position is derived from the terminal
      /// specificity of the referenced modification(s) in the ModificationsDB. Returns -1 if no unique terminal
      /// position can be inferred, i.e. when the modification may also occur internally, or when the @c Modification
      /// names a concrete @c residues amino acid (a residue-specific/internal modification whose exact position is
      /// unknown without @c location).
      SignedSize inferModificationLocation_(MzIdentMLNode* modification_element, Size peptide_length) const;
      void parsePeptideEvidenceElements_(const std::vector<MzIdentMLNode*>& peptideEvidenceElements);
      void parseSpectrumIdentificationElements_(const std::vector<MzIdentMLNode*>& spectrumIdentificationElements);
      void parseSpectrumIdentificationProtocolElements_(const std::vector<MzIdentMLNode*>& spectrumIdentificationProtocolElements);
      void parseInputElements_(const std::vector<MzIdentMLNode*>& inputElements);
      /// Reads a SpectrumIdentificationResult of the SpectrumIdentificationList with the id @p spectrum_identification_list_id
      void parseSpectrumIdentificationResultElement_(MzIdentMLNode* element_res, const std::string& spectrum_identification_list_id);
      /// The index in pro_id_ of the protein identification run of a SpectrumIdentificationList; throws Exception::InvalidParameter if no SpectrumIdentification references the list
      Size runIndexOfList_(const std::string& spectrum_identification_list_ref) const;
      /// The protein identification run of a SpectrumIdentificationList; throws Exception::InvalidParameter if no SpectrumIdentification references the list
      ProteinIdentification& runOfList_(const std::string& spectrum_identification_list_ref);
      void parseSpectrumIdentificationItemSetXLMS(std::set<std::string>::const_iterator set_it, std::multimap<std::string, int> xl_val_map, MzIdentMLNode* element_res, const std::string& spectrumID);
      void parseSpectrumIdentificationItemElement_(MzIdentMLNode* spectrumIdentificationItemElement, PeptideIdentification& spectrum_identification, std::string& spectrumIdentificationList_ref);
      /**
        @brief Reads the protein inference results (ProteinDetectionList) into the identification run they were inferred for

        The run is the one of the SpectrumIdentificationList referenced by the ProteinDetection of the list. Each
        ProteinAmbiguityGroup becomes an indistinguishable protein group (its "leading protein"s; without such
        annotations, the proteins with the largest sets of supporting peptides), clusters ("cluster identifier") become
        general protein groups. See MzIdentMLProteinDetection for the mapping.
      */
      void parseProteinDetectionListElements_(const std::vector<MzIdentMLNode*>& proteinDetectionListElements,
                                              const std::vector<MzIdentMLNode*>& proteinDetectionElements,
                                              const std::vector<MzIdentMLNode*>& proteinDetectionProtocolElements);
      static ProteinIdentification::SearchParameters findSearchParameters_(std::pair<CVTermList, std::map<std::string, DataValue> > as_params);
      //@}

      /**@name Helper functions to build a DOM tree from the internal id structures*/
      void buildCvList_(xmlNode* cvElements);
      void buildAnalysisSoftwareList_(xmlNode* analysisSoftwareElements);
      void buildSequenceCollection_(xmlNode* sequenceCollectionElements);
      void buildAnalysisCollection_(xmlNode* analysisCollectionElements);
      void buildAnalysisProtocolCollection_(xmlNode* protocolElements);
      void buildInputDataCollection_(xmlNode* inputElements);
      void buildEnclosedCV_(xmlNode* parentElement, const std::string& encel, const std::string& acc, const std::string& name, const std::string& cvref);
      void buildAnalysisDataCollection_(xmlNode* analysisElements);
      //@}


private:
      struct SaxReader;

      MzIdentMLDOMHandler();
      MzIdentMLDOMHandler(const MzIdentMLDOMHandler& rhs);
      MzIdentMLDOMHandler& operator=(const MzIdentMLDOMHandler& rhs);

      ///Struct to hold the used analysis software for that file
      struct AnalysisSoftware
      {
        std::string name;
        std::string version;
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
        std::string sequence;
        std::string database_ref;
        std::string accession;
        CVTermList cvs;
      };
      ///Struct to hold the information from the SpectrumIdentification xml tag
      struct SpectrumIdentification
      {
        std::string spectra_data_ref;
        std::string search_database_ref;
        std::string spectrum_identification_protocol_ref;
        std::string spectrum_identification_list_ref;
      };
      ///Struct to hold the information from the ModificationParam xml tag
      struct ModificationParam
      {
        std::string fixed_mod;
        long double mass_delta;
        std::string residues;
        CVTermList modification_param_cvs;
        CVTermList specificities;
      };
      ///Struct to hold the information from the SpectrumIdentificationProtocol xml tag
      struct SpectrumIdentificationProtocol
      {
        CVTerm searchtype;
        std::string enzyme;
        CVTermList parameter_cvs;
        std::map<std::string, DataValue> parameter_ups;
//        std::vector<ModificationParam> modification_parameter;
        CVTermList modification_parameter;
        long double precursor_tolerance;
        long double fragment_tolerance;
        CVTermList threshold_cvs;
        std::map<std::string, DataValue> threshold_ups;
      };
      ///Struct to hold the information from the DatabaseInput xml tag
      struct DatabaseInput
      {
        std::string name;
        std::string location;
        std::string version;
        DateTime date;
      };
      ///Struct to hold the information from the ProteinDetectionHypothesis xml tag
      struct ProteinHypothesis
      {
        std::string id;
        std::string dbsequence_ref;
        std::string accession;
        bool leading = false;         ///< annotated as "leading protein"
        bool non_leading = false;     ///< annotated as "non-leading protein"
        bool representative = false;  ///< annotated as "group representative"
        bool pass_threshold = true;   ///< the passThreshold attribute
        std::vector<std::string> same_set; ///< ids of the ProteinDetectionHypotheses annotated as same-set proteins
        std::set<std::string> peptides; ///< unmodified sequences of the supporting peptides
        std::vector<CVTerm> cv_params; ///< in document order
        std::vector<std::pair<std::string, DataValue>> user_params; ///< in document order
      };
      ///Struct to hold the information from the ProteinAmbiguityGroup xml tag
      struct ProteinAmbiguityGroup
      {
        std::vector<ProteinHypothesis> hypotheses; ///< only those with a known DBSequence
        std::vector<CVTerm> cv_params; ///< in document order
        std::vector<std::pair<std::string, DataValue>> user_params; ///< in document order
      };

      /// Parses a ProteinAmbiguityGroup with its ProteinDetectionHypotheses
      ProteinAmbiguityGroup parseProteinAmbiguityGroupElement_(MzIdentMLNode* proteinAmbiguityGroupElement);
      /// Adds a cvParam or userParam element to @p cv_params or @p user_params (other elements are ignored)
      void collectParam_(MzIdentMLNode* param, std::vector<CVTerm>& cv_params, std::vector<std::pair<std::string, DataValue>>& user_params);
      /// Adds a protein hit for @p db to the run of the SpectrumIdentificationList @p spectrum_identification_list_ref,
      /// unless the run already has a hit with that accession
      void addProteinHit_(const std::string& spectrum_identification_list_ref, const DBSequence& db, bool is_decoy);


      std::unique_ptr<XMLHandler> xml_handler_ = nullptr;

      //from AnalysisSoftware
      std::string search_engine_;
      std::string search_engine_version_;
      //mapping from AnalysisSoftware
      std::map<std::string, AnalysisSoftware> as_map_; ///< mapping AnalysisSoftware id -> AnalysisSoftware

      //mapping from DataCollection Inputs
      std::map<std::string, std::string> sr_map_; ///< mapping sourcefile id -> sourcefile location
      std::map<std::string, std::string> sd_map_; ///< mapping spectradata id -> spectradata location
      std::map<std::string, DatabaseInput> db_map_; ///< mapping database id -> DatabaseInput

      //mapping from SpectrumIdentification - SpectrumIdentification will be the new IdentificationRuns
      std::map<std::string, SpectrumIdentification> si_map_; ///< mapping SpectrumIdentification id -> SpectrumIdentification (id refs)
      std::map<std::string, size_t> si_pro_map_; ///< mapping SpectrumIdentificationList id -> index to ProteinIdentification in pro_id_

      //mapping from SpectrumIdentificationProtocol
      std::map<std::string, SpectrumIdentificationProtocol> sp_map_; ///< mapping SpectrumIdentificationProtocol id -> SpectrumIdentificationProtocol

      //mapping from SequenceCollection
      std::map<std::string, AASequence> pep_map_; ///< mapping Peptide id -> Sequence
      std::map<std::string, PeptideEvidence> pe_ev_map_; ///< mapping PeptideEvidence id -> PeptideEvidence
      std::map<std::string, std::string> pv_db_map_; ///< mapping PeptideEvidence id -> DBSequence id
      std::multimap<std::string, std::string> p_pv_map_; ///< mapping Peptide id -> PeptideEvidence id, multiple PeptideEvidences can have equivalent Peptides.
      std::map<std::string, std::string> pv_p_map_; ///< mapping PeptideEvidence id -> Peptide id
      std::map<std::string, DBSequence> db_sq_map_; ///< mapping DBSequence id -> Sequence
      std::unordered_map<Size, std::unordered_set<std::string>> protein_accessions_of_run_; ///< mapping run index in pro_id_ -> accessions of its protein hits

      std::list<std::list<std::string> > hit_pev_; ///< writing help only

      bool xl_ms_search_ = false; ///< is true when reading a file containing Cross-Linking MS search results
      std::map<std::string, std::string> xl_id_donor_map_; ///< mapping Peptide id -> crosslink donor value
      //std::map<std::string, std::string> xl_id_acceptor_map_; ///< mapping Peptide id -> crosslink acceptor value
      std::map<std::string, std::string> xl_id_acceptor_map_; ///< mapping  peptide id of acceptor peptide -> crosslink acceptor value
      std::map<std::string, SignedSize> xl_donor_pos_map_; ///< mapping donor value -> cross-link modification location
      std::map<std::string, SignedSize> xl_acceptor_pos_map_; ///< mapping acceptor value -> cross-link modification location
      std::map<std::string, double> xl_mass_map_; ///< mapping Peptide id -> cross-link mass
      std::map<std::string, std::string> xl_mod_map_; ///< mapping peptide id -> cross-linking reagent name

      /// cached CV child term sets (computed once, reused per PSM)
      std::set<std::string> q_score_child_terms_;
      std::set<std::string> e_score_child_terms_;
      std::set<std::string> specific_score_child_terms_;
      std::set<std::string> protein_score_terms_; ///< protein- and group-level statistics and search engine specific scores (for ProteinDetectionList)
    };
  } // namespace Internal
} // namespace OpenMS

