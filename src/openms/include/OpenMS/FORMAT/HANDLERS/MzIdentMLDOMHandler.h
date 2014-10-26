// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer $
// $Authors: Mathias Walzer$
// --------------------------------------------------------------------------

#ifndef OPENMS_FORMAT_HANDLERS_MZIDENTMLDOMHANDLER_H
#define OPENMS_FORMAT_HANDLERS_MZIDENTMLDOMHANDLER_H

#include <OpenMS/KERNEL/StandardTypes.h>

#include <OpenMS/FORMAT/ControlledVocabulary.h>
#include <OpenMS/METADATA/ProteinHit.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <xercesc/dom/DOM.hpp>
#include <xercesc/dom/DOMDocument.hpp>
#include <xercesc/dom/DOMDocumentType.hpp>
#include <xercesc/dom/DOMElement.hpp>
#include <xercesc/dom/DOMImplementation.hpp>
#include <xercesc/dom/DOMImplementationLS.hpp>
#include <xercesc/dom/DOMNodeIterator.hpp>
#include <xercesc/dom/DOMNodeList.hpp>
#include <xercesc/dom/DOMText.hpp>
#include <xercesc/util/OutOfMemoryException.hpp>
#include <xercesc/util/XMLString.hpp>
#include <xercesc/util/PlatformUtils.hpp>
#include <xercesc/framework/LocalFileFormatTarget.hpp>

#include <xercesc/parsers/XercesDOMParser.hpp>
#include <xercesc/util/XMLUni.hpp>
#include <xercesc//framework/psvi/XSValue.hpp>

#include <string>
#include <stdexcept>
#include <vector>
#include <map>

#include <OpenMS/CONCEPT/LogStream.h>

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
        @brief XML handler for MzIdentMLFile

        @note Do not use this class. It is only needed in MzIdentMLFile.
    */
    class OPENMS_DLLAPI MzIdentMLDOMHandler
    {
public:
      /**@name Constructors and destructor */
      //@{
      /// Constructor for a write-only handler for internal identification structures
      MzIdentMLDOMHandler(const std::vector<ProteinIdentification> & pro_id, const std::vector<PeptideIdentification> & pep_id, const String & version, const ProgressLogger & logger);

      /// Constructor for a read-only handler for internal identification structures
      MzIdentMLDOMHandler(std::vector<ProteinIdentification> & pro_id, std::vector<PeptideIdentification> & pep_id, const String & version, const ProgressLogger & logger);

      /// Destructor
      virtual ~MzIdentMLDOMHandler();
      //@}

      /// Provides the functionality of reading a mzid with a handler object
      void readMzIdentMLFile(const std::string &mzid_file) throw(std::runtime_error);
      /// Providese the functionality to write a mzid with a handler object
      void writeMzIdentMLFile(const std::string &mzid_file) throw(std::runtime_error);

protected:
      /// Progress logger
      const ProgressLogger & logger_;

      ///Controlled vocabulary (psi-ms from OpenMS/share/OpenMS/CV/psi-ms.obo)
      ControlledVocabulary cv_;
      ///Controlled vocabulary for modifications (unimod from OpenMS/share/OpenMS/CV/unimod.obo)
      ControlledVocabulary unimod_;

      ///internal Identification Item for proteins
      std::vector<ProteinIdentification> * pro_id_;
      ///Identification Item for peptides
      std::vector<PeptideIdentification> * pep_id_;

      const std::vector<ProteinIdentification> * cpro_id_;
      const std::vector<PeptideIdentification> * cpep_id_;

      /// Looks up a child CV term of @p parent_accession with the name @p name. If no such term is found, an empty term is returned.
      ControlledVocabulary::CVTerm getChildWithName_(const String & parent_accession, const String & name) const;

      std::pair<CVTermList, std::map<String,DataValue> > parseParamGroup_( xercesc::DOMNodeList * paramGroup);
      CVTerm parseCvParam_( xercesc::DOMElement* param);
      std::pair<String, DataValue> parseUserParam_( xercesc::DOMElement* param );
      void parseAnalysisSoftwareList_( xercesc::DOMNodeList * analysisSoftwareElements);
      void parseDBSequenceElements_( xercesc::DOMNodeList * dbSequenceElements);
      void parsePeptideElements_( xercesc::DOMNodeList * peptideElements);
      AASequence parsePeptideSiblings_( xercesc::DOMNodeList * peptideSiblings);
      void parsePeptideEvidenceElements_( xercesc::DOMNodeList * peptideEvidenceElements);
      void parseSpectrumIdentificationElements_( xercesc::DOMNodeList * spectrumIdentificationElements);
      void parseSpectrumIdentificationProtocolElements_( xercesc::DOMNodeList * spectrumIdentificationProtocolElements);
      void parseInputElements_( xercesc::DOMNodeList * inputElements);
      void parseSpectrumIdentificationResultElements_( xercesc::DOMNodeList * spectrumIdentificationResultElements);
      void parseSpectrumIdentificationItemElement_(xercesc::DOMElement * spectrumIdentificationItemElement, PeptideIdentification &spectrum_identification, String& spectrumIdentificationList_ref);
      void parseProteinDetectionHypothesisElement_( xercesc::DOMElement * proteinDetectionHypothesisElement, ProteinIdentification& protein_identification);
      void parseProteinAmbiguityGroupElement_(xercesc::DOMElement * proteinAmbiguityGroupElement, ProteinIdentification& protein_identification);
      void parseProteinDetectionListElements_( xercesc::DOMNodeList * proteinDetectionListElements);

      void buildCvList_(xercesc::DOMElement * cvElements);
      void buildAnalysisSoftwareList_(xercesc::DOMElement * analysisSoftwareElements);
      void buildSequenceCollection_(xercesc::DOMElement * sequenceCollectionElements);
      void buildAnalysisCollection_(xercesc::DOMElement * analysisCollectionElements);
      void buildAnalysisProtocolCollection_(xercesc::DOMElement * protocolElements);
      void buildInputDataCollection_(xercesc::DOMElement * inputElements);
      void buildEnclosedCV_(xercesc::DOMElement * parentElement, String encel, String acc, String name, String cvref);
      void buildAnalysisDataCollection_(xercesc::DOMElement * analysisElements);

      ProteinIdentification::SearchParameters findSearchParameters_(std::pair<CVTermList,std::map<String,DataValue> > as_params);

private:
      MzIdentMLDOMHandler();
      MzIdentMLDOMHandler(const MzIdentMLDOMHandler & rhs);
      MzIdentMLDOMHandler & operator=(const MzIdentMLDOMHandler & rhs);

      struct PeptideEvidence
      {
        int start;
        int stop;
        char pre;
        char post;
        bool idec;
      };
      struct DBSequence
      {
        String sequence;
        String database_ref;
        String accession;
        CVTermList cvs;
      };
      struct SpectrumIdentification
      {
        String spectra_data_ref;
        String search_database_ref;
        String spectrum_identification_protocol_ref;
        String spectrum_identification_list_ref;
      };
      struct ModificationParam
      {
        String fixed_mod;
        long double mass_delta;
        String residues;
        CVTermList modification_param_cvs;
        CVTermList specificities;
      };
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
      struct DatabaseInput
      {
        String name;
        String location;
        String version;
        DateTime date;
      };

      XMLCh* TAG_root;
      XMLCh* TAG_CV;
      XMLCh* ATTR_name;

      xercesc::XercesDOMParser mzid_parser_;

      String search_engine_;
      String search_engine_version_;

      std::map<String, AASequence> pep_map_; //holding Peptide
      std::map<String, PeptideEvidence> pe_ev_map_; //holding PeptideEvidence
      std::map<String, String> pv_db_map_; //mapping PeptideEvidence to DBSequence
      std::multimap<String, String> p_pv_map_; //mapping Peptide(reference) -> PeptideEvidence(ID)
      std::map<String, DBSequence> db_sq_map_; //holding DBSequence
      std::map<String, SpectrumIdentification> si_map_; //holding the reference connections of db, hits, spectra - SpectrumIdentification will be the new IdentificationRuns
      std::map<String, ProteinIdentification*> si_pro_map_; //Proteinidentification for each SI list
      std::map<String, SpectrumIdentificationProtocol> sp_map_; //holding DBSequence
      std::map<String, String > input_source_; //holding sourcefiles location
      std::map<String, String > input_spectra_data_; //holding spectradata files location
      std::map<String, DatabaseInput> input_dbs_; //holding used databases location
      std::list<std::list<String> > hit_pev_; //writing help only

    };
  }   // namespace Internal
} // namespace OpenMS

#endif
