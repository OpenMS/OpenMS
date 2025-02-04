// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Mathias Walzer, Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzIdentMLDOMHandler.h>
#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLHelper.h>
#include <OpenMS/CONCEPT/Constants.h>

#include <boost/lexical_cast.hpp>

#include <sys/stat.h>
#include <cerrno>

using namespace std;
using namespace xercesc;

namespace OpenMS::Internal
{

    //TODO remodel CVTermList
    //TODO extend CVTermlist with CVCollection functionality for complete replacement??
    //TODO general id openms struct for overall parameter for one id run
    MzIdentMLDOMHandler::MzIdentMLDOMHandler(const vector<ProteinIdentification>& pro_id, const vector<PeptideIdentification>& pep_id, const String& version, const ProgressLogger& logger) :
      logger_(logger),
      cv_(ControlledVocabulary::getPSIMSCV()),
      //~ ms_exp_(0),
      cpro_id_(&pro_id),
      cpep_id_(&pep_id),
      schema_version_(version),
      mzid_parser_()
    {
      unimod_.loadFromOBO("UNIMOD", File::find("/CV/unimod.obo"));

      try
      {
        XMLPlatformUtils::Initialize(); // Initialize Xerces infrastructure
      }
      catch (XMLException& e)
      {
        char* message = XMLString::transcode(e.getMessage());
        OPENMS_LOG_ERROR << "XML toolkit initialization error: " << message << endl;
        XMLString::release(&message);
        // throw exception here to return ERROR_XERCES_INIT
      }

      // Tags and attributes used in XML file.
      // Can't call transcode till after Xerces Initialize()
      xml_root_tag_ptr_ = XMLString::transcode("MzIdentML");
      xml_cvparam_tag_ptr_ = XMLString::transcode("cvParam");
      xml_name_attr_ptr_ = XMLString::transcode("option_a");

    }

    MzIdentMLDOMHandler::MzIdentMLDOMHandler(vector<ProteinIdentification>& pro_id, vector<PeptideIdentification>& pep_id, const String& version, const ProgressLogger& logger) :
      logger_(logger),
      cv_(ControlledVocabulary::getPSIMSCV()),
      //~ ms_exp_(0),
      pro_id_(&pro_id),
      pep_id_(&pep_id),
      schema_version_(version),
      mzid_parser_(),
      xl_ms_search_(false)
    {
      unimod_.loadFromOBO("UNIMOD", File::find("/CV/unimod.obo"));

      try
      {
        XMLPlatformUtils::Initialize(); // Initialize Xerces infrastructure
      }
      catch (XMLException& e)
      {
        char* message = XMLString::transcode(e.getMessage());
        OPENMS_LOG_ERROR << "XML toolkit initialization error: " << message << endl;
        XMLString::release(&message);
      }

      // Tags and attributes used in XML file.
      // Can't call transcode till after Xerces Initialize()
      xml_root_tag_ptr_ = XMLString::transcode("MzIdentML");
      xml_cvparam_tag_ptr_ = XMLString::transcode("cvParam");
      xml_name_attr_ptr_ = XMLString::transcode("name");
    }

    /*
     *  Class destructor frees memory used to hold the XML tag and
     *  attribute definitions. It also terminates use of the xerces-C
     *  framework.
     */
    MzIdentMLDOMHandler::~MzIdentMLDOMHandler()
    {
      try
      {
        XMLString::release(&xml_root_tag_ptr_);
        XMLString::release(&xml_cvparam_tag_ptr_);
        XMLString::release(&xml_name_attr_ptr_);
//         if(m_name)   XMLString::release( &m_name ); //releasing you here is releasing you twice, dunno yet why?!
      }
      catch (...)
      {
        OPENMS_LOG_ERROR << "Unknown exception encountered in 'TagNames' destructor" << endl;
      }

      // Terminate Xerces
      try
      {
        XMLPlatformUtils::Terminate(); // Terminate after release of memory
      }
      catch (xercesc::XMLException& e)
      {
        char* message = xercesc::XMLString::transcode(e.getMessage());
        OPENMS_LOG_ERROR << "XML toolkit teardown error: " << message << endl;
        XMLString::release(&message);
      }
    }

    /*
     * reads a mzid file into handlers pro_id_ & pep_id_ members
     */
    void MzIdentMLDOMHandler::readMzIdentMLFile(const std::string& mzid_file)
    {
      // Test to see if the file is ok.
      struct stat fileStatus;

      xml_handler_ = make_unique<XMLHandler>(mzid_file, schema_version_);

      errno = 0;
      if (stat(mzid_file.c_str(), &fileStatus) == -1) // ==0 ok; ==-1 error
      {
        if (errno == ENOENT) // errno declared by include file errno.h
        {
          throw (runtime_error("Path file_name does not exist, or path is an empty string."));
        }
        else if (errno == ENOTDIR)
        {
          throw (runtime_error("A component of the path is not a directory."));
        }
        // On MSVC 2008, the ELOOP constant is not declared and thus introduces a compile error
        //else if (errno == ELOOP)
        //  throw (runtime_error("Too many symbolic links encountered while traversing the path."));
        else if (errno == EACCES)
        {
          throw (runtime_error("Permission denied."));
        }
        else if (errno == ENAMETOOLONG)
        {
          throw (runtime_error("File can not be read."));
        }
      }

      // Configure DOM parser.
      mzid_parser_.setValidationScheme(XercesDOMParser::Val_Never);
      mzid_parser_.setDoNamespaces(false);
      mzid_parser_.setDoSchema(false);
      mzid_parser_.setLoadExternalDTD(false);

      try
      {
        mzid_parser_.parse(mzid_file.c_str());
      }
      catch (xercesc::XMLException& e)
      {
        char* message = xercesc::XMLString::transcode(e.getMessage());
        //          ostringstream errBuf;
        //          errBuf << "Error parsing file: " << message << flush;
        OPENMS_LOG_ERROR << "XERCES error parsing file: " << message << flush << endl;
        XMLString::release(&message);
      }

      // we adopt/own the document so we can free it in case this DOMHandler is reused on another file.
      xercesc::DOMDocument* xmlDoc = mzid_parser_.adoptDocument();

      try
      {
        // Catch special case: Cross-Linking MS
        DOMNodeList* additionalSearchParams = xmlDoc->getElementsByTagName(CONST_XMLCH("AdditionalSearchParams"));
        const  XMLSize_t as_node_count = additionalSearchParams->getLength();

        for (XMLSize_t i = 0; i < as_node_count; ++i)
        {
          DOMNode* current_sp = additionalSearchParams->item(i);

          DOMElement* element_SearchParams = dynamic_cast<xercesc::DOMElement*>(current_sp);
          String cross_linking_search = StringManager::convert(element_SearchParams->getAttribute(CONST_XMLCH("id")));
          DOMElement* child = element_SearchParams->getFirstElementChild();

          while (child && !xl_ms_search_)
          {
            String accession = StringManager::convert(child->getAttribute(CONST_XMLCH("accession")));
            if (accession == "MS:1002494") // accession for "crosslinking search"
            {
              xl_ms_search_ = true;
            }
            child = child->getNextElementSibling();
          }
        }

        if (xl_ms_search_)
        {
          OPENMS_LOG_DEBUG << "Reading a Cross-Linking MS file." << endl;
        }

        // 0. AnalysisSoftwareList {0,1}
        DOMNodeList* analysisSoftwareElements = xmlDoc->getElementsByTagName(CONST_XMLCH("AnalysisSoftware"));
        parseAnalysisSoftwareList_(analysisSoftwareElements);

        // 1. DataCollection {1,1}
        DOMNodeList* spectraDataElements = xmlDoc->getElementsByTagName(CONST_XMLCH("SpectraData"));
        if (spectraDataElements->getLength() == 0)
        {
          throw(runtime_error("No SpectraData nodes"));
        }
        parseInputElements_(spectraDataElements);

        // 1.2. SearchDatabase {0,unbounded}
        DOMNodeList* searchDatabaseElements = xmlDoc->getElementsByTagName(CONST_XMLCH("SearchDatabase"));
        parseInputElements_(searchDatabaseElements);

        // 1.1 SourceFile {0,unbounded}
        DOMNodeList* sourceFileElements = xmlDoc->getElementsByTagName(CONST_XMLCH("SourceFile"));
        parseInputElements_(sourceFileElements);

        // 2. SpectrumIdentification  {1,unbounded} ! creates identification runs (or ProteinIdentifications)
        DOMNodeList* spectrumIdentificationElements = xmlDoc->getElementsByTagName(CONST_XMLCH("SpectrumIdentification"));
        if (spectrumIdentificationElements->getLength() == 0)
        {
          throw(runtime_error("No SpectrumIdentification nodes"));
        }
        parseSpectrumIdentificationElements_(spectrumIdentificationElements);

        // 3. AnalysisProtocolCollection {1,1} SpectrumIdentificationProtocol  {1,unbounded} ! identification run parameters
        DOMNodeList* spectrumIdentificationProtocolElements = xmlDoc->getElementsByTagName(CONST_XMLCH("SpectrumIdentificationProtocol"));
        if (spectrumIdentificationProtocolElements->getLength() == 0)
        {
          throw(runtime_error("No SpectrumIdentificationProtocol nodes"));
        }
        parseSpectrumIdentificationProtocolElements_(spectrumIdentificationProtocolElements);

        // 4. SequenceCollection nodes {0,1} DBSequenceElement {1,unbounded} Peptide {0,unbounded} PeptideEvidence {0,unbounded}
        DOMNodeList* dbSequenceElements = xmlDoc->getElementsByTagName(CONST_XMLCH("DBSequence"));
        parseDBSequenceElements_(dbSequenceElements);

        DOMNodeList* peptideElements = xmlDoc->getElementsByTagName(CONST_XMLCH("Peptide"));
        parsePeptideElements_(peptideElements);

        DOMNodeList* peptideEvidenceElements = xmlDoc->getElementsByTagName(CONST_XMLCH("PeptideEvidence"));
        parsePeptideEvidenceElements_(peptideEvidenceElements);
//          mzid_parser_.resetDocumentPool(); //segfault prone: do not use!

        // 5. AnalysisSampleCollection ??? contact stuff

        // 6. AnalysisCollection {1,1} - build final structures PeptideIdentification (and hits)

        // 6.1 SpectrumIdentificationList {0,1}
        DOMNodeList* spectrumIdentificationListElements = xmlDoc->getElementsByTagName(CONST_XMLCH("SpectrumIdentificationList"));
        if (spectrumIdentificationListElements->getLength() == 0)
        {
          throw(runtime_error("No SpectrumIdentificationList nodes"));
        }
        parseSpectrumIdentificationListElements_(spectrumIdentificationListElements);

        // 6.2 ProteinDetection {0,1}
        DOMNodeList* parseProteinDetectionListElements = xmlDoc->getElementsByTagName(CONST_XMLCH("ProteinDetectionList"));
        parseProteinDetectionListElements_(parseProteinDetectionListElements);

        for (vector<ProteinIdentification>::iterator it = pro_id_->begin(); it != pro_id_->end(); ++it)
        {
          it->sort();
        }
        //note: PeptideIdentification sorting here not necessary any more, due to sorting according to cv in SpectrumIdentificationResult

      }
      catch (xercesc::XMLException& e)
      {
        char* message = xercesc::XMLString::transcode(e.getMessage());
//          ostringstream errBuf;
//          errBuf << "Error parsing file: " << message << flush;
        OPENMS_LOG_ERROR << "XERCES error traversing DOM: " << message << flush << endl;
        XMLString::release(&message);
      }
      xmlDoc->release();
      if (xl_ms_search_)
      {
        OPXLHelper::addProteinPositionMetaValues(*this->pep_id_);
        OPXLHelper::addBetaAccessions(*this->pep_id_);
        OPXLHelper::addXLTargetDecoyMV(*this->pep_id_);
        OPXLHelper::removeBetaPeptideHits(*this->pep_id_);
        OPXLHelper::computeDeltaScores(*this->pep_id_);
        OPXLHelper::addPercolatorFeatureList((*this->pro_id_)[0]);
      }
    }

    void MzIdentMLDOMHandler::writeMzIdentMLFile(const std::string& mzid_file)
    {
      xml_handler_ = make_unique<XMLHandler>(mzid_file, schema_version_);

      DOMImplementation* impl =  DOMImplementationRegistry::getDOMImplementation(CONST_XMLCH("XML 1.0")); //XML 3?!
      if (impl != nullptr)
      {
        try
        {
          xercesc::DOMDocument* xmlDoc = impl->createDocument(
            CONST_XMLCH("http://psidev.info/psi/pi/mzIdentML/1.1"),
            CONST_XMLCH("MzIdentML"), // root element name
            nullptr); // document type object (DTD).

          DOMElement* rootElem = xmlDoc->getDocumentElement();
          rootElem->setAttribute(CONST_XMLCH("version"),
                                 StringManager::convertPtr(schema_version_).get());
          rootElem->setAttribute(CONST_XMLCH("xsi:schemaLocation"),
                                 CONST_XMLCH("http://psidev.info/psi/pi/mzIdentML/1.1 ../../schema/mzIdentML1.1.0.xsd"));
          rootElem->setAttribute(CONST_XMLCH("creationDate"),
                                 StringManager::convertPtr(String(DateTime::now().getDate() + "T" + DateTime::now().getTime())).get());

          // * cvList *
          DOMElement* cvl_p = xmlDoc->createElement(CONST_XMLCH("cvList")); // TODO add generically
          buildCvList_(cvl_p);
          rootElem->appendChild(cvl_p);

          // * AnalysisSoftwareList *
          DOMElement* asl_p = xmlDoc->createElement(CONST_XMLCH("AnalysisSoftwareList"));
          for (vector<ProteinIdentification>::const_iterator pi = cpro_id_->begin(); pi != cpro_id_->end(); ++pi)
          {
//                  search_engine_version_ = pi->getSearchEngineVersion();
//                  search_engine_ = pi->getSearchEngine();
          }
          buildAnalysisSoftwareList_(asl_p);
          rootElem->appendChild(asl_p);

//              // * AnalysisSampleCollection *
//              DOMElement* asc_p = xmlDoc->createElement(CONST_XMLCH("AnalysisSampleCollection"));
//              buildAnalysisSampleCollection_(asc_p);
//              rootElem->appendChild(asc_p);

          // * SequenceCollection *
          DOMElement* sc_p = xmlDoc->createElement(CONST_XMLCH("SequenceCollection"));

          for (vector<ProteinIdentification>::const_iterator pi = cpro_id_->begin(); pi != cpro_id_->end(); ++pi)
          {
            String dbref = pi->getSearchParameters().db +  pi->getSearchParameters().db_version + pi->getSearchParameters().taxonomy; //TODO @mths : this needs to be more unique, btw add tax etc. as cv to DBSequence
            for (vector<ProteinHit>::const_iterator ph = pi->getHits().begin(); ph != pi->getHits().end(); ++ph)
            {
              CVTermList cvs;
              DBSequence temp_struct = {ph->getSequence(), dbref, ph->getAccession(), cvs};
              db_sq_map_.insert(make_pair(ph->getAccession(), temp_struct));
            }
          }

          set<AASequence> pepset;
          for (vector<PeptideIdentification>::const_iterator pi = cpep_id_->begin(); pi != cpep_id_->end(); ++pi)
          {
            for (vector<PeptideHit>::const_iterator ph = pi->getHits().begin(); ph != pi->getHits().end(); ++ph)
            {
              list<String> pepevs;
              for (vector<OpenMS::PeptideEvidence>::const_iterator pev = ph->getPeptideEvidences().begin(); pev != ph->getPeptideEvidences().end(); ++pev)
              {
                String pepevref = String("OpenMS") + String(UniqueIdGenerator::getUniqueId());
                pv_db_map_.insert(make_pair(pepevref, pev->getProteinAccession()));
                pepevs.push_back(pepevref);
                bool idec = (String(ph->getMetaValue(Constants::UserParam::TARGET_DECOY))).hasSubstring("decoy");
                PeptideEvidence temp_struct = {pev->getStart(), pev->getEnd(), pev->getAABefore(), pev->getAAAfter(), idec}; //TODO @ mths : completely switch to PeptideEvidence
                pe_ev_map_.insert(make_pair(pepevref, temp_struct)); // TODO @ mths : double check start & end & chars for before & after
              }
              hit_pev_.push_back(pepevs);

              String pepref = String("OpenMS") + String(UniqueIdGenerator::getUniqueId());
              if (pepset.find(ph->getSequence()) != pepset.end())
              {
                pepset.insert(ph->getSequence());
                pep_map_.insert(make_pair(pepref, ph->getSequence()));
                for (list<String>::iterator pepevref = pepevs.begin(); pepevref != pepevs.end(); ++pepevref)
                {
                  p_pv_map_.insert(make_pair(*pepevref, pepref));
                }
              }
            }
          }

          buildSequenceCollection_(sc_p);
          rootElem->appendChild(sc_p);

          // * AnalysisCollection *
          DOMElement* analysis_c_p = xmlDoc->createElement(CONST_XMLCH("AnalysisCollection"));
          buildAnalysisCollection_(analysis_c_p);
          rootElem->appendChild(analysis_c_p);

          // * AnalysisProtocolCollection *
          DOMElement* apc_p = xmlDoc->createElement(CONST_XMLCH("AnalysisProtocolCollection"));
          buildAnalysisCollection_(apc_p);
          rootElem->appendChild(apc_p);

          // * DataCollection *
          DOMElement* dc_p = xmlDoc->createElement(CONST_XMLCH("DataCollection"));
          rootElem->appendChild(dc_p);
          DOMElement* in_p = dc_p->getOwnerDocument()->createElement(CONST_XMLCH("Inputs"));
          DOMElement* ad_p = dc_p->getOwnerDocument()->createElement(CONST_XMLCH("AnalysisData"));
          dc_p->appendChild(in_p);
          dc_p->appendChild(ad_p);

          // * BibliographicReference *
          DOMElement* br_p = xmlDoc->createElement(CONST_XMLCH("BibliographicReference"));
          br_p->setAttribute(CONST_XMLCH("authors"), CONST_XMLCH("all"));
          rootElem->appendChild(br_p);

          // * Serialisation *
          DOMLSSerializer* serializer = ((DOMImplementationLS*)impl)->createLSSerializer();
          // serializer gets prettyprint and stuff
          if (serializer->getDomConfig()->canSetParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true))
          {
            serializer->getDomConfig()->setParameter(XMLUni::fgDOMWRTDiscardDefaultContent, true);
          }
          if (serializer->getDomConfig()->canSetParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true))
          {
            serializer->getDomConfig()->setParameter(XMLUni::fgDOMWRTFormatPrettyPrint, true);
          }
//              // optionally implement DOMErrorHandler (e.g. MyDOMErrorHandler) and set it to the serializer
//              DOMErrorHandler* errHandler = new myDOMErrorHandler();
//              serializer->getDomConfig()->setParameter(XMLUni::fgDOMErrorHandler, myErrorHandler);

          XMLFormatTarget* file_target = new LocalFileFormatTarget(mzid_file.c_str());
          DOMLSOutput* dom_output = ((DOMImplementationLS*)impl)->createLSOutput();
          dom_output->setByteStream(file_target);

          try
          {
            // do the serialization through DOMLSSerializer::write();
            serializer->write(xmlDoc, dom_output);
          }
          catch (const XMLException& toCatch)
          {
            char* message = XMLString::transcode(toCatch.getMessage());
            OPENMS_LOG_ERROR << "Serialisation exception: \n"
                      << message << "\n";
            XMLString::release(&message);
          }
          catch (const DOMException& toCatch)
          {
            char* message = XMLString::transcode(toCatch.msg);
            OPENMS_LOG_ERROR << "Serialisation exception: \n"
                      << message << "\n";
            XMLString::release(&message);
          }
          catch (...)
          {
            OPENMS_LOG_ERROR << "Unexpected exception building the document tree." << endl;
          }

          dom_output->release();
          serializer->release();
//              delete myErrorHandler;
          delete file_target;
        }
        catch (const OutOfMemoryException&)
        {
          OPENMS_LOG_ERROR << "Xerces OutOfMemoryException" << endl;
        }
        catch (const DOMException& e)
        {
          OPENMS_LOG_ERROR << "DOMException code is: " << e.code << endl;
        }
        catch (const exception& e)
        {
          OPENMS_LOG_ERROR << "An error occurred creating the document: " << e.what() << endl;
        }

      } // (inpl != NULL)
      else
      {
        OPENMS_LOG_ERROR << "Requested DOM implementation is not supported" << endl;
      }
    }

    pair<CVTermList, map<String, DataValue> > MzIdentMLDOMHandler::parseParamGroup_(DOMNodeList* paramGroup)
    {
      CVTermList ret_cv;             // cvParam
      map<String, DataValue> ret_up; // userParam
      const  XMLSize_t cv_node_count = paramGroup->getLength();
      for (XMLSize_t cvi = 0; cvi < cv_node_count; ++cvi)
      {
        DOMNode* current_cv = paramGroup->item(cvi);
        if (current_cv->getNodeType() && // true is not NULL
            current_cv->getNodeType() == DOMNode::ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          DOMElement* element_param = dynamic_cast<xercesc::DOMElement*>(current_cv);
          if (XMLString::equals(element_param->getTagName(), CONST_XMLCH("cvParam")))
          {
            ret_cv.addCVTerm(parseCvParam_(element_param));
          }
          else if (XMLString::equals(element_param->getTagName(), CONST_XMLCH("userParam")))
          {
            ret_up.insert(parseUserParam_(element_param));
          }
          else if (XMLString::equals(element_param->getTagName(), CONST_XMLCH("PeptideEvidence"))
                  || XMLString::equals(element_param->getTagName(), CONST_XMLCH("PeptideEvidenceRef"))
                  || XMLString::equals(element_param->getTagName(), CONST_XMLCH("Fragmentation"))
                  || XMLString::equals(element_param->getTagName(), CONST_XMLCH("SpectrumIdentificationItem")))
          {
            //here it's okay to do nothing
          }
          else
          {
            OPENMS_LOG_WARN << "Misplaced elements ignored in 'ParamGroup' in " << StringManager::convert(element_param->getTagName()) << endl;
          }
        }
      }
      return make_pair(ret_cv, ret_up);
    }

    CVTerm MzIdentMLDOMHandler::parseCvParam_(DOMElement* param)
    {
      if (param)
      {
        //      <cvParam accession="MS:1001469" name="taxonomy: scientific name" cvRef="PSI-MS"  value="Drosophila melanogaster"/>
        String accession = StringManager::convert(param->getAttribute(CONST_XMLCH("accession")));
        String name = StringManager::convert(param->getAttribute(CONST_XMLCH("name")));
        String cvRef = StringManager::convert(param->getAttribute(CONST_XMLCH("cvRef")));
        String value = StringManager::convert(param->getAttribute(CONST_XMLCH("value")));

        String unitAcc = StringManager::convert(param->getAttribute(CONST_XMLCH("unitAccession")));
        String unitName = StringManager::convert(param->getAttribute(CONST_XMLCH("unitName")));
        String unitCvRef = StringManager::convert(param->getAttribute(CONST_XMLCH("unitCvRef")));

        CVTerm::Unit u; // TODO @mths : make DataValue usage safe!
        if (!unitAcc.empty() && !unitName.empty())
        {
          u = CVTerm::Unit(unitAcc, unitName, unitCvRef);
          if (unitCvRef.empty())
          {
            OPENMS_LOG_WARN << "This mzid file uses a cv term with units, but without "
                     << "unit cv reference (required)! Please notify the mzid "
                     << "producer of this file. \"" << name << "\" will be read as \""
                     << unitName << "\" but further actions on this unit may fail."
                     << endl;
          }
        }
        return CVTerm(accession, name, cvRef, value, u);
      }
      else
        throw invalid_argument("no cv param here");
    }

    pair<String, DataValue> MzIdentMLDOMHandler::parseUserParam_(DOMElement* param)
    {
      if (param)
      {
        //      <userParam name="Mascot User Comment" value="Example Mascot MS-MS search for PSI mzIdentML"/>
        String name = StringManager::convert(param->getAttribute(CONST_XMLCH("name")));
        String value = StringManager::convert(param->getAttribute(CONST_XMLCH("value")));
        bool has_value = param->hasAttribute(CONST_XMLCH("value"));
        String unitAcc = StringManager::convert(param->getAttribute(CONST_XMLCH("unitAccession")));
        String unitName = StringManager::convert(param->getAttribute(CONST_XMLCH("unitName")));
        String unitCvRef = StringManager::convert(param->getAttribute(CONST_XMLCH("unitCvRef")));
        String type = StringManager::convert(param->getAttribute(CONST_XMLCH("type")));

        DataValue dv = DataValue::EMPTY;
        if (has_value)
        {
          dv = XMLHandler::fromXSDString(type, value);
        }

        // Add unit *after* creating the term
        if (!unitAcc.empty())
        {
          if (unitAcc.hasPrefix("UO:"))
          {
            dv.setUnit(unitAcc.suffix(unitAcc.size() - 3).toInt());
            dv.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
          }
          else if (unitAcc.hasPrefix("MS:"))
          {
            dv.setUnit(unitAcc.suffix(unitAcc.size() - 3).toInt());
            dv.setUnitType(DataValue::UnitType::MS_ONTOLOGY);
          }
          else
          {
            OPENMS_LOG_WARN << String("Unhandled unit '") + unitAcc + "' in tag '" + name + "'." << endl;
          }
        }

        return make_pair(name, dv);
      }
      else
      {
        OPENMS_LOG_ERROR << "No parameters found at given position." << endl;
        throw invalid_argument("no user param here");
      }
    }

    void MzIdentMLDOMHandler::parseAnalysisSoftwareList_(DOMNodeList* analysisSoftwareElements)
    {
      const  XMLSize_t as_node_count = analysisSoftwareElements->getLength();
      for (XMLSize_t swni = 0; swni < as_node_count; ++swni)
      {
        DOMNode* current_as = analysisSoftwareElements->item(swni);
        if (current_as->getNodeType() && // true is not NULL
            current_as->getNodeType() == DOMNode::ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_AnalysisSoftware = dynamic_cast<xercesc::DOMElement*>(current_as);
          String id = StringManager::convert(element_AnalysisSoftware->getAttribute(CONST_XMLCH("id")));
          DOMElement* child = element_AnalysisSoftware->getFirstElementChild();
          String swname, swversion;
          while (child)
          {
            if (XMLString::equals(child->getTagName(),CONST_XMLCH("SoftwareName"))) //must have exactly one SoftwareName
            {
              DOMNodeList* element_pg = child->getChildNodes();

              pair<CVTermList, map<String, DataValue> > swn = parseParamGroup_(element_pg);
              swversion = StringManager::convert(element_AnalysisSoftware->getAttribute(CONST_XMLCH("version")));
              if (!swn.first.getCVTerms().empty())
              {
                set<String> software_terms;
                cv_.getAllChildTerms(software_terms, "MS:1000531");
                for (map<String, vector<CVTerm> >::const_iterator it = swn.first.getCVTerms().begin(); it != swn.first.getCVTerms().end(); ++it)
                {
                  if (software_terms.find(it->first) != software_terms.end())
                  {
                    swname = it->second.front().getName();
                    break;
                  }
                }
              }
              else if (!swn.second.empty())
              {
                for (map<String, DataValue>::const_iterator up = swn.second.begin(); up != swn.second.end(); ++up)
                {
                  if (up->first.hasSubstring("name"))
                  {
                    swname = up->second.toString();
                    break;
                  }
                  else
                  {
                    swname = up->first;
                  }
                }
              }
            }
            child = child->getNextElementSibling();
          }
          if (!swname.empty() && !swversion.empty())
          {
            AnalysisSoftware temp_struct = {swname, swversion};
            as_map_.insert(make_pair(id, temp_struct));
          }
          else
          {
            OPENMS_LOG_ERROR << "No name/version found for 'AnalysisSoftware':" << id << "." << endl;
          }
        }
      }
    }

    void MzIdentMLDOMHandler::parseDBSequenceElements_(DOMNodeList* dbSequenceElements)
    {
      const  XMLSize_t dbs_node_count = dbSequenceElements->getLength();
      for (XMLSize_t c = 0; c < dbs_node_count; ++c)
      {
        DOMNode* current_dbs = dbSequenceElements->item(c);
        if (current_dbs->getNodeType() && // true is not NULL
            current_dbs->getNodeType() == DOMNode::ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_dbs = dynamic_cast<xercesc::DOMElement*>(current_dbs);
          String id = StringManager::convert(element_dbs->getAttribute(CONST_XMLCH("id")));
          String seq = "";
          String dbref = StringManager::convert(element_dbs->getAttribute(CONST_XMLCH("searchDatabase_ref")));
          String acc = StringManager::convert(element_dbs->getAttribute(CONST_XMLCH("accession")));
          CVTermList cvs;

          DOMElement* child = element_dbs->getFirstElementChild();
          while (child)
          {
            if (XMLString::equals(child->getTagName(), CONST_XMLCH("Seq")))
            {
              seq = StringManager::convert(child->getTextContent());
            }
            else if (XMLString::equals(child->getTagName(), CONST_XMLCH("cvParam")))
            {
              cvs.addCVTerm(parseCvParam_(child));
            }
            child = child->getNextElementSibling();
          }
          if (!acc.empty())
          {
            DBSequence temp_struct = {seq, dbref, acc, cvs};
            db_sq_map_.insert(make_pair(id, temp_struct));
          }
        }
      }
    }

    void MzIdentMLDOMHandler::parsePeptideElements_(DOMNodeList* peptideElements)
    {
      const  XMLSize_t pep_node_count = peptideElements->getLength();
      for (XMLSize_t c = 0; c < pep_node_count; ++c)
      {
        DOMNode* current_pep = peptideElements->item(c);
        if (current_pep->getNodeType() && // true is not NULL
            current_pep->getNodeType() == DOMNode::ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_pep = dynamic_cast<xercesc::DOMElement*>(current_pep);
          String id = StringManager::convert(element_pep->getAttribute(CONST_XMLCH("id")));

          //DOMNodeList* pep_sib = element_pep->getChildNodes();
          AASequence aas;
          try
          {
            //aas = parsePeptideSiblings_(pep_sib);
            try
            {
              aas = parsePeptideSiblings_(element_pep);
            }
            catch (Exception::MissingInformation&)
            {
              // We found an unknown modification, we could try to rescue this
              // situation. The "name" attribute, if present, may be parsable:
              //   The potentially ambiguous common identifier, such as a
              //   human-readable name for the instance.
              String name = StringManager::convert(element_pep->getAttribute(CONST_XMLCH("name")));
              if (!name.empty()) aas = AASequence::fromString(name);
            }
          }
          catch (...)
          {
            OPENMS_LOG_ERROR << "No amino acid sequence readable from 'Peptide'" << endl;
          }

          pep_map_.insert(make_pair(id, aas));
        }
      }
    }

    void MzIdentMLDOMHandler::parsePeptideEvidenceElements_(DOMNodeList* peptideEvidenceElements)
    {
      const  XMLSize_t pev_node_count = peptideEvidenceElements->getLength();
      for (XMLSize_t c = 0; c < pev_node_count; ++c)
      {
        DOMNode* current_pev = peptideEvidenceElements->item(c);
        if (current_pev->getNodeType() && // true is not NULL
            current_pev->getNodeType() == DOMNode::ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_pev = dynamic_cast<xercesc::DOMElement*>(current_pev);

//          <PeptideEvidence peptide_ref="peptide_1_1" id="PE_1_1_HSP70_ECHGR_0" start="161" end="172" pre="K" post="I" isDecoy="false" dBSequence_ref="DBSeq_HSP70_ECHGR"/>

          String id = StringManager::convert(element_pev->getAttribute(CONST_XMLCH("id")));
          String peptide_ref = StringManager::convert(element_pev->getAttribute(CONST_XMLCH("peptide_ref")));
          String dBSequence_ref = StringManager::convert(element_pev->getAttribute(CONST_XMLCH("dBSequence_ref")));
          //rest is optional !!
          int start = -1;
          int end = -1;
          try
          {
            start = StringManager::convert(element_pev->getAttribute(CONST_XMLCH("start"))).toInt();
            end = StringManager::convert(element_pev->getAttribute(CONST_XMLCH("end"))).toInt();
          }
          catch (...)
          {
            OPENMS_LOG_WARN << "'PeptideEvidence' without reference to the position in the originating sequence found." << endl;
          }
          char pre = '-';
          char post = '-';
          try
          {
            if (element_pev->hasAttribute(CONST_XMLCH("pre")))
            {
            pre = StringManager::convert(element_pev->getAttribute(CONST_XMLCH("pre")))[0];
            }
            if (element_pev->hasAttribute(CONST_XMLCH("post")))
            {
            post = StringManager::convert(element_pev->getAttribute(CONST_XMLCH("post")))[0];
          }
          }
          catch (...)
          {
            OPENMS_LOG_WARN << "'PeptideEvidence' without reference to the bordering amino acids in the originating sequence found." << endl;
          }
          bool idec = false;
          try
          {
            String d = StringManager::convert(element_pev->getAttribute(CONST_XMLCH("isDecoy")));
            if (d.hasPrefix('t') || d.hasPrefix('1'))
              idec = true;
          }
          catch (...)
          {
            OPENMS_LOG_WARN << "'PeptideEvidence' with unreadable 'isDecoy' status found." << endl;
          }
          PeptideEvidence temp_struct = {start, end, pre, post, idec};
          pe_ev_map_.insert(make_pair(id, temp_struct));
          p_pv_map_.insert(make_pair(peptide_ref, id));
          pv_db_map_.insert(make_pair(id, dBSequence_ref));
        }
      }
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationElements_(DOMNodeList* spectrumIdentificationElements)
    {
      const  XMLSize_t si_node_count = spectrumIdentificationElements->getLength();
      for (XMLSize_t c = 0; c < si_node_count; ++c)
      {
        DOMNode* current_si = spectrumIdentificationElements->item(c);
        if (current_si->getNodeType() && // true is not NULL
            current_si->getNodeType() == DOMNode::ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_si = dynamic_cast<xercesc::DOMElement*>(current_si);
          String id = StringManager::convert(element_si->getAttribute(CONST_XMLCH("id")));
          String spectrumIdentificationProtocol_ref = StringManager::convert(element_si->getAttribute(CONST_XMLCH("spectrumIdentificationProtocol_ref")));
          String spectrumIdentificationList_ref = StringManager::convert(element_si->getAttribute(CONST_XMLCH("spectrumIdentificationList_ref")));
          String spectrumIdentification_date = StringManager::convert(element_si->getAttribute(CONST_XMLCH("activityDate")));

          String searchDatabase_ref = "";
          String spectra_data_ref = "";
          DOMElement* child = element_si->getFirstElementChild();
          while (child)
          {
            if (XMLString::equals(child->getTagName(), CONST_XMLCH("InputSpectra")))
            {
              spectra_data_ref = StringManager::convert(child->getAttribute(CONST_XMLCH("spectraData_ref")));
            }
            else if (XMLString::equals(child->getTagName(), CONST_XMLCH("SearchDatabaseRef")))
            {
              searchDatabase_ref = StringManager::convert(child->getAttribute(CONST_XMLCH("searchDatabase_ref")));
            }
            child = child->getNextElementSibling();
          }
          SpectrumIdentification temp_struct = {spectra_data_ref, searchDatabase_ref, spectrumIdentificationProtocol_ref, spectrumIdentificationList_ref};
          si_map_.insert(make_pair(id, temp_struct));

          pro_id_->push_back(ProteinIdentification());
          ProteinIdentification::SearchParameters sp;
          sp.db = db_map_[searchDatabase_ref].location;
          sp.db_version = db_map_[searchDatabase_ref].version;
          pro_id_->back().setSearchParameters(sp);
          if (xl_ms_search_)
          {
            pro_id_->back().setMetaValue("SpectrumIdentificationProtocol", "MS:1002494"); // XL-MS CV term
          }

          // internally we store a list of files so convert the mzIdentML file String to a StringList
          StringList spectra_data_list;
          spectra_data_list.push_back(sd_map_[spectra_data_ref]);
          pro_id_->back().setMetaValue("spectra_data", spectra_data_list);
          if (!spectrumIdentification_date.empty())
          {
            pro_id_->back().setDateTime(DateTime::fromString(spectrumIdentification_date));
          }
          else
          {
            pro_id_->back().setDateTime(DateTime::now());
          }
          pro_id_->back().setIdentifier(UniqueIdGenerator::getUniqueId()); // no more identification wit engine/date/time!
          //TODO setIdentifier to xml id?
          si_pro_map_.insert(make_pair(spectrumIdentificationList_ref, pro_id_->size() - 1));
        }
      }
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationProtocolElements_(DOMNodeList* spectrumIdentificationProtocolElements)
    {
      const  XMLSize_t si_node_count = spectrumIdentificationProtocolElements->getLength();
      for (XMLSize_t c = 0; c < si_node_count; ++c)
      {
        ProteinIdentification::SearchParameters sp;
        DOMNode* current_sip = spectrumIdentificationProtocolElements->item(c);
        if (current_sip->getNodeType() && // true is not NULL
            current_sip->getNodeType() == DOMNode::ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_sip = dynamic_cast<xercesc::DOMElement*>(current_sip);
          String id = StringManager::convert(element_sip->getAttribute(CONST_XMLCH("id")));
          String swr = StringManager::convert(element_sip->getAttribute(CONST_XMLCH("analysisSoftware_ref")));

          CVTerm searchtype;
          String enzymename;
          CVTermList param_cv;
          map<String, DataValue> param_up;
          CVTermList modparam;
          double p_tol = 0;
          double f_tol = 0;
          CVTermList tcv;
          map<String, DataValue> tup;
          DOMElement* child = element_sip->getFirstElementChild();
          while (child)
          {
            if (XMLString::equals(child->getTagName(), CONST_XMLCH("SearchType")))
            {
              searchtype = parseCvParam_(child->getFirstElementChild());
            }
            else if (XMLString::equals(child->getTagName(), CONST_XMLCH("AdditionalSearchParams")))
            {
              pair<CVTermList, map<String, DataValue> > as_params = parseParamGroup_(child->getChildNodes());
              sp = findSearchParameters_(as_params); // this must be renamed!!
            }
            else if (XMLString::equals(child->getTagName(), CONST_XMLCH("ModificationParams"))) // TODO @all where to store the specificities?
            {
              vector<String> fix, var;
              DOMElement* sm = child->getFirstElementChild();
              while (sm)
              {
                String residues = StringManager::convert(sm->getAttribute(CONST_XMLCH("residues")));
                bool fixedMod = false;
                XSValue::Status status;
                XSValue* val = XSValue::getActualValue(sm->getAttribute(CONST_XMLCH("fixedMod")), XSValue::dt_boolean, status);
                if (status == XSValue::st_Init)
                {
                  fixedMod = val->fData.fValue.f_bool;
                }
                delete val;

//                double massDelta = 0;
//                try
//                {
//                  massDelta = boost::lexical_cast<double>(StringManager::convert(sm->getAttribute(CONST_XMLCH("massDelta"))));
//                }
//                catch (...)
//                {
//                    OPENMS_LOG_ERROR << "Could not cast ModificationParam massDelta from " << StringManager::convert(sm->getAttribute(CONST_XMLCH("massDelta")));
//                }

                String mname;
                CVTermList specificity_rules;
                DOMElement* sub = sm->getFirstElementChild();
                while (sub)
                {
                  if (XMLString::equals(sub->getTagName(), CONST_XMLCH("cvParam")))
                  {
                    mname = StringManager::convert(sub->getAttribute(CONST_XMLCH("name")));
                   if (mname == "unknown modification")
                   {
                     // e.g. <cvParam cvRef="MS" accession="MS:1001460" name="unknown modification" value="N-Glycan"/>
                     mname = StringManager::convert(sub->getAttribute(CONST_XMLCH("value")));
                   }
                  }
                  else if (XMLString::equals(sub->getTagName(), CONST_XMLCH("SpecificityRules")))
                  {
                    specificity_rules.consumeCVTerms(parseParamGroup_(sub->getChildNodes()).first.getCVTerms());
                    // let's press them where all other SearchengineAdapters press them in
                  }
                  else
                  {
                    OPENMS_LOG_ERROR << "Misplaced information in 'ModificationParams' ignored." << endl;
                  }
                  sub = sub->getNextElementSibling();
                }

                if (!mname.empty())
                {
                  String mod;
                  String r = (residues != ".") ? residues : "";

                  if (!specificity_rules.empty())
                  {
                    for (map<String, vector<CVTerm> >::const_iterator spci = specificity_rules.getCVTerms().begin(); spci != specificity_rules.getCVTerms().end(); ++spci)
                    {
                      if (spci->second.front().getAccession() == "MS:1001189")  // nterm
                      {
                        const ResidueModification* m = ModificationsDB::getInstance()->getModification(mname, r, ResidueModification::N_TERM);
                        mod = m->getFullId();
                      }
                      else if (spci->second.front().getAccession() == "MS:1001190")  // cterm
                      {
                        const ResidueModification* m = ModificationsDB::getInstance()->getModification(mname, r, ResidueModification::C_TERM);
                        mod = m->getFullId();
                      }
                      else if (spci->second.front().getAccession() == "MS:1002057")  // protein nterm
                      {
                        const ResidueModification* m = ModificationsDB::getInstance()->getModification(mname,  r, ResidueModification::PROTEIN_N_TERM);
                        mod = m->getFullId();
                      }
                      else if (spci->second.front().getAccession() == "MS:1002058")  // protein cterm
                      {
                        const ResidueModification* m = ModificationsDB::getInstance()->getModification(mname,  r, ResidueModification::PROTEIN_C_TERM);
                        mod = m->getFullId();
                      }
                    }
                  }
                  else  // anywhere
                  {
                    const ResidueModification* m = ModificationsDB::getInstance()->getModification(mname, r);
                    mod = m->getFullId();
                  }

                  if (fixedMod)
                  {
                    fix.push_back(mod);
                  }
                  else
                  {
                    var.push_back(mod);
                  }
                }
                sm = sm->getNextElementSibling();
              }
              sp.fixed_modifications = fix;
              sp.variable_modifications = var;
            }
            else if (XMLString::equals(child->getTagName(), CONST_XMLCH("Enzymes"))) // TODO @all : where store multiple enzymes for one identificationrun?
            {
              DOMElement* enzyme = child->getFirstElementChild(); //Enzyme elements
              while (enzyme)
              {
                int missedCleavages = -1;
                try
                {
                  missedCleavages = StringManager::convert(enzyme->getAttribute(CONST_XMLCH("missedCleavages"))).toInt();
                }
                catch (exception& e)
                {
                  OPENMS_LOG_WARN << "Search engine enzyme settings for 'missedCleavages' unreadable: " << e.what() << StringManager::convert(enzyme->getAttribute(CONST_XMLCH("missedCleavages"))) << endl;
                }
                if (missedCleavages < 0)
                {
                  OPENMS_LOG_WARN << "missedCleavages has a negative value. Assuming unlimited and setting it to 1000." << endl;
                  missedCleavages = 1000;
                }
                sp.missed_cleavages = missedCleavages;

//                String semiSpecific = StringManager::convert(enzyme->getAttribute(CONST_XMLCH("semiSpecific"))); //xsd:boolean
//                String cTermGain = StringManager::convert(enzyme->getAttribute(CONST_XMLCH("cTermGain")));
//                String nTermGain = StringManager::convert(enzyme->getAttribute(CONST_XMLCH("nTermGain")));
//                int minDistance = -1;
//                try
//                {
//                  minDistance = StringManager::convert(enzyme->getAttribute(CONST_XMLCH("minDistance"))).toInt();
//                }
//                catch (...)
//                {
//                    OPENMS_LOG_WARN << "Search engine settings for 'minDistance' unreadable." << endl;
//                }
                enzymename = "UNKNOWN";
                DOMElement* sub = enzyme->getFirstElementChild();
                while (sub)
                {
                  //SiteRegex unstorable just now
                  if (XMLString::equals(sub->getTagName(), CONST_XMLCH("EnzymeName")))
                  {
                    set<String> enzymes_terms;
                    cv_.getAllChildTerms(enzymes_terms, "MS:1001045"); // cleavage agent name
                    pair<CVTermList, map<String, DataValue> > params = parseParamGroup_(sub->getChildNodes());
                    for (map<String, vector<CVTerm> >::const_iterator it = params.first.getCVTerms().begin(); it != params.first.getCVTerms().end(); ++it)
                    {
                      if (enzymes_terms.find(it->first) != enzymes_terms.end())
                      {
                        enzymename = it->second.front().getName();
                      }
                      else
                      {
                        OPENMS_LOG_WARN << "Additional parameters for enzyme settings not readable." << endl;
                      }
                    }
                  }
                  sub = sub->getNextElementSibling();
                }
                if (ProteaseDB::getInstance()->hasEnzyme(enzymename))
                {
                  sp.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzymename));
                }
                enzyme = enzyme->getNextElementSibling();
              }
            }
            else if (XMLString::equals(child->getTagName(), CONST_XMLCH("FragmentTolerance")))
            {
              pair<CVTermList, map<String, DataValue> > params = parseParamGroup_(child->getChildNodes());
              //+- take the numerically greater
              for (map<String, vector<CVTerm> >::const_iterator it = params.first.getCVTerms().begin(); it != params.first.getCVTerms().end(); ++it)
              {
                f_tol = max(f_tol, it->second.front().getValue().toString().toDouble());
                sp.fragment_mass_tolerance = f_tol;
                if (it->second.front().getUnit().name == "parts per million" )
                {
                  sp.fragment_mass_tolerance_ppm = true;
                }
              }
            }
            else if (XMLString::equals(child->getTagName(), CONST_XMLCH("ParentTolerance")))
            {
              pair<CVTermList, map<String, DataValue> > params = parseParamGroup_(child->getChildNodes());
              //+- take the numerically greater
              for (map<String, vector<CVTerm> >::const_iterator it = params.first.getCVTerms().begin(); it != params.first.getCVTerms().end(); ++it)
              {
                p_tol = max(p_tol, it->second.front().getValue().toString().toDouble());
                sp.precursor_mass_tolerance = p_tol;
                if (it->second.front().getUnit().name == "parts per million" )
                {
                  sp.precursor_mass_tolerance_ppm = true;
                }

              }
            }
            else if (XMLString::equals(child->getTagName(), CONST_XMLCH("Threshold")))
            {
              pair<CVTermList, map<String, DataValue> > params = parseParamGroup_(child->getChildNodes());
              tcv = params.first;
              tup = params.second;
            }
            child = child->getNextElementSibling();
            //      <DatabaseFilters> omitted for now, not reflectable by our member structures
            //      <DatabaseTranslation> omitted for now, not reflectable by our member structures
            //      <MassTable> omitted for now, not reflectable by our member structures
          }
          SpectrumIdentificationProtocol temp_struct = {searchtype, enzymename, param_cv, param_up, modparam, p_tol, f_tol, tcv, tup};
          sp_map_.insert(make_pair(id, temp_struct));

          double thresh = 0.0;
          bool use_thresh = false;
          set<String> threshold_terms;
          cv_.getAllChildTerms(threshold_terms, "MS:1002482"); //statistical threshold
          for (map<String, vector<OpenMS::CVTerm> >::const_iterator thit = tcv.getCVTerms().begin(); thit != tcv.getCVTerms().end(); ++thit)
          {
            if (threshold_terms.find(thit->first) != threshold_terms.end())
            {
              if (thit->first != "MS:1001494") // no threshold
              {
                thresh = thit->second.front().getValue().toString().toDouble(); // cast fix needed as DataValue is init with XercesString
                use_thresh = true;
                break;
              }
              else
              {
                break;
              }
            }
          }

          String search_engine, search_engine_version;
          for (map<String, SpectrumIdentification>::const_iterator si_it = si_map_.begin(); si_it != si_map_.end(); ++si_it)
          {
            if (si_it->second.spectrum_identification_protocol_ref == id)
            {
              search_engine = as_map_[swr].name;
              search_engine_version = as_map_[swr].version;
//              String identi = search_engine+"_"+si_pro_map_[si_it->second.spectrum_identification_list_ref]->getDateTime().getDate()+"T"
//                      +si_pro_map_[si_it->second.spectrum_identification_list_ref]->getDateTime().getTime();
//              pro_id_->at(si_pro_map_[si_it->second.spectrum_identification_list_ref]).setIdentifier(identi);
              pro_id_->at(si_pro_map_[si_it->second.spectrum_identification_list_ref]).setSearchEngine(search_engine);
              pro_id_->at(si_pro_map_[si_it->second.spectrum_identification_list_ref]).setSearchEngineVersion(search_engine_version);
              sp.db = pro_id_->at(si_pro_map_[si_it->second.spectrum_identification_list_ref]).getSearchParameters().db; // was previously set, but main parts of sp are set here
              sp.db_version = pro_id_->at(si_pro_map_[si_it->second.spectrum_identification_list_ref]).getSearchParameters().db_version; // was previously set, but main parts of sp are set here
              pro_id_->at(si_pro_map_[si_it->second.spectrum_identification_list_ref]).setSearchParameters(sp);
              if (use_thresh)
              {
                pro_id_->at(si_pro_map_[si_it->second.spectrum_identification_list_ref]).setSignificanceThreshold(thresh);
              }
            }
          }
        }
      }
    }

    void MzIdentMLDOMHandler::parseInputElements_(DOMNodeList* inputElements)
    {
      const  XMLSize_t node_count = inputElements->getLength();
      for (XMLSize_t c = 0; c < node_count; ++c)
      {
        DOMNode* current_in = inputElements->item(c);
        if (current_in->getNodeType() && // true is not NULL
            current_in->getNodeType() == DOMNode::ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_in = dynamic_cast<xercesc::DOMElement*>(current_in);

          String id = StringManager::convert(element_in->getAttribute(CONST_XMLCH("id")));
          String location = StringManager::convert(element_in->getAttribute(CONST_XMLCH("location")));

          if (XMLString::equals(element_in->getTagName(), CONST_XMLCH("SpectraData")))
          {
            //      <FileFormat> omitted for now, not reflectable by our member structures
            //      <SpectrumIDFormat> omitted for now, not reflectable by our member structures
            sd_map_.insert(make_pair(id, location));
          }
          else if (XMLString::equals(element_in->getTagName(), CONST_XMLCH("SourceFile")))
          {
            //      <FileFormat> omitted for now, not reflectable by our member structures
            sr_map_.insert(make_pair(id, location));
          }
          else if (XMLString::equals(element_in->getTagName(), CONST_XMLCH("SearchDatabase")))
          {
            //      <FileFormat> omitted for now, not reflectable by our member structures
            DateTime releaseDate;
//            releaseDate.set(StringManager::convert(element_in->getAttribute(CONST_XMLCH("releaseDate"))));
            String version = StringManager::convert(element_in->getAttribute(CONST_XMLCH("version")));
            String dbname;
            DOMElement* element_dbn = element_in->getFirstElementChild();
            while (element_dbn)
            {
              if (XMLString::equals(element_dbn->getTagName(), CONST_XMLCH("DatabaseName")))
              {
                DOMElement* databasename_param = element_dbn->getFirstElementChild();
                while (databasename_param)
                {
                  if (XMLString::equals(databasename_param->getTagName(), CONST_XMLCH("cvParam")))
                  {
                    CVTerm param = parseCvParam_(databasename_param);
                    dbname = param.getValue();
                  }
                  else if (XMLString::equals(databasename_param->getTagName(), CONST_XMLCH("userParam")))
                  {
                    pair<String, DataValue> param = parseUserParam_(databasename_param);
                    // issue #7099: mzID might have missing "value" for this element
                    // in this case, just use the "name"
                    dbname = param.second.isEmpty() ? dbname = param.first : param.second.toString();
                  }
                  databasename_param = databasename_param->getNextElementSibling();
                }
                //each SearchDatabase element may have one DatabaseName, each DatabaseName only one param
              }
              element_dbn = element_dbn->getNextElementSibling();
            }
            if (dbname.empty())
            {
              OPENMS_LOG_WARN << "No DatabaseName element found, use results at own risk." << endl;
              dbname = "unknown";
            }
            DatabaseInput temp_struct = {dbname, location, version, releaseDate};
            db_map_.insert(make_pair(id, temp_struct));
          }
        }
      }
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationListElements_(DOMNodeList* spectrumIdentificationListElements)
    {
      const  XMLSize_t node_count = spectrumIdentificationListElements->getLength();
      for (XMLSize_t c = 0; c < node_count; ++c)
      {
        DOMNode* current_lis = spectrumIdentificationListElements->item(c);
        if (current_lis->getNodeType() && // true is not NULL
            current_lis->getNodeType() == DOMNode::ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_lis = dynamic_cast<xercesc::DOMElement*>(current_lis);
          String id = StringManager::convert(element_lis->getAttribute(CONST_XMLCH("id")));
//          String name = StringManager::convert(element_res->getAttribute(CONST_XMLCH("name")));

          DOMElement* element_res = element_lis->getFirstElementChild();
          while (element_res)
          {
            if (XMLString::equals(element_res->getTagName(), CONST_XMLCH("SpectrumIdentificationResult")))
            {
              String spectra_data_ref = StringManager::convert(element_res->getAttribute(CONST_XMLCH("spectraData_ref"))); //ref to the sourcefile, could be useful but now nowhere to store
              String spectrumID = StringManager::convert(element_res->getAttribute(CONST_XMLCH("spectrumID")));
              pair<CVTermList, map<String, DataValue> > params = parseParamGroup_(element_res->getChildNodes());

              if (xl_ms_search_) // XL-MS data has a different structure (up to 4 spectrum identification items for the same PSM)
              {
                std::multimap<String, int> xl_val_map;
                std::set<String> xl_val_set;
                int index_counter = 0;
                DOMElement* sii = element_res->getFirstElementChild();

                // loop over all SIIs of a spectrum and group together the SIIs belonging to the same cross-link spectrum match
                while (sii)
                {
                  if (XMLString::equals(sii->getTagName(), CONST_XMLCH("SpectrumIdentificationItem")))
                  {
                    DOMNodeList* sii_cvp = sii->getElementsByTagName(CONST_XMLCH("cvParam"));
                    const  XMLSize_t cv_count = sii_cvp->getLength();
                    for (XMLSize_t i = 0; i < cv_count; ++i)
                    {
                      DOMElement* element_sii_cvp = dynamic_cast<xercesc::DOMElement*>(sii_cvp->item(i));
                      if (XMLString::equals(element_sii_cvp->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1002511"))) // cross-link spectrum identification item
                      {
                        String xl_val = StringManager::convert(element_sii_cvp->getAttribute(CONST_XMLCH("value")));
                        xl_val_map.insert(make_pair(xl_val, index_counter));
                        xl_val_set.insert(xl_val);
                      }
                    }
                  }
                  sii = sii->getNextElementSibling();
                  ++index_counter;
                }

                // fix for label-free mono-links
                // those only have one SII and no "cross-link spectrum identification item" value
                if (xl_val_set.empty())
                {
                  xl_val_set.insert("0");
                  xl_val_map.insert(make_pair("0", 0));
                }
                for (set<String>::const_iterator set_it = xl_val_set.begin(); set_it != xl_val_set.end(); ++set_it)
                {
                  parseSpectrumIdentificationItemSetXLMS(set_it, xl_val_map, element_res, spectrumID);
                }
                pep_id_->back().setIdentifier(pro_id_->at(si_pro_map_[id]).getIdentifier());
              }
              else // general case
              {
                pep_id_->push_back(PeptideIdentification());
                pep_id_->back().setHigherScoreBetter(false); //either a q-value or an e-value, only if neither available there will be another
                pep_id_->back().setSpectrumReference( spectrumID);  // SpectrumIdentificationResult attribute spectrumID is taken from the mz_file and should correspond to MSSpectrum.nativeID, thus spectrum_reference will serve as reference. As the format of the 'reference' widely varies from vendor to vendor, spectrum_reference as string will serve best, indices are not recommended.

                //fill pep_id_->back() with content
                DOMElement* parent = dynamic_cast<xercesc::DOMElement*>(element_res->getParentNode());
                String sil = StringManager::convert(parent->getAttribute(CONST_XMLCH("id")));

                DOMElement* child = element_res->getFirstElementChild();
                while (child)
                {
                  if (XMLString::equals(child->getTagName(), CONST_XMLCH("SpectrumIdentificationItem")))
                  {
                    parseSpectrumIdentificationItemElement_(child, pep_id_->back(), sil);
                  }
                  child = child->getNextElementSibling();
                }

              } // end of "not-XLMS-results"

              // TODO @mths: setSignificanceThreshold, but from where?

              pep_id_->back().setIdentifier(pro_id_->at(si_pro_map_[id]).getIdentifier());

              pep_id_->back().sortByRank();

              //adopt cv s
              for (map<String, vector<CVTerm> >::const_iterator cvit =  params.first.getCVTerms().begin(); cvit != params.first.getCVTerms().end(); ++cvit)
              {
                // check for retention time or scan time entry
                /* N.B.: MzIdentML does not impose the requirement to store
                   'redundant' data (e.g. RT) as the identified spectrum is
                   unambiguously referenceable by the spectrumID (OpenMS
                   internally spectrum_reference) and hence such data can be
                   looked up in the mz file. For convenience, and as OpenMS
                   relies on the smallest common denominator to reference a
                   spectrum (RT/precursor MZ), we provide functionality to amend
                   RT data to identifications and support reading such from mzid
                */
                if (cvit->first == "MS:1000894" || cvit->first == "MS:1000016") //TODO use subordinate terms which define units
                {
                  double rt = cvit->second.front().getValue().toString().toDouble();
                  if (cvit->second.front().getUnit().accession == "UO:0000031")  // minutes
                  {
                    rt *= 60.0;
                  }
                  pep_id_->back().setRT(rt);
                }
                else
                {
                  pep_id_->back().setMetaValue(cvit->first, cvit->second.front().getValue()); // TODO? all DataValues - are there more then one, my guess is this is overdesigned
                }
              }
              //adopt up s
              for (map<String, DataValue>::const_iterator upit = params.second.begin(); upit != params.second.end(); ++upit)
              {
                pep_id_->back().setMetaValue(upit->first, upit->second);
              }
              if (pep_id_->back().getRT() != pep_id_->back().getRT())
              {
                OPENMS_LOG_WARN << "No retention time found for 'SpectrumIdentificationResult'" << endl;
              }
            }
            element_res = element_res->getNextElementSibling();
          }
        }
      }
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationItemSetXLMS(set<String>::const_iterator set_it, std::multimap<String, int> xl_val_map, DOMElement* element_res, const String& spectrumID)
    {
      // each value in the set corresponds to one PeptideIdentification object
      std::pair <std::multimap<String, int>::iterator, std::multimap<String, int>::iterator> range;
      range = xl_val_map.equal_range(*set_it);
      DOMNodeList* siis = element_res->getElementsByTagName(CONST_XMLCH("SpectrumIdentificationItem"));

      DOMElement* parent = dynamic_cast<xercesc::DOMElement*>(element_res->getParentNode());
      String spectrumIdentificationList_ref = StringManager::convert(parent->getAttribute(CONST_XMLCH("id")));

      // initialize all needed values, extract them one by one in vectors, e.g. using max and min to determine which are heavy which light
      // get peptide id, that way determine donor, acceptor (alpha, beta)
      std::vector<String> peptides;
      double score = -1;
      std::vector<double> exp_mzs;
      std::vector<double> RTs;
      int rank = 0;
      int charge = 0;
      vector<PeptideHit::PeakAnnotation> frag_annotations;

      double xcorrx = 0;
      double xcorrc = 0;
      double matchodds = 0;
      double intsum = 0;
      double wTIC = 0;
      vector< vector<String> > userParamNameLists;
      vector< vector<String> > userParamValueLists;
      vector< vector<String> > userParamUnitLists;

      for (std::multimap<String, int>::iterator it=range.first; it!=range.second; ++it)
      {
        DOMElement* cl_sii = dynamic_cast<xercesc::DOMElement*>(siis->item(it->second));
        // Attributes
        String peptide = StringManager::convert(cl_sii->getAttribute(CONST_XMLCH("peptide_ref")));
        peptides.push_back(peptide);
        double exp_mz =StringManager::convert(cl_sii->getAttribute(CONST_XMLCH("experimentalMassToCharge"))).toDouble();
        exp_mzs.push_back(exp_mz);

        if (rank == 0)
        {
          rank = StringManager::convert(cl_sii->getAttribute(CONST_XMLCH("rank"))).toInt();
        }
        if (charge == 0)
        {
          charge = StringManager::convert(cl_sii->getAttribute(CONST_XMLCH("chargeState"))).toInt();
        }

        // CVs
        DOMNodeList* sii_cvp = cl_sii->getElementsByTagName(CONST_XMLCH("cvParam"));
        const  XMLSize_t cv_count = sii_cvp->getLength();
        for (XMLSize_t i = 0; i < cv_count; ++i)
        {
          DOMElement* element_sii_cvp = dynamic_cast<xercesc::DOMElement*>(sii_cvp->item(i));
          if (XMLString::equals(element_sii_cvp->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1002681"))) // OpenXQuest:combined score
          {
            score = StringManager::convert(element_sii_cvp->getAttribute(CONST_XMLCH("value"))).toDouble();
          }
          else if (XMLString::equals(element_sii_cvp->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1002682"))) // OpenXQuest: xcorr common
          {
            xcorrx = StringManager::convert(element_sii_cvp->getAttribute(CONST_XMLCH("value"))).toDouble();
          }
          else if (XMLString::equals(element_sii_cvp->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1002683"))) // OpenXQuest: xcorr xlink
          {
            xcorrc = StringManager::convert(element_sii_cvp->getAttribute(CONST_XMLCH("value"))).toDouble();
          }
          else if (XMLString::equals(element_sii_cvp->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1002684"))) // OpenXQuest: match-odds
          {
            matchodds = StringManager::convert(element_sii_cvp->getAttribute(CONST_XMLCH("value"))).toDouble();
          }
          else if (XMLString::equals(element_sii_cvp->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1002685"))) // OpenXQuest: intsum
          {
            intsum = StringManager::convert(element_sii_cvp->getAttribute(CONST_XMLCH("value"))).toDouble();
          }
          else if (XMLString::equals(element_sii_cvp->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1002686"))) // OpenXQuest: wTIC
          {
            wTIC = StringManager::convert(element_sii_cvp->getAttribute(CONST_XMLCH("value"))).toDouble();
          }
          else if (XMLString::equals(element_sii_cvp->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1003024"))) // OpenPepXL:score
          {
            score = StringManager::convert(element_sii_cvp->getAttribute(CONST_XMLCH("value"))).toDouble();
          }
          else if (XMLString::equals(element_sii_cvp->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1000894"))) // retention time
          {
            double RT = StringManager::convert(element_sii_cvp->getAttribute(CONST_XMLCH("value"))).toDouble();
            RTs.push_back(RT);
          }
        }

        // userParams
        vector<String> userParamNames;
        vector<String> userParamValues;
        vector<String> userParamUnits;

        DOMNodeList* sii_up = cl_sii->getElementsByTagName(CONST_XMLCH("userParam"));
        const  XMLSize_t up_count = sii_up->getLength();
        for (XMLSize_t i = 0; i < up_count; ++i)
        {
          DOMElement* element_sii_up = dynamic_cast<xercesc::DOMElement*>(sii_up->item(i));
          userParamNames.push_back(StringManager::convert(element_sii_up->getAttribute(CONST_XMLCH("name"))));
          userParamValues.push_back(StringManager::convert(element_sii_up->getAttribute(CONST_XMLCH("value"))));
          userParamUnits.push_back(StringManager::convert(element_sii_up->getAttribute(CONST_XMLCH("unitName"))));
        }
        userParamNameLists.push_back(userParamNames);
        userParamValueLists.push_back(userParamValues);
        userParamUnitLists.push_back(userParamUnits);

        // Fragmentation, does not matter where to get them. Look for them as long as the vector is empty
        if (frag_annotations.empty())
        {
          DOMNodeList* frag_element_list = cl_sii->getElementsByTagName(CONST_XMLCH("Fragmentation"));

          if (frag_element_list->getLength() > 0)
          {
            DOMElement* frag_element = dynamic_cast<xercesc::DOMElement*>(frag_element_list->item(0));
            DOMNodeList* ion_types = frag_element->getElementsByTagName(CONST_XMLCH("IonType"));
            const  XMLSize_t ion_type_count = ion_types->getLength();
            for (XMLSize_t i = 0; i < ion_type_count; ++i)
            {
              DOMElement* ion_type_element = dynamic_cast<xercesc::DOMElement*>(ion_types->item(i));
              int ion_charge = StringManager::convert(ion_type_element ->getAttribute(CONST_XMLCH("charge"))).toInt();
              vector<String> indices;
              vector<String> positions;
              vector<String> intensities;
              vector<String> chains;
              vector<String> categories;
              String frag_type;
              String loss = "";

              StringManager::convert(ion_type_element ->getAttribute(CONST_XMLCH("index"))).split(" ", indices);

              DOMNodeList* frag_arrays = ion_type_element ->getElementsByTagName(CONST_XMLCH("FragmentArray"));
              const XMLSize_t frag_array_count = frag_arrays->getLength();
              for (XMLSize_t f = 0; f < frag_array_count; ++f)
              {
                DOMElement* frag_array_element = dynamic_cast<xercesc::DOMElement*>(frag_arrays->item(f));
                if ( XMLString::equals(frag_array_element->getAttribute(CONST_XMLCH("measure_ref")), CONST_XMLCH("Measure_mz")))
                {
                  StringManager::convert(frag_array_element->getAttribute(CONST_XMLCH("values"))).split(" ", positions);
                }
                if ( XMLString::equals(frag_array_element->getAttribute(CONST_XMLCH("measure_ref")), CONST_XMLCH("Measure_int")))
                {
                  StringManager::convert(frag_array_element->getAttribute(CONST_XMLCH("values"))).split(" ", intensities);
                }
              }

              DOMNodeList* userParams = ion_type_element->getElementsByTagName(CONST_XMLCH("userParam"));
              const XMLSize_t userParam_count = userParams->getLength();
              for (XMLSize_t u = 0; u < userParam_count; ++u)
              {
                DOMElement* userParam_element = dynamic_cast<xercesc::DOMElement*>(userParams->item(u));
                if ( XMLString::equals(userParam_element ->getAttribute(CONST_XMLCH("name")), CONST_XMLCH("cross-link_chain")))
                {
                  StringManager::convert(userParam_element ->getAttribute(CONST_XMLCH("value"))).split(" ", chains);
                }
                if ( XMLString::equals(userParam_element ->getAttribute(CONST_XMLCH("name")), CONST_XMLCH("cross-link_ioncategory")))
                {
                  StringManager::convert(userParam_element ->getAttribute(CONST_XMLCH("value"))).split(" ", categories);
                }
              }

              DOMNodeList* cvts = ion_type_element->getElementsByTagName(CONST_XMLCH("cvParam"));
              const XMLSize_t cvt_count = cvts->getLength();
              for (XMLSize_t cvt = 0; cvt < cvt_count; ++cvt)
              {
                DOMElement* cvt_element = dynamic_cast<xercesc::DOMElement*>(cvts->item(cvt));

                // Standard ions
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001229"))) // frag: a ion
                {
                  frag_type = "a";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001224"))) // frag: b ion
                {
                  frag_type = "b";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001231"))) // frag: c ion
                {
                  frag_type = "c";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001228"))) // frag: x ion
                {
                  frag_type = "x";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001220"))) // frag: y ion
                {
                  frag_type = "y";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001230"))) // frag: z ion
                {
                  frag_type = "z";
                }

                // Ions with H2O losses
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001234"))) // frag: a ion - H2O
                {
                  frag_type = "a";
                  loss = "-H2O";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001222"))) // frag: b ion - H20
                {
                  frag_type = "b";
                  loss = "-H2O";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001515"))) // frag: c ion - H20
                {
                  frag_type = "c";
                  loss = "-H2O";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001519"))) // frag: x ion - H20
                {
                  frag_type = "x";
                  loss = "-H2O";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001223"))) // frag: y ion - H20
                {
                  frag_type = "y";
                  loss = "-H2O";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001517"))) // frag: z ion - H20
                {
                  frag_type = "z";
                  loss = "-H2O";
                }

                // Ions with NH3 losses
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001235"))) // frag: a ion - NH3
                {
                  frag_type = "a";
                  loss = "-NH3";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001232"))) // frag: b ion - NH3
                {
                  frag_type = "b";
                  loss = "-NH3";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001516"))) // frag: c ion - NH3
                {
                  frag_type = "c";
                  loss = "-NH3";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001520"))) // frag: x ion - NH3
                {
                  frag_type = "x";
                  loss = "-NH3";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001233"))) // frag: y ion - NH3
                {
                  frag_type = "y";
                  loss = "-NH3";
                }
                if (XMLString::equals(cvt_element->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1001518"))) // frag: z ion - NH3
                {
                  frag_type = "z";
                  loss = "-NH3";
                }
              }

              for (Size s = 0; s < indices.size(); ++s)
              {
                String annotation= "[" + chains[s] + "|" + categories[s]  + "$" + frag_type + indices[s] + loss + "]";

                PeptideHit::PeakAnnotation frag_anno;
                frag_anno.charge = ion_charge;
                frag_anno.mz = positions[s].toDouble();
                frag_anno.intensity = intensities[s].toDouble();
                frag_anno.annotation = annotation;
                frag_annotations.push_back(frag_anno);
              }
            }
          }
        }
      }

      // Generate and fill PeptideIdentification
      vector<Size> light;
      vector<Size> heavy;
      double MZ_light = *std::min_element(exp_mzs.begin(), exp_mzs.end());
      double MZ_heavy = *std::max_element(exp_mzs.begin(), exp_mzs.end());

      // are the cross-links labeled?
      bool labeled = MZ_light != MZ_heavy;

      for (Size i = 0; i < exp_mzs.size(); ++i)
      {
        if (MZ_light == exp_mzs[i])
        {
          light.push_back(i);
        }
        else
        {
          heavy.push_back(i);
        }
      }
      double RT_light = RTs[light[0]];
      double RT_heavy = RT_light;

      if (labeled)
      {
        RT_heavy = RTs[heavy[0]];
      }

      vector<Size> alpha;
      vector<Size> beta;
      for (Size i = 0; i < peptides.size(); ++i)
      {
        try
        {
          String donor_pep = xl_id_donor_map_.at(peptides[i]);      // map::at throws an out-of-range
          alpha.push_back(i);
        }
        catch (const std::out_of_range& /*oor*/)
        {
          beta.push_back(i);
        }
      }
      String xl_type = "mono-link";

      SignedSize alpha_pos = xl_donor_pos_map_.at(xl_id_donor_map_.at(peptides[alpha[0]]));
      vector<String> spectrumIDs;
      spectrumID.split(",", spectrumIDs);

      if (alpha.size() == beta.size()) // if a beta exists at all, it must be cross-link. But they should also each have the same number of SIIs in the mzid
      {
        xl_type = "cross-link";
      }
      else
      {
        try
        {
          String donor_val = xl_id_donor_map_.at(peptides[alpha[0]]);      // map::at throws an out-of-range
          String acceptor_val = xl_id_acceptor_map_.at(peptides[alpha[0]]);
          if (donor_val  == acceptor_val) // if donor and acceptor on the same peptide belong to the same cross-link, it must be a loop-link
          {
            xl_type = "loop-link";
          }
        }
        catch (const std::out_of_range& /*oor*/)
        {
            // do nothing. Must be a mono-link, which is already set
        }
      }

      PeptideIdentification current_pep_id;
      current_pep_id.setRT(RT_light);
      current_pep_id.setMZ(MZ_light);
      current_pep_id.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, spectrumID);
      current_pep_id.setScoreType(Constants::UserParam::OPENPEPXL_SCORE);
      current_pep_id.setHigherScoreBetter(true);

      vector<PeptideHit> phs;
      PeptideHit ph_alpha;
      ph_alpha.setSequence((*pep_map_.find(peptides[alpha[0]])).second);
      ph_alpha.setCharge(charge);
      ph_alpha.setScore(score);
      ph_alpha.setRank(rank);
      ph_alpha.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, spectrumIDs[0]);
      ph_alpha.setMetaValue("xl_chain", "MS:1002509"); // donor

      if (labeled)
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_RT, RT_heavy);
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_MZ, MZ_heavy);
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_HEAVY_SPEC_REF, spectrumIDs[1]);
      }

      ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE, xl_type);
      ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK, rank);

      ph_alpha.setMetaValue("OpenPepXL:xcorr xlink",xcorrx);
      ph_alpha.setMetaValue("OpenPepXL:xcorr common", xcorrc);
      ph_alpha.setMetaValue("OpenPepXL:match-odds", matchodds);
      ph_alpha.setMetaValue("OpenPepXL:intsum", intsum);
      ph_alpha.setMetaValue("OpenPepXL:wTIC", wTIC);

      vector<String> userParamNames_alpha = userParamNameLists[alpha[0]];
      vector<String> userParamValues_alpha = userParamValueLists[alpha[0]];
      vector<String> userParamUnits_alpha = userParamUnitLists[alpha[0]];

      for (Size i = 0; i < userParamNames_alpha.size(); ++i)
      {
        ph_alpha.setMetaValue(userParamNames_alpha[i], XMLHandler::fromXSDString(userParamUnits_alpha[i], userParamValues_alpha[i]));
      }

      ph_alpha.setPeakAnnotations(frag_annotations);

      if (xl_type == "loop-link")
      {
        Size pos2 = xl_acceptor_pos_map_.at(xl_id_acceptor_map_.at(peptides[alpha[0]]));
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2, pos2);
      }

      if (xl_type != "mono-link")
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD, xl_mod_map_.at(peptides[alpha[0]]));
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MASS,DataValue(xl_mass_map_.at(peptides[alpha[0]])));
      }
      else if ( xl_mod_map_.find(peptides[alpha[0]]) != xl_mod_map_.end() )
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_MOD, xl_mod_map_.at(peptides[alpha[0]]));
      }

      // correction for terminal modifications
      if (alpha_pos == -1)
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1, ++alpha_pos);
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA, "N_TERM");
      }
      else if (alpha_pos == static_cast<SignedSize>(ph_alpha.getSequence().size()))
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1, --alpha_pos);
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA, "C_TERM");
      }
      else
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1, alpha_pos);
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA, "ANYWHERE");
      }

      if (xl_type == "cross-link")
      {
        PeptideHit ph_beta;
        SignedSize beta_pos = xl_acceptor_pos_map_.at(xl_id_acceptor_map_.at(peptides[beta[0]]));

        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE, (*pep_map_.find(peptides[beta[0]])).second.toString());
        ph_beta.setSequence((*pep_map_.find(peptides[beta[0]])).second);
        ph_beta.setCharge(charge);
        ph_beta.setScore(score);
        ph_beta.setRank(rank);
        ph_beta.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, spectrumIDs[0]);
        ph_beta.setMetaValue("xl_chain", "MS:1002510"); // receiver

        // correction for terminal modifications
        if (beta_pos == -1)
        {
          ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2, ++beta_pos);
          ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA, "N_TERM");
        }
        else if (beta_pos == static_cast<SignedSize>(ph_beta.getSequence().size()))
        {
          ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2, --beta_pos);
          ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA, "C_TERM");
        }
        else
        {
          ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2, beta_pos);
          ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA, "ANYWHERE");
        }
        phs.push_back(ph_alpha);
        phs.push_back(ph_beta);
      }
      else
      {
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA, "ANYWHERE");
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2, "-");
        ph_alpha.setMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE, "-");
        phs.push_back(ph_alpha);
      }

      std::vector<String> unique_peptides;
      unique_peptides.push_back(peptides[alpha[0]]);
      if (!beta.empty())
      {
        unique_peptides.push_back(peptides[beta[0]]);
      }

      for (Size pep = 0; pep < unique_peptides.size(); ++pep)
      {

        String peptide_ref = unique_peptides[pep];

        // Debug output
//        cout << "peptides: ";
//        for (Size k = 0; k < peptides.size(); ++k)
//        {
//          cout << "nr." << k << ": " << peptides[k] << "\t";
//        }
//        cout << endl << "unique peptides: ";
//        for (Size k = 0; k < unique_peptides.size(); ++k)
//        {
//          cout << "nr." << k << ": " << unique_peptides[k] << "\t";
//        }
//        cout << endl << "alpha: " << alpha[0];
//        if (beta.size() > 0)
//        {
//          cout << "\t beta: " << beta[0];
//        }
//        cout << endl;
//        cout << "xl_type: " << xl_type << "\tphs_size: " << phs.size() << endl;
//        cout << "phs0: " << phs[0].getSequence().toString() << endl;
//        if (phs.size() > 1)
//        {
//          cout << "phs1: " << phs[1].getSequence().toString() << endl;
//        }

        //connect the PeptideHit with PeptideEvidences (for AABefore/After) and subsequently with DBSequence (for ProteinAccession)
        pair<multimap<String, String>::iterator, multimap<String, String>::iterator> pev_its;
        pev_its = p_pv_map_.equal_range(peptide_ref);
        for (multimap<String, String>::iterator pev_it = pev_its.first; pev_it != pev_its.second; ++pev_it)
        {
          bool idec = false;
          OpenMS::PeptideEvidence pev;
          if (pe_ev_map_.find(pev_it->second) != pe_ev_map_.end())
          {
            MzIdentMLDOMHandler::PeptideEvidence& pv = pe_ev_map_[pev_it->second];
            if (pv.pre != '-') pev.setAABefore(pv.pre);
            if (pv.post != '-') pev.setAAAfter(pv.post);

            if (pv.start != OpenMS::PeptideEvidence::UNKNOWN_POSITION && pv.stop != OpenMS::PeptideEvidence::UNKNOWN_POSITION)
            {
              pev.setStart(pv.start);
              pev.setEnd(pv.stop);
            }

            idec = pv.idec;
            if (idec)
            {
              if (phs[pep].metaValueExists(Constants::UserParam::TARGET_DECOY))
              {
                if (phs[pep].getMetaValue(Constants::UserParam::TARGET_DECOY) != "decoy")
                {
                  phs[pep].setMetaValue(Constants::UserParam::TARGET_DECOY, "target+decoy");
                }
                else
                {
                  phs[pep].setMetaValue(Constants::UserParam::TARGET_DECOY, "decoy");
                }
              }
              else
              {
                phs[pep].setMetaValue(Constants::UserParam::TARGET_DECOY, "decoy");
              }
            }
            else
            {
              if (phs[pep].metaValueExists(Constants::UserParam::TARGET_DECOY))
              {
                if (phs[pep].getMetaValue(Constants::UserParam::TARGET_DECOY) != "target")
                {
                  phs[pep].setMetaValue(Constants::UserParam::TARGET_DECOY, "target+decoy");
                }
                else
                {
                  phs[pep].setMetaValue(Constants::UserParam::TARGET_DECOY, "target");
                }
              }
              else
              {
                phs[pep].setMetaValue(Constants::UserParam::TARGET_DECOY, "target");
              }
            }
          }

          if (pv_db_map_.find(pev_it->second) != pv_db_map_.end())
          {
            String& dpv = pv_db_map_[pev_it->second];
            DBSequence& db = db_sq_map_[dpv];
            pev.setProteinAccession(db.accession);

            if (pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).findHit(db.accession)
                == pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).getHits().end())
            { // butt ugly! TODO @ mths for ProteinInference
              pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).insertHit(ProteinHit());
              pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).getHits().back().setSequence(db.sequence);
              pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).getHits().back().setAccession(db.accession);
              if (idec)
              {
                pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).getHits().back().setMetaValue("isDecoy", "true");
              }
              else
              {
                pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).getHits().back().setMetaValue("isDecoy", "false");
              }
            }
          }
          phs[pep].addPeptideEvidence(pev);
        }
      }
      current_pep_id.setHits(phs);
      pep_id_->push_back(current_pep_id);
      pep_id_->back().sortByRank();
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationItemElement_(DOMElement* spectrumIdentificationItemElement, PeptideIdentification& spectrum_identification, String& spectrumIdentificationList_ref)
    {
      String id = StringManager::convert(spectrumIdentificationItemElement->getAttribute(CONST_XMLCH("id")));
      String name = StringManager::convert(spectrumIdentificationItemElement->getAttribute(CONST_XMLCH("name")));

      long double calculatedMassToCharge = StringManager::convert(spectrumIdentificationItemElement->getAttribute(CONST_XMLCH("calculatedMassToCharge"))).toDouble();
//      long double calculatedPI = StringManager::convert(spectrumIdentificationItemElement->getAttribute(CONST_XMLCH("calculatedPI"))).toDouble();
      int chargeState = 0;
      try
      {
        chargeState = StringManager::convert(spectrumIdentificationItemElement->getAttribute(CONST_XMLCH("chargeState"))).toInt();
      }
      catch (...)
      {
        OPENMS_LOG_WARN << "Found unreadable 'chargeState'." << endl;
      }
      long double experimentalMassToCharge = StringManager::convert(spectrumIdentificationItemElement->getAttribute(CONST_XMLCH("experimentalMassToCharge"))).toDouble();
      int rank = 0;
      try
      {
        rank = StringManager::convert(spectrumIdentificationItemElement->getAttribute(CONST_XMLCH("rank"))).toInt();
      }
      catch (...)
      {
        OPENMS_LOG_WARN << "Found unreadable PSM rank." << endl;
      }

      String peptide_ref = StringManager::convert(spectrumIdentificationItemElement->getAttribute(CONST_XMLCH("peptide_ref")));
//      String sample_ref = StringManager::convert(spectrumIdentificationItemElement->getAttribute(CONST_XMLCH("sample_ref")));
//      String massTable_ref = StringManager::convert(spectrumIdentificationItemElement->getAttribute(CONST_XMLCH("massTable_ref")));

      XSValue::Status status;
      std::unique_ptr<XSValue> val(XSValue::getActualValue(spectrumIdentificationItemElement->getAttribute(CONST_XMLCH("passThreshold")), XSValue::dt_boolean, status));
      bool pass = false;
      if (status == XSValue::st_Init)
      {
        pass = val->fData.fValue.f_bool;
      }
      // TODO @all: where to store passThreshold value? set after score type eval in pass_threshold

      long double score = 0;
      const auto& [param_cv, param_user] = parseParamGroup_(spectrumIdentificationItemElement->getChildNodes());
      set<String> q_score_terms;
      set<String> e_score_terms,e_score_tmp;
      set<String> specific_score_terms;
      cv_.getAllChildTerms(q_score_terms, "MS:1002354"); //q-value for peptides
      cv_.getAllChildTerms(e_score_terms, "MS:1001872");
      cv_.getAllChildTerms(e_score_tmp, "MS:1002353");
      e_score_terms.insert(e_score_tmp.begin(),e_score_tmp.end()); //E-value for peptides
      cv_.getAllChildTerms(specific_score_terms, "MS:1001143"); //search engine specific score for PSMs
      bool scoretype = false;
      for (map<String, vector<OpenMS::CVTerm>>::const_iterator scoreit = param_cv.getCVTerms().begin(); scoreit != param_cv.getCVTerms().end(); ++scoreit)
      {
        if (q_score_terms.find(scoreit->first) != q_score_terms.end() || scoreit->first == "MS:1002354")
        {
          if (scoreit->first != "MS:1002055") // do not use peptide-level q-values for now
          {
            score = scoreit->second.front().getValue().toString().toDouble(); // cast fix needed as DataValue is init with XercesString
            spectrum_identification.setHigherScoreBetter(false);
            spectrum_identification.setScoreType("q-value"); //higherIsBetter = false
            scoretype = true;
            break;
          }
        }
        else if (specific_score_terms.find(scoreit->first) != specific_score_terms.end())
        {
          score = scoreit->second.front().getValue().toString().toDouble(); // cast fix needed as DataValue is init with XercesString
          spectrum_identification.setHigherScoreBetter(ControlledVocabulary::CVTerm::isHigherBetterScore(cv_.getTerm(scoreit->first)));
          spectrum_identification.setScoreType(scoreit->second.front().getName());
          scoretype = true;
          break;
        }
        else if (e_score_terms.find(scoreit->first) != e_score_terms.end())
        {
          score = scoreit->second.front().getValue().toString().toDouble(); // cast fix needed as DataValue is init with XercesString
          spectrum_identification.setHigherScoreBetter(false);
          spectrum_identification.setScoreType("E-value"); //higherIsBetter = false
          scoretype = true;
          break;
        }
        else if (scoreit->first == "MS:1001143")
        {
          spectrum_identification.setScoreType("PSM-level search engine specific statistic");
          // TODO this is just an assumption for unknown scores
          spectrum_identification.setHigherScoreBetter(true);
          scoretype = true;
        }
      }
      if (scoretype) //else (i.e. no q/E/raw score or threshold not passed) no hit will be read TODO @mths: yielding no peptide hits will be error prone!!! what to do? remove and warn peptideidentifications with no hits inside?!
      {
        //build the PeptideHit from a SpectrumIdentificationItem
        PeptideHit hit(score, rank, chargeState, pep_map_[peptide_ref]);
        for (const auto& cvs : param_cv.getCVTerms())
        {
          for (const auto& cv : cvs.second) // if the same accession occurred multiple times...
          {
            if (cvs.first == "MS:1002540")
            {
              hit.setMetaValue(cvs.first, cv.getValue().toString());
            }
            else if (cvs.first == "MS:1001143") // this is the CV term "PSM-level search engine specific statistic" and it doesn't have a value
            {
              continue;
            }
            else
            { // value may be empty, e.g. <cvParam accession="MS:1001363" name="peptide unique to one protein" cvRef="PSI-MS" />
              auto value = xml_handler_->cvParamToValue(cv_, cv);
              hit.setMetaValue(cvs.first, value);
            }
          }
        }
        for (map<String, DataValue>::const_iterator up = param_user.begin(); up != param_user.end(); ++up)
        {
          hit.setMetaValue(up->first, up->second);
        }
        hit.setMetaValue("calcMZ", calculatedMassToCharge);
        spectrum_identification.setMZ(experimentalMassToCharge); // TODO @ mths for next PSI meeting: why is this not in SpectrumIdentificationResult in the schema? exp. m/z for one spec should not change from one id for it to the next!
        hit.setMetaValue("pass_threshold", pass); //TODO @ mths do not write metavalue pass_threshold

        //connect the PeptideHit with PeptideEvidences (for AABefore/After) and subsequently with DBSequence (for ProteinAccession)
        pair<multimap<String, String>::iterator, multimap<String, String>::iterator> pev_its;
        pev_its = p_pv_map_.equal_range(peptide_ref);
        for (multimap<String, String>::iterator pev_it = pev_its.first; pev_it != pev_its.second; ++pev_it)
        {
          bool idec = false;
          OpenMS::PeptideEvidence pev;
          if (pe_ev_map_.find(pev_it->second) != pe_ev_map_.end())
          {
            MzIdentMLDOMHandler::PeptideEvidence& pv = pe_ev_map_[pev_it->second];
            if (pv.pre != '-') pev.setAABefore(pv.pre);
            if (pv.post != '-') pev.setAAAfter(pv.post);

            if (pv.start != OpenMS::PeptideEvidence::UNKNOWN_POSITION && pv.stop != OpenMS::PeptideEvidence::UNKNOWN_POSITION)
            {
              hit.setMetaValue("start", pv.start);
              hit.setMetaValue("end", pv.stop);
              pev.setStart(pv.start);
              pev.setEnd(pv.stop);
            }

            idec = pv.idec;
            if (idec)
            {
              if (hit.metaValueExists(Constants::UserParam::TARGET_DECOY))
              {
                if (hit.getMetaValue(Constants::UserParam::TARGET_DECOY) != "decoy")
                {
                  hit.setMetaValue(Constants::UserParam::TARGET_DECOY, "target+decoy");
                }
                else
                {
                  hit.setMetaValue(Constants::UserParam::TARGET_DECOY, "decoy");
                }
              }
              else
              {
                hit.setMetaValue(Constants::UserParam::TARGET_DECOY, "decoy");
              }
            }
            else
            {
              if (hit.metaValueExists(Constants::UserParam::TARGET_DECOY))
              {
                if (hit.getMetaValue(Constants::UserParam::TARGET_DECOY) != "target")
                {
                  hit.setMetaValue(Constants::UserParam::TARGET_DECOY, "target+decoy");
                }
                else
                {
                  hit.setMetaValue(Constants::UserParam::TARGET_DECOY, "target");
                }
              }
              else
              {
                hit.setMetaValue(Constants::UserParam::TARGET_DECOY, "target");
              }
            }
          }

          if (pv_db_map_.find(pev_it->second) != pv_db_map_.end())
          {
            String& dpv = pv_db_map_[pev_it->second];
            DBSequence& db = db_sq_map_[dpv];
            pev.setProteinAccession(db.accession);

            if (pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).findHit(db.accession)
                == pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).getHits().end())
            { // butt ugly! TODO @ mths for ProteinInference
              pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).insertHit(ProteinHit());
              pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).getHits().back().setSequence(db.sequence);
              pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).getHits().back().setAccession(db.accession);
              if (idec)
              {
                pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).getHits().back().setMetaValue("isDecoy", "true");
              }
              else
              {
                pro_id_->at(si_pro_map_[spectrumIdentificationList_ref]).getHits().back().setMetaValue("isDecoy", "false");
              }
            }
          }
          hit.addPeptideEvidence(pev);
        }
        //sort necessary for tests, as DOM children are obviously randomly ordered.
//        vector<String> temp = spectrum_identification.getHits().back().getPeptideEvidences().back().getProteinAccessions();
//        sort(temp.begin(), temp.end());
//        spectrum_identification.getHits().back().setProteinAccessions(temp);
        spectrum_identification.insertHit(hit);

        // due to redundant references this is not needed! peptideref already maps to all those PeptideEvidence elements
        //      DOMElement* child = spectrumIdentificationItemElement->getFirstElementChild();
        //      while ( child )
        //      {
        //        if (XMLString::equals(child->getTagName(), CONST_XMLCH("PeptideEvidenceRef")))
        //        {
        //          ref = StringManager::convert(element_si->getAttribute(CONST_XMLCH("peptideEvidence_ref")));
        //          //...
        //          spectrum_identification.getHits().back().setAABefore(char acid);
        //          spectrum_identification.getHits().back().setAAAfter (char acid);
        //          break;
        //        }
        //        child = child->getNextElementSibling();
        //      }
      }
      // <Fragmentation> omitted for the time being
    }

    void MzIdentMLDOMHandler::parseProteinDetectionListElements_(DOMNodeList* proteinDetectionListElements)
    {
      const  XMLSize_t node_count = proteinDetectionListElements->getLength();
      for (XMLSize_t c = 0; c < node_count; ++c)
      {
        DOMNode* current_pr = proteinDetectionListElements->item(c);
        if (current_pr->getNodeType() && // true is not NULL
            current_pr->getNodeType() == DOMNode::ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          DOMElement* element_pr = dynamic_cast<xercesc::DOMElement*>(current_pr);

//          String id = StringManager::convert(element_pr->getAttribute(CONST_XMLCH("id")));
//          pair<CVTermList, map<String, DataValue> > params = parseParamGroup_(current_pr->getChildNodes());

          // TODO @mths : this needs to be a ProteinIdentification for the ProteinDetectionListElement which is not mandatory and used in downstream analysis ProteinInference etc.
//          pro_id_->push_back(ProteinIdentification());
//          pro_id_->back().setSearchEngine(search_engine_);
//          pro_id_->back().setSearchEngineVersion(search_engine_version_);
//          pro_id_->back().setIdentifier(search_engine_);

          //      SearchParameters  search_parameters_
          //      DateTime  date_
          //      String    protein_score_type_ <- from ProteinDetectionProtocol
          //      DoubleReal    protein_significance_threshold_ <- from ProteinDetectionProtocol

          DOMElement* child = element_pr->getFirstElementChild();
          while (child)
          {
            if (XMLString::equals(child->getTagName(), CONST_XMLCH("ProteinAmbiguityGroup")))
            {
              parseProteinAmbiguityGroupElement_(child, pro_id_->back());
            }
            child = child->getNextElementSibling();
          }
        }
      }
    }

    void MzIdentMLDOMHandler::parseProteinAmbiguityGroupElement_(DOMElement* proteinAmbiguityGroupElement, ProteinIdentification& protein_identification)
    {
//      String id = StringManager::convert(proteinAmbiguityGroupElement->getAttribute(CONST_XMLCH("id")));
//      pair<CVTermList, map<String, DataValue> > params = parseParamGroup_(proteinAmbiguityGroupElement->getChildNodes());

      //fill pro_id_->back() with content,
      DOMElement* child = proteinAmbiguityGroupElement->getFirstElementChild();
      while (child)
      {
        if (XMLString::equals(child->getTagName(), CONST_XMLCH("ProteinDetectionHypothesis")))
        {
          parseProteinDetectionHypothesisElement_(child, protein_identification);
        }
        child = child->getNextElementSibling();
      }
    }

    void MzIdentMLDOMHandler::parseProteinDetectionHypothesisElement_(DOMElement* proteinDetectionHypothesisElement, ProteinIdentification& protein_identification)
    {
      String dBSequence_ref = StringManager::convert(proteinDetectionHypothesisElement->getAttribute(CONST_XMLCH("dBSequence_ref")));

//      pair<CVTermList, map<String, DataValue> > params = parseParamGroup_(proteinDetectionHypothesisElement->getChildNodes());

      DBSequence& db = db_sq_map_[dBSequence_ref];

      protein_identification.insertHit(ProteinHit());
      protein_identification.getHits().back().setSequence(db.sequence);
      protein_identification.getHits().back().setAccession(db.accession);
//      protein_identification.getHits().back().setCoverage((long double)params.first.getCVTerms()["MS:1001093"].front().getValue()); //TODO @ mths: calc percent
//      protein_identification.getHits().back().setScore(boost::lexical_cast<double>(params.first.getCVTerms()["MS:1001171"].front().getValue())); //or any other score
    }

    AASequence MzIdentMLDOMHandler::parsePeptideSiblings_(DOMElement* peptide)
    {
      DOMNodeList* peptideSiblings = peptide->getChildNodes();
      const  XMLSize_t node_count = peptideSiblings->getLength();
      String as;
      //1. Sequence
      for (XMLSize_t c = 0; c < node_count; ++c)
      {
        DOMNode* current_sib = peptideSiblings->item(c);
        if (current_sib->getNodeType() && // true is not NULL
            current_sib->getNodeType() == DOMNode::ELEMENT_NODE)
        {
          DOMElement* element_sib = dynamic_cast<xercesc::DOMElement*>(current_sib);
          if (XMLString::equals(element_sib->getTagName(), CONST_XMLCH("PeptideSequence")))
          {
            DOMNode* tn = element_sib->getFirstChild();
            if (tn->getNodeType() == DOMNode::TEXT_NODE)
            {
              DOMText* data = dynamic_cast<DOMText*>(tn);
              const XMLCh* val = data->getWholeText();
              as = StringManager::convert(val);
            }
            else
            {
              throw std::runtime_error("ERROR : Non Text Node");
            }
          }
        }
      }
      //2. Substitutions
      for (XMLSize_t c = 0; c < node_count; ++c)
      {
        DOMNode* current_sib = peptideSiblings->item(c);
        if (current_sib->getNodeType() && // true is not NULL
            current_sib->getNodeType() == DOMNode::ELEMENT_NODE)
        {
          DOMElement* element_sib = dynamic_cast<xercesc::DOMElement*>(current_sib);
          if (XMLString::equals(element_sib->getTagName(), CONST_XMLCH("SubstitutionModification")))
          {

            String location = StringManager::convert(element_sib->getAttribute(CONST_XMLCH("location")));
            char originalResidue = StringManager::convert(element_sib->getAttribute(CONST_XMLCH("originalResidue")))[0];
            char replacementResidue = StringManager::convert(element_sib->getAttribute(CONST_XMLCH("replacementResidue")))[0];

            if (!location.empty())
            {
              as[location.toInt() - 1] = replacementResidue;
            }
            else if (as.hasSubstring(originalResidue)) //no location - every occurrence will be replaced
            {
              as.substitute(originalResidue, replacementResidue);
            }
            else
            {
              throw std::runtime_error("ERROR : Non Text Node");
            }
          }
        }
      }
      //3. Modifications
      as.trim();
      AASequence aas = AASequence::fromString(as);
      for (XMLSize_t c = 0; c < node_count; ++c)
      {
        DOMNode* current_sib = peptideSiblings->item(c);
        if (current_sib->getNodeType() && // true is not NULL
            current_sib->getNodeType() == DOMNode::ELEMENT_NODE)
        {
          DOMElement* element_sib = dynamic_cast<xercesc::DOMElement*>(current_sib);
          if (XMLString::equals(element_sib->getTagName(), CONST_XMLCH("Modification")))
          {
            SignedSize index = -2;
            try
            {
              index = static_cast<SignedSize>(StringManager::convert(element_sib->getAttribute(CONST_XMLCH("location"))).toInt());
            }
            catch (...)
            {
              OPENMS_LOG_WARN << "Found unreadable modification location." << endl;
            }

            //double monoisotopicMassDelta = StringManager::convert(element_dbs->getAttribute(CONST_XMLCH("monoisotopicMassDelta")));

            if (xl_ms_search_) // special case: XL-MS search results
            {
              String pep_id = StringManager::convert(peptide->getAttribute(CONST_XMLCH("id")));
              //DOMNodeList* cvParams = element_sib->getElementsByTagName(CONST_XMLCH("cvParam"));
              DOMElement* cvp = element_sib->getFirstElementChild();
              //for (XMLSize_t i = 0; i < cvParams.length(); ++i)
              bool donor_acceptor_found = false;
              bool xlink_mod_found = false;

              while (cvp)
              {
                if (XMLString::equals(cvp->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1002509"))) // crosslink donor
                {
                  String donor_val = StringManager::convert(cvp->getAttribute(CONST_XMLCH("value")));
                  xl_id_donor_map_.insert(make_pair(pep_id, donor_val));
                  String massdelta = StringManager::convert(element_sib->getAttribute(CONST_XMLCH("monoisotopicMassDelta")));
                  double monoisotopicMassDelta = massdelta.toDouble();
                  xl_mass_map_.insert(make_pair(pep_id, monoisotopicMassDelta));
                  xl_donor_pos_map_.insert(make_pair(donor_val, index-1));

                  DOMElement* cvp1 = element_sib->getFirstElementChild();
                  String xl_mod_name = StringManager::convert(cvp1->getAttribute(CONST_XMLCH("name")));
                  xl_mod_map_.insert(make_pair(pep_id, xl_mod_name));
                  donor_acceptor_found = true;
                }
                else if (XMLString::equals(cvp->getAttribute(CONST_XMLCH("accession")), CONST_XMLCH("MS:1002510"))) // crosslink acceptor
                {
                  String acceptor_val = StringManager::convert(cvp->getAttribute(CONST_XMLCH("value")));
                  xl_id_acceptor_map_.insert(make_pair(pep_id, acceptor_val));
                  xl_acceptor_pos_map_.insert(make_pair(acceptor_val, index-1));
                  donor_acceptor_found = true;
                }
                else
                {
                  CVTerm cv = parseCvParam_(cvp);
                  const String& cvname = cv.getName();
                  if (cvname.hasPrefix("Xlink") || cv.getAccession().hasPrefix("XLMOD"))
                  {
                    xlink_mod_found = true;
                  }
                  if (cvname.hasSubstring("unknown mono-link"))
                  {
                    xl_mod_map_.insert(make_pair(pep_id, cvname));
                  }
                  else // normal mod, copied from below
                  {
                    if ( (cv.getCVIdentifierRef() != "UNIMOD") && (cv.getCVIdentifierRef() != "XLMOD") )
                    {
                      // e.g.  <cvParam accession="MS:1001524" name="fragment neutral loss" cvRef="PSI-MS" value="0" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO"/>
                      cvp = cvp->getNextElementSibling();
                      continue;
                    }
                    if (index == 0)
                    {
                      try // does not work for cross-links yet, but the information is finally stored as MetaValues of the PeptideHit
                      {
                        if (cvname == "unknown modification")
                        {
                          const String & cvvalue = cv.getValue();
                          if (ModificationsDB::getInstance()->has(cvvalue) && !cvvalue.empty())
                          {
                            aas.setNTerminalModification(cv.getValue());
                          }
                        }
                        else
                        {
                          aas.setNTerminalModification(cvname);
                        }
                        cvp = cvp->getNextElementSibling();
                        continue;
                      }
                      catch (...)
                      {
                        // TODO Residue and AASequence should use CrossLinksDB as well
                      }
                    }
                    else if (index == static_cast<SignedSize>(aas.size() + 1))
                    {
                      // TODO: does not work for cross-links yet, but the information is finally stored as MetaValues of the PeptideHit
                      if (cvname == "unknown modification")
                      {
                        const String & cvvalue = cv.getValue();
                        if (ModificationsDB::getInstance()->has(cvvalue) && !cvvalue.empty())
                        {
                          aas.setCTerminalModification(cvvalue);
                        }
                        else
                        {
                          OPENMS_LOG_WARN << "Modification: " << cvvalue << " not found in ModificationsDB." << endl;
                          // TODO? @enetz look it up in XLDB?
                        }
                      }
                      else
                      {
                        if (ModificationsDB::getInstance()->has(cvname))
                        {
                          aas.setCTerminalModification(cvname);
                        }
                        else
                        {
                          OPENMS_LOG_WARN << "Modification: " << cvname << " not found in ModificationsDB." << endl;
                          // TODO? @enetz look it up in XLDB?
                        }
                      }
                      cvp = cvp->getNextElementSibling();
                      continue;
                    }
                    else
                    {
                      if (cvname == "unknown modification")
                      {
                        const String & cvvalue = cv.getValue();
                        if (ModificationsDB::getInstance()->has(cvvalue) && !cvvalue.empty())
                        {
                          aas.setModification(index - 1, cvvalue);
                        }
                        else
                        {
                          OPENMS_LOG_WARN << "Modification: " << cvvalue << " not found in ModificationsDB." << endl;
                          // TODO? @enetz look it up in XLDB?
                        }
                      }
                      else
                      {
                        if (ModificationsDB::getInstance()->has(cvname))
                        {
                          aas.setModification(index - 1, cvname);
                        }
                        else
                        {
                          OPENMS_LOG_WARN << "Modification: " << cvname << " not found in ModificationsDB." << endl;
                          // TODO? @enetz look it up in XLDB?
                        }
                      }
                      cvp = cvp->getNextElementSibling();
                      continue;
/* TODO enetz: look up in XLDB and remove this ha
                        // this is a bad hack to avoid a long list of warnings in the case of XL-MS data
                        if ( !(String(e.what()).hasSubstring("'DSG'") || String(e.what()).hasSubstring("'DSS'") || String(e.what()).hasSubstring("'EDC'")) || String(e.what()).hasSubstring("'BS3'") || String(e.what()).hasSubstring("'BS2G'") )
                        {
                          OPENMS_LOG_WARN << e.getName() << ": " << e.what() << " Sequence: " << aas.toUnmodifiedString() << ", residue " << aas.getResidue(index - 1).getName() << "@" << String(index) << endl;
                        }
*/
                    }
                  }
                }
                cvp = cvp->getNextElementSibling();
              }
              if ( (!donor_acceptor_found) && (xlink_mod_found) ) // mono-link, here using pep_id also as the CV value, since mono-links don't have a cross-linking CV term
              {
                xl_id_donor_map_.insert(make_pair(pep_id, pep_id));
                xl_donor_pos_map_.insert(make_pair(pep_id, index-1));
              }
            }
            else // general case: no XL-MS result
            {
              DOMElement* cvp = element_sib->getFirstElementChild();
              while (cvp)
              {
                CVTerm cv = parseCvParam_(cvp);
                if (cv.getAccession() == "MS:1001460") // unknown modification
                {
                  const String & cvvalue = cv.getValue();
                  if (cv.hasValue() && ModificationsDB::getInstance()->has(cvvalue) && !cvvalue.empty())  // why do we need to check for empty?
                  {
                    // Case 1: unknown (to e.g., third-party tool) modification known to OpenMS (see value)
                    //  <Modification location="0" monoisotopicMassDelta="17.031558">
                    //  <cvParam cvRef="PSI-MS" accession="MS:1001460" name="unknown modification" value="Methyl:2H(2)13C"/>
                    const String & mname = cvvalue;
                    if (index == 0)
                    {
                      aas.setNTerminalModification(mname);
                    }
                    else if (index == (int)aas.size() + 1)
                    {
                      aas.setCTerminalModification(mname);
                    }
                    else if (index > 0 && index <= (int)aas.size() )
                    {
                      aas.setModification(index - 1, mname);
                    }
                    cvp = cvp->getNextElementSibling();
                    continue;
                  }
                  else
                  {
                    // Case 2: unknown modification (needs to be added to ModificationsDB)
                    // note, this is optional
                    double mass_delta = 0;
                    bool has_mass_delta = false;
                    String mod;

                    // try to parse information, give up if we cannot
                    try
                    {
                      mod = StringManager::convert(element_sib->getAttribute(CONST_XMLCH("monoisotopicMassDelta")));
                      mass_delta = static_cast<double>(mod.toDouble());
                      has_mass_delta = true;
                    }
                    catch (...)
                    {
                      OPENMS_LOG_WARN << "Found unreadable modification location." << endl;
                      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown modification");
                    }
                    if (!has_mass_delta)
                    {
                      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "Unknown modification");
                    }

                    // Parse this and add a new modification of mass "monoisotopicMassDelta" to the AASequence
                    // e.g. <cvParam cvRef="MS" accession="MS:1001460" name="unknown modification" value="N-Glycan"/>

                    // compare with String::ConstIterator AASequence::parseModSquareBrackets_
                    ModificationsDB* mod_db = ModificationsDB::getInstance();
                    if (index == 0)
                    {
                      // n-terminal
                      String residue_name = ".[+" + mod + "]";
                      String residue_id = ".n[" + mod + "]";

                      // Check if it already exists, if not create new modification, transfer
                      // ownership to ModDB
                      if (!mod_db->has(residue_id))
                      {
                        unique_ptr<ResidueModification> new_mod(new ResidueModification);
                        new_mod->setFullId(residue_id); // setting FullId but not Id makes it a user-defined mod
                        new_mod->setFullName(residue_name); // display name
                        new_mod->setDiffMonoMass(mass_delta);
                        new_mod->setMonoMass(mass_delta + Residue::getInternalToNTerm().getMonoWeight());
                        new_mod->setTermSpecificity(ResidueModification::N_TERM);
                        mod_db->addModification(std::move(new_mod));
                      }

                      aas.setNTerminalModification(residue_id);
                      cvp = cvp->getNextElementSibling();
                      continue;
                    }
                    else if (index == (int)aas.size() +1)
                    {
                      // c-terminal
                      String residue_name = ".[" + mod + "]";
                      String residue_id = ".c[" + mod + "]";

                      // Check if it already exists, if not create new modification, transfer
                      // ownership to ModDB
                      if (!mod_db->has(residue_name))
                      {
                        unique_ptr<ResidueModification> new_mod(new ResidueModification);
                        new_mod->setFullId(residue_id); // setting FullId but not Id makes it a user-defined mod
                        new_mod->setFullName(residue_name); // display name
                        new_mod->setDiffMonoMass(mass_delta);
                        new_mod->setMonoMass(mass_delta + Residue::getInternalToCTerm().getMonoWeight());
                        new_mod->setTermSpecificity(ResidueModification::C_TERM);
                        mod_db->addModification(std::move(new_mod));
                      }
                      aas.setCTerminalModification(residue_id);
                      cvp = cvp->getNextElementSibling();
                      continue;
                    }
                    else if (index > 0 && index <= (int)aas.size() )
                    {
                      // internal modification
                      const Residue& residue = aas[index-1];
                      String residue_name = residue.getOneLetterCode() + "[" + mod + "]"; // e.g. N[12345.6]
                      String modification_name = "[" + mod + "]";

                      if (!mod_db->has(residue_name))
                      {
                        // create new modification
                        unique_ptr<ResidueModification> new_mod(new ResidueModification);
                        new_mod->setFullId(residue_name); // setting FullId but not Id makes it a user-defined mod
                        new_mod->setFullName(modification_name); // display name

                        // We will set origin to make sure the same modification will be used
                        // for the same AA
                        new_mod->setOrigin(residue.getOneLetterCode()[0]);

                        // We cannot set origin if we want to use the same modification name
                        // also at other AA (and since we have no information here, it is safer
                        // to assume that this may happen).
                        // new_mod->setOrigin(residue.getOneLetterCode()[0]);

                        new_mod->setMonoMass(mass_delta + residue.getMonoWeight());
                        new_mod->setAverageMass(mass_delta + residue.getAverageWeight());
                        new_mod->setDiffMonoMass(mass_delta);

                        mod_db->addModification(std::move(new_mod));
                      }

                      // now use the new modification
                      Size mod_idx = mod_db->findModificationIndex(residue_name);
                      const ResidueModification* res_mod = mod_db->getModification(mod_idx);

                      // Set a modification on the given AA
                      // Note: this calls setModification_ on a new Residue which changes its
                      // weight to the weight of the modification (set above)
                      //
                      aas.setModification(index-1, res_mod->getFullId());
                      cvp = cvp->getNextElementSibling();
                      continue;
                    }
                  }
                }
                if (cv.getCVIdentifierRef() != "UNIMOD")
                {
                  // e.g.  <cvParam accession="MS:1001524" name="fragment neutral loss" cvRef="PSI-MS" value="0" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO"/>
                  cvp = cvp->getNextElementSibling();
                  continue;
                }

                if (index == 0)
                {
                  if (cv.getName() == "unknown modification")
                  {
                    aas.setNTerminalModification(cv.getValue());
                    cvp = cvp->getNextElementSibling();
                    continue;
                  }
                  else
                  {
                    aas.setNTerminalModification(cv.getName());
                    cvp = cvp->getNextElementSibling();
                    continue;
                  }
                }
                else if (index == static_cast<SignedSize>(aas.size() + 1))
                {
                  aas.setCTerminalModification(cv.getName());
                  cvp = cvp->getNextElementSibling();
                  continue;
                }
                else
                {
                  try
                  {
                     aas.setModification(index - 1, cv.getName()); //TODO @mths,Timo : do this via UNIMOD accessions
                     cvp = cvp->getNextElementSibling();
                     continue;
                  }
                  catch (Exception::BaseException& e)
                  {
                    OPENMS_LOG_WARN << e.getName() << ": " << e.what() << " Sequence: " << aas.toUnmodifiedString() << ", residue " << aas.getResidue(index - 1).getName() << "@" << String(index) << "\n";
                  }
                }

                cvp = cvp->getNextElementSibling();
              }
            }
          }
        }
      }
      return aas;
    }

    void MzIdentMLDOMHandler::buildCvList_(DOMElement* cvElements)
    {
      DOMElement* cv1 = cvElements->getOwnerDocument()->createElement(CONST_XMLCH("cv"));
      cv1->setAttribute(CONST_XMLCH("id"), CONST_XMLCH("PSI-MS"));
      cv1->setAttribute(CONST_XMLCH("fullName"),
                        CONST_XMLCH("Proteomics Standards Initiative Mass Spectrometry Vocabularies"));
      cv1->setAttribute(CONST_XMLCH("uri"),
                        CONST_XMLCH("http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo"));
      cv1->setAttribute(CONST_XMLCH("version"), CONST_XMLCH("2.32.0"));
      cvElements->appendChild(cv1);
      DOMElement* cv2 = cvElements->getOwnerDocument()->createElement(CONST_XMLCH("cv"));
      cv2->setAttribute(CONST_XMLCH("id"), CONST_XMLCH("UNIMOD"));
      cv2->setAttribute(CONST_XMLCH("fullName"),
                        CONST_XMLCH("UNIMOD"));
      cv2->setAttribute(CONST_XMLCH("uri"),
                        CONST_XMLCH("http://www.unimod.org/obo/unimod.obo"));
      cvElements->appendChild(cv2);
      DOMElement* cv3 = cvElements->getOwnerDocument()->createElement(CONST_XMLCH("cv"));
      cv3->setAttribute(CONST_XMLCH("id"), CONST_XMLCH("UO"));
      cv3->setAttribute(CONST_XMLCH("fullName"),
                        CONST_XMLCH("UNIT-ONTOLOGY"));
      cv3->setAttribute(CONST_XMLCH("uri"),
                        CONST_XMLCH("http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo"));
      cvElements->appendChild(cv3);
    }

    void MzIdentMLDOMHandler::buildAnalysisSoftwareList_(DOMElement* analysisSoftwareElements)
    {
      DOMElement* current_as = analysisSoftwareElements->getOwnerDocument()->createElement(CONST_XMLCH("AnalysisSoftware"));
      current_as->setAttribute(CONST_XMLCH("id"), StringManager::convertPtr((String("OpenMS") + UniqueIdGenerator::getUniqueId())).get());
      current_as->setAttribute(CONST_XMLCH("version"), CONST_XMLCH("search_engine_version_"));
      current_as->setAttribute(CONST_XMLCH("name"), CONST_XMLCH("search_engine_"));
      analysisSoftwareElements->appendChild(current_as);
      DOMElement* current_sw = current_as->getOwnerDocument()->createElement(CONST_XMLCH("SoftwareName"));

      //TODO build extract as function and insert cv
      DOMElement* current_cv = current_sw->getOwnerDocument()->createElement(CONST_XMLCH("cvParam"));
      current_cv->setAttribute(CONST_XMLCH("name"), CONST_XMLCH("search_engine_"));
      current_cv->setAttribute(CONST_XMLCH("cvRef"), CONST_XMLCH("PSI-MS"));

      //TODO this needs error handling
      current_cv->setAttribute(CONST_XMLCH("accession"), StringManager::convertPtr(cv_.getTermByName("search_engine_").id).get());
      current_sw->appendChild(current_cv);
      analysisSoftwareElements->appendChild(current_sw);
    }

    void MzIdentMLDOMHandler::buildSequenceCollection_(DOMElement* sequenceCollectionElements)
    {
      for (map<String, DBSequence>::iterator dbs = db_sq_map_.begin(); dbs != db_sq_map_.end(); ++dbs)
      {
        DOMElement* current_dbs = sequenceCollectionElements->getOwnerDocument()->createElement(CONST_XMLCH("DBSequence"));
        current_dbs->setAttribute(CONST_XMLCH("id"), StringManager::convertPtr(dbs->second.accession).get());
        current_dbs->setAttribute(CONST_XMLCH("length"), StringManager::convertPtr(String(dbs->second.sequence.length())).get());
        current_dbs->setAttribute(CONST_XMLCH("accession"), StringManager::convertPtr(dbs->second.accession).get());
        current_dbs->setAttribute(CONST_XMLCH("searchDatabase_ref"), StringManager::convertPtr(dbs->second.database_ref).get()); // This is going to be wrong
        DOMElement* current_seq = current_dbs->getOwnerDocument()->createElement(CONST_XMLCH("Seq"));
        DOMText* current_seqnot = current_seq->getOwnerDocument()->createTextNode(StringManager::convertPtr(dbs->second.sequence).get());
        current_seq->appendChild(current_seqnot);
        current_dbs->appendChild(current_seq);
        sequenceCollectionElements->appendChild(current_dbs);
      }

      for (map<String, AASequence>::iterator peps = pep_map_.begin(); peps != pep_map_.end(); ++peps)
      {
        DOMElement* current_pep = sequenceCollectionElements->getOwnerDocument()->createElement(CONST_XMLCH("Peptide"));
        current_pep->setAttribute(CONST_XMLCH("id"), StringManager::convertPtr(peps->first).get());
        DOMElement* current_seq = current_pep->getOwnerDocument()->createElement(CONST_XMLCH("PeptideSequence"));
        DOMText* current_seqnot = current_seq->getOwnerDocument()->createTextNode(StringManager::convertPtr(peps->second.toUnmodifiedString()).get());
        current_seq->appendChild(current_seqnot);
        current_pep->appendChild(current_seq);
        if (peps->second.hasNTerminalModification())
        {
          const ResidueModification* mod = peps->second.getNTerminalModification();
          DOMElement* current_mod = current_pep->getOwnerDocument()->createElement(CONST_XMLCH("Modification"));
          DOMElement* current_cv = current_pep->getOwnerDocument()->createElement(CONST_XMLCH("cvParam"));
          current_mod->setAttribute(CONST_XMLCH("location"), CONST_XMLCH("0"));
          current_mod->setAttribute(CONST_XMLCH("monoisotopicMassDelta"), StringManager::convertPtr(String(mod->getDiffMonoMass())).get());
          String origin = mod->getOrigin();
          if (origin == "X") origin = ".";
          current_mod->setAttribute(CONST_XMLCH("residues"), StringManager::convertPtr(origin).get());
          current_cv->setAttribute(CONST_XMLCH("name"), StringManager::convertPtr(mod->getName()).get());
          current_cv->setAttribute(CONST_XMLCH("cvRef"), CONST_XMLCH("UNIMOD"));
          current_cv->setAttribute(CONST_XMLCH("accession"), StringManager::convertPtr(mod->getUniModAccession()).get());

          current_mod->appendChild(current_cv);
          current_pep->appendChild(current_mod);
        }
        if (peps->second.hasCTerminalModification())
        {
          const ResidueModification* mod = peps->second.getCTerminalModification();
          DOMElement* current_mod = current_pep->getOwnerDocument()->createElement(CONST_XMLCH("Modification"));
          DOMElement* current_cv = current_mod->getOwnerDocument()->createElement(CONST_XMLCH("cvParam"));
          current_mod->setAttribute(CONST_XMLCH("location"), StringManager::convertPtr(String(peps->second.size() + 1)).get());
          current_mod->setAttribute(CONST_XMLCH("monoisotopicMassDelta"), StringManager::convertPtr(String(mod->getDiffMonoMass())).get());
          String origin = mod->getOrigin();
          if (origin == "X") origin = ".";
          current_mod->setAttribute(CONST_XMLCH("residues"), StringManager::convertPtr(origin).get());
          current_cv->setAttribute(CONST_XMLCH("name"), StringManager::convertPtr(mod->getName()).get());
          current_cv->setAttribute(CONST_XMLCH("cvRef"), CONST_XMLCH("UNIMOD"));
          current_cv->setAttribute(CONST_XMLCH("accession"), StringManager::convertPtr(mod->getUniModAccession()).get());

          current_mod->appendChild(current_cv);
          current_pep->appendChild(current_mod);
        }
        if (peps->second.isModified())
        {
          Size i = 0;
          for (AASequence::ConstIterator res = peps->second.begin(); res != peps->second.end(); ++res, ++i)
          {
            const ResidueModification* mod = res->getModification();
            if (mod == nullptr) continue;
            DOMElement* current_mod = current_pep->getOwnerDocument()->createElement(CONST_XMLCH("Modification"));
            DOMElement* current_cv = current_pep->getOwnerDocument()->createElement(CONST_XMLCH("cvParam"));
            current_mod->setAttribute(CONST_XMLCH("location"), StringManager::convertPtr(String(i)).get());
            current_mod->setAttribute(CONST_XMLCH("monoisotopicMassDelta"), StringManager::convertPtr(String(mod->getDiffMonoMass())).get());
            current_mod->setAttribute(CONST_XMLCH("residues"), StringManager::convertPtr(String(mod->getOrigin())).get());
            current_cv->setAttribute(CONST_XMLCH("name"), StringManager::convertPtr(mod->getName()).get());
            current_cv->setAttribute(CONST_XMLCH("cvRef"), CONST_XMLCH("UNIMOD"));
            current_cv->setAttribute(CONST_XMLCH("accession"), StringManager::convertPtr(mod->getUniModAccession()).get());

            current_mod->appendChild(current_cv);
            current_pep->appendChild(current_mod);
          }
        }
        sequenceCollectionElements->appendChild(current_pep);
      }

      for (map<String, PeptideEvidence>::iterator pevs = pe_ev_map_.begin(); pevs != pe_ev_map_.end(); ++pevs)
      {
        DOMElement* current_pev = sequenceCollectionElements->getOwnerDocument()->createElement(CONST_XMLCH("PeptideEvidence"));
        current_pev->setAttribute(CONST_XMLCH("peptide_ref"), CONST_XMLCH("TBA"));
        current_pev->setAttribute(CONST_XMLCH("id"), StringManager::convertPtr(pevs->first).get());
        current_pev->setAttribute(CONST_XMLCH("start"), StringManager::convertPtr(String(pevs->second.start)).get());
        current_pev->setAttribute(CONST_XMLCH("end"), StringManager::convertPtr(String(pevs->second.stop)).get());
        current_pev->setAttribute(CONST_XMLCH("pre"), StringManager::convertPtr(String(pevs->second.pre)).get());
        current_pev->setAttribute(CONST_XMLCH("post"), StringManager::convertPtr(String(pevs->second.post)).get());
        // do not forget to annotate the decoy
        current_pev->setAttribute(CONST_XMLCH("isDecoy"), CONST_XMLCH("false"));
        sequenceCollectionElements->appendChild(current_pev);
      }
    }

    void MzIdentMLDOMHandler::buildAnalysisCollection_(DOMElement* analysisCollectionElements)
    {
      // for now there is only one search per file
      DOMElement* current_si = analysisCollectionElements->getOwnerDocument()->createElement(CONST_XMLCH("SpectrumIdentification"));
      current_si->setAttribute(CONST_XMLCH("id"), CONST_XMLCH("TBA"));
      current_si->setAttribute(CONST_XMLCH("spectrumIdentificationProtocol_ref"), CONST_XMLCH("SIP"));
      current_si->setAttribute(CONST_XMLCH("spectrumIdentificationList_ref"), CONST_XMLCH("SIL"));
      current_si->setAttribute(CONST_XMLCH("activityDate"), CONST_XMLCH("now"));
      DOMElement* current_is = current_si->getOwnerDocument()->createElement(CONST_XMLCH("InputSpectra"));
      current_is->setAttribute(CONST_XMLCH("spectraData_ref"), CONST_XMLCH("TODO")); // TODO @ mths while DataCollection
      DOMElement* current_sr = current_si->getOwnerDocument()->createElement(CONST_XMLCH("SearchDatabaseRef"));
      current_sr->setAttribute(CONST_XMLCH("searchDatabase_ref"), CONST_XMLCH("TODO")); // TODO @ mths while DataCollection
      current_si->appendChild(current_is);
      current_si->appendChild(current_sr);
      // and no ProteinDetection for now
      analysisCollectionElements->appendChild(current_si);
    }

    void MzIdentMLDOMHandler::buildAnalysisProtocolCollection_(DOMElement* protocolElements)
    {
      // for now there is only one search per file
      DOMElement* current_sp = protocolElements->getOwnerDocument()->createElement(CONST_XMLCH("SpectrumIdentificationProtocol"));
      current_sp->setAttribute(CONST_XMLCH("id"), CONST_XMLCH("SIP"));
      current_sp->setAttribute(CONST_XMLCH("analysisSoftware_ref"), CONST_XMLCH("what now?"));
      protocolElements->appendChild(current_sp);
      DOMElement* current_st = current_sp->getOwnerDocument()->createElement(CONST_XMLCH("SearchType"));
      current_sp->appendChild(current_st);
      DOMElement* current_cv = current_st->getOwnerDocument()->createElement(CONST_XMLCH("cvParam"));
      current_cv->setAttribute(CONST_XMLCH("accession"), CONST_XMLCH("MS:1001083")); // TODO @ mths for now static cv
      current_cv->setAttribute(CONST_XMLCH("name"), CONST_XMLCH("ms-ms search"));
      current_cv->setAttribute(CONST_XMLCH("cvRef"), CONST_XMLCH("PSI-MS"));
      current_st->appendChild(current_cv);

      //for now no <AdditionalSearchParams>, <ModificationParams>, <MassTable id="MT" msLevel="1 2">, <FragmentTolerance>, <ParentTolerance>, <DatabaseFilters>, <DatabaseTranslations>

      DOMElement* current_th = current_sp->getOwnerDocument()->createElement(CONST_XMLCH("Threshold"));
      DOMElement* current_up = current_th->getOwnerDocument()->createElement(CONST_XMLCH("userParam"));
      current_up->setAttribute(CONST_XMLCH("value"), CONST_XMLCH("0.05")); // TODO @ mths for now static cv
      current_up->setAttribute(CONST_XMLCH("name"), CONST_XMLCH("some significance threshold"));
      current_st->appendChild(current_up);
      // and no ProteinDetection for now
      protocolElements->appendChild(current_th);
    }

    void MzIdentMLDOMHandler::buildInputDataCollection_(DOMElement* inputElements)
    {
      DOMElement* current_sf = inputElements->getOwnerDocument()->createElement(CONST_XMLCH("SourceFile"));
      current_sf->setAttribute(CONST_XMLCH("location"), CONST_XMLCH("file:///tmp/test.dat"));
      current_sf->setAttribute(CONST_XMLCH("id"), CONST_XMLCH("SF1"));
      buildEnclosedCV_(current_sf, "FileFormat", "MS:1001199", "Mascot DAT file", "PSI-MS"); // TODO @ mths for now static cv
      inputElements->appendChild(current_sf);

      DOMElement* current_sd = inputElements->getOwnerDocument()->createElement(CONST_XMLCH("SearchDatabase"));
      current_sd->setAttribute(CONST_XMLCH("location"), CONST_XMLCH("file:///tmp/test.fasta"));
      current_sd->setAttribute(CONST_XMLCH("id"), CONST_XMLCH("DB1"));
      current_sd->setAttribute(CONST_XMLCH("name"), CONST_XMLCH("SwissProt"));
      current_sd->setAttribute(CONST_XMLCH("numDatabaseSequences"), CONST_XMLCH("257964"));
      current_sd->setAttribute(CONST_XMLCH("numResidues"), CONST_XMLCH("93947433"));
      current_sd->setAttribute(CONST_XMLCH("releaseDate"), CONST_XMLCH("2011-03-01T21:32:52"));
      current_sd->setAttribute(CONST_XMLCH("version"), CONST_XMLCH("SwissProt_51.6.fasta"));
      buildEnclosedCV_(current_sd, "FileFormat", "MS:1001348", "FASTA format", "PSI-MS"); // TODO @ mths for now static cv
      DOMElement* current_dn = current_sd->getOwnerDocument()->createElement(CONST_XMLCH("DatabaseName"));
      DOMElement* current_up = current_dn->getOwnerDocument()->createElement(CONST_XMLCH("userParam"));
      current_up->setAttribute(CONST_XMLCH("name"), CONST_XMLCH("SwissProt_51.6.fasta")); // TODO @ mths for now static cv
      current_dn->appendChild(current_up);
      current_sd->appendChild(current_dn);

      DOMElement* current_cv = current_sd->getOwnerDocument()->createElement(CONST_XMLCH("cvParam"));
      current_cv->setAttribute(CONST_XMLCH("accession"), CONST_XMLCH("MS:1001073")); // TODO @ mths for now static cv
      current_cv->setAttribute(CONST_XMLCH("name"), CONST_XMLCH("database type amino acid"));
      current_cv->setAttribute(CONST_XMLCH("cvRef"), CONST_XMLCH("PSI-MS"));
      current_sd->appendChild(current_cv);
      inputElements->appendChild(current_sd);

      DOMElement* current_spd = inputElements->getOwnerDocument()->createElement(CONST_XMLCH("SpectraData"));
      current_spd->setAttribute(CONST_XMLCH("location"), CONST_XMLCH("file:///tmp/test.mzML"));
      current_spd->setAttribute(CONST_XMLCH("id"), CONST_XMLCH("SD1"));
      buildEnclosedCV_(current_spd, "FileFormat", "MS:1001062", "Mascot MGF file", "PSI-MS");
      buildEnclosedCV_(current_spd, "SpectrumIDFormat", "MS:1001528", "Mascot query number", "PSI-MS");
      inputElements->appendChild(current_spd);
    }

    void MzIdentMLDOMHandler::buildEnclosedCV_(DOMElement* parentElement, const String& encel, const String& acc, const String& name, const String& cvref)
    {
      DOMElement* current_ff = parentElement->getOwnerDocument()->createElement(StringManager::convertPtr(encel).get());
      DOMElement* current_cv = current_ff->getOwnerDocument()->createElement(CONST_XMLCH("cvParam"));
      current_cv->setAttribute(CONST_XMLCH("accession"), StringManager::convertPtr(acc).get());
      current_cv->setAttribute(CONST_XMLCH("name"), StringManager::convertPtr(name).get());
      current_cv->setAttribute(CONST_XMLCH("cvRef"), StringManager::convertPtr(cvref).get());
      current_ff->appendChild(current_cv);
      parentElement->appendChild(current_ff);
    }

    void MzIdentMLDOMHandler::buildAnalysisDataCollection_(DOMElement* analysisElements)
    {
      DOMElement* current_sil = analysisElements->getOwnerDocument()->createElement(CONST_XMLCH("SpectrumIdentificationList"));
      current_sil->setAttribute(CONST_XMLCH("id"), CONST_XMLCH("SIL1"));
      current_sil->setAttribute(CONST_XMLCH("numSequencesSearched"), CONST_XMLCH("TBA"));
      // for now no FragmentationTable

      for (vector<PeptideIdentification>::iterator pi = pep_id_->begin(); pi != pep_id_->end(); ++pi)
      {
        DOMElement* current_sr = current_sil->getOwnerDocument()->createElement(CONST_XMLCH("SpectrumIdentificationResult"));
        current_sr->setAttribute(CONST_XMLCH("id"), StringManager::convertPtr(String(UniqueIdGenerator::getUniqueId())).get());
        current_sr->setAttribute(CONST_XMLCH("spectrumID"), StringManager::convertPtr(String(UniqueIdGenerator::getUniqueId())).get());
        current_sr->setAttribute(CONST_XMLCH("spectraData_ref"), CONST_XMLCH("SD1"));
        for (vector<PeptideHit>::iterator ph = pi->getHits().begin(); ph != pi->getHits().end(); ++ph)
        {
          DOMElement* current_si = current_sr->getOwnerDocument()->createElement(CONST_XMLCH("SpectrumIdentificationItem"));
          current_si->setAttribute(CONST_XMLCH("id"), StringManager::convertPtr(String(UniqueIdGenerator::getUniqueId())).get());
          current_si->setAttribute(CONST_XMLCH("calculatedMassToCharge"), StringManager::convertPtr(String(ph->getSequence().getMonoWeight(Residue::Full, ph->getCharge()))).get()); //TODO @mths : this is not correct!1elf - these interfaces are BS!
          current_si->setAttribute(CONST_XMLCH("chargeState"), StringManager::convertPtr(String(ph->getCharge())).get());
          current_si->setAttribute(CONST_XMLCH("experimentalMassToCharge"), StringManager::convertPtr(String(ph->getSequence().getMonoWeight(Residue::Full, ph->getCharge()))).get()); //TODO @mths : this is not correct!1elf - these interfaces are BS!
          current_si->setAttribute(CONST_XMLCH("peptide_ref"), CONST_XMLCH("TBA"));
          current_si->setAttribute(CONST_XMLCH("rank"), StringManager::convertPtr(String(ph->getRank())).get());
          current_si->setAttribute(CONST_XMLCH("passThreshold"), CONST_XMLCH("TBA"));
          current_si->setAttribute(CONST_XMLCH("sample_ref"), CONST_XMLCH("TBA"));
          // TODO cvs for score!
          current_sr->appendChild(current_si);
          for (list<String>::iterator pepevref = hit_pev_.front().begin(); pepevref != hit_pev_.front().end(); ++pepevref)
          {
            DOMElement* current_per = current_si->getOwnerDocument()->createElement(CONST_XMLCH("PeptideEvidenceRef"));
            current_per->setAttribute(CONST_XMLCH("peptideEvidence_ref"), StringManager::convertPtr(*pepevref).get());
            current_si->appendChild(current_per);
          }
          hit_pev_.pop_front();
          // and no Fragmentation annotation for now
        }
//          <cvParam accession="MS:1001371" name="Mascot:identity threshold" cvRef="PSI-MS" value="44"/>
//          <cvParam accession="MS:1001370" name="Mascot:homology threshold" cvRef="PSI-MS" value="18"/>
//          <cvParam accession="MS:1001030" name="number of peptide seqs compared to each spectrum" cvRef="PSI-MS" value="26981"/>
//          <cvParam accession="MS:1000796" name="spectrum title" cvRef="PSI-MS" value="dp210198 21-Jan-98 DERIVED SPECTRUM    #9"/>
        current_sil->appendChild(current_sr);
      }

      // and no ProteinDetection for now
    }

    ProteinIdentification::SearchParameters MzIdentMLDOMHandler::findSearchParameters_(pair<CVTermList, map<String, DataValue> > as_params)
    {
      ProteinIdentification::SearchParameters sp = ProteinIdentification::SearchParameters();
      for (std::map<String, vector<CVTerm> >::const_iterator cvs = as_params.first.getCVTerms().begin(); cvs != as_params.first.getCVTerms().end(); ++cvs)
      {
        for (vector<CVTerm>::const_iterator cvit = cvs->second.begin(); cvit != cvs->second.end(); ++cvit)
        {
          sp.setMetaValue(cvs->first, cvit->getValue());
        }
      }
      int minCharge = 0;
      int maxCharge = 0;
      for (map<String, DataValue>::const_iterator upit = as_params.second.begin(); upit != as_params.second.end(); ++upit)
      {
        if (upit->first == "taxonomy")
        {
          sp.taxonomy = upit->second.toString();
        }
        else if (upit->first == "charges")
        {
          sp.charges = upit->second.toString();
        }
        else if (upit->first == "MinCharge")
        {
          minCharge = upit->second.toString().toInt();
        }
        else if (upit->first == "MaxCharge")
        {
          maxCharge = upit->second.toString().toInt();
        }
        else if (upit->first == "NumTolerableTermini")
        {
          sp.enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>(upit->second.toString().toInt());
        }
        else
        {
          sp.setMetaValue(upit->first, upit->second);
        }
      }
      if (minCharge != 0 || maxCharge != 0) // this means "MinCharge" and "MaxCharge" get preference over "charges"
      {
        sp.charges = String(minCharge) + "-" + String(maxCharge);
      }
      return sp;
    }
 
} // namespace OpenMS   //namespace Internal
