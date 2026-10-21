// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Eugen Netz $
// $Authors: Mathias Walzer, Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MzIdentMLDOMHandler.h>
#include <OpenMS/CONCEPT/Init.h>

#include <OpenMS/FORMAT/HANDLERS/XMLHandler.h>
#include <OpenMS/SYSTEM/File.h>

#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueDB.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/ANALYSIS/ID/Scores.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLHelper.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/ProteinRunReferences.h>

#include <boost/lexical_cast.hpp>

#include <sys/stat.h>
#include <algorithm>
#include <cerrno>
#include <limits>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <libxml/xmlsave.h>
#include <memory>

using namespace std;

namespace OpenMS::Internal
{

  namespace
  {
    // The document owns attached nodes. Remember newly allocated writer nodes as
    // well, so an exception before append cannot leak a detached subtree.
    struct NativeDocument
    {
      xmlDoc* doc = xmlNewDoc(BAD_CAST "1.0");
      std::vector<xmlNode*> allocated;
      NativeDocument()
      {
        if (!doc) { throw std::bad_alloc(); }
        doc->_private = this;
      }
      ~NativeDocument()
      {
        for (xmlNode* node : allocated)
        {
          if (!node->parent) { xmlAddChild(reinterpret_cast<xmlNode*>(doc), node); }
        }
        xmlFreeDoc(doc);
      }
      NativeDocument(const NativeDocument&) = delete;
      NativeDocument& operator=(const NativeDocument&) = delete;
      xmlNode* own(xmlNode* node)
      {
        if (!node) { throw std::bad_alloc(); }
        try { allocated.push_back(node); }
        catch (...) { xmlFreeNode(node); throw; }
        return node;
      }
    };

    std::string name_(const xmlNode* node)
    {
      if (!node || !node->name) { return {}; }
      // Prefixes are arbitrary, but extension namespaces must not impersonate
      // mzIdentML elements. Retain support for historical unqualified files.
      if (node->ns && node->ns->href && *node->ns->href &&
          !xmlStrEqual(node->ns->href, BAD_CAST "http://psidev.info/psi/pi/mzIdentML/1.0") &&
          !xmlStrEqual(node->ns->href, BAD_CAST "http://psidev.info/psi/pi/mzIdentML/1.1") &&
          !xmlStrEqual(node->ns->href, BAD_CAST "http://psidev.info/psi/pi/mzIdentML/1.2") &&
          !xmlStrEqual(node->ns->href, BAD_CAST "http://psidev.info/psi/pi/mzIdentML/1.3")) { return {}; }
      return reinterpret_cast<const char*>(node->name);
    }

    xmlNode* nextElement_(xmlNode* node)
    {
      node = xmlNextElementSibling(node);
      while (node && name_(node).empty()) { node = xmlNextElementSibling(node); }
      return node;
    }

    xmlNode* firstElement_(xmlNode* node)
    {
      node = xmlFirstElementChild(node);
      if (node && name_(node).empty()) { return nextElement_(node); }
      return node;
    }

    bool boolean_(std::string value)
    {
      StringUtils::trim(value);
      return value == "true" || value == "1";
    }

    bool equal_(const std::string& left, const std::string& right) { return left == right; }

    std::string attribute_(const xmlNode* node, const char* name)
    {
      std::unique_ptr<xmlChar, decltype(xmlFree)> value(xmlGetNoNsProp(node, BAD_CAST name), xmlFree);
      return value ? reinterpret_cast<const char*>(value.get()) : "";
    }

    bool hasAttribute_(const xmlNode* node, const char* name) { return xmlHasNsProp(node, BAD_CAST name, nullptr) != nullptr; }

    std::string content_(const xmlNode* node)
    {
      std::unique_ptr<xmlChar, decltype(xmlFree)> value(xmlNodeGetContent(node), xmlFree);
      return value ? reinterpret_cast<const char*>(value.get()) : "";
    }

    std::vector<xmlNode*> children_(const xmlNode* node)
    {
      std::vector<xmlNode*> result;
      for (xmlNode* child = node->children; child; child = child->next) { result.push_back(child); }
      return result;
    }

    // Preserve descendant document order while excluding the starting element.
    std::vector<xmlNode*> descendants_(const xmlNode* parent, const std::string& name)
    {
      std::vector<xmlNode*> result;
      xmlNode* node = parent->children;
      while (node)
      {
        if (node->type == XML_ELEMENT_NODE && name_(node) == name) { result.push_back(node); }
        if (node->children) { node = node->children; }
        else
        {
          while (node && node != parent && !node->next) { node = node->parent; }
          if (!node || node == parent) { break; }
          node = node->next;
        }
      }
      return result;
    }

    std::vector<xmlNode*> descendants_(const xmlDoc* doc, const std::string& name)
    {
      // DTD entity declarations can themselves own parsed element trees. They
      // describe replacement content and must not be counted as document data.
      xmlNode* root = xmlDocGetRootElement(doc);
      if (!root) { return {}; }
      auto result = descendants_(root, name);
      if (name_(root) == name) { result.insert(result.begin(), root); }
      return result;
    }

    xmlNode* element_(xmlDoc* doc, const std::string& name)
    {
      xmlNode* root = xmlDocGetRootElement(doc);
      return static_cast<NativeDocument*>(doc->_private)->own(
        xmlNewDocNode(doc, root ? root->ns : nullptr, BAD_CAST name.c_str(), nullptr));
    }

    void setText_(xmlNode* parent, const std::string& value)
    {
      xmlNode* text = xmlNewDocText(parent->doc, BAD_CAST value.c_str());
      if (!text) { throw std::bad_alloc(); }
      if (!xmlAddChild(parent, text)) { xmlFreeNode(text); throw std::bad_alloc(); }
    }

    void append_(xmlNode* parent, xmlNode* child)
    {
      xmlUnlinkNode(child);
      if (!xmlAddChild(parent, child)) { throw std::bad_alloc(); }
    }

    void setAttribute_(xmlNode* node, const std::string& name, const std::string& value)
    {
      xmlAttr* attribute;
      if (name == "xsi:schemaLocation")
      {
        xmlNs* ns = xmlSearchNs(node->doc, node, BAD_CAST "xsi");
        attribute = xmlSetNsProp(node, ns, BAD_CAST "schemaLocation", BAD_CAST value.c_str());
      }
      else { attribute = xmlSetProp(node, BAD_CAST name.c_str(), BAD_CAST value.c_str()); }
      if (!attribute) { throw std::bad_alloc(); }
    }

    // (experimental|calculated)MassToCharge and retention time are optional mzIdentML attributes/values; a
    // missing one is read as an empty string. For these position-like fields convert tolerantly
    // so an absent value becomes NaN instead of throwing Exception::ConversionError while parsing (e.g. on
    // store/load round-trips); NaN is the OpenMS sentinel for "no RT/MZ" (cf. PeptideIdentification::hasRT()).
    inline double toDoubleOrNaN_(const std::string& s)
    {
      return s.empty() ? std::numeric_limits<double>::quiet_NaN() : StringUtils::toDouble(s);
    }

    // For score values an empty value must NOT become NaN: PeptideHit has no "hasScore()" and defaults the
    // score to 0.0, and a NaN score silently breaks sorting/FDR/IDFilter (every NaN comparison is false).
    // Default to 0.0 to match PeptideHit's contract and the other identification parsers.
    inline double toDoubleOrZero_(const std::string& s)
    {
      return s.empty() ? 0.0 : StringUtils::toDouble(s);
    }

    // Search engine specific scores (MS:1001143 and children) count as higher-is-better unless the PSI-MS
    // vocabulary annotates them otherwise (has_order) or, for the many terms without such an annotation,
    // OpenMS' score registry knows them as e-values, PEPs or q-values (e.g. Comet:expectation value).
    bool isLowerBetterEngineScore_(const ControlledVocabulary::CVTerm& term)
    {
      const bool has_order = std::any_of(term.unparsed.begin(), term.unparsed.end(),
        [](const std::string& line) { return StringUtils::hasPrefix(line, "relationship: has_order"); });
      if (has_order)
      {
        return !ControlledVocabulary::CVTerm::isHigherBetterScore(term);
      }
      Scores::IDType type = Scores::IDType::RAW;
      return Scores::findIDTypeByName(term.id, type) && !Scores::isHigherBetter(type);
    }
  }

    //TODO remodel CVTermList
    //TODO extend CVTermlist with CVCollection functionality for complete replacement??
    //TODO general id openms struct for overall parameter for one id run
    MzIdentMLDOMHandler::MzIdentMLDOMHandler(const vector<ProteinIdentification>& pro_id, const PeptideIdentificationList& pep_id, const std::string& version, const ProgressLogger& logger) :
      logger_(logger),
      cv_(ControlledVocabulary::getPSIMSCV()),
      //~ ms_exp_(0),
      cpro_id_(&pro_id),
      cpep_id_(&pep_id),
      schema_version_(version)
    {
      libxmlInitialize();
      unimod_.loadFromOBO("UNIMOD", File::find("/CV/unimod.obo"));

      initScoreTermCaches_();
    }

    MzIdentMLDOMHandler::MzIdentMLDOMHandler(vector<ProteinIdentification>& pro_id, PeptideIdentificationList& pep_id, const std::string& version, const ProgressLogger& logger) :
      logger_(logger),
      cv_(ControlledVocabulary::getPSIMSCV()),
      //~ ms_exp_(0),
      pro_id_(&pro_id),
      pep_id_(&pep_id),
      schema_version_(version),
      xl_ms_search_(false)
    {
      libxmlInitialize();
      unimod_.loadFromOBO("UNIMOD", File::find("/CV/unimod.obo"));

      initScoreTermCaches_();
    }

    void MzIdentMLDOMHandler::initScoreTermCaches_()
    {
      // precompute CV child term sets, used per PSM, but constant across the file
      cv_.addAllChildTerms(q_score_child_terms_, "MS:1002354");
      cv_.addAllChildTerms(e_score_child_terms_, "MS:1001872");
      cv_.addAllChildTerms(e_score_child_terms_, "MS:1002353");
      cv_.addAllChildTerms(specific_score_child_terms_, "MS:1001143");
    }

    MzIdentMLDOMHandler::~MzIdentMLDOMHandler() = default;

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

      std::unique_ptr<xmlParserCtxt, decltype(&xmlFreeParserCtxt)> parser(xmlNewParserCtxt(), xmlFreeParserCtxt);
      if (!parser) { throw std::bad_alloc(); }
      // libxml2 retains the diagnostic on this context. Suppress its default stderr
      // reporter without installing process-global handlers or throwing through C.
      xmlCtxtSetErrorHandler(parser.get(), [](void*, const xmlError*) noexcept {}, nullptr);
      xmlCtxtSetMaxAmplification(parser.get(), 5);
      std::unique_ptr<xmlDoc, decltype(&xmlFreeDoc)> document(
        xmlCtxtReadFile(parser.get(), mzid_file.c_str(), nullptr,
          XML_PARSE_NONET | XML_PARSE_NO_XXE | XML_PARSE_NOENT | XML_PARSE_HUGE), xmlFreeDoc);
      if (!document)
      {
        const xmlError* error = xmlCtxtGetLastError(parser.get());
        std::string message = error && error->message ? error->message : "Failed to parse mzIdentML file.";
        if (error)
        {
          message += " (line " + std::to_string(error->line) + ", column " + std::to_string(error->int2) + ")";
        }
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mzid_file, message);
      }
      auto* xmlDoc = document.get();
      if (!xmlDocGetRootElement(xmlDoc) || name_(xmlDocGetRootElement(xmlDoc)) != "MzIdentML")
      {
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mzid_file, "Expected MzIdentML root element.");
      }
      {
        // Catch special case: Cross-Linking MS
        std::vector<xmlNode*> additionalSearchParams = descendants_(xmlDoc, "AdditionalSearchParams");
        const  Size as_node_count = additionalSearchParams.size();

        for (Size i = 0; i < as_node_count; ++i)
        {
          xmlNode* current_sp = additionalSearchParams[i];

          xmlNode* element_SearchParams = current_sp;
          std::string cross_linking_search = attribute_(element_SearchParams, "id");
          xmlNode* child = firstElement_(element_SearchParams);

          while (child && !xl_ms_search_)
          {
            std::string accession = attribute_(child, "accession");
            if (accession == "MS:1002494") // accession for "crosslinking search"
            {
              xl_ms_search_ = true;
            }
            child = nextElement_(child);
          }
        }

        if (xl_ms_search_)
        {
          OPENMS_LOG_DEBUG << "Reading a Cross-Linking MS file." << endl;
        }

        // 0. AnalysisSoftwareList {0,1}
        std::vector<xmlNode*> analysisSoftwareElements = descendants_(xmlDoc, "AnalysisSoftware");
        parseAnalysisSoftwareList_(analysisSoftwareElements);

        // 1. DataCollection {1,1}
        std::vector<xmlNode*> spectraDataElements = descendants_(xmlDoc, "SpectraData");
        if (spectraDataElements.size() == 0)
        {
          throw(runtime_error("No SpectraData nodes"));
        }
        parseInputElements_(spectraDataElements);

        // 1.2. SearchDatabase {0,unbounded}
        std::vector<xmlNode*> searchDatabaseElements = descendants_(xmlDoc, "SearchDatabase");
        parseInputElements_(searchDatabaseElements);

        // 1.1 SourceFile {0,unbounded}
        std::vector<xmlNode*> sourceFileElements = descendants_(xmlDoc, "SourceFile");
        parseInputElements_(sourceFileElements);

        // 2. SpectrumIdentification  {1,unbounded} ! creates identification runs (or ProteinIdentifications)
        std::vector<xmlNode*> spectrumIdentificationElements = descendants_(xmlDoc, "SpectrumIdentification");
        if (spectrumIdentificationElements.size() == 0)
        {
          throw(runtime_error("No SpectrumIdentification nodes"));
        }
        parseSpectrumIdentificationElements_(spectrumIdentificationElements);

        // 3. AnalysisProtocolCollection {1,1} SpectrumIdentificationProtocol  {1,unbounded} ! identification run parameters
        std::vector<xmlNode*> spectrumIdentificationProtocolElements = descendants_(xmlDoc, "SpectrumIdentificationProtocol");
        if (spectrumIdentificationProtocolElements.size() == 0)
        {
          throw(runtime_error("No SpectrumIdentificationProtocol nodes"));
        }
        parseSpectrumIdentificationProtocolElements_(spectrumIdentificationProtocolElements);

        // 4. SequenceCollection nodes {0,1} DBSequenceElement {1,unbounded} Peptide {0,unbounded} PeptideEvidence {0,unbounded}
        std::vector<xmlNode*> dbSequenceElements = descendants_(xmlDoc, "DBSequence");
        parseDBSequenceElements_(dbSequenceElements);

        std::vector<xmlNode*> peptideElements = descendants_(xmlDoc, "Peptide");
        parsePeptideElements_(peptideElements);

        std::vector<xmlNode*> peptideEvidenceElements = descendants_(xmlDoc, "PeptideEvidence");
        parsePeptideEvidenceElements_(peptideEvidenceElements);


        // 5. AnalysisSampleCollection ??? contact stuff

        // 6. AnalysisCollection {1,1} - build final structures PeptideIdentification (and hits)

        // 6.1 SpectrumIdentificationList {0,1}
        std::vector<xmlNode*> spectrumIdentificationListElements = descendants_(xmlDoc, "SpectrumIdentificationList");
        if (spectrumIdentificationListElements.size() == 0)
        {
          throw(runtime_error("No SpectrumIdentificationList nodes"));
        }
        parseSpectrumIdentificationListElements_(spectrumIdentificationListElements);

        // 6.2 ProteinDetection {0,1}
        std::vector<xmlNode*> parseProteinDetectionListElements = descendants_(xmlDoc, "ProteinDetectionList");
        parseProteinDetectionListElements_(parseProteinDetectionListElements);

        for (vector<ProteinIdentification>::iterator it = pro_id_->begin(); it != pro_id_->end(); ++it)
        {
          it->sort();
        }
        //note: PeptideIdentification sorting here not necessary any more, due to sorting according to cv in SpectrumIdentificationResult

      }
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

      NativeDocument document;
      auto* xmlDoc = document.doc;
      xmlNode* root = xmlNewDocNode(xmlDoc, nullptr, BAD_CAST "MzIdentML", nullptr);
      if (!root) { throw std::bad_alloc(); }
      xmlDocSetRootElement(xmlDoc, root);
      xmlNs* mzid_ns = xmlNewNs(root, BAD_CAST "http://psidev.info/psi/pi/mzIdentML/1.1", nullptr);
      if (!mzid_ns || !xmlNewNs(root, BAD_CAST "http://www.w3.org/2001/XMLSchema-instance", BAD_CAST "xsi")) { throw std::bad_alloc(); }
      xmlSetNs(root, mzid_ns);
          xmlNode* rootElem = xmlDocGetRootElement(xmlDoc);
          setAttribute_(rootElem, "version",
                                 schema_version_);
          setAttribute_(rootElem, "xsi:schemaLocation",
                                 "http://psidev.info/psi/pi/mzIdentML/1.1 ../../schema/mzIdentML1.1.0.xsd");
          setAttribute_(rootElem, "creationDate",
                                 std::string(DateTime::now().getDate() + "T" + DateTime::now().getTime()));

          // * cvList *
          xmlNode* cvl_p = element_(xmlDoc, "cvList"); // TODO add generically
          buildCvList_(cvl_p);
          append_(rootElem, cvl_p);

          // * AnalysisSoftwareList *
          xmlNode* asl_p = element_(xmlDoc, "AnalysisSoftwareList");
          for (vector<ProteinIdentification>::const_iterator pi = cpro_id_->begin(); pi != cpro_id_->end(); ++pi)
          {
//                  search_engine_version_ = pi->getSearchEngineVersion();
//                  search_engine_ = pi->getSearchEngine();
          }
          buildAnalysisSoftwareList_(asl_p);
          append_(rootElem, asl_p);

//              // * AnalysisSampleCollection *
//              xmlNode* asc_p = element_(xmlDoc, "AnalysisSampleCollection");
//              buildAnalysisSampleCollection_(asc_p);
//              append_(rootElem, asc_p);

          // * SequenceCollection *
          xmlNode* sc_p = element_(xmlDoc, "SequenceCollection");

          for (vector<ProteinIdentification>::const_iterator pi = cpro_id_->begin(); pi != cpro_id_->end(); ++pi)
          {
            std::string dbref = pi->getSearchParameters().db +  pi->getSearchParameters().db_version + pi->getSearchParameters().taxonomy; //TODO @mths : this needs to be more unique, btw add tax etc. as cv to DBSequence
            for (vector<ProteinHit>::const_iterator ph = pi->getHits().begin(); ph != pi->getHits().end(); ++ph)
            {
              CVTermList cvs;
              DBSequence temp_struct = {ph->getSequence(), dbref, ph->getAccession(), cvs};
              db_sq_map_.insert(make_pair(ph->getAccession(), temp_struct));
            }
          }

          set<AASequence> pepset;
          for (PeptideIdentificationList::const_iterator pi = cpep_id_->begin(); pi != cpep_id_->end(); ++pi)
          {
            for (vector<PeptideHit>::const_iterator ph = pi->getHits().begin(); ph != pi->getHits().end(); ++ph)
            {
              list<std::string> pepevs;
              for (vector<OpenMS::PeptideEvidence>::const_iterator pev = ph->getPeptideEvidences().begin(); pev != ph->getPeptideEvidences().end(); ++pev)
              {
                std::string pepevref ="OpenMS" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
                pv_db_map_.insert(make_pair(pepevref, pev->getProteinAccession()));
                pepevs.push_back(pepevref);
                bool idec = StringUtils::hasSubstring(StringUtils::toStr(ph->getMetaValue(Constants::UserParam::TARGET_DECOY)), "decoy");
                PeptideEvidence temp_struct = {pev->getStart(), pev->getEnd(), pev->getAABefore(), pev->getAAAfter(), idec}; //TODO @ mths : completely switch to PeptideEvidence
                pe_ev_map_.insert(make_pair(pepevref, temp_struct)); // TODO @ mths : double check start & end & chars for before & after
              }
              hit_pev_.push_back(pepevs);

              std::string pepref ="OpenMS" + StringUtils::toStr(UniqueIdGenerator::getUniqueId());
              if (pepset.contains(ph->getSequence()))
              {
                pepset.insert(ph->getSequence());
                pep_map_.insert(make_pair(pepref, ph->getSequence()));
                for (list<std::string>::iterator pepevref = pepevs.begin(); pepevref != pepevs.end(); ++pepevref)
                {
                  p_pv_map_.insert(make_pair(*pepevref, pepref));
                }
              }
            }
          }

          buildSequenceCollection_(sc_p);
          append_(rootElem, sc_p);

          // * AnalysisCollection *
          xmlNode* analysis_c_p = element_(xmlDoc, "AnalysisCollection");
          buildAnalysisCollection_(analysis_c_p);
          append_(rootElem, analysis_c_p);

          // * AnalysisProtocolCollection *
          xmlNode* apc_p = element_(xmlDoc, "AnalysisProtocolCollection");
          buildAnalysisCollection_(apc_p);
          append_(rootElem, apc_p);

          // * DataCollection *
          xmlNode* dc_p = element_(xmlDoc, "DataCollection");
          append_(rootElem, dc_p);
          xmlNode* in_p = element_(dc_p->doc, "Inputs");
          xmlNode* ad_p = element_(dc_p->doc, "AnalysisData");
          append_(dc_p, in_p);
          append_(dc_p, ad_p);

          // * BibliographicReference *
          xmlNode* br_p = element_(xmlDoc, "BibliographicReference");
          setAttribute_(br_p, "authors", "all");
          append_(rootElem, br_p);

      // Save as UTF8 and propagate open/write/close errors instead of logging success.
      xmlSaveCtxt* output = xmlSaveToFilename(mzid_file.c_str(), "UTF-8", XML_SAVE_FORMAT);
      if (!output) { throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mzid_file); }
      const long saved = xmlSaveDoc(output, xmlDoc);
      const int closed = xmlSaveClose(output);
      if (saved < 0 || closed < 0)
      {
        throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, mzid_file);
      }
    }

    pair<CVTermList, map<std::string, DataValue> > MzIdentMLDOMHandler::parseParamGroup_(const std::vector<xmlNode*>& paramGroup)
    {
      CVTermList ret_cv;             // cvParam
      map<std::string, DataValue> ret_up; // userParam
      const  Size cv_node_count = paramGroup.size();
      for (Size cvi = 0; cvi < cv_node_count; ++cvi)
      {
        xmlNode* current_cv = paramGroup[cvi];
        if (current_cv->type && // true is not NULL
            current_cv->type == XML_ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          xmlNode* element_param = current_cv;
          if (equal_(name_(element_param), "cvParam"))
          {
            ret_cv.addCVTerm(parseCvParam_(element_param));
          }
          else if (equal_(name_(element_param), "userParam"))
          {
            ret_up.insert(parseUserParam_(element_param));
          }
          else if (equal_(name_(element_param), "PeptideEvidence")
                  || equal_(name_(element_param), "PeptideEvidenceRef")
                  || equal_(name_(element_param), "Fragmentation")
                  || equal_(name_(element_param), "SpectrumIdentificationItem"))
          {
            //here it's okay to do nothing
          }
          else
          {
            OPENMS_LOG_WARN << "Misplaced elements ignored in 'ParamGroup' in " << name_(element_param) << endl;
          }
        }
      }
      return make_pair(ret_cv, ret_up);
    }

    CVTerm MzIdentMLDOMHandler::parseCvParam_(xmlNode* param)
    {
      if (param)
      {
        //      <cvParam accession="MS:1001469" name="taxonomy: scientific name" cvRef="PSI-MS"  value="Drosophila melanogaster"/>
        std::string accession = attribute_(param, "accession");
        std::string name = attribute_(param, "name");
        std::string cvRef = attribute_(param, "cvRef");
        std::string value = attribute_(param, "value");

        std::string unitAcc = attribute_(param, "unitAccession");
        std::string unitName = attribute_(param, "unitName");
        std::string unitCvRef = attribute_(param, "unitCvRef");

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

    pair<std::string, DataValue> MzIdentMLDOMHandler::parseUserParam_(xmlNode* param)
    {
      if (param)
      {
        //      <userParam name="Mascot User Comment" value="Example Mascot MS-MS search for PSI mzIdentML"/>
        std::string name = attribute_(param, "name");
        std::string value = attribute_(param, "value");
        bool has_value = hasAttribute_(param, "value");
        std::string unitAcc = attribute_(param, "unitAccession");
        std::string unitName = attribute_(param, "unitName");
        std::string unitCvRef = attribute_(param, "unitCvRef");
        std::string type = attribute_(param, "type");

        DataValue dv = DataValue::EMPTY;
        if (has_value)
        {
          dv = XMLHandler::fromXSDString(type, value);
        }

        // Add unit *after* creating the term
        if (!unitAcc.empty())
        {
          if (StringUtils::hasPrefix(unitAcc, "UO:"))
          {
            dv.setUnit(StringUtils::toInt32(StringUtils::suffix(unitAcc, unitAcc.size() - 3)));
            dv.setUnitType(DataValue::UnitType::UNIT_ONTOLOGY);
          }
          else if (StringUtils::hasPrefix(unitAcc, "MS:"))
          {
            dv.setUnit(StringUtils::toInt32(StringUtils::suffix(unitAcc, unitAcc.size() - 3)));
            dv.setUnitType(DataValue::UnitType::MS_ONTOLOGY);
          }
          else
          {
            OPENMS_LOG_WARN << std::string("Unhandled unit '") + unitAcc + "' in tag '" + name + "'." << endl;
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

    void MzIdentMLDOMHandler::parseAnalysisSoftwareList_(const std::vector<xmlNode*>& analysisSoftwareElements)
    {
      const  Size as_node_count = analysisSoftwareElements.size();
      for (Size swni = 0; swni < as_node_count; ++swni)
      {
        xmlNode* current_as = analysisSoftwareElements[swni];
        if (current_as->type && // true is not NULL
            current_as->type == XML_ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          xmlNode* element_AnalysisSoftware = current_as;
          std::string id = attribute_(element_AnalysisSoftware, "id");
          xmlNode* child = firstElement_(element_AnalysisSoftware);
          std::string swname, swversion;
          while (child)
          {
            if (equal_(name_(child),"SoftwareName")) //must have exactly one SoftwareName
            {
              std::vector<xmlNode*> element_pg = children_(child);

              pair<CVTermList, map<std::string, DataValue> > swn = parseParamGroup_(element_pg);
              swversion = attribute_(element_AnalysisSoftware, "version");
              if (!swn.first.getCVTerms().empty())
              {
                set<std::string> software_terms;
                cv_.getAllChildTerms(software_terms, "MS:1000531");
                for (map<std::string, vector<CVTerm> >::const_iterator it = swn.first.getCVTerms().begin(); it != swn.first.getCVTerms().end(); ++it)
                {
                  if (software_terms.contains(it->first))
                  {
                    swname = it->second.front().getName();
                    break;
                  }
                }
              }
              else if (!swn.second.empty())
              {
                for (map<std::string, DataValue>::const_iterator up = swn.second.begin(); up != swn.second.end(); ++up)
                {
                  if (StringUtils::hasSubstring(up->first, "name"))
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
            child = nextElement_(child);
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

    void MzIdentMLDOMHandler::parseDBSequenceElements_(const std::vector<xmlNode*>& dbSequenceElements)
    {
      const  Size dbs_node_count = dbSequenceElements.size();
      for (Size c = 0; c < dbs_node_count; ++c)
      {
        xmlNode* current_dbs = dbSequenceElements[c];
        if (current_dbs->type && // true is not NULL
            current_dbs->type == XML_ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          xmlNode* element_dbs = current_dbs;
          std::string id = attribute_(element_dbs, "id");
          std::string seq;
          std::string dbref = attribute_(element_dbs, "searchDatabase_ref");
          std::string acc = attribute_(element_dbs, "accession");
          CVTermList cvs;

          xmlNode* child = firstElement_(element_dbs);
          while (child)
          {
            if (equal_(name_(child), "Seq"))
            {
              seq = content_(child);
            }
            else if (equal_(name_(child), "cvParam"))
            {
              cvs.addCVTerm(parseCvParam_(child));
            }
            child = nextElement_(child);
          }
          if (!acc.empty())
          {
            DBSequence temp_struct = {seq, dbref, acc, cvs};
            db_sq_map_.insert(make_pair(id, temp_struct));
          }
        }
      }
    }

    void MzIdentMLDOMHandler::parsePeptideElements_(const std::vector<xmlNode*>& peptideElements)
    {
      const  Size pep_node_count = peptideElements.size();
      for (Size c = 0; c < pep_node_count; ++c)
      {
        xmlNode* current_pep = peptideElements[c];
        if (current_pep->type && // true is not NULL
            current_pep->type == XML_ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          xmlNode* element_pep = current_pep;
          std::string id = attribute_(element_pep, "id");

          //std::vector<xmlNode*> pep_sib = children_(element_pep);
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
              std::string name = attribute_(element_pep, "name");
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

    void MzIdentMLDOMHandler::parsePeptideEvidenceElements_(const std::vector<xmlNode*>& peptideEvidenceElements)
    {
      const  Size pev_node_count = peptideEvidenceElements.size();
      for (Size c = 0; c < pev_node_count; ++c)
      {
        xmlNode* current_pev = peptideEvidenceElements[c];
        if (current_pev->type && // true is not NULL
            current_pev->type == XML_ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          xmlNode* element_pev = current_pev;

//          <PeptideEvidence peptide_ref="peptide_1_1" id="PE_1_1_HSP70_ECHGR_0" start="161" end="172" pre="K" post="I" isDecoy="false" dBSequence_ref="DBSeq_HSP70_ECHGR"/>

          std::string id = attribute_(element_pev, "id");
          std::string peptide_ref = attribute_(element_pev, "peptide_ref");
          std::string dBSequence_ref = attribute_(element_pev, "dBSequence_ref");
          //rest is optional !!
          int start = -1;
          int end = -1;
          try
          {
            start = StringUtils::toInt32(attribute_(element_pev, "start"));
            end = StringUtils::toInt32(attribute_(element_pev, "end"));
          }
          catch (...)
          {
            OPENMS_LOG_WARN << "'PeptideEvidence' without reference to the position in the originating sequence found." << endl;
          }
          char pre = '-';
          char post = '-';
          try
          {
            if (hasAttribute_(element_pev, "pre"))
            {
            pre = attribute_(element_pev, "pre")[0];
            }
            if (hasAttribute_(element_pev, "post"))
            {
            post = attribute_(element_pev, "post")[0];
          }
          }
          catch (...)
          {
            OPENMS_LOG_WARN << "'PeptideEvidence' without reference to the bordering amino acids in the originating sequence found." << endl;
          }
          bool idec = false;
          try
          {
            std::string d = attribute_(element_pev, "isDecoy");
            if (StringUtils::hasPrefix(d, 't') || StringUtils::hasPrefix(d, '1'))
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

    void MzIdentMLDOMHandler::parseSpectrumIdentificationElements_(const std::vector<xmlNode*>& spectrumIdentificationElements)
    {
      const  Size si_node_count = spectrumIdentificationElements.size();
      for (Size c = 0; c < si_node_count; ++c)
      {
        xmlNode* current_si = spectrumIdentificationElements[c];
        if (current_si->type && // true is not NULL
            current_si->type == XML_ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          xmlNode* element_si = current_si;
          std::string id = attribute_(element_si, "id");
          std::string spectrumIdentificationProtocol_ref = attribute_(element_si, "spectrumIdentificationProtocol_ref");
          std::string spectrumIdentificationList_ref = attribute_(element_si, "spectrumIdentificationList_ref");
          std::string spectrumIdentification_date = attribute_(element_si, "activityDate");

          std::string searchDatabase_ref;
          std::string spectra_data_ref;
          xmlNode* child = firstElement_(element_si);
          while (child)
          {
            if (equal_(name_(child), "InputSpectra"))
            {
              spectra_data_ref = attribute_(child, "spectraData_ref");
            }
            else if (equal_(name_(child), "SearchDatabaseRef"))
            {
              searchDatabase_ref = attribute_(child, "searchDatabase_ref");
            }
            child = nextElement_(child);
          }
          SpectrumIdentification temp_struct = {spectra_data_ref, searchDatabase_ref, spectrumIdentificationProtocol_ref, spectrumIdentificationList_ref};
          si_map_.insert(make_pair(id, temp_struct));

          pro_id_->emplace_back();
          ProteinIdentification::SearchParameters sp;
          sp.db = db_map_[searchDatabase_ref].location;
          sp.db_version = db_map_[searchDatabase_ref].version;
          pro_id_->back().setSearchParameters(sp);
          if (xl_ms_search_)
          {
            pro_id_->back().setMetaValue("SpectrumIdentificationProtocol", "MS:1002494"); // XL-MS CV term
          }

          // internally we store a list of files so convert the mzIdentML file std::string to a StringList
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
          pro_id_->back().setIdentifier(StringUtils::toStr(UniqueIdGenerator::getUniqueId())); // no more identification wit engine/date/time!
          //TODO setIdentifier to xml id?
          si_pro_map_.insert(make_pair(spectrumIdentificationList_ref, pro_id_->size() - 1));
        }
      }
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationProtocolElements_(const std::vector<xmlNode*>& spectrumIdentificationProtocolElements)
    {
      const  Size si_node_count = spectrumIdentificationProtocolElements.size();
      for (Size c = 0; c < si_node_count; ++c)
      {
        ProteinIdentification::SearchParameters sp;
        xmlNode* current_sip = spectrumIdentificationProtocolElements[c];
        if (current_sip->type && // true is not NULL
            current_sip->type == XML_ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          xmlNode* element_sip = current_sip;
          std::string id = attribute_(element_sip, "id");
          std::string swr = attribute_(element_sip, "analysisSoftware_ref");

          CVTerm searchtype;
          std::string enzymename;
          CVTermList param_cv;
          map<std::string, DataValue> param_up;
          CVTermList modparam;
          double p_tol = 0;
          double f_tol = 0;
          CVTermList tcv;
          map<std::string, DataValue> tup;
          xmlNode* child = firstElement_(element_sip);
          while (child)
          {
            if (equal_(name_(child), "SearchType"))
            {
              searchtype = parseCvParam_(firstElement_(child));
            }
            else if (equal_(name_(child), "AdditionalSearchParams"))
            {
              pair<CVTermList, map<std::string, DataValue> > as_params = parseParamGroup_(children_(child));
              sp = findSearchParameters_(as_params); // this must be renamed!!
            }
            else if (equal_(name_(child), "ModificationParams")) // TODO @all where to store the specificities?
            {
              vector<std::string> fix, var;
              xmlNode* sm = firstElement_(child);
              while (sm)
              {
                std::string residues = attribute_(sm, "residues");
                const bool fixedMod = boolean_(attribute_(sm, "fixedMod"));

//                double massDelta = 0;
//                try
//                {
//                  massDelta = boost::lexical_cast<double>(attribute_(sm, "massDelta"));
//                }
//                catch (...)
//                {
//                    OPENMS_LOG_ERROR << "Could not cast ModificationParam massDelta from " << attribute_(sm, "massDelta");
//                }

                std::string mname;
                CVTermList specificity_rules;
                xmlNode* sub = firstElement_(sm);
                while (sub)
                {
                  if (equal_(name_(sub), "cvParam"))
                  {
                    // In v1.3, SearchModification may have multiple cvParams: the actual modification
                    // plus PSI-MS annotations (MS:1003392 "search modification id",
                    // MS:1003347 "UNIMOD derivative code", MS:1002509/MS:1002510 donor/acceptor)
                    std::string cvref = attribute_(sub, "cvRef");
                    if (cvref == "UNIMOD" || cvref == "XLMOD")
                    {
                      mname = attribute_(sub, "name");
                      if (mname == "unknown modification")
                      {
                        // e.g. <cvParam cvRef="MS" accession="MS:1001460" name="unknown modification" value="N-Glycan"/>
                        mname = attribute_(sub, "value");
                      }
                    }
                    else if (mname.empty())
                    {
                      // fallback for files that may not specify cvRef
                      mname = attribute_(sub, "name");
                      if (mname == "unknown modification")
                      {
                        mname = attribute_(sub, "value");
                      }
                    }
                  }
                  else if (equal_(name_(sub), "SpecificityRules"))
                  {
                    specificity_rules.consumeCVTerms(parseParamGroup_(children_(sub)).first.getCVTerms());
                    // let's press them where all other SearchengineAdapters press them in
                  }
                  else
                  {
                    OPENMS_LOG_ERROR << "Misplaced information in 'ModificationParams' ignored." << endl;
                  }
                  sub = nextElement_(sub);
                }

                if (!mname.empty())
                {
                  std::string mod;
                  std::string r = (residues != ".") ? residues : "";

                  try
                  {
                    if (!specificity_rules.empty())
                    {
                      for (map<std::string, vector<CVTerm> >::const_iterator spci = specificity_rules.getCVTerms().begin(); spci != specificity_rules.getCVTerms().end(); ++spci)
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
                  }
                  catch (Exception::BaseException& e)
                  {
                    OPENMS_LOG_WARN << "SearchModification '" << mname << "' not found in ModificationsDB"
                                    << e.what() << endl;
                  }

                  if (!mod.empty() && fixedMod)
                  {
                    fix.push_back(mod);
                  }
                  else if (!mod.empty())
                  {
                    var.push_back(mod);
                  }
                }
                sm = nextElement_(sm);
              }
              sp.fixed_modifications = fix;
              sp.variable_modifications = var;
            }
            else if (equal_(name_(child), "Enzymes")) // TODO @all : where store multiple enzymes for one identificationrun?
            {
              xmlNode* enzyme = firstElement_(child); //Enzyme elements
              while (enzyme)
              {
                int missedCleavages = -1;
                try
                {
                  missedCleavages = StringUtils::toInt32(attribute_(enzyme, "missedCleavages"));
                }
                catch (exception& e)
                {
                  OPENMS_LOG_WARN << "Search engine enzyme settings for 'missedCleavages' unreadable: " << e.what() << attribute_(enzyme, "missedCleavages") << endl;
                }
                if (missedCleavages < 0)
                {
                  OPENMS_LOG_WARN << "missedCleavages has a negative value. Assuming unlimited and setting it to 1000." << endl;
                  missedCleavages = 1000;
                }
                sp.missed_cleavages = missedCleavages;

//                std::string semiSpecific = attribute_(enzyme, "semiSpecific"); //xsd:boolean
//                std::string cTermGain = attribute_(enzyme, "cTermGain");
//                std::string nTermGain = attribute_(enzyme, "nTermGain");
//                int minDistance = -1;
//                try
//                {
//                  minDistance = StringUtils::toInt32(attribute_(enzyme, "minDistance"));
//                }
//                catch (...)
//                {
//                    OPENMS_LOG_WARN << "Search engine settings for 'minDistance' unreadable." << endl;
//                }
                enzymename = "UNKNOWN";
                xmlNode* sub = firstElement_(enzyme);
                while (sub)
                {
                  //SiteRegex unstorable just now
                  if (equal_(name_(sub), "EnzymeName"))
                  {
                    set<std::string> enzymes_terms;
                    cv_.getAllChildTerms(enzymes_terms, "MS:1001045"); // cleavage agent name
                    pair<CVTermList, map<std::string, DataValue> > params = parseParamGroup_(children_(sub));
                    for (map<std::string, vector<CVTerm> >::const_iterator it = params.first.getCVTerms().begin(); it != params.first.getCVTerms().end(); ++it)
                    {
                      if (enzymes_terms.contains(it->first))
                      {
                        enzymename = it->second.front().getName();
                      }
                      else
                      {
                        OPENMS_LOG_WARN << "Additional parameters for enzyme settings not readable." << endl;
                      }
                    }
                  }
                  sub = nextElement_(sub);
                }
                if (ProteaseDB::getInstance()->hasEnzyme(enzymename))
                {
                  sp.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzymename));
                }
                enzyme = nextElement_(enzyme);
              }
            }
            else if (equal_(name_(child), "FragmentTolerance"))
            {
              pair<CVTermList, map<std::string, DataValue> > params = parseParamGroup_(children_(child));
              //+- take the numerically greater
              for (map<std::string, vector<CVTerm> >::const_iterator it = params.first.getCVTerms().begin(); it != params.first.getCVTerms().end(); ++it)
              {
                f_tol = max(f_tol, StringUtils::toDouble(it->second.front().getValue().toString()));
                sp.fragment_mass_tolerance = f_tol;
                if (it->second.front().getUnit().name == "parts per million" )
                {
                  sp.fragment_mass_tolerance_ppm = true;
                }
              }
            }
            else if (equal_(name_(child), "ParentTolerance"))
            {
              pair<CVTermList, map<std::string, DataValue> > params = parseParamGroup_(children_(child));
              //+- take the numerically greater
              for (map<std::string, vector<CVTerm> >::const_iterator it = params.first.getCVTerms().begin(); it != params.first.getCVTerms().end(); ++it)
              {
                p_tol = max(p_tol, StringUtils::toDouble(it->second.front().getValue().toString()));
                sp.precursor_mass_tolerance = p_tol;
                if (it->second.front().getUnit().name == "parts per million" )
                {
                  sp.precursor_mass_tolerance_ppm = true;
                }

              }
            }
            else if (equal_(name_(child), "Threshold"))
            {
              pair<CVTermList, map<std::string, DataValue> > params = parseParamGroup_(children_(child));
              tcv = params.first;
              tup = params.second;
            }
            child = nextElement_(child);
            //      <DatabaseFilters> omitted for now, not reflectable by our member structures
            //      <DatabaseTranslation> omitted for now, not reflectable by our member structures
            //      <MassTable> omitted for now, not reflectable by our member structures
          }
          SpectrumIdentificationProtocol temp_struct = {searchtype, enzymename, param_cv, param_up, modparam, p_tol, f_tol, tcv, tup};
          sp_map_.insert(make_pair(id, temp_struct));

          double thresh = 0.0;
          bool use_thresh = false;
          set<std::string> threshold_terms;
          cv_.getAllChildTerms(threshold_terms, "MS:1002482"); //statistical threshold
          for (map<std::string, vector<OpenMS::CVTerm> >::const_iterator thit = tcv.getCVTerms().begin(); thit != tcv.getCVTerms().end(); ++thit)
          {
            if (threshold_terms.contains(thit->first))
            {
              if (thit->first != "MS:1001494") // no threshold
              {
                thresh = StringUtils::toDouble(thit->second.front().getValue().toString());
                use_thresh = true;
                break;
              }
              else
              {
                break;
              }
            }
          }

          std::string search_engine, search_engine_version;
          for (map<std::string, SpectrumIdentification>::const_iterator si_it = si_map_.begin(); si_it != si_map_.end(); ++si_it)
          {
            if (si_it->second.spectrum_identification_protocol_ref == id)
            {
              search_engine = as_map_[swr].name;
              search_engine_version = as_map_[swr].version;
//              std::string identi = search_engine+"_"+si_pro_map_[si_it->second.spectrum_identification_list_ref]->getDateTime().getDate()+"T"
//                      +si_pro_map_[si_it->second.spectrum_identification_list_ref]->getDateTime().getTime();
//              pro_id_->at(si_pro_map_[si_it->second.spectrum_identification_list_ref]).setIdentifier(identi);
              auto& pro = runOfList_(si_it->second.spectrum_identification_list_ref);
              pro.setSearchEngine(search_engine);
              pro.setSearchEngineVersion(search_engine_version);
              sp.db = pro.getSearchParameters().db; // was previously set, but main parts of sp are set here
              sp.db_version = pro.getSearchParameters().db_version; // was previously set, but main parts of sp are set here
              pro.setSearchParameters(sp);
              if (use_thresh)
              {
                pro.setSignificanceThreshold(thresh);
              }
            }
          }
        }
      }
    }

    void MzIdentMLDOMHandler::parseInputElements_(const std::vector<xmlNode*>& inputElements)
    {
      const  Size node_count = inputElements.size();
      for (Size c = 0; c < node_count; ++c)
      {
        xmlNode* current_in = inputElements[c];
        if (current_in->type && // true is not NULL
            current_in->type == XML_ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          xmlNode* element_in = current_in;

          std::string id = attribute_(element_in, "id");
          std::string location = attribute_(element_in, "location");

          if (equal_(name_(element_in), "SpectraData"))
          {
            //      <FileFormat> omitted for now, not reflectable by our member structures
            //      <SpectrumIDFormat> omitted for now, not reflectable by our member structures
            sd_map_.insert(make_pair(id, location));
          }
          else if (equal_(name_(element_in), "SourceFile"))
          {
            //      <FileFormat> omitted for now, not reflectable by our member structures
            sr_map_.insert(make_pair(id, location));
          }
          else if (equal_(name_(element_in), "SearchDatabase"))
          {
            //      <FileFormat> omitted for now, not reflectable by our member structures
            DateTime releaseDate;
//            releaseDate.set(attribute_(element_in, "releaseDate"));
            std::string version = attribute_(element_in, "version");
            std::string dbname;
            xmlNode* element_dbn = firstElement_(element_in);
            while (element_dbn)
            {
              if (equal_(name_(element_dbn), "DatabaseName"))
              {
                xmlNode* databasename_param = firstElement_(element_dbn);
                while (databasename_param)
                {
                  if (equal_(name_(databasename_param), "cvParam"))
                  {
                    CVTerm param = parseCvParam_(databasename_param);
                    dbname = StringUtils::toStr(param.getValue());
                  }
                  else if (equal_(name_(databasename_param), "userParam"))
                  {
                    pair<std::string, DataValue> param = parseUserParam_(databasename_param);
                    // issue #7099: mzID might have missing "value" for this element
                    // in this case, just use the "name"
                    dbname = param.second.isEmpty() ? dbname = param.first : param.second.toString();
                  }
                  databasename_param = nextElement_(databasename_param);
                }
                //each SearchDatabase element may have one DatabaseName, each DatabaseName only one param
              }
              element_dbn = nextElement_(element_dbn);
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

    ProteinIdentification& MzIdentMLDOMHandler::runOfList_(const std::string& spectrum_identification_list_ref)
    {
      const auto run = si_pro_map_.find(spectrum_identification_list_ref);
      if (run == si_pro_map_.end())
      {
        // no SpectrumIdentification (search run) references this list
        throw Exception::InvalidParameter(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
          ProteinRunReferences::missingRunMessage(spectrum_identification_list_ref) + " (SpectrumIdentificationList)");
      }
      return pro_id_->at(run->second);
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationListElements_(const std::vector<xmlNode*>& spectrumIdentificationListElements)
    {
      const  Size node_count = spectrumIdentificationListElements.size();
      for (Size c = 0; c < node_count; ++c)
      {
        xmlNode* current_lis = spectrumIdentificationListElements[c];
        if (current_lis->type && // true is not NULL
            current_lis->type == XML_ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          xmlNode* element_lis = current_lis;
          std::string id = attribute_(element_lis, "id");
//          std::string name = attribute_(element_res, "name");

          xmlNode* element_res = firstElement_(element_lis);
          while (element_res)
          {
            if (equal_(name_(element_res), "SpectrumIdentificationResult"))
            {
              std::string spectra_data_ref = attribute_(element_res, "spectraData_ref"); //ref to the sourcefile, could be useful but now nowhere to store
              std::string spectrumID = attribute_(element_res, "spectrumID");
              pair<CVTermList, map<std::string, DataValue> > params = parseParamGroup_(children_(element_res));

              if (xl_ms_search_) // XL-MS data has a different structure (up to 4 spectrum identification items for the same PSM)
              {
                std::multimap<std::string, int> xl_val_map;
                std::set<std::string> xl_val_set;
                int index_counter = 0;
                xmlNode* sii = firstElement_(element_res);

                // loop over all SIIs of a spectrum and group together the SIIs belonging to the same cross-link spectrum match
                while (sii)
                {
                  if (equal_(name_(sii), "SpectrumIdentificationItem"))
                  {
                    std::vector<xmlNode*> sii_cvp = descendants_(sii, "cvParam");
                    const  Size cv_count = sii_cvp.size();
                    for (Size i = 0; i < cv_count; ++i)
                    {
                      xmlNode* element_sii_cvp = sii_cvp[i];
                      if (equal_(attribute_(element_sii_cvp, "accession"), "MS:1002511")) // cross-link spectrum identification item
                      {
                        std::string xl_val = attribute_(element_sii_cvp, "value");
                        xl_val_map.insert(make_pair(xl_val, index_counter));
                        xl_val_set.insert(xl_val);
                      }
                    }
                  }
                  sii = nextElement_(sii);
                  ++index_counter;
                }

                // fix for label-free mono-links
                // those only have one SII and no "cross-link spectrum identification item" value
                if (xl_val_set.empty())
                {
                  xl_val_set.insert("0");
                  xl_val_map.insert(make_pair("0", 0));
                }
                // (each identification of a cross-link spectrum match gets the run there)
                for (set<std::string>::const_iterator set_it = xl_val_set.begin(); set_it != xl_val_set.end(); ++set_it)
                {
                  parseSpectrumIdentificationItemSetXLMS(set_it, xl_val_map, element_res, spectrumID);
                }
              }
              else // general case
              {
                pep_id_->push_back(PeptideIdentification());
                pep_id_->back().setHigherScoreBetter(false); //either a q-value or an e-value, only if neither available there will be another
                pep_id_->back().setSpectrumReference( spectrumID);  // SpectrumIdentificationResult attribute spectrumID is taken from the mz_file and should correspond to MSSpectrum.nativeID, thus spectrum_reference will serve as reference. As the format of the 'reference' widely varies from vendor to vendor, spectrum_reference as string will serve best, indices are not recommended.

                //fill pep_id_->back() with content
                xmlNode* parent = element_res->parent;
                std::string sil = attribute_(parent, "id");

                xmlNode* child = firstElement_(element_res);
                while (child)
                {
                  if (equal_(name_(child), "SpectrumIdentificationItem"))
                  {
                    parseSpectrumIdentificationItemElement_(child, pep_id_->back(), sil);
                  }
                  child = nextElement_(child);
                }

              } // end of "not-XLMS-results"

              pep_id_->back().setIdentifier(runOfList_(id).getIdentifier());

              pep_id_->back().sort();

              //adopt cv s
              for (map<std::string, vector<CVTerm> >::const_iterator cvit =  params.first.getCVTerms().begin(); cvit != params.first.getCVTerms().end(); ++cvit)
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
                  double rt = StringUtils::toDouble(cvit->second.front().getValue().toString());
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
              for (map<std::string, DataValue>::const_iterator upit = params.second.begin(); upit != params.second.end(); ++upit)
              {
                pep_id_->back().setMetaValue(upit->first, upit->second);
              }
              if (pep_id_->back().getRT() != pep_id_->back().getRT())
              {
                OPENMS_LOG_WARN << "No retention time found for 'SpectrumIdentificationResult'" << endl;
              }
            }
            element_res = nextElement_(element_res);
          }
        }
      }
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationItemSetXLMS(set<std::string>::const_iterator set_it, std::multimap<std::string, int> xl_val_map, xmlNode* element_res, const std::string& spectrumID)
    {
      // each value in the set corresponds to one PeptideIdentification object
      std::pair <std::multimap<std::string, int>::iterator, std::multimap<std::string, int>::iterator> range;
      range = xl_val_map.equal_range(*set_it);
      std::vector<xmlNode*> siis = descendants_(element_res, "SpectrumIdentificationItem");

      xmlNode* parent = element_res->parent;
      std::string spectrumIdentificationList_ref = attribute_(parent, "id");

      // initialize all needed values, extract them one by one in vectors, e.g. using max and min to determine which are heavy which light
      // get peptide id, that way determine donor, acceptor (alpha, beta)
      std::vector<std::string> peptides;
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
      vector< vector<std::string> > userParamNameLists;
      vector< vector<std::string> > userParamValueLists;
      vector< vector<std::string> > userParamUnitLists;

      for (std::multimap<std::string, int>::iterator it=range.first; it!=range.second; ++it)
      {
        xmlNode* cl_sii = siis[it->second];
        // Attributes
        std::string peptide = attribute_(cl_sii, "peptide_ref");
        peptides.push_back(peptide);
        double exp_mz = toDoubleOrNaN_(attribute_(cl_sii, "experimentalMassToCharge"));
        exp_mzs.push_back(exp_mz);

        // collect one RT entry per SII so RTs stays index-aligned with peptides/exp_mzs even when some SIIs lack MS:1000894;
        // a missing RT stays NaN (mirrors the non-XL path and PeptideIdentification::hasRT()) so it is not confused with a real 0 s RT
        double rt = std::numeric_limits<double>::quiet_NaN();
        bool has_rt = false;

        if (rank == 0)
        {
          rank = StringUtils::toInt32(attribute_(cl_sii, "rank"));
        }
        if (charge == 0)
        {
          charge = StringUtils::toInt32(attribute_(cl_sii, "chargeState"));
        }

        // CVs
        std::vector<xmlNode*> sii_cvp = descendants_(cl_sii, "cvParam");
        const  Size cv_count = sii_cvp.size();
        for (Size i = 0; i < cv_count; ++i)
        {
          xmlNode* element_sii_cvp = sii_cvp[i];
          if (equal_(attribute_(element_sii_cvp, "accession"), "MS:1002681")) // OpenXQuest:combined score
          {
            score = toDoubleOrZero_(attribute_(element_sii_cvp, "value"));
          }
          else if (equal_(attribute_(element_sii_cvp, "accession"), "MS:1002682")) // OpenXQuest: xcorr common
          {
            xcorrx = toDoubleOrZero_(attribute_(element_sii_cvp, "value"));
          }
          else if (equal_(attribute_(element_sii_cvp, "accession"), "MS:1002683")) // OpenXQuest: xcorr xlink
          {
            xcorrc = toDoubleOrZero_(attribute_(element_sii_cvp, "value"));
          }
          else if (equal_(attribute_(element_sii_cvp, "accession"), "MS:1002684")) // OpenXQuest: match-odds
          {
            matchodds = toDoubleOrZero_(attribute_(element_sii_cvp, "value"));
          }
          else if (equal_(attribute_(element_sii_cvp, "accession"), "MS:1002685")) // OpenXQuest: intsum
          {
            intsum = toDoubleOrZero_(attribute_(element_sii_cvp, "value"));
          }
          else if (equal_(attribute_(element_sii_cvp, "accession"), "MS:1002686")) // OpenXQuest: wTIC
          {
            wTIC = toDoubleOrZero_(attribute_(element_sii_cvp, "value"));
          }
          else if (equal_(attribute_(element_sii_cvp, "accession"), "MS:1003024")) // OpenPepXL:score
          {
            score = toDoubleOrZero_(attribute_(element_sii_cvp, "value"));
          }
          else if (equal_(attribute_(element_sii_cvp, "accession"), "MS:1000894")) // retention time
          {
            rt = toDoubleOrNaN_(attribute_(element_sii_cvp, "value"));
            if (equal_(attribute_(element_sii_cvp, "unitAccession"), "UO:0000031")) // minutes -> seconds (mirror non-XL path)
            {
              rt *= 60.0;
            }
            has_rt = true;
          }
        }

        // userParams
        vector<std::string> userParamNames;
        vector<std::string> userParamValues;
        vector<std::string> userParamUnits;

        std::vector<xmlNode*> sii_up = descendants_(cl_sii, "userParam");
        const  Size up_count = sii_up.size();
        for (Size i = 0; i < up_count; ++i)
        {
          xmlNode* element_sii_up = sii_up[i];
          userParamNames.push_back(attribute_(element_sii_up, "name"));
          userParamValues.push_back(attribute_(element_sii_up, "value"));
          userParamUnits.push_back(attribute_(element_sii_up, "unitName"));
        }
        userParamNameLists.push_back(userParamNames);
        userParamValueLists.push_back(userParamValues);
        userParamUnitLists.push_back(userParamUnits);

        RTs.push_back(rt);
        if (!has_rt)
        {
          OPENMS_LOG_WARN << "No retention time found for cross-linked SpectrumIdentificationItem in '" << spectrumID << "'." << endl;
        }

        // Fragmentation, does not matter where to get them. Look for them as long as the vector is empty
        if (frag_annotations.empty())
        {
          std::vector<xmlNode*> frag_element_list = descendants_(cl_sii, "Fragmentation");

          if (frag_element_list.size() > 0)
          {
            xmlNode* frag_element = frag_element_list[0];
            std::vector<xmlNode*> ion_types = descendants_(frag_element, "IonType");
            const  Size ion_type_count = ion_types.size();
            for (Size i = 0; i < ion_type_count; ++i)
            {
              xmlNode* ion_type_element = ion_types[i];
              int ion_charge = StringUtils::toInt32(attribute_(ion_type_element, "charge"));
              vector<std::string> indices;
              vector<std::string> positions;
              vector<std::string> intensities;
              vector<std::string> chains;
              vector<std::string> categories;
              std::string frag_type;
              std::string loss;

              StringUtils::split(attribute_(ion_type_element, "index"), " ", indices);

              std::vector<xmlNode*> frag_arrays = descendants_(ion_type_element, "FragmentArray");
              const Size frag_array_count = frag_arrays.size();
              for (Size f = 0; f < frag_array_count; ++f)
              {
                xmlNode* frag_array_element = frag_arrays[f];
                if ( equal_(attribute_(frag_array_element, "measure_ref"), "Measure_mz"))
                {
                  StringUtils::split(attribute_(frag_array_element, "values"), " ", positions);
                }
                if ( equal_(attribute_(frag_array_element, "measure_ref"), "Measure_int"))
                {
                  StringUtils::split(attribute_(frag_array_element, "values"), " ", intensities);
                }
              }

              std::vector<xmlNode*> userParams = descendants_(ion_type_element, "userParam");
              const Size userParam_count = userParams.size();
              for (Size u = 0; u < userParam_count; ++u)
              {
                xmlNode* userParam_element = userParams[u];
                if ( equal_(attribute_(userParam_element, "name"), "cross-link_chain"))
                {
                  StringUtils::split(attribute_(userParam_element, "value"), " ", chains);
                }
                if ( equal_(attribute_(userParam_element, "name"), "cross-link_ioncategory"))
                {
                  StringUtils::split(attribute_(userParam_element, "value"), " ", categories);
                }
              }

              std::vector<xmlNode*> cvts = descendants_(ion_type_element, "cvParam");
              const Size cvt_count = cvts.size();
              for (Size cvt = 0; cvt < cvt_count; ++cvt)
              {
                xmlNode* cvt_element = cvts[cvt];

                // Standard ions
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001229")) // frag: a ion
                {
                  frag_type = "a";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001224")) // frag: b ion
                {
                  frag_type = "b";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001231")) // frag: c ion
                {
                  frag_type = "c";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001228")) // frag: x ion
                {
                  frag_type = "x";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001220")) // frag: y ion
                {
                  frag_type = "y";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001230")) // frag: z ion
                {
                  frag_type = "z";
                }

                // Ions with H2O losses
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001234")) // frag: a ion - H2O
                {
                  frag_type = "a";
                  loss = "-H2O";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001222")) // frag: b ion - H20
                {
                  frag_type = "b";
                  loss = "-H2O";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001515")) // frag: c ion - H20
                {
                  frag_type = "c";
                  loss = "-H2O";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001519")) // frag: x ion - H20
                {
                  frag_type = "x";
                  loss = "-H2O";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001223")) // frag: y ion - H20
                {
                  frag_type = "y";
                  loss = "-H2O";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001517")) // frag: z ion - H20
                {
                  frag_type = "z";
                  loss = "-H2O";
                }

                // Ions with NH3 losses
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001235")) // frag: a ion - NH3
                {
                  frag_type = "a";
                  loss = "-NH3";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001232")) // frag: b ion - NH3
                {
                  frag_type = "b";
                  loss = "-NH3";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001516")) // frag: c ion - NH3
                {
                  frag_type = "c";
                  loss = "-NH3";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001520")) // frag: x ion - NH3
                {
                  frag_type = "x";
                  loss = "-NH3";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001233")) // frag: y ion - NH3
                {
                  frag_type = "y";
                  loss = "-NH3";
                }
                if (equal_(attribute_(cvt_element, "accession"), "MS:1001518")) // frag: z ion - NH3
                {
                  frag_type = "z";
                  loss = "-NH3";
                }
              }

              for (Size s = 0; s < indices.size(); ++s)
              {
                std::string annotation= "[" + chains[s] + "|" + categories[s]  + "$" + frag_type + indices[s] + loss + "]";

                PeptideHit::PeakAnnotation frag_anno;
                frag_anno.charge = ion_charge;
                frag_anno.mz = positions[s].empty() ? 0.0 : StringUtils::toDouble(positions[s]);
                frag_anno.intensity = intensities[s].empty() ? 0.0 : StringUtils::toDouble(intensities[s]);
                frag_anno.annotation = annotation;
                frag_annotations.push_back(frag_anno);
              }
            }
          }
        }
      }

      // Generate and fill PeptideIdentification
      // RTs is now collected one entry per SII (see loop above), so it stays aligned with exp_mzs/peptides.

      // Guard against a SpectrumIdentificationResult without any SpectrumIdentificationItems: the
      // min/max_element and light[0]/heavy[0] accesses below assume at least one experimental m/z.
      if (exp_mzs.empty())
      {
        OPENMS_LOG_WARN << "SpectrumIdentificationResult '" << spectrumID << "' has no SpectrumIdentificationItems; skipping." << endl;
        return;
      }

      vector<std::string> spectrumIDs;
      StringUtils::split(spectrumID, ",", spectrumIDs);

      vector<Size> light;
      vector<Size> heavy;
      double MZ_light = *std::min_element(exp_mzs.begin(), exp_mzs.end());
      double MZ_heavy = *std::max_element(exp_mzs.begin(), exp_mzs.end());

      // Are the cross-links labeled? A labeled (light/heavy) cross-link is reported as a light and a heavy
      // SII whose experimental m/z differ AND that reference two distinct spectra (comma-separated spectrumID).
      // Requiring >1 spectrum avoids misclassifying a single spectrum whose SIIs merely happen to have
      // differing experimental m/z as a labeled pair.
      bool labeled = MZ_light != MZ_heavy && spectrumIDs.size() > 1;

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
          std::string donor_pep = xl_id_donor_map_.at(peptides[i]);      // map::at throws an out-of-range
          alpha.push_back(i);
        }
        catch (const std::out_of_range& /*oor*/)
        {
          beta.push_back(i);
        }
      }
      std::string xl_type = "mono-link";

      // If no crosslink donors were found fall back to parsing each SII as a regular peptide identification
      if (alpha.empty())
      {
        OPENMS_LOG_WARN << "No crosslink donor found for SIIs in '" << spectrumID
                        << "'; parsing as regular (non-XLMS) peptide identifications." << endl;

        xmlNode* parent = element_res->parent;
        std::string sil = attribute_(parent, "id");

        PeptideIdentification current_pep_id;
        current_pep_id.setRT(RTs.empty() ? 0.0 : RTs[0]);
        current_pep_id.setMZ(exp_mzs.empty() ? 0.0 : exp_mzs[0]);
        current_pep_id.setMetaValue(Constants::UserParam::SPECTRUM_REFERENCE, spectrumID);
        current_pep_id.setHigherScoreBetter(false);

        // Reuse the regular SII parser so per-SII score type/userParams/passThreshold/calcMZ handling and
        // pro_id_ protein-hit population behave exactly like the normal non-XL path.
        xmlNode* reg_sii = firstElement_(element_res);
        while (reg_sii)
        {
          if (equal_(name_(reg_sii), "SpectrumIdentificationItem"))
          {
            parseSpectrumIdentificationItemElement_(reg_sii, current_pep_id, sil);
          }
          reg_sii = nextElement_(reg_sii);
        }
        current_pep_id.setIdentifier(runOfList_(sil).getIdentifier());
        current_pep_id.sort();
        pep_id_->push_back(current_pep_id);
        return;
      }

      SignedSize alpha_pos = xl_donor_pos_map_.at(xl_id_donor_map_.at(peptides[alpha[0]]));

      if (alpha.size() == beta.size()) // if a beta exists at all, it must be cross-link. But they should also each have the same number of SIIs in the mzid
      {
        xl_type = "cross-link";
      }
      else
      {
        try
        {
          std::string donor_val = xl_id_donor_map_.at(peptides[alpha[0]]);      // map::at throws an out-of-range
          std::string acceptor_val = xl_id_acceptor_map_.at(peptides[alpha[0]]);
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
      ph_alpha.setRank(rank - 1); // rank is 1-based in mzid, OpenMS uses 0-based
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

      vector<std::string> userParamNames_alpha = userParamNameLists[alpha[0]];
      vector<std::string> userParamValues_alpha = userParamValueLists[alpha[0]];
      vector<std::string> userParamUnits_alpha = userParamUnitLists[alpha[0]];

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
      else if ( xl_mod_map_.contains(peptides[alpha[0]]) )
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
        ph_beta.setRank(rank - 1); // rank is 1-based in mzid, OpenMS uses 0-based
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

      std::vector<std::string> unique_peptides;
      unique_peptides.push_back(peptides[alpha[0]]);
      if (!beta.empty())
      {
        unique_peptides.push_back(peptides[beta[0]]);
      }

      for (Size pep = 0; pep < unique_peptides.size(); ++pep)
      {

        const std::string& peptide_ref = unique_peptides[pep];

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
        pair<multimap<std::string, std::string>::iterator, multimap<std::string, std::string>::iterator> pev_its;
        pev_its = p_pv_map_.equal_range(peptide_ref);
        for (multimap<std::string, std::string>::iterator pev_it = pev_its.first; pev_it != pev_its.second; ++pev_it)
        {
          bool idec = false;
          OpenMS::PeptideEvidence pev;
          if (pe_ev_map_.contains(pev_it->second))
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

          if (pv_db_map_.contains(pev_it->second))
          {
            std::string& dpv = pv_db_map_[pev_it->second];
            DBSequence& db = db_sq_map_[dpv];
            pev.setProteinAccession(db.accession);

            ProteinIdentification& run = runOfList_(spectrumIdentificationList_ref);
            if (run.findHit(db.accession) == run.getHits().end())
            { // butt ugly! TODO @ mths for ProteinInference
              run.insertHit(ProteinHit());
              run.getHits().back().setSequence(db.sequence);
              run.getHits().back().setAccession(db.accession);
              run.getHits().back().setMetaValue("isDecoy", idec ? "true" : "false");
            }
          }
          phs[pep].addPeptideEvidence(pev);
        }
      }
      current_pep_id.setHits(phs);
      // each cross-link spectrum match is a peptide identification of its own and needs the run
      current_pep_id.setIdentifier(runOfList_(spectrumIdentificationList_ref).getIdentifier());
      pep_id_->push_back(current_pep_id);
      pep_id_->back().sort();
    }

    void MzIdentMLDOMHandler::parseSpectrumIdentificationItemElement_(xmlNode* spectrumIdentificationItemElement, PeptideIdentification& spectrum_identification, std::string& spectrumIdentificationList_ref)
    {
      std::string id = attribute_(spectrumIdentificationItemElement, "id");
      std::string name = attribute_(spectrumIdentificationItemElement, "name");

      long double calculatedMassToCharge = toDoubleOrNaN_(attribute_(spectrumIdentificationItemElement, "calculatedMassToCharge"));
//      long double calculatedPI = attribute_(spectrumIdentificationItemElement, "calculatedPI").toDouble();
      int chargeState = 0;
      try
      {
        chargeState = StringUtils::toInt32(attribute_(spectrumIdentificationItemElement, "chargeState"));
      }
      catch (...)
      {
        OPENMS_LOG_WARN << "Found unreadable 'chargeState'." << endl;
      }
      long double experimentalMassToCharge = toDoubleOrNaN_(attribute_(spectrumIdentificationItemElement, "experimentalMassToCharge"));
      int rank = 0;
      try
      {
        rank = StringUtils::toInt32(attribute_(spectrumIdentificationItemElement, "rank"));
        if (rank == 0) 
        {
          // MzIdentML ranks are typically 1-based. For some special data (PMF) it can be 0. 
          // In that case we treat it as 1-based as well otherwise it will underflow if converted to the 0-based OpenMS ranks.
          rank = 1;
          OPENMS_LOG_DEBUG << "Found rank 0. Assuming 1-based rank." << std::endl;
        }
        
      }
      catch (...)
      {
        OPENMS_LOG_WARN << "Found unreadable PSM rank." << endl;
      }

      std::string peptide_ref = attribute_(spectrumIdentificationItemElement, "peptide_ref");
//      std::string sample_ref = attribute_(spectrumIdentificationItemElement, "sample_ref");
//      std::string massTable_ref = attribute_(spectrumIdentificationItemElement, "massTable_ref");

      const bool pass = boolean_(attribute_(spectrumIdentificationItemElement, "passThreshold"));
      // TODO @all: where to store passThreshold value? set after score type eval in pass_threshold

      long double score = 0;
      const auto& [param_cv, param_user] = parseParamGroup_(children_(spectrumIdentificationItemElement));
      // The first recognized score in accession order is used, except that a lower-is-better search engine
      // specific score (e.g. Comet:expectation value or X!Tandem:expect) yields to a PSM-level q-value, which is
      // then used even if other scores lie between the two in accession order.
      const bool has_psm_q_value = std::any_of(param_cv.getCVTerms().begin(), param_cv.getCVTerms().end(),
        [this](const auto& term)
        {
          return (q_score_child_terms_.contains(term.first) || term.first == "MS:1002354") && term.first != "MS:1002055";
        });
      bool scoretype = false;
      bool use_q_value = false;
      for (map<std::string, vector<OpenMS::CVTerm>>::const_iterator scoreit = param_cv.getCVTerms().begin(); scoreit != param_cv.getCVTerms().end(); ++scoreit)
      {
        if (q_score_child_terms_.contains(scoreit->first) || scoreit->first == "MS:1002354")
        {
          if (scoreit->first != "MS:1002055") // do not use peptide-level q-values for now
          {
            score = toDoubleOrZero_(scoreit->second.front().getValue().toString());
            spectrum_identification.setHigherScoreBetter(false);
            spectrum_identification.setScoreType("q-value"); //higherIsBetter = false
            scoretype = true;
            break;
          }
        }
        else if (use_q_value)
        {
          continue;
        }
        else if (scoreit->first != "MS:1001143" && // the parent term itself has no numeric value; handled in the special case below
                 specific_score_child_terms_.contains(scoreit->first))
        {
          const bool lower_better = isLowerBetterEngineScore_(cv_.getTerm(scoreit->first));
          if (lower_better && has_psm_q_value)
          {
            use_q_value = true;
            continue;
          }
          score = toDoubleOrZero_(scoreit->second.front().getValue().toString());
          spectrum_identification.setHigherScoreBetter(!lower_better);
          spectrum_identification.setScoreType(scoreit->second.front().getName());
          scoretype = true;
          break;
        }
        else if (e_score_child_terms_.contains(scoreit->first))
        {
          score = toDoubleOrZero_(scoreit->second.front().getValue().toString());
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
        PeptideHit hit(score, rank - 1, chargeState, pep_map_[peptide_ref]); // rank in OpenMS is 0-based, in mzIdentML it is 1-based
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
        for (map<std::string, DataValue>::const_iterator up = param_user.begin(); up != param_user.end(); ++up)
        {
          hit.setMetaValue(up->first, up->second);
        }
        hit.setMetaValue("calcMZ", calculatedMassToCharge);
        spectrum_identification.setMZ(experimentalMassToCharge); // TODO @ mths for next PSI meeting: why is this not in SpectrumIdentificationResult in the schema? exp. m/z for one spec should not change from one id for it to the next!
        hit.setMetaValue("pass_threshold", pass); //TODO @ mths do not write metavalue pass_threshold

        //connect the PeptideHit with PeptideEvidences (for AABefore/After) and subsequently with DBSequence (for ProteinAccession)
        pair<multimap<std::string, std::string>::iterator, multimap<std::string, std::string>::iterator> pev_its;
        pev_its = p_pv_map_.equal_range(peptide_ref);
        for (multimap<std::string, std::string>::iterator pev_it = pev_its.first; pev_it != pev_its.second; ++pev_it)
        {
          bool idec = false;
          OpenMS::PeptideEvidence pev;
          if (pe_ev_map_.contains(pev_it->second))
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

          if (pv_db_map_.contains(pev_it->second))
          {
            std::string& dpv = pv_db_map_[pev_it->second];
            DBSequence& db = db_sq_map_[dpv];
            pev.setProteinAccession(db.accession);

            ProteinIdentification& run = runOfList_(spectrumIdentificationList_ref);
            if (run.findHit(db.accession) == run.getHits().end())
            { // butt ugly! TODO @ mths for ProteinInference
              run.insertHit(ProteinHit());
              run.getHits().back().setSequence(db.sequence);
              run.getHits().back().setAccession(db.accession);
              run.getHits().back().setMetaValue("isDecoy", idec ? "true" : "false");
            }
          }
          hit.addPeptideEvidence(pev);
        }
        //sort necessary for tests, as DOM children are obviously randomly ordered.
//        vector<std::string> temp = spectrum_identification.getHits().back().getPeptideEvidences().back().getProteinAccessions();
//        sort(temp.begin(), temp.end());
//        spectrum_identification.getHits().back().setProteinAccessions(temp);
        spectrum_identification.insertHit(hit);

        // due to redundant references this is not needed! peptideref already maps to all those PeptideEvidence elements
        //      xmlNode* child = firstElement_(spectrumIdentificationItemElement);
        //      while ( child )
        //      {
        //        if (equal_(name_(child), "PeptideEvidenceRef"))
        //        {
        //          ref = attribute_(element_si, "peptideEvidence_ref");
        //          //...
        //          spectrum_identification.getHits().back().setAABefore(char acid);
        //          spectrum_identification.getHits().back().setAAAfter (char acid);
        //          break;
        //        }
        //        child = nextElement_(child);
        //      }
      }
      // <Fragmentation> omitted for the time being
    }

    void MzIdentMLDOMHandler::parseProteinDetectionListElements_(const std::vector<xmlNode*>& proteinDetectionListElements)
    {
      const  Size node_count = proteinDetectionListElements.size();
      for (Size c = 0; c < node_count; ++c)
      {
        xmlNode* current_pr = proteinDetectionListElements[c];
        if (current_pr->type && // true is not NULL
            current_pr->type == XML_ELEMENT_NODE) // is element - possibly not necessary after getElementsByTagName
        {
          // Found element node: re-cast as element
          xmlNode* element_pr = current_pr;

//          std::string id = attribute_(element_pr, "id");
//          pair<CVTermList, map<std::string, DataValue> > params = parseParamGroup_(children_(current_pr));

          // TODO @mths : this needs to be a ProteinIdentification for the ProteinDetectionListElement which is not mandatory and used in downstream analysis ProteinInference etc.
//          pro_id_->push_back(ProteinIdentification());
//          pro_id_->back().setSearchEngine(search_engine_);
//          pro_id_->back().setSearchEngineVersion(search_engine_version_);
//          pro_id_->back().setIdentifier(search_engine_);

          //      SearchParameters  search_parameters_
          //      DateTime  date_
          //      std::string protein_score_type_ <- from ProteinDetectionProtocol
          //      DoubleReal    protein_significance_threshold_ <- from ProteinDetectionProtocol

          xmlNode* child = firstElement_(element_pr);
          while (child)
          {
            if (equal_(name_(child), "ProteinAmbiguityGroup"))
            {
              parseProteinAmbiguityGroupElement_(child, pro_id_->back());
            }
            child = nextElement_(child);
          }
        }
      }
    }

    void MzIdentMLDOMHandler::parseProteinAmbiguityGroupElement_(xmlNode* proteinAmbiguityGroupElement, ProteinIdentification& protein_identification)
    {
//      std::string id = attribute_(proteinAmbiguityGroupElement, "id");
//      pair<CVTermList, map<std::string, DataValue> > params = parseParamGroup_(children_(proteinAmbiguityGroupElement));

      //fill pro_id_->back() with content,
      xmlNode* child = firstElement_(proteinAmbiguityGroupElement);
      while (child)
      {
        if (equal_(name_(child), "ProteinDetectionHypothesis"))
        {
          parseProteinDetectionHypothesisElement_(child, protein_identification);
        }
        child = nextElement_(child);
      }
    }

    void MzIdentMLDOMHandler::parseProteinDetectionHypothesisElement_(xmlNode* proteinDetectionHypothesisElement, ProteinIdentification& protein_identification)
    {
      std::string dBSequence_ref = attribute_(proteinDetectionHypothesisElement, "dBSequence_ref");

//      pair<CVTermList, map<std::string, DataValue> > params = parseParamGroup_(children_(proteinDetectionHypothesisElement));

      DBSequence& db = db_sq_map_[dBSequence_ref];

      protein_identification.insertHit(ProteinHit());
      protein_identification.getHits().back().setSequence(db.sequence);
      protein_identification.getHits().back().setAccession(db.accession);
//      protein_identification.getHits().back().setCoverage((long double)params.first.getCVTerms()["MS:1001093"].front().getValue()); //TODO @ mths: calc percent
//      protein_identification.getHits().back().setScore(boost::lexical_cast<double>(params.first.getCVTerms()["MS:1001171"].front().getValue())); //or any other score
    }

    SignedSize MzIdentMLDOMHandler::inferModificationLocation_(xmlNode* modification_element, Size peptide_length) const
    {
      // The 'location' attribute of a <Modification> is optional in mzIdentML (schema use="optional").
      // Several search engines omit it for terminal modifications. As a best-effort recovery, look up the
      // referenced modification(s) in the ModificationsDB and derive the position from their terminal
      // specificity: 0 for the N-terminus or peptide_length + 1 for the C-terminus. Returns -1 if no
      // terminal position can be inferred.
      //
      // The N-terminus is preferred because N-terminal modifications (e.g. acetylation, the case reported
      // in #5443) are by far the most common ones for which producers drop the 'location'. The C-terminus
      // is only assigned for modifications that are *exclusively* C-terminal: a modification that may also
      // sit on an internal residue (e.g. Oxidation, which has a rare 'Any C-term' variant next to its common
      // 'Anywhere' ones) has a genuinely unknown position when 'location' is missing and is left unassigned.

      // A location-less modification that names a concrete residue (mzIdentML 'residues' attribute, e.g.
      // "K") is residue-specific, not terminal: even knowing the residue type, its exact position is still
      // unknown without 'location' (the residue may occur several times), so we must not force it onto a
      // terminus. Only the terminus marker "." (or an absent/empty value) is compatible with terminal
      // placement. This matches real data: in mzIdentML from Mascot (and search engines generally) internal
      // modifications carry residues="<AA>", while the terminal modifications for which 'location' is
      // sometimes dropped carry residues="." or no 'residues' attribute at all.
      const std::string residues = attribute_(modification_element, "residues");
      if (!residues.empty() && residues != ".")
      {
        return -1;
      }

      const ModificationsDB* mod_db = ModificationsDB::getInstance();
      for (xmlNode* cvp = firstElement_(modification_element); cvp != nullptr; cvp = nextElement_(cvp))
      {
        const std::string mod_name = attribute_(cvp, "name");
        // guard with has() first: searchModifications() would otherwise log a "not found" warning for
        // non-modification cvParams (e.g. neutral loss terms) that legitimately appear as children here.
        if (mod_name.empty() || !mod_db->has(mod_name))
        {
          continue;
        }

        std::set<const ResidueModification*> mods;
        mod_db->searchModifications(mods, mod_name);

        bool can_be_n_terminal = false; // has an N_TERM / PROTEIN_N_TERM variant
        bool can_be_c_terminal = false; // has a  C_TERM / PROTEIN_C_TERM variant
        bool can_be_internal = false;   // has an ANYWHERE variant (i.e. may sit on an internal residue)
        for (const ResidueModification* mod : mods)
        {
          switch (mod->getTermSpecificity())
          {
            case ResidueModification::N_TERM:
            case ResidueModification::PROTEIN_N_TERM:
              can_be_n_terminal = true;
              break;
            case ResidueModification::C_TERM:
            case ResidueModification::PROTEIN_C_TERM:
              can_be_c_terminal = true;
              break;
            case ResidueModification::ANYWHERE:
              can_be_internal = true;
              break;
            default:
              break;
          }
        }

        if (can_be_n_terminal)
        {
          return 0; // N-terminus
        }
        if (can_be_c_terminal && !can_be_internal)
        {
          return static_cast<SignedSize>(peptide_length) + 1; // C-terminus (only if not also internal)
        }
      }
      return -1; // position could not be inferred
    }

    AASequence MzIdentMLDOMHandler::parsePeptideSiblings_(xmlNode* peptide)
    {
      std::vector<xmlNode*> peptideSiblings = children_(peptide);
      const  Size node_count = peptideSiblings.size();
      std::string as;
      //1. Sequence
      for (Size c = 0; c < node_count; ++c)
      {
        xmlNode* current_sib = peptideSiblings[c];
        if (current_sib->type && // true is not NULL
            current_sib->type == XML_ELEMENT_NODE)
        {
          xmlNode* element_sib = current_sib;
          if (equal_(name_(element_sib), "PeptideSequence"))
          {
            // The text content of the element joins the text of all its children (CDATA sections included,
            // comments excluded) and is an empty string when the element has no children, as in the
            // '<PeptideSequence/>' the writer emits for a hit with an empty sequence. Inspecting the first
            // child node instead dereferenced a null pointer for such elements and discarded sequences
            // that start with a comment or are held in a CDATA section.
            as = content_(element_sib);
          }
        }
      }
      // Trim once, before the substitution positions are applied, so that they count residues only.
      StringUtils::trim(as);
      if (as.empty())
      {
        // The caller catches this, logs the problem and stores an empty sequence for the peptide.
        throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    attribute_(peptide, "id"),
                                    "Peptide has no amino acid sequence.");
      }
      //2. Substitutions
      for (Size c = 0; c < node_count; ++c)
      {
        xmlNode* current_sib = peptideSiblings[c];
        if (current_sib->type && // true is not NULL
            current_sib->type == XML_ELEMENT_NODE)
        {
          xmlNode* element_sib = current_sib;
          if (equal_(name_(element_sib), "SubstitutionModification"))
          {

            std::string location = attribute_(element_sib, "location");
            char originalResidue = attribute_(element_sib, "originalResidue")[0];
            std::string replacementResidueStr = attribute_(element_sib, "replacementResidue");
            if (replacementResidueStr.empty())
            {
              throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          attribute_(peptide, "id"),
                                          "SubstitutionModification is missing the replacementResidue attribute.");
            }
            char replacementResidue = replacementResidueStr[0];

            if (!location.empty())
            {
              // location is a 1-based residue position. A location of 0 (or any
              // value below 1) makes "toInt32(location) - 1" negative, which wraps
              // to a huge size_type and writes one byte before as's character
              // buffer; a location past as.size() writes beyond it. Both are
              // undefined behaviour on std::string::operator[] (CPP-169).
              SignedSize pos = StringUtils::toInt32(location);
              if (pos < 1 || static_cast<Size>(pos) > as.size())
              {
                throw Exception::ParseError(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                            attribute_(peptide, "id"),
                                            "SubstitutionModification location '" + location + "' is out of range for a " + StringUtils::toStr(as.size()) + "-residue sequence.");
              }
              as[pos - 1] = replacementResidue;
            }
            else if (StringUtils::hasSubstring(as, originalResidue)) //no location - every occurrence will be replaced
            {
              StringUtils::substitute(as, originalResidue, replacementResidue);
            }
            else
            {
              throw std::runtime_error("ERROR : Non Text Node");
            }
          }
        }
      }
      //3. Modifications
      AASequence aas = AASequence::fromString(as);
      for (Size c = 0; c < node_count; ++c)
      {
        xmlNode* current_sib = peptideSiblings[c];
        if (current_sib->type && // true is not NULL
            current_sib->type == XML_ELEMENT_NODE)
        {
          xmlNode* element_sib = current_sib;
          if (equal_(name_(element_sib), "Modification"))
          {
            // 'location' is optional in mzIdentML: the 1-based position of the modified residue,
            // with 0 = N-terminus and (peptide length + 1) = C-terminus.
            const std::string location_str = attribute_(element_sib, "location");
            SignedSize index = -2; // sentinel: no valid location determined (yet)
            if (!location_str.empty())
            {
              try
              {
                index = static_cast<SignedSize>(StringUtils::toInt32(location_str));
              }
              catch (...)
              {
                OPENMS_LOG_WARN << "Found unreadable modification location '" << location_str << "'." << endl;
              }
            }

            // Some search engines omit 'location' for terminal modifications (e.g. N-terminal acetylation).
            // Previously the missing/invalid location left 'index' at a negative sentinel, which was then used
            // as an (invalid) residue position and crashed the reader with an IndexOverflow exception (#5443).
            // If the location is missing, unreadable, or out of range, try to infer a terminal position from
            // the modification's terminal specificity; otherwise skip the modification instead of crashing.
            if (index < 0 || index > static_cast<SignedSize>(aas.size()) + 1)
            {
              const SignedSize inferred_index = inferModificationLocation_(element_sib, aas.size());
              if (inferred_index < 0)
              {
                OPENMS_LOG_WARN << "Skipping modification with missing or invalid 'location' attribute; its position could not be inferred." << endl;
                continue;
              }
              index = inferred_index;
            }

            //double monoisotopicMassDelta = attribute_(element_dbs, "monoisotopicMassDelta");

            if (xl_ms_search_) // special case: XL-MS search results
            {
              std::string pep_id = attribute_(peptide, "id");
              //std::vector<xmlNode*> cvParams = descendants_(element_sib, "cvParam");
              xmlNode* cvp = firstElement_(element_sib);
              //for (Size i = 0; i < cvParams.length(); ++i)
              bool donor_acceptor_found = false;
              bool xlink_mod_found = false;

              while (cvp)
              {
                if (equal_(attribute_(cvp, "accession"), "MS:1002509")) // crosslink donor
                {
                  std::string donor_val = attribute_(cvp, "value");
                  xl_id_donor_map_.insert(make_pair(pep_id, donor_val));
                  std::string massdelta = attribute_(element_sib, "monoisotopicMassDelta");
                  double monoisotopicMassDelta = massdelta.empty() ? 0.0 : StringUtils::toDouble(massdelta);
                  xl_mass_map_.insert(make_pair(pep_id, monoisotopicMassDelta));
                  xl_donor_pos_map_.insert(make_pair(donor_val, index-1));

                  // Find the actual modification name from UNIMOD/XLMOD cvParam
                  std::string xl_mod_name;
                  xmlNode* cvp1 = firstElement_(element_sib);
                  while (cvp1)
                  {
                    std::string cvref = attribute_(cvp1, "cvRef");
                    if (cvref == "UNIMOD" || cvref == "XLMOD")
                    {
                      xl_mod_name = attribute_(cvp1, "name");
                      break;
                    }
                    cvp1 = nextElement_(cvp1);
                  }
                  if (xl_mod_name.empty())
                  {
                    // fallback
                    xl_mod_name = attribute_(firstElement_(element_sib), "name");
                  }
                  xl_mod_map_.insert(make_pair(pep_id, xl_mod_name));
                  donor_acceptor_found = true;
                }
                else if (equal_(attribute_(cvp, "accession"), "MS:1002510")) // crosslink acceptor
                {
                  std::string acceptor_val = attribute_(cvp, "value");
                  xl_id_acceptor_map_.insert(make_pair(pep_id, acceptor_val));
                  xl_acceptor_pos_map_.insert(make_pair(acceptor_val, index-1));
                  donor_acceptor_found = true;
                }
                else
                {
                  CVTerm cv = parseCvParam_(cvp);
                  const std::string& cvname = cv.getName();
                  if (StringUtils::hasPrefix(cvname, "Xlink") || StringUtils::hasPrefix(cv.getAccession(), "XLMOD"))
                  {
                    xlink_mod_found = true;
                  }
                  if (StringUtils::hasSubstring(cvname, "unknown mono-link"))
                  {
                    xl_mod_map_.insert(make_pair(pep_id, cvname));
                  }
                  else // normal mod, copied from below
                  {
                    if ( (cv.getCVIdentifierRef() != "UNIMOD") && (cv.getCVIdentifierRef() != "XLMOD") )
                    {
                      // e.g.  <cvParam accession="MS:1001524" name="fragment neutral loss" cvRef="PSI-MS" value="0" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO"/>
                      cvp = nextElement_(cvp);
                      continue;
                    }
                    if (index == 0)
                    {
                      try // does not work for cross-links yet, but the information is finally stored as MetaValues of the PeptideHit
                      {
                        if (cvname == "unknown modification")
                        {
                          const std::string & cvvalue = cv.getValue();
                          if (ModificationsDB::getInstance()->has(cvvalue) && !cvvalue.empty())
                          {
                            aas.setNTerminalModification(cv.getValue());
                          }
                        }
                        else
                        {
                          aas.setNTerminalModification(cvname);
                        }
                        cvp = nextElement_(cvp);
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
                        const std::string & cvvalue = cv.getValue();
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
                      cvp = nextElement_(cvp);
                      continue;
                    }
                    else
                    {
                      if (cvname == "unknown modification")
                      {
                        const std::string & cvvalue = cv.getValue();
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
                      cvp = nextElement_(cvp);
                      continue;
/* TODO enetz: look up in XLDB and remove this ha
                        // this is a bad hack to avoid a long list of warnings in the case of XL-MS data
                        if ( !(StringUtils::toStr(e.what()).hasSubstring("'DSG'") || StringUtils::toStr(e.what()).hasSubstring("'DSS'") || StringUtils::toStr(e.what()).hasSubstring("'EDC'")) || StringUtils::toStr(e.what()).hasSubstring("'BS3'") || StringUtils::toStr(e.what()).hasSubstring("'BS2G'") )
                        {
                          OPENMS_LOG_WARN << e.getName() << ": " << e.what() << " Sequence: " << aas.toUnmodifiedString() << ", residue " << aas.getResidue(index - 1).getName() << "@" << StringUtils::toStr(index) << endl;
                        }
*/
                    }
                  }
                }
                cvp = nextElement_(cvp);
              }
              if ( (!donor_acceptor_found) && (xlink_mod_found) ) // mono-link, here using pep_id also as the CV value, since mono-links don't have a cross-linking CV term
              {
                xl_id_donor_map_.insert(make_pair(pep_id, pep_id));
                xl_donor_pos_map_.insert(make_pair(pep_id, index-1));
              }
            }
            else // general case: no XL-MS result
            {
              xmlNode* cvp = firstElement_(element_sib);
              while (cvp)
              {
                CVTerm cv = parseCvParam_(cvp);
                if (cv.getAccession() == "MS:1001460") // unknown modification
                {
                  const std::string & cvvalue = cv.getValue();
                  if (cv.hasValue() && ModificationsDB::getInstance()->has(cvvalue) && !cvvalue.empty())  // why do we need to check for empty?
                  {
                    // Case 1: unknown (to e.g., third-party tool) modification known to OpenMS (see value)
                    //  <Modification location="0" monoisotopicMassDelta="17.031558">
                    //  <cvParam cvRef="PSI-MS" accession="MS:1001460" name="unknown modification" value="Methyl:2H(2)13C"/>
                    const std::string & mname = cvvalue;
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
                    cvp = nextElement_(cvp);
                    continue;
                  }
                  else
                  {
                    // Case 2: unknown modification (needs to be added to ModificationsDB)
                    // note, this is optional
                    double mass_delta = 0;
                    bool has_mass_delta = false;
                    std::string mod;

                    // try to parse information, give up if we cannot
                    try
                    {
                      mod = attribute_(element_sib, "monoisotopicMassDelta");
                      mass_delta = static_cast<double>(StringUtils::toDouble(mod));
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
                    const ModificationsDB* mod_db = ModificationsDB::getInstance();
                    if (index == 0)
                    {
                      // n-terminal
                      std::string residue_name = ".[+" + mod + "]";
                      std::string residue_id = ".n[" + mod + "]";

                      // Check if it already exists, if not create new modification, transfer
                      // ownership to ModDB
                      if (!mod_db->has(residue_id))
                      {
                        unique_ptr<ResidueModification> new_mod(new ResidueModification);
                        new_mod->setFullId(residue_id); // setting FullId but not Id makes it a user-defined mod
                        new_mod->setProvenance(ResidueModification::MASS_ONLY);
                        new_mod->setFullName(residue_name); // display name
                        new_mod->setDiffMonoMass(mass_delta);
                        new_mod->setMonoMass(mass_delta + Residue::getInternalToNTerm().getMonoWeight());
                        new_mod->setTermSpecificity(ResidueModification::N_TERM);
                        mod_db->addModification(std::move(new_mod));
                      }

                      aas.setNTerminalModification(residue_id);
                      cvp = nextElement_(cvp);
                      continue;
                    }
                    else if (index == (int)aas.size() +1)
                    {
                      // c-terminal
                      std::string residue_name = ".[" + mod + "]";
                      std::string residue_id = ".c[" + mod + "]";

                      // Check if it already exists, if not create new modification, transfer
                      // ownership to ModDB (check residue_id, which is what gets registered)
                      if (!mod_db->has(residue_id))
                      {
                        unique_ptr<ResidueModification> new_mod(new ResidueModification);
                        new_mod->setFullId(residue_id); // setting FullId but not Id makes it a user-defined mod
                        new_mod->setProvenance(ResidueModification::MASS_ONLY);
                        new_mod->setFullName(residue_name); // display name
                        new_mod->setDiffMonoMass(mass_delta);
                        new_mod->setMonoMass(mass_delta + Residue::getInternalToCTerm().getMonoWeight());
                        new_mod->setTermSpecificity(ResidueModification::C_TERM);
                        mod_db->addModification(std::move(new_mod));
                      }
                      aas.setCTerminalModification(residue_id);
                      cvp = nextElement_(cvp);
                      continue;
                    }
                    else if (index > 0 && index <= (int)aas.size() )
                    {
                      // internal modification
                      const Residue& residue = aas[index-1];
                      std::string residue_name = residue.getOneLetterCode() + "[" + mod + "]"; // e.g. N[12345.6]
                      std::string modification_name = "[" + mod + "]";

                      if (!mod_db->has(residue_name))
                      {
                        // create new modification
                        unique_ptr<ResidueModification> new_mod(new ResidueModification);
                        new_mod->setFullId(residue_name); // setting FullId but not Id makes it a user-defined mod
                        new_mod->setProvenance(ResidueModification::MASS_ONLY);
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
                      cvp = nextElement_(cvp);
                      continue;
                    }
                  }
                }
                if (cv.getCVIdentifierRef() != "UNIMOD")
                {
                  // e.g.  <cvParam accession="MS:1001524" name="fragment neutral loss" cvRef="PSI-MS" value="0" unitAccession="UO:0000221" unitName="dalton" unitCvRef="UO"/>
                  cvp = nextElement_(cvp);
                  continue;
                }

                if (index == 0)
                {
                  if (cv.getName() == "unknown modification")
                  {
                    aas.setNTerminalModification(cv.getValue());
                    cvp = nextElement_(cvp);
                    continue;
                  }
                  else
                  {
                    aas.setNTerminalModification(cv.getName());
                    cvp = nextElement_(cvp);
                    continue;
                  }
                }
                else if (index == static_cast<SignedSize>(aas.size() + 1))
                {
                  aas.setCTerminalModification(cv.getName());
                  cvp = nextElement_(cvp);
                  continue;
                }
                else
                {
                  try
                  {
                     aas.setModification(index - 1, cv.getName()); //TODO @mths,Timo : do this via UNIMOD accessions
                     cvp = nextElement_(cvp);
                     continue;
                  }
                  catch (Exception::BaseException& e)
                  {
                    OPENMS_LOG_WARN << e.getName() << ": " << e.what() << " Sequence: " << aas.toUnmodifiedString() << ", residue " << aas.getResidue(index - 1).getName() << "@" << StringUtils::toStr(index) << "\n";
                  }
                }

                cvp = nextElement_(cvp);
              }
            }
          }
        }
      }
      return aas;
    }

    void MzIdentMLDOMHandler::buildCvList_(xmlNode* cvElements)
    {
      xmlNode* cv1 = element_(cvElements->doc, "cv");
      setAttribute_(cv1, "id", "PSI-MS");
      setAttribute_(cv1, "fullName",
                        "Proteomics Standards Initiative Mass Spectrometry Vocabularies");
      setAttribute_(cv1, "uri",
                        "http://psidev.cvs.sourceforge.net/viewvc/*checkout*/psidev/psi/psi-ms/mzML/controlledVocabulary/psi-ms.obo");
      setAttribute_(cv1, "version", "2.32.0");
      append_(cvElements, cv1);
      xmlNode* cv2 = element_(cvElements->doc, "cv");
      setAttribute_(cv2, "id", "UNIMOD");
      setAttribute_(cv2, "fullName",
                        "UNIMOD");
      setAttribute_(cv2, "uri",
                        "http://www.unimod.org/obo/unimod.obo");
      append_(cvElements, cv2);
      xmlNode* cv3 = element_(cvElements->doc, "cv");
      setAttribute_(cv3, "id", "UO");
      setAttribute_(cv3, "fullName",
                        "UNIT-ONTOLOGY");
      setAttribute_(cv3, "uri",
                        "http://obo.cvs.sourceforge.net/*checkout*/obo/obo/ontology/phenotype/unit.obo");
      append_(cvElements, cv3);
    }

    void MzIdentMLDOMHandler::buildAnalysisSoftwareList_(xmlNode* analysisSoftwareElements)
    {
      const std::string engine = cpro_id_->empty() || cpro_id_->front().getSearchEngine().empty()
        ? "TOPP software" : cpro_id_->front().getSearchEngine();
      const std::string version = cpro_id_->empty() || cpro_id_->front().getSearchEngineVersion().empty()
        ? "unknown" : cpro_id_->front().getSearchEngineVersion();
      xmlNode* current_as = element_(analysisSoftwareElements->doc, "AnalysisSoftware");
      setAttribute_(current_as, "id", "OpenMS" + StringUtils::toStr(UniqueIdGenerator::getUniqueId()));
      setAttribute_(current_as, "version", version);
      setAttribute_(current_as, "name", engine);
      append_(analysisSoftwareElements, current_as);
      xmlNode* current_sw = element_(current_as->doc, "SoftwareName");
      xmlNode* parameter = element_(current_as->doc, cv_.hasTermWithName(engine) ? "cvParam" : "userParam");
      setAttribute_(parameter, "name", engine);
      if (cv_.hasTermWithName(engine))
      {
        setAttribute_(parameter, "cvRef", "PSI-MS");
        setAttribute_(parameter, "accession", cv_.getTermByName(engine).id);
      }
      append_(current_sw, parameter);
      append_(current_as, current_sw);
    }

    void MzIdentMLDOMHandler::buildSequenceCollection_(xmlNode* sequenceCollectionElements)
    {
      for (map<std::string, DBSequence>::iterator dbs = db_sq_map_.begin(); dbs != db_sq_map_.end(); ++dbs)
      {
        xmlNode* current_dbs = element_(sequenceCollectionElements->doc, "DBSequence");
        setAttribute_(current_dbs, "id", dbs->second.accession);
        setAttribute_(current_dbs, "length", StringUtils::toStr(dbs->second.sequence.length()));
        setAttribute_(current_dbs, "accession", dbs->second.accession);
        setAttribute_(current_dbs, "searchDatabase_ref", dbs->second.database_ref); // This is going to be wrong
        xmlNode* current_seq = element_(current_dbs->doc, "Seq");
        setText_(current_seq, dbs->second.sequence);
        append_(current_dbs, current_seq);
        append_(sequenceCollectionElements, current_dbs);
      }

      for (map<std::string, AASequence>::iterator peps = pep_map_.begin(); peps != pep_map_.end(); ++peps)
      {
        xmlNode* current_pep = element_(sequenceCollectionElements->doc, "Peptide");
        setAttribute_(current_pep, "id", peps->first);
        xmlNode* current_seq = element_(current_pep->doc, "PeptideSequence");
        setText_(current_seq, peps->second.toUnmodifiedString());
        append_(current_pep, current_seq);
        if (peps->second.hasNTerminalModification())
        {
          const ResidueModification* mod = peps->second.getNTerminalModification();
          xmlNode* current_mod = element_(current_pep->doc, "Modification");
          xmlNode* current_cv = element_(current_pep->doc, "cvParam");
          setAttribute_(current_mod, "location", "0");
          setAttribute_(current_mod, "monoisotopicMassDelta", StringUtils::toStr(mod->getDiffMonoMass()));
          std::string origin(1, mod->getOrigin());
          if (origin == "X") origin = ".";
          setAttribute_(current_mod, "residues", origin);
          setAttribute_(current_cv, "name", mod->getName());
          setAttribute_(current_cv, "cvRef", "UNIMOD");
          setAttribute_(current_cv, "accession", mod->getUniModAccession());

          append_(current_mod, current_cv);
          append_(current_pep, current_mod);
        }
        if (peps->second.hasCTerminalModification())
        {
          const ResidueModification* mod = peps->second.getCTerminalModification();
          xmlNode* current_mod = element_(current_pep->doc, "Modification");
          xmlNode* current_cv = element_(current_mod->doc, "cvParam");
          setAttribute_(current_mod, "location", StringUtils::toStr(peps->second.size() + 1));
          setAttribute_(current_mod, "monoisotopicMassDelta", StringUtils::toStr(mod->getDiffMonoMass()));
          std::string origin(1, mod->getOrigin());
          if (origin == "X") origin = ".";
          setAttribute_(current_mod, "residues", origin);
          setAttribute_(current_cv, "name", mod->getName());
          setAttribute_(current_cv, "cvRef", "UNIMOD");
          setAttribute_(current_cv, "accession", mod->getUniModAccession());

          append_(current_mod, current_cv);
          append_(current_pep, current_mod);
        }
        if (peps->second.isModified())
        {
          Size i = 0;
          for (AASequence::ConstIterator res = peps->second.begin(); res != peps->second.end(); ++res, ++i)
          {
            const ResidueModification* mod = res->getModification();
            if (mod == nullptr) continue;
            xmlNode* current_mod = element_(current_pep->doc, "Modification");
            xmlNode* current_cv = element_(current_pep->doc, "cvParam");
            setAttribute_(current_mod, "location", StringUtils::toStr(i));
            setAttribute_(current_mod, "monoisotopicMassDelta", StringUtils::toStr(mod->getDiffMonoMass()));
            setAttribute_(current_mod, "residues", StringUtils::toStr(mod->getOrigin()));
            setAttribute_(current_cv, "name", mod->getName());
            setAttribute_(current_cv, "cvRef", "UNIMOD");
            setAttribute_(current_cv, "accession", mod->getUniModAccession());

            append_(current_mod, current_cv);
            append_(current_pep, current_mod);
          }
        }
        append_(sequenceCollectionElements, current_pep);
      }

      for (map<std::string, PeptideEvidence>::iterator pevs = pe_ev_map_.begin(); pevs != pe_ev_map_.end(); ++pevs)
      {
        xmlNode* current_pev = element_(sequenceCollectionElements->doc, "PeptideEvidence");
        setAttribute_(current_pev, "peptide_ref", "TBA");
        setAttribute_(current_pev, "id", pevs->first);
        setAttribute_(current_pev, "start", StringUtils::toStr(pevs->second.start));
        setAttribute_(current_pev, "end", StringUtils::toStr(pevs->second.stop));
        setAttribute_(current_pev, "pre", StringUtils::toStr(pevs->second.pre));
        setAttribute_(current_pev, "post", StringUtils::toStr(pevs->second.post));
        // do not forget to annotate the decoy
        setAttribute_(current_pev, "isDecoy", "false");
        append_(sequenceCollectionElements, current_pev);
      }
    }

    void MzIdentMLDOMHandler::buildAnalysisCollection_(xmlNode* analysisCollectionElements)
    {
      // for now there is only one search per file
      xmlNode* current_si = element_(analysisCollectionElements->doc, "SpectrumIdentification");
      setAttribute_(current_si, "id", "TBA");
      setAttribute_(current_si, "spectrumIdentificationProtocol_ref", "SIP");
      setAttribute_(current_si, "spectrumIdentificationList_ref", "SIL");
      setAttribute_(current_si, "activityDate", "now");
      xmlNode* current_is = element_(current_si->doc, "InputSpectra");
      setAttribute_(current_is, "spectraData_ref", "TODO"); // TODO @ mths while DataCollection
      xmlNode* current_sr = element_(current_si->doc, "SearchDatabaseRef");
      setAttribute_(current_sr, "searchDatabase_ref", "TODO"); // TODO @ mths while DataCollection
      append_(current_si, current_is);
      append_(current_si, current_sr);
      // and no ProteinDetection for now
      append_(analysisCollectionElements, current_si);
    }

    void MzIdentMLDOMHandler::buildAnalysisProtocolCollection_(xmlNode* protocolElements)
    {
      // for now there is only one search per file
      xmlNode* current_sp = element_(protocolElements->doc, "SpectrumIdentificationProtocol");
      setAttribute_(current_sp, "id", "SIP");
      setAttribute_(current_sp, "analysisSoftware_ref", "what now?");
      append_(protocolElements, current_sp);
      xmlNode* current_st = element_(current_sp->doc, "SearchType");
      append_(current_sp, current_st);
      xmlNode* current_cv = element_(current_st->doc, "cvParam");
      setAttribute_(current_cv, "accession", "MS:1001083"); // TODO @ mths for now static cv
      setAttribute_(current_cv, "name", "ms-ms search");
      setAttribute_(current_cv, "cvRef", "PSI-MS");
      append_(current_st, current_cv);

      //for now no <AdditionalSearchParams>, <ModificationParams>, <MassTable id="MT" msLevel="1 2">, <FragmentTolerance>, <ParentTolerance>, <DatabaseFilters>, <DatabaseTranslations>

      xmlNode* current_th = element_(current_sp->doc, "Threshold");
      xmlNode* current_up = element_(current_th->doc, "userParam");
      setAttribute_(current_up, "value", "0.05"); // TODO @ mths for now static cv
      setAttribute_(current_up, "name", "some significance threshold");
      append_(current_st, current_up);
      // and no ProteinDetection for now
      append_(protocolElements, current_th);
    }

    void MzIdentMLDOMHandler::buildInputDataCollection_(xmlNode* inputElements)
    {
      xmlNode* current_sf = element_(inputElements->doc, "SourceFile");
      setAttribute_(current_sf, "location", "file:///tmp/test.dat");
      setAttribute_(current_sf, "id", "SF1");
      buildEnclosedCV_(current_sf, "FileFormat", "MS:1001199", "Mascot DAT file", "PSI-MS"); // TODO @ mths for now static cv
      append_(inputElements, current_sf);

      xmlNode* current_sd = element_(inputElements->doc, "SearchDatabase");
      setAttribute_(current_sd, "location", "file:///tmp/test.fasta");
      setAttribute_(current_sd, "id", "DB1");
      setAttribute_(current_sd, "name", "SwissProt");
      setAttribute_(current_sd, "numDatabaseSequences", "257964");
      setAttribute_(current_sd, "numResidues", "93947433");
      setAttribute_(current_sd, "releaseDate", "2011-03-01T21:32:52");
      setAttribute_(current_sd, "version", "SwissProt_51.6.fasta");
      buildEnclosedCV_(current_sd, "FileFormat", "MS:1001348", "FASTA format", "PSI-MS"); // TODO @ mths for now static cv
      xmlNode* current_dn = element_(current_sd->doc, "DatabaseName");
      xmlNode* current_up = element_(current_dn->doc, "userParam");
      setAttribute_(current_up, "name", "SwissProt_51.6.fasta"); // TODO @ mths for now static cv
      append_(current_dn, current_up);
      append_(current_sd, current_dn);

      xmlNode* current_cv = element_(current_sd->doc, "cvParam");
      setAttribute_(current_cv, "accession", "MS:1001073"); // TODO @ mths for now static cv
      setAttribute_(current_cv, "name", "database type amino acid");
      setAttribute_(current_cv, "cvRef", "PSI-MS");
      append_(current_sd, current_cv);
      append_(inputElements, current_sd);

      xmlNode* current_spd = element_(inputElements->doc, "SpectraData");
      setAttribute_(current_spd, "location", "file:///tmp/test.mzML");
      setAttribute_(current_spd, "id", "SD1");
      buildEnclosedCV_(current_spd, "FileFormat", "MS:1001062", "Mascot MGF file", "PSI-MS");
      buildEnclosedCV_(current_spd, "SpectrumIDFormat", "MS:1001528", "Mascot query number", "PSI-MS");
      append_(inputElements, current_spd);
    }

    void MzIdentMLDOMHandler::buildEnclosedCV_(xmlNode* parentElement, const std::string& encel, const std::string& acc, const std::string& name, const std::string& cvref)
    {
      xmlNode* current_ff = element_(parentElement->doc, encel);
      xmlNode* current_cv = element_(current_ff->doc, "cvParam");
      setAttribute_(current_cv, "accession", acc);
      setAttribute_(current_cv, "name", name);
      setAttribute_(current_cv, "cvRef", cvref);
      append_(current_ff, current_cv);
      append_(parentElement, current_ff);
    }

    void MzIdentMLDOMHandler::buildAnalysisDataCollection_(xmlNode* analysisElements)
    {
      xmlNode* current_sil = element_(analysisElements->doc, "SpectrumIdentificationList");
      setAttribute_(current_sil, "id", "SIL1");
      setAttribute_(current_sil, "numSequencesSearched", "TBA");
      // for now no FragmentationTable

      for (PeptideIdentificationList::iterator pi = pep_id_->begin(); pi != pep_id_->end(); ++pi)
      {
        xmlNode* current_sr = element_(current_sil->doc, "SpectrumIdentificationResult");
        setAttribute_(current_sr, "id", StringUtils::toStr(UniqueIdGenerator::getUniqueId()));
        setAttribute_(current_sr, "spectrumID", StringUtils::toStr(UniqueIdGenerator::getUniqueId()));
        setAttribute_(current_sr, "spectraData_ref", "SD1");
        for (vector<PeptideHit>::iterator ph = pi->getHits().begin(); ph != pi->getHits().end(); ++ph)
        {
          xmlNode* current_si = element_(current_sr->doc, "SpectrumIdentificationItem");
          setAttribute_(current_si, "id", StringUtils::toStr(UniqueIdGenerator::getUniqueId()));
          setAttribute_(current_si, "calculatedMassToCharge", StringUtils::toStr(ph->getSequence().getMonoWeight(Residue::Full, ph->getCharge()))); //TODO @mths : this is not correct!1elf - these interfaces are BS!
          setAttribute_(current_si, "chargeState", StringUtils::toStr(ph->getCharge()));
          setAttribute_(current_si, "experimentalMassToCharge", StringUtils::toStr(ph->getSequence().getMonoWeight(Residue::Full, ph->getCharge()))); //TODO @mths : this is not correct!1elf - these interfaces are BS!
          setAttribute_(current_si, "peptide_ref", "TBA");
          setAttribute_(current_si, "rank", StringUtils::toStr(ph->getRank() + 1));
          setAttribute_(current_si, "passThreshold", "TBA");
          setAttribute_(current_si, "sample_ref", "TBA");
          // TODO cvs for score!
          append_(current_sr, current_si);
          for (list<std::string>::iterator pepevref = hit_pev_.front().begin(); pepevref != hit_pev_.front().end(); ++pepevref)
          {
            xmlNode* current_per = element_(current_si->doc, "PeptideEvidenceRef");
            setAttribute_(current_per, "peptideEvidence_ref", *pepevref);
            append_(current_si, current_per);
          }
          hit_pev_.pop_front();
          // and no Fragmentation annotation for now
        }
//          <cvParam accession="MS:1001371" name="Mascot:identity threshold" cvRef="PSI-MS" value="44"/>
//          <cvParam accession="MS:1001370" name="Mascot:homology threshold" cvRef="PSI-MS" value="18"/>
//          <cvParam accession="MS:1001030" name="number of peptide seqs compared to each spectrum" cvRef="PSI-MS" value="26981"/>
//          <cvParam accession="MS:1000796" name="spectrum title" cvRef="PSI-MS" value="dp210198 21-Jan-98 DERIVED SPECTRUM    #9"/>
        append_(current_sil, current_sr);
      }

      // and no ProteinDetection for now
    }

    ProteinIdentification::SearchParameters MzIdentMLDOMHandler::findSearchParameters_(pair<CVTermList, map<std::string, DataValue> > as_params)
    {
      ProteinIdentification::SearchParameters sp = ProteinIdentification::SearchParameters();
      for (std::map<std::string, vector<CVTerm> >::const_iterator cvs = as_params.first.getCVTerms().begin(); cvs != as_params.first.getCVTerms().end(); ++cvs)
      {
        for (vector<CVTerm>::const_iterator cvit = cvs->second.begin(); cvit != cvs->second.end(); ++cvit)
        {
          sp.setMetaValue(cvs->first, cvit->getValue());
        }
      }
      int minCharge = 0;
      int maxCharge = 0;
      for (map<std::string, DataValue>::const_iterator upit = as_params.second.begin(); upit != as_params.second.end(); ++upit)
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
          minCharge = StringUtils::toInt32(upit->second.toString());
        }
        else if (upit->first == "MaxCharge")
        {
          maxCharge = StringUtils::toInt32(upit->second.toString());
        }
        else if (upit->first == "NumTolerableTermini")
        {
          sp.enzyme_term_specificity = static_cast<EnzymaticDigestion::Specificity>(StringUtils::toInt32(upit->second.toString()));
        }
        else
        {
          sp.setMetaValue(upit->first, upit->second);
        }
      }
      if (minCharge != 0 || maxCharge != 0) // this means "MinCharge" and "MaxCharge" get preference over "charges"
      {
        sp.charges =StringUtils::toStr(minCharge) + "-" + StringUtils::toStr(maxCharge);
      }
      return sp;
    }
 
} // namespace OpenMS   //namespace Internal
