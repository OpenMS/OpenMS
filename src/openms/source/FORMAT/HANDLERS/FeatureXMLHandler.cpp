// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/FeatureXMLHandler.h>

#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/DataProcessing.h>

#include <fstream>
#include <map>

using namespace std;

namespace OpenMS::Internal
{
  FeatureXMLHandler::FeatureXMLHandler(FeatureMap& map, const String& filename) :
    Internal::XMLHandler("", "1.9")
  {
    resetMembers_();
    map_ = &map;
    file_ = filename;
  }

  FeatureXMLHandler::FeatureXMLHandler(const FeatureMap& map, const String& filename) :
    Internal::XMLHandler("", "1.9")
  {
    resetMembers_();
    cmap_ = &map;
    file_ = filename;
  }
  FeatureXMLHandler::~FeatureXMLHandler() = default;

  void FeatureXMLHandler::resetMembers_()
  {
    disable_parsing_ = 0;
    current_feature_ = nullptr;
    map_ = nullptr;
    //options_ = FeatureFileOptions(); do NOT reset this, since we need to preserve options!
    size_only_ = false;
    expected_size_ = 0;
    param_ = Param();
    current_chull_ = ConvexHull2D::PointArrayType();
    hull_position_ = DPosition<2>();
    dim_ = 0;
    in_description_ = false;
    subordinate_feature_level_ = 0;
    last_meta_ = nullptr;
    prot_id_ = ProteinIdentification();
    pep_id_ = PeptideIdentification();
    prot_hit_ = ProteinHit();
    pep_hit_ = PeptideHit();
    proteinid_to_accession_.clear();
    accession_to_id_.clear();
    identifier_id_.clear();
    id_identifier_.clear();
    search_param_ = ProteinIdentification::SearchParameters();

  }

  void FeatureXMLHandler::writeTo(std::ostream& os)
  {
    const FeatureMap& feature_map = *(cmap_);

    os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
       << "<featureMap version=\"" << version_ << "\"";
    // file id
    if (!feature_map.getIdentifier().empty())
    {
      os << " document_id=\"" << feature_map.getIdentifier() << "\"";
    }
    // unique id
    if (feature_map.hasValidUniqueId())
    {
      os << " id=\"fm_" << feature_map.getUniqueId() << "\"";
    }
    os << " xsi:noNamespaceSchemaLocation=\"https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/FeatureXML_1_9.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

    // user param
    writeUserParam_("UserParam", os, feature_map, 1);

    // write data processing
    for (Size i = 0; i < feature_map.getDataProcessing().size(); ++i)
    {
      const DataProcessing& processing = feature_map.getDataProcessing()[i];
      os << "\t<dataProcessing completion_time=\"" << processing.getCompletionTime().getDate() << 'T' << processing.getCompletionTime().getTime() << "\">\n";
      os << "\t\t<software name=\"" << processing.getSoftware().getName() << "\" version=\"" << processing.getSoftware().getVersion() << "\" />\n";
      for (set<DataProcessing::ProcessingAction>::const_iterator it = processing.getProcessingActions().begin(); it != processing.getProcessingActions().end(); ++it)
      {
        os << "\t\t<processingAction name=\"" << DataProcessing::NamesOfProcessingAction[*it] << "\" />\n";
      }
      writeUserParam_("UserParam", os, processing, 2);
      os << "\t</dataProcessing>\n";
    }

    // throws if protIDs are not unique, i.e. PeptideIDs will be randomly assigned (bad!)
    checkUniqueIdentifiers_(feature_map.getProteinIdentifications());

    // write identification runs
    Size prot_count = 0;
    for (Size i = 0; i < feature_map.getProteinIdentifications().size(); ++i)
    {
      const ProteinIdentification& current_prot_id = feature_map.getProteinIdentifications()[i];
      os << "\t<IdentificationRun ";
      os << "id=\"PI_" << i << "\" ";
      identifier_id_[current_prot_id.getIdentifier()] = String("PI_") + i;
      os << "date=\"" << current_prot_id.getDateTime().getDate() << "T" << current_prot_id.getDateTime().getTime() << "\" ";
      os << "search_engine=\"" << writeXMLEscape(current_prot_id.getSearchEngine()) << "\" ";
      os << "search_engine_version=\"" << writeXMLEscape(current_prot_id.getSearchEngineVersion()) << "\">\n";

      //write search parameters
      const ProteinIdentification::SearchParameters& search_param = current_prot_id.getSearchParameters();
      os << "\t\t<SearchParameters "
         << "db=\"" << writeXMLEscape(search_param.db) << "\" "
         << "db_version=\"" << writeXMLEscape(search_param.db_version) << "\" "
         << "taxonomy=\"" << writeXMLEscape(search_param.taxonomy) << "\" ";
      if (search_param.mass_type == ProteinIdentification::MONOISOTOPIC)
      {
        os << "mass_type=\"monoisotopic\" ";
      }
      else if (search_param.mass_type == ProteinIdentification::AVERAGE)
      {
        os << "mass_type=\"average\" ";
      }
      os << "charges=\"" << search_param.charges << "\" ";
      String enzyme_name = search_param.digestion_enzyme.getName();
      os << "enzyme=\"" << enzyme_name.toLower() << "\" ";

      String precursor_unit = search_param.precursor_mass_tolerance_ppm ? "true" : "false";
      String peak_unit = search_param.fragment_mass_tolerance_ppm ? "true" : "false";

      os << "missed_cleavages=\"" << search_param.missed_cleavages << "\" "
         << "precursor_peak_tolerance=\"" << search_param.precursor_mass_tolerance << "\" ";
      os << "precursor_peak_tolerance_ppm=\"" << precursor_unit << "\" ";
      os << "peak_mass_tolerance=\"" << search_param.fragment_mass_tolerance << "\" ";
      os << "peak_mass_tolerance_ppm=\"" << peak_unit << "\" ";
      os << ">\n";

      //modifications
      for (Size j = 0; j != search_param.fixed_modifications.size(); ++j)
      {
        os << "\t\t\t<FixedModification name=\"" << writeXMLEscape(search_param.fixed_modifications[j]) << "\" />\n";
        //Add MetaInfo, when modifications has it (Andreas)
      }
      for (Size j = 0; j != search_param.variable_modifications.size(); ++j)
      {
        os << "\t\t\t<VariableModification name=\"" << writeXMLEscape(search_param.variable_modifications[j]) << "\" />\n";
        //Add MetaInfo, when modifications has it (Andreas)
      }

      writeUserParam_("UserParam", os, search_param, 3);

      os << "\t\t</SearchParameters>\n";

      //write protein identifications
      os << "\t\t<ProteinIdentification";
      os << " score_type=\"" << writeXMLEscape(current_prot_id.getScoreType()) << "\"";
      os << " higher_score_better=\"" << (current_prot_id.isHigherScoreBetter() ? "true" : "false") << "\"";
      os << " significance_threshold=\"" << current_prot_id.getSignificanceThreshold() << "\">\n";

      // write protein hits
      for (Size j = 0; j < current_prot_id.getHits().size(); ++j)
      {
        os << "\t\t\t<ProteinHit";

        // prot_count
        os << " id=\"PH_" << prot_count << "\"";
        accession_to_id_[current_prot_id.getIdentifier() + "_" + current_prot_id.getHits()[j].getAccession()] = prot_count;
        ++prot_count;

        os << " accession=\"" << writeXMLEscape(current_prot_id.getHits()[j].getAccession()) << "\"";
        os << " score=\"" << current_prot_id.getHits()[j].getScore() << "\"";
        
        double coverage = current_prot_id.getHits()[j].getCoverage();
        if (coverage != ProteinHit::COVERAGE_UNKNOWN)
        {
          os << " coverage=\"" << coverage << "\"";
        }
        
        os << " sequence=\"" << writeXMLEscape(current_prot_id.getHits()[j].getSequence()) << "\">\n";

        writeUserParam_("UserParam", os, current_prot_id.getHits()[j], 4);

        os << "\t\t\t</ProteinHit>\n";
      }

      writeUserParam_("UserParam", os, current_prot_id, 3);
      os << "\t\t</ProteinIdentification>\n";
      os << "\t</IdentificationRun>\n";
    }

    //write unassigned peptide identifications
    for (Size i = 0; i < feature_map.getUnassignedPeptideIdentifications().size(); ++i)
    {
      writePeptideIdentification_(file_, os, feature_map.getUnassignedPeptideIdentifications()[i], "UnassignedPeptideIdentification", 1);
    }

    // write features with their corresponding attributes
    os << "\t<featureList count=\"" << feature_map.size() << "\">\n";
    startProgress(0, feature_map.size(), "Storing featureXML file");
    for (Size s = 0; s < feature_map.size(); s++)
    {
      writeFeature_(file_, os, feature_map[s], "f_", feature_map[s].getUniqueId(), 0);
      setProgress(s);
      // writeFeature_(file_, os, feature_map[s], "f_", s, 0);
    }
    endProgress();

    os << "\t</featureList>\n";
    os << "</featureMap>\n";

    //Clear members
    accession_to_id_.clear();
    identifier_id_.clear();
  }

  FeatureFileOptions& FeatureXMLHandler::getOptions()
  {
    return options_;
  }

  const FeatureFileOptions& FeatureXMLHandler::getOptions() const
  {
    return options_;
  }

  void FeatureXMLHandler::setOptions(const FeatureFileOptions& options)
  {
    options_ = options;
  }

  void FeatureXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  {
    static const XMLCh* s_dim = xercesc::XMLString::transcode("dim");
    static const XMLCh* s_name = xercesc::XMLString::transcode("name");
    static const XMLCh* s_version = xercesc::XMLString::transcode("version");
    static const XMLCh* s_value = xercesc::XMLString::transcode("value");
    static const XMLCh* s_type = xercesc::XMLString::transcode("type");
    static const XMLCh* s_completion_time = xercesc::XMLString::transcode("completion_time");
    static const XMLCh* s_document_id = xercesc::XMLString::transcode("document_id");
    static const XMLCh* s_id = xercesc::XMLString::transcode("id");

    // TODO The next line should be removed in OpenMS 1.7 or so!
    static const XMLCh* s_unique_id = xercesc::XMLString::transcode("unique_id");

    String tag = sm_.convert(qname);

    // handle skipping of whole sections
    // IMPORTANT: check parent tags first (i.e. tags higher in the tree), since otherwise sections might be enabled/disabled too early/late
    //            disable_parsing_ is an Int, since subordinates might be chained, thus at SO-level 2, the endelement() would switch on parsing again
    //                                      , even though the end of the parent SO was not reached
    if ((!options_.getLoadSubordinates()) && tag == "subordinate")
    {
      ++disable_parsing_;
    }
    else if ((!options_.getLoadConvexHull()) && tag == "convexhull")
    {
      ++disable_parsing_;
    }
    if (disable_parsing_)
    {
      return;
    }
    // do the actual parsing:
    String parent_tag;
    if (!open_tags_.empty())
    {
      parent_tag = open_tags_.back();
    }
    open_tags_.push_back(tag);

    //for downward compatibility, all tags in the old description must be ignored
    if (in_description_)
    {
      return;
    }
    if (tag == "description")
    {
      in_description_ = true;
    }
    else if (tag == "feature")
    {
      // create new feature at appropriate level
      updateCurrentFeature_(true);
      current_feature_->setUniqueId(attributeAsString_(attributes, s_id));
    }
    else if (tag == "subordinate") // this is not safe towards malformed xml!
    {
      ++subordinate_feature_level_;
    }
    else if (tag == "featureList")
    {
      if (options_.getMetadataOnly())
      {
        throw EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
      Size count = attributeAsInt_(attributes, "count");
      if (size_only_) // true if loadSize() was used instead of load()
      {
        expected_size_ = count;
        throw EndParsingSoftly(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION);
      }
      map_->reserve(std::min(Size(1e5), count)); // reserve vector for faster push_back, but with upper boundary of 1e5 (as >1e5 is most likely an invalid feature count)
      startProgress(0, count, "Loading featureXML file");
    }
    else if (tag == "quality" || tag == "hposition" || tag == "position")
    {
      dim_ = attributeAsInt_(attributes, s_dim);
    }
    else if (tag == "pt")
    {
      hull_position_[0] = attributeAsDouble_(attributes, "x");
      hull_position_[1] = attributeAsDouble_(attributes, "y");
    }
    else if (tag == "convexhull")
    {
      current_chull_.clear();
    }
    else if (tag == "hullpoint")
    {
      hull_position_ = DPosition<2>::zero();
    }
    else if (tag == "param")
    {
      String name = attributeAsString_(attributes, s_name);
      String value = attributeAsString_(attributes, s_value);
      if (!name.empty() && !value.empty())
        param_.setValue(name, value);
    }
    else if (tag == "userParam" || tag == "UserParam") // correct: "UserParam". Test for backwards compatibility.
    {
      if (last_meta_ == nullptr)
      {
        fatalError(LOAD, String("Unexpected UserParam in tag '") + parent_tag + "'");
      }

      String name = attributeAsString_(attributes, s_name);
      String type = attributeAsString_(attributes, s_type);

      if (type == "int")
      {
        last_meta_->setMetaValue(name, attributeAsInt_(attributes, s_value));
      }
      else if (type == "float")
      {
        last_meta_->setMetaValue(name, attributeAsDouble_(attributes, s_value));
      }
      else if (type == "string")
      {
        last_meta_->setMetaValue(name, (String)attributeAsString_(attributes, s_value));
      }
      else if (type == "intList")
      {
        last_meta_->setMetaValue(name, attributeAsIntList_(attributes, "value"));
      }
      else if (type == "floatList")
      {
        last_meta_->setMetaValue(name, attributeAsDoubleList_(attributes, "value"));
      }
      else if (type == "stringList")
      {
        last_meta_->setMetaValue(name, attributeAsStringList_(attributes, "value"));
      }
      else
      {
        fatalError(LOAD, String("Invalid UserParam type '") + type + "'");
      }
    }
    else if (tag == "featureMap")
    {
      //check file version against schema version
      String file_version = "";
      optionalAttributeAsString_(file_version, attributes, s_version);
      if (file_version.empty())
      {
        file_version = "1.0"; //default version is 1.0
      }
      if (file_version.toDouble() > version_.toDouble())
      {
        warning(LOAD, String("The XML file (") + file_version + ") is newer than the parser (" + version_ + "). This might lead to undefined program behavior.");
      }
      //handle document id
      String document_id;
      if (optionalAttributeAsString_(document_id, attributes, s_document_id))
      {
        map_->setIdentifier(document_id);
      }
      //handle unique id
      String unique_id;
      if (optionalAttributeAsString_(unique_id, attributes, s_id))
      {
        map_->setUniqueId(unique_id);
      }
      // TODO The next four lines should be removed in OpenMS 1.7 or so!
      if (optionalAttributeAsString_(unique_id, attributes, s_unique_id))
      {
        map_->setUniqueId(unique_id);
      }
      last_meta_ = map_;
    }
    else if (tag == "dataProcessing")
    {
      DataProcessing tmp;
      tmp.setCompletionTime(asDateTime_(attributeAsString_(attributes, s_completion_time)));
      map_->getDataProcessing().push_back(tmp);
      last_meta_ = &(map_->getDataProcessing().back());
    }
    else if (tag == "software" && parent_tag == "dataProcessing")
    {
      map_->getDataProcessing().back().getSoftware().setName(attributeAsString_(attributes, s_name));
      map_->getDataProcessing().back().getSoftware().setVersion(attributeAsString_(attributes, s_version));
    }
    else if (tag == "processingAction" && parent_tag == "dataProcessing")
    {
      String name = attributeAsString_(attributes, s_name);
      for (Size i = 0; i < DataProcessing::SIZE_OF_PROCESSINGACTION; ++i)
      {
        if (name == DataProcessing::NamesOfProcessingAction[i])
        {
          map_->getDataProcessing().back().getProcessingActions().insert((DataProcessing::ProcessingAction)i);
        }
      }
    }
    else if (tag == "IdentificationRun")
    {
      prot_id_.setSearchEngine(attributeAsString_(attributes, "search_engine"));
      prot_id_.setSearchEngineVersion(attributeAsString_(attributes, "search_engine_version"));
      prot_id_.setDateTime(DateTime::fromString(attributeAsString_(attributes, "date")));
      // set identifier
      // always generate a unique id to link a ProteinIdentification and the corresponding PeptideIdentifications
      // , since any FeatureLinker might just carelessly concatenate PepIDs from different FeatureMaps.
      // If these FeatureMaps have identical identifiers (SearchEngine time + type match exactly), then ALL PepIDs would be falsely attributed
      // to a single ProtID...

      String id = attributeAsString_(attributes, "id");
      while (true)
      { // loop until the identifier is unique (should be on the first iteration -- very(!) unlikely it will not be unique)
        // Note: technically, it would be preferable to prefix the UID for faster string comparison, but this results in random write-orderings during file store (breaks tests)
        String identifier = prot_id_.getSearchEngine() + '_' + attributeAsString_(attributes, "date") + '_' + String(UniqueIdGenerator::getUniqueId());

        if (id_identifier_.find(id) == id_identifier_.end())
        {
          prot_id_.setIdentifier(identifier);
          id_identifier_[id] = identifier;
          break;
        }
      }
    }
    else if (tag == "SearchParameters")
    {
      //load parameters
      search_param_.db = attributeAsString_(attributes, "db");
      search_param_.db_version = attributeAsString_(attributes, "db_version");
      optionalAttributeAsString_(search_param_.taxonomy, attributes, "taxonomy");
      search_param_.charges = attributeAsString_(attributes, "charges");
      optionalAttributeAsUInt_(search_param_.missed_cleavages, attributes, "missed_cleavages");
      search_param_.fragment_mass_tolerance = attributeAsDouble_(attributes, "peak_mass_tolerance");
      String peak_unit;
      optionalAttributeAsString_(peak_unit, attributes, "peak_mass_tolerance_ppm");
      search_param_.fragment_mass_tolerance_ppm = peak_unit == "true" ? true : false;
      search_param_.precursor_mass_tolerance = attributeAsDouble_(attributes, "precursor_peak_tolerance");
      String precursor_unit;
      optionalAttributeAsString_(precursor_unit, attributes, "precursor_peak_tolerance_ppm");
      search_param_.precursor_mass_tolerance_ppm = precursor_unit == "true" ? true : false;
      //mass type
      String mass_type = attributeAsString_(attributes, "mass_type");
      if (mass_type == "monoisotopic")
      {
        search_param_.mass_type = ProteinIdentification::MONOISOTOPIC;
      }
      else if (mass_type == "average")
      {
        search_param_.mass_type = ProteinIdentification::AVERAGE;
      }
      //enzyme
      String enzyme;
      optionalAttributeAsString_(enzyme, attributes, "enzyme");
      if (ProteaseDB::getInstance()->hasEnzyme(enzyme))
      {
        search_param_.digestion_enzyme = *(ProteaseDB::getInstance()->getEnzyme(enzyme));
      }
      last_meta_ = &search_param_;
    }
    else if (tag == "FixedModification")
    {
      search_param_.fixed_modifications.push_back(attributeAsString_(attributes, "name"));
      //change this line as soon as there is a MetaInfoInterface for modifications (Andreas)
      last_meta_ = nullptr;
    }
    else if (tag == "VariableModification")
    {
      search_param_.variable_modifications.push_back(attributeAsString_(attributes, "name"));
      //change this line as soon as there is a MetaInfoInterface for modifications (Andreas)
      last_meta_ = nullptr;
    }
    else if (tag == "ProteinIdentification")
    {
      prot_id_.setScoreType(attributeAsString_(attributes, "score_type"));

      //optional significance threshold
      double tmp = 0.0;
      optionalAttributeAsDouble_(tmp, attributes, "significance_threshold");
      if (tmp != 0.0)
      {
        prot_id_.setSignificanceThreshold(tmp);
      }

      //score orientation
      prot_id_.setHigherScoreBetter(asBool_(attributeAsString_(attributes, "higher_score_better")));

      last_meta_ = &prot_id_;
    }
    else if (tag == "ProteinHit")
    {
      prot_hit_ = ProteinHit();
      String accession = attributeAsString_(attributes, "accession");
      prot_hit_.setAccession(accession);
      prot_hit_.setScore(attributeAsDouble_(attributes, "score"));

      // coverage
      double coverage = -std::numeric_limits<double>::max();
      optionalAttributeAsDouble_(coverage, attributes, "coverage");
      if (coverage != -std::numeric_limits<double>::max())
      {
        prot_hit_.setCoverage(coverage);
      }

      //sequence
      String tmp = "";
      optionalAttributeAsString_(tmp, attributes, "sequence");
      prot_hit_.setSequence(tmp);

      last_meta_ = &prot_hit_;

      //insert id and accession to map
      proteinid_to_accession_[attributeAsString_(attributes, "id")] = accession;
    }
    else if (tag == "PeptideIdentification" || tag == "UnassignedPeptideIdentification")
    {
      String id = attributeAsString_(attributes, "identification_run_ref");
      if (id_identifier_.find(id) == id_identifier_.end())
      {
        warning(LOAD, String("Peptide identification without ProteinIdentification found (id: '") + id + "')!");
      }
      pep_id_.setIdentifier(id_identifier_[id]);

      pep_id_.setScoreType(attributeAsString_(attributes, "score_type"));

      //optional significance threshold
      double tmp = 0.0;
      optionalAttributeAsDouble_(tmp, attributes, "significance_threshold");
      if (tmp != 0.0)
      {
        pep_id_.setSignificanceThreshold(tmp);
      }

      //score orientation
      pep_id_.setHigherScoreBetter(asBool_(attributeAsString_(attributes, "higher_score_better")));

      //MZ
      double tmp2 = -numeric_limits<double>::max();
      optionalAttributeAsDouble_(tmp2, attributes, "MZ");
      if (tmp2 != -numeric_limits<double>::max())
      {
        pep_id_.setMZ(tmp2);
      }
      //RT
      tmp2 = -numeric_limits<double>::max();
      optionalAttributeAsDouble_(tmp2, attributes, "RT");
      if (tmp2 != -numeric_limits<double>::max())
      {
        pep_id_.setRT(tmp2);
      }
      String tmp3;
      optionalAttributeAsString_(tmp3, attributes, "spectrum_reference");
      if (!tmp3.empty())
      {
        pep_id_.setSpectrumReference( tmp3);
      }

      last_meta_ = &pep_id_;
    }
    else if (tag == "PeptideHit")
    {
      pep_hit_ = PeptideHit();
      vector<PeptideEvidence> peptide_evidences_;

      pep_hit_.setCharge(attributeAsInt_(attributes, "charge"));
      pep_hit_.setScore(attributeAsDouble_(attributes, "score"));
      pep_hit_.setSequence(AASequence::fromString(String(attributeAsString_(attributes, "sequence"))));

      //parse optional protein ids to determine accessions
      const XMLCh* refs = attributes.getValue(sm_.convert("protein_refs").c_str());
      if (refs != nullptr)
      {
        String accession_string = sm_.convert(refs);
        accession_string.trim();
        vector<String> accessions;
        accession_string.split(' ', accessions);
        if (!accession_string.empty() && accessions.empty())
        {
          accessions.push_back(accession_string);
        }
        for (vector<String>::const_iterator it = accessions.begin(); it != accessions.end(); ++it)
        {
          std::map<String, String>::const_iterator it2 = proteinid_to_accession_.find(*it);
          if (it2 != proteinid_to_accession_.end())
          {
            PeptideEvidence pe;
            pe.setProteinAccession(it2->second);
            peptide_evidences_.push_back(pe);
          }
          else
          {
            fatalError(LOAD, String("Invalid protein reference '") + *it + "'");
          }
        }
      }

      //aa_before
      String tmp = "";
      optionalAttributeAsString_(tmp, attributes, "aa_before");
      if (!tmp.empty())
      {
        std::vector<String> splitted;
        tmp.split(' ', splitted);
        for (Size i = 0; i != splitted.size(); ++i)
        { 
          if (peptide_evidences_.size() < i + 1) 
          {
            peptide_evidences_.emplace_back();
          }
          peptide_evidences_[i].setAABefore(splitted[i][0]);
        }
      }

      //aa_after
      tmp = "";
      optionalAttributeAsString_(tmp, attributes, "aa_after");
      if (!tmp.empty())
      {
        std::vector<String> splitted;
        tmp.split(' ', splitted);
        for (Size i = 0; i != splitted.size(); ++i)
        { 
          if (peptide_evidences_.size() < i + 1) 
          {
            peptide_evidences_.emplace_back();
          }
          peptide_evidences_[i].setAAAfter(splitted[i][0]);
        }
      }

      //start
      tmp = "";
      optionalAttributeAsString_(tmp, attributes, "start");

      if (!tmp.empty())
      {
        std::vector<String> splitted;
        tmp.split(' ', splitted);
        for (Size i = 0; i != splitted.size(); ++i)
        { 
          if (peptide_evidences_.size() < i + 1) 
          {
            peptide_evidences_.emplace_back();
          }
          peptide_evidences_[i].setStart(splitted[i].toInt());
        }
      }

      //end
      tmp = "";
      optionalAttributeAsString_(tmp, attributes, "end");
      if (!tmp.empty())
      {
        std::vector<String> splitted;
        tmp.split(' ', splitted);
        for (Size i = 0; i != splitted.size(); ++i)
        { 
          if (peptide_evidences_.size() < i + 1) 
          {
            peptide_evidences_.emplace_back();
          }
          peptide_evidences_[i].setEnd(splitted[i].toInt());
        }
      }

      pep_hit_.setPeptideEvidences(peptide_evidences_);
      last_meta_ = &pep_hit_;
    }
  }

  void FeatureXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    String tag = sm_.convert(qname);

    // handle skipping of whole sections
    // IMPORTANT: check parent tags first (i.e. tags higher in the tree), since otherwise sections might be enabled/disabled too early/late
    if (((!options_.getLoadSubordinates()) && tag == "subordinate")
       || ((!options_.getLoadConvexHull()) && tag == "convexhull"))
    {
      --disable_parsing_;
      return; // even if disable_parsing is false now, we still exit (since this endelement() should be ignored)
    }

    if (disable_parsing_)
    {
      return;
    }
    // do the actual parsing:
    open_tags_.pop_back();

    //for downward compatibility, all tags in the old description must be ignored
    if (tag == "description")
    {
      in_description_ = false;
    }
    if (in_description_)
    {
      return;
    }
    if (tag == "feature")
    {
      if ((!options_.hasRTRange() || options_.getRTRange().encloses(current_feature_->getRT()))
         &&  (!options_.hasMZRange() || options_.getMZRange().encloses(current_feature_->getMZ()))
         &&  (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(current_feature_->getIntensity())))
      {
      }
      else
      {
        // this feature does not pass the restrictions --> remove it
        if (subordinate_feature_level_ == 0)
        {
          map_->pop_back();
        }
        else
        {
          Feature* f1(nullptr);
          if (!map_->empty())
          {
            f1 = &(map_->back());
          }
          else
          {
            fatalError(LOAD, "Feature with unexpected location.");
          }

          for (Int level = 1; level < subordinate_feature_level_; ++level)
          {
            f1 = &(f1->getSubordinates().back());
          }
          // delete the offending feature
          f1->getSubordinates().pop_back();
        }
      }
      updateCurrentFeature_(false);
    }
    else if (tag == "model")
    {
      warning(LOAD, String("The featureXML file contains a 'model' description, but the internal datastructure has no model support since OpenMS 1.12. Model will be ignored!"));
    }
    else if (tag == "hullpoint" || tag == "pt")
    {
      current_chull_.push_back(hull_position_);
    }
    else if (tag == "convexhull")
    {
      ConvexHull2D hull;
      hull.setHullPoints(current_chull_);
      current_feature_->getConvexHulls().push_back(hull);
    }
    else if (tag == "subordinate")
    {
      --subordinate_feature_level_;
      // reset current_feature
      updateCurrentFeature_(false);
    }
    else if (tag == "IdentificationRun")
    {
      map_->getProteinIdentifications().push_back(prot_id_);
      prot_id_ = ProteinIdentification();
      last_meta_  = nullptr;
    }
    else if (tag == "SearchParameters")
    {
      prot_id_.setSearchParameters(search_param_);
      search_param_ = ProteinIdentification::SearchParameters();
    }
    else if (tag == "FixedModification")
    {
      last_meta_ = &search_param_;
    }
    else if (tag == "VariableModification")
    {
      last_meta_ = &search_param_;
    }

    else if (tag == "ProteinHit")
    {
      prot_id_.insertHit(prot_hit_);
      last_meta_ = &prot_id_;
    }
    else if (tag == "PeptideIdentification")
    {
      current_feature_->getPeptideIdentifications().push_back(pep_id_);
      pep_id_ = PeptideIdentification();
      last_meta_  = &map_->back();
    }
    else if (tag == "UnassignedPeptideIdentification")
    {
      map_->getUnassignedPeptideIdentifications().push_back(pep_id_);
      pep_id_ = PeptideIdentification();
      last_meta_  = nullptr;
    }
    else if (tag == "PeptideHit")
    {
      pep_id_.insertHit(pep_hit_);
      last_meta_ = &pep_id_;
    }
    else if (tag == "featureList")
    {
      endProgress();
    }
  }

  void FeatureXMLHandler::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
  {
    // handle skipping of whole sections
    if (disable_parsing_)
    {
      return;
    }
    // do the actual parsing:

    //for downward compatibility, all tags in the old description must be ignored
    if (in_description_)
    {
      return;
    }
    // we are before first tag or beyond last tag
    if (open_tags_.empty())
    {
      return;
    }
    String& current_tag = open_tags_.back();
    if (current_tag == "intensity")
    {
      current_feature_->setIntensity(asDouble_(sm_.convert(chars)));
    }
    else if (current_tag == "position")
    {
      current_feature_->getPosition()[dim_] = asDouble_(sm_.convert(chars));
    }
    else if (current_tag == "quality")
    {
      current_feature_->setQuality(dim_, asDouble_(sm_.convert(chars)));
    }
    else if (current_tag == "overallquality")
    {
      current_feature_->setOverallQuality(asDouble_(sm_.convert(chars)));
    }
    else if (current_tag == "charge")
    {
      current_feature_->setCharge(asInt_(chars));
    }
    else if (current_tag == "hposition")
    {
      hull_position_[dim_] = asDouble_(sm_.convert(chars));
    }
  }

  void FeatureXMLHandler::writeFeature_(const String& filename, ostream& os, const Feature& feat, const String& identifier_prefix, UInt64 identifier, UInt indentation_level)
  {
    String indent = String(indentation_level, '\t');

    os << indent << "\t\t<feature id=\"" << identifier_prefix << identifier << "\">\n";
    for (Size i = 0; i < 2; i++)
    {
      os << indent << "\t\t\t<position dim=\"" << i << "\">" << precisionWrapper(feat.getPosition()[i]) << "</position>\n";
    }
    os << indent << "\t\t\t<intensity>" << precisionWrapper(feat.getIntensity()) << "</intensity>\n";
    for (Size i = 0; i < 2; i++)
    {
      os << indent << "\t\t\t<quality dim=\"" << i << "\">" << precisionWrapper(feat.getQuality(i)) << "</quality>\n";
    }
    os << indent << "\t\t\t<overallquality>" << precisionWrapper(feat.getOverallQuality()) << "</overallquality>\n";
    os << indent << "\t\t\t<charge>" << feat.getCharge() << "</charge>\n";

    // write convex hull
    vector<ConvexHull2D> hulls = feat.getConvexHulls();

    Size hulls_count = hulls.size();

    for (Size i = 0; i < hulls_count; i++)
    {
      os << indent << "\t\t\t<convexhull nr=\"" << i << "\">\n";

      ConvexHull2D current_hull = hulls[i];
      current_hull.compress();
      Size hull_size = current_hull.getHullPoints().size();

      for (Size j = 0; j < hull_size; j++)
      {
        DPosition<2> pos = current_hull.getHullPoints()[j];
        /*Size pos_size = pos.size();
            os << indent << "\t\t\t\t<hullpoint>\n";
    for (Size k=0; k<pos_size; k++)
            {
                os << indent << "\t\t\t\t\t<hposition dim=\"" << k << "\">" << precisionWrapper(pos[k]) << "</hposition>\n";
            }
            os << indent << "\t\t\t\t</hullpoint>\n";*/
        os << indent << "\t\t\t\t<pt x=\"" << precisionWrapper(pos[0]) << "\" y=\"" << precisionWrapper(pos[1]) << "\" />\n";
      }

      os << indent << "\t\t\t</convexhull>\n";
    }

    if (!feat.getSubordinates().empty())
    {
      os << indent << "\t\t\t<subordinate>\n";
      for (size_t i = 0; i < feat.getSubordinates().size(); ++i)
      {
        // These subordinate identifiers are a bit long, but who cares about subordinates anyway?  :-P
        // This way the parent stands out clearly.  However,
        // note that only the portion after the last '_' is parsed when this is read back.
        writeFeature_(filename, os, feat.getSubordinates()[i], identifier_prefix + identifier + "_", feat.getSubordinates()[i].getUniqueId(), indentation_level + 2);
      }
      os << indent << "\t\t\t</subordinate>\n";
    }

    // write PeptideIdentification
    for (Size i = 0; i < feat.getPeptideIdentifications().size(); ++i)
    {
      writePeptideIdentification_(filename, os, feat.getPeptideIdentifications()[i], "PeptideIdentification", 3);
    }

    writeUserParam_("UserParam", os, feat, indentation_level + 3);

    os << indent << "\t\t</feature>\n";
  }

  void FeatureXMLHandler::writePeptideIdentification_(const String& filename, std::ostream& os, const PeptideIdentification& id, const String& tag_name, UInt indentation_level)
  {
    String indent = String(indentation_level, '\t');

    if (identifier_id_.find(id.getIdentifier()) == identifier_id_.end())
    {
      warning(STORE, String("Omitting peptide identification because of missing ProteinIdentification with identifier '") + id.getIdentifier() + "' while writing '" + filename + "'!");
      return;
    }
    os << indent << "<" << tag_name << " ";
    os << "identification_run_ref=\"" << identifier_id_[id.getIdentifier()] << "\" ";
    os << "score_type=\"" << writeXMLEscape(id.getScoreType()) << "\" ";
    os << "higher_score_better=\"" << (id.isHigherScoreBetter() ? "true" : "false") << "\" ";
    os << "significance_threshold=\"" << id.getSignificanceThreshold() << "\" ";
    //mz
    if (id.hasMZ())
    {
      os << "MZ=\"" << id.getMZ() << "\" ";
    }
    // rt
    if (id.hasRT())
    {
      os << "RT=\"" << id.getRT() << "\" ";
    }
    // spectrum_reference
    DataValue dv = id.getMetaValue("spectrum_reference");
    if (dv != DataValue::EMPTY)
    {
      os << "spectrum_reference=\"" << writeXMLEscape(dv.toString()) << "\" ";
    }
    os << ">\n";

    // write peptide hits
    for (Size j = 0; j < id.getHits().size(); ++j)
    {
      const PeptideHit& h = id.getHits()[j];
      os << indent << "\t<PeptideHit";
      os << " score=\"" << h.getScore() << "\"";
      os << " sequence=\"" << writeXMLEscape(h.getSequence().toString()) << "\"";
      os << " charge=\"" << h.getCharge() << "\"";

      const vector<PeptideEvidence>& pes = id.getHits()[j].getPeptideEvidences();

      IdXMLFile::createFlankingAAXMLString_(pes, os);
      IdXMLFile::createPositionXMLString_(pes, os);

      String accs;
      for (vector<PeptideEvidence>::const_iterator pe = pes.begin(); pe != pes.end(); ++pe)
      {
        if (!accs.empty())
        {
          accs += " ";
        }
        String protein_accession = pe->getProteinAccession();

        // empty accessions are not written out (legacy code)
        if (!protein_accession.empty())
        {
          accs += "PH_";
          accs += String(accession_to_id_[id.getIdentifier() + "_" + protein_accession]);
        }
      }

      // don't write protein_refs if no peptide evidences present
      if (!accs.empty())
      {
        os << " protein_refs=\"" << accs << "\"";
      }

      os << ">\n";

      writeUserParam_("UserParam", os, id.getHits()[j], indentation_level + 2);
      os << indent << "\t</PeptideHit>\n";
    }

    //do not write "spectrum_reference" since it is written as attribute already
    MetaInfoInterface tmp = id;
    tmp.removeMetaValue("spectrum_reference");
    writeUserParam_("UserParam", os, tmp, indentation_level + 1);
    os << indent << "</" << tag_name << ">\n";
  }

  void FeatureXMLHandler::updateCurrentFeature_(bool create)
  {
    if (subordinate_feature_level_ == 0)
    {
      if (create)
      {
        setProgress(map_->size());
        map_->push_back(Feature());
        current_feature_ = &map_->back();
        last_meta_ =  &map_->back();
      }
      else
      {
        if (map_->empty())
        {
          current_feature_ = nullptr;
          last_meta_ =  nullptr;
        }
        else
        {
          current_feature_ = &map_->back();
          last_meta_ =  &map_->back();
        }
      }
      return;
    }

    Feature* f1 = nullptr;
    if (map_->empty())
    {
      // do NOT throw an exception here. this is a valid case! e.g. the
      // only one feature in a map was discarded during endElement(), thus
      // the map_ is empty() now and we cannot assign a current_feature,
      // because there is none!
      current_feature_ = nullptr;
      last_meta_ = nullptr;
      return;
    }
    else
    {
      f1 = &map_->back();
    }

    for (Int level = 1; level < subordinate_feature_level_; ++level)
    {
      // if all features of the current level are discarded (due to
      // range-restrictions etc), then the current feature is the one which
      // is one level up
      if (f1->getSubordinates().empty())
      {
        current_feature_ = f1;
        last_meta_ = f1;
        return;
      }
      f1 = &f1->getSubordinates().back();
    }
    if (create)
    {
      f1->getSubordinates().emplace_back();
      current_feature_ = &f1->getSubordinates().back();
      last_meta_ = &f1->getSubordinates().back();
      return;
    }
    else
    {
      if (f1->getSubordinates().empty())
      {
        current_feature_ = nullptr;
        last_meta_ = nullptr;
        return;
      }
      else
      {
        current_feature_ = &f1->getSubordinates().back();
        last_meta_ = &f1->getSubordinates().back();
        return;
      }
    }
  }

}
