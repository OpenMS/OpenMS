// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm, Mathias Walzer $
// --------------------------------------------------------------------------


#include <OpenMS/FORMAT/HANDLERS/ConsensusXMLHandler.h>

#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/SYSTEM/File.h>

#include <map>
#include <fstream>

using namespace std;

namespace OpenMS::Internal
{
  ConsensusXMLHandler::ConsensusXMLHandler(ConsensusMap& map, const String& filename) :
    XMLHandler("", "1.7"),
    ProgressLogger(),
    act_cons_element_(),
    last_meta_(nullptr)
  {
    consensus_map_ = &map;
    file_ = filename;
  }

  ConsensusXMLHandler::ConsensusXMLHandler(const ConsensusMap& map, const String& filename) :
    XMLHandler("", "1.7"),
    ProgressLogger(),
    act_cons_element_(),
    last_meta_(nullptr)
  {
    cconsensus_map_ = &map;
    file_ = filename;
  }

  ConsensusXMLHandler::~ConsensusXMLHandler() = default;

  PeakFileOptions& ConsensusXMLHandler::getOptions()
  {
    return options_;
  }

  void ConsensusXMLHandler::setOptions(const PeakFileOptions& options)
  {
    options_ = options;
  }

  const PeakFileOptions& ConsensusXMLHandler::getOptions() const
  {
    return options_;
  }

  void ConsensusXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    String tag = sm_.convert(qname);
    open_tags_.pop_back();

    if (tag == "consensusElement")
    {
      if ((!options_.hasRTRange() || options_.getRTRange().encloses(act_cons_element_.getRT())) && (!options_.hasMZRange() || options_.getMZRange().encloses(
                                                                                                      act_cons_element_.getMZ())) && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(act_cons_element_.getIntensity())))
      {
        consensus_map_->push_back(act_cons_element_);
        act_cons_element_.getPeptideIdentifications().clear();
      }
      last_meta_ = nullptr;
    }
    else if (tag == "IdentificationRun")
    {
      // post processing of ProteinGroups (hack)
      getProteinGroups_(prot_id_.getProteinGroups(), "protein_group");
      getProteinGroups_(prot_id_.getIndistinguishableProteins(),
                        "indistinguishable_proteins");
      consensus_map_->getProteinIdentifications().emplace_back(std::move(prot_id_));
      prot_id_ = ProteinIdentification();
      last_meta_ = nullptr;
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
      act_cons_element_.getPeptideIdentifications().emplace_back(std::move(pep_id_));
      pep_id_ = PeptideIdentification();
      last_meta_ = &act_cons_element_;
    }
    else if (tag == "UnassignedPeptideIdentification")
    {
      consensus_map_->getUnassignedPeptideIdentifications().emplace_back(std::move(pep_id_));
      pep_id_ = PeptideIdentification();
      last_meta_ = consensus_map_;
    }
    else if (tag == "PeptideHit")
    {
      pep_hit_.setPeptideEvidences(peptide_evidences_);
      pep_id_.insertHit(pep_hit_);
      last_meta_ = &pep_id_;
    }
    else if (tag == "consensusXML")
    {
      endProgress();
    }
  }

  void ConsensusXMLHandler::characters(const XMLCh* const /*chars*/, const XMLSize_t /*length*/)
  {
  }

  void ConsensusXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  {
    const String& parent_tag = (open_tags_.empty() ? "" : open_tags_.back());
    open_tags_.push_back(sm_.convert(qname));
    const String& tag = open_tags_.back();

    String tmp_str;
    if (tag == "map")
    {
      setProgress(++progress_);
      Size last_map = attributeAsInt_(attributes, "id");
      last_meta_ = &consensus_map_->getColumnHeaders()[last_map];
      consensus_map_->getColumnHeaders()[last_map].filename = attributeAsString_(attributes, "name");
      String unique_id;
      if (XMLHandler::optionalAttributeAsString_(unique_id, attributes, "unique_id"))
      {
        UniqueIdInterface tmp;
        tmp.setUniqueId(unique_id);
        consensus_map_->getColumnHeaders()[last_map].unique_id = tmp.getUniqueId();
      }
      String label;
      if (XMLHandler::optionalAttributeAsString_(label, attributes, "label"))
      {
        consensus_map_->getColumnHeaders()[last_map].label = label;
      }
      UInt size;
      if (XMLHandler::optionalAttributeAsUInt_(size, attributes, "size"))
      {
        consensus_map_->getColumnHeaders()[last_map].size = size;
      }
    }
    else if (tag == "consensusElement")
    {
      setProgress(++progress_);
      act_cons_element_ = ConsensusFeature();
      last_meta_ = &act_cons_element_;
      // quality
      if (double quality; optionalAttributeAsDouble_(quality, attributes, "quality"))
      {
        act_cons_element_.setQuality(quality);
      }
      // charge
      if (Int charge; optionalAttributeAsInt_(charge, attributes, "charge"))
      {
        act_cons_element_.setCharge(charge);
      }
      // unique id
      act_cons_element_.setUniqueId(attributeAsString_(attributes, "id"));
      last_meta_ = &act_cons_element_;
    }
    else if (tag == "centroid")
    {
      tmp_str = attributeAsString_(attributes, "rt");
      if (!tmp_str.empty())
      {
        pos_[Peak2D::RT] = asDouble_(tmp_str);
      }

      tmp_str = attributeAsString_(attributes, "mz");
      if (!tmp_str.empty())
      {
        pos_[Peak2D::MZ] = asDouble_(tmp_str);
      }

      tmp_str = attributeAsString_(attributes, "it");
      if (!tmp_str.empty())
      {
        it_ = asDouble_(tmp_str);
      }

    }
    else if (tag == "element")
    {
      FeatureHandle act_index_tuple;
      UniqueIdInterface tmp_unique_id_interface;

      tmp_str = attributeAsString_(attributes, "map");
      if (!tmp_str.empty())
      {
        tmp_unique_id_interface.setUniqueId(tmp_str);
        UInt64 map_index = tmp_unique_id_interface.getUniqueId();

        tmp_str = attributeAsString_(attributes, "id");
        if (!tmp_str.empty())
        {
          tmp_unique_id_interface.setUniqueId(tmp_str);
          UInt64 unique_id = tmp_unique_id_interface.getUniqueId();

          act_index_tuple.setMapIndex(map_index);
          act_index_tuple.setUniqueId(unique_id);

          tmp_str = attributeAsString_(attributes, "rt");
          DPosition<2> pos;
          pos[0] = asDouble_(tmp_str);
          tmp_str = attributeAsString_(attributes, "mz");
          pos[1] = asDouble_(tmp_str);

          act_index_tuple.setPosition(pos);
          act_index_tuple.setIntensity(attributeAsDouble_(attributes, "it"));

          Int charge = 0;
          if (optionalAttributeAsInt_(charge, attributes, "charge"))
          {
            act_index_tuple.setCharge(charge);
          }

          act_cons_element_.insert(std::move(act_index_tuple));
        }
      }
      act_cons_element_.getPosition() = pos_;
      act_cons_element_.setIntensity(it_);
    }
    else if (tag == "consensusXML")
    {
      startProgress(0, 0, "loading consensusXML file");
      progress_ = 0;
      setProgress(++progress_);
      //check file version against schema version
      String file_version = "";
      optionalAttributeAsString_(file_version, attributes, "version");
      if (file_version.empty())
      {
        file_version = "1.0"; //default version is 1.0
      }
      if (file_version.toDouble() > version_.toDouble())
      {
        warning(LOAD, "The XML file (" + file_version + ") is newer than the parser (" + version_ + "). This might lead to undefined program behavior.");
      }
      // handle document id
      String document_id;
      if (optionalAttributeAsString_(document_id, attributes, "document_id"))
      {
        consensus_map_->setIdentifier(document_id);
      }
      // handle unique id
      String unique_id;
      if (optionalAttributeAsString_(unique_id, attributes, "id"))
      {
        consensus_map_->setUniqueId(unique_id);
      }
      // TODO The next four lines should be removed in OpenMS 1.7 or so!
      if (optionalAttributeAsString_(unique_id, attributes, "unique_id"))
      {
        consensus_map_->setUniqueId(unique_id);
      }
      //handle experiment type
      String experiment_type;
      if (optionalAttributeAsString_(experiment_type, attributes, "experiment_type"))
      {
        consensus_map_->setExperimentType(experiment_type);
      }
      last_meta_ = consensus_map_;
    }
    else if (tag == "userParam" || tag == "UserParam") // remain backwards compatible. Correct is "UserParam"
    {
      if (last_meta_ == nullptr)
      {
        fatalError(LOAD, String("Unexpected UserParam in tag '") + parent_tag + "'");
      }

      String name = attributeAsString_(attributes, "name");
      String type = attributeAsString_(attributes, "type");

      if (type == "int")
      {
        last_meta_->setMetaValue(name, attributeAsInt_(attributes, "value"));
      }
      else if (type == "float")
      {
        last_meta_->setMetaValue(name, attributeAsDouble_(attributes, "value"));
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
      else if (type == "string")
      {
        last_meta_->setMetaValue(name, (String) attributeAsString_(attributes, "value"));
      }
      else
      {
        fatalError(LOAD, String("Invalid UserParam type '") + type + "'");
      }
    }
    else if (tag == "IdentificationRun")
    {
      setProgress(++progress_);
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
      setProgress(++progress_);
      prot_hit_ = ProteinHit();
      String accession = attributeAsString_(attributes, "accession");
      prot_hit_.setAccession(accession);
      prot_hit_.setScore(attributeAsDouble_(attributes, "score"));

      // coverage
      if (double coverage; optionalAttributeAsDouble_(coverage, attributes, "coverage"))
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
      if (double thresh; optionalAttributeAsDouble_(thresh, attributes, "significance_threshold"))
      {
        pep_id_.setSignificanceThreshold(thresh);
      }

      //score orientation
      pep_id_.setHigherScoreBetter(asBool_(attributeAsString_(attributes, "higher_score_better")));

      //MZ
      if (double mz; optionalAttributeAsDouble_(mz, attributes, "MZ"))
      {
        pep_id_.setMZ(mz);
      }
      //RT
      if (double rt; optionalAttributeAsDouble_(rt, attributes, "RT"))
      {
        pep_id_.setRT(rt);
      }

      if (String ref; optionalAttributeAsString_(ref, attributes, "spectrum_reference"))
      {
        pep_id_.setSpectrumReference( ref);
      }

      last_meta_ = &pep_id_;
    }
    else if (tag == "PeptideHit")
    {
      setProgress(++progress_);
      pep_hit_ = PeptideHit();
      peptide_evidences_ = vector<PeptideEvidence>();
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
          accessions.push_back(std::move(accession_string));
        }

        for (vector<String>::const_iterator it = accessions.begin(); it != accessions.end(); ++it)
        {
          std::map<String, String>::const_iterator it2 = proteinid_to_accession_.find(*it);
          if (it2 != proteinid_to_accession_.end())
          {
            PeptideEvidence pe;
            pe.setProteinAccession(it2->second);
            peptide_evidences_.push_back(std::move(pe));
          }
          else
          {
            fatalError(LOAD, String("Invalid protein reference '") + *it + "'");
          }
        }
      }

      //aa_before
      String tmp; 
      std::vector<String> splitted;
      if (optionalAttributeAsString_(tmp, attributes, "aa_before"))
      {
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
      if (optionalAttributeAsString_(tmp, attributes, "aa_after"))
      {
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
      if (optionalAttributeAsString_(tmp, attributes, "start"))
      {
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
      if (optionalAttributeAsString_(tmp, attributes, "end"))
      {
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

      last_meta_ = &pep_hit_;
    }
    else if (tag == "dataProcessing")
    {
      setProgress(++progress_);
      DataProcessing tmp;
      tmp.setCompletionTime(asDateTime_(attributeAsString_(attributes, "completion_time")));
      consensus_map_->getDataProcessing().push_back(std::move(tmp));
      last_meta_ = &(consensus_map_->getDataProcessing().back());
    }
    else if (tag == "software" && parent_tag == "dataProcessing")
    {
      consensus_map_->getDataProcessing().back().getSoftware().setName(attributeAsString_(attributes, "name"));
      consensus_map_->getDataProcessing().back().getSoftware().setVersion(attributeAsString_(attributes, "version"));
    }
    else if (tag == "processingAction" && parent_tag == "dataProcessing")
    {
      String name = attributeAsString_(attributes, "name");
      for (Size i = 0; i < DataProcessing::SIZE_OF_PROCESSINGACTION; ++i)
      {
        if (name == DataProcessing::NamesOfProcessingAction[i])
        {
          consensus_map_->getDataProcessing().back().getProcessingActions().insert((DataProcessing::ProcessingAction) i);
        }
      }
    }
  }

  void ConsensusXMLHandler::writeTo(std::ostream& os)
  {
    const ConsensusMap& consensus_map = *(cconsensus_map_);

    startProgress(0, 0, "storing consensusXML file");
    progress_ = 0;
    setProgress(++progress_);

    setProgress(++progress_);
    os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    os << "<?xml-stylesheet type=\"text/xsl\" href=\"https://www.openms.de/xml-stylesheet/ConsensusXML.xsl\" ?>\n";

    setProgress(++progress_);
    os << "<consensusXML version=\"" << version_ << "\"";
    // file id
    if (!consensus_map.getIdentifier().empty())
    {
      os << " document_id=\"" << consensus_map.getIdentifier() << "\"";
    }
    // unique id
    if (consensus_map.hasValidUniqueId())
    {
      os << " id=\"cm_" << consensus_map.getUniqueId() << "\"";
    }
    if (!consensus_map.getExperimentType().empty())
    {
      os << " experiment_type=\"" << consensus_map.getExperimentType() << "\"";
    }
    os
      << " xsi:noNamespaceSchemaLocation=\"https://raw.githubusercontent.com/OpenMS/OpenMS/develop/share/OpenMS/SCHEMAS/ConsensusXML_1_7.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

    // user param
    writeUserParam_("UserParam", os, consensus_map, 1);
    setProgress(++progress_);

    // write data processing
    for (Size i = 0; i < consensus_map.getDataProcessing().size(); ++i)
    {
      const DataProcessing& processing = consensus_map.getDataProcessing()[i];
      os << "\t<dataProcessing completion_time=\"" << processing.getCompletionTime().getDate() << 'T' << processing.getCompletionTime().getTime() << "\">\n";
      os << "\t\t<software name=\"" << processing.getSoftware().getName() << "\" version=\"" << processing.getSoftware().getVersion() << "\" />\n";
      for (set<DataProcessing::ProcessingAction>::const_iterator it = processing.getProcessingActions().begin(); it != processing.getProcessingActions().end(); ++it)
      {
        os << "\t\t<processingAction name=\"" << DataProcessing::NamesOfProcessingAction[*it] << "\" />\n";
      }
      writeUserParam_("UserParam", os, processing, 2);
      os << "\t</dataProcessing>\n";
    }
    setProgress(++progress_);

    // write identification run
    UInt prot_count = 0;

    // throws if protIDs are not unique, i.e. PeptideIDs will be randomly assigned (bad!)
    checkUniqueIdentifiers_(consensus_map.getProteinIdentifications());

    for (UInt i = 0; i < consensus_map.getProteinIdentifications().size(); ++i)
    {
      setProgress(++progress_);
      const ProteinIdentification& current_prot_id = consensus_map.getProteinIdentifications()[i];
      os << "\t<IdentificationRun ";
      os << "id=\"PI_" << i << "\" ";
      identifier_id_[current_prot_id.getIdentifier()] = String("PI_") + i;
      os << "date=\"" << current_prot_id.getDateTime().getDate() << "T" << current_prot_id.getDateTime().getTime() << "\" ";
      os << "search_engine=\"" << writeXMLEscape(current_prot_id.getSearchEngine()) << "\" ";
      os << "search_engine_version=\"" << writeXMLEscape(current_prot_id.getSearchEngineVersion()) << "\">\n";

      //write search parameters
      const ProteinIdentification::SearchParameters& search_param = current_prot_id.getSearchParameters();
      os << "\t\t<SearchParameters " << "db=\"" << search_param.db << "\" " << "db_version=\"" << search_param.db_version << "\" " << "taxonomy=\""
         << search_param.taxonomy << "\" ";
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
      }
      for (Size j = 0; j != search_param.variable_modifications.size(); ++j)
      {
        os << "\t\t\t<VariableModification name=\"" << writeXMLEscape(search_param.variable_modifications[j]) << "\" />\n";
      }

      writeUserParam_("UserParam", os, search_param, 4);

      os << "\t\t</SearchParameters>\n";

      //write protein identifications
      os << "\t\t<ProteinIdentification";
      os << " score_type=\"" << writeXMLEscape(current_prot_id.getScoreType()) << "\"";
      os << " higher_score_better=\"" << (current_prot_id.isHigherScoreBetter() ? "true" : "false") << "\"";
      os << " significance_threshold=\"" << current_prot_id.getSignificanceThreshold() << "\">\n";

      //TODO @julianus @timo IMPLEMENT PROTEIN GROUP SUPPORT!!
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

      // add ProteinGroup info to metavalues (hack)
      MetaInfoInterface meta = current_prot_id;
      addProteinGroups_(meta, current_prot_id.getProteinGroups(),
                        "protein_group", accession_to_id_, current_prot_id.getIdentifier(), STORE);
      addProteinGroups_(meta, current_prot_id.getIndistinguishableProteins(),
                        "indistinguishable_proteins", accession_to_id_, current_prot_id.getIdentifier(), STORE);
      writeUserParam_("UserParam", os, meta, 3);
      os << "\t\t</ProteinIdentification>\n";
      os << "\t</IdentificationRun>\n";
    }

    //write unassigned peptide identifications
    for (UInt i = 0; i < consensus_map.getUnassignedPeptideIdentifications().size(); ++i)
    {
      writePeptideIdentification_(file_, os, consensus_map.getUnassignedPeptideIdentifications()[i], "UnassignedPeptideIdentification", 1);
    }

    //file descriptions
    const ConsensusMap::ColumnHeaders& description_vector = consensus_map.getColumnHeaders();
    os << "\t<mapList count=\"" << description_vector.size() << "\">\n";
    for (ConsensusMap::ColumnHeaders::const_iterator it = description_vector.begin(); it != description_vector.end(); ++it)
    {
      setProgress(++progress_);
      os << "\t\t<map id=\"" << it->first;
      os << "\" name=\"" << it->second.filename;
      if (UniqueIdInterface::isValid(it->second.unique_id))
      {
        os << "\" unique_id=\"" << it->second.unique_id;
      }
      os << "\" label=\"" << it->second.label;
      os << "\" size=\"" << it->second.size << "\">\n";
      writeUserParam_("UserParam", os, it->second, 3);
      os << "\t\t</map>\n";
    }
    os << "\t</mapList>\n";

    // write all consensus elements
    os << "\t<consensusElementList>\n";
    for (Size i = 0; i < consensus_map.size(); ++i)
    {
      setProgress(++progress_);
      // write a consensusElement
      const ConsensusFeature& elem = consensus_map[i];
      os << "\t\t<consensusElement id=\"e_" << elem.getUniqueId() << "\" quality=\"" << precisionWrapper(elem.getQuality()) << "\"";
      if (elem.getCharge() != 0)
      {
        os << " charge=\"" << elem.getCharge() << "\"";
      }
      os << ">\n";
      // write centroid
      os << "\t\t\t<centroid rt=\"" << precisionWrapper(elem.getRT()) << "\" mz=\"" << precisionWrapper(elem.getMZ()) << "\" it=\"" << precisionWrapper(
        elem.getIntensity()) << "\"/>\n";
      // write groupedElementList
      os << "\t\t\t<groupedElementList>\n";
      for (ConsensusFeature::HandleSetType::const_iterator it = elem.begin(); it != elem.end(); ++it)
      {
        os << "\t\t\t\t<element"
              " map=\"" << it->getMapIndex() << "\""
                                                " id=\"" << it->getUniqueId() << "\""
                                                                                 " rt=\"" << precisionWrapper(it->getRT()) << "\""
                                                                                                                              " mz=\"" << precisionWrapper(it->getMZ()) << "\""
                                                                                                                                                                           " it=\"" << precisionWrapper(it->getIntensity()) << "\"";
        if (it->getCharge() != 0)
        {
          os << " charge=\"" << it->getCharge() << "\"";
        }
        os << "/>\n";
      }
      os << "\t\t\t</groupedElementList>\n";

      // write PeptideIdentification
      for (UInt j = 0; j < elem.getPeptideIdentifications().size(); ++j)
      {
        writePeptideIdentification_(file_, os, elem.getPeptideIdentifications()[j], "PeptideIdentification", 3);
      }

      writeUserParam_("UserParam", os, elem, 3);
      os << "\t\t</consensusElement>\n";
    }
    os << "\t</consensusElementList>\n";

    os << "</consensusXML>\n";

    //Clear members
    identifier_id_.clear();
    accession_to_id_.clear();
    endProgress();
  }

  void ConsensusXMLHandler::writePeptideIdentification_(const String& filename, std::ostream& os, const PeptideIdentification& id, const String& tag_name,
                                                UInt indentation_level)
  {
    String indent = String(indentation_level, '\t');

    if (identifier_id_.find(id.getIdentifier()) == identifier_id_.end())
    {
      warning(STORE, String("Omitting peptide identification because of missing ProteinIdentification with identifier '") + id.getIdentifier()
              + "' while writing '" + filename + "'!");
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
      os << indent << "\t<PeptideHit";
      os << " score=\"" << id.getHits()[j].getScore() << "\"";
      os << " sequence=\"" << writeXMLEscape(id.getHits()[j].getSequence().toString()) << "\"";
      os << " charge=\"" << id.getHits()[j].getCharge() << "\"";

      vector<PeptideEvidence> pes = id.getHits()[j].getPeptideEvidences();

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

    // do not write "spectrum_reference" since it is written as attribute already
    MetaInfoInterface tmp = id;
    tmp.removeMetaValue("spectrum_reference");
    writeUserParam_("UserParam", os, tmp, indentation_level + 1);
    os << indent << "</" << tag_name << ">\n";
  }

  void ConsensusXMLHandler::addProteinGroups_(
      MetaInfoInterface& meta, const std::vector<ProteinIdentification::ProteinGroup>& groups,
      const String& group_name, const std::unordered_map<string, UInt>& accession_to_id, const String& runid,
      XMLHandler::ActionMode mode)
  {
    for (Size g = 0; g < groups.size(); ++g)
    {
      String name = group_name + "_" + String(g);
      if (meta.metaValueExists(name))
      {
        warning(mode, String("Metavalue '") + name + "' already exists. Overwriting...");
      }
      String accessions;
      for (StringList::const_iterator acc_it = groups[g].accessions.begin();
           acc_it != groups[g].accessions.end(); ++acc_it)
      {
        if (acc_it != groups[g].accessions.begin())
          accessions += ",";
        const auto pos = accession_to_id.find(runid + "_" + *acc_it);
        if (pos != accession_to_id.end())
        {
          accessions += "PH_" + String(pos->second);
        }
        else
        {
          fatalError(mode, String("Invalid protein reference '") + *acc_it + "'");
        }
      }
      String value = String(groups[g].probability) + "," + accessions;
      meta.setMetaValue(name, value);
    }
  }

  void ConsensusXMLHandler::getProteinGroups_(std::vector<ProteinIdentification::ProteinGroup>&
  groups, const String& group_name)
  {
    groups.clear();
    Size g_id = 0;
    String current_meta = group_name + "_" + String(g_id);
    StringList values;
    while (last_meta_->metaValueExists(current_meta)) // assumes groups have incremental g_IDs
    {
      // convert to proper ProteinGroup
      ProteinIdentification::ProteinGroup g;
      String(last_meta_->getMetaValue(current_meta)).split(',', values);
      if (values.size() < 2)
      {
        fatalError(LOAD, String("Invalid UserParam for ProteinGroups (not enough values)'"));
      }
      g.probability = values[0].toDouble();
      for (Size i_ind = 1; i_ind < values.size(); ++i_ind)
      {
        g.accessions.push_back(proteinid_to_accession_[values[i_ind]]);
      }
      groups.push_back(std::move(g));
      last_meta_->removeMetaValue(current_meta);
      current_meta = group_name + "_" + String(++g_id);
    }
  }

} // namespace OpenMS
