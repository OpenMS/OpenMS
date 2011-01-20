// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Clemens Groepl $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <fstream>

using namespace std;

namespace OpenMS
{
  ConsensusXMLFile::ConsensusXMLFile() :
    XMLHandler("", "1.4"), XMLFile("/SCHEMAS/ConsensusXML_1_4.xsd", "1.4"), ProgressLogger(), consensus_map_(0), act_cons_element_(), last_meta_(0)
  {
  }

  ConsensusXMLFile::~ConsensusXMLFile()
  {
  }

  PeakFileOptions&
  ConsensusXMLFile::getOptions()
  {
    return options_;
  }

  const PeakFileOptions&
  ConsensusXMLFile::getOptions() const
  {
    return options_;
  }

  void
  ConsensusXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
  {
    String tag = sm_.convert(qname);
    open_tags_.pop_back();

    if ( tag == "consensusElement" )
    {
      if ( (!options_.hasRTRange() || options_.getRTRange().encloses(act_cons_element_.getRT())) && (!options_.hasMZRange() || options_.getMZRange().encloses(
          act_cons_element_.getMZ())) && (!options_.hasIntensityRange() || options_.getIntensityRange().encloses(act_cons_element_.getIntensity())) )
      {
        consensus_map_->push_back(act_cons_element_);
        act_cons_element_.getPeptideIdentifications().clear();
      }
      last_meta_ = 0;
    }
    else if ( tag == "IdentificationRun" )
    {
      consensus_map_->getProteinIdentifications().push_back(prot_id_);
      prot_id_ = ProteinIdentification();
      last_meta_ = 0;
    }
    else if ( tag == "SearchParameters" )
    {
      prot_id_.setSearchParameters(search_param_);
    }
    else if ( tag == "ProteinHit" )
    {
      prot_id_.insertHit(prot_hit_);
      last_meta_ = &prot_id_;
    }
    else if ( tag == "PeptideIdentification" )
    {
      act_cons_element_.getPeptideIdentifications().push_back(pep_id_);
      pep_id_ = PeptideIdentification();
      last_meta_ = &act_cons_element_;
    }
    else if ( tag == "UnassignedPeptideIdentification" )
    {
      consensus_map_->getUnassignedPeptideIdentifications().push_back(pep_id_);
      pep_id_ = PeptideIdentification();
      last_meta_ = consensus_map_;
    }
    else if ( tag == "PeptideHit" )
    {
      pep_id_.insertHit(pep_hit_);
      last_meta_ = &pep_id_;
    }
    else if ( tag == "consensusXML" )
    {
      endProgress();
    }
  }

  void
  ConsensusXMLFile::characters(const XMLCh* const /*chars*/, const XMLSize_t /*length*/)
  {
  }

  void
  ConsensusXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  {
    String tag = sm_.convert(qname);
    String parent_tag;
    if ( open_tags_.size() != 0 )
      parent_tag = open_tags_.back();
    open_tags_.push_back(tag);

    String tmp_str;
    if ( tag == "map" )
    {
      setProgress(++progress_);
      last_map_ = attributeAsInt_(attributes, "id");
      last_meta_ = &consensus_map_->getFileDescriptions()[last_map_];
      consensus_map_->getFileDescriptions()[last_map_].filename = attributeAsString_(attributes, "name");
      String unique_id;
      if ( XMLHandler::optionalAttributeAsString_(unique_id, attributes, "unique_id") )
      {
        UniqueIdInterface tmp;
        tmp.setUniqueId(unique_id);
        consensus_map_->getFileDescriptions()[last_map_].unique_id = tmp.getUniqueId();
      }
      String label;
      if ( XMLHandler::optionalAttributeAsString_(label, attributes, "label") )
      {
        consensus_map_->getFileDescriptions()[last_map_].label = label;
      }
      UInt size;
      if ( XMLHandler::optionalAttributeAsUInt_(size, attributes, "size") )
      {
        consensus_map_->getFileDescriptions()[last_map_].size = size;
      }
    }
    else if ( tag == "consensusElement" )
    {
      setProgress(++progress_);
      act_cons_element_ = ConsensusFeature();
      last_meta_ = &act_cons_element_;
      // quality
      DoubleReal quality = 0.0;
      if ( optionalAttributeAsDouble_(quality, attributes, "quality") )
      {
        act_cons_element_.setQuality(quality);
      }
      // charge
      Int charge = 0;
      if ( optionalAttributeAsInt_(charge, attributes, "charge") )
      {
        act_cons_element_.setCharge(charge);
      }
      // unique id
      act_cons_element_.setUniqueId(attributeAsString_(attributes,"id"));
      last_meta_ = &act_cons_element_;
    }
    else if ( tag == "centroid" )
    {
      tmp_str = attributeAsString_(attributes, "rt");
      if ( tmp_str != "" )
      {
        pos_[Peak2D::RT] = asDouble_(tmp_str);
      }

      tmp_str = attributeAsString_(attributes, "mz");
      if ( tmp_str != "" )
      {
        pos_[Peak2D::MZ] = asDouble_(tmp_str);
      }

      tmp_str = attributeAsString_(attributes, "it");
      if ( tmp_str != "" )
      {
        it_ = asDouble_(tmp_str);
      }

    }
    else if ( tag == "element" )
    {
      FeatureHandle act_index_tuple;
      UniqueIdInterface tmp_unique_id_interface;

      tmp_str = attributeAsString_(attributes, "map");
      if ( tmp_str != "" )
      {
        tmp_unique_id_interface.setUniqueId(tmp_str);
        UInt64 map_index = tmp_unique_id_interface.getUniqueId();

        tmp_str = attributeAsString_(attributes, "id");
        if ( tmp_str != "" )
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
          if ( optionalAttributeAsInt_(charge, attributes, "charge") )
          {
            act_index_tuple.setCharge(charge);
          }

          act_cons_element_.insert(act_index_tuple);
        }
      }
      act_cons_element_.getPosition() = pos_;
      act_cons_element_.setIntensity(it_);
    }
    else if ( tag == "consensusXML" )
    {
      startProgress(0, 0, "loading consensusXML file");
      progress_ = 0;
      setProgress(++progress_);
      //check file version against schema version
      String file_version = "";
      optionalAttributeAsString_(file_version, attributes, "version");
      if ( file_version == "" )
        file_version = "1.0"; //default version is 1.0
      if ( file_version.toDouble() > version_.toDouble() )
      {
        warning(LOAD, "The XML file (" + file_version + ") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
      }
      // handle document id
      String document_id;
      if ( optionalAttributeAsString_(document_id, attributes, "document_id") )
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
      if ( optionalAttributeAsString_(experiment_type, attributes, "experiment_type") )
      {
        consensus_map_->setExperimentType(experiment_type);
      }
      last_meta_ = consensus_map_;
    }
    else if ( tag == "userParam" )
    {
      if ( last_meta_ == 0 )
      {
        fatalError(LOAD, String("Unexpected userParam in tag '") + parent_tag + "'");
      }

      String name = attributeAsString_(attributes, "name");
      String type = attributeAsString_(attributes, "type");

      if ( type == "int" )
      {
        last_meta_->setMetaValue(name, attributeAsInt_(attributes, "value"));
      }
      else if ( type == "float" )
      {
        last_meta_->setMetaValue(name, attributeAsDouble_(attributes, "value"));
      }
      else if ( type == "string" )
      {
        last_meta_->setMetaValue(name, (String) attributeAsString_(attributes, "value"));
      }
      else
      {
        fatalError(LOAD, String("Invalid userParam type '") + type + "'");
      }
    }
    else if ( tag == "IdentificationRun" )
    {
      setProgress(++progress_);
      prot_id_.setSearchEngine(attributeAsString_(attributes, "search_engine"));
      prot_id_.setSearchEngineVersion(attributeAsString_(attributes, "search_engine_version"));
      prot_id_.setDateTime(DateTime::fromString(String(attributeAsString_(attributes, "date")).toQString(), "yyyy-MM-ddThh:mm:ss"));
      //set identifier
      String identifier = prot_id_.getSearchEngine() + '_' + attributeAsString_(attributes, "date");
      prot_id_.setIdentifier(identifier);
      id_identifier_[attributeAsString_(attributes, "id")] = identifier;
    }
    else if ( tag == "SearchParameters" )
    {
      //load parameters
      search_param_.db = attributeAsString_(attributes, "db");
      search_param_.db_version = attributeAsString_(attributes, "db_version");
      optionalAttributeAsString_(search_param_.taxonomy, attributes, "taxonomy");
      search_param_.charges = attributeAsString_(attributes, "charges");
      optionalAttributeAsUInt_(search_param_.missed_cleavages, attributes, "missed_cleavages");
      search_param_.peak_mass_tolerance = attributeAsDouble_(attributes, "peak_mass_tolerance");
      search_param_.precursor_tolerance = attributeAsDouble_(attributes, "precursor_peak_tolerance");
      //mass type
      String mass_type = attributeAsString_(attributes, "mass_type");
      if ( mass_type == "monoisotopic" )
      {
        search_param_.mass_type = ProteinIdentification::MONOISOTOPIC;
      }
      else if ( mass_type == "average" )
      {
        search_param_.mass_type = ProteinIdentification::AVERAGE;
      }
      //enzyme
      String enzyme;
      optionalAttributeAsString_(enzyme, attributes, "enzyme");
      if ( enzyme == "trypsin" )
      {
        search_param_.enzyme = ProteinIdentification::TRYPSIN;
      }
      else if ( enzyme == "pepsin_a" )
      {
        search_param_.enzyme = ProteinIdentification::PEPSIN_A;
      }
      else if ( enzyme == "protease_k" )
      {
        search_param_.enzyme = ProteinIdentification::PROTEASE_K;
      }
      else if ( enzyme == "chymotrypsin" )
      {
        search_param_.enzyme = ProteinIdentification::CHYMOTRYPSIN;
      }
      else if ( enzyme == "no_enzyme" )
      {
        search_param_.enzyme = ProteinIdentification::NO_ENZYME;
      }
      else if ( enzyme == "unknown_enzyme" )
      {
        search_param_.enzyme = ProteinIdentification::UNKNOWN_ENZYME;
      }
      last_meta_ = &search_param_;
    }
    else if ( tag == "FixedModification" )
    {
      search_param_.fixed_modifications.push_back(attributeAsString_(attributes, "name"));
      //change this line as soon as there is a MetaInfoInterface for modifications (Andreas)
      last_meta_ = 0;
    }
    else if ( tag == "VariableModification" )
    {
      search_param_.variable_modifications.push_back(attributeAsString_(attributes, "name"));
      //change this line as soon as there is a MetaInfoInterface for modifications (Andreas)
      last_meta_ = 0;
    }
    else if ( tag == "ProteinIdentification" )
    {
      prot_id_.setScoreType(attributeAsString_(attributes, "score_type"));

      //optional significance threshold
      DoubleReal tmp = 0.0;
      optionalAttributeAsDouble_(tmp, attributes, "significance_threshold");
      if ( tmp != 0.0 )
      {
        prot_id_.setSignificanceThreshold(tmp);
      }

      //score orientation
      prot_id_.setHigherScoreBetter(asBool_(attributeAsString_(attributes, "higher_score_better")));

      last_meta_ = &prot_id_;
    }
    else if ( tag == "ProteinHit" )
    {
      setProgress(++progress_);
      prot_hit_ = ProteinHit();
      String accession = attributeAsString_(attributes, "accession");
      prot_hit_.setAccession(accession);
      prot_hit_.setScore(attributeAsDouble_(attributes, "score"));

      //sequence
      String tmp = "";
      optionalAttributeAsString_(tmp, attributes, "sequence");
      prot_hit_.setSequence(tmp);

      last_meta_ = &prot_hit_;

      //insert id and accession to map
      proteinid_to_accession_[attributeAsString_(attributes, "id")] = accession;
    }
    else if ( tag == "PeptideIdentification" || tag == "UnassignedPeptideIdentification" )
    {
      String id = attributeAsString_(attributes, "identification_run_ref");
      if ( !id_identifier_.has(id) )
      {
        warning(LOAD, String("Peptide identification without ProteinIdentification found (id: '") + id + "')!");
      }
      pep_id_.setIdentifier(id_identifier_[id]);

      pep_id_.setScoreType(attributeAsString_(attributes, "score_type"));

      //optional significance threshold
      DoubleReal tmp = 0.0;
      optionalAttributeAsDouble_(tmp, attributes, "significance_threshold");
      if ( tmp != 0.0 )
      {
        pep_id_.setSignificanceThreshold(tmp);
      }

      //score orientation
      pep_id_.setHigherScoreBetter(asBool_(attributeAsString_(attributes, "higher_score_better")));

      //MZ
      DoubleReal tmp2 = -numeric_limits<DoubleReal>::max();
      optionalAttributeAsDouble_(tmp2, attributes, "MZ");
      if ( tmp2 != -numeric_limits<DoubleReal>::max() )
      {
        pep_id_.setMetaValue("MZ", tmp2);
      }
      //RT
      tmp2 = -numeric_limits<DoubleReal>::max();
      optionalAttributeAsDouble_(tmp2, attributes, "RT");
      if ( tmp2 != -numeric_limits<DoubleReal>::max() )
      {
        pep_id_.setMetaValue("RT", tmp2);
      }
      Int tmp3 = -numeric_limits<Int>::max();
      optionalAttributeAsInt_(tmp3, attributes, "spectrum_reference");
      if ( tmp3 != -numeric_limits<Int>::max() )
      {
        pep_id_.setMetaValue("spectrum_reference", tmp3);
      }

      last_meta_ = &pep_id_;
    }
    else if ( tag == "PeptideHit" )
    {
      setProgress(++progress_);
      pep_hit_ = PeptideHit();
      pep_hit_.setCharge(attributeAsInt_(attributes, "charge"));
      pep_hit_.setScore(attributeAsDouble_(attributes, "score"));
      pep_hit_.setSequence(attributeAsString_(attributes, "sequence"));

      //aa_before
      String tmp = "";
      optionalAttributeAsString_(tmp, attributes, "aa_before");
      if ( !tmp.empty() )
      {
        pep_hit_.setAABefore(tmp[0]);
      }
      //aa_after
      tmp = "";
      optionalAttributeAsString_(tmp, attributes, "aa_after");
      if ( !tmp.empty() )
      {
        pep_hit_.setAAAfter(tmp[0]);
      }

      //parse optional protein ids to determine accessions
      const XMLCh* refs = attributes.getValue(sm_.convert("protein_refs"));
      if ( refs != 0 )
      {
        String accession_string = sm_.convert(refs);
        accession_string.trim();
        vector<String> accessions;
        accession_string.split(' ', accessions);
        if ( accession_string != "" && accessions.size() == 0 )
        {
          accessions.push_back(accession_string);
        }
        for ( vector<String>::const_iterator it = accessions.begin(); it != accessions.end(); ++it )
        {
          Map<String, String>::const_iterator it2 = proteinid_to_accession_.find(*it);
          if ( it2 != proteinid_to_accession_.end() )
          {
            pep_hit_.addProteinAccession(it2->second);
          }
          else
          {
            fatalError(LOAD, String("Invalid protein reference '") + *it + "'");
          }
        }
      }
      last_meta_ = &pep_hit_;
    }
    else if ( tag == "dataProcessing" )
    {
      setProgress(++progress_);
      DataProcessing tmp;
      tmp.setCompletionTime(asDateTime_(attributeAsString_(attributes, "completion_time")));
      consensus_map_->getDataProcessing().push_back(tmp);
      last_meta_ = &(consensus_map_->getDataProcessing().back());
    }
    else if ( tag == "software" && parent_tag == "dataProcessing" )
    {
      consensus_map_->getDataProcessing().back().getSoftware().setName(attributeAsString_(attributes, "name"));
      consensus_map_->getDataProcessing().back().getSoftware().setVersion(attributeAsString_(attributes, "version"));
    }
    else if ( tag == "processingAction" && parent_tag == "dataProcessing" )
    {
      String name = attributeAsString_(attributes, "name");
      for ( Size i = 0; i < DataProcessing::SIZE_OF_PROCESSINGACTION; ++i )
      {
        if ( name == DataProcessing::NamesOfProcessingAction[i] )
        {
          consensus_map_->getDataProcessing().back().getProcessingActions().insert((DataProcessing::ProcessingAction) i);
        }
      }
    }
  }

  void
  ConsensusXMLFile::store(const String& filename, const ConsensusMap& consensus_map)
  {
    startProgress(0, 0, "storing consensusXML file");
    progress_ = 0;
    setProgress(++progress_);

    if ( Size invalid_unique_ids = consensus_map.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId) )
    {
      // TODO Take care *outside* that this does not happen.
      // We can detect this here but it is too late to fix the problem;
      // there is no straightforward action to be taken in all cases.
      // Note also that we are given a const reference.
      LOG_INFO << String("ConsensusXMLFile::store():  found ") + invalid_unique_ids + " invalid unique ids" << std::endl;
    }

    // This will throw if the unique ids are not unique,
    // so we never create bad files in this respect.
    try
    {
      consensus_map.updateUniqueIdToIndex();
    }
    catch ( Exception::Postcondition& e )
    {
      LOG_FATAL_ERROR << e.getName() << ' ' << e.getMessage() << std::endl;
      throw;
    }

    //open stream
    ofstream os(filename.c_str());
    if ( !os )
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
    }

    os.precision(writtenDigits<DoubleReal> ());

    setProgress(++progress_);
    os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    //add XSLT file if it can be found
    try
    {
      String xslt_file = File::find("XSL/ConsensusXML.xsl");
      os << "<?xml-stylesheet type=\"text/xsl\" href=\"file:///" << xslt_file << "\"?>\n";
    }
    catch(Exception::FileNotFound&)
    {
    }

    setProgress(++progress_);
    os << "<consensusXML version=\"" << version_ << "\"";
    // file id
    if ( consensus_map.getIdentifier() != "" )
    {
      os << " document_id=\"" << consensus_map.getIdentifier() << "\"";
    }
    // unique id
    if (consensus_map.hasValidUniqueId())
    {
      os << " id=\"cm_" << consensus_map.getUniqueId() << "\"";
    }
    if ( consensus_map.getExperimentType() != "" )
    {
      os << " experiment_type=\"" << consensus_map.getExperimentType() << "\"";
    }
    os
        << " xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/ConsensusXML_1_4.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

    //user param
    writeUserParam_("userParam", os, consensus_map, 1);
    setProgress(++progress_);

    //write data processing
    for ( Size i = 0; i < consensus_map.getDataProcessing().size(); ++i )
    {
      const DataProcessing& processing = consensus_map.getDataProcessing()[i];
      os << "\t<dataProcessing completion_time=\"" << processing.getCompletionTime().getDate() << 'T' << processing.getCompletionTime().getTime() << "\">\n";
      os << "\t\t<software name=\"" << processing.getSoftware().getName() << "\" version=\"" << processing.getSoftware().getVersion() << "\" />\n";
      for ( set<DataProcessing::ProcessingAction>::const_iterator it = processing.getProcessingActions().begin(); it != processing.getProcessingActions().end(); ++it )
      {
        os << "\t\t<processingAction name=\"" << DataProcessing::NamesOfProcessingAction[*it] << "\" />\n";
      }
      writeUserParam_("userParam", os, processing, 2);
      os << "\t</dataProcessing>\n";
    }
    setProgress(++progress_);

    // write identification run
    UInt prot_count = 0;

    for ( UInt i = 0; i < consensus_map.getProteinIdentifications().size(); ++i )
    {
      setProgress(++progress_);
      const ProteinIdentification& current_prot_id = consensus_map.getProteinIdentifications()[i];
      os << "\t<IdentificationRun ";
      os << "id=\"PI_" << i << "\" ";
      identifier_id_[current_prot_id.getIdentifier()] = String("PI_") + i;
      os << "date=\"" << current_prot_id.getDateTime().getDate() << "T" << current_prot_id.getDateTime().getTime() << "\" ";
      os << "search_engine=\"" << current_prot_id.getSearchEngine() << "\" ";
      os << "search_engine_version=\"" << current_prot_id.getSearchEngineVersion() << "\">\n";

      //write search parameters
      const ProteinIdentification::SearchParameters& search_param = current_prot_id.getSearchParameters();
      os << "\t\t<SearchParameters " << "db=\"" << search_param.db << "\" " << "db_version=\"" << search_param.db_version << "\" " << "taxonomy=\""
          << search_param.taxonomy << "\" ";
      if ( search_param.mass_type == ProteinIdentification::MONOISOTOPIC )
      {
        os << "mass_type=\"monoisotopic\" ";
      }
      else if ( search_param.mass_type == ProteinIdentification::AVERAGE )
      {
        os << "mass_type=\"average\" ";
      }
      os << "charges=\"" << search_param.charges << "\" ";
      if ( search_param.enzyme == ProteinIdentification::TRYPSIN )
      {
        os << "enzyme=\"trypsin\" ";
      }
      if ( search_param.enzyme == ProteinIdentification::PEPSIN_A )
      {
        os << "enzyme=\"pepsin_a\" ";
      }
      if ( search_param.enzyme == ProteinIdentification::PROTEASE_K )
      {
        os << "enzyme=\"protease_k\" ";
      }
      if ( search_param.enzyme == ProteinIdentification::CHYMOTRYPSIN )
      {
        os << "enzyme=\"chymotrypsin\" ";
      }
      else if ( search_param.enzyme == ProteinIdentification::NO_ENZYME )
      {
        os << "enzyme=\"no_enzyme\" ";
      }
      else if ( search_param.enzyme == ProteinIdentification::UNKNOWN_ENZYME )
      {
        os << "enzyme=\"unknown_enzyme\" ";
      }
      os << "missed_cleavages=\"" << search_param.missed_cleavages << "\" " << "precursor_peak_tolerance=\"" << search_param.precursor_tolerance << "\" "
          << "peak_mass_tolerance=\"" << search_param.peak_mass_tolerance << "\" " << ">\n";

      //modifications
      for ( Size j = 0; j != search_param.fixed_modifications.size(); ++j )
      {
        os << "\t\t\t<FixedModification name=\"" << search_param.fixed_modifications[j] << "\" />\n";
        //Add MetaInfo, when modifications has it (Andreas)
      }
      for ( Size j = 0; j != search_param.variable_modifications.size(); ++j )
      {
        os << "\t\t\t<VariableModification name=\"" << search_param.variable_modifications[j] << "\" />\n";
        //Add MetaInfo, when modifications has it (Andreas)
      }

      writeUserParam_("UserParam", os, search_param, 4);

      os << "\t\t</SearchParameters>\n";

      //write protein identifications
      os << "\t\t<ProteinIdentification";
      os << " score_type=\"" << current_prot_id.getScoreType() << "\"";
      os << " higher_score_better=\"" << (current_prot_id.isHigherScoreBetter() ? "true" : "false") << "\"";
      os << " significance_threshold=\"" << current_prot_id.getSignificanceThreshold() << "\">\n";

      // write protein hits
      for ( Size j = 0; j < current_prot_id.getHits().size(); ++j )
      {
        os << "\t\t\t<ProteinHit";

        // prot_count
        os << " id=\"PH_" << prot_count << "\"";
        accession_to_id_[current_prot_id.getIdentifier() + "_" + current_prot_id.getHits()[j].getAccession()] = prot_count;
        ++prot_count;

        os << " accession=\"" << current_prot_id.getHits()[j].getAccession() << "\"";
        os << " score=\"" << current_prot_id.getHits()[j].getScore() << "\"";
        os << " sequence=\"" << current_prot_id.getHits()[j].getSequence() << "\">\n";

        writeUserParam_("userParam", os, current_prot_id.getHits()[j], 4);

        os << "\t\t\t</ProteinHit>\n";
      }

      writeUserParam_("userParam", os, current_prot_id, 3);
      os << "\t\t</ProteinIdentification>\n";
      os << "\t</IdentificationRun>\n";
    }

    //write unassigned peptide identifications
    for ( UInt i = 0; i < consensus_map.getUnassignedPeptideIdentifications().size(); ++i )
    {
      writePeptideIdentification_(filename, os, consensus_map.getUnassignedPeptideIdentifications()[i], "UnassignedPeptideIdentification", 1);
    }

    //file descriptions
    const ConsensusMap::FileDescriptions& description_vector = consensus_map.getFileDescriptions();
    os << "\t<mapList count=\"" << description_vector.size() << "\">\n";
    for ( ConsensusMap::FileDescriptions::const_iterator it = description_vector.begin(); it != description_vector.end(); ++it )
    {
      setProgress(++progress_);
      os << "\t\t<map id=\"" << it->first;
      os << "\" name=\"" << it->second.filename;
      if ( UniqueIdInterface::isValid(it->second.unique_id) )
      {
        os << "\" unique_id=\"" << it->second.unique_id;
      }
      os << "\" label=\"" << it->second.label;
      os << "\" size=\"" << it->second.size << "\">\n";
      writeUserParam_("userParam", os, it->second, 3);
      os << "\t\t</map>\n";
    }
    os << "\t</mapList>\n";

    // write all consensus elements
    os << "\t<consensusElementList>\n";
    for ( Size i = 0; i < consensus_map.size(); ++i )
    {
      setProgress(++progress_);
      // write a consensusElement
      const ConsensusFeature& elem = consensus_map[i];
      os << "\t\t<consensusElement id=\"e_" << elem.getUniqueId() << "\" quality=\"" << precisionWrapper(elem.getQuality()) << "\"";
      if ( elem.getCharge() != 0 )
      {
        os << " charge=\"" << elem.getCharge() << "\"";
      }
      os << ">\n";
      // write centroid
      os << "\t\t\t<centroid rt=\"" << precisionWrapper(elem.getRT()) << "\" mz=\"" << precisionWrapper(elem.getMZ()) << "\" it=\"" << precisionWrapper(
          elem.getIntensity()) << "\"/>\n";
      // write groupedElementList
      os << "\t\t\t<groupedElementList>\n";
      for ( ConsensusFeature::HandleSetType::const_iterator it = elem.begin(); it != elem.end(); ++it )
      {
        os << "\t\t\t\t<element"
          " map=\"" << it->getMapIndex() << "\""
          " id=\"" << it->getUniqueId() << "\""
          " rt=\"" << precisionWrapper(it->getRT()) << "\""
          " mz=\"" << precisionWrapper(it->getMZ()) << "\""
          " it=\"" << precisionWrapper(it->getIntensity()) << "\"";
        if ( it->getCharge() != 0 )
        {
          os << " charge=\"" << it->getCharge() << "\"";
        }
        os << "/>\n";
      }
      os << "\t\t\t</groupedElementList>\n";

      // write PeptideIdentification
      for ( UInt i = 0; i < elem.getPeptideIdentifications().size(); ++i )
      {
        writePeptideIdentification_(filename, os, elem.getPeptideIdentifications()[i], "PeptideIdentification", 3);
      }

      writeUserParam_("userParam", os, elem, 3);
      os << "\t\t</consensusElement>\n";
    }
    os << "\t</consensusElementList>\n";

    os << "</consensusXML>\n";
    ;

    //Clear members
    identifier_id_.clear();
    accession_to_id_.clear();
    endProgress();
  }

  void
  ConsensusXMLFile::load(const String& filename, ConsensusMap& map)
  {
    //Filename for error messages in XMLHandler
    file_ = filename;

    map.clear(true); // clear map
    consensus_map_ = &map;

    //set DocumentIdentifier
    consensus_map_->setLoadedFileType(file_);
    consensus_map_->setLoadedFilePath(file_);

    parse_(filename, this);

    //reset members
    consensus_map_ = 0;
    act_cons_element_ = ConsensusFeature();
    pos_.clear();
    it_ = 0;
    last_map_ = 0;
    last_meta_ = 0;
    prot_id_ = ProteinIdentification();
    pep_id_ = PeptideIdentification();
    prot_hit_ = ProteinHit();
    pep_hit_ = PeptideHit();
    proteinid_to_accession_.clear();
    accession_to_id_.clear();
    identifier_id_.clear();
    id_identifier_.clear();
    search_param_ = ProteinIdentification::SearchParameters();
    progress_ = 0;
  }

  void
  ConsensusXMLFile::writePeptideIdentification_(const String& filename, std::ostream& os, const PeptideIdentification& id, const String& tag_name,
      UInt indentation_level)
  {
    String indent = String(indentation_level, '\t');

    if ( !identifier_id_.has(id.getIdentifier()) )
    {
      warning(STORE, String("Omitting peptide identification because of missing ProteinIdentification with identifier '") + id.getIdentifier()
          + "' while writing '" + filename + "'!");
      return;
    }
    os << indent << "<" << tag_name << " ";
    os << "identification_run_ref=\"" << identifier_id_[id.getIdentifier()] << "\" ";
    os << "score_type=\"" << id.getScoreType() << "\" ";
    os << "higher_score_better=\"" << (id.isHigherScoreBetter() ? "true" : "false") << "\" ";
    os << "significance_threshold=\"" << id.getSignificanceThreshold() << "\" ";
    //mz
    DataValue dv = id.getMetaValue("MZ");
    if ( dv != DataValue::EMPTY )
    {
      os << "MZ=\"" << dv.toString() << "\" ";
    }
    // rt
    dv = id.getMetaValue("RT");
    if ( dv != DataValue::EMPTY )
    {
      os << "RT=\"" << dv.toString() << "\" ";
    }
    // spectrum_reference
    dv = id.getMetaValue("spectrum_reference");
    if ( dv != DataValue::EMPTY )
    {
      os << "spectrum_reference=\"" << dv.toString() << "\" ";
    }
    os << ">\n";

    // write peptide hits
    for ( Size j = 0; j < id.getHits().size(); ++j )
    {
      os << indent << "\t<PeptideHit";
      os << " score=\"" << id.getHits()[j].getScore() << "\"";
      os << " sequence=\"" << id.getHits()[j].getSequence() << "\"";
      os << " charge=\"" << id.getHits()[j].getCharge() << "\"";
      if ( id.getHits()[j].getAABefore() != ' ' )
      {
        os << " aa_before=\"" << id.getHits()[j].getAABefore() << "\"";
      }
      if ( id.getHits()[j].getAAAfter() != ' ' )
      {
        os << " aa_after=\"" << id.getHits()[j].getAAAfter() << "\"";
      }
      if ( id.getHits()[j].getProteinAccessions().size() != 0 )
      {
        String accs = "";
        for ( Size m = 0; m < id.getHits()[j].getProteinAccessions().size(); ++m )
        {
          if ( m )
            accs += " ";
          accs += "PH_";
          accs += String(accession_to_id_[id.getIdentifier() + "_" + id.getHits()[j].getProteinAccessions()[m]]);
        }
        os << " protein_refs=\"" << accs << "\"";
      }
      os << ">\n";
      writeUserParam_("userParam", os, id.getHits()[j], indentation_level + 2);
      os << indent << "\t</PeptideHit>\n";
    }

    //do not write "RT", "MZ" and "spectrum_reference" as they are written as attributes already
    MetaInfoInterface tmp = id;
    tmp.removeMetaValue("RT");
    tmp.removeMetaValue("MZ");
    tmp.removeMetaValue("spectrum_reference");
    writeUserParam_("userParam", os, tmp, indentation_level + 1);
    os << indent << "</" << tag_name << ">\n";
  }

}// namespace OpenMS

