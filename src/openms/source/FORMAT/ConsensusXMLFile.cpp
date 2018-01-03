// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm, Mathias Walzer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
#include <OpenMS/SYSTEM/File.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/PrecisionWrapper.h>
#include <OpenMS/METADATA/DataProcessing.h>
#include <OpenMS/CONCEPT/UniqueIdGenerator.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <fstream>

using namespace std;

namespace OpenMS
{
  ConsensusXMLFile::ConsensusXMLFile() :
    XMLHandler("", "1.7"),
    XMLFile("/SCHEMAS/ConsensusXML_1_7.xsd", "1.7"),
    ProgressLogger(),
    consensus_map_(nullptr),
    act_cons_element_(),
    last_meta_(nullptr)
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
      consensus_map_->getProteinIdentifications().push_back(prot_id_);
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
      act_cons_element_.getPeptideIdentifications().push_back(pep_id_);
      pep_id_ = PeptideIdentification();
      last_meta_ = &act_cons_element_;
    }
    else if (tag == "UnassignedPeptideIdentification")
    {
      consensus_map_->getUnassignedPeptideIdentifications().push_back(pep_id_);
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

  void
  ConsensusXMLFile::characters(const XMLCh* const /*chars*/, const XMLSize_t /*length*/)
  {
  }

  void
  ConsensusXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
  {
    String tag = sm_.convert(qname);
    String parent_tag;
    if (!open_tags_.empty())
    {
      parent_tag = open_tags_.back();
    }
    open_tags_.push_back(tag);

    String tmp_str;
    if (tag == "map")
    {
      setProgress(++progress_);
      Size last_map = attributeAsInt_(attributes, "id");
      last_meta_ = &consensus_map_->getFileDescriptions()[last_map];
      consensus_map_->getFileDescriptions()[last_map].filename = attributeAsString_(attributes, "name");
      String unique_id;
      if (XMLHandler::optionalAttributeAsString_(unique_id, attributes, "unique_id"))
      {
        UniqueIdInterface tmp;
        tmp.setUniqueId(unique_id);
        consensus_map_->getFileDescriptions()[last_map].unique_id = tmp.getUniqueId();
      }
      String label;
      if (XMLHandler::optionalAttributeAsString_(label, attributes, "label"))
      {
        consensus_map_->getFileDescriptions()[last_map].label = label;
      }
      UInt size;
      if (XMLHandler::optionalAttributeAsUInt_(size, attributes, "size"))
      {
        consensus_map_->getFileDescriptions()[last_map].size = size;
      }
    }
    else if (tag == "consensusElement")
    {
      setProgress(++progress_);
      act_cons_element_ = ConsensusFeature();
      last_meta_ = &act_cons_element_;
      // quality
      double quality = 0.0;
      if (optionalAttributeAsDouble_(quality, attributes, "quality"))
      {
        act_cons_element_.setQuality(quality);
      }
      // charge
      Int charge = 0;
      if (optionalAttributeAsInt_(charge, attributes, "charge"))
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
      if (tmp_str != "")
      {
        pos_[Peak2D::RT] = asDouble_(tmp_str);
      }

      tmp_str = attributeAsString_(attributes, "mz");
      if (tmp_str != "")
      {
        pos_[Peak2D::MZ] = asDouble_(tmp_str);
      }

      tmp_str = attributeAsString_(attributes, "it");
      if (tmp_str != "")
      {
        it_ = asDouble_(tmp_str);
      }

    }
    else if (tag == "element")
    {
      FeatureHandle act_index_tuple;
      UniqueIdInterface tmp_unique_id_interface;

      tmp_str = attributeAsString_(attributes, "map");
      if (tmp_str != "")
      {
        tmp_unique_id_interface.setUniqueId(tmp_str);
        UInt64 map_index = tmp_unique_id_interface.getUniqueId();

        tmp_str = attributeAsString_(attributes, "id");
        if (tmp_str != "")
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

          act_cons_element_.insert(act_index_tuple);
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
      if (file_version == "")
        file_version = "1.0"; //default version is 1.0
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
      prot_id_.setDateTime(DateTime::fromString(String(attributeAsString_(attributes, "date")).toQString(), "yyyy-MM-ddThh:mm:ss"));
      //set identifier
      String identifier = prot_id_.getSearchEngine() + '_' + attributeAsString_(attributes, "date");
      String id = attributeAsString_(attributes, "id");

      if (!id_identifier_.has(id))
      {
        prot_id_.setIdentifier(identifier);
        id_identifier_[id] = identifier;
      }
      else
      {
        warning(LOAD, "Non-unique identifier for IdentificationRun encountered '" + identifier + "'. Generating a unique one.");
        UInt64 uid = UniqueIdGenerator::getUniqueId();
        identifier = identifier + String(uid);
        prot_id_.setIdentifier(identifier);
        id_identifier_[id] = identifier;
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
      if (!id_identifier_.has(id))
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
        pep_id_.setMetaValue("spectrum_reference", tmp3);
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
        if (accession_string != "" && accessions.empty())
        {
          accessions.push_back(accession_string);
        }

        for (vector<String>::const_iterator it = accessions.begin(); it != accessions.end(); ++it)
        {
          Map<String, String>::const_iterator it2 = proteinid_to_accession_.find(*it);
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
            peptide_evidences_.push_back(PeptideEvidence());
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
            peptide_evidences_.push_back(PeptideEvidence());
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
            peptide_evidences_.push_back(PeptideEvidence());
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
            peptide_evidences_.push_back(PeptideEvidence());
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
      consensus_map_->getDataProcessing().push_back(tmp);
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

  void
  ConsensusXMLFile::store(const String& filename, const ConsensusMap& consensus_map)
  {
    if (!FileHandler::hasValidExtension(filename, FileTypes::CONSENSUSXML))
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename, "invalid file extension, expected '" + FileTypes::typeToName(FileTypes::CONSENSUSXML) + "'");
    }

    if (!consensus_map.isMapConsistent(&LOG_WARN))
    {
      // Currently it is possible that FeatureLinkerUnlabeledQT triggers this exception
      // throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The ConsensusXML file contains invalid maps or references thereof. No data was written! Please fix the file or notify the maintainer of this tool if you did not provide a consensusXML file!");
      std::cerr << "The ConsensusXML file contains invalid maps or references thereof. Please fix the file or notify the maintainer of this tool if you did not provide a consensusXML file! Note that this warning will be a fatal error in the next version of OpenMS!" << std::endl;
    }

    startProgress(0, 0, "storing consensusXML file");
    progress_ = 0;
    setProgress(++progress_);

    if (Size invalid_unique_ids = consensus_map.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId))
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
    catch (Exception::Postcondition& e)
    {
      LOG_FATAL_ERROR << e.getName() << ' ' << e.getMessage() << std::endl;
      throw;
    }

    //open stream
    ofstream os(filename.c_str());
    if (!os)
    {
      throw Exception::UnableToCreateFile(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, filename);
    }

    os.precision(writtenDigits<double>(0.0));

    setProgress(++progress_);
    os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n";
    //add XSLT file if it can be found
    try
    {
      String xslt_file = File::find("XSL/ConsensusXML.xsl");
      os << "<?xml-stylesheet type=\"text/xsl\" href=\"file:///" << xslt_file << "\"?>\n";
    }
    catch (Exception::FileNotFound&)
    {
    }

    setProgress(++progress_);
    os << "<consensusXML version=\"" << version_ << "\"";
    // file id
    if (consensus_map.getIdentifier() != "")
    {
      os << " document_id=\"" << consensus_map.getIdentifier() << "\"";
    }
    // unique id
    if (consensus_map.hasValidUniqueId())
    {
      os << " id=\"cm_" << consensus_map.getUniqueId() << "\"";
    }
    if (consensus_map.getExperimentType() != "")
    {
      os << " experiment_type=\"" << consensus_map.getExperimentType() << "\"";
    }
    os
      << " xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/ConsensusXML_1_7.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

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
    for (UInt i = 0; i < consensus_map.getUnassignedPeptideIdentifications().size(); ++i)
    {
      writePeptideIdentification_(filename, os, consensus_map.getUnassignedPeptideIdentifications()[i], "UnassignedPeptideIdentification", 1);
    }

    //file descriptions
    const ConsensusMap::FileDescriptions& description_vector = consensus_map.getFileDescriptions();
    os << "\t<mapList count=\"" << description_vector.size() << "\">\n";
    for (ConsensusMap::FileDescriptions::const_iterator it = description_vector.begin(); it != description_vector.end(); ++it)
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
        writePeptideIdentification_(filename, os, elem.getPeptideIdentifications()[j], "PeptideIdentification", 3);
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

    if (!map.isMapConsistent(&LOG_WARN)) // a warning is printed to LOG_WARN during isMapConsistent()
    {
      // don't throw exception for now, since this would prevent us from reading old files...
      // throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION, "The ConsensusXML file contains invalid maps or references thereof. Please fix the file!");

    }

    //reset members
    consensus_map_ = nullptr;
    act_cons_element_ = ConsensusFeature();
    pos_.clear();
    it_ = 0;
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
    progress_ = 0;
    map.updateRanges();
  }

  void
  ConsensusXMLFile::writePeptideIdentification_(const String& filename, std::ostream& os, const PeptideIdentification& id, const String& tag_name,
                                                UInt indentation_level)
  {
    String indent = String(indentation_level, '\t');

    if (!identifier_id_.has(id.getIdentifier()))
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

      os << IdXMLFile::createFlankingAAXMLString_(pes);
      os << IdXMLFile::createPositionXMLString_(pes);

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

} // namespace OpenMS
