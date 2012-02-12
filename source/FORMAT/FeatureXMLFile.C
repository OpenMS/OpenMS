// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow, Clemens Groepl $
// $Authors: Marc Sturm, Chris Bielow, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/CONCEPT/LogStream.h>

#include <fstream>

using namespace std;

namespace OpenMS
{
	FeatureXMLFile::FeatureXMLFile()
		: Internal::XMLHandler("","1.6"),
			Internal::XMLFile("/SCHEMAS/FeatureXML_1_6.xsd","1.6"),
		 	map_(0),
		 	in_description_(false),
			subordinate_feature_level_(0),
			last_meta_(0)
	{
	}

	FeatureXMLFile::~FeatureXMLFile()
	{
	}

	void FeatureXMLFile::load(String filename, FeatureMap<>& feature_map)
	{
  	//Filename for error messages in XMLHandler
  	file_ = filename;

  	feature_map.clear(true);
  	map_ = &feature_map;

		//set DocumentIdentifier
		map_->setLoadedFileType(file_);
		map_->setLoadedFilePath(file_);

		parse_(filename, this);

		// !!! Hack: set feature FWHM from meta info entries as 
    // long as featureXML doesn't support a width entry.
    // See also hack in BaseFeature::setWidth().
		for (FeatureMap<>::Iterator it = map_->begin(); it != map_->end(); ++it)
		{
			if (it->metaValueExists("FWHM"))
			{
				it->setWidth((double)it->getMetaValue("FWHM"));
			}
		}

    // reset members
		current_feature_ = 0;
    map_ = 0;
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
		
		return;
	}

	void FeatureXMLFile::store(String filename, const FeatureMap<>& feature_map)
	{
  	//open stream
		ofstream os(filename.c_str());
		if (!os)
		{
			throw Exception::UnableToCreateFile(__FILE__, __LINE__, __PRETTY_FUNCTION__, filename);
		}

		if ( Size invalid_unique_ids = feature_map.applyMemberFunction(&UniqueIdInterface::hasInvalidUniqueId) )
		{

		  // TODO Take care *outside* that this does not happen.
		  // We can detect this here but it is too late to fix the problem;
		  // there is no straightforward action to be taken in all cases.
		  // Note also that we are given a const reference.
		  LOG_INFO << String("FeatureXMLFile::store():  found ") + invalid_unique_ids + " invalid unique ids" << std::endl;
		}

		// This will throw if the unique ids are not unique,
		// so we never create bad files in this respect.
		try
		{
		  feature_map.updateUniqueIdToIndex();
		}
		catch ( Exception::Postcondition& e )
		{
		  LOG_FATAL_ERROR << e.getName() << ' ' << e.getMessage() << std::endl;
		  throw;
		}

		os.precision(writtenDigits<DoubleReal>());

		os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\"?>\n"
			 << "<featureMap version=\"" << version_ << "\"";
		// file id
		if (feature_map.getIdentifier()!="")
		{
			os << " document_id=\"" << feature_map.getIdentifier() << "\"";
		}
		// unique id
    if (feature_map.hasValidUniqueId())
    {
      os << " id=\"fm_" << feature_map.getUniqueId() << "\"";
    }
		os << " xsi:noNamespaceSchemaLocation=\"http://open-ms.sourceforge.net/schemas/FeatureXML_1_4.xsd\" xmlns:xsi=\"http://www.w3.org/2001/XMLSchema-instance\">\n";

		//write data processing
		for (Size i=0; i< feature_map.getDataProcessing().size(); ++i)
		{
			const DataProcessing& processing = feature_map.getDataProcessing()[i];
			os << "\t<dataProcessing completion_time=\"" << processing.getCompletionTime().getDate() << 'T' << processing.getCompletionTime().getTime() << "\">\n";
			os << "\t\t<software name=\"" << processing.getSoftware().getName() << "\" version=\"" << processing.getSoftware().getVersion() << "\" />\n";
			for (set<DataProcessing::ProcessingAction>::const_iterator it = processing.getProcessingActions().begin(); it!=processing.getProcessingActions().end(); ++it)
			{
				os << "\t\t<processingAction name=\"" << DataProcessing::NamesOfProcessingAction[*it] << "\" />\n";
			}
			writeUserParam_ ("userParam", os, processing, 2);
			os << "\t</dataProcessing>\n";
		}

		// write identification runs
		Size prot_count = 0;
		for ( Size i = 0; i < feature_map.getProteinIdentifications().size(); ++i )
		{
			const ProteinIdentification& current_prot_id = feature_map.getProteinIdentifications()[i];
			os << "\t<IdentificationRun ";
			os << "id=\"PI_" << i << "\" ";
			identifier_id_[current_prot_id.getIdentifier()] = String("PI_") + i;
			os << "date=\"" << current_prot_id.getDateTime().getDate() << "T" << current_prot_id.getDateTime().getTime() << "\" ";
			os << "search_engine=\"" << current_prot_id.getSearchEngine() << "\" ";
			os << "search_engine_version=\"" << current_prot_id.getSearchEngineVersion() << "\">\n";

			//write search parameters
			const ProteinIdentification::SearchParameters& search_param = current_prot_id.getSearchParameters();
			os << "\t\t<SearchParameters "
				 << "db=\"" << search_param.db << "\" "
				 << "db_version=\"" << search_param.db_version << "\" "
				 << "taxonomy=\"" << search_param.taxonomy << "\" ";
			if (search_param.mass_type == ProteinIdentification::MONOISOTOPIC)
			{
				os << "mass_type=\"monoisotopic\" ";
			}
			else if (search_param.mass_type == ProteinIdentification::AVERAGE)
			{
				os << "mass_type=\"average\" ";
			}
			os << "charges=\"" << search_param.charges << "\" ";
			if (search_param.enzyme == ProteinIdentification::TRYPSIN)
			{
				os << "enzyme=\"trypsin\" ";
			}
			if (search_param.enzyme == ProteinIdentification::PEPSIN_A)
			{
				os << "enzyme=\"pepsin_a\" ";
			}
			if (search_param.enzyme == ProteinIdentification::PROTEASE_K)
			{
				os << "enzyme=\"protease_k\" ";
			}
			if (search_param.enzyme == ProteinIdentification::CHYMOTRYPSIN)
			{
				os << "enzyme=\"chymotrypsin\" ";
			}
			else if (search_param.enzyme == ProteinIdentification::NO_ENZYME)
			{
				os << "enzyme=\"no_enzyme\" ";
			}
			else if (search_param.enzyme == ProteinIdentification::UNKNOWN_ENZYME)
			{
				os << "enzyme=\"unknown_enzyme\" ";
			}
			os << "missed_cleavages=\"" << search_param.missed_cleavages << "\" "
				 << "precursor_peak_tolerance=\"" << search_param.precursor_tolerance << "\" "
				 << "peak_mass_tolerance=\"" << search_param.peak_mass_tolerance << "\" "
				 << ">\n";

			//modifications
			for (Size j=0; j!=search_param.fixed_modifications.size(); ++j)
			{
				os << "\t\t\t<FixedModification name=\"" << search_param.fixed_modifications[j] << "\" />\n";
				//Add MetaInfo, when modifications has it (Andreas)
			}
			for (Size j=0; j!=search_param.variable_modifications.size(); ++j)
			{
				os << "\t\t\t<VariableModification name=\"" << search_param.variable_modifications[j] << "\" />\n";
				//Add MetaInfo, when modifications has it (Andreas)
			}

			writeUserParam_("UserParam", os, search_param, 4);

			os << "\t\t</SearchParameters>\n";

			//write protein identifications
			os << "\t\t<ProteinIdentification";
			os << " score_type=\"" << current_prot_id.getScoreType() << "\"";
			os << " higher_score_better=\"" << ( current_prot_id.isHigherScoreBetter() ? "true" : "false" ) << "\"";
			os << " significance_threshold=\"" << current_prot_id.getSignificanceThreshold() << "\">\n";

			// write protein hits
			for (Size j=0; j<current_prot_id.getHits().size(); ++j)
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
		for ( Size i = 0; i < feature_map.getUnassignedPeptideIdentifications().size(); ++i )
		{
			writePeptideIdentification_(filename, os, feature_map.getUnassignedPeptideIdentifications()[i], "UnassignedPeptideIdentification", 1);
		}

		// write features with their corresponding attributes
		os << "\t<featureList count=\"" << feature_map.size() << "\">\n";
		for (Size s=0; s<feature_map.size(); s++)
		{
      writeFeature_(filename, os, feature_map[s], "f_", feature_map[s].getUniqueId(), 0);
      // writeFeature_(filename, os, feature_map[s], "f_", s, 0);
		}

		os << "\t</featureList>\n";
		os << "</featureMap>\n";

		//Clear members
		accession_to_id_.clear();
		identifier_id_.clear();
	}

	PeakFileOptions& FeatureXMLFile::getOptions()
	{
		return options_;
	}

  const PeakFileOptions& FeatureXMLFile::getOptions() const
  {
  	return options_;
  }

	void FeatureXMLFile::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
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
		String parent_tag;
		if (open_tags_.size()!=0) parent_tag = open_tags_.back();
		open_tags_.push_back(tag);

		//for downward compatibility, all tags in the old description must be ignored
		if (in_description_) return;

		if (tag=="description")
		{
			in_description_ = true;
		}
		else if (tag=="feature")
		{
			// create new feature at appropriate level
			updateCurrentFeature_(true);
			current_feature_->setUniqueId(attributeAsString_(attributes,s_id));
		}
		else if (tag=="subordinate")
		{ // this is not safe towards malformed xml!
			++subordinate_feature_level_;
		}
		else if (tag=="featureList")
		{
			if (options_.getMetadataOnly()) throw EndParsingSoftly(__FILE__,__LINE__,__PRETTY_FUNCTION__);
		}
		else if (tag=="quality" || tag=="hposition" || tag=="position")
		{
			dim_ = attributeAsInt_(attributes,s_dim);
		}
    else if (tag=="pt")
		{
      hull_position_[0] = attributeAsDouble_(attributes,"x");
      hull_position_[1] = attributeAsDouble_(attributes,"y");
    }
		else if (tag=="convexhull")
		{
			current_chull_.clear();
		}
		else if (tag=="hullpoint")
		{
			hull_position_ = DPosition<2>::zero();
		}
		else if (tag=="model")
		{
			model_desc_ = new ModelDescription<2>();
			param_.clear();
			model_desc_->setName(attributeAsString_(attributes,s_name));
		}
		else if (tag=="param")
		{
			String name = attributeAsString_(attributes,s_name);
			String value = attributeAsString_(attributes,s_value);
			if (name != "" && value != "") param_.setValue(name, value);
		}
		else if (tag == "userParam")
		{
			if (last_meta_ == 0)
			{
				fatalError(LOAD, String("Unexpected userParam in tag '") + parent_tag + "'" );
			}

			String name = attributeAsString_(attributes,s_name);
			String type = attributeAsString_(attributes,s_type);

			if (type == "int")
			{
				last_meta_->setMetaValue(name, attributeAsInt_(attributes,s_value));
			}
			else if (type == "float")
			{
				last_meta_->setMetaValue(name, attributeAsDouble_(attributes,s_value));
			}
			else if (type == "string")
			{
				last_meta_->setMetaValue(name, (String)attributeAsString_(attributes,s_value));
			}
		        else if ( type == "intList" )
		        {
		          last_meta_->setMetaValue(name, attributeAsIntList_(attributes, "value"));
		        }
		        else if ( type == "floatList" )
		        {
		          last_meta_->setMetaValue(name, attributeAsDoubleList_(attributes, "value"));
		        }
		        else if ( type == "stringList" )
		        {
		          last_meta_->setMetaValue(name, attributeAsStringList_(attributes, "value"));
		        }
			else
			{
				fatalError(LOAD, String("Invalid userParam type '") + type + "'" );
			}
		}
		else if (tag=="featureMap")
		{
			//check file version against schema version
			String file_version="";
			optionalAttributeAsString_(file_version,attributes,s_version);
			if (file_version == "") 
			{
				file_version="1.0"; //default version is 1.0
			}
			if (file_version.toDouble() > version_.toDouble())
			{
				warning(LOAD, String("The XML file (") + file_version +") is newer than the parser (" + version_ + "). This might lead to undefinded program behaviour.");
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
		}
		else if (tag=="dataProcessing")
		{
			DataProcessing tmp;
			tmp.setCompletionTime(asDateTime_(attributeAsString_(attributes, s_completion_time)));
			map_->getDataProcessing().push_back(tmp);
			last_meta_ = &(map_->getDataProcessing().back());
		}
		else if (tag=="software" && parent_tag=="dataProcessing")
		{
			map_->getDataProcessing().back().getSoftware().setName(attributeAsString_(attributes, s_name));
			map_->getDataProcessing().back().getSoftware().setVersion(attributeAsString_(attributes, s_version));
		}
		else if (tag=="processingAction" && parent_tag=="dataProcessing")
		{
			String name = attributeAsString_(attributes, s_name);
			for (Size i=0; i< DataProcessing::SIZE_OF_PROCESSINGACTION; ++i)
			{
				if (name == DataProcessing::NamesOfProcessingAction[i])
				{
					map_->getDataProcessing().back().getProcessingActions().insert((DataProcessing::ProcessingAction)i);
				}
			}
		}
		else if (tag == "IdentificationRun")
		{
			prot_id_.setSearchEngine(attributeAsString_(attributes,"search_engine"));
			prot_id_.setSearchEngineVersion(attributeAsString_(attributes,"search_engine_version"));
			prot_id_.setDateTime(DateTime::fromString(String(attributeAsString_(attributes,"date")).toQString(),"yyyy-MM-ddThh:mm:ss"));
			//set identifier
			String identifier = prot_id_.getSearchEngine() + '_' + attributeAsString_(attributes,"date");
			prot_id_.setIdentifier(identifier);
			id_identifier_[attributeAsString_(attributes,"id")] = identifier;
		}
		else if (tag =="SearchParameters")
		{
			//load parameters
			search_param_.db = attributeAsString_(attributes,"db");
			search_param_.db_version = attributeAsString_(attributes,"db_version");
			optionalAttributeAsString_(search_param_.taxonomy, attributes,"taxonomy");
			search_param_.charges = attributeAsString_(attributes,"charges");
			optionalAttributeAsUInt_(search_param_.missed_cleavages, attributes,"missed_cleavages");
			search_param_.peak_mass_tolerance = attributeAsDouble_(attributes,"peak_mass_tolerance");
			search_param_.precursor_tolerance = attributeAsDouble_(attributes,"precursor_peak_tolerance");
			//mass type
			String mass_type = attributeAsString_(attributes,"mass_type");
			if (mass_type=="monoisotopic")
			{
				search_param_.mass_type = ProteinIdentification::MONOISOTOPIC;
			}
			else if (mass_type=="average")
			{
				search_param_.mass_type = ProteinIdentification::AVERAGE;
			}
			//enzyme
			String enzyme;
			optionalAttributeAsString_(enzyme,attributes,"enzyme");
			if (enzyme == "trypsin")
			{
				search_param_.enzyme = ProteinIdentification::TRYPSIN;
			}
			else if (enzyme == "pepsin_a")
			{
				search_param_.enzyme = ProteinIdentification::PEPSIN_A;
			}
			else if (enzyme == "protease_k")
			{
				search_param_.enzyme = ProteinIdentification::PROTEASE_K;
			}
			else if (enzyme == "chymotrypsin")
			{
				search_param_.enzyme = ProteinIdentification::CHYMOTRYPSIN;
			}
			else if (enzyme == "no_enzyme")
			{
				search_param_.enzyme = ProteinIdentification::NO_ENZYME;
			}
			else if (enzyme == "unknown_enzyme")
			{
				search_param_.enzyme = ProteinIdentification::UNKNOWN_ENZYME;
			}
			last_meta_ = &search_param_;
		}
		else if (tag =="FixedModification")
		{
			search_param_.fixed_modifications.push_back(attributeAsString_(attributes,"name"));
			//change this line as soon as there is a MetaInfoInterface for modifications (Andreas)
			last_meta_ = 0;
		}
		else if (tag =="VariableModification")
		{
			search_param_.variable_modifications.push_back(attributeAsString_(attributes,"name"));
			//change this line as soon as there is a MetaInfoInterface for modifications (Andreas)
			last_meta_ = 0;
		}
		else if ( tag == "ProteinIdentification" )
		{
			prot_id_.setScoreType(attributeAsString_(attributes,"score_type"));

			//optional significance threshold
			DoubleReal tmp=0.0;
			optionalAttributeAsDouble_(tmp,attributes,"significance_threshold");
			if (tmp!=0.0)
			{
				prot_id_.setSignificanceThreshold(tmp);
			}

			//score orientation
			prot_id_.setHigherScoreBetter(asBool_(attributeAsString_(attributes,"higher_score_better")));

			last_meta_ = &prot_id_;
		}
		else if (tag == "ProteinHit")
		{
			prot_hit_ = ProteinHit();
			String accession = attributeAsString_(attributes,"accession");
			prot_hit_.setAccession(accession);
			prot_hit_.setScore(attributeAsDouble_(attributes,"score"));

			//sequence
			String tmp="";
			optionalAttributeAsString_(tmp,attributes,"sequence");
			prot_hit_.setSequence(tmp);

			last_meta_ = &prot_hit_;

			//insert id and accession to map
			proteinid_to_accession_[attributeAsString_(attributes,"id")] = accession;
		}
		else if (tag == "PeptideIdentification" || tag == "UnassignedPeptideIdentification")
		{
			String id = attributeAsString_(attributes,"identification_run_ref");
			if (!id_identifier_.has(id))
			{
				warning(LOAD, String("Peptide identification without ProteinIdentification found (id: '") + id + "')!");
			}
			pep_id_.setIdentifier(id_identifier_[id]);

			pep_id_.setScoreType(attributeAsString_(attributes,"score_type"));

			//optional significance threshold
			DoubleReal tmp=0.0;
			optionalAttributeAsDouble_(tmp,attributes,"significance_threshold");
			if (tmp!=0.0)
			{
				pep_id_.setSignificanceThreshold(tmp);
			}

			//score orientation
			pep_id_.setHigherScoreBetter(asBool_(attributeAsString_(attributes,"higher_score_better")));

			//MZ
			DoubleReal tmp2 = - numeric_limits<DoubleReal>::max();
			optionalAttributeAsDouble_(tmp2, attributes,"MZ");
			if ( tmp2 != - numeric_limits<DoubleReal>::max() )
			{
				pep_id_.setMetaValue("MZ", tmp2);
			}
			//RT
			tmp2 = - numeric_limits<DoubleReal>::max();
			optionalAttributeAsDouble_(tmp2, attributes,"RT");
			if ( tmp2 != - numeric_limits<DoubleReal>::max())
			{
				pep_id_.setMetaValue("RT", tmp2);
			}
			Int tmp3 = - numeric_limits<Int>::max();
			optionalAttributeAsInt_(tmp3, attributes,"spectrum_reference");
			if (tmp3 != - numeric_limits<Int>::max())
			{
				pep_id_.setMetaValue("spectrum_reference", tmp3);
			}

			last_meta_ = &pep_id_;
		}
		else if (tag == "PeptideHit")
		{
			pep_hit_ = PeptideHit();

			pep_hit_.setCharge(attributeAsInt_(attributes,"charge"));
			pep_hit_.setScore(attributeAsDouble_(attributes,"score"));
			pep_hit_.setSequence(attributeAsString_(attributes,"sequence"));

			//aa_before
			String tmp="";
			optionalAttributeAsString_(tmp,attributes,"aa_before");
			if (!tmp.empty())
			{
				pep_hit_.setAABefore(tmp[0]);
			}
			//aa_after
			tmp="";
			optionalAttributeAsString_(tmp,attributes,"aa_after");
			if (!tmp.empty())
			{
				pep_hit_.setAAAfter(tmp[0]);
			}

			//parse optional protein ids to determine accessions
			const XMLCh* refs = attributes.getValue(sm_.convert("protein_refs"));
			if (refs!=0)
			{
				String accession_string = sm_.convert(refs);
				accession_string.trim();
				vector<String> accessions;
				accession_string.split(' ', accessions);
				if (accession_string!="" && accessions.empty())
				{
					accessions.push_back(accession_string);
				}
				for(vector<String>::const_iterator it = accessions.begin(); it!=accessions.end(); ++it)
				{
					Map<String,String>::const_iterator it2 = proteinid_to_accession_.find(*it);
					if (it2!=proteinid_to_accession_.end())
					{
						pep_hit_.addProteinAccession(it2->second);
					}
					else
					{
						fatalError(LOAD, String("Invalid protein reference '") + *it + "'" );
					}
				}
			}
			last_meta_ = &pep_hit_;
		}
	}

	void FeatureXMLFile::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
	{
		String tag = sm_.convert(qname);
		open_tags_.pop_back();

		//for downward compatibility, all tags in the old description must be ignored
		if (tag=="description")
		{
			in_description_ = false;
		}
		if (in_description_) return;

		if (tag=="feature")
		{
			if ((!options_.hasRTRange() || options_.getRTRange().encloses(current_feature_->getRT()))
					&&	(!options_.hasMZRange() || options_.getMZRange().encloses(current_feature_->getMZ()))
					&&	(!options_.hasIntensityRange() || options_.getIntensityRange().encloses(current_feature_->getIntensity())))
			{
			}
			else
			{
				// this feature does not pass the restrictions --> remove it
				if (subordinate_feature_level_==0)
				{
					map_->pop_back();
				}
				else
				{
					Feature* f1(0);
					if (!map_->empty())
					{
						f1 = &(map_->back());
					}
					else
					{
						fatalError(LOAD, "Feature with unexpected location.");
					}

					for (Int level = 1; level<subordinate_feature_level_;++level)
					{
						f1 = &(f1->getSubordinates().back());
					}
					// delete the offending feature
					f1->getSubordinates().pop_back();
				}
			}
			updateCurrentFeature_(false);
		}
		else if (tag=="model")
		{
			model_desc_->setParam(param_);
			current_feature_->setModelDescription(*model_desc_);
			delete model_desc_;
		}
		else if (tag=="hullpoint" || tag=="pt")
		{
			current_chull_.push_back(hull_position_);
		}
		else if (tag=="convexhull")
		{
			ConvexHull2D hull;
			hull.setHullPoints(current_chull_);
			current_feature_->getConvexHulls().push_back(hull);
		}
		else if (tag=="subordinate")
		{
			--subordinate_feature_level_;
			// reset current_feature
			updateCurrentFeature_(false);
		}
		else if (tag == "IdentificationRun")
		{
			map_->getProteinIdentifications().push_back(prot_id_);
			prot_id_ = ProteinIdentification();
			last_meta_  = 0;
		}
		else if (tag == "SearchParameters")
		{
			prot_id_.setSearchParameters(search_param_);
			search_param_ = ProteinIdentification::SearchParameters();
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
			last_meta_  = 0;
		}
		else if (tag == "PeptideHit")
		{
			pep_id_.insertHit(pep_hit_);
			last_meta_ = &pep_id_;
		}
	}

	void FeatureXMLFile::characters(const XMLCh* const chars, const XMLSize_t /*length*/)
	{
		//for downward compatibility, all tags in the old description must be ignored
		if (in_description_) return;

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

	void FeatureXMLFile::writeFeature_(const String& filename, ostream& os, const Feature& feat, const String& identifier_prefix, UInt64 identifier, UInt indentation_level)
	{
		String indent = String(indentation_level,'\t');

		os << indent << "\t\t<feature id=\"" << identifier_prefix << identifier << "\">\n";
		for (Size i=0; i<2;i++)
		{
			os << indent <<	"\t\t\t<position dim=\"" << i << "\">" << precisionWrapper(feat.getPosition()[i]) << "</position>\n";
		}
		os << indent << "\t\t\t<intensity>" << precisionWrapper(feat.getIntensity()) << "</intensity>\n";
		for (Size i = 0; i < 2; i++)
		{
			os << indent << "\t\t\t<quality dim=\"" << i << "\">" << precisionWrapper(feat.getQuality(i)) << "</quality>\n";
		}
		os << indent << "\t\t\t<overallquality>" << precisionWrapper(feat.getOverallQuality()) << "</overallquality>\n";
		os << indent << "\t\t\t<charge>" << feat.getCharge() << "</charge>\n";

		// write model description
		ModelDescription<2> desc = feat.getModelDescription();
		if (!desc.getName().empty() || !desc.getParam().empty())
		{
			os << indent << "\t\t\t<model name=\"" << desc.getName() << "\">\n";
			Param modelp = desc.getParam();
			Param::ParamIterator piter = modelp.begin();
			while (piter != modelp.end())
			{
				os << indent << "\t\t\t\t<param name=\"" << piter.getName() << "\" value=\"" << piter->value << "\"/>\n";
				piter++;
			}
			os << indent << "\t\t\t</model>\n";
		}

		// write convex hull
		vector<ConvexHull2D> hulls = feat.getConvexHulls();

		Size hulls_count = hulls.size();

		for (Size i = 0;i < hulls_count; i++)
		{
			os << indent << "\t\t\t<convexhull nr=\"" << i << "\">\n";

			ConvexHull2D current_hull = hulls[i];
      current_hull.compress();
			Size hull_size	= current_hull.getHullPoints().size();

			for (Size j=0;j<hull_size;j++)
			{
				DPosition<2> pos = current_hull.getHullPoints()[j];
			/*Size pos_size = pos.size();
				os << indent << "\t\t\t\t<hullpoint>\n";
        for (Size k=0; k<pos_size; k++)
				{
					os << indent << "\t\t\t\t\t<hposition dim=\"" << k << "\">" << precisionWrapper(pos[k]) << "</hposition>\n";
				}
				os << indent << "\t\t\t\t</hullpoint>\n";*/
        os << indent << "\t\t\t\t<pt x=\"" << precisionWrapper(pos[0]) << "\" y=\""<< precisionWrapper(pos[1]) << "\" />\n";
			}

			os << indent << "\t\t\t</convexhull>\n";
		}

		if (!feat.getSubordinates().empty())
		{
			os << indent << "\t\t\t<subordinate>\n";
			for (size_t i=0;i<feat.getSubordinates().size();++i)
			{
			  // These subordinate identifiers are a bit long, but who cares about subordinates anyway?  :-P
			  // This way the parent stands out clearly.  However,
			  // note that only the portion after the last '_' is parsed when this is read back.
			  writeFeature_(filename, os, feat.getSubordinates()[i], identifier_prefix+identifier+"_", feat.getSubordinates()[i].getUniqueId(), indentation_level+2);
			}
			os << indent << "\t\t\t</subordinate>\n";
		}

		// write PeptideIdentification
		for ( Size i = 0; i < feat.getPeptideIdentifications().size(); ++i )
		{
			writePeptideIdentification_(filename, os, feat.getPeptideIdentifications()[i], "PeptideIdentification", 3);
		}

		writeUserParam_("userParam", os, feat, indentation_level + 3);

		os << indent << "\t\t</feature>\n";
	}

	void FeatureXMLFile::writePeptideIdentification_(const String& filename, std::ostream& os, const PeptideIdentification& id, const String& tag_name, UInt indentation_level)
	{
		String indent = String(indentation_level,'\t');

		if (!identifier_id_.has(id.getIdentifier()))
		{
			warning(STORE, String("Omitting peptide identification because of missing ProteinIdentification with identifier '") + id.getIdentifier() + "' while writing '" + filename + "'!");
			return;
		}
		os << indent << "<" << tag_name << " ";
		os << "identification_run_ref=\"" << identifier_id_[id.getIdentifier()] << "\" ";
		os << "score_type=\"" << id.getScoreType() << "\" ";
		os << "higher_score_better=\"" << ( id.isHigherScoreBetter() ? "true" : "false" ) << "\" ";
		os << "significance_threshold=\"" << id.getSignificanceThreshold() << "\" ";
		//mz
		DataValue dv = id.getMetaValue("MZ");
		if (dv!=DataValue::EMPTY)
		{
			os << "MZ=\"" << dv.toString() << "\" ";
		}
		// rt
		dv = id.getMetaValue("RT");
		if (dv!=DataValue::EMPTY)
		{
			os << "RT=\"" << dv.toString() << "\" ";
		}
		// spectrum_reference
		dv = id.getMetaValue("spectrum_reference");
		if (dv!=DataValue::EMPTY)
		{
			os << "spectrum_reference=\"" << dv.toString() << "\" ";
		}
		os << ">\n";

		// write peptide hits
		for (Size j=0; j<id.getHits().size(); ++j)
		{
			os << indent << "\t<PeptideHit";
			os << " score=\"" << id.getHits()[j].getScore() << "\"";
			os << " sequence=\"" << id.getHits()[j].getSequence() << "\"";
			os << " charge=\"" << id.getHits()[j].getCharge() << "\"";
			if (id.getHits()[j].getAABefore()!=' ')
			{
				os << " aa_before=\"" << id.getHits()[j].getAABefore() << "\"";
			}
			if (id.getHits()[j].getAAAfter()!=' ')
			{
				os << " aa_after=\"" << id.getHits()[j].getAAAfter() << "\"";
			}
			if(id.getHits()[j].getProteinAccessions().size()!=0)
			{
				String accs = "";
				for (Size m=0; m<id.getHits()[j].getProteinAccessions().size(); ++m)
				{
					if (m) accs += " ";
					accs += "PH_";
					accs += String(accession_to_id_[id.getIdentifier() + "_" + id.getHits()[j].getProteinAccessions()[m]]);
				}
				os << " protein_refs=\"" << accs << "\"";
			}
			os << ">\n";
			writeUserParam_("userParam", os, id.getHits()[j], indentation_level+2);
			os << indent << "\t</PeptideHit>\n";
		}

		//do not write "RT", "MZ" and "spectrum_reference" as they are written as attributes already
		MetaInfoInterface tmp = id;
		tmp.removeMetaValue("RT");
		tmp.removeMetaValue("MZ");
		tmp.removeMetaValue("spectrum_reference");
		writeUserParam_("userParam", os, tmp, indentation_level+1);
		os << indent << "</" << tag_name << ">\n";
	}

	void FeatureXMLFile::updateCurrentFeature_(bool create)
	{
		if (subordinate_feature_level_==0)
		{
			if (create)
			{
				map_->push_back(Feature());
				current_feature_ = &map_->back();
				last_meta_ =  &map_->back();
			}
			else
			{
				if (map_->empty())
				{
					current_feature_ = 0;
					last_meta_ =  0;
				}
				else
				{
					current_feature_ = &map_->back();
					last_meta_ =  &map_->back();
				}
			}
			return;
		}

		Feature* f1 = 0;
		if (map_->empty())
		{
			// do NOT throw an exception here. this is a valid case!	e.g. the
			// only one feature in a map was discarded during endElement(), thus
			// the map_ is empty() now and we cannot assign a current_feature,
			// because there is none!
			current_feature_ = 0;
			last_meta_ = 0;
			return;
		}
		else
		{
			f1 = &map_->back();
		}

		for (Int level = 1; level<subordinate_feature_level_;++level)
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
			f1->getSubordinates().push_back(Feature());
			current_feature_ = &f1->getSubordinates().back();
			last_meta_ = &f1->getSubordinates().back();
			return;
		}
		else
		{
			if (f1->getSubordinates().empty())
			{
				current_feature_ = 0;
				last_meta_ = 0;
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
