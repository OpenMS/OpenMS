// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/AnalysisXMLHandler.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <xercesc/sax2/Attributes.hpp>

#include <iostream>
#include <map>

using namespace std;
using namespace xercesc;

namespace OpenMS
{
	namespace Internal
	{
  AnalysisXMLHandler::AnalysisXMLHandler(vector<ProteinIdentification>& protein_identifications,
  									 std::vector<IdentificationData>& id_data,
   									 const String& filename ) :
    XMLHandler(filename),
    protein_identifications_(&protein_identifications),
    id_data_(&id_data),
    actual_protein_hit_(),
    actual_protein_hits_(),
    actual_peptide_hit_(),
    actual_peptide_hits_(),
    peptide_identification_index_(0),
    protein_identification_index_(0),
    inside_peptide_(false),
    const_protein_identifications_(),
    const_id_data_(),
		const_predicted_retention_times_(),
		tag_(),        
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    predicted_retention_times_(),
    predicted_sigma_(0),
    const_predicted_sigma_(0),
    date_times_temp_(),
    date_times_counter_(0),
    actual_date_time_("0000-00-00 00:00:00")
  {
  }

  AnalysisXMLHandler::AnalysisXMLHandler(vector<ProteinIdentification>& protein_identifications, 
      									 std::vector<IdentificationData>& id_data,
      									 map<String, double>& predicted_retention_times,
      									 DoubleReal& predicted_sigma,
   									 		 const String& filename) :
    XMLHandler(filename),
    protein_identifications_(&protein_identifications),
    id_data_(&id_data),
    actual_protein_hit_(),
    actual_protein_hits_(),
    actual_peptide_hit_(),
    actual_peptide_hits_(),
    peptide_identification_index_(0),
    protein_identification_index_(0),
    inside_peptide_(false),
    const_protein_identifications_(),
    const_id_data_(),
		const_predicted_retention_times_(),
		tag_(),
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    predicted_retention_times_(&predicted_retention_times),
    predicted_sigma_(&predicted_sigma),
    const_predicted_sigma_(0),
    date_times_temp_(),
    date_times_counter_(0),
    actual_date_time_("0000-00-00 00:00:00")
  {
  }

  AnalysisXMLHandler::AnalysisXMLHandler(const vector<ProteinIdentification>& protein_identifications,
  									 const std::vector<IdentificationData>& id_data,
   									 const String& filename) :
    XMLHandler(filename),
    protein_identifications_(),
    id_data_(0),
    actual_protein_hit_(),
    actual_protein_hits_(),
    actual_peptide_hit_(),
    actual_peptide_hits_(),
    peptide_identification_index_(0),
    protein_identification_index_(0),
    inside_peptide_(false),
    const_protein_identifications_(protein_identifications),
    const_id_data_(id_data),
		const_predicted_retention_times_(),
		tag_(),        
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    predicted_retention_times_(),
    predicted_sigma_(0),
    const_predicted_sigma_(0),
    date_times_temp_(),
    date_times_counter_(0),
    actual_date_time_("0000-00-00 00:00:00")    
  {
  }
      									 	      									 	  
  AnalysisXMLHandler::AnalysisXMLHandler(const vector<ProteinIdentification>& protein_identifications, 
      									 const std::vector<IdentificationData>& id_data,
      									 const map<String, double>& const_predicted_retention_times,
      									 DoubleReal predicted_sigma,
   									 const String& filename) :
    XMLHandler(filename),
    protein_identifications_(0),
    id_data_(0),
    actual_protein_hit_(),
    actual_protein_hits_(),
    actual_peptide_hit_(),
    actual_peptide_hits_(),
    peptide_identification_index_(0),
    protein_identification_index_(0),
    inside_peptide_(false),
    const_protein_identifications_(protein_identifications),
    const_id_data_(id_data),
		const_predicted_retention_times_(const_predicted_retention_times),
		tag_(),
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    predicted_retention_times_(),
    predicted_sigma_(),
    const_predicted_sigma_(predicted_sigma),
    date_times_temp_(),
    date_times_counter_(0),
    actual_date_time_("0000-00-00 00:00:00")
  {
  	
  }
   
  AnalysisXMLHandler::~AnalysisXMLHandler()
  {
    
  }
   
  void AnalysisXMLHandler::writeTo(std::ostream& os)
  {
  	vector<ProteinHit> 										temp_protein_hits;
  	vector<PeptideHit> 										temp_peptide_hits;
  	vector<PeptideHit>* 									referencing_peptide_hits;
  	vector<PeptideHit>* 									non_referencing_peptide_hits;
  	vector<PeptideHit>  									all_peptide_hits;
  	vector<UnsignedInt> 									indices;
  	multimap< String, ProteinHit > 				all_protein_hits;
  	
  	bool 																	protein_hits_present = false;
  	String 																temp_peptide_sequence;
  	map< String , UnsignedInt>						date_times;
  	UnsignedInt 													date_times_counter = 0;
  	String 																date_time_string;
  	DateTime															date_time;
  	String																actual_date_time_;
  	map< String , UnsignedInt>::iterator	date_times_iterator;
		UnsignedInt 													group_count = 0;
		
		date_time.now();
		date_time.get(actual_date_time_);
		date_time.clear();
		
		os << "<?xml version=\"1.0\" encoding=\"ISO-8859-1\" standalone=\"yes\"?>\n"
		 //<< "<?xml-stylesheet href=\"identifications.xsl\" type=\"text/xsl\"?>\n"
		 << "<mzIdent>\n"
		 << "\t<description>\n"
		 << "\t\t<admin>\n"
		 << "\t\t\t<sampleName>";
		if (true)
		{
			os << " unknown ";
		}
		else
		{
			os << "";			
		}
		os << "</sampleName>\n"
		 << "\t\t</admin>\n"
		 << "\t\t<settings>\n"
		 << "\t\t\t<userParam name=\"number_of_identifications\" value=\""
  	 << const_id_data_.size() << "\" />\n"
		 << "\t\t\t<userParam name=\"number_of_protein_identifications\" value=\""
  	 << const_protein_identifications_.size() << "\" />\n"
		 << "\t\t\t<userParam name=\"retention_time_type\" value=\""
  	 << "float\" />\n"  	 
		 << "\t\t\t<userParam name=\"mz_type\" value=\""
  	 << "float\" />\n"  	 
		 << "\t\t\t<userParam name=\"retention_time_unit\" value=\""
  	 << "seconds\" />\n"  	 
		 << "\t\t\t<userParam name=\"mz_unit\" value=\""
  	 << "Thomson\" />\n";
  	if (const_predicted_sigma_ != 0)
  	{	 
			os << "\t\t\t<userParam name=\"predicted_sigma\" value=\""
  	 		<< const_predicted_sigma_ << "\" />\n";  	 
  	}
  	for(vector<IdentificationData>::const_iterator it = const_id_data_.begin();
  		  it != const_id_data_.end();
  		  it++)
  	{
  		date_time = it->id.getDateTime();
  		date_time.get(date_time_string);
  		if (date_time_string == "0000-00-00 00:00:00")
  		{
  			date_time_string = actual_date_time_;
  		}
  		if (date_times.find(date_time_string) == date_times.end())
  		{
  			date_times.insert(make_pair(date_time_string, date_times_counter));
  			os << "\t\t\t<userParam name=\"date_group_" << date_times_counter << "\" value=\""
  				 << date_time_string << "\"/>\n";
  			date_times_counter++;
  		}
  	}
  	for(vector<ProteinIdentification>::const_iterator it = const_protein_identifications_.begin();
  		  it != const_protein_identifications_.end();
  		  it++)
  	{
  		date_time = it->getDateTime();
  		date_time.get(date_time_string);
  		if (date_time_string == "0000-00-00 00:00:00")
  		{
  			date_time_string = actual_date_time_;
  		}
  		if (date_times.find(date_time_string) == date_times.end())
  		{
  			date_times.insert(make_pair(date_time_string, date_times_counter));
  			os << "\t\t\t<userParam name=\"date_group_" << date_times_counter << "\" value=\""
  				 << date_time_string << "\"/>\n";
  			date_times_counter++;
  		}
  	}
		os << "\t\t</settings>\n"				
		 << "\t</description>\n"
		 << "\t<ident>\n";

		if (const_protein_identifications_.size() > 0)
		{
			group_count = 0;
			for(UnsignedInt j = 0; j < const_protein_identifications_.size(); j++)
			{
				date_time = const_protein_identifications_[j].getDateTime();
				date_time.get(date_time_string);
	  		if (date_time_string == "0000-00-00 00:00:00")
	  		{
	  			date_time_string = actual_date_time_;
	  		}
									
				os << "\t\t<proteinGroup count=\"" << group_count << "\">\n";
				temp_protein_hits = const_protein_identifications_[j].getProteinHits();
				for(vector<ProteinHit>::const_iterator protein_hits_it = temp_protein_hits.begin();
						protein_hits_it != temp_protein_hits.end();
						protein_hits_it++)
				{
					os << "\t\t\t<protein>\n"
					<< "\t\t\t\t<dbName>";
					if (protein_hits_it->getAccessionType() != "")
					{
						os  << protein_hits_it->getAccessionType();
					}
					os << "</dbName>\n"				 
						<< "\t\t\t\t<dbID>";
					if (protein_hits_it->getAccession() != "")
					{
						os << protein_hits_it->getAccession();
					}
					os << "</dbID>\n";
					for(UnsignedInt i = 0; i < const_id_data_.size(); i++)
					{
						referencing_peptide_hits = const_id_data_[i].id.getReferencingHits(date_time_string, 
																																	 								  protein_hits_it->getAccession());
						for(vector<PeptideHit>::iterator peptide_hits_it = referencing_peptide_hits->begin();
								peptide_hits_it != referencing_peptide_hits->end();
								peptide_hits_it++)
						{							
							writePeptideHit(os, 
  														String("\t\t\t\t"),
  														*peptide_hits_it,
  														const_id_data_[i].id.getPeptideSignificanceThreshold(),
  														i,
  														peptide_hits_it->getCharge(), 
  														const_id_data_[i].rt,
  														const_id_data_[i].mz,
  														const_id_data_[i].id.getDateTime(),
															date_times);

						} // peptide_hits_it
						delete referencing_peptide_hits;
					} // identifications
					os	<< "\t\t\t\t<userParam name=\"protein_identification_index\" value=\"" 
					   << j << "\"/>\n"
					   << "\t\t\t\t<userParam name=\"date_group_index\" value=\"" 
					   << getDateGroupIndex(date_time, date_times) << "\"/>\n"
					   << "\t\t\t\t<userParam name=\"score\" value=\"" 
					   << protein_hits_it->getScore() << "\"/>\n"
					   << "\t\t\t\t<userParam name=\"score_type\" value=\"" 
					   << protein_hits_it->getScoreType() << "\"/>\n"
					    << "\t\t\t</protein>\n";
				} // protein hits
				os << "\t\t\t<userParam name=\"proteins\" value=\"all\"/>\n"
					 << "\t\t</proteinGroup>\n";
			} // protein identifications
		}

		for(UnsignedInt i = 0; i < const_id_data_.size(); i++)
		{
			if (const_id_data_[i].id.getProteinHits().size() > 0)
			{
				protein_hits_present = true;
			}
		}

		if (protein_hits_present)
		{
			os << "\t\t<proteinGroup count=\"" << group_count << "\">\n";
			for(UnsignedInt j = 0; j < const_id_data_.size(); j++)
			{
				date_time = const_id_data_[j].id.getDateTime();
				date_time.get(date_time_string);
	  		if (date_time_string == "0000-00-00 00:00:00")
	  		{
	  			date_time_string = actual_date_time_;
	  		}				
					
				temp_protein_hits = const_id_data_[j].id.getProteinHits();
				for(vector<ProteinHit>::const_iterator protein_hits_it = temp_protein_hits.begin();
						protein_hits_it != temp_protein_hits.end();
						protein_hits_it++)
				{
					os << "\t\t\t<protein>\n"
					<< "\t\t\t\t<dbName>";
					if (protein_hits_it->getAccessionType() != "")
					{
						os  << protein_hits_it->getAccessionType();
					}
					os << "</dbName>\n"				 
						<< "\t\t\t\t<dbID>";
					if (protein_hits_it->getAccession() != "")
					{
						os << protein_hits_it->getAccession();
					}
					os << "</dbID>\n";
					for(UnsignedInt i = 0; i < const_id_data_.size(); i++)
					{
						referencing_peptide_hits = const_id_data_[i].id.getReferencingHits(date_time_string, 
																																						 protein_hits_it->getAccession());
						for(vector<PeptideHit>::iterator peptide_hits_it = referencing_peptide_hits->begin();
								peptide_hits_it != referencing_peptide_hits->end();
								peptide_hits_it++)
						{							
							writePeptideHit(os, 
  														String("\t\t\t\t"),
  														*peptide_hits_it,
  														const_id_data_[i].id.getPeptideSignificanceThreshold(),
  														i,
  														peptide_hits_it->getCharge(),
  														const_id_data_[i].rt,
  														const_id_data_[i].mz,
  														const_id_data_[i].id.getDateTime(),
															date_times);

						} // peptide_hits_it
						delete referencing_peptide_hits;
					} // identifications
					os  << "\t\t\t\t<userParam name=\"identification_index\" value=\""
							<< j << "\"/>\n"
					   << "\t\t\t\t<userParam name=\"date_group_index\" value=\"" 
					   << getDateGroupIndex(date_time, date_times) << "\"/>\n"
					   << "\t\t\t\t<userParam name=\"score\" value=\"" 
					   << protein_hits_it->getScore() << "\"/>\n"
					   << "\t\t\t\t<userParam name=\"score_type\" value=\"" 
					   << protein_hits_it->getScoreType() << "\"/>\n"
					   << "\t\t\t</protein>\n";
				} // protein hits
			} // identifications
			os << "\t\t\t<userParam name=\"proteins\" value=\"all\"/>\n"
				 << "\t\t</proteinGroup>\n";
		}

		// All protein hits with the corresponding date and time are stored.
		// They are used afterwards to determine the peptide hits that do
		// not reference any protein hits.
		for(UnsignedInt i = 0; i < const_id_data_.size(); i++)
		{
			const_id_data_[i].id.getDateTime().get(date_time_string);
  		if (date_time_string == "0000-00-00 00:00:00")
  		{
  			date_time_string = actual_date_time_;
  		}
			
			for(vector<ProteinHit>::const_iterator it = const_id_data_[i].id.getProteinHits().begin();
					it != const_id_data_[i].id.getProteinHits().end();
					it++)
			{
				all_protein_hits.insert(make_pair(date_time_string, *it));				
			}
		}
		for(UnsignedInt i = 0; i < const_protein_identifications_.size(); i++)
		{
			const_protein_identifications_[i].getDateTime().get(date_time_string);
  		if (date_time_string == "0000-00-00 00:00:00")
  		{
  			date_time_string = actual_date_time_;
  		}
			
			for(vector<ProteinHit>::const_iterator it = const_protein_identifications_[i].getProteinHits().begin();
					it != const_protein_identifications_[i].getProteinHits().end();
					it++)
			{
				all_protein_hits.insert(make_pair(date_time_string, *it));				
			}
		}
		for(UnsignedInt i = 0; i < const_id_data_.size(); i++)
		{
			non_referencing_peptide_hits = 
				const_id_data_[i].id.getNonReferencingHits(all_protein_hits);
				 								  														  
																  														 
			for(vector<PeptideHit>::const_iterator peptide_hits_it = non_referencing_peptide_hits->begin();
					peptide_hits_it != 	non_referencing_peptide_hits->end();
					peptide_hits_it++)
			{	
				writePeptideHit(os, 
												String("\t\t"),
												*peptide_hits_it,
												const_id_data_[i].id.getPeptideSignificanceThreshold(),
												i,
												peptide_hits_it->getCharge(), 
												const_id_data_[i].rt,
												const_id_data_[i].mz,
												const_id_data_[i].id.getDateTime(),
												date_times);
			}
			delete non_referencing_peptide_hits;
		}

		os << "\t</ident>\n"
		<< "</mzIdent>";				 
		 
  }

  void AnalysisXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		tag_ = String(xercesc::XMLString::transcode(qname));
		
		String attribute_value;

		if (tag_ == "userParam")
		{
			attribute_value = String(XMLString::transcode(attributes.getValue(0u))).trim();
			
			if (attribute_value == "peptide_significance_threshold")
			{
				(*id_data_)[peptide_identification_index_].id.setPeptideSignificanceThreshold(
					((String) XMLString::transcode(attributes.getValue(1u))).toFloat());
			}
			else if (attribute_value == "precursor_charge")
			{
				if (inside_peptide_)
				{
					actual_peptide_hit_.setCharge(((String) XMLString::transcode(attributes.getValue(1u))).toInt());
				}
			}
			else if (attribute_value == "precursor_retention_time")
			{
				(*id_data_)[peptide_identification_index_].rt = 
					((String) XMLString::transcode(attributes.getValue(1u))).toFloat();
			}
			else if (attribute_value == "predicted_sigma")
			{
				if (predicted_sigma_ != 0)
				{
					*predicted_sigma_ = ((String) XMLString::transcode(attributes.getValue(1u))).toFloat();
				}
			}
			else if (attribute_value == "predicted_precursor_retention_time")
			{
				if (predicted_sigma_ != 0)
				{
					predicted_retention_times_->insert(make_pair(actual_peptide_hit_.getSequence(), 
						((String) XMLString::transcode(attributes.getValue(1u))).toFloat()));
				}
			}
			else if (attribute_value == "number_of_identifications" || attribute_value == "number_of_db_searches")
			{
				id_data_->resize(((String) XMLString::transcode(attributes.getValue(1u))).toInt());
			}
			else if (attribute_value == "number_of_protein_identifications")
			{
				for(int i = 0; i < ((String) XMLString::transcode(attributes.getValue(1u))).toInt(); i++)
				{
					ProteinIdentification temp_protein_identification;
					protein_identifications_->push_back(temp_protein_identification);
				}
			}
			else if (attribute_value == "precursor_mz")
			{
				(*id_data_)[peptide_identification_index_].mz = ((String) XMLString::transcode(attributes.getValue(1u))).toFloat();
			}
			else if (attribute_value == "score")
			{
				if (inside_peptide_)
				{
					actual_peptide_hit_.setScore(((String) XMLString::transcode(attributes.getValue(1u))).toFloat());
				}
				else
				{
					actual_protein_hit_.setScore(((String) XMLString::transcode(attributes.getValue(1u))).toFloat());
				}
			}
			else if (attribute_value == "score_type")
			{
				if (inside_peptide_)
				{
					actual_peptide_hit_.setScoreType(((String) XMLString::transcode(attributes.getValue(1u))));										
				}
				else
				{
					actual_protein_hit_.setScoreType(((String) XMLString::transcode(attributes.getValue(1u))));					
				}
			}
			else if (attribute_value == "identification_index")
			{
				if (inside_peptide_)
				{
					peptide_identification_index_ = ((String) XMLString::transcode(attributes.getValue(1u))).toInt();
					if (peptide_identification_index_ >= id_data_->size())
					{ 
						throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,String(peptide_identification_index_),"Peptide identification_index larger then the overall number of peptides");
					}
				}
				else
				{
					protein_identification_index_ = ((String) XMLString::transcode(attributes.getValue(1u))).toInt();
					if (protein_identification_index_ >= protein_identifications_->size())
					{ 
						throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,String(protein_identification_index_),"Protein identification_index larger then the overall number of proteins");
					}
				}
			}
			else if (attribute_value == "protein_identification_index")
			{
				protein_identification_index_ = ((String) XMLString::transcode(attributes.getValue(1u))).toInt();
				inside_global_protein_ = true;
			}
			else if (attribute_value.hasPrefix("date_group_") 
							 && attribute_value != "date_group_index")
			{
				date_times_temp_.push_back((String) XMLString::transcode(attributes.getValue(1u)));
				date_times_counter_++;
			}
			else if (attribute_value == "date_group_index")
			{		
				UnsignedInt index = ((String) XMLString::transcode(attributes.getValue(1u))).toInt();
				
				if (index>=date_times_temp_.size())
				{
					throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,String(index),"Undefined date_group index");
				}

				String date_time_string = date_times_temp_.at(index);
					if (inside_peptide_)
					{
						(*id_data_)[peptide_identification_index_].id.getDateTime().set(date_time_string);
					}
					else
					{
						if (inside_global_protein_)
						{
							(*protein_identifications_)[protein_identification_index_].getDateTime().set(date_time_string);
						}
						else
						{
							(*id_data_)[protein_identification_index_].id.getDateTime().set(date_time_string);
						}
					}
			}						
		}
		else if (tag_ == "peptide")
		{
			inside_peptide_ = true;			
		}
		else if (tag_ == "protein")
		{
                        inside_global_protein_ = true;
			inside_protein_ = true;
			inside_peptide_ = false;			
		}
	}

  void AnalysisXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
 	{
 		tag_ = String(xercesc::XMLString::transcode(qname)).trim();
		
 		if (tag_ == "protein")
 		{	
			if (inside_global_protein_)
			{
	 			(*protein_identifications_)[protein_identification_index_].insertProteinHit(actual_protein_hit_);
			}
			else
			{				
	 			(*id_data_)[protein_identification_index_].id.insertProteinHit(actual_protein_hit_);
			}
	 			actual_protein_hit_.clear();
	 			actual_peptide_indices_.clear();
				inside_protein_ = false;
				inside_global_protein_ = false;
 		}
 		else if (tag_ == "peptide")
 		{
			bool													already_stored = false;
			vector<PeptideHit>::iterator  it;
 			
			vector<PeptideHit>& temp_peptide_hits = 
				(*id_data_)[peptide_identification_index_].id.getPeptideHits();
				
			it = temp_peptide_hits.begin();
			while(it != temp_peptide_hits.end() && !already_stored)
			{
				if (it->getSequence() == actual_peptide_hit_.getSequence())
				{
					already_stored = true;
					if (inside_protein_)
					{
		 				String 							date_time;
		 			
		 				(*id_data_)[peptide_identification_index_].id.getDateTime().get(date_time);
						it->addProteinIndex(make_pair(date_time, actual_protein_hit_.getAccession()));						
					}
				}
				it++;
			}
			if (!already_stored)
			{
				if (inside_protein_)
				{
	 				String 							date_time;
	 			
	 				(*id_data_)[peptide_identification_index_].id.getDateTime().get(date_time);
					actual_peptide_hit_.addProteinIndex(make_pair(date_time, actual_protein_hit_.getAccession()));
				}
	 			(*id_data_)[peptide_identification_index_].id.insertPeptideHit(actual_peptide_hit_); 			
			}
 			actual_peptide_hit_.clear();
 			inside_peptide_ = false;
 		}

		tag_ = "";
 	} 

  void AnalysisXMLHandler::characters(const XMLCh* const chars, const unsigned int /*length*/)
  {
		if (tag_ == "dbName")
		{
			actual_protein_hit_.setAccessionType(((String) xercesc::XMLString::transcode(chars)).trim());
			tag_ = "";
		}
		else if (tag_ == "dbID")
		{
			actual_protein_hit_.setAccession(((String) xercesc::XMLString::transcode(chars)).trim());			
			tag_ = "";
		}
		else if (tag_ == "sequence")
		{
			actual_protein_hit_.setSequence(((String) xercesc::XMLString::transcode(chars)).trim());			
			tag_ = "";
		}
		else if (tag_ == "seq")
		{
			actual_peptide_hit_.setSequence(((String) xercesc::XMLString::transcode(chars)).trim());			
			tag_ = "";
		}
  }
  
  UnsignedInt AnalysisXMLHandler::getDateGroupIndex(DateTime 											date_time,
			  																						map< String , UnsignedInt> 		date_times)
	{
		String 																date_time_string 					= "";
		UnsignedInt 													date_group_index					= 0;
		map< String , UnsignedInt>::iterator 	date_times_iterator;
		
		// determining the date group
		date_time.get(date_time_string);
 		if (date_time_string == "0000-00-00 00:00:00")
 		{
 			date_time_string = actual_date_time_;
 		}
		
		date_times_iterator = date_times.find(date_time_string);
		if (date_times_iterator != date_times.end())
		{
			date_group_index = date_times_iterator->second;
		}
		else
		{
			date_group_index = 0;
		}
  	
		return date_group_index;		
	}
  
  void AnalysisXMLHandler::writePeptideHit(ostream& os, 
  																				String shift,
			  																	PeptideHit hit,
			  																	Real significance_threshold,
			  																	UnsignedInt identification_index,
			  																	SignedInt charge, 
			  																	Real precursor_retention_time,
			  																	Real precursor_mz,
			  																	DateTime date_time,
			  																	map< String , UnsignedInt> date_times)
  {
  	String 																temp_peptide_sequence 		= "";
		Real 																	predicted_retention_time 	= -1;
  	map<String, double>::const_iterator 	predictions_iterator;
		UnsignedInt 													date_group_index					= 0;
		
		// determining the predicted retention time (default: -1)
		if (const_predicted_retention_times_.size() > 0)
		{
			predictions_iterator = 
				const_predicted_retention_times_.find(hit.getSequence());
			if (predictions_iterator != const_predicted_retention_times_.end())
			{
				predicted_retention_time = predictions_iterator->second;
			}
			else
			{
				predicted_retention_time = -1;
			}
		}
		
		date_group_index = getDateGroupIndex(date_time, date_times);
		
		os <<  shift << "<peptide>\n"
		<< shift << "\t<seq>";
		temp_peptide_sequence = hit.getSequence();
		if (temp_peptide_sequence != "")
		{
			os << hit.getSequence();
		}
		os << "</seq>\n"
		<< shift << "\t<userParam name=\"identification_index\" value=\""
		<< identification_index << "\" />\n"
		<<  shift << "\t<userParam name=\"peptide_significance_threshold\" value=\"" 
 		<< significance_threshold << "\"/>\n"
		<<  shift << "\t<userParam name=\"precursor_charge\" value=\"" 
		<< charge << "\" />\n"
		<< shift << "\t<userParam name=\"precursor_retention_time\" value=\""
		<< precursor_retention_time << "\" />\n";						
		if (predicted_retention_time >= 0)
		{
			os << shift << "\t<userParam name=\"predicted_precursor_retention_time\" "
				 << "value=\"" 
				 << predicted_retention_time
				 << "\" />\n";
		}
		os 	<< shift << "\t<userParam name=\"precursor_mz\" value=\""
				<< precursor_mz << "\" />\n"
				<< shift << "\t<userParam name=\"score\" value=\""
				<< hit.getScore() << "\" />\n"
				<< shift << "\t<userParam name=\"score_type\" value=\""
				<< hit.getScoreType() << "\"/>\n"
				<< shift << "\t<userParam name=\"date_group_index\" value=\""
				<< date_group_index << "\"/>\n"
				<<  shift << "</peptide>\n";
	}
  

	} // namespace Internal
} // namespace OpenMS
