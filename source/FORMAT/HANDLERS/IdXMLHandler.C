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

#include <OpenMS/FORMAT/HANDLERS/IdXMLHandler.h>
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
  IdXMLHandler::IdXMLHandler(vector<Identification>& protein_identifications,
  									 std::vector<PeptideIdentification>& id_data,
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
		tag_(),        
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    date_times_temp_(),
    date_times_counter_(0),
    actual_date_time_("0000-00-00 00:00:00"),
    actual_search_parameters_()
  {
  }

  IdXMLHandler::IdXMLHandler(const vector<Identification>& protein_identifications,
  									 const std::vector<PeptideIdentification>& id_data,
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
		tag_(),        
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    date_times_temp_(),
    date_times_counter_(0),
    actual_date_time_("0000-00-00 00:00:00"),
    actual_search_parameters_()    
  {
  }
   
  IdXMLHandler::~IdXMLHandler()
  {
    
  }
   
  void IdXMLHandler::writeTo(std::ostream& os)
  {
  	vector<ProteinHit> 										temp_protein_hits;
  	vector<PeptideHit> 										temp_peptide_hits;
  	vector<PeptideHit>  									all_peptide_hits;
  	vector<UInt> 													indices;
  	vector< String >							 				all_protein_hits;
  	
  	String 																temp_peptide_sequence;
  	map< String , UInt>										date_times;
  	String 																date_time_string;
  	DateTime															date_time;
  	String																actual_date_time_;
  	map< String , UInt>::iterator					date_times_iterator;
		UInt 																	group_count = 0;
		Identification::SearchParameters 			search_parameters;
		
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
  	 << "Thomson\" />\n"
  	 << "\t\t</settings>\n"				
		 << "\t</description>\n"
		 << "\t<ident>\n";

		if (const_protein_identifications_.size() > 0)
		{
			group_count = 0;
			for(UInt j = 0; j < const_protein_identifications_.size(); j++)
			{
				date_time = const_protein_identifications_[j].getDateTime();
	  		if (!date_time.isValid())
	  		{
	  			date_time_string = actual_date_time_;
	  		}
	  		else
	  		{
					date_time.get(date_time_string);
	  		}
									
				os << "\t\t<proteinGroup count=\"" << j << "\">\n";
				temp_protein_hits = const_protein_identifications_[j].getHits();
				String temp_id = const_protein_identifications_[j].getIdentifier();
				vector<PeptideHit> referencing_peptide_hits;
				DoubleReal predicted_rt;
				DoubleReal predicted_rt_p_value;

				for(vector<ProteinHit>::const_iterator protein_hits_it = temp_protein_hits.begin();
						protein_hits_it != temp_protein_hits.end();
						protein_hits_it++)
				{
					os << "\t\t\t<protein>\n"
						<< "\t\t\t\t<dbID>";
					if (protein_hits_it->getAccession() != "")
					{
						os << protein_hits_it->getAccession();
					}
					os << "</dbID>\n";

					for(UInt i = 0; i < const_id_data_.size(); i++)
					{
						if (temp_id == const_id_data_[i].getIdentifier())
						{
							const_id_data_[i].getReferencingHits(protein_hits_it->getAccession(),
																									 referencing_peptide_hits);
							for(vector<PeptideHit>::iterator peptide_hits_it = referencing_peptide_hits.begin();
									peptide_hits_it != referencing_peptide_hits.end();
									peptide_hits_it++)
							{
								if (const_id_data_[i].metaValueExists("predicted_RT"))
								{
									predicted_rt = DoubleReal(const_id_data_[i].getMetaValue("predicted_RT"));
								}
								else
								{
									predicted_rt = -1;
								}							
								if (const_id_data_[i].metaValueExists("predicted_RT_p_value"))
								{
									predicted_rt_p_value = DoubleReal(const_id_data_[i].getMetaValue("predicted_RT_p_value"));
								}
								else
								{
									predicted_rt_p_value = -1;
								}							
								writePeptideHit(os, 
	  														String("\t\t\t\t"),
	  														*peptide_hits_it,
	  														const_id_data_[i].getSignificanceThreshold(),
	  														i,
	  														const_id_data_[i].getMetaValue("RT"),
	  														const_id_data_[i].getMetaValue("MZ"),
																const_id_data_[i].getIdentifier(),
																const_id_data_[i].getScoreType(),
																const_id_data_[i].isHigherScoreBetter(),
	  														predicted_rt,
	  														predicted_rt_p_value);
	
							} // peptide_hits_it
							referencing_peptide_hits.clear();
						}
					} // identifications
					os	<< "\t\t\t\t<userParam name=\"protein_identification_index\" value=\"" 
					   << j << "\"/>\n"
					   << "\t\t\t\t<userParam name=\"score\" value=\"" 
					   << protein_hits_it->getScore() << "\"/>\n"
					    << "\t\t\t</protein>\n";
				} // protein hits

				search_parameters = const_protein_identifications_[j].getSearchParameters();				
				os << "\t\t\t<userParam name=\"db\" value=\"" << search_parameters.db << "\"/>\n"
					 << "\t\t\t<userParam name=\"db_version\" value=\"" << search_parameters.db_version << "\"/>\n"
					 << "\t\t\t<userParam name=\"taxonomy\" value=\"" << search_parameters.taxonomy << "\"/>\n"
					 << "\t\t\t<userParam name=\"charges\" value=\"" << search_parameters.charges << "\"/>\n"
					 << "\t\t\t<userParam name=\"mass_type\" value=\"" << search_parameters.mass_type << "\"/>\n";
				for(vector<String>::iterator param_it = search_parameters.fixed_modifications.begin();
						param_it != search_parameters.fixed_modifications.end();
						++param_it)
				{
					 os << "\t\t\t<userParam name=\"fixed_modifications\" value=\"" << *param_it << "\"/>\n";
				}
				for(vector<String>::iterator param_it = search_parameters.variable_modifications.begin();
						param_it != search_parameters.variable_modifications.end();
						++param_it)
				{				
					 os << "\t\t\t<userParam name=\"variable_modifications\" value=\"" << *param_it << "\"/>\n";
				}
				const_protein_identifications_[j].getDateTime().get(date_time_string);
  			if (date_time_string == "0000-00-00 00:00:00")
  			{
  				date_time_string = actual_date_time_;
  			}
  			
  			if (const_protein_identifications_[j].getIdentifier() == "")
  			{
//  				throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"No identifier set for identification ",String(j));
  			}  			

				os << "\t\t\t<userParam name=\"enzyme\" value=\"" << search_parameters.enzyme << "\"/>\n"
					 << "\t\t\t<userParam name=\"missed_cleavages\" value=\"" << search_parameters.missed_cleavages << "\"/>\n"
					 << "\t\t\t<userParam name=\"peak_mass_tolerance\" value=\"" << search_parameters.peak_mass_tolerance << "\"/>\n"
					 << "\t\t\t<userParam name=\"precursor_tolerance\" value=\"" << search_parameters.precursor_tolerance << "\"/>\n"
					 << "\t\t\t<userParam name=\"date\" value=\"" << date_time_string << "\"/>\n"
					 << "\t\t\t<userParam name=\"search_engine\" value=\"" << const_protein_identifications_[j].getSearchEngine() << "\"/>\n"
					 << "\t\t\t<userParam name=\"search_engine_version\" value=\"" << const_protein_identifications_[j].getSearchEngineVersion() << "\"/>\n"
					 << "\t\t\t<userParam name=\"identifier\" value=\"" << const_protein_identifications_[j].getIdentifier() << "\"/>\n"
					 << "\t\t\t<userParam name=\"score_type\" value=\"" << const_protein_identifications_[j].getScoreType() << "\"/>\n"
					 << "\t\t\t<userParam name=\"higher_score_better\" value=\"" << const_protein_identifications_[j].isHigherScoreBetter() << "\"/>\n"
					 <<	"\t\t\t<userParam name=\"proteins\" value=\"all\"/>\n"
					 << "\t\t</proteinGroup>\n";
			} // protein identifications
		}

		// All protein hits with the corresponding date and time are stored.
		// They are used afterwards to determine the peptide hits that do
		// not reference any protein hits.
		for(UInt i = 0; i < const_protein_identifications_.size(); i++)
		{
			const_protein_identifications_[i].getDateTime().get(date_time_string);
  		if (date_time_string == "0000-00-00 00:00:00")
  		{
  			date_time_string = actual_date_time_;
  		}
			
			for(vector<ProteinHit>::const_iterator it = const_protein_identifications_[i].getHits().begin();
					it != const_protein_identifications_[i].getHits().end();
					it++)
			{
				all_protein_hits.push_back(it->getAccession());				
			}
		}

		vector<PeptideHit> non_referencing_peptide_hits;
		for(UInt i = 0; i < const_id_data_.size(); i++)
		{			
			const_id_data_[i].getNonReferencingHits(all_protein_hits, non_referencing_peptide_hits);				 								  														  
																  														 
			DoubleReal predicted_rt;
			DoubleReal predicted_rt_p_value;
			for(vector<PeptideHit>::const_iterator peptide_hits_it = non_referencing_peptide_hits.begin();
					peptide_hits_it != 	non_referencing_peptide_hits.end();
					peptide_hits_it++)
			{	
				if (const_id_data_[i].metaValueExists("predicted_RT"))
				{
					predicted_rt = DoubleReal(const_id_data_[i].getMetaValue("predicted_RT"));
				}
				else
				{
					predicted_rt = -1;
				}							
				if (const_id_data_[i].metaValueExists("predicted_RT_p_value"))
				{
					predicted_rt_p_value = DoubleReal(const_id_data_[i].getMetaValue("predicted_RT_p_value"));
				}
				else
				{
					predicted_rt_p_value = -1;
				}							

				writePeptideHit(os, 
	  										String("\t\t"),
												*peptide_hits_it,
												const_id_data_[i].getSignificanceThreshold(),
												i,
												const_id_data_[i].getMetaValue("RT"),
												const_id_data_[i].getMetaValue("MZ"),
												const_id_data_[i].getIdentifier(),
												const_id_data_[i].getScoreType(),
												const_id_data_[i].isHigherScoreBetter(),
												predicted_rt,
												predicted_rt_p_value);
			}
		}

		os << "\t</ident>\n"
		<< "</mzIdent>";				 
		 
  }

  void IdXMLHandler::startElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname, const xercesc::Attributes& attributes)
	{
		tag_ = String(xercesc::XMLString::transcode(qname));
		
		String attribute_value;

		if (tag_ == "userParam")
		{
			attribute_value = String(XMLString::transcode(attributes.getValue(0u))).trim();
			
			if (attribute_value == "peptide_significance_threshold" 
				|| (attribute_value == "significance_threshold" && inside_peptide_))
			{
				(*id_data_)[peptide_identification_index_].setSignificanceThreshold(
					((String) XMLString::transcode(attributes.getValue(1u))).toFloat());
			}
			else if (attribute_value == "significance_threshold" && !inside_peptide_)
			{
				(*protein_identifications_)[protein_identification_index_].setSignificanceThreshold(
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
				(*id_data_)[peptide_identification_index_].setMetaValue("RT", 
					((String) XMLString::transcode(attributes.getValue(1u))).toFloat());				
			}
			else if (attribute_value == "predicted_precursor_retention_time")
			{
				actual_peptide_hit_.setMetaValue("predicted_RT", 
					((String) XMLString::transcode(attributes.getValue(1u))).toFloat());				
			}
			else if (attribute_value == "number_of_identifications" || attribute_value == "number_of_db_searches")
			{
				id_data_->resize(((String) XMLString::transcode(attributes.getValue(1u))).toInt());
			}
			else if (attribute_value == "number_of_protein_identifications")
			{
				for(Int i = 0; i < ((String) XMLString::transcode(attributes.getValue(1u))).toInt(); i++)
				{
					Identification temp_protein_identification;
					protein_identifications_->push_back(temp_protein_identification);
				}
			}
			else if (attribute_value == "precursor_mz")
			{
				(*id_data_)[peptide_identification_index_].setMetaValue("MZ", 
					((String) XMLString::transcode(attributes.getValue(1u))).toFloat());				
			}
			else if (attribute_value == "predicted_rt_p_value")
			{			
				actual_peptide_hit_.setMetaValue("predicted_RT_p_value", 
					((String) XMLString::transcode(attributes.getValue(1u))).toFloat());				
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
					(*id_data_)[peptide_identification_index_].setScoreType(((String) XMLString::transcode(attributes.getValue(1u))));										
				}
				else
				{
					(*protein_identifications_)[protein_identification_index_].setScoreType(((String) XMLString::transcode(attributes.getValue(1u))));					
				}
			}
			else if (attribute_value == "identifier")
			{
				if (inside_peptide_)
				{
					(*id_data_)[peptide_identification_index_].setIdentifier(((String) XMLString::transcode(attributes.getValue(1u))));										
				}
				else
				{
					(*protein_identifications_)[protein_identification_index_].setIdentifier(((String) XMLString::transcode(attributes.getValue(1u))));					
				}
			}
			else if (attribute_value == "higher_score_better")
			{
				String temp_string;
				bool temp_bool;
				
				temp_string = (String) XMLString::transcode(attributes.getValue(1u));
				temp_bool = (temp_string == "true");

				if (inside_peptide_)
				{
					(*id_data_)[peptide_identification_index_].setHigherScoreBetter(temp_bool);
				}
				else
				{
					(*protein_identifications_)[protein_identification_index_].setHigherScoreBetter(temp_bool);
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
				//just for compatibility to older xml version begin
				
						
				UInt index = ((String) XMLString::transcode(attributes.getValue(1u))).toInt();
				
				if (index>=date_times_temp_.size())
				{
					throw Exception::ParseError(__FILE__,__LINE__,__PRETTY_FUNCTION__,String(index),"Undefined date_group index");
				}

				String date_time_string = date_times_temp_.at(index);
				DateTime temp_date_time;
				temp_date_time.set(date_time_string);
				(*protein_identifications_)[protein_identification_index_].setDateTime(temp_date_time);
				(*protein_identifications_)[protein_identification_index_].setIdentifier(date_time_string);
				if (inside_peptide_)
				{
					(*id_data_)[peptide_identification_index_].setIdentifier(date_time_string);
				}
				// end
			}						
			else if (attribute_value == "db")
			{
				actual_search_parameters_.db = ((String) XMLString::transcode(attributes.getValue(1u)));
			}
			else if (attribute_value == "db_version")
			{
				actual_search_parameters_.db_version = ((String) XMLString::transcode(attributes.getValue(1u)));
			}
			else if (attribute_value == "taxonomy")
			{
				actual_search_parameters_.taxonomy = ((String) XMLString::transcode(attributes.getValue(1u)));
			}
			else if (attribute_value == "charges")
			{
				actual_search_parameters_.charges = ((String) XMLString::transcode(attributes.getValue(1u)));
			}
			else if (attribute_value == "mass_type")
			{
				actual_search_parameters_.mass_type = (Identification::PeakMassType)((String) XMLString::transcode(attributes.getValue(1u))).toInt();
			}
			else if (attribute_value == "fixed_modifications")
			{							
				actual_search_parameters_.fixed_modifications.push_back(((String) XMLString::transcode(attributes.getValue(1u))));
			}
			else if (attribute_value == "variable_modifications")
			{			
				actual_search_parameters_.variable_modifications.push_back(((String) XMLString::transcode(attributes.getValue(1u))));
			}
			else if (attribute_value == "enzyme")
			{
				actual_search_parameters_.enzyme = (Identification::DigestionEnzyme)((String) XMLString::transcode(attributes.getValue(1u))).toInt();
			}
			else if (attribute_value == "missed_cleavages")
			{
				actual_search_parameters_.missed_cleavages = ((String) XMLString::transcode(attributes.getValue(1u))).toInt();
			}
			else if (attribute_value == "peak_mass_tolerance")
			{
				actual_search_parameters_.peak_mass_tolerance = ((String) XMLString::transcode(attributes.getValue(1u))).toFloat();
			}
			else if (attribute_value == "precursor_tolerance")
			{
				actual_search_parameters_.precursor_tolerance = ((String) XMLString::transcode(attributes.getValue(1u))).toFloat();
			}
			else if (attribute_value == "date")
			{
				DateTime temp_date;
				temp_date.set(((String) XMLString::transcode(attributes.getValue(1u))));
				(*protein_identifications_)[protein_identification_index_].setDateTime(temp_date);
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

  void IdXMLHandler::endElement(const XMLCh* const /*uri*/, const XMLCh* const /*local_name*/, const XMLCh* const qname)
 	{
 		tag_ = String(xercesc::XMLString::transcode(qname)).trim();
		
 		if (tag_ == "protein")
 		{	
 			(*protein_identifications_)[protein_identification_index_].insertHit(actual_protein_hit_);

	 		actual_protein_hit_ = ProteinHit();
	 		actual_peptide_indices_.clear();
			inside_protein_ = false;
			inside_global_protein_ = false;
 		}
 		else if (tag_ == "peptide")
 		{
			bool already_stored = false;
			vector<PeptideHit>::iterator  it;
 			
			vector<PeptideHit> temp_peptide_hits = 
				(*id_data_)[peptide_identification_index_].getHits();
				
			it = temp_peptide_hits.begin();
			while(it != temp_peptide_hits.end() && !already_stored)
			{
				if (it->getSequence() == actual_peptide_hit_.getSequence())
				{
					already_stored = true;
					if (inside_protein_)
					{
						// If the peptide is already stored and references 
						// a protein then add the reference.
						it->addProteinAccession(actual_protein_hit_.getAccession());						
					}
				}
				it++;
			}
			(*id_data_)[peptide_identification_index_].setHits(temp_peptide_hits);
			if (!already_stored)
			{
				if (inside_protein_)
				{
	 				String 							date_time;
	 			
					actual_peptide_hit_.addProteinAccession(actual_protein_hit_.getAccession());
				}
	 			(*id_data_)[peptide_identification_index_].insertHit(actual_peptide_hit_); 			
			}
 			actual_peptide_hit_ = PeptideHit();
 			inside_peptide_ = false;
 		}
 		
 		if (tag_ == "proteinGroup")
 		{
 			(*protein_identifications_)[protein_identification_index_].setSearchParameters(actual_search_parameters_);
 			actual_search_parameters_ = Identification::SearchParameters();
 		}

		tag_ = "";
 	} 

  void IdXMLHandler::characters(const XMLCh* const chars, unsigned int /*length*/)
  {
		if (tag_ == "dbID")
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
  
  void IdXMLHandler::writePeptideHit(ostream& os, 
  																				String shift,
			  																	PeptideHit hit,
			  																	Real significance_threshold,
			  																	UInt identification_index,
			  																	DataValue precursor_retention_time,
			  																	DataValue precursor_mz,
																					String identifier,
																					String score_type,
																					bool higher_score_better,
			  																	DoubleReal predicted_retention_time,
																					DoubleReal predicted_rt_p_value)
  {
  	String 																temp_peptide_sequence 		= "";
		
		if (identifier == "")
		{
//			throw Exception::InvalidValue(__FILE__, __LINE__, __PRETTY_FUNCTION__,"No identifier set for identification ",String(identification_index));
		}  			

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
		<<  shift << "\t<userParam name=\"significance_threshold\" value=\"" 
 		<< significance_threshold << "\"/>\n"
		<<  shift << "\t<userParam name=\"precursor_charge\" value=\"" 
		<< hit.getCharge() << "\" />\n"
		<< shift << "\t<userParam name=\"precursor_retention_time\" value=\""
		<< precursor_retention_time << "\" />\n";						
		if (predicted_rt_p_value != -1)
		{
			os << shift << "\t<userParam name=\"predicted_rt_p_value\" "
				 << "value=\"" 
				 << predicted_rt_p_value
				 << "\" />\n"
				 << shift << "\t<userParam name=\"predicted_precursor_retention_time\" "
				 << "value=\"" 
				 << predicted_retention_time
				 << "\" />\n";
		}

		os 	<< shift << "\t<userParam name=\"precursor_mz\" value=\""
				<< precursor_mz << "\" />\n"
				<< shift << "\t<userParam name=\"score\" value=\""
				<< hit.getScore() << "\" />\n"
				<< shift << "\t<userParam name=\"score_type\" value=\""
				<< score_type << "\" />\n"
				<< shift << "\t<userParam name=\"identifier\" value=\""
				<< identifier << "\" />\n"
				<< shift << "\t<userParam name=\"higher_score_better\" value=\""
				<< higher_score_better << "\" />\n"
				<<  shift << "</peptide>\n";
	}
  

	} // namespace Internal
} // namespace OpenMS
