// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2006 -- Oliver Kohlbacher, Knut Reinert
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
// $Id: AnalysisXMLHandler.C,v 1.14 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/AnalysisXMLHandler.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <iostream>
#include <map>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
  
  AnalysisXMLHandler::AnalysisXMLHandler(vector<ProteinIdentification>& protein_identifications,
  									 vector<Identification>& identifications, 
   									 vector<float>& precursor_retention_times, 
   									 vector<float>& precursor_mz_values) :
    XMLHandler(),
    protein_identifications_(&protein_identifications),
    identifications_(&identifications),
    precursor_retention_times_(&precursor_retention_times),
    precursor_mz_values_(&precursor_mz_values),
    actual_protein_hit_(),
    actual_protein_hits_(),
    actual_peptide_hit_(),
    actual_peptide_hits_(),
    peptide_identification_index_(0),
    protein_identification_index_(0),
    inside_peptide_(false),
    contact_person_(0),
    const_protein_identifications_(),
    const_identifications_(),
    const_precursor_retention_times_(),
    const_precursor_mz_values_(),
		const_contact_person_(),
		const_predicted_retention_times_(),
		tag_(),        
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    predicted_retention_times_(),
    predicted_sigma_(0),
    const_predicted_sigma_(0),
    date_times_temp_(),
    date_times_counter_(0)
  {
  }
   									 
  AnalysisXMLHandler::AnalysisXMLHandler(const vector<ProteinIdentification>& protein_identifications,
  									 const vector<Identification>& identifications, 
   									 const vector<float>& precursor_retention_times, 
   									 const vector<float>& precursor_mz_values) :
    XMLHandler(),
    protein_identifications_(),
    identifications_(0),
    precursor_retention_times_(0),
    precursor_mz_values_(0),
    actual_protein_hit_(),
    actual_protein_hits_(),
    actual_peptide_hit_(),
    actual_peptide_hits_(),
    peptide_identification_index_(0),
    protein_identification_index_(0),
    inside_peptide_(false),
    contact_person_(0),
    const_protein_identifications_(protein_identifications),
    const_identifications_(identifications),
    const_precursor_retention_times_(precursor_retention_times),
    const_precursor_mz_values_(precursor_mz_values),
		const_contact_person_(),
		const_predicted_retention_times_(),
		tag_(),        
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    predicted_retention_times_(),
    predicted_sigma_(0),
    const_predicted_sigma_(0),
    date_times_temp_(),
    date_times_counter_(0)    
  {
  }
   									 
  AnalysisXMLHandler::AnalysisXMLHandler(vector<ProteinIdentification>* protein_identifications, 
      									 vector<Identification>* identifications, 
      									 vector<float>* precursor_retention_times, 
      									 vector<float>* precursor_mz_values,
      									 ContactPerson* contact_person) :
    XMLHandler(),
    protein_identifications_(protein_identifications),
    identifications_(identifications),
    precursor_retention_times_(precursor_retention_times),
    precursor_mz_values_(precursor_mz_values),
    actual_protein_hit_(),
    actual_protein_hits_(),
    actual_peptide_hit_(),
    actual_peptide_hits_(),
    peptide_identification_index_(0),
    protein_identification_index_(0),
    inside_peptide_(false),
    contact_person_(contact_person),
    const_protein_identifications_(),
    const_identifications_(),
    const_precursor_retention_times_(),
    const_precursor_mz_values_(),
		const_contact_person_(),
		const_predicted_retention_times_(),
		tag_(),
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    predicted_retention_times_(),
    predicted_sigma_(0),
    const_predicted_sigma_(0),
    date_times_temp_(),
    date_times_counter_(0)
  {
  }

  AnalysisXMLHandler::AnalysisXMLHandler(vector<ProteinIdentification>* protein_identifications, 
      									 vector<Identification>* identifications, 
      									 vector<float>* precursor_retention_times, 
      									 vector<float>* precursor_mz_values,
      									 ContactPerson* contact_person,
      									 map<String, double>* predicted_retention_times,
      									 DoubleReal* predicted_sigma) :
    XMLHandler(),
    protein_identifications_(protein_identifications),
    identifications_(identifications),
    precursor_retention_times_(precursor_retention_times),
    precursor_mz_values_(precursor_mz_values),
    actual_protein_hit_(),
    actual_protein_hits_(),
    actual_peptide_hit_(),
    actual_peptide_hits_(),
    peptide_identification_index_(0),
    protein_identification_index_(0),
    inside_peptide_(false),
    contact_person_(contact_person),
    const_protein_identifications_(),
    const_identifications_(),
    const_precursor_retention_times_(),
    const_precursor_mz_values_(),
		const_contact_person_(),
		const_predicted_retention_times_(),
		tag_(),
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    predicted_retention_times_(predicted_retention_times),
    predicted_sigma_(predicted_sigma),
    const_predicted_sigma_(0),
    date_times_temp_(),
    date_times_counter_(0)
  {
  }

  AnalysisXMLHandler::AnalysisXMLHandler(const vector<ProteinIdentification>& protein_identifications, 
      									 const vector<Identification>& identifications, 
      									 const vector<float>& precursor_retention_times, 
      									 const vector<float>& precursor_mz_values,
      									 const ContactPerson& contact_person) :
    XMLHandler(),
    protein_identifications_(0),
    identifications_(0),
    precursor_retention_times_(0),
    precursor_mz_values_(0),
    actual_protein_hit_(),
    actual_protein_hits_(),
    actual_peptide_hit_(),
    actual_peptide_hits_(),
    peptide_identification_index_(0),
    protein_identification_index_(0),
    inside_peptide_(false),
    contact_person_(0),
    const_protein_identifications_(protein_identifications),
    const_identifications_(identifications),
    const_precursor_retention_times_(precursor_retention_times),
    const_precursor_mz_values_(precursor_mz_values),
		const_contact_person_(contact_person),
		const_predicted_retention_times_(),
		tag_(),
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    predicted_retention_times_(),
    predicted_sigma_(0),
		const_predicted_sigma_(0),
    date_times_temp_(),
    date_times_counter_(0)
  {
  }
      									 	      									 	  
  AnalysisXMLHandler::AnalysisXMLHandler(const vector<ProteinIdentification>& protein_identifications, 
      									 const vector<Identification>& identifications, 
      									 const vector<float>& precursor_retention_times, 
      									 const vector<float>& precursor_mz_values,
      									 const ContactPerson& contact_person,
      									 const map<String, double>& const_predicted_retention_times,
      									 DoubleReal predicted_sigma) :
    XMLHandler(),
    protein_identifications_(0),
    identifications_(0),
    precursor_retention_times_(0),
    precursor_mz_values_(0),
    actual_protein_hit_(),
    actual_protein_hits_(),
    actual_peptide_hit_(),
    actual_peptide_hits_(),
    peptide_identification_index_(0),
    protein_identification_index_(0),
    inside_peptide_(false),
    contact_person_(0),
    const_protein_identifications_(protein_identifications),
    const_identifications_(identifications),
    const_precursor_retention_times_(precursor_retention_times),
    const_precursor_mz_values_(precursor_mz_values),
		const_contact_person_(contact_person),
		const_predicted_retention_times_(const_predicted_retention_times),
		tag_(),
		charge_identification_index_(0),        
    inside_protein_(false),
    inside_global_protein_(false),
    predicted_retention_times_(),
    predicted_sigma_(),
    const_predicted_sigma_(predicted_sigma),
    date_times_temp_(),
    date_times_counter_(0)
  {
  }
      									 	      									 	  
  AnalysisXMLHandler::AnalysisXMLHandler(const AnalysisXMLHandler& source) :
    XMLHandler(source),
    protein_identifications_(source.protein_identifications_),
    identifications_(source.identifications_),
    precursor_retention_times_(source.precursor_retention_times_),
    precursor_mz_values_(source.precursor_mz_values_),
    actual_protein_hit_(source.actual_protein_hit_),
    actual_protein_hits_(source.actual_protein_hits_),
    actual_peptide_hit_(source.actual_peptide_hit_),
    actual_peptide_hits_(source.actual_peptide_hits_),
    peptide_identification_index_(source.peptide_identification_index_),
    protein_identification_index_(source.protein_identification_index_),
    inside_peptide_(source.inside_peptide_),
    contact_person_(source.contact_person_),
    const_protein_identifications_(source.const_protein_identifications_),
    const_identifications_(source.const_identifications_),
    const_precursor_retention_times_(source.const_precursor_retention_times_),
    const_precursor_mz_values_(source.const_precursor_mz_values_),
		const_contact_person_(source.const_contact_person_),
		const_predicted_retention_times_(source.const_predicted_retention_times_),
		tag_(),
		charge_identification_index_(source.charge_identification_index_),        
    inside_protein_(false),
    inside_global_protein_(false),
    predicted_retention_times_(source.predicted_retention_times_),
    predicted_sigma_(source.predicted_sigma_),                                    
    const_predicted_sigma_(source.const_predicted_sigma_),
    date_times_temp_(),
    date_times_counter_(0)                                    
  {

  }
   
  AnalysisXMLHandler::~AnalysisXMLHandler()
  {
    
  }
  
  AnalysisXMLHandler& AnalysisXMLHandler::operator = (const AnalysisXMLHandler& source)
  {
    if (&source == this) return *this;
    
    XMLHandler::operator=(source);
    identifications_ = source.identifications_;
    precursor_retention_times_ = source.precursor_retention_times_;
    precursor_mz_values_ = source.precursor_mz_values_;
    contact_person_ = source.contact_person_;
    predicted_retention_times_ = source.predicted_retention_times_;
    
    return *this;
  }

  bool AnalysisXMLHandler::operator == (const AnalysisXMLHandler& rhs) const
  {
    return 
      ( identifications_ == rhs.identifications_ ) &&
      ( precursor_retention_times_ == rhs.precursor_retention_times_ ) &&
      ( precursor_mz_values_ == rhs.precursor_mz_values_ ) &&
      ( contact_person_ == rhs.contact_person_ ) && 
      ( predicted_retention_times_ == rhs.predicted_retention_times_ );
  }

  bool AnalysisXMLHandler::operator != (const AnalysisXMLHandler& rhs) const
  {
    return !(operator == (rhs));
  }
   
  void AnalysisXMLHandler::writeTo(std::ostream& os)
  {
  	vector<ProteinHit> 										temp_protein_hits;
  	vector<PeptideHit> 										temp_peptide_hits;
  	vector<PeptideHit>* 									referencing_peptide_hits;
  	vector<PeptideHit>* 									non_referencing_peptide_hits;
  	vector<PeptideHit>  									all_peptide_hits;
  	vector<UnsignedInt> 									indices;
  	vector<String>												date_time_strings;
  	multimap< String, ProteinHit > 				all_protein_hits;
  	
  	bool 																	protein_hits_present = false;
  	String 																temp_peptide_sequence;
  	map< String , UnsignedInt>						date_times;
  	UnsignedInt 													date_times_counter = 0;
  	String 																date_time_string;
  	DateTime															date_time;
  	map< String , UnsignedInt>::iterator	date_times_iterator;
		UnsignedInt 													group_count = 0;
		
		os << "<!-- -*- Mode: XML; tab-width: 2; -*- -->\n"
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
		os << "</sampleName>\n\t\t\t<contact>\n"
		 << "\t\t\t\t<name>";
		if (const_contact_person_.getName() != "")
		{
			os << const_contact_person_.getName();			
		}
		else
		{
			os << "unknown";
		}
		os << "</name>\n"		
		 << "\t\t\t\t<institution>";
		if (const_contact_person_.getInstitution() != "")
		{
		  os << const_contact_person_.getInstitution();
		}
		else
		{
			os << "unknown";			
		}		 
		os << "</institution>\n";	
		if (const_contact_person_.getContactInfo() != "")
		{
			os << "\t\t\t\t<contactInfo>"
				 << const_contact_person_.getContactInfo()
				 << "</contactInfo>\n";
		}
		os << "\t\t\t</contact>\n"	
		 << "\t\t</admin>\n"
		 << "\t\t<settings>\n"
		 << "\t\t\t<userParam name=\"number_of_identifications\" value=\""
  	 << const_identifications_.size() << "\" />\n"
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
  	for(vector<Identification>::const_iterator it = const_identifications_.begin();
  		  it != const_identifications_.end();
  		  it++)
  	{
  		date_time = it->getDateTime();
  		date_time.get(date_time_string);
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
					for(UnsignedInt i = 0; i < const_identifications_.size(); i++)
					{
						referencing_peptide_hits = const_identifications_[i].getReferencingHits(date_time_string, 
																																	 								  protein_hits_it->getAccession());
						for(vector<PeptideHit>::iterator peptide_hits_it = referencing_peptide_hits->begin();
								peptide_hits_it != referencing_peptide_hits->end();
								peptide_hits_it++)
						{							
							writePeptideHit(os, 
  														String("\t\t\t\t"),
  														*peptide_hits_it,
  														const_identifications_[i].getPeptideSignificanceThreshold(),
  														i,
  														const_identifications_[i].getCharge(), 
  														const_precursor_retention_times_[i],
  														const_precursor_mz_values_[i],
  														const_identifications_[i].getDateTime(),
															date_times);

						} // peptide_hits_it
					} // identifications
					os	<< "\t\t\t\t<userParam name=\"protein_identification_index\" value=\"" 
					   << j << "\"/>\n"
					   << "\t\t\t\t<userParam name=\"date_group_index\" value=\"" 
					   << getDateGroupIndex(date_time, date_times) << "\"/>\n"
					   << "\t\t\t\t<userParam name=\"score\" value=\"" 
					   << protein_hits_it->getScore() << "\"/>\n"
					    << "\t\t\t</protein>\n";
				} // protein hits
				os << "\t\t\t<userParam name=\"proteins\" value=\"all\"/>\n"
					 << "\t\t</proteinGroup>\n";
			} // protein identifications
		}

		for(UnsignedInt i = 0; i < const_identifications_.size(); i++)
		{
			if (const_identifications_[i].getProteinHits().size() > 0)
			{
				protein_hits_present = true;
			}
		}

		if (protein_hits_present)
		{
			os << "\t\t<proteinGroup count=\"" << group_count << "\">\n";
			for(UnsignedInt j = 0; j < const_identifications_.size(); j++)
			{
				date_time = const_identifications_[j].getDateTime();
				date_time.get(date_time_string);
					
				temp_protein_hits = const_identifications_[j].getProteinHits();
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
					for(UnsignedInt i = 0; i < const_identifications_.size(); i++)
					{
						referencing_peptide_hits = const_identifications_[i].getReferencingHits(date_time_string, 
																																						 protein_hits_it->getAccession());
						for(vector<PeptideHit>::iterator peptide_hits_it = referencing_peptide_hits->begin();
								peptide_hits_it != referencing_peptide_hits->end();
								peptide_hits_it++)
						{							
							writePeptideHit(os, 
  														String("\t\t\t\t"),
  														*peptide_hits_it,
  														const_identifications_[i].getPeptideSignificanceThreshold(),
  														i,
  														const_identifications_[i].getCharge(), 
  														const_precursor_retention_times_[i],
  														const_precursor_mz_values_[i],
  														const_identifications_[i].getDateTime(),
															date_times);

						} // peptide_hits_it
					} // identifications
					os  << "\t\t\t<userParam name=\"identification_index\" value=\""
							<< j << "\"/>\n"
					   << "\t\t\t<userParam name=\"date_group_index\" value=\"" 
					   << getDateGroupIndex(date_time, date_times) << "\"/>\n"
					   << "\t\t\t\t<userParam name=\"score\" value=\"" 
					   << protein_hits_it->getScore() << "\"/>\n"
					   << "\t\t\t</protein>\n";
				} // protein hits
			} // identifications
			os << "\t\t\t<userParam name=\"proteins\" value=\"all\"/>\n"
				 << "\t\t</proteinGroup>\n";
		}

		// All protein hits with the corresponding date and time are stored.
		// They are used afterwards to determine the peptide hits that do
		// not reference any protein hits.
		for(UnsignedInt i = 0; i < const_identifications_.size(); i++)
		{
			const_identifications_[i].getDateTime().get(date_time_string);
			for(vector<ProteinHit>::const_iterator it = const_identifications_[i].getProteinHits().begin();
					it != const_identifications_[i].getProteinHits().end();
					it++)
			{
				all_protein_hits.insert(make_pair(date_time_string, *it));				
			}
			// if the date time is new add it to the date times vector
			if (find(date_time_strings.begin(), date_time_strings.end(), date_time_string) == date_time_strings.end())
			{
				date_time_strings.push_back(date_time_string);
			}
		}
		for(UnsignedInt i = 0; i < const_protein_identifications_.size(); i++)
		{
			const_protein_identifications_[i].getDateTime().get(date_time_string);
			for(vector<ProteinHit>::const_iterator it = const_protein_identifications_[i].getProteinHits().begin();
					it != const_protein_identifications_[i].getProteinHits().end();
					it++)
			{
				all_protein_hits.insert(make_pair(date_time_string, *it));				
			}
			// if the date time is new add it to the date times vector
			if (find(date_time_strings.begin(), date_time_strings.end(), date_time_string) == date_time_strings.end())
			{
				date_time_strings.push_back(date_time_string);
			}
		}
		for(UnsignedInt i = 0; i < const_identifications_.size(); i++)
		{
			non_referencing_peptide_hits = 
				const_identifications_[i].getNonReferencingHits(all_protein_hits, 
				 								  														  date_time_strings);
				 								  														  
																  														 
			for(vector<PeptideHit>::const_iterator peptide_hits_it = non_referencing_peptide_hits->begin();
					peptide_hits_it != 	non_referencing_peptide_hits->end();
					peptide_hits_it++)
			{	
				writePeptideHit(os, 
												String("\t\t"),
												*peptide_hits_it,
												const_identifications_[i].getPeptideSignificanceThreshold(),
												i,
												const_identifications_[i].getCharge(), 
												const_precursor_retention_times_[i],
												const_precursor_mz_values_[i],
												const_identifications_[i].getDateTime(),
												date_times);
			}
		}

		os << "\t</ident>\n"
		<< "</mzIdent>";				 
		 
  }

  bool AnalysisXMLHandler::startElement(const QString & /* uri */, 
  																			const QString & /* local_name */,
																				const QString & qname, 
																				const QXmlAttributes & attributes )
	{
		tag_ = String(qname.ascii());

		String attribute_value;

		if (tag_ == "userParam")
		{
			attribute_value = String(attributes.value(0).ascii()).trim();
			if (attribute_value == "peptide_significance_threshold")
			{
				(*identifications_)[peptide_identification_index_].setPeptideSignificanceThreshold(
					((String) attributes.value(1).ascii()).toFloat());
			}
			else if (attribute_value == "precursor_charge")
			{
				(*identifications_)[peptide_identification_index_].setCharge(((String) attributes.value(1).ascii()).toInt());
			}
			else if (attribute_value == "precursor_retention_time")
			{
				(*precursor_retention_times_)[peptide_identification_index_] = 
					((String) attributes.value(1).ascii()).toFloat();
			}
			else if (attribute_value == "predicted_sigma")
			{
				if (predicted_sigma_ != 0)
				{
					*predicted_sigma_ = ((String) attributes.value(1).ascii()).toFloat();
				}
			}
			else if (attribute_value == "predicted_precursor_retention_time")
			{
				if (predicted_sigma_ != 0)
				{
					predicted_retention_times_->insert(make_pair(actual_peptide_hit_.getSequence(), 
						((String) attributes.value(1).ascii()).toFloat()));
				}
			}
			else if (attribute_value == "number_of_identifications")
			{
				for(int i = 0; i < ((String) attributes.value(1).ascii()).toInt(); i++)
				{
					Identification temp_identification;
					identifications_->push_back(temp_identification);
					(*precursor_retention_times_).push_back(0);
					(*precursor_mz_values_).push_back(0);
				}
			}
			else if (attribute_value == "number_of_protein_identifications")
			{
				for(int i = 0; i < ((String) attributes.value(1).ascii()).toInt(); i++)
				{
					ProteinIdentification temp_protein_identification;
					protein_identifications_->push_back(temp_protein_identification);
				}
			}
			else if (attribute_value == "precursor_mz")
			{
				(*precursor_mz_values_)[peptide_identification_index_] = ((String) attributes.value(1).ascii()).toFloat();
			}
			else if (attribute_value == "score")
			{
				if (inside_peptide_)
				{
					actual_peptide_hit_.setScore(((String) attributes.value(1).ascii()).toFloat());
				}
				else
				{
					actual_protein_hit_.setScore(((String) attributes.value(1).ascii()).toFloat());					
				}
			}
			else if (attribute_value == "identification_index")
			{
				if (inside_peptide_)
				{
					peptide_identification_index_ = ((String) attributes.value(1).ascii()).toInt();
				}
				else
				{
					protein_identification_index_ = ((String) attributes.value(1).ascii()).toInt();
				}
			}
			else if (attribute_value == "protein_identification_index")
			{
				protein_identification_index_ = ((String) attributes.value(1).ascii()).toInt();
				inside_global_protein_ = true;
			}
			else if (attribute_value.hasPrefix("date_group_") 
							 && attribute_value != "date_group_index")
			{
				date_times_temp_.push_back((String) attributes.value(1).ascii());
				date_times_counter_++;
			}
			else if (attribute_value == "date_group_index")
			{		
					String   		date_time_string;
					
					date_time_string = date_times_temp_.at(((String) attributes.value(1).ascii()).toInt());
					if (inside_peptide_)
					{
						(*identifications_)[peptide_identification_index_].getDateTime().set(date_time_string);
					}
					else
					{
						if (inside_global_protein_)
						{
							(*protein_identifications_)[protein_identification_index_].getDateTime().set(date_time_string);
						}
						else
						{
							(*identifications_)[protein_identification_index_].getDateTime().set(date_time_string);
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
			inside_protein_ = true;
			inside_peptide_ = false;			
		}
		
		return true;

	}
	  
  bool AnalysisXMLHandler::endElement(const QString & /* uri */, 
  																		const QString & /* local_name */,
 								  										const QString & qname)
 	{
 		tag_ = String(qname.ascii()).trim();
 		
 		if (tag_ == "protein")
 		{	
			if (inside_global_protein_)
			{
	 			(*protein_identifications_)[protein_identification_index_].insertProteinHit(actual_protein_hit_);
			}
			else
			{				
	 			(*identifications_)[protein_identification_index_].insertProteinHit(actual_protein_hit_);
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
				(*identifications_)[peptide_identification_index_].getPeptideHits();
				
			it = temp_peptide_hits.begin();
			while(it != temp_peptide_hits.end() && !already_stored)
			{
				if (it->getSequence() == actual_peptide_hit_.getSequence())
				{
					already_stored = true;
					if (inside_protein_)
					{
		 				String 							date_time;
		 			
		 				(*identifications_)[peptide_identification_index_].getDateTime().get(date_time);
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
	 			
	 				(*identifications_)[peptide_identification_index_].getDateTime().get(date_time);
					actual_peptide_hit_.addProteinIndex(make_pair(date_time, actual_protein_hit_.getAccession()));
				}
	 			(*identifications_)[peptide_identification_index_].insertPeptideHit(actual_peptide_hit_); 			
			}
 			actual_peptide_hit_.clear();
 			inside_peptide_ = false;
 		}

		tag_ = "";
 		
 		return true;
 		
 	} 

  bool AnalysisXMLHandler::characters( const QString & chars )
  {
		if (tag_ == "name")
		{
			contact_person_->setName(((String) chars.ascii()).trim());
			tag_ = "";
		}
		else if (tag_ == "institution")
		{
			contact_person_->setInstitution(((String) chars.ascii()).trim());
			tag_ = "";
		}
		else if (tag_ == "contactInfo")
		{
			contact_person_->setContactInfo(((String) chars.ascii()).trim());
			tag_ = "";
		}
		else if (tag_ == "dbName")
		{
			actual_protein_hit_.setAccessionType(((String) chars.ascii()).trim());
			tag_ = "";
		}
		else if (tag_ == "dbID")
		{
			actual_protein_hit_.setAccession(((String) chars.ascii()).trim());			
			tag_ = "";
		}
		else if (tag_ == "sequence")
		{
			actual_protein_hit_.setSequence(((String) chars.ascii()).trim());			
			tag_ = "";
		}
		else if (tag_ == "seq")
		{
			actual_peptide_hit_.setSequence(((String) chars.ascii()).trim());			
			tag_ = "";
		}
  	
		return true;
  }
  
  UnsignedInt AnalysisXMLHandler::getDateGroupIndex(DateTime 											date_time,
			  																						map< String , UnsignedInt> 		date_times)
	{
		String 																date_time_string 					= "";
		UnsignedInt 													date_group_index					= 0;
		map< String , UnsignedInt>::iterator 	date_times_iterator;
		
		// determining the date group
		date_time.get(date_time_string);
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
				<< shift << "\t<userParam name=\"date_group_index\" value=\""
				<< date_group_index << "\"/>\n"
				<<  shift << "</peptide>\n";
	}
  

	} // namespace Internal
} // namespace OpenMS
