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
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/HANDLERS/MascotXMLHandler.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	namespace Internal
	{
  
  MascotXMLHandler::MascotXMLHandler(ProteinIdentification* protein_identification,
								  									 vector<Identification>* identifications, 
								   									 vector<float>* precursor_retention_times, 
								   									 vector<float>* precursor_mz_values) :
    XMLHandler(),
    protein_identification_(protein_identification),
    identifications_(identifications),
    precursor_retention_times_(precursor_retention_times),
    precursor_mz_values_(precursor_mz_values),
    actual_protein_hit_(),
    actual_protein_hits_(),
    actual_peptide_hit_(),
    actual_peptide_hits_(),
    peptide_identification_index_(0),
    protein_identification_index_(0),
    const_protein_identification_(),
    const_identifications_(),
    const_precursor_retention_times_(),
    const_precursor_mz_values_(),
		tag_(),
		date_()        
  {
  }
   									       									 	      									 	  
  MascotXMLHandler::MascotXMLHandler(const MascotXMLHandler& source) :
    XMLHandler(source),
    protein_identification_(source.protein_identification_),
    identifications_(source.identifications_),
    precursor_retention_times_(source.precursor_retention_times_),
    precursor_mz_values_(source.precursor_mz_values_),
    actual_protein_hit_(source.actual_protein_hit_),
    actual_protein_hits_(source.actual_protein_hits_),
    actual_peptide_hit_(source.actual_peptide_hit_),
    actual_peptide_hits_(source.actual_peptide_hits_),
    peptide_identification_index_(source.peptide_identification_index_),
    protein_identification_index_(source.protein_identification_index_),
		tag_(),
		date_()
  {

  }
   
  MascotXMLHandler::~MascotXMLHandler()
  {
    
  }
  
  MascotXMLHandler& MascotXMLHandler::operator = (const MascotXMLHandler& source)
  {
    if (&source == this) return *this;
    
    XMLHandler::operator=(source);
    identifications_ = source.identifications_;
    precursor_retention_times_ = source.precursor_retention_times_;
    precursor_mz_values_ = source.precursor_mz_values_;
    
    return *this;
  }

  bool MascotXMLHandler::operator == (const MascotXMLHandler& rhs) const
  {
    return 
      ( identifications_ == rhs.identifications_ ) &&
      ( precursor_retention_times_ == rhs.precursor_retention_times_ ) &&
      ( precursor_mz_values_ == rhs.precursor_mz_values_ );
  }

  bool MascotXMLHandler::operator != (const MascotXMLHandler& rhs) const
  {
    return !(operator == (rhs));
  }
   
  bool MascotXMLHandler::startElement(const QString & /* uri */, 
  																			const QString & /* local_name */,
																				const QString & qname, 
																				const QXmlAttributes & attributes )
	{

		tag_ = String(qname.ascii());
		
		if (tag_ == "protein")
		{
			String attribute_value;
			attribute_value = String(attributes.value(0).ascii()).trim();
  		actual_protein_hit_.setAccession(attribute_value);
		}		
		else 
		{
			if (tag_ == "peptide")
			{
				String attribute_value;
				attribute_value = String(attributes.value(0).ascii()).trim();
	  		peptide_identification_index_ = attribute_value.toInt() - 1;
			}
			else if (tag_ == "u_peptide")
			{
				String attribute_value;
				attribute_value = String(attributes.value(0).ascii()).trim();
	  		peptide_identification_index_ = attribute_value.toInt() - 1;
			}
  		if (peptide_identification_index_ > identifications_->size())
  		{
				throw Exception::ParseError(__FILE__, __LINE__
																		, __PRETTY_FUNCTION__, ".mascotXML", 
																		"No header information present: use "
																		"the show_header=1 option in the "
																		"./export_dat.pl script");  			
  		}			
		}
		
		return true;

	}
	  
  bool MascotXMLHandler::endElement(const QString & /* uri */, 
  																	const QString & /* local_name */,
 								  									const QString & qname)
 	{
 		tag_ = String(qname.ascii()).trim();
 		 
 		if (tag_ == "protein")
 		{	
 			/// since Mascot uses SwissProt IDs we set this type here
			actual_protein_hit_.setAccessionType("SwissProt");
			actual_protein_hit_.setScoreType("Mascot");
 			protein_identification_->insertProteinHit(actual_protein_hit_);
 			actual_protein_hit_.clear();
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
				}
				it++;
			}
			if (!already_stored)
			{
				actual_peptide_hit_.setScoreType("Mascot");
				actual_peptide_hit_.addProteinIndex(make_pair(date_time_string_, actual_protein_hit_.getAccession()));
	 			(*identifications_)[peptide_identification_index_].insertPeptideHit(actual_peptide_hit_); 			
			}
			else
			{
				it--;
				it->addProteinIndex(make_pair(date_time_string_, actual_protein_hit_.getAccession()));
			}
 			actual_peptide_hit_.clear();
 		}
 		else if (tag_ == "u_peptide")
 		{
 			(*identifications_)[peptide_identification_index_].insertPeptideHit(actual_peptide_hit_); 			
 			actual_peptide_hit_.clear();
 		}
		tag_ = "";
 		
 		return true;
 		
 	} 

  bool MascotXMLHandler::characters( const QString & chars )
  {

		if (tag_ == "NumQueries")
		{
			Identification temp_identification;
			
			temp_identification.setDateTime(date_);
			for(int i = 0; i < ((String) chars.ascii()).trim().toInt(); i++)
			{
				identifications_->push_back(temp_identification);
				(*precursor_retention_times_).push_back(0);
				(*precursor_mz_values_).push_back(0);
			}
			tag_ = "";
		}
		else if (tag_ == "prot_score")
		{
			actual_protein_hit_.setScore(((String) chars.ascii()).trim().toInt());
		}
		else if (tag_ == "pep_exp_mz")
		{
			(*precursor_mz_values_)[peptide_identification_index_] = ((String) chars.ascii()).trim().toFloat();
			tag_ = "";
		}
		else if (tag_ == "pep_exp_z")
		{
			(*identifications_)[peptide_identification_index_].setCharge(((String) chars.ascii()).trim().toInt());
			tag_ = "";
		}
		else if (tag_ == "pep_score")
		{
			actual_peptide_hit_.setScore(((String) chars.ascii()).trim().toFloat());
			tag_ = "";
		}
		else if (tag_ == "pep_homol")
		{			
			(*identifications_)[peptide_identification_index_].setPeptideSignificanceThreshold(
					((String) chars.ascii()).trim().toFloat());
			tag_ = "";
		}
		else if (tag_ == "pep_ident")
		{
			DoubleReal temp_homology = 0;
			DoubleReal temp_identity = 0;
			
			/// According to matrixscience the homology threshold is only used if it exists and is
			/// smaller than the identity threshold.
			temp_homology = 
				(*identifications_)[peptide_identification_index_].getPeptideSignificanceThreshold();
			temp_identity = ((String) chars.ascii()).trim().toFloat();
			if (temp_homology > temp_identity || temp_homology == 0)
			{
				(*identifications_)[peptide_identification_index_].setPeptideSignificanceThreshold(
					temp_identity);				
			}
			tag_ = "";
		}
		else if (tag_ == "pep_seq")
		{
			actual_peptide_hit_.setSequence(((String) chars.ascii()).trim());
			tag_ = "";
		}
		else if (tag_ == "Date")
		{	
			vector< String > parts;
			
			((String) chars.ascii()).trim().split('T', parts);
			if (parts.size() == 2)
			{
				date_.set(parts[0] + ' ' + parts[1].prefix('Z'));
				date_time_string_ = parts[0] + ' ' + parts[1].prefix('Z');
			}
			protein_identification_->setDateTime(date_);
		}
		return true;
  }

	} // namespace Internal
} // namespace OpenMS
