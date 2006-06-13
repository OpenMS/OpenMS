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
// $Id: MascotOutfile.C,v 1.22 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/TimeStamp.h>
#include <OpenMS/FORMAT/MascotOutfile.h>
#include <OpenMS/FORMAT/TextFile.h>

#include <cassert>
#include <string>
#include <sstream>
#include <algorithm>
#include <cmath>

using namespace std;

namespace OpenMS 
{

  MascotOutfile::MascotOutfile(const string& filename, Real p)  
  	throw (Exception::ParseError)
    : db_searches_(), 
    	peptide_hits_(), 
    	protein_hits_(),
    	precursor_retention_times_(),
    	precursor_mz_values_(), 
    	ok_(false)
  {
  	TextFile f(filename);
  	vector<PeptideHit>::iterator peptide_hit_iterator;
		vector<String> parts;
		UnsignedInt number_of_queries = 0;
		map<UnsignedInt, UnsignedInt> indices;
		map<UnsignedInt, UnsignedInt>::iterator indices_iterator;
		int temp_int;
		Identification temp_db_search;
		SignedInt temp_charge = 0;
		String temp_identifier = "";
		vector<Real> temp_scores;
		map<String, vector<Real> > protein_map;
		Real temp_score = 0;
		Real temp_value = 0;
		Real temp_significance_threshold = 0;
		UnsignedInt index = 0;
		String::SizeType tag_start;
		String::SizeType tag_end;
	
		if (f.size() == 0)
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
				"File is empty!", filename);
		}
	
		///Mascot search was not successfull
		if (f.size()<5) return;
		  	
  	TextFile::iterator it;

  	/// (1.0) parse for retention time
  	it = f.search("_RETENTION_TIME=");
  	if (it==f.end())
  	{
	  	it = f.search("sequences=");			
  	}
  	else
  	{
 			precursor_retention_times_.push_back(it->suffix('=').trim().toFloat());
  	}
  	/// (1.0) parse for date

  	it = f.search(it, "date=");
  	if (it==f.end())
  	{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
  			"date in header section not found!" ,filename);
  	}
		PreciseTime precise_date(it->suffix('=').trim().toInt(),0);
		stringstream ss;
		ss << precise_date;
		
		DateTime date;
		it = f.search(it, "time=");
  	if (it==f.end())
  	{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
  			"time in header section not found!" ,filename);
  	}
		
		
		date.set(ss.str().substr(6,2) + "." + ss.str().substr(4,2) + "." + 
			ss.str().substr(0,4) + " " + it->suffix('=').trim());
		temp_db_search.setDateTime(date);
		
  	/// (1.0.1) parse for number of queries
  	it = f.search(it, "queries=");
  	if (it==f.end())
  	{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
  			"number of queries not found!" ,filename);
  	}
		number_of_queries = it->suffix('=').trim().toInt();

		/// (1.0.2) Searching for query indices for which peptides are present
		if (number_of_queries > 1)
		{
			for(UnsignedInt i = 1; i <= number_of_queries; i++)
			{
	  		it = f.search(it, "q" + String(i) + "_p1=");
	 		 	if (it!=f.end())
	  		{				  		
	  			temp_int = it->suffix('=').trim().toInt();
	  			if (temp_int != -1)
	  			{
	  				indices.insert(make_pair(i, indices.size()));
	  			}
	  		}
	  		else
	  		{
	  			i = number_of_queries;
	  		}
	  	}
	  }
	  else
	  {
	  	indices.insert(make_pair(1, 0));
	  }
	  
		/// (1.1) parse for precursor values
		for(indices_iterator = indices.begin(); 
				indices_iterator != indices.end();
				indices_iterator++)
		{
			it = f.search("qexp" + String(indices_iterator->first) + "=");
			if (it == f.end())
			{
 				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
 				"precursor information for query " + String(indices_iterator->first)
 				+ " not found!" ,filename);			
 			}
			it->suffix('=').split(',',parts);
			precursor_mz_values_.push_back(parts[0].toFloat());
			temp_charge = String(parts[1].trim()[0]).toInt();
			if (String(parts[1].trim()[1]) == "+")
			{
				temp_db_search.setCharge(temp_charge);
			}
			else
			{
				temp_db_search.setCharge((-1 * temp_charge));				
			}
			parts.clear();
			db_searches_.push_back(temp_db_search);
		}
		
		/// (1.2) parse for peptide significance threshold
		for(indices_iterator = indices.begin(); 
				indices_iterator != indices.end();
				indices_iterator++)
		{
			it = f.search("qplughole" + String(indices_iterator->first) + "=");
			if (it == f.end())
			{
	  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
	  			"significance threshold for query " + String(indices_iterator->first)
 				+ " in summary section not found!" ,filename);			
			}
			db_searches_[indices_iterator->second].setPeptideSignificanceThreshold(
				it->suffix('=').trim().toFloat());
		}
		for(indices_iterator = indices.begin(); 
				indices_iterator != indices.end();
				indices_iterator++)
		{
			it = f.search("qmatch" + String(indices_iterator->first) + "=");
			if (it == f.end())
			{
	  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
	  			"number of possible matches in the precursor window for query " + String(indices_iterator->first)
 				+ " in summary section not found!" ,filename);			
			}
			
			temp_value = it->suffix('=').trim().toFloat();
			temp_value = 10 * log10(temp_value / p / 20);
			if (temp_value 
						< db_searches_[indices_iterator->second].getPeptideSignificanceThreshold())
			{
				db_searches_[indices_iterator->second].setPeptideSignificanceThreshold(temp_value);
			}
		}
		
		
		/// The protein significance threshold is left zero because it is not
		/// stored directly in the Mascot outfile and information about the 
		/// calculation is, to our knowledge, not publicly available.				
		
		/// (2.1) parse for ProteinHit information (MudPIT scoring)
		
		if (number_of_queries > 1000)
		{
	  	it = f.searchSuffix("\"proteins\"", true);
	  	
			if (it == f.end())
			{
				cout << "no \"proteins\" tag found " << endl;
			}

	  	///Go to first protein hit entry
			if (it != f.end() && (it + 1) != f.end())
			{
				it += 2;
			}
	  	while(it != f.end())
	  	{
	  		/// search for first "
	  		if ((tag_start = (*it).find('"')) != string::npos)
	  		{
	  			/// search for second "
		  		if ((tag_end = (*it).find('"', tag_start + 1)) != string::npos)
		  		{
		  			temp_identifier = (*it).substr(tag_start + 1, tag_end - tag_start - 1);
						temp_scores.clear();
						for(UnsignedInt k = 0; k < 3; k++)
						{
							temp_scores.push_back(0);
						}
						protein_map.insert(make_pair(temp_identifier, temp_scores));
						it++;
		  		}	  			
					else
					{
						it = f.end();
					}
	  		}
				else
				{
					it = f.end();
				}
			}							
		}

	  UnsignedInt i = 1; /// first index of the peptidehits
	  UnsignedInt j = 1; /// second index of the peptidehits
  	/// (2.2) parse for PeptideHit information
		for(indices_iterator = indices.begin(); 
				indices_iterator != indices.end();
				indices_iterator++)
		{
			i = indices_iterator->first;
			j = 1;
			UnsignedInt counter = 1; ///counter of the peptidehits
	  	it = f.search(String("q")+String(i)+"_p" + String(j) + "=");
	  	while(it != f.end())
	  	{
	  		PeptideHit hit;
				
	  		///(2.1) parse for peptide sequence
				it->suffix('=').split(',',parts);
	  		hit.setSequence(parts[4]);
				temp_score = parts[7].toFloat();			
				hit.setScore(temp_score);
				hit.setScoreType("Mascot");
	  		
	  		hit.setRank(counter);
	   		///(2.3) insert into hits vector
	  		if (temp_score > 0)
	  		{
	  			db_searches_[indices_iterator->second].insertPeptideHit(hit);
					counter++;
	  		}
	  		
	  		if (number_of_queries > 1000)
	  		{
	  			temp_significance_threshold = 
	  				db_searches_[indices_iterator->second].getPeptideSignificanceThreshold();
	
		  		if (temp_score > temp_significance_threshold)
		  		{
		  			
						it->suffix('=').split('"',parts);
			  		index = 1;
			  		while(index < parts.size() - 1)
			  		{
			  			temp_scores = protein_map[parts[index]];
			  					  			
			  			/// the score of the protein hit
			  			temp_scores[0] += (temp_score - temp_significance_threshold);
			  			/// sum of the used significance thresholds
			  			temp_scores[1] += temp_significance_threshold;
			  			/// number of significance thresholds used in sum
			  			temp_scores[2] = temp_scores[2] + 1;
			  			
			  			protein_map[parts[index]] = temp_scores;
			  					  			
			  			index += 2;
			  		}
		  		}
	  		}
	  		
	  		///search for the next hit
	  		++j;
	  		it = f.search(it,String("q")+String(i)+"_p" + String(j) + "=");
				parts.clear();	
	  	}
		}

  	///(3) search for protein hit information
 		i = 1;
 		j = 1;
 		if (number_of_queries == 1)
 		{  	
	  	it = f.search(String("h")+String(i) + "=");
	  	peptide_hits_ = db_searches_[0].getPeptideHits();
	  	while(it != f.end())
	  	{
				ProteinHit protein_hit;
				String temp_peptide_sequence;
				int peptide_index = -1;
	
				protein_hit.setAccession(it->suffix('=').prefix(','));
				protein_hit.setAccessionType("SwissProt");
				protein_hit.setScore(it->substr(it->find(',')+1).prefix(',').toFloat());
				protein_hit.setScoreType("Mascot");
				protein_hit.setRank(i);
				
	  		it = f.search(it,String("h")+String(i)+"_q" + String(j) + "=");	
				if (it==f.end() && j == 1)
		  	{
		  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
		  			String("Line starting with 'h")+String(i)+"_q1=' not found!" ,
		  			filename);
		  	}
		  	/// search for the peptide hits belonging to the actual protein hit
		  	while(it != f.end())
		  	{ 	
					vector<String> parts;
					it->suffix('=').split(',',parts);
					temp_peptide_sequence = parts[6];
					
					for(uint index = 0; index < peptide_hits_.size(); index++)
					{
						if (peptide_hits_[index].getSequence() == temp_peptide_sequence)
						{
							peptide_index = index;
						}
					}
					/// setting of the indices to store the relational information 
					if (peptide_index != -1)
					{
//						peptide_hits_[peptide_index].addProteinIndex(i - 1);
//						protein_hit.addPeptideIndex(peptide_index);
						peptide_index = -1;
					}
									
					j++;
		  		it = f.search(it,String("h")+String(i)+"_q" + String(j) + "=");	
				}
				protein_hits_.push_back(protein_hit);
				i++;
				j = 1;
		  	it = f.search(String("h")+String(i) + "=");
			}						  	
			curr_peptide_hit_ = peptide_hits_.begin();
			curr_protein_hit_ = protein_hits_.begin();
			db_searches_[0].setPeptideAndProteinHits(peptide_hits_, protein_hits_);
		}
		
		for(indices_iterator = indices.begin(); 
				indices_iterator != indices.end();
				indices_iterator++)
		{
  		it = f.searchSuffix(String("\"query")+String(indices_iterator->first) + String("\""), true);	
			if (it==f.end())
	  	{
				precursor_retention_times_.push_back(0.f);
	  	}
			else
			{
				it = f.search(it, String("rtinseconds="));	
				if (it==f.end())
		  	{
		  		it = f.search(String("\"query")+String(indices_iterator->first) + String("\""), true);
				}
				else
				{
					precursor_retention_times_.push_back(it->suffix('=').trim().toFloat());
				}
			}				
		}
		
		for(map<String, vector<Real> >::iterator protein_map_iterator = protein_map.begin();
				protein_map_iterator != protein_map.end();
				protein_map_iterator++)
		{

			if (protein_map_iterator->second[2] > 0.9)
			{
				ProteinHit protein_hit;
				
				/// the protein score is the score + the average of the used thresholds
				protein_hit.setScore(protein_map_iterator->second[0] + 
														 (protein_map_iterator->second[1] /
														 protein_map_iterator->second[2]));
				
				protein_hit.setAccession(protein_map_iterator->first);
				protein_hit.setAccessionType("SwissProt");
				protein_hit.setScoreType("Mascot");
	
				db_searches_[0].insertProteinHit(protein_hit);
			}
		}

 		ok_ = true;
  }

  MascotOutfile::MascotOutfile(const MascotOutfile& source)
    : db_searches_(source.db_searches_), 
    	peptide_hits_(source.peptide_hits_), 
    	protein_hits_(source.protein_hits_), 
    	precursor_retention_times_(source.precursor_retention_times_), 
    	precursor_mz_values_(source.precursor_mz_values_), 
    	ok_(source.ok_)
  {
  }

  MascotOutfile::~MascotOutfile()
  {
  	peptide_hits_.clear();
  	protein_hits_.clear();
  }

  bool MascotOutfile::ok() const
  {
  	return ok_;
  }
  
  MascotOutfile& MascotOutfile::operator>>(Identification& db_search)
  {
  	
  	db_search.clear();
		db_search = Identification(db_searches_[0]);
    
    return *this;
  }
  

  MascotOutfile& MascotOutfile::operator>>(PeptideHit& peptide_hit)
  {
    peptide_hit.clear();
		
		if (curr_peptide_hit_ == peptide_hits_.end())
		{
      return *this;			
		}
		
		/// copy values from current hit
		peptide_hit = PeptideHit(*curr_peptide_hit_);
    ++curr_peptide_hit_;
    
    return *this;
  }

  MascotOutfile& MascotOutfile::operator>>(ProteinHit& protein_hit)
  {
    protein_hit.clear();
		
		if (curr_protein_hit_ == protein_hits_.end())
		{
      return *this;			
		}
		/// copy values from current hit
		protein_hit = ProteinHit(*curr_protein_hit_);
    ++curr_protein_hit_;
    
    return *this;
  }

	/// Assignment operator
  MascotOutfile& MascotOutfile::operator=(const MascotOutfile& source)
  {
  	if (this == &source)
  	{
  		return *this;
  	}
  	precursor_retention_times_ = source.precursor_retention_times_;  			
  	precursor_mz_values_ = source.precursor_mz_values_;  			
  	db_searches_ = source.db_searches_;
  	peptide_hits_ = source.peptide_hits_;
  	protein_hits_ = source.protein_hits_;

		return *this;
  }
  
  /// returns the retention time of the Mascot search
  const vector<float>& MascotOutfile::getPrecursorRetentionTimes() const
  {
  	return precursor_retention_times_;
  }

  /// sets the retention time of the Mascot search
  void MascotOutfile::setPrecursorRetentionTimes(const vector<float>& precursor_retention_times)
	{
		precursor_retention_times_ = precursor_retention_times;
	}
	
  /// returns the m/z of the precursor peak of the Mascot search
  const vector<float>& MascotOutfile::getPrecursorMZValues() const
  {
  	return precursor_mz_values_;
  }

  /// sets the m/z of the precursor peak of the Mascot search
  void MascotOutfile::setPrecursorMZValues(const vector<float>& precursor_mz_values)
  {
  	precursor_mz_values_ = precursor_mz_values;
  } 
  
  const vector<Identification>& MascotOutfile::getIdentifications() const
	{
		return db_searches_; 
	}
	
  void MascotOutfile::setIdentifications(const vector<Identification>& db_searches)
  {
  	db_searches_ = db_searches;
  }

} ///namespace OpenMS
