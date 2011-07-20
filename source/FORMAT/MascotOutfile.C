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
// $Maintainer: Nico Pfeifer $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/FORMAT/MascotOutfile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/DATASTRUCTURES/DateTime.h>

#include <cmath>

using namespace std;

namespace OpenMS 
{

  MascotOutfile::MascotOutfile()  
  {
  }

	void MascotOutfile::load(String filename,	ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& peptide_identifications, Real p)
	{
  	TextFile f(filename);
  	vector<PeptideHit>::iterator peptide_hit_iterator;
  	vector<PeptideHit> peptide_hits;
  	vector<ProteinHit> protein_hits;
		vector<String> parts;
		UInt number_of_queries = 0;
		map<UInt, UInt> indices;
		map<UInt, UInt>::iterator indices_iterator;
		PeptideIdentification temp_identification;
		vector<Int> charges;
		Int temp_charge = 0;
		String temp_identifier = "";
		vector<Real> temp_scores;
		map<String, vector<Real> > protein_map;
		Real temp_score = 0;
		Real temp_value = 0;
		Real temp_significance_threshold = 0;
		UInt index = 0;
		String::SizeType tag_start = 0;
		String::SizeType tag_end = 0;

		peptide_identifications.clear();
	
		if (f.size() == 0)
		{
			throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
				"File is empty!", filename);
		}
	
		//Mascot search was not successfull
		if (f.size()<5) return;
		  	
  	TextFile::iterator it;

  	// (1.0) parse for date

  	it = f.search("date=");
  	if (it==f.end())
  	{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
  			"date in header section not found!" ,filename);
  	}

		//PreciseTime precise_date(it->suffix('=').trim().toInt(),0);
		DateTime precise_date;
		precise_date.setTime_t(it->suffix('=').trim().toInt());
		
		DateTime date;
		it = f.search(it, "time=");
  	if (it==f.end())
  	{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
  			"time in header section not found!" ,filename);
  	}
		
		//date.set(ss.str().substr(6,2) + "." + ss.str().substr(4,2) + "." + 
		//	ss.str().substr(0,4) + " " + it->suffix('=').trim());
		date.setTime(it->suffix('=').trim());
		// now add the date
		// @todo fix this after cmake installation (Andreas, Chris)
		//date.setDate(precise_date);		
		
		//temp_identification.id.setDateTime(date);
		
		// TODO
		//temp_identification.setDateTime(date);
		
  	// (1.0.1) parse for number of queries
  	it = f.search(it, "queries=");
  	if (it==f.end())
  	{
  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
  			"number of queries not found!" ,filename);
  	}
		number_of_queries = it->suffix('=').trim().toInt();

		// (1.0.2) Searching for query indices for which peptides are present
		if (number_of_queries > 1)
		{
			for (UInt i = 1; i <= number_of_queries; i++)
			{
				//if no peptide is found for a certain query, there is just a -1 in the qi_p1 line
	  		it = f.search(it, "q" + String(i) + "_p1=");
	 		 	if (it!=f.end())
	  		{				  		
	  			if (it->suffix('=').at(0) != '-')
	  			{
	  				indices.insert(make_pair(i, (UInt)indices.size()));
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
	  
		// (1.1) parse for precursor values
		for(indices_iterator = indices.begin(); 
				indices_iterator != indices.end();
				++indices_iterator)
		{
			it = f.search("qexp" + String(indices_iterator->first) + "=");
			if (it == f.end())
			{
 				throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
 				"precursor information for query " + String(indices_iterator->first)
 				+ " not found!" ,filename);			
 			}
			it->suffix('=').split(',',parts);

			temp_charge = String(parts[1].trim()[0]).toInt();
			if (String(parts[1].trim()[1]) == "+")
			{
				charges.push_back(temp_charge);
			}
			else
			{
				charges.push_back((-1 * temp_charge));				
			}
			//temp_identification.mz = parts[0].toFloat();
			temp_identification.setMetaValue("MZ", parts[0].toFloat());
			peptide_identifications.push_back(temp_identification);
			parts.clear();
		}
		
		// (1.2) parse for peptide significance threshold
		for(indices_iterator = indices.begin(); 
				indices_iterator != indices.end();
				++indices_iterator)
		{
			it = f.search("qplughole" + String(indices_iterator->first) + "=");
			if (it == f.end())
			{
	  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
	  			"significance threshold for query " + String(indices_iterator->first)
 				+ " in summary section not found!" ,filename);			
			}
			peptide_identifications[indices_iterator->second].setSignificanceThreshold(
				it->suffix('=').trim().toFloat());
		}
		for(indices_iterator = indices.begin(); 
				indices_iterator != indices.end();
				++indices_iterator)
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
						< peptide_identifications[indices_iterator->second].getSignificanceThreshold())
			{
				peptide_identifications[indices_iterator->second].setSignificanceThreshold(temp_value);
			}
		}
		
		
		// The protein significance threshold is left zero because it is not
		// stored directly in the Mascot outfile and information about the 
		// calculation is, to our knowledge, not publicly available.				
		
		// (2.1) parse for ProteinHit information (MudPIT scoring)
		
		if (number_of_queries > 1000)
		{
	  	it = f.searchSuffix("\"proteins\"", true);
	  	
			if (it == f.end())
			{
				cout << "no \"proteins\" tag found " << "\n";
			}

	  	//Go to first protein hit entry
			if (it != f.end() && (it + 1) != f.end())
			{
				it += 2;
			}
	  	while(it != f.end())
	  	{
	  		// search for first "
	  		if ((tag_start = (*it).find('"')) != String::npos)
	  		{
	  			// search for second "
		  		if ((tag_end = (*it).find('"', tag_start + 1)) != String::npos)
		  		{
		  			temp_identifier = (*it).substr(tag_start + 1, tag_end - tag_start - 1);
						temp_scores.clear();
						for (Size k = 0; k < 3; k++)
						{
							temp_scores.push_back(0);
						}
						protein_map.insert(make_pair(temp_identifier, temp_scores));
						++it;
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

	  UInt i = 1; // first index of the peptidehits
	  UInt j = 1; // second index of the peptidehits
  	// (2.2) parse for PeptideHit information
		for(indices_iterator = indices.begin(); 
				indices_iterator != indices.end();
				++indices_iterator)
		{
			i = indices_iterator->first;
			j = 1;
			UInt counter = 1; //counter of the peptidehits
	  	it = f.search(String("q")+String(i)+"_p" + String(j) + "=");
			peptide_identifications[indices_iterator->second].setScoreType("Mascot");
	  	while(it != f.end())
	  	{
	  		PeptideHit hit;
				
	  		//(2.1) parse for peptide sequence
				it->suffix('=').split(',',parts);
	  		hit.setSequence(parts[4]);
				temp_score = parts[7].toFloat();			
				hit.setScore(temp_score);
				//hit.setScoreType("Mascot");
				hit.setCharge(charges[i - 1]);
	  		
	  		hit.setRank(counter);
	   		//(2.3) insert into hits vector
	  		if (temp_score > 0)
	  		{
	  			peptide_identifications[indices_iterator->second].insertHit(hit);
					counter++;
	  		}
	  		
	  		if (number_of_queries > 1000)
	  		{
	  			temp_significance_threshold = 
	  				peptide_identifications[indices_iterator->second].getSignificanceThreshold();
	
		  		if (temp_score > temp_significance_threshold)
		  		{
		  			
						it->suffix('=').split('"',parts);
			  		index = 1;
			  		while(index < parts.size() - 1)
			  		{
			  			temp_scores = protein_map[parts[index]];
			  					  			
			  			// the score of the protein hit
			  			temp_scores[0] += (temp_score - temp_significance_threshold);
			  			// sum of the used significance thresholds
			  			temp_scores[1] += temp_significance_threshold;
			  			// number of significance thresholds used in sum
			  			temp_scores[2] = temp_scores[2] + 1;
			  			
			  			protein_map[parts[index]] = temp_scores;
			  					  			
			  			index += 2;
			  		}
		  		}
	  		}
	  		
	  		//search for the next hit
	  		++j;
	  		it = f.search(it,String("q")+String(i)+"_p" + String(j) + "=");
				parts.clear();	
	  	}
		}

  	//(3) search for protein hit information
 		i = 1;
 		j = 1;
 		if (number_of_queries == 1)
 		{  	
	  	it = f.search(String("h")+String(i) + "=");
	  	peptide_hits = peptide_identifications[0].getHits();
	  	while(it != f.end())
	  	{
				ProteinHit protein_hit;
				String temp_peptide_sequence;
				SignedSize peptide_index = -1;
	
				protein_hit.setAccession(it->suffix('=').prefix(','));
				//protein_hit.setAccessionType("SwissProt");
				protein_hit.setScore(it->substr(it->find(',')+1).prefix(',').toFloat());
				//protein_hit.setScoreType("Mascot");
				protein_hit.setRank(i);
				
	  		it = f.search(it,String("h")+String(i)+"_q" + String(j) + "=");	
				if (it==f.end() && j == 1)
		  	{
		  		throw Exception::ParseError(__FILE__, __LINE__, __PRETTY_FUNCTION__, 
		  			String("Line starting with 'h")+String(i)+"_q1=' not found!" ,
		  			filename);
		  	}
		  	// search for the peptide hits belonging to the actual protein hit
		  	while(it != f.end())
		  	{ 	
					vector<String> parts;
					it->suffix('=').split(',',parts);
					temp_peptide_sequence = parts[6];
					
					for (Size index = 0; index < peptide_hits.size(); index++)
					{
						if (peptide_hits[index].getSequence() == temp_peptide_sequence)
						{
							peptide_index = index;
						}
					}
					// setting of the indices to store the relational information 
					if (peptide_index != -1)
					{
//						peptide_hits[peptide_index].addProteinIndex(i - 1);
//						protein_hit.addPeptideIndex(peptide_index);
						peptide_index = -1;
					}
									
					j++;
		  		it = f.search(it,String("h")+String(i)+"_q" + String(j) + "=");	
				}
				protein_hits.push_back(protein_hit);
				i++;
				j = 1;
		  	it = f.search(String("h")+String(i) + "=");
			}						  	
			//peptide_identifications[0].setPeptideAndProteinHits(peptide_hits, protein_hits);
			protein_identification.setHits(protein_hits);
		}
		
		
		UInt count = 0;
		for(indices_iterator = indices.begin(); 
				indices_iterator != indices.end();
				++indices_iterator)
		{
  		it = f.searchSuffix(String("\"query") + indices_iterator->first + "\"", true);	
			if (it!=f.end())
			{
				it = f.search(it, String("rtinseconds="));	
				if (it==f.end())
		  	{
		  		it = f.search(String("\"query") + indices_iterator->first + "\"", true);
				}
				else
				{
					peptide_identifications[count].setMetaValue("RT", it->suffix('=').trim().toFloat());
				}
			}				
			++count;
		}
	
		//protein_identification.setProteinAccessionType("SwissProt");
		protein_identification.setScoreType("Mascot");
		for(map<String, vector<Real> >::iterator protein_map_iterator = protein_map.begin();
				protein_map_iterator != protein_map.end();
				++protein_map_iterator)
		{

			if (protein_map_iterator->second[2] > 0.9)
			{
				ProteinHit protein_hit;
				
				// the protein score is the score + the average of the used thresholds
				protein_hit.setScore(protein_map_iterator->second[0] + 
														 (protein_map_iterator->second[1] /
														 protein_map_iterator->second[2]));
				
				protein_hit.setAccession(protein_map_iterator->first);
	
				protein_identification.insertHit(protein_hit);
			}
		}
	}

} //namespace OpenMS
