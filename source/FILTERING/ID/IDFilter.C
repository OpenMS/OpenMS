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
// $Id: IDFilter.C,v 1.15 2006/06/09 23:47:35 nicopfeifer Exp $
// $Author: nicopfeifer $
// $Maintainer: Nico Pfeifer $
// --------------------------------------------------------------------------

#include <OpenMS/FILTERING/ID/IDFilter.h>

#include <cmath>

using namespace std;

namespace OpenMS 
{
  IDFilter::IDFilter(): 
    peptide_threshold_fraction_(1),
    protein_threshold_fraction_(1),
    proteins_()
  {
    
  }
  
  IDFilter::IDFilter(const IDFilter& source):
    peptide_threshold_fraction_(source.peptide_threshold_fraction_),
    protein_threshold_fraction_(source.protein_threshold_fraction_),
    proteins_(source.proteins_)
  {
    
  }
   
  IDFilter::~IDFilter()
  {
    
  }
  
  IDFilter& IDFilter::operator = (const IDFilter& source)
  {
    if (this == &source)
    { 
    	return *this;
    }
    
    peptide_threshold_fraction_ = source.peptide_threshold_fraction_;
    protein_threshold_fraction_ = source.protein_threshold_fraction_;
    proteins_ = source.proteins_;
    
    return *this;
  }

  const double& IDFilter::getPeptideThresholdFraction() const
  {
    return peptide_threshold_fraction_;
  }
 
  void IDFilter::setPeptideThresholdFraction(const double& peptide_threshold_fraction)
  {
    peptide_threshold_fraction_ = peptide_threshold_fraction;
  }
  
  const double& IDFilter::getProteinThresholdFraction() const
  {
    return protein_threshold_fraction_;
  }
 
  void IDFilter::setProteinThresholdFraction(const double& protein_threshold_fraction)
  {
    protein_threshold_fraction_ = protein_threshold_fraction;
  }
 
  const vector< pair<String, String> >& IDFilter::getProteins() const
  {
    return proteins_;
  }
 
  void IDFilter::setProteins(const vector< pair<String, String> >& proteins)
  {
    proteins_ = proteins;
  }
	
	const Identification& IDFilter::filterIdentificationsByThresholds(const Identification& identification,
																															const double& peptide_threshold_fraction,
																															const double& protein_threshold_fraction,
																															bool strict)
	{
		peptide_threshold_fraction_ = peptide_threshold_fraction;
		protein_threshold_fraction_ = protein_threshold_fraction;
		
		return filterIdentificationsByThresholds(identification, strict);
	}

	const Identification& IDFilter::filterIdentificationsByThresholds(const Identification& identification, bool strict)
	{
		Identification* filtered_identification = 0;
		vector<PeptideHit> temp_peptide_hits;
		vector<PeptideHit> filtered_peptide_hits;
		PeptideHit temp_peptide_hit;
		vector<ProteinHit> temp_protein_hits;
		vector<ProteinHit> filtered_protein_hits;
		ProteinHit temp_protein_hit;
		float max_value = 0;
		float temp_score = 0;
		vector< UnsignedInt > new_peptide_indices;		
		vector< UnsignedInt > new_protein_indices;		
		DateTime date;
		
		date = identification.getDateTime();

		temp_protein_hits = identification.getProteinHits();
		temp_peptide_hits = identification.getPeptideHits();
		
		if (temp_peptide_hits.size() > 0)
		{
			max_value = temp_peptide_hits[0].getScore();
			new_peptide_indices.push_back(0);
			
			/// searching for peptide(s) with maximal score			
			for(UnsignedInt i = 1; i < temp_peptide_hits.size(); i++)
			{
				temp_score = temp_peptide_hits[i].getScore();
				if (temp_score > max_value)
				{
					max_value = temp_score;
					new_peptide_indices.clear();
					new_peptide_indices.push_back(i);
				}				
				else if (temp_score == max_value)
				{
					new_peptide_indices.push_back(i);
				}
			}
						
			/// testing whether best peptide fullfills filter constraint
			if (max_value >= 
					peptide_threshold_fraction_ * 
					identification.getPeptideSignificanceThreshold())
			{	
				if (!strict || new_peptide_indices.size() == 1)
				{
					for(UnsignedInt i = 0; i < new_peptide_indices.size(); i++)
					{
						temp_peptide_hit = PeptideHit(temp_peptide_hits[new_peptide_indices[i]]);
						temp_peptide_hit.setRank((i + 1));
  					filtered_peptide_hits.push_back(temp_peptide_hit);
  				}
  			}
			}	
			else
			{
				new_peptide_indices.clear();	
			}	
		}

		if (temp_protein_hits.size() > 0)
		{
			max_value = temp_protein_hits[0].getScore();
			new_protein_indices.push_back(0);

			/// searching for protein(s) with maximal score			
			for(UnsignedInt i = 1; i < temp_protein_hits.size(); i++)
			{
				temp_score = temp_protein_hits[i].getScore();
				if (temp_score > max_value)
				{
					max_value = temp_score;
					new_protein_indices.clear();
					new_protein_indices.push_back(i);
				}				
				else if (temp_score == max_value)
				{
					new_protein_indices.push_back(i);
				}
			}
						
			/// testing whether best protein(s) fullfill(s) filter constraint
			if (max_value >= 
					protein_threshold_fraction_ * 
					identification.getProteinSignificanceThreshold())
			{		
				if (!strict || new_protein_indices.size() == 1)
				{
					for(UnsignedInt i = 0; i < new_protein_indices.size(); i++)
					{
						temp_protein_hit = ProteinHit(temp_protein_hits[new_protein_indices[i]]);
						temp_protein_hit.setRank((i + 1));
 	 					filtered_protein_hits.push_back(temp_protein_hit);
 	 				}
 	 			}
			}
			else
			{
				new_protein_indices.clear();	
			}	
		}
		
		filtered_identification = new Identification();
		if (filtered_peptide_hits.size() > 0 || filtered_protein_hits.size() > 0)
		{
  		filtered_identification->setPeptideAndProteinHits(filtered_peptide_hits, 
  																											filtered_protein_hits);
			filtered_identification->setPeptideSignificanceThreshold(identification.getPeptideSignificanceThreshold());
			filtered_identification->setProteinSignificanceThreshold(identification.getProteinSignificanceThreshold());  																								
			filtered_identification->setCharge(identification.getCharge());  																								
			filtered_identification->setDateTime(date);  																								
		}
		return *filtered_identification;		 
	}

	const Identification& IDFilter::filterIdentificationsByProteins(const Identification& identification, 
																														vector< pair<String, String> >proteins)
	{
		proteins_ = proteins;
		
		return filterIdentificationsByProteins(identification);
	}

	const Identification& IDFilter::filterIdentificationsByProteins(const Identification& identification)
	{
		String protein_sequences;
		String accession_sequences;
		vector< UnsignedInt > new_peptide_indices;		
		vector< UnsignedInt > new_protein_indices;
		vector<PeptideHit> temp_peptide_hits;
		vector<PeptideHit> filtered_peptide_hits;
		vector<ProteinHit> temp_protein_hits;
		vector<ProteinHit> filtered_protein_hits;
		ProteinHit temp_protein_hit;
		PeptideHit temp_peptide_hit;
		Identification* filtered_identification = 0;		
		DateTime date;		
		
		date = identification.getDateTime();
		temp_peptide_hits = identification.getPeptideHits();
		temp_protein_hits = identification.getProteinHits();
		for(UnsignedInt i = 0; i < proteins_.size(); i++)
		{
			accession_sequences.append("*" + proteins_[i].first);
			protein_sequences.append("*" + proteins_[i].second);
		}
		accession_sequences.append("*");
		
		for(UnsignedInt i = 0; i < temp_peptide_hits.size(); i++)
		{
	  	if (protein_sequences.find(temp_peptide_hits[i].getSequence().c_str()) 
	  			!= string::npos)
	  	{
	  		new_peptide_indices.push_back(i);
	  	}
		}
		for(UnsignedInt i = 0; i < temp_protein_hits.size(); i++)
		{
	  	if (accession_sequences.find("*" 
	  			+ temp_protein_hits[i].getAccession() 
					+ "*") != string::npos)
	  	{
	  		new_protein_indices.push_back(i);
	  	}
		}
		for(UnsignedInt i = 0; i < new_peptide_indices.size(); i++)
		{
			temp_peptide_hit = PeptideHit(temp_peptide_hits[new_peptide_indices[i]]);
			temp_peptide_hit.setRank((i + 1));
			filtered_peptide_hits.push_back(temp_peptide_hit);
		}
		for(UnsignedInt i = 0; i < new_protein_indices.size(); i++)
		{
			temp_protein_hit = ProteinHit(temp_protein_hits[new_protein_indices[i]]);
			temp_protein_hit.setRank((i + 1));
			filtered_protein_hits.push_back(temp_protein_hit);
		}
		filtered_identification = new Identification();
		if (filtered_peptide_hits.size() > 0 || filtered_protein_hits.size() > 0)
		{
			/*
			for(UnsignedInt i = 0; i < filtered_peptide_hits.size(); i++)
			{
				filtered_peptide_hits[i].updateProteinIndices(new_protein_indices);
			}
			for(UnsignedInt i = 0; i < filtered_protein_hits.size(); i++)
			{
				filtered_protein_hits[i].updatePeptideIndices(new_peptide_indices);
			}
			*/					
  		filtered_identification->setPeptideAndProteinHits(filtered_peptide_hits, 
  																								filtered_protein_hits);
			filtered_identification->setPeptideSignificanceThreshold(identification.getPeptideSignificanceThreshold());
			filtered_identification->setProteinSignificanceThreshold(identification.getProteinSignificanceThreshold());  																								
			filtered_identification->setCharge(identification.getCharge());  																								
			filtered_identification->setDateTime(date);  																								
		}
		return *filtered_identification;		 				
	}
	
	const Identification& IDFilter::filterIdentificationsByRetentionTimes(const Identification& identification,
																																	const map<String, double>& predicted_retention_times,
																																	double measured_retention_time,
																																	double predicted_sigma,
																																	double allowed_deviation_factor,
																																	double total_gradient_time)
	{
		Identification* filtered_identification = 0;
		vector<PeptideHit> temp_peptide_hits;
		vector<PeptideHit>::iterator it;
		PeptideHit temp_peptide_hit;
		String temp_sequence;
		double temp_retention_time = 0;
		map<String, double>::const_iterator const_it;
		double difference = 0;
		double allowed_deviation = 0;
		DateTime date;		
		
		allowed_deviation = allowed_deviation_factor * sqrt(2.0) * predicted_sigma;
		date = identification.getDateTime();

		filtered_identification = new Identification();
		
		temp_peptide_hits = identification.getPeptideHits();
		for(it = temp_peptide_hits.begin();
				it != temp_peptide_hits.end();
				it++)
		{
			temp_sequence = it->getSequence();
			if ((const_it = predicted_retention_times.find(temp_sequence)) 
					!= predicted_retention_times.end())
			{
				temp_retention_time = const_it->second;
				difference = measured_retention_time - temp_retention_time;
				difference /= total_gradient_time;
				if (difference < 0)
				{
					difference *= -1;
				}
				if (difference <= allowed_deviation)
				{
					filtered_identification->insertPeptideHit(*it);
				}
			}
		}
		if (filtered_identification->getPeptideHits().size() > 0)
		{
			filtered_identification->setPeptideSignificanceThreshold(identification.getPeptideSignificanceThreshold());
			filtered_identification->setProteinSignificanceThreshold(identification.getProteinSignificanceThreshold());  																								
			filtered_identification->setCharge(identification.getCharge());  																								
			filtered_identification->setDateTime(date);  																								
		}
		
		return *filtered_identification;
	}
	
	const Identification& IDFilter::filterIdentificationsByExclusionPeptides(const Identification& identification,
																																		 vector<String> peptides)
	{
		String protein_sequences;
		String accession_sequences;
		vector< UnsignedInt > new_peptide_indices;		
		vector< UnsignedInt > new_protein_indices;
		vector<PeptideHit> temp_peptide_hits;
		vector<PeptideHit> filtered_peptide_hits;
		vector<ProteinHit> temp_protein_hits;
		vector<ProteinHit> filtered_protein_hits;
		ProteinHit temp_protein_hit;
		PeptideHit temp_peptide_hit;
		Identification* filtered_identification = 0;
		DateTime date;		
		
		temp_peptide_hits = identification.getPeptideHits();
		date = identification.getDateTime();
		
		for(UnsignedInt i = 0; i < temp_peptide_hits.size(); i++)
		{
	  	if (find(peptides.begin(), peptides.end(), 
	  			temp_peptide_hits[i].getSequence().c_str()) 
	  			== peptides.end())
	  	{
	  		new_peptide_indices.push_back(i);
	  	}
		}
		for(UnsignedInt i = 0; i < temp_protein_hits.size(); i++)
		{
	  		new_protein_indices.push_back(i);
		}
		for(UnsignedInt i = 0; i < new_peptide_indices.size(); i++)
		{
			temp_peptide_hit = PeptideHit(temp_peptide_hits[new_peptide_indices[i]]);
			temp_peptide_hit.setRank((i + 1));
			filtered_peptide_hits.push_back(temp_peptide_hit);
		}
		for(UnsignedInt i = 0; i < new_protein_indices.size(); i++)
		{
			temp_protein_hit = ProteinHit(temp_protein_hits[new_protein_indices[i]]);
			temp_protein_hit.setRank((i + 1));
			filtered_protein_hits.push_back(temp_protein_hit);
		}
		filtered_identification = new Identification();
		if (filtered_peptide_hits.size() > 0 || filtered_protein_hits.size() > 0)
		{
			/*
			for(UnsignedInt i = 0; i < filtered_peptide_hits.size(); i++)
			{
				filtered_peptide_hits[i].updateProteinIndices(new_protein_indices);
			}
			for(UnsignedInt i = 0; i < filtered_protein_hits.size(); i++)
			{
				filtered_protein_hits[i].updatePeptideIndices(new_peptide_indices);
			}
			*/					
  		filtered_identification->setPeptideAndProteinHits(filtered_peptide_hits, 
  																								filtered_protein_hits);
			filtered_identification->setPeptideSignificanceThreshold(identification.getPeptideSignificanceThreshold());
			filtered_identification->setProteinSignificanceThreshold(identification.getProteinSignificanceThreshold());  																								
			filtered_identification->setCharge(identification.getCharge());
			filtered_identification->setDateTime(date);  																								
		}
		return *filtered_identification;		 				
	}
		
} // namespace OpenMS
