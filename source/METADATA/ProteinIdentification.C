// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/METADATA/PeptideHit.h>

#include <sstream>
#include <algorithm>

using namespace std;

namespace OpenMS
{

	const std::string ProteinIdentification::NamesOfPeakMassType[] = {"Monoisotopic","Average"};
	const std::string ProteinIdentification::NamesOfDigestionEnzyme[] = {"Trypsin", "Pepsin A", "Protease K", "Chymotrypsin", "No enzyme", "Unknown"};

  ProteinIdentification::ProteinIdentification()
    : MetaInfoInterface(),
			id_(),
			search_engine_(),
			search_engine_version_(),
			search_parameters_(),
			date_(),
			protein_score_type_(),
			higher_score_better_(true),
			protein_hits_(),
			protein_groups_(),
			protein_significance_threshold_(0.0)
  {
  }

  ProteinIdentification::ProteinIdentification(const ProteinIdentification& source)
		: MetaInfoInterface(source),
			id_(source.id_),
			search_engine_(source.search_engine_),
			search_engine_version_(source.search_engine_version_),
			search_parameters_(source.search_parameters_),
			date_(source.date_),
			protein_score_type_(source.protein_score_type_),
			higher_score_better_(source.higher_score_better_),
			protein_hits_(source.protein_hits_),
			protein_groups_(source.protein_groups_),
	  	protein_significance_threshold_(source.protein_significance_threshold_)
  {
  }

  ProteinIdentification::~ProteinIdentification()
  {
  }

	void ProteinIdentification::setDateTime(const DateTime& date)
	{
		date_ = date;
	}

	const DateTime& ProteinIdentification::getDateTime() const
	{
		return date_;
	}

	const vector<ProteinHit>& ProteinIdentification::getHits() const
	{
		return protein_hits_;
	}

	vector<ProteinHit>& ProteinIdentification::getHits()
	{
		return protein_hits_;
	}

	void ProteinIdentification::setHits(const vector<ProteinHit>& protein_hits)
	{
		// groups might become invalid by this operation
		if (protein_groups_.size()>0)
		{
			LOG_ERROR << "New protein hits set while proteins groups are non-empty! This might invalidate groups! Delete groups before setting new hits!\n";
		}
		protein_hits_ = protein_hits;
	}
	
	const vector<ProteinIdentification::ProteinGroup>& ProteinIdentification::getProteinGroups() const
	{
		return protein_groups_;
	}
	vector<ProteinIdentification::ProteinGroup>& ProteinIdentification::getProteinGroups()
	{
		return protein_groups_;
	}
	void ProteinIdentification::insertGroup(const ProteinIdentification::ProteinGroup& group)
	{
		protein_groups_.push_back(group);
	}


	// retrival of the peptide significance threshold value
  DoubleReal ProteinIdentification::getSignificanceThreshold() const
  {
  	return protein_significance_threshold_;
  }

	// setting of the peptide significance threshold value
	void ProteinIdentification::setSignificanceThreshold(DoubleReal value)
	{
		protein_significance_threshold_ = value;
	}

	void ProteinIdentification::setScoreType(const String& type)
	{
		protein_score_type_ = type;
	}

	const String& ProteinIdentification::getScoreType() const
	{
		return protein_score_type_;
	}

	void ProteinIdentification::insertHit(const ProteinHit& protein_hit)
	{
		protein_hits_.push_back(protein_hit);
	}

  ProteinIdentification& ProteinIdentification::operator=(const ProteinIdentification& source)
  {
  	if (this == &source)
  	{
  		return *this;
  	}
    MetaInfoInterface::operator=(source);
		id_ = source.id_;
		search_engine_ = source.search_engine_;
		search_engine_version_ = source.search_engine_version_;
		search_parameters_ = source.search_parameters_;
		date_ = source.date_;
    protein_hits_ = source.protein_hits_;
    protein_groups_ = source.protein_groups_;
		protein_score_type_ = source.protein_score_type_;
    protein_significance_threshold_ = source.protein_significance_threshold_;
		higher_score_better_ = source.higher_score_better_;
    return *this;
  }

	// Equality operator
	bool ProteinIdentification::operator == (const ProteinIdentification& rhs) const
	{
		return 	MetaInfoInterface::operator==(rhs) &&
						id_ == rhs.id_ &&
						search_engine_ == rhs.search_engine_ &&
				    search_engine_version_ == rhs.search_engine_version_ &&
						search_parameters_ == rhs.search_parameters_ &&
						date_ == rhs.date_ &&
						protein_hits_ == rhs.protein_hits_ &&
						protein_groups_ == rhs.protein_groups_ &&
						protein_score_type_ == rhs.protein_score_type_ &&
						protein_significance_threshold_ == rhs.protein_significance_threshold_ &&
						higher_score_better_ == rhs.higher_score_better_;

	}

	// Inequality operator
	bool ProteinIdentification::operator != (const ProteinIdentification& rhs) const
	{
		return !operator==(rhs);
	}


  void ProteinIdentification::sort()
  {
		// while sorting proteins, we need to update the groups as well
		
		// create index vector ..
		vector<size_t> index;
		for (size_t i = 0; i < protein_hits_.size(); ++i)	index.push_back(i);

		// .. and sort it
 		if (higher_score_better_)
  	{
			std::sort(index.begin(), index.end(), ProteinCmp<std::vector<ProteinHit>&, ProteinHit::ScoreMore>(protein_hits_, ProteinHit::ScoreMore()) );
			// here we sort the actual proteins (could be done using sorted indices as well, but the list is short)
			std::sort(protein_hits_.begin(), protein_hits_.end(), ProteinHit::ScoreMore());
  	}
  	else
  	{
  		std::sort(index.begin(), index.end(), ProteinCmp<std::vector<ProteinHit>&, ProteinHit::ScoreLess>(protein_hits_, ProteinHit::ScoreLess()) );
  		// here we sort the actual proteins (could be done using sorted indices as well, but the list is short)
  		std::sort(protein_hits_.begin(), protein_hits_.end(), ProteinHit::ScoreLess());
  	}
  	
  	// create mapping (from old index to new one)
  	Map<size_t, size_t> index_to_index;
  	for (size_t i = 0; i < protein_hits_.size(); ++i)	index_to_index[index[i]] = i;
  	
  	// update groups
  	for (std::vector<ProteinGroup>::iterator it=protein_groups_.begin(); it!=protein_groups_.end();++it)
  	{
  		std::set<Size> indices_new;
  		for (std::set<Size>::const_iterator it_set=it->indices.begin();it_set!=it->indices.end();++it_set)
  		{
  			indices_new.insert(index_to_index[*it_set]);
  		}
  		it->indices=indices_new;
  	}
  	
  }
  
  
  void ProteinIdentification::assignRanks()
  {
		UInt rank = 1;
    sort();
    vector<ProteinHit>::iterator lit = protein_hits_.begin();
    Real tmpscore = lit->getScore();
    while (  lit != protein_hits_.end() )
    {
     	lit->setRank(rank);
     	++lit;
     	if ( lit != protein_hits_.end() && lit->getScore() != tmpscore )
     	{
       	++rank;
       	tmpscore = lit->getScore();
     	}
		}
  }


	bool ProteinIdentification::isHigherScoreBetter() const
	{
		return higher_score_better_;
	}

	void ProteinIdentification::setHigherScoreBetter(bool value)
	{
		higher_score_better_ = value;
	}

	const String& ProteinIdentification::getIdentifier() const
	{
		return id_;
	}

	void ProteinIdentification::setIdentifier(const String& id)
	{
		id_ = id;
	}

	void ProteinIdentification::setSearchEngine(const String& search_engine)
	{
		search_engine_ = search_engine;
	}

	const String& ProteinIdentification::getSearchEngine() const
	{
		return search_engine_;
	}

	void ProteinIdentification::setSearchEngineVersion(const String& search_engine_version)
	{
		search_engine_version_ = search_engine_version;
	}

	const String& ProteinIdentification::getSearchEngineVersion() const
	{
		return search_engine_version_;
	}

	void ProteinIdentification::setSearchParameters(const SearchParameters& search_parameters)
	{
		search_parameters_ = search_parameters;
	}

	const ProteinIdentification::SearchParameters& ProteinIdentification::getSearchParameters() const
	{
		return search_parameters_;
	}

}// namespace OpenMS
