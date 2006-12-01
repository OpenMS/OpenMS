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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------


#include <OpenMS/ANALYSIS/ID/PILISSequenceDB.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <fstream>
#include <iostream>

using namespace std;

namespace OpenMS
{
	
	PILISSequenceDB::PILISSequenceDB()
		:	factor_(10.0),
			replace_X_and_L_(true)
	{
	}

	PILISSequenceDB::PILISSequenceDB(const PILISSequenceDB& db)
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}
	
	PILISSequenceDB::~PILISSequenceDB()
	{
	}

	void PILISSequenceDB::addPeptidesFromFile(const String& filename)
	{
		ifstream seq_db_in(filename.c_str());
		char line[1000];
		while (seq_db_in.getline(line, 1000))
  	{
    	String s = String(line).trim();
    	vector<String> split;
    	s.split(' ', split);
			String peptide = split[0];
			if (replace_X_and_L_)
			{
				replace(peptide.begin(), peptide.end(), 'X', 'I');
				replace(peptide.begin(), peptide.end(), 'L', 'I');
			}
    	PepStruct new_peptide;
    	new_peptide.peptide = peptide;
    	new_peptide.weight = split[1].toFloat();
    	new_peptide.charge = split[2].toInt();
    	peptides_[(Size)(new_peptide.weight * factor_)].push_back(new_peptide);
  	}
		return;
	}

	void PILISSequenceDB::addFASTAFile(const String& filename)
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}

	void PILISSequenceDB::digestSequencesTryptic(Size missed_cleavages)
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}

	bool PILISSequenceDB::has(const String& peptide)
	{
		// TODO replace Z by K,Q
		for (Size i = 0; i != peptide.size(); ++i)
		{
			if (peptide[i] == 'Z')
			{
				String p1(peptide), p2(peptide);
				p1[i] = 'K';
				p2[i] = 'Q';
				return has(p1) || has(p2); 
			}
		}
		
		try 
		{
			AASequence peptide_sequence(peptide);
			double weight = peptide_sequence.getAverageWeight(Residue::Full, 2);

			for (vector<PepStruct>::const_iterator it = peptides_[(Size)(weight * factor_)].begin(); it != peptides_[(Size)(weight * factor_)].end(); ++it)
			{
				if (it->peptide == peptide)
				{
					return true;
				}
			}

			weight = peptide_sequence.getAverageWeight(Residue::Full, 1);
			for (vector<PepStruct>::const_iterator it = peptides_[(Size)(weight * factor_)].begin(); it != peptides_[(Size)(weight * factor_)].end(); ++it)
 	   	{
 	    	if (it->peptide == peptide)
 	     	{
 	       	return true;
 	     	}
 	   	}
		}
		catch (Exception::ParseError e)
		{
			return false;
		}

		return false;
	}

	void PILISSequenceDB::getPeptides(std::vector<PepStruct>& peptides, double range_start, double range_stop)
	{
		if (range_start > range_stop)
		{
			cerr << "PILISSequenceDB: nonsense ranges, start=" << range_start << ", stop=" << range_stop << endl;
			return;
		}

		for (Size i = (Size)(range_start * factor_); i <= (Size)(range_stop * factor_); ++i)
		{
			for (vector<PepStruct>::const_iterator it = peptides_[i].begin(); it != peptides_[i].end(); ++it)
			{
				peptides.push_back(*it);
			}
		}

		return;
	}

	void PILISSequenceDB::clearProteins()
	{
		proteins_.clear();
	}

	void PILISSequenceDB::clearPeptides()
	{
		peptides_.clear();
	}

	void PILISSequenceDB::setFactor(double factor)
	{
		factor_ = factor;
	}

	double PILISSequenceDB::getFactor() const
	{
		return factor_;
	}

	void PILISSequenceDB::setReplaceXandL(bool replace)
	{
		replace_X_and_L_ = replace;
	}

	bool PILISSequenceDB::isReplaceXandL() const
	{
		return replace_X_and_L_;
	}

	unsigned int PILISSequenceDB::countPeptides() const
	{
		unsigned int count(0);
		for (HashMap<Size, std::vector<PepStruct> >::ConstIterator it = peptides_.begin(); it != peptides_.end(); ++it)
		{
			count += it->second.size();
		}
		return count;
	}
}

