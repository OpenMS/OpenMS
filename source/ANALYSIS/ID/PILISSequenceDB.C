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
		: proteins_(db.proteins_),
			peptides_(db.peptides_),
			factor_(db.factor_),
			replace_X_and_L_(db.replace_X_and_L_)
	{
	}

	PILISSequenceDB& PILISSequenceDB::operator = (const PILISSequenceDB& rhs)
	{
		if (this != &rhs)
		{
			proteins_ = rhs.proteins_;
			peptides_ = rhs.peptides_;
			factor_ = rhs.factor_;
			replace_X_and_L_ = rhs.replace_X_and_L_;
		}
		return *this;
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
			if (split.size() == 3)
			{
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

			if (split.size() == 0)
			{
				// contains only peptide
				String peptide(s);
				PepStruct new_peptide;

				// add peptide with charge 1
				new_peptide.peptide = peptide;
				new_peptide.weight = AASequence(peptide).getAverageWeight(Residue::Full, 1);
				new_peptide.charge = 1;
				peptides_[(Size)(new_peptide.weight * factor_)].push_back(new_peptide);

				// add peptide with charge 2
				new_peptide.weight = (new_peptide.weight + 1.0)/2.0;
				new_peptide.charge = 2;
				peptides_[(Size)(new_peptide.weight * factor_)].push_back(new_peptide);
			}
  	}
		return;
	}

	void PILISSequenceDB::addFASTAFile(const String& /*filename*/)
	{
		throw Exception::NotImplemented(__FILE__, __LINE__, __PRETTY_FUNCTION__);
	}

	void PILISSequenceDB::digestProteinsTryptic(Size missed_cleavages)
	{
		for (vector<pair<String, String> >::const_iterator it = proteins_.begin(); it != proteins_.end(); ++it)
		{
			vector<String> peptides;
			digestTryptic_(it->second, peptides, missed_cleavages);

			for (vector<String>::const_iterator it1 = peptides.begin(); it1 != peptides.end(); ++it1)
			{
				addPeptide_(*it1);
			}
		}
	}

	void PILISSequenceDB::addPeptide_(const String& peptide)
	{
		if (!has(peptide))
		{
			AASequence peptide_sequence(peptide);
			double weight = peptide_sequence.getAverageWeight(Residue::Full, 2)/2.0;
			PepStruct p;
			p.peptide = peptide;
			p.weight = weight;
			p.charge = 2;
			peptides_[(Size)(weight * factor_)].push_back(p);

			weight = peptide_sequence.getAverageWeight(Residue::Full, 1);
			p.weight = weight;
			p.charge = 1;
			peptides_[(Size)(weight * factor_)].push_back(p);
		}
	}
	
	bool PILISSequenceDB::has(const String& peptide) const
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
			double weight = peptide_sequence.getAverageWeight(Residue::Full, 2)/2.0;
			if (!peptides_.has((Size)(weight * factor_)))
			{
				return false;
			}

			for (vector<PepStruct>::const_iterator it = peptides_[(Size)(weight * factor_)].begin(); it != peptides_[(Size)(weight * factor_)].end(); ++it)
			{
				if (it->peptide == peptide)
				{
					return true;
				}
			}

			weight = peptide_sequence.getAverageWeight(Residue::Full, 1);
			if (!peptides_.has((Size)(weight * factor_)))
			{
				return false;
			}
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

	unsigned int PILISSequenceDB::countProteins() const
	{
		return proteins_.size();
	}

	void PILISSequenceDB::digestTryptic_(const String& seq, vector<String>& peptides, Size missedcleavages)
	{
		// trypsin cuts C-terminal to K,R, and not before proline
    for (Size i=0; i != seq.size(); ++i)
    {
     	if (seq[i] == 'K' || seq[i] == 'R')
     	{
       	Size j = i + 1;
       	for (Size k = 0; k <= missedcleavages; ++k)
       	{
         	while (j < seq.size() && !(seq[j] == 'K' || seq[j] == 'R'))
         	{
           	++j;
         	}
         	if (seq[j] == 'K' || seq[j] == 'R')
         	{
           	if (j - i > 2)
           	{
             	if (j + 1 != seq.size() && seq[j + 1] != 'P')
             	{
               	peptides.push_back(seq.substr(i + 1, j - i));
             	}
             	else
             	{
               	if (j + 1 == seq.size())
               	{
                 	peptides.push_back(seq.substr(i + 1, j - i));
               	}
             	}
           	}
         	}
         	++j;
       	}
     	}
    }
		return;
	}
}

