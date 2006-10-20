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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CHEMISTRY/EnzymaticDigestion.h>

#include <iostream>

using namespace std;

namespace OpenMS
{
	const std::string EnzymaticDigestion::NamesOfEnzymes[] = {"Trypsin"};
	
	EnzymaticDigestion::EnzymaticDigestion()
		: missed_cleavages_(0),
			enzyme_(TRYPSIN)
	{
		
	}
	
	UnsignedInt EnzymaticDigestion::getMissedCleavages() const
	{
		return missed_cleavages_;
	}

	void EnzymaticDigestion::setMissedCleavages(UnsignedInt missed_cleavages)
	{
		missed_cleavages_ = missed_cleavages;
	}

	EnzymaticDigestion::Enzyme EnzymaticDigestion::getEnzyme() const
	{
		return enzyme_;
	}
	
	void EnzymaticDigestion::setEnzyme(Enzyme enzyme)
	{
		if (enzyme < SIZE_OF_ENZYMES) enzyme_ = enzyme;
	}

	void EnzymaticDigestion::nextCleavageSite_(const PeptideSequence& protein, PeptideSequence::ConstIterator& iterator)
	{
		//cout << "nextCleavageSite_" << endl;
		switch (enzyme_)
		{
			case TRYPSIN:	 
				while (iterator != protein.end())
				{
					//R or K at the end or not P afterwards
					if ((**iterator=='R' || **iterator=='K') && ((iterator+1)==protein.end() || **(iterator+1)!='P'))
					{
						++iterator;
						return;
					}
					//cout << "it: "<< **iterator << endl;
					++iterator;
				}
				break;
			default:
				return;
		};
	}
	
	UnsignedInt EnzymaticDigestion::peptideCount(const PeptideSequence& protein)
	{
		UnsignedInt count = 1;
		PeptideSequence::ConstIterator iterator = protein.begin();
		while(nextCleavageSite_(protein,iterator), iterator != protein.end())
		{
			++count;
		}
		
		//missed cleavages
		UnsignedInt sum = count;
		for (UnsignedInt i=1 ; ((i<=missed_cleavages_) && (count > i)); ++i)
		{
			sum += count - i;
		}
		
		return sum;
	}

	void EnzymaticDigestion::digest(const PeptideSequence& protein, std::vector<PeptideSequence>& output)
	{
		//initialization
		UnsignedInt count = 1;
		output.clear();
		
		//missed cleavage iterators
		vector<PeptideSequence::ConstIterator> mc_iterators;
		if (missed_cleavages_ != 0) mc_iterators.push_back(protein.begin());
		
		PeptideSequence::ConstIterator begin = protein.begin();
		PeptideSequence::ConstIterator end = protein.begin();
		while(nextCleavageSite_(protein,end), end != protein.end())
		{
			++count;
			if (missed_cleavages_ != 0) mc_iterators.push_back(end);
			output.push_back(PeptideSequence(begin,end));
			begin = end;
		}
		output.push_back(PeptideSequence(begin,end));
		if (missed_cleavages_ != 0) mc_iterators.push_back(end);
		
		//missed cleavages
		if (mc_iterators.size()>2) //there is at least one cleavage site!
		{
			//resize to number of fragments
			UnsignedInt sum = count;
			for (UnsignedInt i=1 ; ((i<=missed_cleavages_) && (count > i)); ++i)
			{
				sum += count - i;
			}

			output.resize(sum);
			
			//generate fragments with missed cleavages
			UnsignedInt pos = count;
			for (UnsignedInt i=1 ; ((i<=missed_cleavages_) && (count > i)) ; ++i)
			{
				vector<PeptideSequence::ConstIterator>::const_iterator b = mc_iterators.begin();
				vector<PeptideSequence::ConstIterator>::const_iterator e = b+(i+1);
				while (e != mc_iterators.end())
				{
					output[pos] = PeptideSequence(*b,*e);
					++b;
					++e;
					++pos;
				}
			}
		}
	}

}	//namespace
