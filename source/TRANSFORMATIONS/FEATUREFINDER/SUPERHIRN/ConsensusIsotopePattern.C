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
// $Maintainer: Florian Zeller $
// $Authors: Lukas Mueller, Markus Mueller $
// --------------------------------------------------------------------------
//
///////////////////////////////////////////////////////////////////////////
//
//  written by Lukas N Mueller, 30.3.05
//  Lukas.Mueller@imsb.biol.ethz.ch
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
// 


#include <map>
#include <vector>
#include <math.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ConsensusIsotopePattern.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnUtil.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>

namespace OpenMS
{

	using namespace std;

////////////////////////////////////////////////
// constructor for the object ConsensusIsotopePattern:
	ConsensusIsotopePattern::ConsensusIsotopePattern()
	{
	}

//////////////////////////////////////////////////
// class desctructor of ConsensusIsotopePattern
	ConsensusIsotopePattern::~ConsensusIsotopePattern()
	{

		isotopesTrace_.clear();
		mzIsotopesStDev_.clear();
		intensIsotopesStDev_.clear();
		rawIsotopes_.clear();

	}

//////////////////////////////////////////////////
// class copy constructor of ConsensusIsotopePattern
	ConsensusIsotopePattern::ConsensusIsotopePattern(const ConsensusIsotopePattern& tmp)
	{
		isotopesTrace_ = tmp.isotopesTrace_;
		mzIsotopesStDev_ = tmp.mzIsotopesStDev_;
		intensIsotopesStDev_ = tmp.intensIsotopesStDev_;
		rawIsotopes_ = tmp.rawIsotopes_;
	}

//////////////////////////////////////////////////
// class copy constructor of ConsensusIsotopePattern
	ConsensusIsotopePattern::ConsensusIsotopePattern(const ConsensusIsotopePattern* tmp)
	{
		isotopesTrace_ = tmp->isotopesTrace_;
		mzIsotopesStDev_ = tmp->mzIsotopesStDev_;
		intensIsotopesStDev_ = tmp->intensIsotopesStDev_;
		rawIsotopes_ = tmp->rawIsotopes_;
	}

//////////////////////////////////////////////////
// copy constructor:
	ConsensusIsotopePattern& ConsensusIsotopePattern::operator=(const ConsensusIsotopePattern& tmp)
	{
		isotopesTrace_ = tmp.isotopesTrace_;
		mzIsotopesStDev_ = tmp.mzIsotopesStDev_;
		intensIsotopesStDev_ = tmp.intensIsotopesStDev_;
		rawIsotopes_ = tmp.rawIsotopes_;
		return *this;
	}


/////////////////////////////////////////////////
// order an isotope trace in the correct cluster:
	void ConsensusIsotopePattern::addIsotopeTrace(double mz, double intens)
	{

		map<double, pair<vector<double>, vector<double> > >::iterator F = rawIsotopes_.lower_bound(mz);
		bool match = false;
		if (F != rawIsotopes_.end())
		{

			// compute the delta:
			if (SuperHirnUtil::compareMassValuesAtPPMLevel(mz, (*F).first, SuperHirnParameters::instance()->getMzTolPpm()))
			{
				(*F).second.first.push_back(mz);
				(*F).second.second.push_back(mz);
				match = true;
			}
			else if (F != rawIsotopes_.begin())
			{
				F--;
				if (SuperHirnUtil::compareMassValuesAtPPMLevel(mz, (*F).first, SuperHirnParameters::instance()->getMzTolPpm()))
				{
					(*F).second.first.push_back(mz);
					(*F).second.second.push_back(mz);
					match = true;
				}

			}

		}

		if (!match)
		{
			vector<double> mzTmp;
			mzTmp.push_back(mz);
			vector<double> intensTmp;
			intensTmp.push_back(intens);
			rawIsotopes_.insert(make_pair(mz, make_pair(mzTmp, intensTmp)));
		}

	}

/////////////////////////////////////////////////
// construc the consenus pattern:
	void ConsensusIsotopePattern::constructConsusPattern()
	{

		map<double, pair<vector<double>, vector<double> > >::iterator I = rawIsotopes_.begin();
		while (I != rawIsotopes_.end())
		{
			// condens a isotope peak trace:
			condensIsotopePattern(&(*I).second);
			I++;
		}

	}

// copied from simple_math
	pair<double, double> simple_math_AVERAGE_and_STDEV(vector<double>* IN)
	{

		double AVERAGE = 0;
		double STDEV = 0;

		if (IN->empty())
		{
			return make_pair(AVERAGE, STDEV);
		}

		if (IN->size() > 1)
		{
			vector<double>::iterator START = IN->begin();
			while (START != IN->end())
			{
				AVERAGE += (*START);
				START++;
			}
			AVERAGE /= double(IN->size());

			START = IN->begin();
			while (START != IN->end())
			{
				STDEV += ((AVERAGE - (*START)) * (AVERAGE - (*START)) );
				START++;
			}
			STDEV /= double(IN->size());
			STDEV = sqrt(STDEV);
			return make_pair(AVERAGE, STDEV);
		}
		else
		{
			return make_pair((*IN->begin()), 0.0);
		}
	}

//////////////////////////////////////////////////
// condens the pattern, make averge peaks from the traces:
	void ConsensusIsotopePattern::condensIsotopePattern(pair<vector<double>, vector<double> >* in)
	{

		// mz
		pair<double, double> mz = simple_math_AVERAGE_and_STDEV(&(in->first));
		// intens:
		pair<double, double> intens = simple_math_AVERAGE_and_STDEV(&(in->second));

		isotopesTrace_.insert(make_pair(mz.first, intens.first));
		mzIsotopesStDev_.push_back(mz.second);
		intensIsotopesStDev_.push_back(intens.second);

	}

}
