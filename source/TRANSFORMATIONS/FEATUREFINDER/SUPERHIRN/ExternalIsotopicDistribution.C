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
//
//  by Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch
//  October 2005
//  
//  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
//  December 2010
//
//  Group of Prof. Ruedi Aebersold, IMSB, ETH Hoenggerberg, Zurich
// 
//

#include <string>
#include <vector>
#include <map>
#include <math.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnUtil.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/PeptideIsotopeDistribution.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ExternalIsotopicDistribution.h>

namespace OpenMS
{

	using namespace std;

	bool ExternalIsotopicDistribution::EXTERNAL_ISOTOPIC_PROFILES = false;
	string ExternalIsotopicDistribution::XMLInputFile;
	double ExternalIsotopicDistribution::EXTERNAL_DISTRIBUTION_MONO_ISOTOPE_PPM_TOLERANCE;
	multimap<double, PeptideIsotopeDisribution> ExternalIsotopicDistribution::allExternalPepdistributions;

////////////////////////////////////////////////
// constructor for the object ExternalIsotopicDistribution:
	ExternalIsotopicDistribution::ExternalIsotopicDistribution()
	{
	}

//////////////////////////////////////////////////
// class desctructor of ExternalIsotopicDistribution
	ExternalIsotopicDistribution::~ExternalIsotopicDistribution()
	{
	}


//////////////////////////////////////////////////
// function to extract external isotopic profiles
// by an input monoisotopic mass
	PeptideIsotopeDisribution* ExternalIsotopicDistribution::extractExternalIsotopicProfile(double monoMass, int charge,
			double RT)
	{

		multimap<double, PeptideIsotopeDisribution>::iterator F = allExternalPepdistributions.lower_bound(monoMass);
		multimap<double, PeptideIsotopeDisribution>::iterator P = F;

		// search down:
		do
		{

			if (checkMonoIsotopOverlap(monoMass, charge, RT, &(P->second)))
			{
				return &(P->second);
			}

			P--;
		}
		while (P != allExternalPepdistributions.begin());

		// search up:
		// PK: does not make sense, F iterator is never used, just P iterator which is fixed to begin().
		while (F != allExternalPepdistributions.end())
		{
			if (checkMonoIsotopOverlap(monoMass, charge, RT, &(P->second)))
			{
				return &(P->second);
			}
			F++;
		}
		return NULL;

	}

//////////////////////////////////////////////////
//  check if two masses to be the same
	bool ExternalIsotopicDistribution::checkMonoIsotopOverlap(double MeasMonoMass, int TestCharge, double TR,
			PeptideIsotopeDisribution* dist)
	{

		// check the charge state:
		if (TestCharge != dist->getChargeState())
		{
			return false;
		}

		// check the mass tolerance:
		if (!SuperHirnUtil::compareMassValuesAtPPMLevel(dist->getMonoMass(), MeasMonoMass,
				ExternalIsotopicDistribution::EXTERNAL_DISTRIBUTION_MONO_ISOTOPE_PPM_TOLERANCE))
		{
			return false;
		}

		if ((TR < dist->getRTStart()) || (TR > dist->getRTEnd()))
		{
			return false;
		}

		return true;
	}

//////////////////////////////////////////////////
// init the retention time segments
	void ExternalIsotopicDistribution::initRetentionTimeSegments(double start, double end)
	{

		if (!allExternalPepdistributions.empty())
		{

			double maxS = 0;
			multimap<double, PeptideIsotopeDisribution>::iterator P = allExternalPepdistributions.begin();
			while (P != allExternalPepdistributions.end())
			{

				if (P->second.getRTSegment() > maxS)
				{
					maxS = P->second.getRTSegment();
				}

				P++;
			}

			// divide into segments:
			double Ssize = (end - start) / maxS;

			P = allExternalPepdistributions.begin();
			while (P != allExternalPepdistributions.end())
			{

				double s = P->second.getRTSegment();

				double lStart = (s - 1.0) * Ssize;
				P->second.setRTStart(lStart);
				double lEnd = (s) * Ssize;
				P->second.setRTEnd(lEnd);

				P++;
			}
		}

	}
}
