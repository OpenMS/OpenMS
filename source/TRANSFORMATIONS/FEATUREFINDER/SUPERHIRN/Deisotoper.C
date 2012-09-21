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
/*
 *  Deisotoper.cpp
 *  PeakDetection
 *
 *  Created by Markus Mueller on 10/19/06.
 *  
 *  Ported to OpenMS by Florian Zeller, florian.zeller@bsse.ethz.ch
 *  December 2010
 *
 */

#include <list>
#include <iostream>
#include <map>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/Deisotoper.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidPeak.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/CentroidData.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/IsotopicDist.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/SuperHirnParameters.h>

namespace OpenMS
{

	using namespace std;

//	int Deisotoper::sfMinCharge = 1; // minimum charge considered in peak detection
//	int Deisotoper::sfMaxCharge = 5; // maximum charge considered in peak detection

// Default constructor
	Deisotoper::Deisotoper()
	{
		fMinPeakGroupSize = 0;
		fTheta = 0.0;
		fScanNumber = -1;
	}

// Constructor. Takes centroide values and deisotopes them 
	Deisotoper::Deisotoper(CentroidData& pCentroidData) // Data objects containing centroid values
	{
		go(pCentroidData);
	}

// destructor
	Deisotoper::~Deisotoper()
	{
		fDeconvPeaks.clear();
	}

// Operators

// Writes data to out stream using the << operator
	ostream& operator<<(ostream& pOut, // output stream
			Deisotoper& pDeisotoper) //
	{
		list<DeconvPeak> p;
		list<DeconvPeak>::iterator pi;

		p = pDeisotoper.getDeconvPeaks();

		for (pi = p.begin(); pi != p.end(); ++pi)
		{
			if (pDeisotoper.getShortReportFlag())
			{
				pOut << *((CentroidPeak*) &(*pi)) << endl;
			}
			else
			{
				pOut << *pi << " " << pDeisotoper.getScanNumber() << endl;
			}
		}

		return pOut;
	}

// Takes centroide values and deisotopes them 
	void Deisotoper::go(CentroidData& pCentroidData) // Data objects containing centroid values
	{
		int cnt, charge;
		double alpha;
		bool matched;
		list<CentroidPeak> centroidPeaks;
		list<CentroidPeak>::iterator start, end, pi;
		list<list<CentroidPeak>::iterator> matchedPeaks;

		pCentroidData.get(centroidPeaks);

		fMinPeakGroupSize = 2;

		if (SuperHirnParameters::instance()->getMinIntensity() < SuperHirnParameters::instance()->getIntensityFloor())
		{
			pCentroidData.setNoise(30.0); // set noise level at 30 prcentile
			fTheta = pCentroidData.getNoise();
		}
		else
		{
			fTheta = SuperHirnParameters::instance()->getMinIntensity();
		}

		pCentroidData.resetPeakGroupIter();
		while (pCentroidData.getNextPeakGroup(start, end))
		{ // isotopic patterns are withing the same peak group
			for (cnt = 0, pi = start; pi != end; ++pi, ++cnt)
			{
			};

			if (cnt >= fMinPeakGroupSize)
			{ // Discard peak groups with only one peak
				for (pi = start; pi != end; ++pi, --cnt)
				{
					if (pi->getIntensity() < fTheta || cnt < fMinPeakGroupSize)
						continue; // skip small peaks
					/*
					 if( CentroidData::MonoIsoDebugging ){
					 if( ( CentroidData::DebugMonoIsoMassMin <= pi->getMass()) && ( CentroidData::DebugMonoIsoMassMax >= pi->getMass()) ){
					 cout<<"To deisotope: "<<pi->getMass()<<endl;
					 }
					 }
					 */

					for (charge = SuperHirnParameters::instance()->getMaxFeatureChrg(); charge >= SuperHirnParameters::instance()->getMinFeatureChrg(); --charge)
					{

						matched = IsotopicDist::getMatchingPeaks(pi, end, charge, alpha, fTheta, matchedPeaks); // get peak that match isotopic pattern of charge
						if (matched && pi->getIntensity() >= fTheta)
						{ // subtract isotopic match from peaks if match is significant

							/*
							 if( CentroidData::MonoIsoDebugging ){
							 if( ( CentroidData::DebugMonoIsoMassMin <= pi->getMass()) && ( CentroidData::DebugMonoIsoMassMax >= pi->getMass()) ){
							 cout<<"matched to +"<<charge<<endl;
							 }
							 }
							 */

							DeconvPeak mono(pi->getMass(), 0.0, charge, 0, 0.0, 0.0);
							if (!pi->getExtraPeakInfo().empty())
							{
								mono.setExtraPeakInfo(pi->getExtraPeakInfo());
							}

							IsotopicDist::subtractMatchingPeaks(matchedPeaks, charge, alpha, mono);
							fDeconvPeaks.push_back(mono);

							/*
							 if( CentroidData::MonoIsoDebugging ){
							 if( ( CentroidData::DebugMonoIsoMassMin <= pi->getMass()) && ( CentroidData::DebugMonoIsoMassMax >= pi->getMass()) ){
							 mono.show_info();
							 }
							 }
							 */
						}
						matchedPeaks.clear();
					}
				}
			}
		}
	}

// removes spooky or very small monoisotopic peaks
	void Deisotoper::cleanDeconvPeaks()
	{
		double tol, mass;
		list<DeconvPeak>::iterator pi, beg, end, most_intense;

		for (pi = fDeconvPeaks.begin(); pi != fDeconvPeaks.end(); ++pi)
		{
			beg = pi;
			mass = pi->getMass();
			most_intense = pi;
			tol = SuperHirnParameters::instance()->getMassTolPpm() * mass / 1.0e6
					+ SuperHirnParameters::instance()->getMassTolDa();
			++pi;
			for (; pi != fDeconvPeaks.end(); ++pi)
			{ // cluster peaks and define max intensity within cluster
				if (pi->getMass() > mass + 2.0 * tol)
				{
					break;
				}

				if (most_intense->getIntensity() < pi->getIntensity())
				{
					most_intense = pi; // store most intense peak
				}
			}
			end = pi;

			for (pi = beg; pi != fDeconvPeaks.end() && pi != end; ++pi)
			{ // remove all 'very' small peak within cluster
//			cout << "remove: " << pi->getMass() << " " << pi->getIntensity() << " " << pi->getCharge() << " | " << most_intense->getMass() << " " << most_intense->getIntensity() << endl;
				if (most_intense->getIntensity() > 2.0 * pi->getIntensity())
				{
//				cout << "remove: " << pi->getMass() << " " << pi->getIntensity() << endl;
					pi = fDeconvPeaks.erase(pi);
					if (pi != fDeconvPeaks.begin()) // FLO: Fix windows error (crash "could not decrement")
						--pi;
				}
			}

			if (pi != fDeconvPeaks.end())
			{
				--pi;
			}
			else
			{
				break;
			}
		}
	}
}
