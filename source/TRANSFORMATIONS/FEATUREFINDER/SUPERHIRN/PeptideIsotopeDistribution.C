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
//  PEAK DETECTION OF FOURIER TRANSFORME MS INSTRUMENT DATA
//
//  written by Markus Mueller, markus.mueller@imsb.biol.ethz.ch
//  ( and Lukas Mueller, Lukas.Mueller@imsb.biol.ethz.ch)
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
#include <string.h>
#include <stdio.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/PeptideIsotopeDistribution.h>

namespace OpenMS
{

	using namespace std;

////////////////////////////////////////////////
// constructor for the object PeptideIsotopeDisribution:
	PeptideIsotopeDisribution::PeptideIsotopeDisribution(vector<double> iMZ, vector<double> iIntens)
	{
		this->mass = iMZ;
		this->intens = iIntens;
		chargeState = 2;
		id = -1;
		intensArray = NULL;
	}

////////////////////////////////////////////////
// constructor for the object PeptideIsotopeDisribution:
	PeptideIsotopeDisribution::PeptideIsotopeDisribution(vector<double> iMZ, vector<double> iIntens, int iZ, string iName,
			string iSQ, int iID)
	{
		this->mass = iMZ;
		this->intens = iIntens;
		this->name = iName;
		this->sq = iSQ;
		this->chargeState = iZ;
		this->id = iID;
		intensArray = NULL;

		this->constructSummaryString();
	}

////////////////////////////////////////////////
// constructor for the object PeptideIsotopeDisribution:
	PeptideIsotopeDisribution::PeptideIsotopeDisribution(vector<double> iMZ, vector<double> iIntens, int iZ, string iName,
			string iSQ, int iID, double rtSeg)
	{
		this->mass = iMZ;
		this->intens = iIntens;
		this->name = iName;
		this->sq = iSQ;
		this->chargeState = iZ;
		this->id = iID;
		this->RtSegment = rtSeg;
		intensArray = NULL;

		this->constructSummaryString();
	}

//////////////////////////////////////////////////
// class desctructor of PeptideIsotopeDisribution
	PeptideIsotopeDisribution::~PeptideIsotopeDisribution()
	{
		this->mass.clear();
		this->intens.clear();

		if (intensArray != NULL)
		{
			delete intensArray;
			intensArray = NULL;
		}

	}

//////////////////////////////////////////////////
// class copy constructor of PeptideIsotopeDisribution
	PeptideIsotopeDisribution::PeptideIsotopeDisribution(const PeptideIsotopeDisribution& tmp)
	{
		this->mass = tmp.mass;
		this->intens = tmp.intens;
		this->name = tmp.name;
		this->id = tmp.id;
		this->sq = tmp.sq;
		this->summary = tmp.summary;
		this->chargeState = tmp.chargeState;
		this->RtSegment = tmp.RtSegment;
		this->RtEnd = tmp.RtEnd;
		this->RtStart = tmp.RtStart;

		this->intensArray = NULL;
		if (tmp.intensArray != NULL)
		{
			intensArray = new double[sizeof(tmp.intensArray)];
			memcpy(intensArray, tmp.intensArray, sizeof(tmp.intensArray) * sizeof(double));
		}
	}

//////////////////////////////////////////////////
// class copy constructor of PeptideIsotopeDisribution
	PeptideIsotopeDisribution::PeptideIsotopeDisribution(const PeptideIsotopeDisribution* tmp)
	{
		this->mass = tmp->mass;
		this->intens = tmp->intens;
		this->name = tmp->name;
		this->sq = tmp->sq;
		this->id = tmp->id;
		this->summary = tmp->summary;
		this->chargeState = tmp->chargeState;
		this->RtSegment = tmp->RtSegment;
		this->RtEnd = tmp->RtEnd;
		this->RtStart = tmp->RtStart;

		this->intensArray = NULL;
		if (tmp->intensArray != NULL)
		{
			intensArray = new double[sizeof(tmp->intensArray)];
			memcpy(intensArray, tmp->intensArray, sizeof(tmp->intensArray) * sizeof(double));
		}
	}

//////////////////////////////////////////////////
// copy constructor:
	PeptideIsotopeDisribution& PeptideIsotopeDisribution::operator=(const PeptideIsotopeDisribution& tmp)
	{
		this->mass = tmp.mass;
		this->intens = tmp.intens;
		this->name = tmp.name;
		this->id = tmp.id;
		this->sq = tmp.sq;
		this->summary = tmp.summary;
		this->RtSegment = tmp.RtSegment;
		this->RtEnd = tmp.RtEnd;
		this->RtStart = tmp.RtStart;
		this->chargeState = tmp.chargeState;
		if (intensArray != NULL)
		{
//			delete intensArray; // PK this is passed outside, if it is deleted here it might cause failure later
			intensArray = NULL;
		}

		if (tmp.intensArray != NULL)
		{
			intensArray = new double[sizeof(tmp.intensArray)];
			memcpy(intensArray, tmp.intensArray, sizeof(tmp.intensArray) * sizeof(double));
		}

		return *this;
	}

//////////////////////////////////////////////////
// construct an information string:
	void PeptideIsotopeDisribution::constructSummaryString()
	{
		char tmp[100];
		sprintf(tmp, "INFO:%s_%d;%s", name.c_str(), id, sq.c_str());
		this->summary = tmp;
	}

//////////////////////////////////////////////////
// get the intensity values in form of an array
	double* PeptideIsotopeDisribution::getIntensityArray()
	{

		if (intensArray == NULL && !intens.empty())
		{

			intensArray = new double[intens.size()];int i = 0;
			vector< double>::iterator P = intens.begin();
			while( P != intens.end() )
			{
				intensArray[i] = *P;
				// cout<<intensArray[i]<<endl;
					i++;
					P++;
				}
			}

			return intensArray;
		}

//////////////////////////////////////////////////
// show info:
	void PeptideIsotopeDisribution::show_info()
	{

		//printf( "\n \t External Isotope Distribution:  %d isotopes", mass.size() );
		printf("\n \t Info: Name=%s, SQ=%s, Id=%d, z=+%d, rtSeg=%.0f \n", name.c_str(), sq.c_str(), id, chargeState,
				RtSegment);

		vector<double>::iterator M = mass.begin();
		vector<double>::iterator I = intens.begin();

		while (M != mass.end())
		{

			printf("\t \t isotope m/z %0.2f - %0.2f \n", *M, *I);

			M++;
			I++;
		}

	}

////////////////////////////////////////////////
// return info in form or a string about PeptideIsotopeDisribution:
	string PeptideIsotopeDisribution::getIsotopeDistInfo()
	{
		string info;
		info = "Name: " + this->name;
		info += ", SQ: " + this->sq;
		return info;
	}
}
