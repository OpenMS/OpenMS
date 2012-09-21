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

#include <stdio.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Fragment.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/ClusteredMS2ConsensusSpectrum.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/SUPERHIRN/MS2Feature.h>

namespace OpenMS
{

////////////////////////////////////////////////
// constructor for the object MS2Feature:
	MS2Feature::MS2Feature(MS2Fragment* in) :
			ClusteredMS2ConsensusSpectrum(in)
	{
		ID = -1;
	}

////////////////////////////////////////////////
// constructor for the object MS2ConsensusSpectrum:
	MS2Feature::MS2Feature(double iPrecursorMZ, double iTR, int iChrg, int iApexScan) :
			ClusteredMS2ConsensusSpectrum(iPrecursorMZ, iTR, iChrg, iApexScan)
	{
		ID = -1;
	}

//////////////////////////////////////////////////
// class desctructor of MS2Feature
	MS2Feature::~MS2Feature()
	{
	}

//////////////////////////////////////////////////
// class copy constructor of MS2Feature
	MS2Feature::MS2Feature(const MS2Feature& tmp) :
			ClusteredMS2ConsensusSpectrum(tmp)
	{
		ID = tmp.ID;
	}

//////////////////////////////////////////////////
// class copy constructor of MS2Feature
	MS2Feature::MS2Feature(const MS2Feature* tmp) :
			ClusteredMS2ConsensusSpectrum(tmp)
	{
		ID = tmp->ID;
	}

/////////////////////////////////////////////
// show info 
	void MS2Feature::show_info()
	{
		//printf("DELETED");
	}

}
