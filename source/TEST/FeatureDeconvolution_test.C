// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DECHARGING/FeatureDeconvolution.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <CoinMessageHandler.hpp>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureDeconvolution, "$Id:$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureDeconvolution* ptr = 0;
START_SECTION(FeatureDeconvolution())
	ptr = new FeatureDeconvolution();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~FeatureDeconvolution())
	delete ptr;
END_SECTION

START_SECTION(void compute)
//_CrtSetDbgFlag(_CrtSetDbgFlag(0)|_CRTDBG_CHECK_ALWAYS_DF);


	CoinMessageHandler a;
	{
	CoinMessageHandler bb;
	}

	/*CoinMessageHandler * handler_ =	new CoinMessageHandler();
	delete handler_;
	handler_ =	new CoinMessageHandler();
	delete handler_;
*/

	FeatureDeconvolution fd;
	FeatureMap<> fm;
	ConsensusMap cm, cm2;
	FeatureXMLFile fl;
	fl.load("/home/chris/svn2/SIM_EASY_INPUT.FEATUREXML",fm);
	fd.compute(fm, cm, cm2);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


