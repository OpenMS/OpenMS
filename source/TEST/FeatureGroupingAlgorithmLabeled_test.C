// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/FeatureGroupingAlgorithmLabeled.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(FeatureGroupingAlgorithmLabeled, "$Id FeatureFinder_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureGroupingAlgorithmLabeled* ptr = 0;
FeatureGroupingAlgorithmLabeled* nullPointer = 0;
START_SECTION((FeatureGroupingAlgorithmLabeled()))
	ptr = new FeatureGroupingAlgorithmLabeled();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureGroupingAlgorithmLabeled()))
	delete ptr;
END_SECTION

START_SECTION((static FeatureGroupingAlgorithm* create()))
	FeatureGroupingAlgorithm* ptr2 = 0;
  FeatureGroupingAlgorithm* base_NullPointer = 0;
	ptr2 = FeatureGroupingAlgorithmLabeled::create();
  TEST_NOT_EQUAL(ptr2, base_NullPointer)
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(FeatureGroupingAlgorithmLabeled::getProductName(),"labeled")
END_SECTION

START_SECTION((virtual void group(const std::vector< FeatureMap<> > &maps, ConsensusMap &out)))
	TOLERANCE_ABSOLUTE(0.001)
	
	FeatureGroupingAlgorithmLabeled fga;
	std::vector< FeatureMap<> > in;
	ConsensusMap out; 
	
	//test exception (no input)
	TEST_EXCEPTION(Exception::IllegalArgument, fga.group(in,out));
		
	//real test
	in.resize(1);
	in[0].resize(10);
	//start
	in[0][0].setRT(1.0f);
	in[0][0].setMZ(1.0f);
	in[0][0].setCharge(1);
	in[0][0].setOverallQuality(1);
	in[0][0].setIntensity(4.0f);
	//best
	in[0][1].setRT(1.5f);
	in[0][1].setMZ(5.0f);
	in[0][1].setCharge(1);
	in[0][1].setOverallQuality(1);
	in[0][1].setIntensity(2.0f);
	//inside (down, up, left, right)
	in[0][2].setRT(1.0f);
	in[0][2].setMZ(5.0f);
	in[0][2].setCharge(1);
	in[0][2].setOverallQuality(1);
	
	in[0][3].setRT(3.0f);
	in[0][3].setMZ(5.0f);
	in[0][3].setCharge(1);
	in[0][3].setOverallQuality(1);
	
	in[0][4].setRT(1.5f);
	in[0][4].setMZ(4.8f);
	in[0][4].setCharge(1);
	in[0][4].setOverallQuality(1);
	
	in[0][5].setRT(1.5f);
	in[0][5].setMZ(5.2f);
	in[0][5].setCharge(1);
	in[0][5].setOverallQuality(1);
	
	//outside (down, up, left, right)
	in[0][6].setRT(0.0f);
	in[0][6].setMZ(5.0f);
	in[0][6].setCharge(1);
	in[0][6].setOverallQuality(1);
	
	in[0][7].setRT(4.0f);
	in[0][7].setMZ(5.0f);
	in[0][7].setCharge(1);
	in[0][7].setOverallQuality(1);
	
	in[0][8].setRT(1.5f);
	in[0][8].setMZ(4.0f);
	in[0][8].setCharge(1);
	in[0][8].setOverallQuality(1);
	
	in[0][9].setRT(1.5f);
	in[0][9].setMZ(6.0f);
	in[0][9].setCharge(1);
	in[0][9].setOverallQuality(1);
	
	Param p;
	p.setValue("rt_estimate","false");
	p.setValue("rt_pair_dist",0.4);
	p.setValue("rt_dev_low",1.0);
	p.setValue("rt_dev_high",2.0);
	p.setValue("mz_pair_dists",DoubleList::create(4.0));
	p.setValue("mz_dev",0.6);
	fga.setParameters(p);


	//test exception (no file name set in out)
	TEST_EXCEPTION(Exception::IllegalArgument, fga.group(in,out));
	
	out.getFileDescriptions()[5].label = "light";
	out.getFileDescriptions()[5].filename = "filename";
	out.getFileDescriptions()[8] = out.getFileDescriptions()[5];
	out.getFileDescriptions()[8].label = "heavy";
	fga.group(in,out);
	
	TEST_EQUAL(out.size(),1)
	TEST_REAL_SIMILAR(out[0].getQuality(),0.959346);
	TEST_EQUAL(out[0].size(),2)
	ConsensusFeature::HandleSetType::const_iterator it = out[0].begin();
	TEST_REAL_SIMILAR(it->getMZ(),1.0f);
	TEST_REAL_SIMILAR(it->getRT(),1.0f);
	++it;
	TEST_REAL_SIMILAR(it->getMZ(),5.0f);
	TEST_REAL_SIMILAR(it->getRT(),1.5f);
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



