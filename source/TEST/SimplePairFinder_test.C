// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>


///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/SimplePairFinder.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DPosition <2> PositionType;


START_TEST(SimplePairFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SimplePairFinder* ptr = 0;
CHECK((SimplePairFinder()))
	ptr = new SimplePairFinder();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~SimplePairFinder()))
	delete ptr;
RESULT

CHECK((double getDiffExponent(UInt dim)))
  SimplePairFinder spf;
  
  TEST_REAL_EQUAL(spf.getDiffExponent(0),2)
  TEST_REAL_EQUAL(spf.getDiffExponent(1),1)
RESULT

CHECK((double getPairMinQuality()))
  SimplePairFinder spf;
  
  TEST_REAL_EQUAL(spf.getPairMinQuality(),0.01)
RESULT

CHECK((void setDiffExponent(UInt dim, DoubleReal exponent)))
  SimplePairFinder spf;
  spf.setDiffExponent(0,20);
  spf.setDiffExponent(1,25);
  
  TEST_REAL_EQUAL(spf.getDiffExponent(0),20)
  TEST_REAL_EQUAL(spf.getDiffExponent(1),25)
RESULT

CHECK((void setPairMinQuality(DoubleReal quality)))
  SimplePairFinder spf;
  spf.setPairMinQuality(0.9);
  
  TEST_REAL_EQUAL(spf.getPairMinQuality(),0.9)
RESULT

CHECK((static BasePairFinder<ConsensusMap>* create()))
	BasePairFinder* base_ptr = 0;
	base_ptr = SimplePairFinder::create();
	TEST_NOT_EQUAL(base_ptr, 0)
RESULT

CHECK((static const String getProductName()))
  SimplePairFinder spf;
  
  TEST_EQUAL(spf.getName() == "simple",true)
RESULT

CHECK((virtual void findElementPairs()))
  FeatureMap<> scene;
  Feature feat1;
  Feature feat2;
  Feature feat3;
  PositionType pos1(0,0);
  PositionType pos2(200,300);
  PositionType pos3(400,500);
  feat1.setPosition(pos1);
  feat1.setIntensity(100);
  feat2.setPosition(pos2);
  feat2.setIntensity(300);
  feat3.setPosition(pos3);
  feat3.setIntensity(400);
  scene.push_back(feat1);
  scene.push_back(feat2);
  scene.push_back(feat3);
  
  FeatureMap<> model;
  Feature feat4;
  Feature feat5;
  Feature feat6;
  PositionType pos4(4,4);
  PositionType pos5(204,304);
  PositionType pos6(404,504);
  feat4.setPosition(pos4);
  feat4.setIntensity(100);
  feat5.setPosition(pos5);
  feat5.setIntensity(300);
  feat6.setPosition(pos6);
  feat6.setIntensity(400);
  model.push_back(feat4);
  model.push_back(feat5);
  model.push_back(feat6);
  
  SimplePairFinder spf;
	ConsensusMap model2;
	SimplePairFinder::convert(0,model,model2);
	spf.setModelMap(0,model2);
	ConsensusMap scene2;
	SimplePairFinder::convert(1,scene,scene2);
	spf.setSceneMap(1,scene2);
	SimplePairFinder::ElementPairVectorType pairs;
  spf.setElementPairs(pairs);
  spf.findElementPairs();
    
  TEST_EQUAL((pairs.begin())->first == feat1, true)
  TEST_EQUAL((pairs.begin())->second == feat4, true)
  TEST_EQUAL((pairs.begin()+1)->first == feat2, true)
  TEST_EQUAL((pairs.begin()+1)->second == feat5, true)
  TEST_EQUAL((pairs.begin()+2)->first == feat3,true)
  TEST_EQUAL((pairs.begin()+2)->second == feat6,true)
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



