// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/KERNEL/ConsensusFeature.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusFeature, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusFeature* ptr = 0;
CHECK((ConsensusFeature()))
	ptr = new ConsensusFeature();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~ConsensusFeature()))
	delete ptr;
RESULT

Feature tmp_feature;
tmp_feature.setRT(1);
tmp_feature.setMZ(2);
tmp_feature.setIntensity(200);

Feature tmp_feature2;
tmp_feature2.setRT(2);
tmp_feature2.setMZ(3);
tmp_feature2.setIntensity(300);

Feature tmp_feature3;
tmp_feature3.setRT(3);
tmp_feature3.setMZ(4);
tmp_feature3.setIntensity(400);

CHECK((ConsensusFeature& operator=(const ConsensusFeature &rhs)))
  ConsensusFeature cons(tmp_feature);
  cons.insert(1,3,tmp_feature);
  
  ConsensusFeature cons_copy;
  cons_copy = cons;
  
  TEST_REAL_EQUAL(cons_copy.getRT(),1)
  TEST_REAL_EQUAL(cons_copy.getMZ(),2)
  TEST_REAL_EQUAL(cons_copy.getIntensity(),200)
  TEST_REAL_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_REAL_EQUAL((cons_copy.begin())->getElementIndex(),3)
  TEST_REAL_EQUAL((cons_copy.begin())->getIntensity(),200)
RESULT

CHECK((ConsensusFeature(const ConsensusFeature &rhs)))
  
  ConsensusFeature cons(tmp_feature);
  cons.insert(1,3,tmp_feature);
  ConsensusFeature cons_copy(cons);
  
  TEST_REAL_EQUAL(cons_copy.getRT(),1)
  TEST_REAL_EQUAL(cons_copy.getMZ(),2)
  TEST_REAL_EQUAL(cons_copy.getIntensity(),200)
  TEST_REAL_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_REAL_EQUAL((cons_copy.begin())->getElementIndex(),3)
  TEST_REAL_EQUAL((cons_copy.begin())->getIntensity(),200)
RESULT

CHECK((ConsensusFeature(const Peak2D &point)))
  
  ConsensusFeature cons(tmp_feature);
  TEST_REAL_EQUAL(cons.getRT(),1)
  TEST_REAL_EQUAL(cons.getMZ(),2)
  TEST_REAL_EQUAL(cons.getIntensity(),200)
  TEST_EQUAL(cons.empty(), true)
RESULT

CHECK((ConsensusFeature(UInt map_index, UInt element_index, const Feature &element)))
 	ConsensusFeature cons(1,3,tmp_feature);
  DRange<2> pos_range(1,2,1,2);
  DRange<1> int_range(200,200);
    
  TEST_REAL_EQUAL(cons.getRT(),1)
  TEST_REAL_EQUAL(cons.getMZ(),2)
  TEST_REAL_EQUAL(cons.getIntensity(),200)
  ConsensusFeature::HandleSetType::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),1)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getIntensity(),200)
RESULT

CHECK((DRange<1> getIntensityRange() const))
  ConsensusFeature cons;
  Feature f;
  f.setIntensity(0);
  cons.insert(0,0,f);
  f.setIntensity(200);
  cons.insert(0,1,f);
  
  TEST_REAL_EQUAL(cons.getIntensityRange().minX(),0.0)
  TEST_REAL_EQUAL(cons.getIntensityRange().maxX(),200.0)
RESULT

CHECK((DRange<2> getPositionRange() const))
  ConsensusFeature cons;
  Feature f;
  f.setRT(1.0);
  f.setMZ(500.0);  
  cons.insert(0,0,f);
  f.setRT(1000.0);
  f.setMZ(1500.0);  
  cons.insert(0,1,f);
  
  TEST_REAL_EQUAL(cons.getPositionRange().minX(),1.0)
  TEST_REAL_EQUAL(cons.getPositionRange().maxX(),1000.0)
  TEST_REAL_EQUAL(cons.getPositionRange().minY(),500.0)
  TEST_REAL_EQUAL(cons.getPositionRange().maxY(),1500.0)
RESULT

CHECK((const HandleSetType& getFeatures() const))
  ConsensusFeature cons;
  cons.insert(2,3,tmp_feature);
  const ConsensusFeature cons_copy(cons);
  
  ConsensusFeature::HandleSetType group = cons_copy.getFeatures();
    
  ConsensusFeature::HandleSetType::const_iterator it = group.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getIntensity(),200)
RESULT


CHECK((void insert(FeatureHandle const &handle)))
  ConsensusFeature cons;
  FeatureHandle h1(2,3,tmp_feature);
  FeatureHandle h2(4,5,tmp_feature);
  cons.insert(h1);
  cons.insert(h2);
      
  ConsensusFeature::HandleSetType::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getIntensity(),200)
  ++it;
  TEST_REAL_EQUAL(it->getMapIndex(),4)
  TEST_REAL_EQUAL(it->getElementIndex(),5)
  TEST_REAL_EQUAL(it->getIntensity(),200)
  ++it;
  TEST_EQUAL(it==cons.end(), true)
RESULT

CHECK((void insert(UInt map_index, UInt element_index, const Feature &element)))
  ConsensusFeature cons;
  cons.insert(2,3,tmp_feature);
      
  ConsensusFeature::HandleSetType::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getIntensity(),200)
  ++it;
  TEST_EQUAL(it==cons.end(),true)
RESULT

CHECK((DoubleReal getQuality() const))
	ConsensusFeature cons;
	TEST_REAL_EQUAL(cons.getQuality(),0.0)
RESULT

CHECK((void setQuality(DoubleReal quality)))
	ConsensusFeature cons;
	cons.setQuality(4.5);
	TEST_REAL_EQUAL(cons.getQuality(),4.5)
RESULT

CHECK((void computeConsensus()))
  ConsensusFeature cons;
  //one point
  cons.insert(2,3,tmp_feature);
	cons.computeConsensus();
	TEST_REAL_EQUAL(cons.getIntensity(),200)
	TEST_REAL_EQUAL(cons.getRT(),1)
	TEST_REAL_EQUAL(cons.getMZ(),2)
	//two points
  cons.insert(4,5,tmp_feature2);
	cons.computeConsensus();
	TEST_REAL_EQUAL(cons.getIntensity(),250)
	TEST_REAL_EQUAL(cons.getRT(),1.5)
	TEST_REAL_EQUAL(cons.getMZ(),2.5)	
	//three points
  cons.insert(6,7,tmp_feature3);
	cons.computeConsensus();
	TEST_REAL_EQUAL(cons.getIntensity(),300)
	TEST_REAL_EQUAL(cons.getRT(),2)
	TEST_REAL_EQUAL(cons.getMZ(),3)	
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



