// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

///////////////////////////
#include <OpenMS/KERNEL/ConsensusFeature.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusFeature, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusFeature<>* ptr = 0;
CHECK(ConsensusFeature())
	ptr = new ConsensusFeature<>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~ConsensusFeature())
	delete ptr;
RESULT

CHECK(ConsensusFeature& operator=(const ConsensusFeature& source))
  DPosition<2> pos(1,2);
  ConsensusFeature<> cons(pos,200);
  DFeature<2> feat;
  feat.setPosition(pos);
  feat.setIntensity(200);
  
  IndexTuple<> ind(1,3,feat);
  cons.insert(ind);
  
  ConsensusFeature<> cons_copy;
  cons_copy = cons;
  
  TEST_REAL_EQUAL(cons_copy.getPosition()[0],1)
  TEST_REAL_EQUAL(cons_copy.getPosition()[1],2)
  TEST_REAL_EQUAL(cons_copy.getIntensity(),200)
  TEST_EQUAL(cons_copy.getPositionRange() == cons.getPositionRange(), true)
  TEST_EQUAL(cons_copy.getIntensityRange() == cons.getIntensityRange(), true)
  TEST_REAL_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_REAL_EQUAL((cons_copy.begin())->getElementIndex(),3)
  TEST_REAL_EQUAL((cons_copy.begin())->getElement().getIntensity(),200)
RESULT

CHECK((ConsensusFeature(const ConsensusFeature& c_feature_1, const ConsensusFeature& c_feature_2)))
  DPosition<2> pos(1,2);
  ConsensusFeature<> cons1(pos,200);
  DFeature<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  IndexTuple<> ind1(1,3,feat1);
  cons1.insert(ind1);
  
  pos[0]=2;
  pos[1]=3;
  ConsensusFeature<> cons2(pos,200);
  DFeature<2> feat2;
  feat2.setPosition(pos);
  feat2.setIntensity(200);
  IndexTuple<> ind2(2,3,feat2);
  cons2.insert(ind2);
  
  ConsensusFeature<> cons3(cons1,cons2);
  DRange<2> pos_range(1,2,2,3);
  DRange<1> int_range(200,200);
    
  TEST_REAL_EQUAL(cons3.getPosition()[0],1.5)
  TEST_REAL_EQUAL(cons3.getPosition()[1],2.5)
  TEST_REAL_EQUAL(cons3.getIntensity(),200)
  TEST_EQUAL(cons3.getPositionRange() == pos_range, true)
  TEST_EQUAL(cons3.getIntensityRange() == int_range, true)
  ConsensusFeature<>::Group::const_iterator it = cons3.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),1)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
  ++it;
  TEST_REAL_EQUAL(it->getMapIndex(),2)
 	TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(ConsensusFeature(const ConsensusFeature& source))
  DPosition<2> pos(1,2);
  ConsensusFeature<> cons(pos,200);
  DFeature<2> feat;
  feat.setPosition(pos);
  feat.setIntensity(200);
  IndexTuple<> ind(1,3,feat);
  cons.insert(ind);
  ConsensusFeature<> cons_copy(cons);
  
  TEST_REAL_EQUAL(cons_copy.getPosition()[0],1)
  TEST_REAL_EQUAL(cons_copy.getPosition()[1],2)
  TEST_REAL_EQUAL(cons_copy.getIntensity(),200)
  TEST_EQUAL(cons_copy.getPositionRange() == cons.getPositionRange(), true)
  TEST_EQUAL(cons_copy.getIntensityRange() == cons.getIntensityRange(), true)
  TEST_REAL_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_REAL_EQUAL((cons_copy.begin())->getElementIndex(),3)
  TEST_REAL_EQUAL((cons_copy.begin())->getElement().getIntensity(),200)
RESULT

CHECK((ConsensusFeature(const PositionType& pos, const IntensityType& i)))
  DPosition<2> pos(1,2);
  ConsensusFeature<> cons(pos,200);
    
  DRange<2> pos_range;
  DRange<1> int_range;
  TEST_REAL_EQUAL(cons.getPosition()[0],1)
  TEST_REAL_EQUAL(cons.getPosition()[1],2)
  TEST_REAL_EQUAL(cons.getIntensity(),200)
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
  TEST_EQUAL(cons.isEmpty(), true)
RESULT

CHECK((ConsensusFeature(const UnsignedInt& map_1_index, const UnsignedInt& feature_index_1, const ElementType& feature_1, const UnsignedInt& map_2_index, const UnsignedInt& feature_index_2, const ElementType& feature_2)))
  DPosition<2> pos(1,2);
  DFeature<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  
  pos[0]=2;
  pos[1]=3;
  DFeature<2> feat2;
  feat2.setPosition(pos);
  feat2.setIntensity(200);
  
  ConsensusFeature<> cons3(1,3,feat1,2,3,feat2);
  DRange<2> pos_range(1,2,2,3);
  DRange<1> int_range(200,200);
    
  TEST_REAL_EQUAL(cons3.getPosition()[0],1.5)
  TEST_REAL_EQUAL(cons3.getPosition()[1],2.5)
  TEST_REAL_EQUAL(cons3.getIntensity(),200)
  //std::cout << "cons3.getPositionRange() " << cons3.getPositionRange() << std::endl;
  //std::cout << "pos_range " << pos_range << std::endl;
  TEST_EQUAL(cons3.getPositionRange() == pos_range, true)
  TEST_EQUAL(cons3.getIntensityRange() == int_range, true)
  ConsensusFeature<>::Group::const_iterator it = cons3.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),1)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
  ++it;
  TEST_REAL_EQUAL(it->getMapIndex(),2)
 	TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK((ConsensusFeature(const UnsignedInt& map_index, const UnsignedInt& feature_index, const ElementType& feature)))
  DPosition<2> pos(1,2);
  DFeature<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  
 	ConsensusFeature<> cons(1,3,feat1);
  DRange<2> pos_range(1,2,1,2);
  DRange<1> int_range(200,200);
    
  TEST_REAL_EQUAL(cons.getPosition()[0],1)
  TEST_REAL_EQUAL(cons.getPosition()[1],2)
  TEST_REAL_EQUAL(cons.getIntensity(),200)
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
  ConsensusFeature<>::Group::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),1)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK((ConsensusFeature(const UnsignedInt& map_index, const UnsignedInt& feature_index, const ElementType& feature, const ConsensusFeature& c_feature)))
  DPosition<2> pos(1,2);
  DFeature<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
    
  pos[0]=2;
  pos[1]=3;
  ConsensusFeature<> cons2(pos,200);
  DFeature<2> feat2;
  feat2.setPosition(pos);
  feat2.setIntensity(200);
  IndexTuple<> ind2(2,3,feat2);
  cons2.insert(ind2);
  
  ConsensusFeature<> cons3(1,3,feat1,cons2);
  DRange<2> pos_range(1,2,2,3);
  DRange<1> int_range(200,200);
    
  TEST_REAL_EQUAL(cons3.getPosition()[0],1.5)
  TEST_REAL_EQUAL(cons3.getPosition()[1],2.5)
  TEST_REAL_EQUAL(cons3.getIntensity(),200)
  TEST_EQUAL(cons3.getPositionRange() == pos_range, true)
  TEST_EQUAL(cons3.getIntensityRange() == int_range, true)
  ConsensusFeature<>::Group::const_iterator it = cons3.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),1)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
  ++it;
  TEST_REAL_EQUAL(it->getMapIndex(),2)
 	TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(Group& getFeatures())
	DPosition<2> pos(1,2);
  DFeature<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  IndexTuple<> ind(2,3,feat1);
  
  ConsensusFeature<>::Group group;
  group.insert(ind);
  
  ConsensusFeature<> cons;
  cons.getFeatures() = group;
    
  ConsensusFeature<>::Group::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(IntensityBoundingBoxType& getIntensityRange())
  DRange<1> int_range(0,200);
  ConsensusFeature<> cons;
  cons.getIntensityRange() = int_range;
  
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
RESULT

CHECK(PositionBoundingBoxType& getPositionRange())
  DRange<2> pos_range(0,1,100,200);
  ConsensusFeature<> cons;
  cons.getPositionRange() = pos_range;
  
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
RESULT

CHECK(const Group& getFeatures() const)
  DPosition<2> pos(1,2);
  DFeature<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  IndexTuple<> ind(2,3,feat1);
  ConsensusFeature<> cons;
  cons.insert(ind);
  const ConsensusFeature<> cons_copy(cons);
  
  ConsensusFeature<>::Group group = cons_copy.getFeatures();
    
  ConsensusFeature<>::Group::const_iterator it = group.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(const IntensityBoundingBoxType& getIntensityRange() const)
  DRange<1> int_range;
  const ConsensusFeature<> cons;
  
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
RESULT

CHECK(const PositionBoundingBoxType& getPositionRange() const)
  DRange<2> pos_range;
  const ConsensusFeature<> cons;
  
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
RESULT

CHECK(void insert(const IndexTuple& tuple))
  DPosition<2> pos(1,2);
  DFeature<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  IndexTuple<> ind(2,3,feat1);
  
  ConsensusFeature<> cons;
  cons.insert(ind);
      
  ConsensusFeature<>::Group::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(void setFeatures(const Group& g))
  DPosition<2> pos(1,2);
  DFeature<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  IndexTuple<> ind(2,3,feat1);
  
  ConsensusFeature<>::Group group;
  group.insert(ind);
  
  ConsensusFeature<> cons;
  cons.setFeatures(group);
    
  ConsensusFeature<>::Group::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(void setIntensityRange(const IntensityBoundingBoxType& i))
  DRange<1> int_range(0,200);
  ConsensusFeature<> cons;
  cons.setIntensityRange(int_range);
  
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
RESULT

CHECK(void setPositionRange(const PositionBoundingBoxType& p))
  DRange<2> pos_range(0,1,100,200);
  ConsensusFeature<> cons;
  cons.setPositionRange(pos_range);
  
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



