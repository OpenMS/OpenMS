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
// $Maintainer: Eva Lange$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/KERNEL/ConsensusPeak.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusPeak, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusPeak<>* ptr = 0;
CHECK(ConsensusPeak())
	ptr = new ConsensusPeak<>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~ConsensusPeak())
	delete ptr;
RESULT

CHECK(ConsensusPeak& operator=(const ConsensusPeak& source))
  DPosition<2> pos(1,2);
  ConsensusPeak<> cons(pos,200);
  DPeak<2> feat;
  feat.setPosition(pos);
  feat.setIntensity(200);
  
  IndexTuple<DPeakArray<2, Peak2D> > ind(1,3,feat);
  cons.insert(ind);
  
  ConsensusPeak<> cons_copy;
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

CHECK((ConsensusPeak(const ConsensusPeak& c_peak_1, const ConsensusPeak& c_peak_2)))
  DPosition<2> pos(1,2);
  ConsensusPeak<> cons1(pos,200);
  DPeak<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  IndexTuple<DPeakArray<2, Peak2D> > ind1(1,3,feat1);
  cons1.insert(ind1);
  
  pos[0]=2;
  pos[1]=3;
  ConsensusPeak<> cons2(pos,200);
  DPeak<2> feat2;
  feat2.setPosition(pos);
  feat2.setIntensity(200);
  IndexTuple<DPeakArray<2, Peak2D> > ind2(2,3,feat2);
  cons2.insert(ind2);
  
  ConsensusPeak<> cons3(cons1,cons2);
  DRange<2> pos_range(1,2,2,3);
  DRange<1> int_range(200,200);
    
  TEST_REAL_EQUAL(cons3.getPosition()[0],1.5)
  TEST_REAL_EQUAL(cons3.getPosition()[1],2.5)
  TEST_REAL_EQUAL(cons3.getIntensity(),200)
  TEST_EQUAL(cons3.getPositionRange() == pos_range, true)
  TEST_EQUAL(cons3.getIntensityRange() == int_range, true)
  ConsensusPeak<>::Group::const_iterator it = cons3.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),1)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
  ++it;
  TEST_REAL_EQUAL(it->getMapIndex(),2)
 	TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(ConsensusPeak(const ConsensusPeak& source))
  DPosition<2> pos(1,2);
  ConsensusPeak<> cons(pos,200);
  DPeak<2> feat;
  feat.setPosition(pos);
  feat.setIntensity(200);
  IndexTuple<DPeakArray<2, Peak2D> > ind(1,3,feat);
  cons.insert(ind);
  ConsensusPeak<> cons_copy(cons);
  
  TEST_REAL_EQUAL(cons_copy.getPosition()[0],1)
  TEST_REAL_EQUAL(cons_copy.getPosition()[1],2)
  TEST_REAL_EQUAL(cons_copy.getIntensity(),200)
  TEST_EQUAL(cons_copy.getPositionRange() == cons.getPositionRange(), true)
  TEST_EQUAL(cons_copy.getIntensityRange() == cons.getIntensityRange(), true)
  TEST_REAL_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_REAL_EQUAL((cons_copy.begin())->getElementIndex(),3)
  TEST_REAL_EQUAL((cons_copy.begin())->getElement().getIntensity(),200)
RESULT

CHECK((ConsensusPeak(const PositionType& pos, const IntensityType& i)))
  DPosition<2> pos(1,2);
  ConsensusPeak<> cons(pos,200);
    
  DRange<2> pos_range;
  DRange<1> int_range;
  TEST_REAL_EQUAL(cons.getPosition()[0],1)
  TEST_REAL_EQUAL(cons.getPosition()[1],2)
  TEST_REAL_EQUAL(cons.getIntensity(),200)
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
  TEST_EQUAL(cons.isEmpty(), true)
RESULT

CHECK((ConsensusPeak(const UnsignedInt& map_1_index, const UnsignedInt& peak_index_1, const ElementType& peak_1, const UnsignedInt& map_2_index, const UnsignedInt& peak_index_2, const ElementType& peak_2)))
  DPosition<2> pos(1,2);
  DPeak<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  
  pos[0]=2;
  pos[1]=3;
  DPeak<2> feat2;
  feat2.setPosition(pos);
  feat2.setIntensity(200);
  
  ConsensusPeak<> cons3(1,3,feat1,2,3,feat2);
  DRange<2> pos_range(1,2,2,3);
  DRange<1> int_range(200,200);
    
  TEST_REAL_EQUAL(cons3.getPosition()[0],1.5)
  TEST_REAL_EQUAL(cons3.getPosition()[1],2.5)
  TEST_REAL_EQUAL(cons3.getIntensity(),200)
  TEST_EQUAL(cons3.getPositionRange() == pos_range, true)
  TEST_EQUAL(cons3.getIntensityRange() == int_range, true)
  ConsensusPeak<>::Group::const_iterator it = cons3.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),1)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
  ++it;
  TEST_REAL_EQUAL(it->getMapIndex(),2)
 	TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK((ConsensusPeak(const UnsignedInt& map_index, const UnsignedInt& peak_index, const ElementType& peak)))
  DPosition<2> pos(1,2);
  DPeak<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  
 	ConsensusPeak<> cons(1,3,feat1);
  DRange<2> pos_range(1,2,1,2);
  DRange<1> int_range(200,200);
    
  TEST_REAL_EQUAL(cons.getPosition()[0],1)
  TEST_REAL_EQUAL(cons.getPosition()[1],2)
  TEST_REAL_EQUAL(cons.getIntensity(),200)
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
  ConsensusPeak<>::Group::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),1)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK((ConsensusPeak(const UnsignedInt& map_index, const UnsignedInt& peak_index, const ElementType& peak, const ConsensusPeak& c_peak)))
  DPosition<2> pos(1,2);
  DPeak<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
    
  pos[0]=2;
  pos[1]=3;
  ConsensusPeak<> cons2(pos,200);
  DPeak<2> feat2;
  feat2.setPosition(pos);
  feat2.setIntensity(200);
  IndexTuple<DPeakArray<2, Peak2D> > ind2(2,3,feat2);
  cons2.insert(ind2);
  
  ConsensusPeak<> cons3(1,3,feat1,cons2);
  DRange<2> pos_range(1,2,2,3);
  DRange<1> int_range(200,200);
    
  TEST_REAL_EQUAL(cons3.getPosition()[0],1.5)
  TEST_REAL_EQUAL(cons3.getPosition()[1],2.5)
  TEST_REAL_EQUAL(cons3.getIntensity(),200)
  TEST_EQUAL(cons3.getPositionRange() == pos_range, true)
  TEST_EQUAL(cons3.getIntensityRange() == int_range, true)
  ConsensusPeak<>::Group::const_iterator it = cons3.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),1)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
  ++it;
  TEST_REAL_EQUAL(it->getMapIndex(),2)
 	TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(Group& getPeaks())
  DPosition<2> pos(1,2);
  DPeak<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  IndexTuple<DPeakArray<2, Peak2D> > ind(2,3,feat1);
  
  ConsensusPeak<>::Group group;
  group.insert(ind);
  
  ConsensusPeak<> cons;
  cons.getPeaks() = group;
    
  ConsensusPeak<>::Group::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(IntensityBoundingBoxType& getIntensityRange())
  DRange<1> int_range(0,200);
  ConsensusPeak<> cons;
  cons.getIntensityRange() = int_range;
  
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
RESULT

CHECK(PositionBoundingBoxType& getPositionRange())
  DRange<2> pos_range(0,1,100,200);
  ConsensusPeak<> cons;
  cons.getPositionRange() = pos_range;
  
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
RESULT

CHECK(const Group& getPeaks() const)
  DPosition<2> pos(1,2);
  DPeak<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  IndexTuple<DPeakArray<2, Peak2D> > ind(2,3,feat1);
  ConsensusPeak<> cons;
  cons.insert(ind);
  const ConsensusPeak<> cons_copy(cons);
  
  ConsensusPeak<>::Group group = cons_copy.getPeaks();
    
  ConsensusPeak<>::Group::const_iterator it = group.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(const IntensityBoundingBoxType& getIntensityRange() const)
  DRange<1> int_range;
  const ConsensusPeak<> cons;
  
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
RESULT

CHECK(const PositionBoundingBoxType& getPositionRange() const)
  DRange<2> pos_range;
  const ConsensusPeak<> cons;
  
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
RESULT

CHECK(void insert(const IndexTuple& tuple))
  DPosition<2> pos(1,2);
  DPeak<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  IndexTuple<DPeakArray<2, Peak2D> > ind(2,3,feat1);
  
  ConsensusPeak<> cons;
  cons.insert(ind);
      
  ConsensusPeak<>::Group::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(void setIntensityRange(const IntensityBoundingBoxType& i))
  DPosition<2> pos(1,2);
  DPeak<2> feat1;
  feat1.setPosition(pos);
  feat1.setIntensity(200);
  IndexTuple<DPeakArray<2, Peak2D> > ind(2,3,feat1);
  
  ConsensusPeak<>::Group group;
  group.insert(ind);
  
  ConsensusPeak<> cons;
  cons.setPeaks(group);
    
  ConsensusPeak<>::Group::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getElement().getIntensity(),200)
RESULT

CHECK(void setPositionRange(const PositionBoundingBoxType& p))
  DRange<1> int_range(0,200);
  ConsensusPeak<> cons;
  cons.setIntensityRange(int_range);
  
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
RESULT

CHECK(void setPeaks(const Group& g))
  DRange<2> pos_range(0,1,100,200);
  ConsensusPeak<> cons;
  cons.setPositionRange(pos_range);
  
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



