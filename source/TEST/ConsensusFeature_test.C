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

///////////////////////////
#include <OpenMS/KERNEL/ConsensusFeature.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusFeature, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusFeature* ptr = 0;
CHECK(ConsensusFeature())
	ptr = new ConsensusFeature();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~ConsensusFeature())
	delete ptr;
RESULT

Feature tmp_feature;
tmp_feature.setRT(1);
tmp_feature.setMZ(2);
tmp_feature.setIntensity(200);

CHECK(ConsensusFeature& operator=(const ConsensusFeature& source))
  ConsensusFeature cons(tmp_feature);
  cons.insert(1,3,tmp_feature);
  
  ConsensusFeature cons_copy;
  cons_copy = cons;
  
  TEST_REAL_EQUAL(cons_copy.getRT(),1)
  TEST_REAL_EQUAL(cons_copy.getMZ(),2)
  TEST_REAL_EQUAL(cons_copy.getIntensity(),200)
  TEST_EQUAL(cons_copy.getPositionRange() == cons.getPositionRange(), true)
  TEST_EQUAL(cons_copy.getIntensityRange() == cons.getIntensityRange(), true)
  TEST_REAL_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_REAL_EQUAL((cons_copy.begin())->getElementIndex(),3)
  TEST_REAL_EQUAL((cons_copy.begin())->getIntensity(),200)
RESULT

CHECK(ConsensusFeature(const ConsensusFeature& source))
  
  ConsensusFeature cons(tmp_feature);
  cons.insert(1,3,tmp_feature);
  ConsensusFeature cons_copy(cons);
  
  TEST_REAL_EQUAL(cons_copy.getRT(),1)
  TEST_REAL_EQUAL(cons_copy.getMZ(),2)
  TEST_REAL_EQUAL(cons_copy.getIntensity(),200)
  TEST_EQUAL(cons_copy.getPositionRange() == cons.getPositionRange(), true)
  TEST_EQUAL(cons_copy.getIntensityRange() == cons.getIntensityRange(), true)
  TEST_REAL_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_REAL_EQUAL((cons_copy.begin())->getElementIndex(),3)
  TEST_REAL_EQUAL((cons_copy.begin())->getIntensity(),200)
RESULT

CHECK((ConsensusFeature(const RawDataPoint2D& pos)))
  
  ConsensusFeature cons(tmp_feature);
  TEST_REAL_EQUAL(cons.getRT(),1)
  TEST_REAL_EQUAL(cons.getMZ(),2)
  TEST_REAL_EQUAL(cons.getIntensity(),200)
  TEST_EQUAL(cons.empty(), true)
RESULT

CHECK((ConsensusFeature(UInt map_index, UInt feature_index, const ElementType& feature)))
 	ConsensusFeature cons(1,3,tmp_feature);
  DRange<2> pos_range(1,2,1,2);
  DRange<1> int_range(200,200);
    
  TEST_REAL_EQUAL(cons.getRT(),1)
  TEST_REAL_EQUAL(cons.getMZ(),2)
  TEST_REAL_EQUAL(cons.getIntensity(),200)
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
  ConsensusFeature::Group::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),1)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getIntensity(),200)
RESULT

CHECK(IntensityBoundingBoxType& getIntensityRange())
  DRange<1> int_range(0,200);
  ConsensusFeature cons;
  cons.getIntensityRange() = int_range;
  
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
RESULT

CHECK(PositionBoundingBoxType& getPositionRange())
  DRange<2> pos_range(0,1,100,200);
  ConsensusFeature cons;
  cons.getPositionRange() = pos_range;
  
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
RESULT

CHECK(const Group& getFeatures() const)
  ConsensusFeature cons;
  cons.insert(2,3,tmp_feature);
  const ConsensusFeature cons_copy(cons);
  
  ConsensusFeature::Group group = cons_copy.getFeatures();
    
  ConsensusFeature::Group::const_iterator it = group.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getIntensity(),200)
RESULT

CHECK(const IntensityBoundingBoxType& getIntensityRange() const)
  DRange<1> int_range;
  const ConsensusFeature cons;
  
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
RESULT

CHECK(const PositionBoundingBoxType& getPositionRange() const)
  DRange<2> pos_range;
  const ConsensusFeature cons;
  
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
RESULT

CHECK(void insert(const IndexTuple& tuple))
  ConsensusFeature cons;
  cons.insert(2,3,tmp_feature);
      
  ConsensusFeature::Group::const_iterator it = cons.begin();
  TEST_REAL_EQUAL(it->getMapIndex(),2)
  TEST_REAL_EQUAL(it->getElementIndex(),3)
  TEST_REAL_EQUAL(it->getIntensity(),200)
RESULT

CHECK(void setIntensityRange(const IntensityBoundingBoxType& i))
  DRange<1> int_range(0,200);
  ConsensusFeature cons;
  cons.setIntensityRange(int_range);
  
  TEST_EQUAL(cons.getIntensityRange() == int_range, true)
RESULT

CHECK(void setPositionRange(const PositionBoundingBoxType& p))
  DRange<2> pos_range(0,1,100,200);
  ConsensusFeature cons;
  cons.setPositionRange(pos_range);
  
  TEST_EQUAL(cons.getPositionRange() == pos_range, true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



