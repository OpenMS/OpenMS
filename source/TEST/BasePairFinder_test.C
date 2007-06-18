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
#include <OpenMS/ANALYSIS/MAPMATCHING/BasePairFinder.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef LinearMapping TransformationType;
typedef Feature ElementType;
typedef FeatureMap< ElementType> ElementMapType;
typedef ElementPair < Feature > ElementPairType;
typedef vector< ElementPairType >  ElementPairVectorType;

class TestPairFinder 
	: public BasePairFinder<ElementMapType>
{
  public:
	TestPairFinder() 
		: BasePairFinder<ElementMapType>()
	{
		check_defaults_ = false; 
	}
	
	virtual void findElementPairs() 
	{
		
	}

};

START_TEST(BasePairFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestPairFinder* ptr = 0;
CHECK((BasePairFinder()))
	ptr = new TestPairFinder();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~BasePairFinder()))
	delete ptr;
RESULT


CHECK((const ElementPairVectorType& getElementPairs() const))
	ElementPairVectorType pairs;
	TestPairFinder bpf;
	bpf.setElementPairs(pairs);
  const TestPairFinder bpf_copy(bpf);
  
  TEST_EQUAL(&(bpf_copy.getElementPairs()) == &pairs,true)
RESULT

CHECK((const PointMapType& getElementMap(UInt index) const))
  ElementMapType map;
  TestPairFinder bpf;
  bpf.setElementMap(0,map);
  const TestPairFinder bpf_copy(bpf);
  TEST_EQUAL(&(bpf_copy.getElementMap(0)) == &map,true)
RESULT

CHECK((static void registerChildren()))
  
RESULT

CHECK((void findElementPairs()))
  
RESULT

CHECK((void setElementMap(UInt const index, const PointMapType& element_map)))
  ElementMapType map;
  TestPairFinder bpf;
  bpf.setElementMap(0,map);
  TEST_EQUAL(&(bpf.getElementMap(0)) == &map,true)
RESULT

CHECK((void setElementPairs(ElementPairVectorType& element_pairs)))
  ElementPairVectorType pairs;
  TestPairFinder bpf;
  bpf.setElementPairs(pairs);
  const TestPairFinder bpf_copy(bpf);
  TEST_EQUAL(&(bpf_copy.getElementPairs()) == &pairs,true)
RESULT

CHECK((void setTransformation(UInt dim, const TransformationType& trafo)))
  TransformationType trafo(1.,2.);
  TestPairFinder bpf;
  bpf.setTransformation(0,trafo);
    
  TEST_REAL_EQUAL((bpf.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bpf.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT

CHECK((const TransformationType& getTransformation(UInt dim) const))
  TransformationType trafo(1.,2.);
  TestPairFinder bpf;
  bpf.setTransformation(0,trafo);
  const TestPairFinder bpf_copy(bpf);
  
  TEST_REAL_EQUAL((bpf.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bpf.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT


CHECK((int dumpElementPairs(const String &filename)))

RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



