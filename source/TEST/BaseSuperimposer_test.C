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
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef LinearMapping TransformationType;
typedef Feature ElementType;
typedef FeatureMap< ElementType> ElementMapType;

class TestSuperimposer : public BaseSuperimposer<ElementMapType>
{
  public:
	TestSuperimposer() 
		: BaseSuperimposer<ElementMapType>()
	{ 
		check_defaults_ = false; 
	}

	virtual void run() 
	{
	}
};

START_TEST(BaseSuperimposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestSuperimposer* ptr = 0;
CHECK((BaseSuperimposer()))
	ptr = new TestSuperimposer();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~BaseSuperimposer()))
	delete ptr;
RESULT

CHECK((const PointMapType& getElementMap(UInt index)))
  
RESULT

CHECK((const PointMapType& getElementMap(UInt index) const))
  ElementMapType map;
  TestSuperimposer bpf;
  bpf.setElementMap(0,map);
  const TestSuperimposer bpf_copy(bpf);
  TEST_EQUAL(&(bpf_copy.getElementMap(0)) == &map,true)
RESULT

CHECK((const TransformationType& getTransformation(UInt dim) const))
  TransformationType trafo(1.,2.);
  TestSuperimposer bsi;
  bsi.setTransformation(0,trafo);
  const TestSuperimposer bsi_copy(bsi);
  
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT

CHECK((virtual void run()=0))
  
RESULT

CHECK((void registerChildren()))
  
RESULT

CHECK((void setElementMap(UInt const index, const PointMapType &element_map)))
  ElementMapType map;
  TestSuperimposer bsi;
  bsi.setElementMap(0,map);
  TEST_EQUAL(&(bsi.getElementMap(0)) == &map,true)
RESULT

CHECK((void setTransformation(UInt dim, const TransformationType& trafo)))
  TransformationType trafo(1.,2.);
  TestSuperimposer bsi;
  bsi.setTransformation(0,trafo);
    
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getSlope(),trafo.getSlope())
  TEST_REAL_EQUAL((bsi.getTransformation(0)).getIntercept(),trafo.getIntercept())
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



