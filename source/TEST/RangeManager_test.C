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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/DPeak.h>
#include <OpenMS/DATASTRUCTURES/RangeManager.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

class RM
	: public RangeManager<2>
{
	public:
		RM()
			: RangeManager<2>()
		{
			
		}
	
		RM(const RM& rhs)
			: RangeManager<2>(rhs)
		{
			
		}
	
		RM& operator = (const RM& rhs)
		{
			if (this==&rhs) return *this;
			
			RangeManager<2>::operator=(rhs);
			
			return *this;
		}
	
		bool operator == (const RM& rhs) const
		{
			return
				RangeManager<2>::operator==(rhs);
				;				
		}
		
		bool operator != (const RM& rhs) const
		{
			return !(operator==(rhs));
		}
	
		virtual void updateRanges()
		{
			std::vector<DPeak<2> > vec;
			DPeak<2> tmp;
			
			tmp.getPosition()[0] = 2.0;
			tmp.getPosition()[1] = 500.0;
			tmp.setIntensity(1.0);
			vec.push_back(tmp);
			
			tmp.getPosition()[0] = 100.0;
			tmp.getPosition()[1] = 1300.0;
			tmp.setIntensity(47110.0);
			vec.push_back(tmp);

			tmp.getPosition()[0] = 2.0;
			tmp.getPosition()[1] = 500.0;
			tmp.setIntensity(1.0);
			vec.push_back(tmp);
			
			clearRanges();
			updateRanges_(vec.begin(), vec.end());
		}

		virtual void updateRanges2()
		{
			std::vector<DPeak<2> > vec;
			DPeak<2> tmp;
			
			tmp.getPosition()[0] = 2.0;
			tmp.getPosition()[1] = 500.0;
			tmp.setIntensity(1.0);
			vec.push_back(tmp);
			
			clearRanges();
			updateRanges_(vec.begin(), vec.end());
		}
		
}; // class RM

START_TEST(RangeManager, "RangeManager")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RM* ptr;
CHECK(RangeManager())
	ptr = new RM();
RESULT

CHECK(~RangeManager())
	delete ptr;
RESULT

CHECK(const PositionType& getMin() const)
	TEST_EQUAL(RM().getMin(), RM::PositionType::max)
RESULT

CHECK(const PositionType& getMax() const)
	TEST_EQUAL(RM().getMax(), RM::PositionType::min_negative)
RESULT

CHECK(const IntensityType getMinInt() const)
	TEST_REAL_EQUAL(RM().getMinInt(), numeric_limits<RM::IntensityType>::max())
RESULT

CHECK(const IntensityType getMaxInt() const)
	TEST_REAL_EQUAL(RM().getMaxInt(), -numeric_limits<RM::IntensityType>::max())
RESULT

CHECK(void updateRanges())
	RM rm;
	
	rm.updateRanges();
	rm.updateRanges(); //second time to check the initialization
	
	TEST_REAL_EQUAL(rm.getMin()[0], 2.0)
	TEST_REAL_EQUAL(rm.getMin()[1], 500.0)
	TEST_REAL_EQUAL(rm.getMax()[0], 100.0)
	TEST_REAL_EQUAL(rm.getMax()[1], 1300.0)
	TEST_REAL_EQUAL(rm.getMinInt(), 1.0)
	TEST_REAL_EQUAL(rm.getMaxInt(), 47110.0)	

	//test with only one point
	rm.updateRanges2(); //second time to check the initialization
	
	TEST_REAL_EQUAL(rm.getMin()[0], 2.0)
	TEST_REAL_EQUAL(rm.getMin()[1], 500.0)
	TEST_REAL_EQUAL(rm.getMax()[0], 2.0)
	TEST_REAL_EQUAL(rm.getMax()[1], 500.0)
	TEST_REAL_EQUAL(rm.getMinInt(), 1.0)
	TEST_REAL_EQUAL(rm.getMaxInt(), 1.0)	
RESULT

CHECK(void clearRanges())
	RM rm;
	rm.updateRanges();
	TEST_REAL_EQUAL(rm.getMin()[0], 2.0)
	TEST_REAL_EQUAL(rm.getMin()[1], 500.0)
	TEST_REAL_EQUAL(rm.getMax()[0], 100.0)
	TEST_REAL_EQUAL(rm.getMax()[1], 1300.0)
	TEST_REAL_EQUAL(rm.getMinInt(), 1.0)
	TEST_REAL_EQUAL(rm.getMaxInt(), 47110.0)
	
	rm.clearRanges();
	TEST_EQUAL(RM().getMin(), RM::PositionType::max)
	TEST_EQUAL(RM().getMax(), RM::PositionType::min_negative)
	TEST_REAL_EQUAL(RM().getMinInt(), numeric_limits<RM::IntensityType>::max())
	TEST_REAL_EQUAL(RM().getMaxInt(), -numeric_limits<RM::IntensityType>::max())
RESULT

CHECK(RangeManager(const RangeManager& rhs))
	RM rm0;
	rm0.updateRanges();
	RM rm(rm0);
	TEST_REAL_EQUAL(rm.getMin()[0], 2.0)
	TEST_REAL_EQUAL(rm.getMin()[1], 500.0)
	TEST_REAL_EQUAL(rm.getMax()[0], 100.0)
	TEST_REAL_EQUAL(rm.getMax()[1], 1300.0)
	TEST_REAL_EQUAL(rm.getMinInt(), 1.0)
	TEST_REAL_EQUAL(rm.getMaxInt(), 47110.0)		
RESULT

CHECK(RangeManager& operator = (const RangeManager& rhs))
	RM rm0;
	rm0.updateRanges();
	RM rm;
	rm = rm0;
	TEST_REAL_EQUAL(rm.getMin()[0], 2.0)
	TEST_REAL_EQUAL(rm.getMin()[1], 500.0)
	TEST_REAL_EQUAL(rm.getMax()[0], 100.0)
	TEST_REAL_EQUAL(rm.getMax()[1], 1300.0)
	TEST_REAL_EQUAL(rm.getMinInt(), 1.0)
	TEST_REAL_EQUAL(rm.getMaxInt(), 47110.0)	
RESULT

CHECK(bool operator == (const RangeManager& rhs) const)
	RM rm0 , rm;
	TEST_EQUAL(rm==rm0, true);
	rm0.updateRanges();
	TEST_EQUAL(rm==rm0, false);
RESULT

CHECK(bool operator != (const RangeManager& rhs) const)
	RM rm0 , rm;
	TEST_EQUAL(rm!=rm0, false);
	rm0.updateRanges();
	TEST_EQUAL(rm!=rm0, true);
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
