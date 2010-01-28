// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2010 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/RangeManager.h>

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
			std::vector<Peak2D > vec;
			Peak2D tmp;

			tmp.getPosition()[0] = 2.0;
			tmp.getPosition()[1] = 500.0;
			tmp.setIntensity(1.0f);
			vec.push_back(tmp);

			tmp.getPosition()[0] = 100.0;
			tmp.getPosition()[1] = 1300.0;
			tmp.setIntensity(47110.0);
			vec.push_back(tmp);

			tmp.getPosition()[0] = 2.0;
			tmp.getPosition()[1] = 500.0;
			tmp.setIntensity(1.0f);
			vec.push_back(tmp);

			clearRanges();
			updateRanges_(vec.begin(), vec.end());
		}

		virtual void updateRanges2()
		{
			std::vector<Peak2D > vec;
			Peak2D tmp;

			tmp.getPosition()[0] = 2.0;
			tmp.getPosition()[1] = 500.0;
			tmp.setIntensity(1.0f);
			vec.push_back(tmp);

			clearRanges();
			updateRanges_(vec.begin(), vec.end());
		}

}; // class RM

START_TEST(RangeManager, "RangeManager")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RM* ptr;
START_SECTION((RangeManager()))
	ptr = new RM();
  TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~RangeManager()))
	delete ptr;
END_SECTION

START_SECTION((const PositionType& getMin() const))
	TEST_EQUAL(RM().getMin(), RM::PositionType::maxPositive())
END_SECTION

START_SECTION((const PositionType& getMax() const))
	TEST_EQUAL(RM().getMax(), RM::PositionType::minNegative())
END_SECTION

START_SECTION((DoubleReal getMinInt() const ))
	TEST_REAL_SIMILAR(RM().getMinInt(), numeric_limits<DoubleReal>::max())
END_SECTION

START_SECTION((DoubleReal getMaxInt() const ))
	TEST_REAL_SIMILAR(RM().getMaxInt(), -numeric_limits<DoubleReal>::max())
END_SECTION

START_SECTION((virtual void updateRanges()=0))
	RM rm;

	rm.updateRanges();
	rm.updateRanges(); //second time to check the initialization

	TEST_REAL_SIMILAR(rm.getMin()[0], 2.0)
	TEST_REAL_SIMILAR(rm.getMin()[1], 500.0)
	TEST_REAL_SIMILAR(rm.getMax()[0], 100.0)
	TEST_REAL_SIMILAR(rm.getMax()[1], 1300.0)
	TEST_REAL_SIMILAR(rm.getMinInt(), 1.0)
	TEST_REAL_SIMILAR(rm.getMaxInt(), 47110.0)

	//test with only one point
	rm.updateRanges2(); //second time to check the initialization

	TEST_REAL_SIMILAR(rm.getMin()[0], 2.0)
	TEST_REAL_SIMILAR(rm.getMin()[1], 500.0)
	TEST_REAL_SIMILAR(rm.getMax()[0], 2.0)
	TEST_REAL_SIMILAR(rm.getMax()[1], 500.0)
	TEST_REAL_SIMILAR(rm.getMinInt(), 1.0)
	TEST_REAL_SIMILAR(rm.getMaxInt(), 1.0)
END_SECTION

START_SECTION((void clearRanges()))
	RM rm;
	rm.updateRanges();
	TEST_REAL_SIMILAR(rm.getMin()[0], 2.0)
	TEST_REAL_SIMILAR(rm.getMin()[1], 500.0)
	TEST_REAL_SIMILAR(rm.getMax()[0], 100.0)
	TEST_REAL_SIMILAR(rm.getMax()[1], 1300.0)
	TEST_REAL_SIMILAR(rm.getMinInt(), 1.0)
	TEST_REAL_SIMILAR(rm.getMaxInt(), 47110.0)

	rm.clearRanges();
	TEST_EQUAL(RM().getMin(), RM::PositionType::maxPositive())
	TEST_EQUAL(RM().getMax(), RM::PositionType::minNegative())
	TEST_REAL_SIMILAR(RM().getMinInt(), numeric_limits<DoubleReal>::max())
	TEST_REAL_SIMILAR(RM().getMaxInt(), -numeric_limits<DoubleReal>::max())
END_SECTION

START_SECTION((RangeManager(const RangeManager& rhs)))
	RM rm0;
	rm0.updateRanges();
	RM rm(rm0);
	TEST_REAL_SIMILAR(rm.getMin()[0], 2.0)
	TEST_REAL_SIMILAR(rm.getMin()[1], 500.0)
	TEST_REAL_SIMILAR(rm.getMax()[0], 100.0)
	TEST_REAL_SIMILAR(rm.getMax()[1], 1300.0)
	TEST_REAL_SIMILAR(rm.getMinInt(), 1.0)
	TEST_REAL_SIMILAR(rm.getMaxInt(), 47110.0)
END_SECTION

START_SECTION((RangeManager& operator = (const RangeManager& rhs)))
	RM rm0;
	rm0.updateRanges();
	RM rm;
	rm = rm0;
	TEST_REAL_SIMILAR(rm.getMin()[0], 2.0)
	TEST_REAL_SIMILAR(rm.getMin()[1], 500.0)
	TEST_REAL_SIMILAR(rm.getMax()[0], 100.0)
	TEST_REAL_SIMILAR(rm.getMax()[1], 1300.0)
	TEST_REAL_SIMILAR(rm.getMinInt(), 1.0)
	TEST_REAL_SIMILAR(rm.getMaxInt(), 47110.0)
END_SECTION

START_SECTION((bool operator == (const RangeManager& rhs) const))
	RM rm0 , rm;
	TEST_EQUAL(rm==rm0, true);
	rm0.updateRanges();
	TEST_EQUAL(rm==rm0, false);
END_SECTION

START_SECTION((bool operator != (const RangeManager& rhs) const))
	RM rm0 , rm;
	TEST_EQUAL(rm!=rm0, false);
	rm0.updateRanges();
	TEST_EQUAL(rm!=rm0, true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
