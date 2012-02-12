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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar, Holger Franken $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FILTERING/SMOOTHING/LowessSmoothing.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(LowessSmoothing, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

LowessSmoothing* ptr = 0;
LowessSmoothing* null_ptr = 0;
START_SECTION(LowessSmoothing())
{
	ptr = new LowessSmoothing();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~LowessSmoothing())
{
	delete ptr;
}
END_SECTION

/////

std::vector<DoubleReal> x, y, out;

DoubleReal exp1[20] = {10, 20, 30, 40, 50, 60, 70.75057075, 81.91410776, 90.49938042, 93.67188213, 90.49938042, 81.91410776, 70.75057075, 60, 50, 40, 30, 20, 10, 7.638334409e-14};
DoubleReal exp2[20] = {4.940778184, 19.1953138, 32.45871201, 44.62566121, 55.59150285, 65.28588352, 73.78027456, 81.64413917, 87.38364167, 89.36964666, 87.38364167, 81.64413917, 72.86539444, 63.49165214, 53.94643243, 43.76172539, 32.89091229, 21.38760603, 9.323517923, -3.233540303};


for (DoubleReal i = 1.0; i <= 20.0; i += 1.0)
{
    x.push_back(i);
}

for (DoubleReal i = 1.0; i <= 10.0; i += 1.0)
{
    y.push_back(i*10);
}

for (DoubleReal i = 1.0; i <= 10.0; i += 1.0)
{
    y.push_back(100.0 - i*10);
}

y.push_back(10.0);

LowessSmoothing lowsmooth;
Param lowpar;
lowpar.setValue("window_size", 15);

START_SECTION(void smoothData(const DoubleVector&, const DoubleVector&, DoubleVector&))
{
    TEST_EXCEPTION(Exception::InvalidValue, lowsmooth.smoothData(x, y, out));

    y.pop_back();
    out.clear();

    lowsmooth.smoothData(x, y, out);

    for (Size i = 0; i < out.size(); ++i)
    {
        TEST_REAL_SIMILAR(out[i], exp1[i]);
    }

    out.clear();
    lowsmooth.setParameters(lowpar);
    lowsmooth.smoothData(x, y, out);

    for (Size i = 0; i < out.size(); ++i)
    {
        TEST_REAL_SIMILAR(out[i], exp2[i]);
    }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



