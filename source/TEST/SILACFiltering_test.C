// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Lars Nilse $
// $Authors: Lars Nilse, Holger Plattfaut $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFiltering.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SILACFiltering, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SILACFiltering* ptr = 0;
START_SECTION(SILACFiltering())
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(~SILACFiltering())
{
	delete ptr;
}
END_SECTION

START_SECTION((SILACFiltering(MSExperiment< Peak1D > &exp, const DoubleReal mz_stepwidth, const DoubleReal intensity_cutoff, const DoubleReal intensity_correlation, const bool allow_missing_peaks)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void addFilter(SILACFilter &filter)))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((void filterDataPoints()))
{
  NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

