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
#include <OpenMS/FILTERING/DATAREDUCTION/SILACFilter.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SILACFilter, "$Id$")

std::vector<DoubleReal> mass_separations;
mass_separations.push_back(4);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((SILACFilter(std::vector< DoubleReal > mass_separations, Int charge, DoubleReal model_deviation, Int isotopes_per_peptide, DoubleReal intensity_cutoff, DoubleReal intensity_correlation, bool allow_missing_peaks)))
{
  SILACFilter f(mass_separations, 1, 1, 3, 0, 0, false);
  TEST_EQUAL(f.getCharge(), 1);
}
END_SECTION

START_SECTION((std::vector<DoubleReal> getPeakPositions()))
{
  SILACFilter f(mass_separations, 1, 1, 3, 0, 0, false);
  // XXX: Segfaults
  // TEST_EQUAL(f.getPeakPositions().size(), 0);
}
END_SECTION

START_SECTION((const std::vector<DoubleReal>& getExpectedMzShifts()))
{
  const UInt peaks_per_peptide = 3;
  SILACFilter f(mass_separations, 1, 1, peaks_per_peptide, 0, 0, false);
  TEST_EQUAL(f.getExpectedMzShifts().size(), (mass_separations.size() + 1) * peaks_per_peptide);
}
END_SECTION

START_SECTION((std::vector<SILACPattern>& getElements()))
{
  SILACFilter f(mass_separations, 1, 1, 3, 0, 0, false);
  TEST_EQUAL(f.getElements().size(), 0);
}
END_SECTION

START_SECTION((Int getCharge()))
{
  SILACFilter f(mass_separations, 1, 1, 3, 0, 0, false);
  TEST_EQUAL(f.getCharge(), 1);
}
END_SECTION

START_SECTION((std::vector<DoubleReal>& getMassSeparations()))
{
  SILACFilter f(mass_separations, 1, 1, 3, 0, 0, false);
  TEST_EQUAL(f.getMassSeparations() == mass_separations, true);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

