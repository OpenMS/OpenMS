// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/DENOVO/CompNovoIonScoring.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(CompNovoIonScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

CompNovoIonScoring* ptr = 0;
START_SECTION(CompNovoIonScoring())
{
	ptr = new CompNovoIonScoring();
	TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((CompNovoIonScoring(const CompNovoIonScoring &source)))
{
  // TODO
}
END_SECTION

START_SECTION((virtual ~CompNovoIonScoring()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void scoreSpectra(Map< DoubleReal, IonScore > &CID_ion_scores, PeakSpectrum &CID_spec, PeakSpectrum &ETD_spec, DoubleReal precursor_weight, Size charge)))
{
  // TODO
}
END_SECTION

START_SECTION((CompNovoIonScoring& operator=(const CompNovoIonScoring &source)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



