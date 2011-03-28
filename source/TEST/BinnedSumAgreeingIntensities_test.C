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
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/SPECTRA/BinnedSumAgreeingIntensities.h>
#include <OpenMS/COMPARISON/SPECTRA/BinnedSpectrum.h>
#include <OpenMS/FORMAT/DTAFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(BinnedSumAgreeingIntensities, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

BinnedSumAgreeingIntensities* ptr = 0;
BinnedSumAgreeingIntensities* nullPointer = 0;
START_SECTION(BinnedSumAgreeingIntensities())
{
	ptr = new BinnedSumAgreeingIntensities();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~BinnedSumAgreeingIntensities())
{
	delete ptr;
}
END_SECTION

ptr = new BinnedSumAgreeingIntensities();

START_SECTION((BinnedSumAgreeingIntensities(const BinnedSumAgreeingIntensities &source)))
{
	BinnedSumAgreeingIntensities copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((BinnedSumAgreeingIntensities& operator=(const BinnedSumAgreeingIntensities &source)))
{
	BinnedSumAgreeingIntensities copy;
	copy = *ptr;
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((double operator()(const BinnedSpectrum &spec1, const BinnedSpectrum &spec2) const))
{
  PeakSpectrum s1, s2;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s2);
  s2.pop_back();
  BinnedSpectrum bs1 (1.5,2,s1);
  BinnedSpectrum bs2 (1.5,2,s2);  

  double score = (*ptr)(bs1, bs2);
  TEST_REAL_SIMILAR(score, 0.997576)
}
END_SECTION

START_SECTION((double operator()(const BinnedSpectrum &spec) const ))
{
  PeakSpectrum s1;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
  BinnedSpectrum bs1 (1.5,2,s1);
  double score = (*ptr)(bs1);
  TEST_REAL_SIMILAR(score, 1);
}
END_SECTION

START_SECTION((static BinnedSpectrumCompareFunctor* create()))
{
	BinnedSpectrumCompareFunctor* bsf = BinnedSumAgreeingIntensities::create();
	BinnedSumAgreeingIntensities bsp;
	TEST_EQUAL(bsf->getParameters(), bsp.getParameters())
	TEST_EQUAL(bsf->getName(), bsp.getName())
}
END_SECTION

START_SECTION((static const String getProductName()))
{
	TEST_EQUAL(ptr->getProductName(), "BinnedSumAgreeingIntensities")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



