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
// $Maintainer: Mathias Walzer$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/COMPARISON/SPECTRA/PeakAlignment.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakAlignment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakAlignment* ptr = 0;
CHECK(PeakAlignment())
{
	ptr = new PeakAlignment();
	TEST_NOT_EQUAL(ptr, 0)
}
RESULT

CHECK(~PeakAlignment())
{
	delete ptr;
}
RESULT

ptr = new PeakAlignment();

CHECK((PeakAlignment(const PeakAlignment &source)))
{
	PeakAlignment copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
RESULT

CHECK((PeakAlignment& operator=(const PeakAlignment &source)))
{
	PeakAlignment copy;
	copy = *ptr;
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
RESULT

CHECK((double operator()(const PeakSpectrum &spec1, const PeakSpectrum &spec2) const ))
{
  PeakSpectrum s1, s2;
  DTAFile().load("data/PILISSequenceDB_DFPIANGER_1.dta", s1);
  DTAFile().load("data/PILISSequenceDB_DFPIANGER_1.dta", s2);
  s2.getContainer().pop_back();
  double score = (*ptr)(s1, s2);
  TEST_REAL_EQUAL(score, 0.997477)
}
RESULT

CHECK((double operator()(const PeakSpectrum &spec) const ))
{
  PeakSpectrum s1;
  DTAFile().load("data/PILISSequenceDB_DFPIANGER_1.dta", s1);
  double score = (*ptr)(s1);
  TEST_REAL_EQUAL(score, 1);
}
RESULT

CHECK((vector< pair<UInt,UInt> > getAlignmentTraceback(const PeakSpectrum &spec1, const PeakSpectrum &spec2) const ))
{
  // TODO
}
RESULT

CHECK((static PeakSpectrumCompareFunctor* create()))
{
	PeakSpectrumCompareFunctor* psf = PeakAlignment::create();
	PeakAlignment pa;
	TEST_EQUAL(psf->getParameters(), pa.getParameters())
	TEST_EQUAL(psf->getName(), pa.getName())
}
RESULT

CHECK((static const String getProductName()))
{
	TEST_EQUAL(ptr->getProductName(), "PeakAlignment")
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



