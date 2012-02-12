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
// $Maintainer: Mathias Walzer$
// $Authors: $
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
PeakAlignment* nullPointer = 0;
START_SECTION(PeakAlignment())
{
	ptr = new PeakAlignment();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~PeakAlignment())
{
	delete ptr;
}
END_SECTION

START_SECTION((PeakAlignment(const PeakAlignment &source)))
{
	ptr = new PeakAlignment();
	PeakAlignment copy(*ptr);
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((PeakAlignment& operator=(const PeakAlignment &source)))
{
	PeakAlignment copy;
	copy = *ptr;
	TEST_EQUAL(copy.getName(), ptr->getName());
	TEST_EQUAL(copy.getParameters(), ptr->getParameters());
}
END_SECTION

START_SECTION((double operator()(const PeakSpectrum &spec1, const PeakSpectrum &spec2) const ))
{
	PeakAlignment pa;
	PeakSpectrum s1, s2;
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s2);
	s2.pop_back();
	double score = pa(s1, s2);
	TEST_REAL_SIMILAR(score, 0.997477)
}
END_SECTION

START_SECTION((double operator()(const PeakSpectrum &spec) const ))
{
  PeakSpectrum s1;
  DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
  double score = (*ptr)(s1);
  TEST_REAL_SIMILAR(score, 1);
}
END_SECTION

START_SECTION((vector< pair<Size,Size> > getAlignmentTraceback(const PeakSpectrum &spec1, const PeakSpectrum &spec2) const ))
{
	PeakAlignment pa;
	PeakSpectrum s1, s2;
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s1);
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("PILISSequenceDB_DFPIANGER_1.dta"), s2);
	vector< pair<Size,Size> > result, tester;
	result = pa.getAlignmentTraceback(s1,s2);
	for (Size i = 0; i < 127; ++i)
	{
		tester.push_back(pair<Size,Size>(i,i));
	}
	TEST_EQUAL(tester.size(),result.size())
	for (Size i = 0; i < tester.size(); ++i)
	{
		TEST_EQUAL(tester.at(i).first,result.at(i).first)
	}
}
END_SECTION

START_SECTION((static PeakSpectrumCompareFunctor* create()))
{
	PeakSpectrumCompareFunctor* psf = PeakAlignment::create();
	PeakAlignment pa;
	TEST_EQUAL(psf->getParameters(), pa.getParameters())
	TEST_EQUAL(psf->getName(), pa.getName())
}
END_SECTION

START_SECTION((static const String getProductName()))
{
	TEST_EQUAL(ptr->getProductName(), "PeakAlignment")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



