// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Mathias Walzer$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/COMPARISON/PeakAlignment.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <vector>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeakAlignment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeakAlignment* ptr = nullptr;
PeakAlignment* nullPointer = nullptr;
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

  // Test empty spectra - they should return zero
  PeakSpectrum empty_spectrum;
	score = pa(empty_spectrum, s2);
	TEST_REAL_SIMILAR(score, 0.0)

	score = pa(s1, empty_spectrum);
	TEST_REAL_SIMILAR(score, 0.0)
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

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



