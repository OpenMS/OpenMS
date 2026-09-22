// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Volker Mosthaf, Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/COMPARISON/SpectrumCheapDPCorr.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/FORMAT/DTAFile.h>

#include <cstring>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(SpectrumCheapDPCorr, "$Id$")

/////////////////////////////////////////////////////////////

SpectrumCheapDPCorr* e_ptr = nullptr;
SpectrumCheapDPCorr* e_nullPointer = nullptr;
START_SECTION(SpectrumCheapDPCorr())
	e_ptr = new SpectrumCheapDPCorr;
	TEST_NOT_EQUAL(e_ptr, e_nullPointer)
END_SECTION

START_SECTION(~SpectrumCheapDPCorr())
	delete e_ptr;
END_SECTION

e_ptr = new SpectrumCheapDPCorr();

START_SECTION(SpectrumCheapDPCorr(const SpectrumCheapDPCorr& source))
	SpectrumCheapDPCorr copy(*e_ptr);
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION(SpectrumCheapDPCorr& operator = (const SpectrumCheapDPCorr& source))
	SpectrumCheapDPCorr copy;
	copy = *e_ptr;
	TEST_EQUAL(copy.getParameters(), e_ptr->getParameters())
	TEST_EQUAL(copy.getName(), e_ptr->getName())
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a, const PeakSpectrum& b) const)
	DTAFile dta_file;
	PeakSpectrum spec1;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec1);

	DTAFile dta_file2;
	PeakSpectrum spec2;
	dta_file2.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests_2.dta"), spec2);

	double score = (*e_ptr)(spec1, spec2);

	TOLERANCE_ABSOLUTE(0.1)
	TEST_REAL_SIMILAR(score, 10145.4)

	score = (*e_ptr)(spec1, spec1);

	TEST_REAL_SIMILAR(score, 12295.5)
	
	SpectrumCheapDPCorr corr;
	score = corr(spec1, spec2);
	TEST_REAL_SIMILAR(score, 10145.4)

	score = corr(spec1, spec1);

	TEST_REAL_SIMILAR(score, 12295.5)
	
END_SECTION

START_SECTION(const PeakSpectrum& lastconsensus() const)
	TEST_EQUAL(e_ptr->lastconsensus().size(), 121)
END_SECTION

START_SECTION((Map<UInt, UInt> getPeakMap() const))
	TEST_EQUAL(e_ptr->getPeakMap().size(), 121)
END_SECTION

START_SECTION([EXTRA] the 'keeppeaks' parameter reaches the consensus spectrum)
{
	DTAFile dta_file;
	PeakSpectrum spec1;
	dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec1);
	PeakSpectrum spec2;
	DTAFile().load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests_2.dta"), spec2);

	// keeppeaks = 0 (default): peaks without an alignment partner are dropped
	SpectrumCheapDPCorr corr;
	corr(spec1, spec2);
	const Size dropped = corr.lastconsensus().size();

	// keeppeaks = 1: they are kept, so the consensus has more peaks
	Param p = corr.getParameters();
	p.setValue("keeppeaks", 1);
	corr.setParameters(p);
	corr(spec1, spec2);
	const Size kept = corr.lastconsensus().size();

	TEST_EQUAL(dropped, 9)
	TEST_EQUAL(kept, 199)

	// dynprog_() -- the O(n^2) dynamic-programming step that resolves an ambiguous, many-to-many
	// block of peaks -- reads keeppeaks_ directly (not through the local above, which only covers
	// operator()'s own single-peak "no partner" branches). A wide 'variation' makes many peaks of
	// these two real spectra mutually ambiguous, so dynprog_() runs repeatedly. To prove the class
	// member is genuinely read (not just consistently zero-initialized by chance on this platform),
	// construct into memory explicitly poisoned to a non-zero pattern beforehand: with the member
	// correctly (re-)assigned in the constructor, the poison must make no difference to the result.
	auto count_ambiguous = [&](bool poison, bool keeppeaks) -> Size
	{
		alignas(SpectrumCheapDPCorr) unsigned char raw[sizeof(SpectrumCheapDPCorr)];
		std::memset(raw, poison ? 0xFF : 0x00, sizeof(raw));
		SpectrumCheapDPCorr* poisoned_corr = new (raw) SpectrumCheapDPCorr();
		Param poisoned_param = poisoned_corr->getParameters();
		poisoned_param.setValue("keeppeaks", keeppeaks ? 1 : 0);
		poisoned_param.setValue("variation", 0.02);
		poisoned_corr->setParameters(poisoned_param);
		(*poisoned_corr)(spec1, spec2);
		const Size result = poisoned_corr->lastconsensus().size();
		poisoned_corr->~SpectrumCheapDPCorr();
		return result;
	};
	TEST_EQUAL(count_ambiguous(false, false), count_ambiguous(true, false))
	TEST_EQUAL(count_ambiguous(false, true), count_ambiguous(true, true))

	// a copy keeps the setting
	SpectrumCheapDPCorr copy(corr);
	copy(spec1, spec2);
	TEST_EQUAL(copy.lastconsensus().size(), kept)

	SpectrumCheapDPCorr assigned;
	assigned = corr;
	assigned(spec1, spec2);
	TEST_EQUAL(assigned.lastconsensus().size(), kept)
}
END_SECTION

START_SECTION(double operator () (const PeakSpectrum& a) const)
  DTAFile dta_file;
  PeakSpectrum spec1;
  dta_file.load(OPENMS_GET_TEST_DATA_PATH("Transformers_tests.dta"), spec1);

  double score = (*e_ptr)(spec1);

  TEST_REAL_SIMILAR(score, 12295.5)

END_SECTION

START_SECTION(void setFactor(double f))
	e_ptr->setFactor(0.3);

	TEST_EXCEPTION(Exception::OutOfRange, e_ptr->setFactor(1.1))
END_SECTION

delete e_ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
