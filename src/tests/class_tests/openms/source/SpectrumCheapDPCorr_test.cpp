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

	// A wide 'variation' gives many peaks of these two spectra several possible partners, so
	// operator() hands those blocks to dynprog_(), which has to honor keeppeaks as well.
	auto count_ambiguous = [&](bool poison, bool keeppeaks) -> Size
	{
		alignas(SpectrumCheapDPCorr) unsigned char raw[sizeof(SpectrumCheapDPCorr)];
		// Construct the object with placement new in storage filled with 0xFF (poison) or 0x00
		// (reference). In the poisoned run, a member the constructor forgets to initialize holds 0xFF
		// ('true' for a bool) instead of the zero that fresh memory usually happens to contain, so a
		// missing initialization changes the result instead of passing by luck. keeppeaks_ used to be
		// such a member (#10237): poisoned, dynprog_() kept 72 instead of 55 peaks.
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
	// the result must not depend on what was in memory before construction ...
	TEST_EQUAL(count_ambiguous(false, false), count_ambiguous(true, false))
	TEST_EQUAL(count_ambiguous(false, true), count_ambiguous(true, true))
	// ... and dynprog_() must honor keeppeaks: a keeppeaks_ that is initialized to false but never
	// set from the parameter passes the checks above, yet drops dynprog_()'s unaligned peaks even
	// with keeppeaks=1
	TEST_EQUAL(count_ambiguous(false, false), 55)
	TEST_EQUAL(count_ambiguous(false, true), 146)

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
