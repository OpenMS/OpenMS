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

#include <iostream>

#include <OpenMS/CHEMISTRY/TheoreticalSpectrumGenerator.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////

START_TEST(TheoreticalSpectrumGenerator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

TheoreticalSpectrumGenerator* ptr = 0;

START_SECTION(TheoreticalSpectrumGenerator())
	ptr = new TheoreticalSpectrumGenerator();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator& source))
	TheoreticalSpectrumGenerator copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(~TheoreticalSpectrumGenerator())
	delete ptr;
END_SECTION

ptr = new TheoreticalSpectrumGenerator();
AASequence peptide("IFSQVGK");

START_SECTION(TheoreticalSpectrumGenerator& operator = (const TheoreticalSpectrumGenerator& tsg))
	TheoreticalSpectrumGenerator copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
END_SECTION

START_SECTION(void addPeaks(RichPeakSpectrum& spectrum, const AASequence& peptide, Residue::ResidueType res_type, Int charge = 1))
	RichPeakSpectrum y_spec, b_spec, a_spec;
	ptr->addPeaks(y_spec, peptide, Residue::YIon, 1);
	ptr->addPeaks(b_spec, peptide, Residue::BIon, 1);
	ptr->addPeaks(a_spec, peptide, Residue::AIon, 1);
	TOLERANCE_ABSOLUTE(0.001)
	double y_result[] = {147.113, 204.135, 303.203, 431.262, 518.294, 665.362};
	for (Size i = 0; i != y_spec.size(); ++i)
	{
		TEST_REAL_SIMILAR(y_spec[i].getPosition()[0], y_result[i])
	}
	double b_result[] = {115.1, 261.16, 348.192, 476.251, 575.319, 632.341};
	for (Size i = 0; i != b_spec.size(); ++i)
	{
		TEST_REAL_SIMILAR(b_spec[i].getPosition()[0], b_result[i])
	}

	double a_result[] = {87.1048, 233.165, 320.197, 448.256, 547.324, 604.346};
	for (Size i = 0; i != a_spec.size(); ++i)
	{
		TEST_REAL_SIMILAR(a_spec[i].getPosition()[0], a_result[i])
	}

	RichPeakSpectrum y_spec2;
	ptr->addPeaks(y_spec2, peptide, Residue::YIon, 2);
	TOLERANCE_ABSOLUTE(0.01)
	for (Size i = 0; i != y_spec2.size(); ++i)
	{
		TEST_REAL_SIMILAR(y_spec2[i].getPosition()[0], (y_result[i]+1.0)/2.0)
	}
END_SECTION

START_SECTION(void addPrecursorPeaks(RichPeakSpectrum& spec, const AASequence& peptide, Int charge = 1))
	RichPeakSpectrum spec;
	ptr->addPrecursorPeaks(spec, peptide, 1);
	double result[] = {778.4457, 760.4352, 761.4192};
	for (Size i = 0; i != spec.size(); ++i)
	{
		TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
	}

	RichPeakSpectrum spec2;
	ptr->addPrecursorPeaks(spec2, peptide, 2);
	double result2[] = {389.7265, 380.7212, 381.2132};
	for (Size i = 0; i != spec2.size(); ++i)
	{
		TEST_REAL_SIMILAR(spec2[i].getPosition()[0], result2[i])
	}
	
END_SECTION

START_SECTION(void getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, Int charge = 1))
	RichPeakSpectrum spec;
	ptr->getSpectrum(spec, peptide, 1);
	TEST_EQUAL(spec.size(), 12)

	TOLERANCE_ABSOLUTE(0.001)

	double result[] = {115.1, 147.113, 204.135, 261.16, 303.203, 348.192, 431.262, 476.251, 518.294, 575.319, 632.341, 665.362};
	for (Size i = 0; i != spec.size(); ++i)
	{
		TEST_REAL_SIMILAR(spec[i].getPosition()[0], result[i])
	}

	spec.clear();
	ptr->getSpectrum(spec, peptide, 2);
	TEST_EQUAL(spec.size(), 24)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
