// -*- mode: C++; tab-width: 2; -*-
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
// $Maintainer: Andreas Bertsch $
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

CHECK(TheoreticalSpectrumGenerator())
	ptr = new TheoreticalSpectrumGenerator();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(TheoreticalSpectrumGenerator(const TheoreticalSpectrumGenerator& source))
	TheoreticalSpectrumGenerator copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
RESULT

CHECK(~TheoreticalSpectrumGenerator())
	delete ptr;
RESULT

ptr = new TheoreticalSpectrumGenerator();
AASequence peptide("IFSQVGK");

CHECK(TheoreticalSpectrumGenerator& operator = (const TheoreticalSpectrumGenerator& tsg))
	TheoreticalSpectrumGenerator copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
RESULT

CHECK(void addPeaks(RichPeakSpectrum& spectrum, const AASequence& peptide, Residue::ResidueType res_type, Int charge = 1))
	RichPeakSpectrum y_spec, b_spec, a_spec;
	ptr->addPeaks(y_spec, peptide, Residue::YIon, 1);
	ptr->addPeaks(b_spec, peptide, Residue::BIon, 1);
	ptr->addPeaks(a_spec, peptide, Residue::AIon, 1);
	PRECISION(0.001)
	double y_result[] = {147.113, 204.135, 303.203, 431.262, 518.294, 665.362};
	for (unsigned int i = 0; i != y_spec.size(); ++i)
	{
		TEST_REAL_EQUAL(y_spec.getContainer()[i].getPosition()[0], y_result[i])
	}
	double b_result[] = {115.1, 261.16, 348.192, 476.251, 575.319, 632.341};
	for (unsigned int i = 0; i != b_spec.size(); ++i)
	{
		TEST_REAL_EQUAL(b_spec.getContainer()[i].getPosition()[0], b_result[i])
	}

	double a_result[] = {87.1048, 233.165, 320.197, 448.256, 547.324, 604.346};
	for (unsigned int i = 0; i != a_spec.size(); ++i)
	{
		TEST_REAL_EQUAL(a_spec.getContainer()[i].getPosition()[0], a_result[i])
	}

	RichPeakSpectrum y_spec2;
	ptr->addPeaks(y_spec2, peptide, Residue::YIon, 2);
	PRECISION(0.01)
	for (unsigned int i = 0; i != y_spec2.size(); ++i)
	{
		TEST_REAL_EQUAL(y_spec2.getContainer()[i].getPosition()[0], (y_result[i]+1.0)/2.0)
	}
RESULT

CHECK(void addPrecursorPeaks(RichPeakSpectrum& spec, const AASequence& peptide, Int charge = 1))
	RichPeakSpectrum spec;
	ptr->addPrecursorPeaks(spec, peptide, 1);
	double result[] = {778.916, 760.901, 761.885};
	for (unsigned int i = 0; i != spec.size(); ++i)
	{
		TEST_REAL_EQUAL(spec.getContainer()[i].getPosition()[0], result[i])
	}

	RichPeakSpectrum spec2;
	ptr->addPrecursorPeaks(spec2, peptide, 2);
	double result2[] = {389.962, 380.954, 381.447};
	for (unsigned int i = 0; i != spec2.size(); ++i)
	{
		TEST_REAL_EQUAL(spec2.getContainer()[i].getPosition()[0], result2[i])
	}
	
RESULT

CHECK(void getSpectrum(RichPeakSpectrum& spec, const AASequence& peptide, Int charge = 1))
	RichPeakSpectrum spec;
	ptr->getSpectrum(spec, peptide, 1);
	TEST_EQUAL(spec.size(), 12)

	PRECISION(0.001)

	double result[] = {115.1, 147.113, 204.135, 261.16, 303.203, 348.192, 431.262, 476.251, 518.294, 575.319, 632.341, 665.362};
	for (unsigned int i = 0; i != spec.size(); ++i)
	{
		TEST_REAL_EQUAL(spec.getContainer()[i].getPosition()[0], result[i])
	}

	spec.clear();
	ptr->getSpectrum(spec, peptide, 2);
	TEST_EQUAL(spec.size(), 24)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
