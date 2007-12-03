// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/ANALYSIS/ID/PILISModel.h>
#include <OpenMS/FORMAT/DTAFile.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////

START_TEST(PILISModel_test.C, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

PILISModel* ptr = 0;
const AASequence peptide("DFPIANGER");
CHECK(PILISModel())
	ptr = new PILISModel();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~PILISModel())
	delete ptr;
RESULT

ptr = new PILISModel();

CHECK(PILISModel(const PILISModel& model))
	PILISModel copy(*ptr);
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
RESULT

CHECK(PILISModel& operator = (const PILISModel& mode))
	PILISModel copy;
	copy = *ptr;
	TEST_EQUAL(copy.getParameters(), ptr->getParameters())
RESULT

CHECK(void writeGraphMLFile(const String& filename))
RESULT

CHECK(void readFromFile(const String& filename))
	ptr->readFromFile("PILIS/PILIS_default_model.dat");
RESULT

CHECK(void writeToFile(const String& filename))
	String temp_filename("data/PILISModel_model.dat");
	NEW_TMP_FILE(temp_filename)
	ptr->writeToFile(temp_filename);
RESULT

CHECK(void getSpectrum(PeakSpectrum& spec, const AASequence& peptide, UInt charge))
	PeakSpectrum spec;
	ptr->getSpectrum(spec, peptide, 1);
	TEST_EQUAL(spec.size(), 90)
RESULT

CHECK(void train(const PeakSpectrum&, const AASequence& peptide, UInt charge))
	PeakSpectrum spec;
	DTAFile().load("data/PILISSequenceDB_DFPIANGER_1.dta", spec);
	ptr->train(spec, peptide, 1);
RESULT

CHECK(void evaluate())
	ptr->evaluate();
RESULT
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

END_TEST
