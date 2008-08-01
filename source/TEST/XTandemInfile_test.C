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

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/XTandemInfile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

#include <vector>

///////////////////////////

START_TEST(XTandemInfile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

XTandemInfile xml_file;
XTandemInfile* ptr;
ProteinIdentification protein_identification;
vector<PeptideIdentification> peptide_identifications; 
vector<PeptideIdentification> peptide_identifications2; 
String date_string_1;
String date_string_2;
PeptideHit peptide_hit;

CHECK((XTandemInfile()))
	ptr = new XTandemInfile();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~XTandemInfile())
	delete ptr;
RESULT

ptr = new XTandemInfile();

CHECK(void setFragmentMassTolerance(double tolerance))
	ptr->setFragmentMassTolerance(13.0);
	TEST_REAL_EQUAL(ptr->getFragmentMassTolerance(), 13.0)
RESULT
 
CHECK(double getFragmentMassTolerance() const)
	NOT_TESTABLE
RESULT

CHECK(void setPrecursorMassTolerancePlus(double tol))
	ptr->setPrecursorMassTolerancePlus(14.0);
	TEST_REAL_EQUAL(ptr->getPrecursorMassTolerancePlus(), 14.0)
RESULT

CHECK(double getPrecursorMassTolerancePlus() const)
	NOT_TESTABLE
RESULT

CHECK(void setPrecursorMassToleranceMinus(double tol))
	ptr->setPrecursorMassToleranceMinus(15.0);
	TEST_REAL_EQUAL(ptr->getPrecursorMassToleranceMinus(), 15.0)
RESULT

CHECK(double getPrecursorMassToleranceMinus() const)
	NOT_TESTABLE
RESULT

CHECK(void setPrecursorMassErrorUnit(ErrorUnit unit))
	ptr->setPrecursorMassErrorUnit(XTandemInfile::DALTONS);
	TEST_EQUAL(ptr->getPrecursorMassErrorUnit(), XTandemInfile::DALTONS)
	ptr->setPrecursorMassErrorUnit(XTandemInfile::PPM);
	TEST_EQUAL(ptr->getPrecursorMassErrorUnit(), XTandemInfile::PPM)
RESULT

CHECK(ErrorUnit getPrecursorMassErrorUnit() const)
	NOT_TESTABLE
RESULT
  
CHECK(void setNumberOfThreads(UInt threads))
	ptr->setNumberOfThreads(16);
	TEST_EQUAL(ptr->getNumberOfThreads(), 16)
RESULT

CHECK(UInt getNumberOfThreads() const)
	NOT_TESTABLE
RESULT

CHECK(void setModifications(const ModificationDefinitionsSet &mods))
	ModificationDefinitionsSet sets("MOD:00720,MOD:00719", "MOD:01061,MOD:01060");
	ptr->setModifications(sets);
	TEST_EQUAL(ptr->getModifications() == sets, true)
RESULT

CHECK(void setOutputFilename(const String &output))
	ptr->setOutputFilename("blubb_new_outputfilename");
	TEST_STRING_EQUAL(ptr->getOutputFilename(), "blubb_new_outputfilename")
RESULT

CHECK(const String& getOutputFilename() const)
	NOT_TESTABLE
RESULT

CHECK(void setInputFilename(const String &input_file))
	ptr->setInputFilename("blubb_new_inputfilename");
	TEST_STRING_EQUAL(ptr->getInputFilename(), "blubb_new_inputfilename")
RESULT

CHECK(const String& getInputFilename() const)
	NOT_TESTABLE
RESULT

CHECK(void setTaxonomyFilename(const String &filename))
	ptr->setTaxonomyFilename("blubb_new_taxonomy_file");
	TEST_STRING_EQUAL(ptr->getTaxonomyFilename(), "blubb_new_taxonomy_file")
RESULT

CHECK(const String& getTaxonomyFilename() const)
	NOT_TESTABLE
RESULT

CHECK(void setDefaultParametersFilename(const String &filename))
	ptr->setDefaultParametersFilename("blubb_new_default_parameters_file");
	TEST_STRING_EQUAL(ptr->getDefaultParametersFilename(), "blubb_new_default_parameters_file")
RESULT

CHECK(const String& getDefaultParametersFilename() const)
	NOT_TESTABLE
RESULT

CHECK(void setTaxon(const String &taxon))
	ptr->setTaxon("blubb_taxon");
	TEST_STRING_EQUAL(ptr->getTaxon(), "blubb_taxon")
RESULT

CHECK(const String& getTaxon() const)
	NOT_TESTABLE
RESULT

CHECK(void setMaxPrecursorCharge(Int max_charge))
	ptr->setMaxPrecursorCharge(17);
	TEST_EQUAL(ptr->getMaxPrecursorCharge(), 17)
RESULT

CHECK(Int getMaxPrecursorCharge() const)
	NOT_TESTABLE
RESULT

CHECK(void setNumberOfMissedCleavages(UInt missed_cleavages))
	ptr->setNumberOfMissedCleavages(18);
	TEST_EQUAL(ptr->getNumberOfMissedCleavages(), 18)
RESULT

CHECK(UInt getNumberOfMissedCleavages() const)
	NOT_TESTABLE
RESULT

CHECK(void setMaxValidEValue(double value))
	ptr->setMaxValidEValue(19.0);
	TEST_REAL_EQUAL(ptr->getMaxValidEValue(), 19.0)
RESULT

CHECK(double getMaxValidEValue() const)
	NOT_TESTABLE
RESULT

CHECK(void setPrecursorErrorType(MassType mono_isotopic))
	XTandemInfile::MassType mono = XTandemInfile::MONOISOTOPIC;
	XTandemInfile::MassType average = XTandemInfile::AVERAGE;
	ptr->setPrecursorErrorType(mono);
	TEST_EQUAL(ptr->getPrecursorErrorType(), mono)
	ptr->setPrecursorErrorType(average);
	TEST_EQUAL(ptr->getPrecursorErrorType(), average)
RESULT

CHECK(MassType getPrecursorErrorType() const)
	NOT_TESTABLE
RESULT

CHECK(void setFragmentMassErrorUnit(ErrorUnit unit))
	XTandemInfile::ErrorUnit daltons = XTandemInfile::DALTONS;
	XTandemInfile::ErrorUnit ppm = XTandemInfile::PPM;
	ptr->setFragmentMassErrorUnit(daltons);
	TEST_EQUAL(ptr->getFragmentMassErrorUnit(), daltons)
	ptr->setFragmentMassErrorUnit(ppm);
	TEST_EQUAL(ptr->getFragmentMassErrorUnit(), ppm)
RESULT

CHECK(ErrorUnit getFragmentMassErrorUnit() const)
	NOT_TESTABLE
RESULT

CHECK(const ModificationDefinitionsSet& getModifications() const)
	NOT_TESTABLE
RESULT

CHECK(void write(const String &filename))
	string filename("XTandemInfile_test.tmp");
	NEW_TMP_FILE(filename);
	ptr->write(filename);
	XTandemInfile file;
	file.load(filename);
	NOT_TESTABLE
RESULT

CHECK(void load(const String &filename))
	ptr->load("data/XTandemInfile_test.xml");
	NOT_TESTABLE
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
