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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
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

START_SECTION((XTandemInfile()))
	ptr = new XTandemInfile();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION(~XTandemInfile())
	delete ptr;
END_SECTION

ptr = new XTandemInfile();

START_SECTION(void setFragmentMassTolerance(double tolerance))
	ptr->setFragmentMassTolerance(13.0);
	TEST_REAL_SIMILAR(ptr->getFragmentMassTolerance(), 13.0)
END_SECTION
 
START_SECTION(double getFragmentMassTolerance() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setPrecursorMassTolerancePlus(double tol))
	ptr->setPrecursorMassTolerancePlus(14.0);
	TEST_REAL_SIMILAR(ptr->getPrecursorMassTolerancePlus(), 14.0)
END_SECTION

START_SECTION(double getPrecursorMassTolerancePlus() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setPrecursorMassToleranceMinus(double tol))
	ptr->setPrecursorMassToleranceMinus(15.0);
	TEST_REAL_SIMILAR(ptr->getPrecursorMassToleranceMinus(), 15.0)
END_SECTION

START_SECTION(double getPrecursorMassToleranceMinus() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setPrecursorMassErrorUnit(ErrorUnit unit))
	ptr->setPrecursorMassErrorUnit(XTandemInfile::DALTONS);
	TEST_EQUAL(ptr->getPrecursorMassErrorUnit(), XTandemInfile::DALTONS)
	ptr->setPrecursorMassErrorUnit(XTandemInfile::PPM);
	TEST_EQUAL(ptr->getPrecursorMassErrorUnit(), XTandemInfile::PPM)
END_SECTION

START_SECTION(ErrorUnit getPrecursorMassErrorUnit() const)
	NOT_TESTABLE
END_SECTION
  
START_SECTION(void setNumberOfThreads(UInt threads))
	ptr->setNumberOfThreads(16);
	TEST_EQUAL(ptr->getNumberOfThreads(), 16)
END_SECTION

START_SECTION(UInt getNumberOfThreads() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setModifications(const ModificationDefinitionsSet &mods))
	ModificationDefinitionsSet sets("Oxidation (M)", "Carboxymethyl (C)");
	ptr->setModifications(sets);
	TEST_EQUAL(ptr->getModifications() == sets, true)
END_SECTION

START_SECTION(void setOutputFilename(const String &output))
	ptr->setOutputFilename("blubb_new_outputfilename");
	TEST_STRING_EQUAL(ptr->getOutputFilename(), "blubb_new_outputfilename")
END_SECTION

START_SECTION(const String& getOutputFilename() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setInputFilename(const String &input_file))
	ptr->setInputFilename("blubb_new_inputfilename");
	TEST_STRING_EQUAL(ptr->getInputFilename(), "blubb_new_inputfilename")
END_SECTION

START_SECTION(const String& getInputFilename() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setTaxonomyFilename(const String &filename))
	ptr->setTaxonomyFilename("blubb_new_taxonomy_file");
	TEST_STRING_EQUAL(ptr->getTaxonomyFilename(), "blubb_new_taxonomy_file")
END_SECTION

START_SECTION(const String& getTaxonomyFilename() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setDefaultParametersFilename(const String &filename))
	ptr->setDefaultParametersFilename("blubb_new_default_parameters_file");
	TEST_STRING_EQUAL(ptr->getDefaultParametersFilename(), "blubb_new_default_parameters_file")
END_SECTION

START_SECTION(const String& getDefaultParametersFilename() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setTaxon(const String &taxon))
	ptr->setTaxon("blubb_taxon");
	TEST_STRING_EQUAL(ptr->getTaxon(), "blubb_taxon")
END_SECTION

START_SECTION(const String& getTaxon() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setMaxPrecursorCharge(Int max_charge))
	ptr->setMaxPrecursorCharge(17);
	TEST_EQUAL(ptr->getMaxPrecursorCharge(), 17)
END_SECTION

START_SECTION(Int getMaxPrecursorCharge() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setNumberOfMissedCleavages(UInt missed_cleavages))
	ptr->setNumberOfMissedCleavages(18);
	TEST_EQUAL(ptr->getNumberOfMissedCleavages(), 18)
END_SECTION

START_SECTION(UInt getNumberOfMissedCleavages() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setMaxValidEValue(double value))
	ptr->setMaxValidEValue(19.0);
	TEST_REAL_SIMILAR(ptr->getMaxValidEValue(), 19.0)
END_SECTION

START_SECTION(double getMaxValidEValue() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setPrecursorErrorType(MassType mono_isotopic))
	XTandemInfile::MassType mono = XTandemInfile::MONOISOTOPIC;
	XTandemInfile::MassType average = XTandemInfile::AVERAGE;
	ptr->setPrecursorErrorType(mono);
	TEST_EQUAL(ptr->getPrecursorErrorType(), mono)
	ptr->setPrecursorErrorType(average);
	TEST_EQUAL(ptr->getPrecursorErrorType(), average)
END_SECTION

START_SECTION(MassType getPrecursorErrorType() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void setFragmentMassErrorUnit(ErrorUnit unit))
	XTandemInfile::ErrorUnit daltons = XTandemInfile::DALTONS;
	XTandemInfile::ErrorUnit ppm = XTandemInfile::PPM;
	ptr->setFragmentMassErrorUnit(daltons);
	TEST_EQUAL(ptr->getFragmentMassErrorUnit(), daltons)
	ptr->setFragmentMassErrorUnit(ppm);
	TEST_EQUAL(ptr->getFragmentMassErrorUnit(), ppm)
END_SECTION

START_SECTION(ErrorUnit getFragmentMassErrorUnit() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(const ModificationDefinitionsSet& getModifications() const)
	NOT_TESTABLE
END_SECTION

START_SECTION(void write(const String &filename))
	string filename("XTandemInfile_test.tmp");
	NEW_TMP_FILE(filename);
	ptr->write(filename);
	XTandemInfile file;
	file.load(filename);
	NOT_TESTABLE
END_SECTION

START_SECTION(void load(const String &filename))
	ptr->load(OPENMS_GET_TEST_DATA_PATH("XTandemInfile_test.xml"));
	NOT_TESTABLE
END_SECTION

START_SECTION(bool isRefining() const )
  XTandemInfile file;
  TEST_EQUAL(file.isRefining()==true, true)
END_SECTION

START_SECTION(void setRefine(const bool refine))
  XTandemInfile file;
  file.setRefine(false);
  TEST_EQUAL(file.isRefining()==false, true)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
