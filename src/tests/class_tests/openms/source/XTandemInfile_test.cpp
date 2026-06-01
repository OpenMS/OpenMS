// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/XTandemInfile.h>

///////////////////////////

START_TEST(XTandemInfile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

XTandemInfile* ptr;
XTandemInfile* nullPointer = nullptr;

START_SECTION((XTandemInfile()))
	ptr = new XTandemInfile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~XTandemInfile())
	delete ptr;
END_SECTION

XTandemInfile xml_file;

START_SECTION(void setFragmentMassTolerance(double tolerance))
	xml_file.setFragmentMassTolerance(13.0);
	TEST_REAL_SIMILAR(xml_file.getFragmentMassTolerance(), 13.0)
END_SECTION

START_SECTION(double getFragmentMassTolerance() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setPrecursorMassTolerancePlus(double tol))
	xml_file.setPrecursorMassTolerancePlus(14.0);
	TEST_REAL_SIMILAR(xml_file.getPrecursorMassTolerancePlus(), 14.0)
END_SECTION

START_SECTION(double getPrecursorMassTolerancePlus() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setPrecursorMassToleranceMinus(double tol))
	xml_file.setPrecursorMassToleranceMinus(15.0);
	TEST_REAL_SIMILAR(xml_file.getPrecursorMassToleranceMinus(), 15.0)
END_SECTION

START_SECTION(double getPrecursorMassToleranceMinus() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setPrecursorMassErrorUnit(ErrorUnit unit))
	xml_file.setPrecursorMassErrorUnit(XTandemInfile::DALTONS);
	TEST_EQUAL(xml_file.getPrecursorMassErrorUnit(), XTandemInfile::DALTONS)
	xml_file.setPrecursorMassErrorUnit(XTandemInfile::PPM);
	TEST_EQUAL(xml_file.getPrecursorMassErrorUnit(), XTandemInfile::PPM)
END_SECTION

START_SECTION(ErrorUnit getPrecursorMassErrorUnit() const)
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setNumberOfThreads(UInt threads))
	xml_file.setNumberOfThreads(16);
	TEST_EQUAL(xml_file.getNumberOfThreads(), 16)
END_SECTION

START_SECTION(UInt getNumberOfThreads() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setModifications(const ModificationDefinitionsSet& mods))
	ModificationDefinitionsSet sets(ListUtils::create<String>("Oxidation (M)"), ListUtils::create<String>("Carboxymethyl (C)"));
	xml_file.setModifications(sets);
	TEST_EQUAL(xml_file.getModifications() == sets, true)
END_SECTION

START_SECTION(const ModificationDefinitionsSet& getModifications() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setOutputFilename(const String& output))
	xml_file.setOutputFilename("blubb_new_outputfilename");
	TEST_STRING_EQUAL(xml_file.getOutputFilename(), "blubb_new_outputfilename")
END_SECTION

START_SECTION(const String& getOutputFilename() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setInputFilename(const String& input_file))
	xml_file.setInputFilename("blubb_new_inputfilename");
	TEST_STRING_EQUAL(xml_file.getInputFilename(), "blubb_new_inputfilename")
END_SECTION

START_SECTION(const String& getInputFilename() const)
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setTaxonomyFilename(const String& filename))
	xml_file.setTaxonomyFilename("blubb_new_taxonomy_file");
	TEST_STRING_EQUAL(xml_file.getTaxonomyFilename(), "blubb_new_taxonomy_file")
END_SECTION

START_SECTION(const String& getTaxonomyFilename() const)
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setDefaultParametersFilename(const String& filename))
	xml_file.setDefaultParametersFilename("blubb_new_default_parameters_file");
	TEST_STRING_EQUAL(xml_file.getDefaultParametersFilename(), "blubb_new_default_parameters_file")
END_SECTION

START_SECTION(const String& getDefaultParametersFilename() const)
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setTaxon(const String& taxon))
	xml_file.setTaxon("blubb_taxon");
	TEST_STRING_EQUAL(xml_file.getTaxon(), "blubb_taxon")
END_SECTION

START_SECTION(const String& getTaxon() const)
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setMaxPrecursorCharge(Int max_charge))
	xml_file.setMaxPrecursorCharge(17);
	TEST_EQUAL(xml_file.getMaxPrecursorCharge(), 17)
END_SECTION

START_SECTION(Int getMaxPrecursorCharge() const)
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setNumberOfMissedCleavages(UInt missed_cleavages))
	xml_file.setNumberOfMissedCleavages(18);
	TEST_EQUAL(xml_file.getNumberOfMissedCleavages(), 18)
END_SECTION

START_SECTION(UInt getNumberOfMissedCleavages() const)
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setMaxValidEValue(double value))
	xml_file.setMaxValidEValue(19.0);
	TEST_REAL_SIMILAR(xml_file.getMaxValidEValue(), 19.0)
END_SECTION

START_SECTION(double getMaxValidEValue() const)
	NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setPrecursorErrorType(MassType mono_isotopic))
	XTandemInfile::MassType mono = XTandemInfile::MONOISOTOPIC;
	XTandemInfile::MassType average = XTandemInfile::AVERAGE;
	xml_file.setPrecursorErrorType(mono);
	TEST_EQUAL(xml_file.getPrecursorErrorType(), mono)
	xml_file.setPrecursorErrorType(average);
	TEST_EQUAL(xml_file.getPrecursorErrorType(), average)
END_SECTION

START_SECTION(MassType getPrecursorErrorType() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void setFragmentMassErrorUnit(ErrorUnit unit))
	XTandemInfile::ErrorUnit daltons = XTandemInfile::DALTONS;
	XTandemInfile::ErrorUnit ppm = XTandemInfile::PPM;
	xml_file.setFragmentMassErrorUnit(daltons);
	TEST_EQUAL(xml_file.getFragmentMassErrorUnit(), daltons)
	xml_file.setFragmentMassErrorUnit(ppm);
	TEST_EQUAL(xml_file.getFragmentMassErrorUnit(), ppm)
END_SECTION

START_SECTION(ErrorUnit getFragmentMassErrorUnit() const)
  NOT_TESTABLE // tested above
END_SECTION

  START_SECTION(void write(const String& filename, bool ignore_member_parameters = false, bool force_default_mods = false))
	String filename("XTandemInfile_test.tmp");
	NEW_TMP_FILE(filename);
  ModificationDefinitionsSet sets(ListUtils::create<String>("Oxidation (M),Dimethyl (N-term),Carboxymethyl (C)"), ListUtils::create<String>("Ammonium (C-term),ICDID (C)"));
  xml_file.setModifications(sets);
  xml_file.setAllowIsotopeError(true);
  xml_file.setSemiCleavage(true);
  xml_file.setOutputResults("all");
	xml_file.write(filename);
  TEST_FILE_SIMILAR(filename.c_str(), OPENMS_GET_TEST_DATA_PATH("XTandemInfile_test_write.xml"))
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
