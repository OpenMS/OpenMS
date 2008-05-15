// -*- Mode: C++; tab-width: 2; -*-
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

CHECK(void setPrecursorMassErrorUnit(ERROR_UNIT unit))
	ptr->setPrecursorMassErrorUnit(XTandemInfile::DALTONS);
	TEST_EQUAL(ptr->getPrecursorMassErrorUnit(), XTandemInfile::DALTONS)
	ptr->setPrecursorMassErrorUnit(XTandemInfile::PPM);
	TEST_EQUAL(ptr->getPrecursorMassErrorUnit(), XTandemInfile::PPM)
RESULT

CHECK(ERROR_UNIT getPrecursorMassErrorUnit() const)
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

CHECK(const ModificationDefinitionsSet& getFixedModifications() const)
	NOT_TESTABLE
RESULT

/*
    - 'void setOutputFilename(const String &output)'
    - 'const String& getOutputFilename() const '
    - 'void setInputFilename(const String &input_file)'
    - 'const String& getInputFilename() const '
    - 'void setTaxonomyFilename(const String &filename)'
    - 'const String& getTaxonomyFilename() const '
    - 'void setDefaultParametersFilename(const String &filename)'
    - 'const String& getDefaultParametersFilename() const '
    - 'void setTaxon(const String &taxon)'
    - 'const String& getTaxon() const '
    - 'void setMaxPrecursorCharge(Int max_charge)'
    - 'Int getMaxPrecursorCharge() const '
    - 'void setNumberOfMissedCleavages(UInt missed_cleavages)'
    - 'UInt getNumberOfMissedCleavages() const '
    - 'void setMaxValidEValue(double value)'
    - 'double getMaxValidEValue() const '
    - 'void write(const String &filename) throw (Exception::UnableToCreateFile)'
    - 'void load(const String &filename) throw (Exception::FileNotFound, Exception::ParseError)
*/

CHECK(void load(const String& filename, ProteinIdentification& protein_identification, std::vector<PeptideIdentification>& id_data) const throw(Exception::FileNotFound, Exception::ParseError))
	//ptr->load("/home/andreas/DATA/OpenMS/share/OpenMS/FORMAT/XTandem_default_input.xml");
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
