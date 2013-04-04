// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2013.
// 
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution 
//    may be used to endorse or promote products derived from this software 
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS. 
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING 
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, 
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, 
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; 
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, 
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR 
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF 
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
XTandemInfile* nullPointer = 0;
ProteinIdentification protein_identification;
vector<PeptideIdentification> peptide_identifications; 
vector<PeptideIdentification> peptide_identifications2; 
String date_string_1;
String date_string_2;
PeptideHit peptide_hit;

START_SECTION((XTandemInfile()))
	ptr = new XTandemInfile();
	TEST_NOT_EQUAL(ptr, nullPointer)
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
