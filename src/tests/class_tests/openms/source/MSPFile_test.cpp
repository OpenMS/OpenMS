// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
///////////////////////////

#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/FileTypes.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MSPFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSPFile* ptr = nullptr;
MSPFile* nullPointer = nullptr;
START_SECTION((MSPFile()))
	ptr = new MSPFile;
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MSPFile()))
	delete ptr;
END_SECTION

START_SECTION(MSPFile(const MSPFile &rhs))
	MSPFile f1, f2;
	Param p = f1.getParameters();
	p.setValue("instrument", "it");
	f1.setParameters(p);
	TEST_EQUAL(f1.getParameters() == f2.getParameters(), false)
	MSPFile f3(f1);
	TEST_EQUAL(f1.getParameters() == f3.getParameters(), true)
END_SECTION

START_SECTION(MSPFile& operator=(const MSPFile &rhs))
	MSPFile f1, f2;
	Param p = f1.getParameters();
	p.setValue("instrument", "it");
	f1.setParameters(p);
	TEST_EQUAL(f1.getParameters() == f2.getParameters(), false)
	f2 = f1;
	TEST_EQUAL(f1.getParameters() == f2.getParameters(), true)
END_SECTION

START_SECTION(void load(const String &filename, std::vector< PeptideIdentification > &ids, PeakMap &exp))
	MSPFile msp_file;
	vector<PeptideIdentification> ids;
	PeakMap exp;
	msp_file.load(OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"), ids, exp);
	TEST_EQUAL(exp.size(), 7)
	TEST_EQUAL(ids.size(), 7)

		//test DocumentIdentifier addition
	TEST_STRING_EQUAL(exp.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"));
  TEST_STRING_EQUAL(FileTypes::typeToName(exp.getLoadedFileType()),"msp");


	TEST_STRING_EQUAL(exp[0].getNativeID(), "index=0")
	TEST_STRING_EQUAL(exp[1].getNativeID(), "index=1")
	TEST_STRING_EQUAL(exp[2].getNativeID(), "index=2")
	TEST_STRING_EQUAL(exp[3].getNativeID(), "index=3")
	TEST_STRING_EQUAL(exp[4].getNativeID(), "index=4")
  TEST_STRING_EQUAL(exp[5].getNativeID(), "index=5")
  TEST_STRING_EQUAL(exp[6].getNativeID(), "index=6")
  TEST_STRING_EQUAL(ids[5].getHits()[0].getSequence().toString(), ".(Acetyl)AAAAAAGAGPEM(Oxidation)VR")
  TEST_STRING_EQUAL( ids[5].getHits()[0].getPeakAnnotations()[0].annotation, "a3")
  TEST_STRING_EQUAL( ids[5].getHits()[0].getPeakAnnotations()[1].annotation, "b3")
  TEST_EQUAL( ids[5].getHits()[0].getPeakAnnotations()[1].charge, 1)
  // next only with parse_firstonly = false
  //TEST_STRING_EQUAL( ids[5].getHits()[0].getPeakAnnotations()[2].annotation, "y2-H2O")
  //TEST_EQUAL( ids[5].getHits()[0].getPeakAnnotations()[2].charge, 1)
  TEST_STRING_EQUAL( ids[5].getHits()[0].getPeakAnnotations()[2].annotation, "?")
  TEST_EQUAL( ids[5].getHits()[0].getPeakAnnotations()[2].charge, 0)
  TEST_STRING_EQUAL( ids[5].getHits()[0].getPeakAnnotations()[3].annotation, "y4")
  TEST_EQUAL( ids[5].getHits()[0].getPeakAnnotations()[3].charge, 2)
  TEST_STRING_EQUAL(ids[6].getHits()[0].getSequence().toString(), ".(Acetyl)AAAAAAVGPGAGGAGSAVPGGAGPC(Carbamidomethyl)ATVSVFPGAR")


	Param p(msp_file.getParameters());
	p.setValue("instrument", "qtof");
	msp_file.setParameters(p);
	ids.clear();
	exp.clear(true);
	msp_file.load(OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"), ids, exp);
	TEST_EQUAL(exp.size(), 2)
	TEST_EQUAL(ids.size(), 2)

	TEST_STRING_EQUAL(exp[0].getNativeID(), "index=0")
	TEST_STRING_EQUAL(exp[1].getNativeID(), "index=1")

	p.setValue("instrument", "it");
	msp_file.setParameters(p);
	ids.clear();
	exp.clear(true);
	msp_file.load(OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"), ids, exp);
	TEST_EQUAL(exp.size(), 5)
	TEST_EQUAL(ids.size(), 5)

	TEST_STRING_EQUAL(exp[0].getNativeID(), "index=2")
	TEST_STRING_EQUAL(exp[1].getNativeID(), "index=3")
	TEST_STRING_EQUAL(exp[2].getNativeID(), "index=4")

END_SECTION

START_SECTION(void store(const String& filename, const PeakMap& exp) const)
	MSPFile msp_file;
	vector<PeptideIdentification> ids;
  PeakMap exp;
  msp_file.load(OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"), ids, exp);
	for (Size i = 0; i != ids.size(); ++i)
	{
		exp[i].getPeptideIdentifications().push_back(ids[i]);
	}
	String filename;
	NEW_TMP_FILE(filename)
	msp_file.store(filename, exp);

	exp.clear(true);
	ids.clear();
	msp_file.load(filename, ids, exp);
	TEST_EQUAL(ids.size(), 7)
	TEST_EQUAL(exp.size(), 7)

	TEST_EQUAL(ids[0].getHits().size(), 1)
	TEST_EQUAL(ids[1].getHits().size(), 1)
	TEST_EQUAL(ids[2].getHits().size(), 1)
	TEST_EQUAL(ids[3].getHits().size(), 1)
	TEST_EQUAL(ids[4].getHits().size(), 1)
	TEST_EQUAL(ids[0].getHits().begin()->getSequence().isModified(), false)
	TEST_EQUAL(ids[1].getHits().begin()->getSequence().isModified(), false)
	TEST_EQUAL(ids[2].getHits().begin()->getSequence().isModified(), false)
	TEST_EQUAL(ids[3].getHits().begin()->getSequence().isModified(), true)
	TEST_EQUAL(ids[4].getHits().begin()->getSequence().isModified(), false)
	TEST_EQUAL(ids[0].getHits().begin()->getCharge(), 2)
	TEST_EQUAL(ids[1].getHits().begin()->getCharge(), 2)
	TEST_EQUAL(ids[2].getHits().begin()->getCharge(), 2)
	TEST_EQUAL(ids[3].getHits().begin()->getCharge(), 2)
	TEST_EQUAL(ids[4].getHits().begin()->getCharge(), 3)
	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
