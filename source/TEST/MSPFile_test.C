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

#include <OpenMS/FORMAT/MSPFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/KERNEL/StandardTypes.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(MSPFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MSPFile* ptr = 0;
START_SECTION((MSPFile()))
	ptr = new MSPFile;
	TEST_NOT_EQUAL(ptr, 0)
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

START_SECTION(void load(const String &filename, std::vector< PeptideIdentification > &ids, RichPeakMap &exp))
	MSPFile msp_file;
	vector<PeptideIdentification> ids;
	RichPeakMap exp;
	msp_file.load(OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"), ids, exp);
	TEST_EQUAL(exp.size(), 5)
	TEST_EQUAL(ids.size(), 5)

		//test DocumentIdentifier addition
	TEST_STRING_EQUAL(exp.getLoadedFilePath(), OPENMS_GET_TEST_DATA_PATH("MSPFile_test.msp"));
	TEST_STRING_EQUAL(FileHandler::typeToName(exp.getLoadedFileType()),"MSP");


	TEST_STRING_EQUAL(exp[0].getNativeID(), "index=0")
	TEST_STRING_EQUAL(exp[1].getNativeID(), "index=1")
	TEST_STRING_EQUAL(exp[2].getNativeID(), "index=2")
	TEST_STRING_EQUAL(exp[3].getNativeID(), "index=3")
	TEST_STRING_EQUAL(exp[4].getNativeID(), "index=4")

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
	TEST_EQUAL(exp.size(), 3)
	TEST_EQUAL(ids.size(), 3)

	TEST_STRING_EQUAL(exp[0].getNativeID(), "index=2")
	TEST_STRING_EQUAL(exp[1].getNativeID(), "index=3")
	TEST_STRING_EQUAL(exp[2].getNativeID(), "index=4")

END_SECTION

START_SECTION(void store(const String& filename, const RichPeakMap& exp) const)
	MSPFile msp_file;
	vector<PeptideIdentification> ids;
  RichPeakMap exp;
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
	TEST_EQUAL(ids.size(), 5)
	TEST_EQUAL(exp.size(), 5)

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
