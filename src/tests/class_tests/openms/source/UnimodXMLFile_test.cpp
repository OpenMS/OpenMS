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

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/FORMAT/UnimodXMLFile.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>

#include <vector>

///////////////////////////

START_TEST(UnimodXMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

UnimodXMLFile xml_file;
UnimodXMLFile* ptr;
UnimodXMLFile* nullPointer = nullptr;

START_SECTION((UnimodXMLFile()))
	ptr = new UnimodXMLFile();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~UnimodXMLFile())
	delete ptr;
END_SECTION

ptr = new UnimodXMLFile();

START_SECTION(void load(const String& filename, vector<ResidueModification*>& modifications))
	vector<ResidueModification*> modifications;
	ptr->load("CHEMISTRY/unimod.xml", modifications);

	//cerr << "#modifications read: " << modifications.size() << endl;
	//for (vector<ResidueModification>::const_iterator it = modifications.begin(); it != modifications.end(); ++it)
	//{
	//	cerr << it->getTitle() << "\t" << it->getFullName() << "\t" << it->getAllowedPositionName() << "\t" << it->getSite() << "\t" << it->getClassification() << "\t" << it->getComposition() << "\t" << it->getMonoMass() << endl;
	//}

	TEST_EQUAL(modifications.size() > 1, true)

  // cleanup
  for (Size k = 0; k < modifications.size(); k++)
  {
    delete modifications[k];
  }
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
