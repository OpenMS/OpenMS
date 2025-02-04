// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/SYSTEM/FileWatcher.h>

/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ControlledVocabulary, "$Id$")

/////////////////////////////////////////////////////////////

FileWatcher* ptr = nullptr;
FileWatcher* nullPointer = nullptr;
START_SECTION(FileWatcher(QObject *parent=0))
	ptr = new FileWatcher();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~FileWatcher())
	delete ptr;
END_SECTION

START_SECTION(void setDelayInSeconds(double delay))
	NOT_TESTABLE
END_SECTION

START_SECTION(void addFile(const String& path))
	NOT_TESTABLE
END_SECTION

START_SECTION(void removeFile(const String& path))
	NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
