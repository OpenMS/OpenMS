// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Andreas Bertsch, Daniel Jameson, Chris Bielow$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MascotRemoteQuery.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MascotRemoteQuery, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MascotRemoteQuery* ptr = nullptr;
MascotRemoteQuery* nullPointer = nullptr;
START_SECTION(MascotRemoteQuery(QObject *parent=0))
{
	ptr = new MascotRemoteQuery();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~MascotRemoteQuery())
{
	delete ptr;
}
END_SECTION

START_SECTION((void setQuerySpectra(const String &exp)))
{
  MascotRemoteQuery query;
	query.setQuerySpectra("BEGIN IONS\n1 1\n1 1\nEND IONS");
	NOT_TESTABLE
}
END_SECTION

START_SECTION((const QByteArray& getMascotXMLResponse() const ))
{
  MascotRemoteQuery query;
	TEST_EQUAL(query.getMascotXMLResponse().size(), 0)
}
END_SECTION

START_SECTION((bool hasError() const ))
{
  MascotRemoteQuery query;
	TEST_EQUAL(query.hasError(), false)
}
END_SECTION

START_SECTION((const String& getErrorMessage() const ))
{
  MascotRemoteQuery query;
	TEST_STRING_EQUAL(query.getErrorMessage(), "")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



