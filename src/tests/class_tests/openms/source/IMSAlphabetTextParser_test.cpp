// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetTextParser.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

START_TEST(IMSAlphabetTextParser, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IMSAlphabetTextParser* ptr = nullptr;
IMSAlphabetTextParser* null_ptr = nullptr;
START_SECTION(IMSAlphabetTextParser())
{
	ptr = new IMSAlphabetTextParser();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IMSAlphabetTextParser())
{
	delete ptr;
}
END_SECTION

IMSAlphabetParser<> * parser = new IMSAlphabetTextParser();

START_SECTION((virtual void parse(std::istream &is)))
{
  String filename;
  NEW_TMP_FILE(filename)

  // just create the file
  ofstream of;
  of.open(filename.c_str());
  of << "# a comment which should be ignored" << std::endl;
  of << "A\t71.03711" << std::endl;
  of << "R\t156.10111" << std::endl;
  of.close();

  ifstream ifs;
  ifs.open(filename.c_str());

  ABORT_IF(!ifs)

  parser->parse(ifs);

  ifs.close();

  TEST_EQUAL(parser->getElements().empty(), false)
}
END_SECTION

START_SECTION((virtual ContainerType& getElements()))
{
  TEST_EQUAL(parser->getElements().size(), 2)
  TEST_REAL_SIMILAR(parser->getElements()["A"], 71.03711)
  TEST_REAL_SIMILAR(parser->getElements()["R"], 156.10111)
}
END_SECTION

delete parser;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



