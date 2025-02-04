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
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetParser.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

class IMSAlphabetParserImpl
    : public IMSAlphabetParser<>
{
private:
  ContainerType elements_;
public:
  ContainerType& getElements() override
  {
    return elements_;
  }

  void parse(std::istream& ) override
  {
    // ignore istream, just enter something into the map
    elements_.insert(std::make_pair("A", 71.03711));
    elements_.insert(std::make_pair("R", 156.10111));
  }
};

START_TEST(IMSAlphabetParser, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// this is an abstract class, that only provides the load method
// it cannot be instanciated so it cannot be tested, therefor we
// test the implementation from above


IMSAlphabetParser<>* ptr = nullptr;
IMSAlphabetParser<>* null_ptr = nullptr;

START_SECTION(IMSAlphabetParser())
{
  ptr = new IMSAlphabetParserImpl();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~IMSAlphabetParser())
{
  delete ptr;
}
END_SECTION

IMSAlphabetParser<> * parser = new IMSAlphabetParserImpl();

START_SECTION((void load(const std::string &fname)))
{
  TEST_EXCEPTION(Exception::IOException ,parser->load(""))

  String filename;
  NEW_TMP_FILE(filename)

  // just create the file
  ofstream of;
  of.open(filename.c_str());
  of << "just text" << std::endl;
  of.close();

  parser->load(filename);

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

START_SECTION((virtual void parse(InputSource &is)))
{
  // already tested by load
  NOT_TESTABLE
}
END_SECTION

delete parser;


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



