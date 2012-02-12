// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

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
  virtual ContainerType& getElements()
  {
    return elements_;
  }

  virtual void parse(std::istream& )
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


IMSAlphabetParser<>* ptr = 0;
IMSAlphabetParser<>* null_ptr = 0;

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



