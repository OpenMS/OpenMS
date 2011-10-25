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
// $Maintainer: Stephan Aiche $
// $Authors: Stephan Aiche $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/IMS/IMSAlphabetTextParser.h>
///////////////////////////

using namespace OpenMS;
using namespace ims;
using namespace std;

START_TEST(IMSAlphabetTextParser, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

IMSAlphabetTextParser* ptr = 0;
IMSAlphabetTextParser* null_ptr = 0;
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



