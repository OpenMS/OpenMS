// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Stephan Aiche$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

using namespace OpenMS;

START_TEST(Exception::Base, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Exception::BaseException* e_ptr = 0;
START_SECTION(Base() )
	e_ptr = new Exception::BaseException;
	TEST_NOT_EQUAL(e_ptr, 0)
END_SECTION

START_SECTION(~Base() )
	delete e_ptr;
END_SECTION


START_SECTION(Base(const Base& exception) )
  // ???
END_SECTION

START_SECTION((Base(const char* file, int line, const char* function) ))
  // ???
END_SECTION

START_SECTION((Base(const char* file, int line, const char* function, const std::string& name, const std::string& message) ))
  // ???
END_SECTION

START_SECTION(const char* getFile() const )
  // ???
END_SECTION

START_SECTION(const char* getName() const )
  // ???
END_SECTION

START_SECTION(const char* what() const )
  // ???
END_SECTION

START_SECTION(int getLine() const )
  // ???
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
