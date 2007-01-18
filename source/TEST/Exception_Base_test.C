// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/KERNEL/DPeak.h>

///////////////////////////

START_TEST(DPeak<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

Exception::Base* e_ptr = 0;
CHECK(Base() throw())
	e_ptr = new Exception::Base;
	TEST_NOT_EQUAL(e_ptr, 0)
RESULT

CHECK(~Base() throw())
	delete e_ptr;
RESULT


CHECK(Base(const Base& exception) throw())
  // ???
RESULT

CHECK((Base(const char* file, int line, const char* function) throw()))
  // ???
RESULT

CHECK((Base(const char* file, int line, const char* function, const std::string& name, const std::string& message) throw()))
  // ???
RESULT

CHECK(const char* getFile() const throw())
  // ???
RESULT

CHECK(const char* getName() const throw())
  // ???
RESULT

CHECK(const char* what() const throw())
  // ???
RESULT

CHECK(int getLine() const throw())
  // ???
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
