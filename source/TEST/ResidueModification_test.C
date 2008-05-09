// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CHEMISTRY/Residue.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(Residue, "$Id$")

/////////////////////////////////////////////////////////////

// Modification tests
ResidueModification* m_ptr = 0;
CHECK(ResidueModification())
  m_ptr = new ResidueModification();
	TEST_NOT_EQUAL(m_ptr, 0)
RESULT

CHECK(~ResidueModification())
	delete m_ptr;
RESULT

m_ptr = new ResidueModification();

CHECK(ResidueModification(const ResidueModification& modification))
  ResidueModification m(*m_ptr);
	TEST_EQUAL(m == *m_ptr, true)
RESULT

CHECK(ResidueModification& operator = (const ResidueModification& modification))
	ResidueModification m;
	m = *m_ptr;
	TEST_EQUAL(m == *m_ptr, true)
RESULT

END_TEST
