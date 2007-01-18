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
// $Maintainer: Andreas Bertsch $
// --------------------------------------------------------------------------
//

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/SparseVector.h>

///////////////////////////


///////////////////////////
START_TEST(SparseVector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

SparseVector* svp;

CHECK(SparseVector::SparseVector())
  svp = new SparseVector();
  TEST_NOT_EQUAL(svp, 0)
RESULT

CHECK(SparseVector::operator[]())
  TEST_EQUAL(svp->nonzero_size(),0)
  (*svp)[3] = 1.2;
  TEST_REAL_EQUAL((*svp)[3],1.2)
  (*svp)[5] = 1.1;
  TEST_REAL_EQUAL((*svp)[5],1.1)
  (*svp)[1] = 1.3;
  TEST_REAL_EQUAL((*svp)[1],1.3)
RESULT


CHECK(SparseVector::~SparseVector())
  delete svp;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
