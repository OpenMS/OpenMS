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

#include <OpenMS/DATASTRUCTURES/BinnedSparseVector.h>

///////////////////////////


///////////////////////////
START_TEST(BinnedSparseVector, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

BinnedSparseVector* sv2p;

CHECK(BinnedSparseVector())
  sv2p = new BinnedSparseVector();
  TEST_NOT_EQUAL(sv2p, 0)
RESULT

CHECK(void resize())
  TEST_EQUAL(sv2p->size(),0)
  sv2p->resize(10);
  TEST_EQUAL(sv2p->size(),10)
RESULT

CHECK(operator[]())
  TEST_EQUAL(sv2p->nonzero_size(),0)
  (*sv2p)[3] = 1.2;
  TEST_EQUAL(sv2p->nonzero_size(),1)
  (*sv2p)[9] = 1.2;
  TEST_EQUAL(sv2p->nonzero_size(),2)
RESULT
 
CHECK(push_back())
  sv2p->push_back(1.3);
  TEST_EQUAL(sv2p->size(),11)
  TEST_EQUAL(sv2p->nonzero_size(),3)
RESULT

CHECK(contents + size)
  uint i;
  for ( i = 0; i < sv2p->size(); ++i )
  {
    if ( i == 3 || i == 9 )
    {
      TEST_REAL_EQUAL(1.2,(*sv2p)[i])
    }
    else if ( i == 10)
    {
      TEST_REAL_EQUAL(1.3,(*sv2p)[i])
    }
    else
    {
      TEST_REAL_EQUAL(0,(*sv2p)[i])
    }
  }
RESULT

CHECK(void clear())
  sv2p->clear();
  TEST_EQUAL(sv2p->size(),0)
  TEST_EQUAL(sv2p->nonzero_size(),0)
RESULT

CHECK(iterator)
  sv2p->resize(10);
  uint i = 0;
  for ( BinnedSparseVector::iterator vit = sv2p->begin(); vit != sv2p->end();++vit)
  {
    if ( i == 2 || i == 5 || i == 7 ) *vit = 1.1;
    else *vit = 0;
    ++i;
  }
  TEST_EQUAL(i,10)
  TEST_EQUAL(sv2p->size(),10)
  TEST_EQUAL(sv2p->nonzero_size(),3)
  i = 0;
  for ( BinnedSparseVector::iterator vit = sv2p->begin(); vit != sv2p->end();++vit)
  {
    if ( i == 2 || i == 5 || i == 7 ) 
    {
      TEST_REAL_EQUAL(*vit,1.1)
    }
    else 
    {
      TEST_REAL_EQUAL(*vit,0)
    }
    ++i;
  }
  TEST_EQUAL(3,sv2p->nonzero_size())
RESULT

CHECK(const_iterator)
  uint i = 0;
  for ( BinnedSparseVector::const_iterator cvit = sv2p->begin(); cvit != sv2p->end();++cvit)
  {
    if ( i == 2 || i == 5 || i == 7 ) 
    {
      TEST_REAL_EQUAL(*cvit,1.1)
    }
    else 
    {
      TEST_REAL_EQUAL(*cvit,0)
    }
    ++i;
  }
  TEST_EQUAL(3,sv2p->nonzero_size())
  TEST_EQUAL(i,10)
RESULT

CHECK(~BinnedSparseVector())
  delete sv2p;
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
