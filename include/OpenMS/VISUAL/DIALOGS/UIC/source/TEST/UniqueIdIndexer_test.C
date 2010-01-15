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
// $Maintainer: Clemens Groepl $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/CONCEPT/UniqueIdIndexer.h>

#include <vector>

///////////////////////////

using namespace OpenMS;
using namespace std;

namespace OpenMS
{

struct Dummy
  : UniqueIdInterface
{
    Dummy () : dummy(0) {};
    Size dummy;
};

class DummyVectorIndexed
: public vector<Dummy>,
  public UniqueIdIndexer<DummyVectorIndexed>
{
};

// this is used for testing purposes only
template < typename T >
struct CanAccessTheUniqueIdMap : T
{
  DummyVectorIndexed::UniqueIdMap & getUniqueIdMap() { return this->uniqueid_to_index_; }
};

// this saves us some typing
template < typename T >
CanAccessTheUniqueIdMap<T>& canAccessTheUniqueIdMap(T &rhs)
{
  return static_cast<CanAccessTheUniqueIdMap<T>&>(rhs);
}

}

START_TEST(UniqueIdIndexer, "$Id: UniqueIdIndexer_test.C 6446 2009-11-20 16:21:41Z andreas_bertsch $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DummyVectorIndexed* ptr = 0;
START_SECTION((UniqueIdIndexer()))
{
  ptr = new DummyVectorIndexed();
  TEST_NOT_EQUAL(ptr, 0)
}
END_SECTION

START_SECTION((~UniqueIdIndexer()))
{
  delete ptr;
}
END_SECTION

START_SECTION((Size uniqueIdToIndex(UInt64 unique_id) const))
{
  DummyVectorIndexed dvi;
  Size num_uii = 10;
  dvi.resize(num_uii);
  for ( Size i = 0; i < num_uii; ++i )
  {
    dvi[i].dummy = i;
    dvi[i].setUniqueId(10*i+1000);
  }

  for ( Size i = 0; i < num_uii; ++i )
  {
    TEST_EQUAL(dvi.uniqueIdToIndex(10*i+1000),i);
  }

  STATUS("shuffling ...");
  std::random_shuffle(dvi.begin(), dvi.end());

  for ( Size i = 0; i < num_uii; ++i )
  {
    Dummy const & current_dummy = dvi[i];
    TEST_EQUAL(dvi.uniqueIdToIndex(current_dummy.getUniqueId()),i);
    // TEST_EQUAL(dvidvi[current_dummy.dummy].getUniqueId()),current_dummy.dummy);
  }

  dvi.pop_back();
  dvi.pop_back();

  dvi.push_back(Dummy());
  dvi.back().setUniqueId(12345678);

  dvi.push_back(Dummy());
  dvi.push_back(Dummy());
  dvi.back().setUniqueId(12345678);
  dvi.push_back(Dummy());

  STATUS("shuffling ...");
  std::random_shuffle(dvi.begin(), dvi.end());

  TEST_EXCEPTION_WITH_MESSAGE(Exception::Postcondition,dvi.updateUniqueIdToIndex(),"Duplicate valid unique ids detected!   RandomAccessContainer has size()==12, num_valid_unique_id==10, uniqueid_to_index_.size()==9");

}
END_SECTION

START_SECTION((void updateUniqueIdToIndex() const))
{
  // see uniqueIdToIndex()
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((void swap(UniqueIdIndexer &rhs)))
{

  DummyVectorIndexed dvi;
  Size num_uii = 10;

  dvi.resize(num_uii);

  for ( Size i = 0; i < num_uii; ++i )
  {
    dvi[i].dummy = i;
    dvi[i].setUniqueId(10*i+1000);
  }

  dvi.updateUniqueIdToIndex();

  DummyVectorIndexed dvi2;

  TEST_EQUAL(canAccessTheUniqueIdMap(dvi).getUniqueIdMap().size(),num_uii);
  TEST_EQUAL(canAccessTheUniqueIdMap(dvi2).getUniqueIdMap().size(),0);

  std::swap(dvi, dvi2);

  TEST_EQUAL(canAccessTheUniqueIdMap(dvi).getUniqueIdMap().size(),0);
  TEST_EQUAL(canAccessTheUniqueIdMap(dvi2).getUniqueIdMap().size(),num_uii);

  dvi = dvi2;

  TEST_EQUAL(canAccessTheUniqueIdMap(dvi).getUniqueIdMap().size(),num_uii);
  canAccessTheUniqueIdMap(dvi).getUniqueIdMap().clear();
  TEST_EQUAL(canAccessTheUniqueIdMap(dvi).getUniqueIdMap().size(),0);

  TEST_EQUAL(dvi.uniqueIdToIndex(4321234324124ull),Size(-1));

  TEST_EQUAL(canAccessTheUniqueIdMap(dvi).getUniqueIdMap().size(),num_uii);

}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
