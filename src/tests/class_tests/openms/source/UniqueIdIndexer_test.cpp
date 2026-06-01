// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/UniqueIdIndexer.h>
#include <OpenMS/DATASTRUCTURES/ExposedVector.h>
#include <OpenMS/MATH/MathFunctions.h>

#include <vector>

///////////////////////////

using namespace OpenMS;
using namespace std;

namespace OpenMS
{

struct Dummy
  : UniqueIdInterface
{
    Dummy () : 
      dummy(0) 
    {}

    Size dummy;
};


class DummyVectorIndexed :
  public ExposedVector<Dummy>,
  public UniqueIdIndexer<DummyVectorIndexed>
{
public:
  EXPOSED_VECTOR_INTERFACE(Dummy)
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

START_TEST(UniqueIdIndexer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DummyVectorIndexed* ptr = nullptr;
DummyVectorIndexed* nullPointer = nullptr;
START_SECTION((UniqueIdIndexer()))
{
  ptr = new DummyVectorIndexed();
	TEST_NOT_EQUAL(ptr, nullPointer)
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
  Math::RandomShuffler r{0};
  r.portable_random_shuffle(dvi.begin(), dvi.end());

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
  r.portable_random_shuffle(dvi.begin(), dvi.end());

  TEST_EXCEPTION_WITH_MESSAGE(Exception::Postcondition,dvi.updateUniqueIdToIndex(),"Duplicate valid unique ids detected!   RandomAccessContainer has size()==12, num_valid_unique_id==10, uniqueid_to_index_.size()==9");

}
END_SECTION

START_SECTION((void updateUniqueIdToIndex() const))
{
  // see uniqueIdToIndex()
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((Size resolveUniqueIdConflicts()))
  DummyVectorIndexed dvi;
  Size num_uii = 10;
  dvi.resize(num_uii);
  for ( Size i = 0; i < num_uii; ++i )
  {
    dvi[i].dummy = i;
    dvi[i].setUniqueId(10*i+1000);
  }
  
  TEST_EQUAL(dvi.resolveUniqueIdConflicts(), 0);
  
  // introduce three doubles
  Dummy a,b;
  a.setUniqueId(1000);
  b.setUniqueId(1000+30);
  dvi.push_back(a);
  dvi.push_back(b);
	dvi.push_back(b);
	TEST_EXCEPTION(Exception::Postcondition,dvi.updateUniqueIdToIndex());
	TEST_EQUAL(dvi.resolveUniqueIdConflicts(), 3);

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
