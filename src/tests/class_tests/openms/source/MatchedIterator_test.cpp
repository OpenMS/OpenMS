// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/DATASTRUCTURES/MatchedIterator.h>
///////////////////////////
#include <OpenMS/KERNEL/MSSpectrum.h>

#include <vector>
#include <ostream>

using namespace OpenMS;
using namespace std;

using MIV = MatchedIterator<vector<double>, ValueTrait, true>;

std::ostream& operator<<(std::ostream& os, const MIV& m)
{
  os << (*m) << "\n";
  return os;
}

START_TEST(MatchedIterator, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MIV* ptr = nullptr;
MIV* null_ptr = nullptr;

START_SECTION(MatchedIterator())
{
  ptr = new MIV();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MatchedIterator())
{
  delete ptr;
}
END_SECTION

vector<double> ref = { 0, 1, 2, 3, 4, 5, 6, 7, 13 };
vector<double> target = { -0.01, 2.5, 3.5, 7, 11 };
vector<double> empty;

START_SECTION((explicit MatchedIterator(const CONT& ref, const CONT& target, float tolerance)))
{
  { // empty reference container
    MIV mi(empty, target, 0.001);
    TEST_EQUAL(mi == mi.end(), true);
  }

  { // empty target container
    MIV mi(ref, empty, 0.001);
    TEST_EQUAL(mi == mi.end(), true);
  }

  { // actual data
    MIV mi(ref, target, 0.001);
    // only a single hit
    TEST_EQUAL(mi.ref(), 7);
    TEST_EQUAL(*mi, 7);
    TEST_EQUAL(mi.refIdx(), 7)
    TEST_EQUAL(mi.tgtIdx(), 3)
    ++mi; // advance to end
    TEST_EQUAL(mi == mi.end(), true);
  }

  { // actual data
    MIV mi_src;
    MIV mi(ref, target, 0.5);
    TEST_EQUAL(mi.ref(), 0);
    TEST_EQUAL(*mi, -0.01);
    TEST_EQUAL(mi.refIdx(), 0)
    TEST_EQUAL(mi.tgtIdx(), 0)
    ++mi;

    TEST_EQUAL(mi.ref(), 2);
    TEST_EQUAL(*mi, 2.5);
    TEST_EQUAL(mi.refIdx(), 2)
    TEST_EQUAL(mi.tgtIdx(), 1)

    mi_src = mi++; // throw in some post-increment
    TEST_EQUAL(mi.ref(), 3);
    TEST_EQUAL(*mi, 2.5);
    TEST_EQUAL(mi.refIdx(), 3)
    TEST_EQUAL(mi.tgtIdx(), 1)

    mi_src = mi++; // throw in some post-increment
    TEST_EQUAL(mi.ref(), 4);
    TEST_EQUAL(*mi, 3.5);
    TEST_EQUAL(mi.refIdx(), 4)
    TEST_EQUAL(mi.tgtIdx(), 2)
    *mi_src; // just to use it once;

    ++mi;    
    TEST_EQUAL(mi.ref(), 7);
    TEST_EQUAL(*mi, 7);
    TEST_EQUAL(mi.refIdx(), 7)
    TEST_EQUAL(mi.tgtIdx(), 3)

    ++mi;
    TEST_EQUAL(mi == mi.end(), true);
    TEST_EQUAL(mi.refIdx(), 9) // points to ref.end()
    TEST_EQUAL(mi.tgtIdx(), 4) // points to last element of target

  }
  // test ppm
  MSSpectrum s, s2;
  s.emplace_back(90.0, 0.0);
  s.emplace_back(100.0, 0.0);
  s.emplace_back(200.0, 0.0);
  s.emplace_back(300.0, 0.0);
  s.emplace_back(400.0, 0.0);
  s2 = s;
  // add a constant to all peaks
  for (auto& p : s2) p.setMZ(p.getMZ() + Math::ppmToMass(2.5, 250.0));
  MatchedIterator<MSSpectrum, PpmTrait, true> it(s, s2, 2.5);
  // the first 3 peaks (90, 100, 200) should not match (since 2.5ppm is smaller than the offset we added)
  TEST_EQUAL(it.refIdx(), 3) // match 300.x
  TEST_EQUAL(it.tgtIdx(), 3)
  ++it;
  TEST_EQUAL(it.refIdx(), 4) // match 400.x
  TEST_EQUAL(it.tgtIdx(), 4)
  ++it;
  TEST_EQUAL(it.refIdx(), 5) // end
  TEST_EQUAL(it == it.end(), true)
  TEST_EQUAL(it.tgtIdx(), 4) // target is always valid

}
END_SECTION

START_SECTION(explicit MatchedIterator())
{
  MIV it;
  TEST_EQUAL(it != it.end(), true)
  it = MIV(empty, empty, 1.0); // assigment is the only valid thing...
  TEST_EQUAL(it == it.end(), true)
}
END_SECTION

START_SECTION(bool operator==(const MatchedIterator& rhs) const)
{
  MIV mi(ref, target, 0.5);
  TEST_EQUAL(mi.ref(), 0);
  TEST_EQUAL(*mi, -0.01);
  TEST_EQUAL(mi.refIdx(), 0)
  ++mi;
  MIV mi2(mi);
  TEST_TRUE(mi == mi2)
  
  TEST_EQUAL(mi.ref(), mi2.ref());
  TEST_EQUAL(*mi, *mi2);
  TEST_EQUAL(mi.refIdx(), mi2.refIdx())
  TEST_EQUAL(mi.ref(), 2);
  TEST_EQUAL(*mi, 2.5);
  TEST_EQUAL(mi.refIdx(), 2)
}
END_SECTION

START_SECTION(bool operator!=(const MatchedIterator& rhs) const)
{
  MIV mi(ref, target, 0.5);
  MIV mi2(mi);
  TEST_EQUAL(mi != mi2, false)
  ++mi;
  TEST_FALSE(mi == mi2)
  MIV mi3(mi);
  TEST_EQUAL(mi != mi3, false)
}
END_SECTION

START_SECTION(const value_type& operator*() const)
{
  NOT_TESTABLE // tested above
}
END_SECTION


typedef pair<double, bool> PDB;

/// Trait for MatchedIterator to find pairs with a certain Th/Da distance in m/z
struct PairTrait
{
  static float allowedTol(float tol, const PDB& /*mz_ref*/)
  {
    return tol;
  }
  /// just use fabs on the value directly
  static float getDiffAbsolute(const PDB& elem_ref, const PDB& elem_tgt)
  {
    return fabs(elem_ref.first - elem_tgt.first);
  }
};

START_SECTION(const value_type& operator->() const)
{
  vector<PDB> ref = { {1, true}, {2, false} };
  
  MatchedIterator<vector<pair<double, bool>>, PairTrait> it(ref, ref, 0.1);
  ++it;
  TEST_EQUAL(it->first, 2.0)

}
END_SECTION

START_SECTION(const value_type& ref() const)
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION(size_t refIdx() const)
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION(size_t tgtIdx() const)
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION(MatchedIterator& operator++())
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION(MatchedIterator operator++(int) const)
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION(static MatchedIterator end())
{
  NOT_TESTABLE // tested above
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

