// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------


#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/KERNEL/Mobilogram.h>
///////////////////////////

#include <sstream>

using namespace OpenMS;
using namespace std;

static_assert(OpenMS::Test::fulfills_rule_of_5<Mobilogram>(), "Must fulfill rule of 5");
static_assert(OpenMS::Test::fulfills_rule_of_6<Mobilogram>(), "Must fulfill rule of 6");
static_assert(OpenMS::Test::fulfills_fast_vector<Mobilogram>(), "Must have fast vector semantics");
static_assert(std::is_nothrow_move_constructible_v<Mobilogram>, "Must have nothrow move constructible");

START_TEST(Mobilogram, "$Id$")

/////////////////////////////////////////////////////////////
// Dummy peak data

MobilityPeak1D p1;
p1.setIntensity(1.0f);
p1.setMobility(2.0);

MobilityPeak1D p2;
p2.setIntensity(2.0f);
p2.setMobility(10.0);

MobilityPeak1D p3;
p3.setIntensity(3.0f);
p3.setMobility(30.0);

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Mobilogram* ptr = nullptr;
Mobilogram* nullPointer = nullptr;
START_SECTION((Mobilogram()))
{
  ptr = new Mobilogram();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~Mobilogram()))
{
  delete ptr;
}
END_SECTION

START_SECTION(([EXTRA] Mobilogram()))
{
  Mobilogram tmp;
  MobilityPeak1D peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);
  TEST_EQUAL(tmp.size(), 1);
  TEST_REAL_SIMILAR(tmp[0].getMobility(), 47.11);
}
END_SECTION

/////////////////////////////////////////////////////////////
// Member accessors

START_SECTION((double getRT() const))
{
  Mobilogram s;
  TEST_REAL_SIMILAR(s.getRT(), -1.0)
}
END_SECTION

START_SECTION((void setRT(double rt)))
{
  Mobilogram s;
  s.setRT(0.451);
  TEST_REAL_SIMILAR(s.getRT(), 0.451)
}
END_SECTION


START_SECTION((double getDriftTimeUnit() const))
{
  Mobilogram s;
  TEST_EQUAL(s.getDriftTimeUnit() == DriftTimeUnit::NONE, true);
}
END_SECTION

START_SECTION((double getDriftTimeUnitAsString() const))
{
  Mobilogram s;
  TEST_EQUAL(s.getDriftTimeUnitAsString(), "<NONE>");
}
END_SECTION

START_SECTION((void setDriftTimeUnit(double dt)))
{
  Mobilogram s;
  s.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  TEST_EQUAL(s.getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true);
  TEST_EQUAL(s.getDriftTimeUnitAsString(), "ms");
}
END_SECTION


/////////////////////////////////////////////////////////////
// RangeManager

START_SECTION((virtual void updateRanges()))
{
  Mobilogram s;
  s.push_back(p1);
  s.push_back(p2);
  s.push_back(p1);

  s.updateRanges();
  s.updateRanges(); // second time to check the initialization

  TEST_REAL_SIMILAR(s.getMaxIntensity(), 2)
  TEST_REAL_SIMILAR(s.getMinIntensity(), 1)
  TEST_REAL_SIMILAR(s.getMaxMobility(), 10)
  TEST_REAL_SIMILAR(s.getMinMobility(), 2)

  // test with only one peak

  s.clear();
  s.push_back(p1);
  s.updateRanges();
  TEST_REAL_SIMILAR(s.getMaxIntensity(), 1)
  TEST_REAL_SIMILAR(s.getMinIntensity(), 1)
  TEST_REAL_SIMILAR(s.getMaxMobility(), 2)
  TEST_REAL_SIMILAR(s.getMinMobility(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Copy constructor, move constructor, assignment operator, move assignment operator, equality

START_SECTION((Mobilogram(const Mobilogram& source)))
{
  Mobilogram tmp;
  tmp.setRT(7.0);
  tmp.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  // peaks
  Mobilogram::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);

  Mobilogram tmp2(tmp);
  TEST_REAL_SIMILAR(tmp2.getRT(), 7.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true)
  // peaks
  TEST_EQUAL(tmp2.size(), 1)
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0], 47.11)
}
END_SECTION

START_SECTION((Mobilogram(const Mobilogram&& source)))
{
  // Ensure that Mobilogram has a no-except move constructor (otherwise
  // std::vector<Mobilogram> is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(Mobilogram(std::declval<Mobilogram&&>())), true)

  Mobilogram tmp;
  tmp.setRT(9.0);
  tmp.setDriftTimeUnit(DriftTimeUnit::VSSC);
  // peaks
  Mobilogram::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);
  peak.getPosition()[0] = 48.11;
  tmp.push_back(peak);

  // copy tmp so we can move one of them
  Mobilogram orig = tmp;
  Mobilogram tmp2(std::move(tmp));

  TEST_EQUAL(tmp2, orig) // should be equal to the original

  TEST_REAL_SIMILAR(tmp2.getRT(), 9.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::VSSC, true)
  TEST_EQUAL(tmp2.size(), 2)
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0], 47.11)
  TEST_REAL_SIMILAR(tmp2[1].getPosition()[0], 48.11)

  // test move -- if this fails, then the move-operator did a copy, not a move... so this better not fail
  TEST_EQUAL(tmp.size(), 0)
}
END_SECTION

START_SECTION((Mobilogram & operator=(const Mobilogram& source)))
{
  Mobilogram tmp;
  tmp.setRT(7.0);
  tmp.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  // peaks
  Mobilogram::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);

  // normal assignment
  Mobilogram tmp2;
  tmp2 = tmp;
  TEST_REAL_SIMILAR(tmp2.getRT(), 7.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::MILLISECOND, true)
  TEST_EQUAL(tmp2.size(), 1);
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0], 47.11);

  // Assignment of empty object
  // normal assignment
  tmp2 = Mobilogram();
  TEST_REAL_SIMILAR(tmp2.getRT(), -1.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::NONE, true)
  TEST_EQUAL(tmp2.size(), 0)
}
END_SECTION

START_SECTION((Mobilogram & operator=(const Mobilogram&& source)))
{
  Mobilogram tmp;
  tmp.setRT(9.0);
  tmp.setDriftTimeUnit(DriftTimeUnit::VSSC);
  // peaks
  Mobilogram::PeakType peak;
  peak.getPosition()[0] = 47.11;
  tmp.push_back(peak);
  peak.getPosition()[0] = 48.11;
  tmp.push_back(peak);

  // copy tmp so we can move one of them
  Mobilogram orig = tmp;

  // move assignment
  Mobilogram tmp2;
  tmp2 = std::move(tmp);

  TEST_EQUAL(tmp2, orig) // should be equal to the original

  TEST_REAL_SIMILAR(tmp2.getRT(), 9.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::VSSC, true)
  TEST_EQUAL(tmp2.size(), 2)
  TEST_REAL_SIMILAR(tmp2[0].getPosition()[0], 47.11)
  TEST_REAL_SIMILAR(tmp2[1].getPosition()[0], 48.11)

  // test move -- if this fails, then the move-operator did a copy, not a move... so this better not fail
  TEST_EQUAL(tmp.size(), 0)

  // Assignment of empty object
  // normal assignment
#ifndef OPENMS_WINDOWSPLATFORM
  #pragma clang diagnostic push
  // Ignore -Wpessimizing-move, because we want to test the move assignment operator.
  #pragma clang diagnostic ignored "-Wpessimizing-move"
#endif
  tmp2 = std::move(Mobilogram());
#ifndef OPENMS_WINDOWSPLATFORM
  #pragma clang diagnostic pop
#endif
  TEST_REAL_SIMILAR(tmp2.getRT(), -1.0)
  TEST_EQUAL(tmp2.getDriftTimeUnit() == DriftTimeUnit::NONE, true)
  TEST_EQUAL(tmp2.size(), 0)
}
END_SECTION

START_SECTION((bool operator==(const Mobilogram& rhs) const))
{
  Mobilogram edit, empty;

  TEST_TRUE(edit == empty);

  edit = empty;
  edit.resize(1);
  TEST_EQUAL(edit == empty, false);

  edit = empty;
  edit.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  TEST_EQUAL(empty == edit, false);

  edit = empty;
  edit.setRT(5);
  TEST_EQUAL(empty == edit, false);

  edit = empty;
  edit.push_back(p1);
  edit.push_back(p2);
  edit.updateRanges();
  edit.clear();
  TEST_TRUE(empty == edit);
}
END_SECTION

START_SECTION((bool operator!=(const Mobilogram& rhs) const))
{
  Mobilogram edit, empty;

  TEST_EQUAL(edit != empty, false);

  edit = empty;
  edit.resize(1);
  TEST_FALSE(edit == empty);

  edit = empty;
  edit.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  TEST_FALSE(edit == empty);

  edit = empty;
  edit.setRT(5);
  TEST_FALSE(edit == empty);

  edit = empty;
  edit.push_back(p1);
  edit.push_back(p2);
  edit.updateRanges();
  edit.clear();
  TEST_TRUE(edit == empty);
}
END_SECTION


/////////////////////////////////////////////////////////////
// Sorting

START_SECTION((void sortByIntensity(bool reverse = false)))
{
  Mobilogram ds;
  MobilityPeak1D p;
  std::vector<double> mzs, intensities;
  intensities.push_back(201);
  mzs.push_back(420.130);
  intensities.push_back(60);
  mzs.push_back(412.824);
  intensities.push_back(56);
  mzs.push_back(423.269);
  intensities.push_back(37);
  mzs.push_back(415.287);
  intensities.push_back(34);
  mzs.push_back(413.800);
  intensities.push_back(31);
  mzs.push_back(419.113);
  intensities.push_back(31);
  mzs.push_back(416.293);
  intensities.push_back(31);
  mzs.push_back(418.232);
  intensities.push_back(29);
  mzs.push_back(414.301);
  intensities.push_back(29);
  mzs.push_back(412.321);

  for (Size i = 0; i < mzs.size(); ++i)
  {
    p.setIntensity(intensities[i]);
    p.setMobility(mzs[i]);
    ds.push_back(p);
  }
  ds.sortByIntensity();
  std::vector<double> intensities_copy(intensities);
  std::sort(intensities_copy.begin(), intensities_copy.end());
  Mobilogram::iterator it_ds = ds.begin();
  ABORT_IF(ds.size() != intensities_copy.size())
  for (std::vector<double>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
  {
    TEST_EQUAL(it_ds->getIntensity(), *it);
    ++it_ds;
  }
  ds.clear();
  for (Size i = 0; i < mzs.size(); ++i)
  {
    p.setIntensity(intensities[i]);
    p.setMobility(mzs[i]);
    ds.push_back(p);
  }

  ds.sortByIntensity();

  Mobilogram::iterator it1 = ds.begin();
  TOLERANCE_ABSOLUTE(0.0001)
  for (std::vector<double>::iterator it = intensities_copy.begin(); it != intensities_copy.end(); ++it)
  {
    if (it1 != ds.end())
    {
      // metadataarray values == mz values
      TEST_REAL_SIMILAR(it1->getIntensity(), *it);
      ++it1;
    }
    else
    {
      TEST_EQUAL(true, false)
    }
  }
}
END_SECTION

START_SECTION((void sortByPosition()))
{
  Mobilogram ds;
  std::vector<double> mzs {423.269, 420.130, 419.113, 418.232, 416.293, 415.287, 414.301, 413.800, 412.824, 412.321};
  std::vector<double> intensities {56, 201, 31, 31, 31, 37, 29, 34, 60, 29};

  for (Size i = 0; i < mzs.size(); ++i)
  {
    ds.emplace_back(mzs[i], intensities[i]);
  }
  ds.sortByPosition();
  Mobilogram::iterator it = ds.begin();
  for (std::vector<double>::reverse_iterator rit = intensities.rbegin(); rit != intensities.rend(); ++rit)
  {
    if (it == ds.end())
    {
      TEST_EQUAL(true, false)
    }
    TEST_EQUAL(it->getIntensity(), *rit);
    ++it;
  }
  ds.clear();
  for (Size i = 0; i < mzs.size(); ++i)
  {
    ds.emplace_back(mzs[i], intensities[i]);
  }
  ds.sortByPosition();

  Size size = intensities.size();
  ABORT_IF(ds.size() != size);
  Mobilogram::iterator it1 = ds.begin();
  for (std::vector<double>::reverse_iterator rit = intensities.rbegin(); rit != intensities.rend(); ++rit)
  {
    // metadataarray values == intensity values
    TEST_REAL_SIMILAR(it1->getIntensity(), *rit);
    ++it1;
  }
}
END_SECTION

START_SECTION(([EXTRA] sorting reorders all data arrays alongside the peaks))
{
  // Regression test: sortByPosition()/sortByIntensity() used to be a bare stable_sort
  // over the peaks and never touched the float/string/integer data arrays, so any
  // per-peak metadata was silently mis-associated after a sort.

  // peaks in mobility order 3,1,2 with intensities 10,30,20; arrays mirror their peak
  auto make = []() {
    Mobilogram m;
    m.emplace_back(3.0, 10.0f);
    m.emplace_back(1.0, 30.0f);
    m.emplace_back(2.0, 20.0f);

    Mobilogram::FloatDataArray fda;
    fda.push_back(3.5f); fda.push_back(1.5f); fda.push_back(2.5f);
    fda.setName("f1");
    m.getFloatDataArrays().push_back(fda);

    Mobilogram::StringDataArray sda;
    sda.push_back("mb3"); sda.push_back("mb1"); sda.push_back("mb2");
    sda.setName("s1");
    m.getStringDataArrays().push_back(sda);

    Mobilogram::IntegerDataArray ida;
    ida.push_back(3); ida.push_back(1); ida.push_back(2);
    ida.setName("i1");
    m.getIntegerDataArrays().push_back(ida);
    return m;
  };

  // --- sortByPosition: mobility 1,2,3
  Mobilogram m = make();
  m.sortByPosition();
  TEST_REAL_SIMILAR(m[0].getMobility(), 1.0)
  TEST_REAL_SIMILAR(m[1].getMobility(), 2.0)
  TEST_REAL_SIMILAR(m[2].getMobility(), 3.0)
  TEST_REAL_SIMILAR(m.getFloatDataArrays()[0][0], 1.5)
  TEST_REAL_SIMILAR(m.getFloatDataArrays()[0][1], 2.5)
  TEST_REAL_SIMILAR(m.getFloatDataArrays()[0][2], 3.5)
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][0], "mb1")
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][1], "mb2")
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][2], "mb3")
  TEST_EQUAL(m.getIntegerDataArrays()[0][0], 1)
  TEST_EQUAL(m.getIntegerDataArrays()[0][1], 2)
  TEST_EQUAL(m.getIntegerDataArrays()[0][2], 3)
  // array names must survive the permutation
  TEST_STRING_EQUAL(m.getFloatDataArrays()[0].getName(), "f1")
  TEST_STRING_EQUAL(m.getStringDataArrays()[0].getName(), "s1")
  TEST_STRING_EQUAL(m.getIntegerDataArrays()[0].getName(), "i1")

  // --- sortByIntensity ascending: 10,20,30 -> mobility 3,2,1
  m = make();
  m.sortByIntensity();
  TEST_REAL_SIMILAR(m[0].getIntensity(), 10.0)
  TEST_REAL_SIMILAR(m[2].getIntensity(), 30.0)
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][0], "mb3")
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][2], "mb1")
  TEST_EQUAL(m.getIntegerDataArrays()[0][0], 3)
  TEST_EQUAL(m.getIntegerDataArrays()[0][2], 1)

  // --- sortByIntensity reverse: 30,20,10 -> mobility 1,2,3
  m = make();
  m.sortByIntensity(true);
  TEST_REAL_SIMILAR(m[0].getIntensity(), 30.0)
  TEST_REAL_SIMILAR(m[2].getIntensity(), 10.0)
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][0], "mb1")
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][2], "mb3")
  TEST_EQUAL(m.getIntegerDataArrays()[0][0], 1)
  TEST_EQUAL(m.getIntegerDataArrays()[0][2], 3)

  // --- a mis-sized data array must be rejected rather than silently corrupting
  m = make();
  m.getIntegerDataArrays()[0].push_back(99); // now 4 entries for 3 peaks
  TEST_EXCEPTION(Exception::Precondition, m.sortByPosition())

  // --- stable order across *tied* keys: the data arrays must stay aligned with the peaks
  //     through the tie, and equal keys must keep their original relative order. Four peaks
  //     with tied mobilities and tied intensities; string/integer arrays carry each peak's
  //     original identity (A..D / 1..4). The starting order is deliberately unsorted by
  //     mobility, by ascending intensity and by descending intensity, so none of the three
  //     sorts can short-circuit -- each really permutes.
  auto make_ties = []() {
    Mobilogram m;
    m.emplace_back(1.0, 10.0f); // A / 1
    m.emplace_back(2.0, 20.0f); // B / 2
    m.emplace_back(1.0, 20.0f); // C / 3
    m.emplace_back(2.0, 10.0f); // D / 4

    Mobilogram::StringDataArray sda;
    sda.push_back("A"); sda.push_back("B"); sda.push_back("C"); sda.push_back("D");
    m.getStringDataArrays().push_back(sda);

    Mobilogram::IntegerDataArray ida;
    ida.push_back(1); ida.push_back(2); ida.push_back(3); ida.push_back(4);
    m.getIntegerDataArrays().push_back(ida);
    return m;
  };

  auto check_order = [](const Mobilogram& m, const std::vector<std::string>& s, const std::vector<Int>& i) {
    for (Size k = 0; k < s.size(); ++k)
    {
      TEST_STRING_EQUAL(m.getStringDataArrays()[0][k], s[k])
      TEST_EQUAL(m.getIntegerDataArrays()[0][k], i[k])
    }
  };

  // sortByPosition: mobility ascending; ties (mob=1: A,C) and (mob=2: B,D) keep original order
  Mobilogram mt = make_ties();
  mt.sortByPosition();
  TEST_REAL_SIMILAR(mt[0].getMobility(), 1.0) TEST_REAL_SIMILAR(mt[1].getMobility(), 1.0)
  TEST_REAL_SIMILAR(mt[2].getMobility(), 2.0) TEST_REAL_SIMILAR(mt[3].getMobility(), 2.0)
  check_order(mt, {"A", "C", "B", "D"}, {1, 3, 2, 4});

  // sortByIntensity ascending; ties (I=10: A,D) then (I=20: B,C) keep original order
  mt = make_ties();
  mt.sortByIntensity();
  check_order(mt, {"A", "D", "B", "C"}, {1, 4, 2, 3});

  // sortByIntensity reverse (descending); ties (I=20: B,C) then (I=10: A,D) keep original order
  mt = make_ties();
  mt.sortByIntensity(true);
  check_order(mt, {"B", "C", "A", "D"}, {2, 3, 1, 4});
}
END_SECTION

START_SECTION(bool isSorted() const)
{
  // make test dataset
  Mobilogram spec;
  MobilityPeak1D p;
  p.setIntensity(1.0);
  p.setMobility(1000.0);
  spec.push_back(p);

  p.setIntensity(1.0);
  p.setMobility(1001.0);
  spec.push_back(p);

  p.setIntensity(1.0);
  p.setMobility(1002.0);
  spec.push_back(p);

  TEST_EQUAL(spec.isSorted(), true)

  reverse(spec.begin(), spec.end());
  TEST_EQUAL(spec.isSorted(), false)
}
END_SECTION

START_SECTION(template<class Predicate> bool isSorted(const Predicate& lamdba) const)
{
  Mobilogram ds;
  std::vector<double> mzs {423.269, 420.130, 419.113, 418.232, 416.293, 415.287, 414.301, 413.800, 412.824, 412.321};
  std::vector<double> intensities {56, 201, 31, 31, 31, 37, 29, 34, 60, 29};

  for (Size i = 0; i < mzs.size(); ++i)
  {
    ds.emplace_back(mzs[i], intensities[i]);
  }
  ds.sortByPosition();

  // more expensive than isSorted(), but just to make sure
  TEST_EQUAL(ds.isSorted([&ds](Size a, Size b) { return ds[a].getMobility() < ds[b].getMobility(); }), true)
  TEST_EQUAL(ds.isSorted(), true) // call other method. Should give the same result

  ds.sortByIntensity();
  TEST_EQUAL(ds.isSorted([&ds](Size a, Size b) { return ds[a].getIntensity() < ds[b].getIntensity(); }), true)
  TEST_EQUAL(ds.isSorted([&ds](Size a, Size b) { return ds[a].getMobility() < ds[b].getMobility(); }), false)
  TEST_EQUAL(ds.isSorted(), false) // call other method. Should give the same result
}
END_SECTION

START_SECTION(template<class Predicate> void sort(const Predicate& lambda))
{
  // sort by a data array rather than by the peaks themselves; peaks and *all*
  // data arrays must follow the same permutation
  Mobilogram m;
  m.emplace_back(1.0, 10.0f);
  m.emplace_back(2.0, 20.0f);
  m.emplace_back(3.0, 30.0f);

  Mobilogram::IntegerDataArray rank;
  rank.push_back(2); rank.push_back(0); rank.push_back(1);
  rank.setName("rank");
  m.getIntegerDataArrays().push_back(rank);

  Mobilogram::StringDataArray tag;
  tag.push_back("a"); tag.push_back("b"); tag.push_back("c");
  m.getStringDataArrays().push_back(tag);

  const auto& r = m.getIntegerDataArrays()[0];
  m.sort([&r](const Size i1, const Size i2) { return r[i1] < r[i2]; });

  // ascending 'rank' means original order 1,2,0
  TEST_EQUAL(m.getIntegerDataArrays()[0][0], 0)
  TEST_EQUAL(m.getIntegerDataArrays()[0][1], 1)
  TEST_EQUAL(m.getIntegerDataArrays()[0][2], 2)
  TEST_REAL_SIMILAR(m[0].getMobility(), 2.0)
  TEST_REAL_SIMILAR(m[1].getMobility(), 3.0)
  TEST_REAL_SIMILAR(m[2].getMobility(), 1.0)
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][0], "b")
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][1], "c")
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][2], "a")
}
END_SECTION

START_SECTION((Mobilogram& select(const std::vector<Size>& indices)))
{
  Mobilogram m;
  m.emplace_back(1.0, 10.0f);
  m.emplace_back(2.0, 20.0f);
  m.emplace_back(3.0, 30.0f);

  Mobilogram::FloatDataArray fda;
  fda.push_back(1.5f); fda.push_back(2.5f); fda.push_back(3.5f);
  fda.setName("f1");
  m.getFloatDataArrays().push_back(fda);

  Mobilogram::StringDataArray sda;
  sda.push_back("a"); sda.push_back("b"); sda.push_back("c");
  m.getStringDataArrays().push_back(sda);

  Mobilogram::IntegerDataArray ida;
  ida.push_back(1); ida.push_back(2); ida.push_back(3);
  m.getIntegerDataArrays().push_back(ida);

  // subset + reorder in one go
  m.select(std::vector<Size>{2, 0});
  TEST_EQUAL(m.size(), 2)
  TEST_REAL_SIMILAR(m[0].getMobility(), 3.0)
  TEST_REAL_SIMILAR(m[1].getMobility(), 1.0)
  TEST_REAL_SIMILAR(m.getFloatDataArrays()[0][0], 3.5)
  TEST_REAL_SIMILAR(m.getFloatDataArrays()[0][1], 1.5)
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][0], "c")
  TEST_STRING_EQUAL(m.getStringDataArrays()[0][1], "a")
  TEST_EQUAL(m.getIntegerDataArrays()[0][0], 3)
  TEST_EQUAL(m.getIntegerDataArrays()[0][1], 1)
  TEST_STRING_EQUAL(m.getFloatDataArrays()[0].getName(), "f1")

  // a valid (unique) reorder keeps the data arrays aligned
  Mobilogram rev;
  rev.emplace_back(1.0, 10.0f);
  rev.emplace_back(2.0, 20.0f);
  Mobilogram::StringDataArray rev_sda;
  rev_sda.push_back("a"); rev_sda.push_back("b");
  rev.getStringDataArrays().push_back(rev_sda);
  rev.select(std::vector<Size>{1, 0});
  TEST_EQUAL(rev.size(), 2)
  TEST_STRING_EQUAL(rev.getStringDataArrays()[0][0], "b")
  TEST_STRING_EQUAL(rev.getStringDataArrays()[0][1], "a")

  // size mismatch between peaks and a data array is an error
  Mobilogram bad;
  bad.emplace_back(1.0, 10.0f);
  bad.emplace_back(2.0, 20.0f);
  Mobilogram::IntegerDataArray too_long;
  too_long.push_back(1); too_long.push_back(2); too_long.push_back(3);
  bad.getIntegerDataArrays().push_back(too_long);
  TEST_EXCEPTION(Exception::Precondition, bad.select(std::vector<Size>{0, 1}))
  // ... and the rejected call must leave the mobilogram untouched
  TEST_EQUAL(bad.size(), 2)
  TEST_REAL_SIMILAR(bad[0].getMobility(), 1.0)
  TEST_REAL_SIMILAR(bad[1].getMobility(), 2.0)
  TEST_EQUAL(bad.getIntegerDataArrays()[0].size(), 3)

  // out-of-range indices are rejected rather than being undefined behaviour
  Mobilogram oor;
  oor.emplace_back(1.0, 10.0f);
  oor.emplace_back(2.0, 20.0f);
  TEST_EXCEPTION(Exception::Precondition, oor.select(std::vector<Size>{0, 5}))
  TEST_EQUAL(oor.size(), 2)
}
END_SECTION


  /////////////////////////////////////////////////////////////
  // Finding peaks or peak ranges

  const Mobilogram spec_find = []() {
    Mobilogram spec;
    spec.push_back({1.0, 29.0f});
    spec.push_back({2.0, 60.0f});
    spec.push_back({3.0, 34.0f});
    spec.push_back({4.0, 29.0f});
    spec.push_back({5.0, 37.0f});
    spec.push_back({6.0, 31.0f});
    return spec;
  }();

START_SECTION((Iterator MBEnd(CoordinateType mb)))
{
  Mobilogram::Iterator it;
  auto tmp = spec_find;
  it = tmp.MBEnd(4.5);
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = tmp.MBEnd(5.0);
  TEST_EQUAL(it->getPosition()[0], 6.0)
  it = tmp.MBEnd(5.5);
  TEST_EQUAL(it->getPosition()[0], 6.0)
}
END_SECTION

START_SECTION((Iterator MBBegin(CoordinateType mb)))
{
  Mobilogram::Iterator it;
  auto tmp = spec_find;

  it = tmp.MBBegin(4.5);
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = tmp.MBBegin(5.0);
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = tmp.MBBegin(5.5);
  TEST_EQUAL(it->getPosition()[0], 6.0)
}
END_SECTION

START_SECTION((Iterator MBBegin(Iterator begin, CoordinateType mb, Iterator end)))
{
  Mobilogram::Iterator it;
  auto tmp = spec_find;

  it = tmp.MBBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = tmp.MBBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = tmp.MBBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0], tmp.begin()->getPosition()[0])
}
END_SECTION

START_SECTION((ConstIterator MBBegin(ConstIterator begin, CoordinateType mb, ConstIterator end) const))
{
  Mobilogram::Iterator it;
  auto tmp = spec_find;

  it = tmp.MBBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = tmp.MBBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = tmp.MBBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0], tmp.begin()->getPosition()[0])
}
END_SECTION

START_SECTION((Iterator MBEnd(Iterator begin, CoordinateType mb, Iterator end)))
{
  Mobilogram::Iterator it;
  auto tmp = spec_find;

  it = tmp.MBEnd(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = tmp.MBEnd(tmp.begin(), 5, tmp.end());
  TEST_EQUAL(it->getPosition()[0], 6.0)
  it = tmp.MBEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPosition()[0], tmp.begin()->getPosition()[0])
}
END_SECTION

START_SECTION((ConstIterator MBEnd(ConstIterator begin, CoordinateType mb, ConstIterator end) const))
{
  Mobilogram::ConstIterator it;

  it = spec_find.MBEnd(spec_find.begin(), 4.5, spec_find.end());
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = spec_find.MBEnd(spec_find.begin(), 5, spec_find.end());
  TEST_EQUAL(it->getPosition()[0], 6.0)
  it = spec_find.MBEnd(spec_find.begin(), 4.5, spec_find.begin());
  TEST_EQUAL(it->getPosition()[0], spec_find.begin()->getPosition()[0])
}
END_SECTION

START_SECTION((ConstIterator MBEnd(CoordinateType mb) const))
{
  Mobilogram::ConstIterator it;

  it = spec_find.MBEnd(4.5);
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = spec_find.MBEnd(5.0);
  TEST_EQUAL(it->getPosition()[0], 6.0)
  it = spec_find.MBEnd(5.5);
  TEST_EQUAL(it->getPosition()[0], 6.0)
}
END_SECTION

START_SECTION((ConstIterator MBBegin(CoordinateType mb) const))
{
  Mobilogram::ConstIterator it;

  it = spec_find.MBBegin(4.5);
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = spec_find.MBBegin(5.0);
  TEST_EQUAL(it->getPosition()[0], 5.0)
  it = spec_find.MBBegin(5.5);
  TEST_EQUAL(it->getPosition()[0], 6.0)
}
END_SECTION

auto tmp = spec_find;

START_SECTION((Iterator PosBegin(CoordinateType mb)))
{
  Mobilogram::Iterator it;
  it = tmp.PosBegin(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.0);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((Iterator PosBegin(Iterator begin, CoordinateType mb, Iterator end)))
{
  Mobilogram::Iterator it;
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(tmp.begin(), 5.5, tmp.end());
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it - 1)->getPos(), (tmp.end() - 1)->getPos())
}
END_SECTION

START_SECTION((ConstIterator PosBegin(CoordinateType mb) const))
{
  Mobilogram::ConstIterator it;
  it = tmp.PosBegin(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.0);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((ConstIterator PosBegin(ConstIterator begin, CoordinateType mb, ConstIterator end) const))
{
  Mobilogram::ConstIterator it;
  it = tmp.PosBegin(tmp.begin(), 3.5, tmp.end());
  TEST_EQUAL(it->getPos(), 4.0)
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosBegin(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it - 1)->getPos(), (tmp.end() - 1)->getPos())
}
END_SECTION

START_SECTION((Iterator PosEnd(CoordinateType mb)))
{
  Mobilogram::Iterator it;
  it = tmp.PosEnd(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(5.0);
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosEnd(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((Iterator PosEnd(Iterator begin, CoordinateType mb, Iterator end)))
{
  Mobilogram::Iterator it;
  it = tmp.PosEnd(tmp.begin(), 3.5, tmp.end());
  TEST_EQUAL(it->getPos(), 4.0)
  it = tmp.PosEnd(tmp.begin(), 4.0, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it - 1)->getPos(), (tmp.end() - 1)->getPos())
}
END_SECTION

START_SECTION((ConstIterator PosEnd(CoordinateType mb) const))
{
  Mobilogram::ConstIterator it;
  it = tmp.PosEnd(4.5);
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(5.0);
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosEnd(5.5);
  TEST_EQUAL(it->getPos(), 6.0)
}
END_SECTION

START_SECTION((ConstIterator PosEnd(ConstIterator begin, CoordinateType mb, ConstIterator end) const))
{
  Mobilogram::ConstIterator it;
  it = tmp.PosEnd(tmp.begin(), 4.5, tmp.end());
  TEST_EQUAL(it->getPos(), 5.0)
  it = tmp.PosEnd(tmp.begin(), 5.0, tmp.end());
  TEST_EQUAL(it->getPos(), 6.0)
  it = tmp.PosEnd(tmp.begin(), 4.5, tmp.begin());
  TEST_EQUAL(it->getPos(), tmp.begin()->getPos())
  it = tmp.PosBegin(tmp.begin(), 8.0, tmp.end());
  TEST_EQUAL((it - 1)->getPos(), (tmp.end() - 1)->getPos())
}
END_SECTION

const Mobilogram spec_test = []() {
  Mobilogram spec_test;
  spec_test.push_back({412.321, 29.0f});
  spec_test.push_back({412.824, 60.0f});
  spec_test.push_back({413.8, 34.0f});
  spec_test.push_back({414.301, 29.0f});
  spec_test.push_back({415.287, 37.0f});
  spec_test.push_back({416.293, 31.0f});
  spec_test.push_back({418.232, 31.0f});
  spec_test.push_back({419.113, 31.0f});
  spec_test.push_back({420.13, 201.0f});
  spec_test.push_back({423.269, 56.0f});
  spec_test.push_back({426.292, 34.0f});
  spec_test.push_back({427.28, 82.0f});
  spec_test.push_back({428.322, 87.0f});
  spec_test.push_back({430.269, 30.0f});
  spec_test.push_back({431.246, 29.0f});
  spec_test.push_back({432.289, 42.0f});
  spec_test.push_back({436.161, 32.0f});
  spec_test.push_back({437.219, 54.0f});
  spec_test.push_back({439.186, 40.0f});
  spec_test.push_back({440.27, 40});
  spec_test.push_back({441.224, 23.0f});
  return spec_test;
}();

START_SECTION((Size findNearest(CoordinateType mb) const))
{
  Mobilogram tmp = spec_test;

  // test outside mass range
  TEST_EQUAL(tmp.findNearest(400.0), 0);
  TEST_EQUAL(tmp.findNearest(500.0), 20);
  // test mass range borders
  TEST_EQUAL(tmp.findNearest(412.4), 0);
  TEST_EQUAL(tmp.findNearest(441.224), 20);
  // test inside scan
  TEST_EQUAL(tmp.findNearest(426.29), 10);
  TEST_EQUAL(tmp.findNearest(426.3), 10);
  TEST_EQUAL(tmp.findNearest(427.2), 11);
  TEST_EQUAL(tmp.findNearest(427.3), 11);

  // empty spectrum
  Mobilogram tmp2;
  TEST_PRECONDITION_VIOLATED(tmp2.findNearest(427.3));
}
END_SECTION

START_SECTION((Size findNearest(CoordinateType mb, CoordinateType tolerance) const))
{
  // test outside mass range
  TEST_EQUAL(spec_test.findNearest(400.0, 1.0), -1);
  TEST_EQUAL(spec_test.findNearest(500.0, 1.0), -1);

  // test mass range borders
  TEST_EQUAL(spec_test.findNearest(412.4, 0.01), -1);
  TEST_EQUAL(spec_test.findNearest(412.4, 0.1), 0);
  TEST_EQUAL(spec_test.findNearest(441.3, 0.01), -1);
  TEST_EQUAL(spec_test.findNearest(441.3, 0.1), 20);

  // test inside scan
  TEST_EQUAL(spec_test.findNearest(426.29, 0.1), 10);
  TEST_EQUAL(spec_test.findNearest(426.3, 0.1), 10);
  TEST_EQUAL(spec_test.findNearest(427.2, 0.1), 11);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.1), 11);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.001), -1);

  // empty spectrum
  Mobilogram spec_test2;
  TEST_EQUAL(spec_test2.findNearest(427.3, 1.0, 1.0), -1);
}
END_SECTION
START_SECTION((Size findNearest(CoordinateType mb, CoordinateType left_tolerance, CoordinateType right_tolerance) const))
{
  // test outside mass range
  TEST_EQUAL(spec_test.findNearest(400.0, 1.0, 1.0), -1);
  TEST_EQUAL(spec_test.findNearest(500.0, 1.0, 1.0), -1);

  // test mass range borders
  TEST_EQUAL(spec_test.findNearest(412.4, 0.01, 0.01), -1);
  TEST_EQUAL(spec_test.findNearest(412.4, 0.1, 0.1), 0);
  TEST_EQUAL(spec_test.findNearest(441.3, 0.01, 0.01), -1);
  TEST_EQUAL(spec_test.findNearest(441.3, 0.1, 0.1), 20);

  // test inside scan
  TEST_EQUAL(spec_test.findNearest(426.29, 0.1, 0.1), 10);
  TEST_EQUAL(spec_test.findNearest(426.3, 0.1, 0.1), 10);
  TEST_EQUAL(spec_test.findNearest(427.2, 0.1, 0.1), 11);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.1, 0.1), 11);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.001, 0.001), -1);

  TEST_EQUAL(spec_test.findNearest(427.3, 0.1, 0.001), 11);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.001, 1.01), -1);
  TEST_EQUAL(spec_test.findNearest(427.3, 0.001, 1.1), 12);

  // empty spectrum
  Mobilogram spec_test2;
  TEST_EQUAL(spec_test2.findNearest(427.3, 1.0, 1.0), -1);
}
END_SECTION
START_SECTION((Size findHighestInWindow(CoordinateType mb, CoordinateType tolerance_left, CoordinateType tolerance_righ) const))
{
  // test outside mass range
  TEST_EQUAL(spec_test.findHighestInWindow(400.0, 1.0, 1.0), -1);
  TEST_EQUAL(spec_test.findHighestInWindow(500.0, 1.0, 1.0), -1);

  // test mass range borders
  TEST_EQUAL(spec_test.findHighestInWindow(412.4, 0.01, 0.01), -1);
  TEST_EQUAL(spec_test.findHighestInWindow(412.4, 0.1, 0.1), 0);
  TEST_EQUAL(spec_test.findHighestInWindow(441.3, 0.01, 0.01), -1);
  TEST_EQUAL(spec_test.findHighestInWindow(441.3, 0.1, 0.1), 20);

  // test inside scan
  TEST_EQUAL(spec_test.findHighestInWindow(426.29, 0.1, 0.1), 10);
  TEST_EQUAL(spec_test.findHighestInWindow(426.3, 0.1, 0.1), 10);
  TEST_EQUAL(spec_test.findHighestInWindow(427.2, 0.1, 0.1), 11);
  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 0.1, 0.1), 11);
  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 0.001, 0.001), -1);

  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 0.1, 0.001), 11);
  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 0.001, 1.01), -1);
  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 0.001, 1.1), 12);

  TEST_EQUAL(spec_test.findHighestInWindow(427.3, 9.0, 4.0), 8);
  TEST_EQUAL(spec_test.findHighestInWindow(430.25, 1.9, 1.01), 13);

  // empty spectrum
  Mobilogram spec_test2;
  TEST_EQUAL(spec_test2.findHighestInWindow(427.3, 1.0, 1.0), -1);
}
END_SECTION

START_SECTION(ConstIterator getBasePeak() const)
{
  const auto it = spec_test.getBasePeak();
  TEST_REAL_SIMILAR(it->getIntensity(), 201.0)
  TEST_EQUAL(std::distance(spec_test.begin(), it), 8);
  Mobilogram empty;
  TEST_EQUAL(empty.getBasePeak() == empty.end(), true);
}
END_SECTION


START_SECTION(Iterator getBasePeak())
{
  Mobilogram test = spec_test;
  auto it = test.getBasePeak();
  it->setIntensity(it->getIntensity() + 0.0);
  TEST_REAL_SIMILAR(it->getIntensity(), 201.0)
  TEST_EQUAL(std::distance(test.begin(), it), 8);
}
END_SECTION


START_SECTION(PeakType::IntensityType calculateTIC() const)
{
  auto r = spec_test.calculateTIC();
  TEST_REAL_SIMILAR(r, 1032.0)
  TEST_EQUAL(Mobilogram().calculateTIC(), 0.0)
}
END_SECTION


START_SECTION(void clear())
{
  Mobilogram edit;
  edit.resize(1);
  edit.setRT(5);
  edit.setDriftTimeUnit(DriftTimeUnit::MILLISECOND);
  edit.getFloatDataArrays().resize(3);
  edit.getIntegerDataArrays().resize(2);
  edit.getStringDataArrays().resize(1);

  edit.clear();
  TEST_EQUAL(edit.size(), 0)
  TEST_EQUAL(edit == Mobilogram(), false)
  TEST_EQUAL(edit.empty(), true)
  // the data arrays are parallel to the peaks, so dropping the peaks must drop them too
  // (previously they survived clear(), so a subsequent fill appended to stale entries)
  TEST_EQUAL(edit.getFloatDataArrays().empty(), true)
  TEST_EQUAL(edit.getIntegerDataArrays().empty(), true)
  TEST_EQUAL(edit.getStringDataArrays().empty(), true)
}
END_SECTION

START_SECTION(([Mobilogram::RTLess] bool operator()(const Mobilogram& a, const Mobilogram& b) const))
{
  vector<Mobilogram> v;

  Mobilogram sp1;
  sp1.setRT(3.0f);
  v.push_back(sp1);

  Mobilogram sp2;
  sp2.setRT(2.0f);
  v.push_back(sp2);

  Mobilogram sp3;
  sp3.setRT(1.0f);
  v.push_back(sp3);

  std::sort(v.begin(), v.end(), Mobilogram::RTLess());

  TEST_REAL_SIMILAR(v[0].getRT(), 1.0);
  TEST_REAL_SIMILAR(v[1].getRT(), 2.0);
  TEST_REAL_SIMILAR(v[2].getRT(), 3.0);

  ///
  Mobilogram s1;
  s1.setRT(0.451);

  Mobilogram s2;
  s2.setRT(0.5);

  TEST_EQUAL(Mobilogram::RTLess()(s1, s2), true);
  TEST_EQUAL(Mobilogram::RTLess()(s2, s1), false);
  TEST_EQUAL(Mobilogram::RTLess()(s2, s2), false);
}
END_SECTION

START_SECTION(([EXTRA] std::ostream & operator<<(std::ostream& os, const Mobilogram& spec)))
{
  Mobilogram spec;
  MobilityPeak1D p;
  p.setIntensity(29.0f);
  p.setMobility(412.321);
  spec.push_back(p); // 0
  p.setIntensity(60.0f);
  p.setMobility(412.824);
  spec.push_back(p); // 1
  p.setIntensity(34.0f);
  p.setMobility(413.8);
  spec.push_back(p); // 2
  p.setIntensity(29.0f);
  p.setMobility(414.301);
  spec.push_back(p); // 3
  p.setIntensity(37.0f);
  p.setMobility(415.287);
  spec.push_back(p); // 4
  p.setIntensity(31.0f);
  p.setMobility(416.293);
  spec.push_back(p); // 5
  p.setIntensity(31.0f);
  p.setMobility(418.232);
  spec.push_back(p); // 6
  p.setIntensity(31.0f);
  p.setMobility(419.113);
  spec.push_back(p); // 7
  p.setIntensity(201.0f);
  p.setMobility(420.13);
  spec.push_back(p); // 8
  p.setIntensity(56.0f);
  p.setMobility(423.269);
  spec.push_back(p); // 9
  p.setIntensity(34.0f);
  p.setMobility(426.292);
  spec.push_back(p); // 10

  spec.setRT(7.0);

  ostringstream test_stream;
  test_stream << spec;

  TEST_EQUAL(test_stream.str(), "-- MOBILOGRAM BEGIN --\n"
                                "POS: 412.321 INT: 29\n"
                                "POS: 412.824 INT: 60\n"
                                "POS: 413.8 INT: 34\n"
                                "POS: 414.301 INT: 29\n"
                                "POS: 415.287 INT: 37\n"
                                "POS: 416.293 INT: 31\n"
                                "POS: 418.232 INT: 31\n"
                                "POS: 419.113 INT: 31\n"
                                "POS: 420.13 INT: 201\n"
                                "POS: 423.269 INT: 56\n"
                                "POS: 426.292 INT: 34\n"
                                "-- MOBILOGRAM END --\n")
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
