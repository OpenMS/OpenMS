// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/RangeManager.h>

///////////////////////////

#include <sstream>

using namespace OpenMS;
using namespace std;

// test with additional Mobility (should always be empty)
using RangeMType = RangeManagerContainer<RangeRT, RangeMZ, RangeIntensity, RangeMobility>;
using RangeMTypeInt = RangeManager<RangeIntensity>;
using RangeMTypeMzInt = RangeManager<RangeMZ, RangeIntensity>;
using RangeMTypeRT = RangeManager<RangeRT>;

class RM : public RangeMType
{
  public:  
    // avoid compiler warning, but in production code virtual should not be done on a RM to save
    // space and time
    virtual ~RM() = default;

    bool operator == (const RM& rhs) const
    {
      return RangeMType::operator==(rhs);
    }

    bool operator != (const RM& rhs) const
    {
      return !(operator==(rhs));
    }

    void updateRanges() override
    {
      std::vector<Peak2D > vec;
      Peak2D tmp;

      tmp.getPosition()[0] = 2.0;
      tmp.getPosition()[1] = 500.0;
      tmp.setIntensity(1.0f);
      vec.push_back(tmp);

      tmp.getPosition()[0] = 100.0;
      tmp.getPosition()[1] = 1300.0;
      tmp.setIntensity(47110.0);
      vec.push_back(tmp);

      tmp.getPosition()[0] = 2.0;
      tmp.getPosition()[1] = 500.0;
      tmp.setIntensity(1.0f);
      vec.push_back(tmp);

      clearRanges();
      for (const auto& peak : vec)
      {
        extendRT(peak.getRT());
        extendMZ(peak.getMZ());
        extendIntensity(peak.getIntensity());
      }
    }

    virtual void updateRanges2()
    {
      std::vector<Peak2D > vec;
      Peak2D tmp;

      tmp.getPosition()[0] = 2.0;
      tmp.getPosition()[1] = 500.0;
      tmp.setIntensity(1.0f);
      vec.push_back(tmp);

      clearRanges();
      for (const auto& peak : vec)
      {
        extendRT(peak.getRT());
        extendMZ(peak.getMZ());
        extendIntensity(peak.getIntensity());
      }
    }
}; // class RM

START_TEST(RangeManager, "RangeManager")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

// tests for RangeBase

START_SECTION(RangeBase())
  RangeBase b;
  TEST_EQUAL(b.isEmpty(), true)
END_SECTION

START_SECTION(RangeBase(const double min, const double max))
  RangeBase b(4, 6);
  TEST_EQUAL(b.isEmpty(), false)
  TEST_EQUAL(b.getMin(), 4)
  TEST_EQUAL(b.getMax(), 6)

  TEST_EXCEPTION(Exception::InvalidRange, RangeBase(6, 3))
END_SECTION

START_SECTION(const RangeBase& rhs)
  RangeBase b_(4, 6);
  auto b(b_);
  TEST_EQUAL(b.isEmpty(), false)
  TEST_EQUAL(b.getMin(), 4)
  TEST_EQUAL(b.getMax(), 6)
END_SECTION

START_SECTION(RangeBase& operator=(const RangeBase& rhs))
  RangeBase b_(4, 6);
  RangeBase b;
  b = b_;
  TEST_EQUAL(b.isEmpty(), false)
  TEST_EQUAL(b.getMin(), 4)
  TEST_EQUAL(b.getMax(), 6)
END_SECTION

START_SECTION(void clear())
  RangeBase b(4, 6);
  TEST_EQUAL(b.isEmpty(), false)
  b.clear();
  TEST_EQUAL(b.isEmpty(), true)
END_SECTION

START_SECTION(bool isEmpty() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(bool contains(const double value) const)
  RangeBase b(4, 6);
  TEST_EQUAL(b.contains(5), true)
  TEST_EQUAL(b.contains(3), false)
  TEST_EQUAL(b.contains(7), false)
  RangeBase empty;
  TEST_EQUAL(empty.contains(5), false)
END_SECTION

START_SECTION(bool contains(const RangeBase& inner_range) const)
  RangeBase b(2, 6), inner1(2,4), inner2(3,4), inner3(4,6), over1(1,4), over2(3,7), outer(1,7);
  TEST_EQUAL(b.contains(inner1), true)
  TEST_EQUAL(b.contains(inner2), true)
  TEST_EQUAL(b.contains(inner3), true)
  TEST_EQUAL(b.contains(over1), false)
  TEST_EQUAL(b.contains(over2), false)
  TEST_EQUAL(b.contains(outer), false)
  TEST_EQUAL(outer.contains(b), true)
END_SECTION

START_SECTION(void setMin(const double min))
  RangeBase b(4, 6);
  b.setMin(5);
  TEST_EQUAL(b.getMin(), 5)
  b.setMin(7);// also increases max
  TEST_EQUAL(b.getMin(), 7)
  TEST_EQUAL(b.getMax(), 7);
END_SECTION

START_SECTION(void setMax(const double max))
  RangeBase b(4, 6);
  b.setMax(5);
  TEST_EQUAL(b.getMax(), 5)
  b.setMax(2);// also decreases min
  TEST_EQUAL(b.getMin(), 2)
  TEST_EQUAL(b.getMax(), 2);
END_SECTION

START_SECTION(double getMin() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(double getMax() const)
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION(void extend(const RangeBase& other))
  RangeBase b(4, 6);
  RangeBase other(1, 8);  
  b.extend(other);
  TEST_EQUAL(b.getMin(), 1)
  TEST_EQUAL(b.getMax(), 8)
END_SECTION

/// extend the range such that it includes the given @p value
START_SECTION(void extend(const double value))
  RangeBase b(4, 6);
  b.extend(1);
  TEST_EQUAL(b.getMin(), 1)
  TEST_EQUAL(b.getMax(), 6)
  RangeBase b2(4, 6);
  b2.extend(8);
  TEST_EQUAL(b2.getMin(), 4)
  TEST_EQUAL(b2.getMax(), 8)
  RangeBase b3(4, 6);
  b3.extend(5);
  TEST_EQUAL(b3.getMin(), 4)
  TEST_EQUAL(b3.getMax(), 6)
END_SECTION

START_SECTION(void extendLeftRight(const double by))
  RangeBase b(4, 6);
  b.extendLeftRight(1);
  TEST_EQUAL(b.getMin(), 3)
  TEST_EQUAL(b.getMax(), 7)
  RangeBase b2(2, 8);
  b2.extendLeftRight(-2);
  TEST_EQUAL(b2.getMin(), 4)
  TEST_EQUAL(b2.getMax(), 6)
  b2.extendLeftRight(-19);
  TEST_TRUE(b2.isEmpty())
  RangeBase empty;
  empty.extendLeftRight(100);
  TEST_TRUE(empty.isEmpty())
END_SECTION

START_SECTION(void clampTo(const RangeBase& other))
  RangeBase b(-4, 6);
  b.clampTo({-2, 7});
  TEST_EQUAL(b.getMin(), -2)
  TEST_EQUAL(b.getMax(), 6)
  RangeBase b2(4, 6);
  b2.clampTo({1, 5});
  TEST_EQUAL(b2.getMin(), 4)
  TEST_EQUAL(b2.getMax(), 5)
  b2.clampTo({4.5, 4.5});
  TEST_EQUAL(b2.getMin(), 4.5)
  TEST_EQUAL(b2.getMax(), 4.5)
  RangeBase b3(4, 6);
  b3.clampTo({10, 11});
  TEST_TRUE(b3.isEmpty())
  RangeBase b4(4, 6);
  TEST_EXCEPTION(Exception::InvalidRange, b4.clampTo(RangeBase()))
END_SECTION

START_SECTION(void pushInto(const RangeBase& sandbox))
  RangeBase b(-4, 6);
  b.pushInto({-2, 7});        // moves and clips
  TEST_EQUAL(b.getMin(), -2)
  TEST_EQUAL(b.getMax(), 7)
  RangeBase b2(4, 6);
  b2.pushInto({1, 15});       // does nothing
  TEST_EQUAL(b2.getMin(), 4)
  TEST_EQUAL(b2.getMax(), 6)
  b2.pushInto({4.5, 4.5});    // hard clip inside old range
  TEST_EQUAL(b2.getMin(), 4.5)
  TEST_EQUAL(b2.getMax(), 4.5)
  RangeBase b3(4, 6);         // move left, no clip
  b3.pushInto({-10, 5});
  TEST_EQUAL(b3.getMin(), 3)
  TEST_EQUAL(b3.getMax(), 5)
  b3.pushInto({4, 10});       // move right, no clip
  TEST_EQUAL(b3.getMin(), 4)
  TEST_EQUAL(b3.getMax(), 6)
  RangeBase b4(4, 6);         // hard clip outside old range
  b4.pushInto({-10, -10});
  TEST_EQUAL(b4.getMin(), -10)
  TEST_EQUAL(b4.getMax(), -10)
  RangeBase b5(4, 6);
  TEST_EXCEPTION(Exception::InvalidRange, b5.pushInto(RangeBase()))
  END_SECTION

START_SECTION(void scaleBy(const double factor))
  RangeBase b(4, 6);
  b.scaleBy(10); // diff is 2, so extend distance to 20, by increase of 9 on each side
  TEST_EQUAL(b.getMin(), 4-9)
  TEST_EQUAL(b.getMax(), 6+9)

  // scaling empty ranges does nothing
  RangeBase empty1, empty2;
  empty1.scaleBy(10);
  TEST_EQUAL(empty1, empty2)
END_SECTION

START_SECTION(void shift(const double distance))
  RangeBase b(4, 6);
  b.shift(10);
  TEST_EQUAL(b.getMin(), 14)
  TEST_EQUAL(b.getMax(), 16)

  // shifting empty ranges does nothing
  RangeBase empty1, empty2;
  empty1.shift(10);
  TEST_EQUAL(empty1, empty2)
END_SECTION

START_SECTION(double center() const)
  RangeBase b(4, 6);
  TEST_EQUAL(b.center(), 5)

  RangeBase empty;
  TEST_TRUE(std::isnan(empty.center()));
END_SECTION

START_SECTION(double getSpan() const)
  RangeBase b(4, 6);
  TEST_EQUAL(b.getSpan(), 2)

  RangeBase empty;
  TEST_TRUE(std::isnan(empty.getSpan()));
END_SECTION


START_SECTION(bool operator==(const RangeBase& rhs) const)
  RangeBase b(4, 6), b2(4, 6), empty1, empty2;
  TEST_NOT_EQUAL(b, empty1)
  TEST_EQUAL(b, b2)
  TEST_EQUAL(empty1, empty2)
END_SECTION


RM* ptr;
RM* nullPointer = nullptr;
START_SECTION((RangeMType()))
  ptr = new RM();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~RangeMType()))
  delete ptr;
END_SECTION

START_SECTION((double getMinMZ() const))
  TEST_EQUAL(RM().getMinMZ(), std::numeric_limits<double>::max())
END_SECTION

START_SECTION((double getMaxMZ() const))
  TEST_EQUAL(RM().getMaxMZ(), -std::numeric_limits<double>::max())
END_SECTION

START_SECTION((double getMinIntensity() const))
TEST_EQUAL(RM().getMinIntensity(), std::numeric_limits<double>::max())
END_SECTION

START_SECTION((double getMaxIntensity() const))
TEST_EQUAL(RM().getMaxIntensity(), -std::numeric_limits<double>::max())
END_SECTION

START_SECTION((double getMinMobility() const))
TEST_EQUAL(RM().getMinMobility(), std::numeric_limits<double>::max())
END_SECTION

START_SECTION((double getMaxMobility() const))
TEST_EQUAL(RM().getMaxMobility(), -std::numeric_limits<double>::max())
END_SECTION



START_SECTION((RangeManager(const RangeManager& rhs)))
  RM rm0;
  rm0.updateRanges();
  RM rm(rm0);
  TEST_REAL_SIMILAR(rm.getMinRT(), 2.0)
  TEST_REAL_SIMILAR(rm.getMinMZ(), 500.0)
  TEST_REAL_SIMILAR(rm.getMaxRT(), 100.0)
  TEST_REAL_SIMILAR(rm.getMaxMZ(), 1300.0)
  TEST_REAL_SIMILAR(rm.getMinIntensity(), 1.0)
  TEST_REAL_SIMILAR(rm.getMaxIntensity(), 47110.0)
END_SECTION

START_SECTION((RangeManager& operator=(const RangeManager& rhs)))
  RM rm0;
  rm0.updateRanges();
  RM rm;
  rm = rm0;
  TEST_REAL_SIMILAR(rm.getMinRT(), 2.0)
  TEST_REAL_SIMILAR(rm.getMinMZ(), 500.0)
  TEST_REAL_SIMILAR(rm.getMaxRT(), 100.0)
  TEST_REAL_SIMILAR(rm.getMaxMZ(), 1300.0)
  TEST_REAL_SIMILAR(rm.getMinIntensity(), 1.0)
  TEST_REAL_SIMILAR(rm.getMaxIntensity(), 47110.0)
END_SECTION

START_SECTION((bool operator==(const RangeManager& rhs) const))
  RM rm0 , rm;
  TEST_TRUE(rm == rm0);
  rm0.updateRanges();
  TEST_EQUAL(rm==rm0, false);
END_SECTION

START_SECTION((bool operator!=(const RangeManager& rhs) const))
  RM rm0 , rm;
  TEST_EQUAL(rm!=rm0, false);
  rm0.updateRanges();
  TEST_FALSE(rm == rm0);
END_SECTION

START_SECTION((virtual void updateRanges()=0))
  RM rm;

  rm.updateRanges();
  rm.updateRanges(); //second time to check the initialization

  TEST_REAL_SIMILAR(rm.getMinRT(), 2.0)
  TEST_REAL_SIMILAR(rm.getMinMZ(), 500.0)
  TEST_EQUAL(rm.RangeRT::isEmpty(), false)
  TEST_REAL_SIMILAR(rm.getMaxRT(), 100.0)
  TEST_REAL_SIMILAR(rm.getMaxMZ(), 1300.0)
  TEST_EQUAL(rm.RangeMZ::isEmpty(), false)
  TEST_REAL_SIMILAR(rm.getMinIntensity(), 1.0)
  TEST_REAL_SIMILAR(rm.getMaxIntensity(), 47110.0)
  TEST_EQUAL(rm.RangeIntensity::isEmpty(), false)
  TEST_EQUAL(rm.getMinMobility(), std::numeric_limits<double>::max())
  TEST_EQUAL(rm.getMaxMobility(), -std::numeric_limits<double>::max())
  TEST_EQUAL(rm.RangeMobility::isEmpty(), true)

  //test with only one point
  rm.updateRanges2(); //second time to check the initialization

  TEST_REAL_SIMILAR(rm.getMinRT(), 2.0)
  TEST_REAL_SIMILAR(rm.getMinMZ(), 500.0)
  TEST_REAL_SIMILAR(rm.getMaxRT(), 2.0)
  TEST_REAL_SIMILAR(rm.getMaxMZ(), 500.0)
  TEST_REAL_SIMILAR(rm.getMinIntensity(), 1.0)
  TEST_REAL_SIMILAR(rm.getMaxIntensity(), 1.0)
END_SECTION


START_SECTION(HasRangeType hasRange() const)
  RM rm;
  TEST_EQUAL(rm.hasRange() == HasRangeType::NONE, true);
  rm.updateRanges();
  TEST_EQUAL(rm.hasRange() == HasRangeType::SOME, true);
  rm.extendMobility(56.4);
  TEST_EQUAL(rm.hasRange() == HasRangeType::ALL, true);
END_SECTION

START_SECTION(template<typename... RangeBasesOther>
              bool containsAll(const RangeManager<RangeBasesOther...>& rhs) const)
  RM rm;
  rm.updateRanges();
  RM outer = rm;
  TEST_EQUAL(rm.containsAll(outer), true);
  TEST_EQUAL(outer.containsAll(rm), true);
  outer.scaleBy(1.1);
  TEST_EQUAL(rm.containsAll(outer), false);
  TEST_EQUAL(outer.containsAll(rm), true);
  outer.scaleBy(0.5);
  TEST_EQUAL(rm.containsAll(outer), true);
  TEST_EQUAL(outer.containsAll(rm), false);
  
  outer = rm;
  // empty dimensions in the rhs are considered contained
  outer.extendMobility(56.4); // rm.mobility is empty
  TEST_EQUAL(rm.containsAll(outer), false);
  TEST_EQUAL(outer.containsAll(rm), true);
  // empty dimensions do not count
  outer.RangeMZ::scaleBy(0.5); // mz range is smaller
  rm.RangeMZ::clear();         // but now does not count anymore
  TEST_EQUAL(rm.containsAll(outer), false); // due to mobility from above
  TEST_EQUAL(outer.containsAll(rm), true);

  // no ranges overlap...
  RangeManager<RangeRT, RangeMZ> rmz;
  RangeManager<RangeIntensity, RangeMobility> im;
  TEST_EXCEPTION(Exception::InvalidRange, rmz.containsAll(im))

END_SECTION

START_SECTION(template<typename... RangeBasesOther>
              void extend(const RangeManager<RangeBasesOther...>& rhs))
  RM rm;
  rm.updateRanges();
  RangeMTypeMzInt mid;
  mid.assign(rm); // assigns only overlapping dimensions
  TEST_REAL_SIMILAR(mid.getMinMZ(), 500.0)
  TEST_REAL_SIMILAR(mid.getMaxMZ(), 1300.0)
  TEST_REAL_SIMILAR(mid.getMinIntensity(), 1.0)
  TEST_REAL_SIMILAR(mid.getMaxIntensity(), 47110.0)

  RangeMTypeInt small;
  small.extendIntensity(123456.7);
  mid.extend(small);
  TEST_REAL_SIMILAR(mid.getMinMZ(), 500.0)
  TEST_REAL_SIMILAR(mid.getMaxMZ(), 1300.0)
  TEST_REAL_SIMILAR(mid.getMinIntensity(), 1.0)
  TEST_REAL_SIMILAR(mid.getMaxIntensity(), 123456.7)
END_SECTION

START_SECTION(void scaleBy(const double factor))
  RM rm;
  rm.updateRanges();
  rm.scaleBy(2);
  TEST_REAL_SIMILAR(rm.getMinRT(), 2.0-49)
  TEST_REAL_SIMILAR(rm.getMaxRT(), 100.0+49)
  TEST_REAL_SIMILAR(rm.getMinMZ(), 500.0-400)
  TEST_REAL_SIMILAR(rm.getMaxMZ(), 1300.0+400)
  TEST_REAL_SIMILAR(rm.getMinIntensity(), 1.0 - (47109.0/2))
  TEST_REAL_SIMILAR(rm.getMaxIntensity(), 47110.0 + (47109.0/2))
  TEST_TRUE(rm.RangeMobility::isEmpty())

  // scaling a dimension where min == max does nothing
  RangeManager<RangeRT, RangeMZ> rtmz;
  rtmz.extendMZ(100);
  rtmz.extendRT(50);
  auto copy = rtmz;
  rtmz.scaleBy(2.0);
  TEST_EQUAL(rtmz, copy)

  // scaling empty dimensions does nothing
  RM rm_empty, rm_empty2;
  rm_empty.scaleBy(4);
  TEST_EQUAL(rm_empty, rm_empty2)
END_SECTION

START_SECTION(template<typename... RangeBasesOther> void pushInto(const RangeManager<RangeBasesOther...>& sandbox))
  RM rm;
  rm.updateRanges();
  RangeMTypeMzInt rmi;
  rmi.extendMZ(700);     // shift
  rmi.extendMZ(2000);
  rmi.extendIntensity(500);  // shift and clamp
  rmi.extendIntensity(600);
  rm.pushInto(rmi);
  TEST_REAL_SIMILAR(rm.getMinRT(), 2.0)
  TEST_REAL_SIMILAR(rm.getMaxRT(), 100.0)
  TEST_REAL_SIMILAR(rm.getMinMZ(), 500.0 + 200)
  TEST_REAL_SIMILAR(rm.getMaxMZ(), 1300.0 + 200)
  TEST_REAL_SIMILAR(rm.getMinIntensity(), 1.0 + 499)
  TEST_REAL_SIMILAR(rm.getMaxIntensity(), 600) // was 47110.0
  TEST_TRUE(rm.RangeMobility::isEmpty())

  // if no dimensions overlap...
  RangeMTypeRT rt;
  TEST_EXCEPTION(Exception::InvalidRange, rmi.pushInto(rt))
END_SECTION


START_SECTION(template<typename... RangeBasesOther> void clampTo(const RangeManager<RangeBasesOther...>& rhs))
  RM rm;
  rm.updateRanges();
  RangeMTypeMzInt rmi;
  rmi.extendMZ(700);         // clamp left
  rmi.extendMZ(2000);
  rmi.extendIntensity(-10);  // clamp to empty
  rmi.extendIntensity(-9);
  rm.clampTo(rmi);
  TEST_REAL_SIMILAR(rm.getMinRT(), 2.0)              // should be untouched (since rmi.RT is empty)
  TEST_REAL_SIMILAR(rm.getMaxRT(), 100.0)
  TEST_REAL_SIMILAR(rm.getMinMZ(), 500.0 + 200)
  TEST_REAL_SIMILAR(rm.getMaxMZ(), 1300.0 + 0)
  TEST_TRUE(rm.getRangeForDim(MSDim::INT).isEmpty()) 
  TEST_TRUE(rm.RangeMobility::isEmpty())
  
  // if no dimensions overlap...
  RangeMTypeRT rt;
  TEST_EXCEPTION(Exception::InvalidRange, rmi.clampTo(rt))
END_SECTION

START_SECTION(RangeBase& getRangeForDim(MSDim dim))
  RM rm;
  rm.updateRanges();
  auto rt = rm.getRangeForDim(MSDim::RT);
  auto mz = rm.getRangeForDim(MSDim::MZ);
  auto in = rm.getRangeForDim(MSDim::INT);
  auto im = rm.getRangeForDim(MSDim::IM);
  TEST_REAL_SIMILAR(rt.getMin(), 2.0)
  TEST_REAL_SIMILAR(mz.getMin(), 500.0)
  TEST_REAL_SIMILAR(rt.getMax(), 100.0)
  TEST_REAL_SIMILAR(mz.getMax(), 1300.0)
  TEST_REAL_SIMILAR(in.getMin(), 1.0)
  TEST_REAL_SIMILAR(in.getMax(), 47110.0)
  TEST_FALSE(rt.isEmpty())
  TEST_TRUE(im.isEmpty())
END_SECTION

START_SECTION((void clearRanges()))
  RM rm;
  rm.updateRanges();
  TEST_REAL_SIMILAR(rm.getMinRT(), 2.0)
  TEST_REAL_SIMILAR(rm.getMinMZ(), 500.0)
  TEST_REAL_SIMILAR(rm.getMaxRT(), 100.0)
  TEST_REAL_SIMILAR(rm.getMaxMZ(), 1300.0)
  TEST_REAL_SIMILAR(rm.getMinIntensity(), 1.0)
  TEST_REAL_SIMILAR(rm.getMaxIntensity(), 47110.0)
  TEST_EQUAL(rm.RangeRT::isEmpty(), false)
  TEST_EQUAL(rm.RangeMZ::isEmpty(), false)
  TEST_EQUAL(rm.RangeIntensity::isEmpty(), false)
  TEST_EQUAL(rm.RangeMobility::isEmpty(), true)

  rm.clearRanges();
  TEST_EQUAL(rm.getMinRT(), std::numeric_limits<double>::max())
  TEST_EQUAL(rm.getMaxRT(), -std::numeric_limits<double>::max())
  TEST_REAL_SIMILAR(rm.getMinIntensity(), numeric_limits<double>::max())
  TEST_REAL_SIMILAR(rm.getMaxIntensity(), -numeric_limits<double>::max())
  TEST_EQUAL(rm.RangeRT::isEmpty(), true)
  TEST_EQUAL(rm.RangeMZ::isEmpty(), true)
  TEST_EQUAL(rm.RangeIntensity::isEmpty(), true)
  TEST_EQUAL(rm.RangeMobility::isEmpty(), true)
END_SECTION


START_SECTION(void printRange(std::ostream& out) const)
  RM rm;
  rm.extendRT(1.0);
  rm.extendMZ(2.0);
  rm.extendIntensity(3.0);
  rm.extendMobility(4.0);
  stringstream ss;
  rm.printRange(ss);
  TEST_EQUAL(ss.str(), "rt: [1, 1]\n"
                       "mz: [2, 2]\n"
                       "intensity: [3, 3]\n"
                       "mobility: [4, 4]\n");
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
