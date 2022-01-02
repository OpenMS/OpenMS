// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

START_SECTION(void assign(const RangeBase& rhs))
  RangeBase b(4, 6), empty;
  empty.assign(b);
  TEST_EQUAL(empty.getMin(), 4)
  TEST_EQUAL(empty.getMax(), 6)
END_SECTION

START_SECTION(bool operator==(const RangeBase& rhs) const)
  RangeBase b(4, 6), empty;
  TEST_EQUAL(b == empty, false)
  TEST_EQUAL(b == b, true)
  TEST_EQUAL(empty == empty, true)
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
  TEST_EQUAL(rm.RangeMobility::isEmpty(), true)

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
  TEST_EQUAL(rt.isEmpty(), false)
  TEST_EQUAL(im.isEmpty(), true)
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

START_SECTION((RangeManager& operator = (const RangeManager& rhs)))
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

START_SECTION((bool operator == (const RangeManager& rhs) const))
  RM rm0 , rm;
  TEST_EQUAL(rm==rm0, true);
  rm0.updateRanges();
  TEST_EQUAL(rm==rm0, false);
END_SECTION

START_SECTION((bool operator != (const RangeManager& rhs) const))
  RM rm0 , rm;
  TEST_EQUAL(rm!=rm0, false);
  rm0.updateRanges();
  TEST_EQUAL(rm!=rm0, true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
