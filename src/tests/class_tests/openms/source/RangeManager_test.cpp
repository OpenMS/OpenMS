// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/KERNEL/RangeManager.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

class RM
  : public RangeManager<2>
{
  public:
    RM()
      : RangeManager<2>()
    {

    }

    RM(const RM& rhs)
      : RangeManager<2>(rhs)
    {

    }

    RM& operator = (const RM& rhs)
    {
      if (this==&rhs) return *this;

      RangeManager<2>::operator=(rhs);

      return *this;
    }

    bool operator == (const RM& rhs) const
    {
      return
        RangeManager<2>::operator==(rhs);
        ;
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
      updateRanges_(vec.begin(), vec.end());
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
      updateRanges_(vec.begin(), vec.end());
    }

}; // class RM

START_TEST(RangeManager, "RangeManager")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RM* ptr;
RM* nullPointer = nullptr;
START_SECTION((RangeManager()))
  ptr = new RM();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~RangeManager()))
  delete ptr;
END_SECTION

START_SECTION((const PositionType& getMin() const))
  TEST_EQUAL(RM().getMin(), RM::PositionType::maxPositive())
END_SECTION

START_SECTION((const PositionType& getMax() const))
  TEST_EQUAL(RM().getMax(), RM::PositionType::minNegative())
END_SECTION

START_SECTION((double getMinInt() const ))
  TEST_REAL_SIMILAR(RM().getMinInt(), numeric_limits<double>::max())
END_SECTION

START_SECTION((double getMaxInt() const ))
  TEST_REAL_SIMILAR(RM().getMaxInt(), -numeric_limits<double>::max())
END_SECTION

START_SECTION((virtual void updateRanges()=0))
  RM rm;

  rm.updateRanges();
  rm.updateRanges(); //second time to check the initialization

  TEST_REAL_SIMILAR(rm.getMin()[0], 2.0)
  TEST_REAL_SIMILAR(rm.getMin()[1], 500.0)
  TEST_REAL_SIMILAR(rm.getMax()[0], 100.0)
  TEST_REAL_SIMILAR(rm.getMax()[1], 1300.0)
  TEST_REAL_SIMILAR(rm.getMinInt(), 1.0)
  TEST_REAL_SIMILAR(rm.getMaxInt(), 47110.0)

  //test with only one point
  rm.updateRanges2(); //second time to check the initialization

  TEST_REAL_SIMILAR(rm.getMin()[0], 2.0)
  TEST_REAL_SIMILAR(rm.getMin()[1], 500.0)
  TEST_REAL_SIMILAR(rm.getMax()[0], 2.0)
  TEST_REAL_SIMILAR(rm.getMax()[1], 500.0)
  TEST_REAL_SIMILAR(rm.getMinInt(), 1.0)
  TEST_REAL_SIMILAR(rm.getMaxInt(), 1.0)
END_SECTION

START_SECTION((void clearRanges()))
  RM rm;
  rm.updateRanges();
  TEST_REAL_SIMILAR(rm.getMin()[0], 2.0)
  TEST_REAL_SIMILAR(rm.getMin()[1], 500.0)
  TEST_REAL_SIMILAR(rm.getMax()[0], 100.0)
  TEST_REAL_SIMILAR(rm.getMax()[1], 1300.0)
  TEST_REAL_SIMILAR(rm.getMinInt(), 1.0)
  TEST_REAL_SIMILAR(rm.getMaxInt(), 47110.0)

  rm.clearRanges();
  TEST_EQUAL(RM().getMin(), RM::PositionType::maxPositive())
  TEST_EQUAL(RM().getMax(), RM::PositionType::minNegative())
  TEST_REAL_SIMILAR(RM().getMinInt(), numeric_limits<double>::max())
  TEST_REAL_SIMILAR(RM().getMaxInt(), -numeric_limits<double>::max())
END_SECTION

START_SECTION((RangeManager(const RangeManager& rhs)))
  RM rm0;
  rm0.updateRanges();
  RM rm(rm0);
  TEST_REAL_SIMILAR(rm.getMin()[0], 2.0)
  TEST_REAL_SIMILAR(rm.getMin()[1], 500.0)
  TEST_REAL_SIMILAR(rm.getMax()[0], 100.0)
  TEST_REAL_SIMILAR(rm.getMax()[1], 1300.0)
  TEST_REAL_SIMILAR(rm.getMinInt(), 1.0)
  TEST_REAL_SIMILAR(rm.getMaxInt(), 47110.0)
END_SECTION

START_SECTION((RangeManager& operator = (const RangeManager& rhs)))
  RM rm0;
  rm0.updateRanges();
  RM rm;
  rm = rm0;
  TEST_REAL_SIMILAR(rm.getMin()[0], 2.0)
  TEST_REAL_SIMILAR(rm.getMin()[1], 500.0)
  TEST_REAL_SIMILAR(rm.getMax()[0], 100.0)
  TEST_REAL_SIMILAR(rm.getMax()[1], 1300.0)
  TEST_REAL_SIMILAR(rm.getMinInt(), 1.0)
  TEST_REAL_SIMILAR(rm.getMaxInt(), 47110.0)
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
