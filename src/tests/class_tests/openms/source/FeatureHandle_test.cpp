// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/FeatureHandle.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef FeatureMap ContainerType;
typedef ContainerType::value_type ElementType;
typedef Feature::PositionType PositionType;

START_TEST(FeatureHandle, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureHandle* ptr = nullptr;
FeatureHandle* nullPointer = nullptr;
START_SECTION((FeatureHandle()))
	ptr = new FeatureHandle();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureHandle()))
	delete ptr;
END_SECTION

START_SECTION((FeatureHandle& operator=(const FeatureHandle &rhs)))
  ElementType e;
  e.setUniqueId(2);
  FeatureHandle it(1,e);

  FeatureHandle it_copy;
  it_copy = it;

  TEST_EQUAL(it.getUniqueId() == it_copy.getUniqueId(), true)
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)
  TEST_EQUAL(it.getIntensity() == it_copy.getIntensity(), true)
  TEST_EQUAL(it.getPosition() == it_copy.getPosition(), true)
END_SECTION

START_SECTION((FeatureHandle(const FeatureHandle &rhs)))
  ElementType e;
  e.setUniqueId(2);
  FeatureHandle it(1,e);

  FeatureHandle it_copy(it);

  TEST_EQUAL(it.getUniqueId() == it_copy.getUniqueId(), true)
  TEST_EQUAL(it.getMapIndex() == it_copy.getMapIndex(), true)
  TEST_EQUAL(it.getIntensity() == it_copy.getIntensity(), true)
  TEST_EQUAL(it.getPosition() == it_copy.getPosition(), true)
END_SECTION

START_SECTION((void setCharge(ChargeType charge)))
{
  FeatureHandle fh;
  fh.setCharge(-17);
  TEST_EQUAL(fh.getCharge(),-17);
  fh.setCharge(-1717);
  TEST_EQUAL(fh.getCharge(),-1717);
}
END_SECTION

START_SECTION((ChargeType getCharge() const))
{
  NOT_TESTABLE; // see setCharge()
}
END_SECTION


START_SECTION((void setWidth(WidthType width)))
{
    FeatureHandle fh_tmp;
    fh_tmp.setWidth(10.7);
    TEST_REAL_SIMILAR(fh_tmp.getWidth(), 10.7);
    fh_tmp.setWidth(-8.9);
    TEST_REAL_SIMILAR(fh_tmp.getWidth(), -8.9);
}
END_SECTION

START_SECTION((WidthType getWidth() const ))
{
    NOT_TESTABLE;
}
END_SECTION


START_SECTION((FeatureHandle(UInt64 map_index, const Peak2D &point, UInt64 element_index)))
  ElementType e;
  FeatureHandle it(1,e,2);

  TEST_EQUAL(it.getUniqueId() == 2, true)
  TEST_EQUAL(it.getMapIndex() == 1, true)
  TEST_EQUAL(it.getPosition() == e.getPosition(), true)
END_SECTION

START_SECTION((FeatureHandle(UInt64 map_index, const BaseFeature& feature)))

  Feature f;
  f.setCharge(-17);
  f.setRT(44324.6);
  f.setMZ(867.4);
  f.setUniqueId(23);
  const Feature& f_cref = f;
  FeatureHandle fh(99,f_cref);

  TEST_EQUAL(fh.getMapIndex(),99);
  TEST_EQUAL(fh.getUniqueId(),23);
  TEST_EQUAL(fh.getRT(),44324.6);
  TEST_EQUAL(fh.getMZ(),867.4);
  TEST_EQUAL(fh.getCharge(),-17);

END_SECTION

START_SECTION((FeatureHandleMutable_ & asMutable() const))
  ConsensusFeature f;
  f.setCharge(-17);
  f.setRT(44324.6);
  f.setMZ(867.4);
  f.setUniqueId(23);
 const ConsensusFeature& f_cref = f;
  FeatureHandle fh(99, f_cref);

  const FeatureHandle& fh_cref = fh;
  // fh_cref.setRT(-64544.3); // compile time error
  fh_cref.asMutable().setRT(-64544.3); // ok

  TEST_EQUAL(fh.getMapIndex(),99);
  TEST_EQUAL(fh.getUniqueId(),23);
  TEST_EQUAL(fh.getRT(),-64544.3);
  TEST_EQUAL(fh.getMZ(),867.4);
  TEST_EQUAL(fh.getCharge(),-17);

END_SECTION


START_SECTION((bool operator!=(const FeatureHandle &i) const))
  ElementType e;
  e.setUniqueId(2);
  FeatureHandle it1(1,e);
  FeatureHandle it2(2,e);

  TEST_EQUAL(it1 != it2, true)
END_SECTION

START_SECTION((bool operator==(const FeatureHandle &i) const))
  ElementType e;
  e.setUniqueId(2);
  FeatureHandle it1(2,e);
  FeatureHandle it2(2,e);

  TEST_EQUAL(it1 == it2, true)
END_SECTION

START_SECTION((UInt64 getMapIndex() const))
  ElementType e;
  e.setUniqueId(2);
  FeatureHandle it(1,e);

  TEST_EQUAL(it.getMapIndex() == 1, true)
END_SECTION

START_SECTION((void setMapIndex(UInt64 i)))
  FeatureHandle it;
  it.setMapIndex(2);
  it.setUniqueId(77);

  TEST_EQUAL(it.getMapIndex() == 2, true)
END_SECTION

START_SECTION(([FeatureHandle::IndexLess] bool operator()(FeatureHandle const &left, FeatureHandle const &right) const))
  FeatureHandle lhs, rhs;
  lhs.setMapIndex(2);
  lhs.setUniqueId(77);
  rhs.setMapIndex(4);
  lhs.setUniqueId(29);

  FeatureHandle::IndexLess il;

  TEST_EQUAL(il(lhs, rhs), 1);
  TEST_EQUAL(il(rhs, lhs), 0);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



