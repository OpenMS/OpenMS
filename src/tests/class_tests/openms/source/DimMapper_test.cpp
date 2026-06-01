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
#include <OpenMS/KERNEL/DimMapper.h>
///////////////////////////

#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/ConsensusFeature.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/KERNEL/Peak2D.h>

#include <limits>

using namespace OpenMS;
using namespace std;

START_TEST(DimMapper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


START_SECTION(DimRT())
{
  DimRT rt;
  TEST_TRUE(rt.clone()->getUnit() == DIM_UNIT::RT);
}
END_SECTION

START_SECTION(std::unique_ptr<DimBase> clone() const override)
{
  DimRT rt;
  TEST_TRUE(rt.clone()->getUnit() == DIM_UNIT::RT);
}
END_SECTION

START_SECTION(ValueType map(const Peak1D& p) const override)
{
  DimRT rt;
  TEST_EXCEPTION(Exception::InvalidRange, rt.map(Peak1D{1, 2}));
}
END_SECTION

START_SECTION(ValueType map(const Peak2D& p) const override)
{
  DimRT rt;
  TEST_EQUAL(rt.map(Peak2D({1, 2}, 3)), 1.0)
}
END_SECTION

START_SECTION(ValueTypes map(const MSSpectrum& spec) const override)
{
  DimRT rt;
  TEST_EXCEPTION(Exception::InvalidRange, rt.map(MSSpectrum()))
}
END_SECTION

START_SECTION(ValueType map(MSExperiment::ConstAreaIterator it) const override)
{
  DimRT rt;
  MSSpectrum spec;
  spec.push_back({1, 2});
  spec.setRT(5);
  MSExperiment exp;
  exp.addSpectrum(spec);
  TEST_EQUAL(rt.map(exp.areaBeginConst(4, 6, 0, 2)), 5)
}
END_SECTION

START_SECTION(ValueType map(const BaseFeature& bf) const override)
{
  DimRT rt;
  TEST_EQUAL(rt.map(BaseFeature{Peak2D({1, 2}, 3)}), 1)
}
END_SECTION

START_SECTION(ValueType map(const PeptideIdentification& pi) const override)
{
  DimRT rt;
  PeptideIdentification pi;
  pi.setRT(1);
  TEST_EQUAL(rt.map(pi), 1)
}
END_SECTION

START_SECTION(RangeBase map(const RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>& rm) const override)
{
  DimRT rt;
  RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility> rm;
  rm.extendRT(1);
  rm.extendRT(1.1);
  rm.extendMZ(2);
  rm.extendIntensity(3);
  TEST_EQUAL(rt.map(rm), RangeBase(1, 1.1));
}
END_SECTION

START_SECTION(void setRange(const RangeBase& in, RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>& out) const)
{
  DimRT rt;
  RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility> rm;
  rm.extendRT(1);
  rm.extendRT(1.1);
  rm.extendMZ(2);
  rm.extendIntensity(3);
  auto rm_old = rm;
  rt.setRange(RangeBase{10, 10.1}, rm);
  rm_old.RangeRT::operator=(RangeBase{10, 10.1});
  TEST_TRUE(rm == rm_old);
}
END_SECTION


using DimMapper3 = DimMapper<3>;

DimMapper3* ptr = nullptr;
DimMapper3* null_ptr = nullptr;

const DIM_UNIT unitsIMR[3] {DIM_UNIT::INT, DIM_UNIT::MZ, DIM_UNIT::RT};
const DIM_UNIT unitsRMI[3] {DIM_UNIT::RT, DIM_UNIT::MZ, DIM_UNIT::INT};

START_SECTION(DimMapper(const DIM_UNIT (&units)[N_DIM]))
{
  ptr = new DimMapper3(unitsIMR);
  TEST_NOT_EQUAL(ptr, null_ptr)
  delete ptr;
}
END_SECTION

START_SECTION(DimMapper(const DimMapper& rhs))
{
  DimMapper3 d1(unitsIMR);
  auto d2(d1);
  TEST_TRUE(d2 == d1);
}
END_SECTION


START_SECTION(DimMapper& operator=(const DimMapper& rhs))
{
  DimMapper3 d1(unitsIMR);
  DimMapper3 d2(unitsRMI);
  TEST_EQUAL(d2 == d1, false);
  d1 = d2;
  TEST_TRUE(d2 == d1);
}
END_SECTION

START_SECTION(bool operator==(const DimMapper& rhs) const)
{
  DimMapper3 d1(unitsIMR);
  DimMapper3 d2(unitsIMR);
  TEST_TRUE(d2 == d1);
  DimMapper3 d3(unitsRMI);
  TEST_TRUE(d3 != d1);
}
END_SECTION

START_SECTION(bool operator!=(const DimMapper& rhs) const)
{
  DimMapper3 d1(unitsIMR);
  DimMapper3 d2(unitsIMR);
  TEST_FALSE(d2 != d1);
  DimMapper3 d3(unitsRMI);
  TEST_FALSE(d3 == d1);
}
END_SECTION


START_SECTION(template<typename T> Point map(const T& data))
{
  DimMapper3 d1(unitsIMR);
  Feature f1;
  f1.setRT(1);
  f1.setMZ(2);
  f1.setIntensity(3);

  TEST_EQUAL(d1.map(f1), DimMapper3::Point(3, 2, 1))
}
END_SECTION

using FullRange = RangeManager<RangeRT, RangeMZ, RangeIntensity, RangeMobility>;

//auto na = nan("");
//constexpr auto dmax = std::numeric_limits<double>::max();

START_SECTION(template<typename... Ranges> DRange<N_DIM> mapRange(const RangeManager<Ranges...>& ranges) const)
{
  FullRange fr;
  fr.extendMobility(4); // not considered
  fr.extendRT(1);
  fr.extendRT(1.1);
  DimMapper3 d1(unitsIMR);
  auto areaXY = d1.mapRange(fr);
  // RT is Z-dimension:
  DRange<3> resXY;
  resXY.setDimMinMax(2, {1, 1.1});
  TEST_EQUAL(areaXY, resXY)
}
END_SECTION

START_SECTION(template<typename... Ranges> void fromXY(const DRange<N_DIM>& in, const RangeManager<Ranges...>& output) const)
{
  FullRange fr;
  fr.extendMobility(-4); // not considered
  fr.extendRT(12134);
  DimMapper3 d1(unitsIMR);
  // RT is Z-dimension:
  DRange<3> areaXY(DPosition<3>(77, 99, 1), DPosition<3>(777, 999, 1.1));
  d1.fromXY(areaXY, fr);
  TEST_EQUAL(fr.getMinRT(), 1)         // overwritten
  TEST_EQUAL(fr.getMaxRT(), 1.1)
  TEST_EQUAL(fr.getMinMZ(), 99)        // overwritten
  TEST_EQUAL(fr.getMaxMZ(), 999)
  TEST_EQUAL(fr.getMinIntensity(), 77)  // overwritten
  TEST_EQUAL(fr.getMaxIntensity(), 777)
  TEST_EQUAL(fr.getMinMobility(), -4)   // not modified
  TEST_EQUAL(fr.getMaxMobility(), -4)
}
END_SECTION

START_SECTION(template<typename... Ranges> void fromXY(const Point& in, RangeManager<Ranges...>& output) const)
{
  FullRange fr;
  fr.extendMobility(-4); // not considered
  fr.extendRT(12134);
  DimMapper3 d1(unitsIMR);
  // RT is Z-dimension:
  DRange<3> areaXY(DPosition<3>(77, 99, 1), DPosition<3>(777, 999, 1.1));
  d1.fromXY(DimMapper3::Point{2, 3, 1}, fr);
  TEST_EQUAL(fr.getMinRT(), 1) // overwritten
  TEST_EQUAL(fr.getMaxRT(), 1)
  TEST_EQUAL(fr.getMinMZ(), 3) // overwritten
  TEST_EQUAL(fr.getMaxMZ(), 3)
  TEST_EQUAL(fr.getMinIntensity(), 2) // overwritten
  TEST_EQUAL(fr.getMaxIntensity(), 2)
  TEST_EQUAL(fr.getMinMobility(), -4) // not modified
  TEST_EQUAL(fr.getMaxMobility(), -4)
}
END_SECTION


START_SECTION(const DimBase& getDim(DIM d) const)
{
  DimMapper3 d1(unitsIMR);
  TEST_TRUE(d1.getDim(DIM::X).getUnit() == DIM_UNIT::INT)
  TEST_TRUE(d1.getDim(DIM::Y).getUnit() == DIM_UNIT::MZ)
  TEST_TRUE(d1.getDim(DIM::Z).getUnit() == DIM_UNIT::RT)
}
END_SECTION

DimMapper3 dm_IMR(unitsIMR);
DimMapper3 dm_RMI(unitsRMI);
using Area3 = Area<3>;

/////// TEST for Area class
START_SECTION(Area(const DimMapper<N_DIM>* const dims))
{
  Area3 a(&dm_IMR);
  NOT_TESTABLE // tested below
}
END_SECTION

START_SECTION(Area(const Area& range) = default)
{
  Area3 a(&dm_IMR);
  Area3 o(a);
  TEST_TRUE(a == o)
}
END_SECTION

START_SECTION(Area& operator=(const Area& rhs) = default)
{
  Area3 a(&dm_IMR);
  auto ar = DRange<3>({1, 1, 1}, {2, 2, 2});
  a.setArea(ar);
  Area3 o(&dm_IMR);
  TEST_TRUE(a != o)
  o = a;
  TEST_TRUE(a == o)
  TEST_EQUAL(o.getAreaXY(), ar);
}
END_SECTION

START_SECTION(bool operator==(const Area& rhs) const)
{
  FullRange fr;
  fr.extendRT(1);
  Area3 a(&dm_IMR);
  Area3 o(&dm_IMR);
  TEST_TRUE(a == o)
  o = a;
  TEST_TRUE(a == o)
  a.setArea(fr);
  TEST_TRUE(a != o)
  o = a;
  TEST_TRUE(a == o)
  DRange<3> areaXY(DPosition<3>(77, 99, 1), DPosition<3>(777, 999, 1.1));
  a.setArea(areaXY);
  TEST_TRUE(a != o)
  TEST_FALSE(a == o)
}
END_SECTION

START_SECTION(bool operator!=(const Area& rhs) const)
{
  NOT_TESTABLE // tested above
}
END_SECTION

//constexpr auto min_value = -std::numeric_limits<DRange<3>::CoordinateType>::max();
//constexpr auto max_value = std::numeric_limits<DRange<3>::CoordinateType>::max();

START_SECTION(const Area& setArea(const UnitRange& data))
{
  FullRange fr;
  fr.RangeRT::operator=(RangeBase{1, 1.1});
  fr.RangeMobility::operator=(RangeBase{4, 4.4}); // not considered by DimMapper
  fr.RangeIntensity::operator=(RangeBase{2, 2.2});
  Area3 a(&dm_IMR);
  a.setArea(fr);
  TEST_EQUAL(fr, a.getAreaUnit()) // unchanged; just what we put in
  DRange<3> areaXY;
  areaXY.setDimMinMax(2, {1, 1.1}); // RT is mapped to dim2
  areaXY.setDimMinMax(0, {2, 2.2}); // Intensity is mapped to dim0
  TEST_EQUAL(a.getAreaXY(), areaXY)
}
END_SECTION


START_SECTION(const Area& setArea(const AreaXYType& data))
{
  DRange<3> areaXY;
  areaXY.setDimMinMax(2, {1, 1.1}); // RT is mapped to dim2
  areaXY.setDimMinMax(0, {2, 2.2}); // Intensity is mapped to dim0
  Area3 a(&dm_IMR);
  a.setArea(areaXY);
  TEST_EQUAL(a.getAreaXY(), areaXY) // unchanged; just what we put in

  FullRange fr;
  fr.RangeRT::operator=(RangeBase{1, 1.1});
  fr.RangeMobility::operator=(RangeBase{4, 4.4}); // not considered by DimMapper
  fr.RangeIntensity::operator=(RangeBase{2, 2.2});
  TEST_EQUAL((RangeRT)fr, (RangeRT)a.getAreaUnit())
  TEST_EQUAL((RangeIntensity)fr, (RangeIntensity)a.getAreaUnit())
  TEST_NOT_EQUAL(fr, a.getAreaUnit()) // due to mobility
}
END_SECTION


START_SECTION(const AreaXYType& getAreaXY() const)
{
  NOT_TESTABLE // tested above
}
END_SECTION


START_SECTION(const UnitRange& getAreaUnit() const)
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION(Area cloneWith(const AreaXYType& data) const)
{
  FullRange fr;
  fr.RangeRT::operator=(RangeBase{1, 1.1});
  fr.RangeMobility::operator=(RangeBase{4, 4.4}); // not considered by DimMapper
  fr.RangeIntensity::operator=(RangeBase{2, 2.2});
  Area3 a_old(&dm_IMR);
  auto a = a_old.cloneWith(fr);
  TEST_EQUAL(fr, a.getAreaUnit()) // unchanged; just what we put in
  DRange<3> areaXY;
  areaXY.setDimMinMax(2, {1, 1.1}); // RT is mapped to dim2
  areaXY.setDimMinMax(0, {2, 2.2}); // Intensity is mapped to dim0
  TEST_EQUAL(a.getAreaXY(), areaXY)
}
END_SECTION


START_SECTION(Area cloneWith(const UnitRange& data) const)
{
  DRange<3> areaXY;
  areaXY.setDimMinMax(2, {1, 1.1}); // RT is mapped to dim2
  areaXY.setDimMinMax(0, {2, 2.2}); // Intensity is mapped to dim0
  Area3 a_old(&dm_IMR);
  auto a = a_old.cloneWith(areaXY);
  TEST_EQUAL(a.getAreaXY(), areaXY) // unchanged; just what we put in

  FullRange fr;
  fr.RangeRT::operator=(RangeBase{1, 1.1});
  fr.RangeMobility::operator=(RangeBase{4, 4.4}); // not considered by DimMapper
  fr.RangeIntensity::operator=(RangeBase{2, 2.2});
  TEST_EQUAL((RangeRT)fr, (RangeRT)a.getAreaUnit())
  TEST_EQUAL((RangeIntensity)fr, (RangeIntensity)a.getAreaUnit())
  TEST_NOT_EQUAL(fr, a.getAreaUnit()) // due to mobility
}
END_SECTION
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



