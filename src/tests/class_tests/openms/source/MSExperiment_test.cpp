// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Marc Sturm, Tom Waschischeck $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/Peak2D.h>
#include <OpenMS/METADATA/DataProcessing.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


MSExperiment createPeakMapWithRTs(std::vector<double> RTs)
{
  MSExperiment map;
  MSSpectrum s;
  for (auto rt : RTs)
  {
    s.setRT(rt);
    map.addSpectrum(s);
  }
  return map;
}

MSExperiment setMSLevel(MSExperiment exp, std::vector<int> ms_levels)
{
  for (size_t i = 0; i < ms_levels.size(); ++i)
  {
    exp[i].setMSLevel(ms_levels[i]);
  }
  return exp;
}

START_TEST(MSExperiment, "$Id$");

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////



PeakMap* ptr = nullptr;
PeakMap* nullPointer = nullptr;
START_SECTION((MSExperiment()))
{
  ptr = new PeakMap;
  TEST_NOT_EQUAL(ptr, nullPointer);
}
END_SECTION

START_SECTION(([EXTRA]~MSExperiment()))
{
  delete ptr;
}
END_SECTION

/////////////////////////////////////////////////////////////
// Copy constructor, move constructor, assignment operator, move assignment operator, equality

START_SECTION((MSExperiment(const MSExperiment& source)))
{
  PeakMap tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setFirstName("Name");
  tmp.resize(1);

  PeakMap tmp2(tmp);
  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL(tmp2.getContacts()[0].getFirstName(),"Name");
  TEST_EQUAL(tmp2.size(),1);
}
END_SECTION

START_SECTION((MSExperiment(const MSExperiment&& source)))
{
  // Ensure that MSExperiment has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  // TODO: wont work for MSVS (yet)
#ifndef OPENMS_COMPILER_MSVC
  TEST_EQUAL(noexcept(MSExperiment(std::declval<MSExperiment&&>())), true)
#endif
  PeakMap tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setFirstName("Name");
  tmp.resize(1);

  //copy tmp so we can move one of them
  PeakMap orig = tmp;
  PeakMap tmp2(std::move(tmp));

  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL(tmp2.getContacts()[0].getFirstName(),"Name");
  TEST_EQUAL(tmp2.size(),1);

  // test move
  TEST_EQUAL(tmp.size(),0);
}
END_SECTION

START_SECTION((MSExperiment& operator= (const MSExperiment& source)))
{
  PeakMap tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setFirstName("Name");
  tmp.resize(1);
  Peak1D p;
  p.setMZ(5.0);
  tmp[0].push_back(p);
  p.setMZ(10.0);
  tmp[0].push_back(p);
  tmp.updateRanges();

  PeakMap tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL(tmp2.getContacts()[0].getFirstName(),"Name");
  TEST_EQUAL(tmp2.size(),1);
  TEST_REAL_SIMILAR(tmp2.getMinMZ(),5.0);
  TEST_REAL_SIMILAR(tmp2.getMaxMZ(),10.0);

  tmp2 = PeakMap();
  TEST_EQUAL(tmp2.getContacts().size(),0);
  TEST_EQUAL(tmp2.size(),0);
}
END_SECTION

START_SECTION((MSExperiment& operator= (const MSExperiment&& source)))
{
  PeakMap tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setFirstName("Name");
  tmp.resize(1);
  Peak1D p;
  p.setMZ(5.0);
  tmp[0].push_back(p);
  p.setMZ(10.0);
  tmp[0].push_back(p);
  tmp.updateRanges();

  //copy tmp so we can move one of them
  PeakMap orig = tmp;
  PeakMap tmp2 = std::move(tmp);

  TEST_EQUAL(tmp2, orig); // should be equal to the original

  TEST_EQUAL(tmp2.getContacts().size(),1);
  TEST_EQUAL(tmp2.getContacts()[0].getFirstName(),"Name");
  TEST_EQUAL(tmp2.size(),1);
  TEST_REAL_SIMILAR(tmp2.getMinMZ(),5.0);
  TEST_REAL_SIMILAR(tmp2.getMaxMZ(),10.0);

  // test move
  TEST_EQUAL(tmp.size(),0);

  tmp2 = PeakMap(); // use rvalue assignment
  TEST_EQUAL(tmp2.getContacts().size(),0);
  TEST_EQUAL(tmp2.size(),0);
}
END_SECTION

START_SECTION((bool operator== (const MSExperiment& rhs) const))
{
  PeakMap edit,empty;

  TEST_TRUE(edit == empty);

  edit.getContacts().resize(1);
  TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.resize(1);
  TEST_EQUAL(edit==empty, false);
}
END_SECTION

START_SECTION((bool operator!= (const MSExperiment& rhs) const))
{
  PeakMap edit,empty;

  TEST_EQUAL(edit!=empty, false);

  edit.getContacts().resize(1);
  TEST_FALSE(edit == empty);

  edit = empty;
  edit.resize(1);
  TEST_FALSE(edit == empty);
}
END_SECTION

START_SECTION((template<class Container> void get2DData(Container& cont) const))
{
  PeakMap exp;
  PeakMap::SpectrumType spec;
  PeakMap::PeakType peak;

  // first spectrum (MS)
  spec.setRT(11.1);
  spec.setMSLevel(1);
  peak.getPosition()[0] = 5;
  peak.setIntensity(47.11f);
  spec.push_back(peak);
  peak.getPosition()[0] = 10;
  peak.setIntensity(48.11f);
  spec.push_back(peak);
  peak.getPosition()[0] = 15;
  spec.push_back(peak);
  exp.addSpectrum(spec);

  // second spectrum (MS/MS)
  spec.clear(true);
  spec.setRT(11.5);
  spec.setMSLevel(2);
  peak.getPosition()[0] = 6;
  spec.push_back(peak);
  peak.getPosition()[0] = 11;
  spec.push_back(peak);
  exp.addSpectrum(spec);

  // third spectrum (MS)
  spec.clear(true);
  spec.setRT(12.2);
  spec.setMSLevel(1);
  peak.getPosition()[0] = 20;
  spec.push_back(peak);
  peak.getPosition()[0] = 25;
  spec.push_back(peak);
  exp.addSpectrum(spec);

  // forth spectrum (MS/MS)
  spec.clear(true);
  spec.setRT(12.5);
  spec.setMSLevel(2);
  peak.getPosition()[0] = 21;
  spec.push_back(peak);
  peak.getPosition()[0] = 26;
  spec.push_back(peak);
  peak.getPosition()[0] = 31;
  spec.push_back(peak);
  exp.addSpectrum(spec);

  //Convert
  std::vector<Peak2D> a;
  exp.get2DData(a);

  //Tests
  TEST_EQUAL(a.size(),5);
  TEST_REAL_SIMILAR(a[0].getRT(),11.1);
  TEST_REAL_SIMILAR(a[0].getMZ(),5);
  TEST_REAL_SIMILAR(a[0].getIntensity(),47.11);
  TEST_REAL_SIMILAR(a[1].getRT(),11.1);
  TEST_REAL_SIMILAR(a[1].getMZ(),10);
  TEST_REAL_SIMILAR(a[1].getIntensity(),48.11);
  TEST_REAL_SIMILAR(a[2].getRT(),11.1);
  TEST_REAL_SIMILAR(a[2].getMZ(),15);
  TEST_REAL_SIMILAR(a[3].getRT(),12.2);
  TEST_REAL_SIMILAR(a[3].getMZ(),20);
  TEST_REAL_SIMILAR(a[4].getRT(),12.2);
  TEST_REAL_SIMILAR(a[4].getMZ(),25);

  //Convert
  std::vector<Peak2D> list;
  exp.get2DData(list);

  //Tests
  TEST_EQUAL(list.size(),5);
  std::vector<Peak2D>::const_iterator it = list.begin();
  TEST_REAL_SIMILAR(it->getRT(),11.1);
  TEST_REAL_SIMILAR(it->getMZ(),5);
  TEST_REAL_SIMILAR(it->getIntensity(),47.11);
  ++it;
  TEST_REAL_SIMILAR(it->getRT(),11.1);
  TEST_REAL_SIMILAR(it->getMZ(),10);
  TEST_REAL_SIMILAR(it->getIntensity(),48.11);
  ++it;
  TEST_REAL_SIMILAR(it->getRT(),11.1);
  TEST_REAL_SIMILAR(it->getMZ(),15);
  ++it;
  TEST_REAL_SIMILAR(it->getRT(),12.2);
  TEST_REAL_SIMILAR(it->getMZ(),20);
  ++it;
  TEST_REAL_SIMILAR(it->getRT(),12.2);
  TEST_REAL_SIMILAR(it->getMZ(),25);
}
END_SECTION

START_SECTION((template <class Container> void set2DData(const Container& cont, const StringList& store_metadata_names = StringList())))
{
  NOT_TESTABLE // tested below
}
END_SECTION

START_SECTION((template <bool add_mass_traces, class Container> void set2DData(const Container& cont, const StringList& store_metadata_names = StringList())))
{
  PeakMap exp;

  // create sample data
  std::vector<Peak2D> input;

  Peak2D p1(Peak2D::PositionType(2.0, 3.0), 1.0);
  input.push_back(p1);

  Peak2D p2(Peak2D::PositionType(5.0, 6.0), 4.0);
  input.push_back(p2);

  Peak2D p3(Peak2D::PositionType(8.5, 9.5), 7.5);
  input.push_back(p3);

  exp.set2DData(input);

  // retrieve data again and check for changes
  std::vector<Peak2D> output;

  exp.get2DData(output);
  TEST_EQUAL(output==input,true);

  //////////////////////////////////////////////////////////////////////////
  // test if meta values are added as floatDataArrays in MSSpectra
  std::vector<RichPeak2D> inputr;

  RichPeak2D pr1(RichPeak2D::PositionType(2.0, 3.0), 1.0);
  pr1.setMetaValue("meta1", 111.1);
  inputr.push_back(pr1);
  RichPeak2D pr2(RichPeak2D::PositionType(5.0, 6.0), 4.0);
  inputr.push_back(pr2);
  RichPeak2D pr3(RichPeak2D::PositionType(8.5, 9.5), 7.5);
  pr3.setMetaValue("meta3", 333.3);
  inputr.push_back(pr3);

  // create float data arrays for these two meta values (missing values in data will be set to NaN)
  exp.set2DData(inputr, ListUtils::create<String>("meta1,meta3"));
  TEST_EQUAL(exp.getNrSpectra(), 3);
  // retrieve data again and check for changes
  std::vector<Peak2D> outputr;
  exp.get2DData(outputr);
  TEST_TRUE(outputr == input); // we compare to non-meta output, since floatdata is not converted back to metavalues
  // check for meta data
  TEST_EQUAL(exp[0].getFloatDataArrays().size(), 2);
  TEST_EQUAL(exp[0].getFloatDataArrays()[0][0], 111.1);
  TEST_EQUAL(exp[1].getFloatDataArrays().size(), 2); // present but all NaN
  TEST_EQUAL(exp[2].getFloatDataArrays().size(), 2);
  TEST_EQUAL(exp[2].getFloatDataArrays()[1][0], 333.3);

  ///////////////////////////////////////
  // test adding of mass traces
  FeatureMap fm, fm2, fm_out;
  Feature f1;
  f1.setIntensity(7.5f);
  f1.setRT(8.5);
  f1.setMZ(9.5);
  Feature f2;
  f2.setIntensity(17.5f);
  f2.setRT(18.5);
  f2.setMZ(19.5);
  fm.push_back(f1);
  fm.push_back(f2);  
  fm2 = fm; // copy without meta values (get2DData will not have them)
  fm.back().setMetaValue(Constants::UserParam::NUM_OF_MASSTRACES, 2);
  fm.back().setMetaValue("masstrace_intensity_0", 11.0f);
  fm.back().setMetaValue("masstrace_intensity_1", 12.0f);
  fm.back().setCharge(2);
  exp.set2DData<true>(fm);
  exp.get2DData(fm_out);
  Feature f2_x = f2;
  f2_x.setIntensity(11.0f);
  f2_x.setMZ(f2_x.getMZ() + Constants::C13C12_MASSDIFF_U/2*0);
  fm2.back() = (f2_x);   // replace
  f2_x.setIntensity(12.0f);
  f2_x.setMZ(f2_x.getMZ() + Constants::C13C12_MASSDIFF_U/2*1);
  fm2.push_back(f2_x);  // add +1Th trace
  TEST_EQUAL(fm_out.size(), fm2.size());
  //std::cout << fm_out[1] << "\n\n" << fm2[1] << std::endl;
  //std::cout << fm_out.back() << "\n\n" << fm2.back() << std::endl;
  TEST_EQUAL(fm_out==fm2,true);

  // test precondition (input sorted by RT)
  input.push_back(p1);
  TEST_PRECONDITION_VIOLATED(exp.set2DData(input));

}
END_SECTION

START_SECTION(([EXTRA] PeakMap()))
{
  PeakMap tmp;
  tmp.resize(1);
  tmp[0].resize(1);
  tmp[0][0].getPosition()[0] = 47.11;
  TEST_REAL_SIMILAR(tmp[0][0].getPosition()[0],47.11)
}
END_SECTION

START_SECTION((CoordinateType getMinMZ() const))
{
  PeakMap tmp;
  TEST_REAL_SIMILAR(tmp.getMinMZ(),numeric_limits<DPosition<2>::CoordinateType>::max())
}
END_SECTION

START_SECTION((CoordinateType getMaxMZ() const))
{
  PeakMap tmp;
  TEST_REAL_SIMILAR(tmp.getMaxMZ(),-numeric_limits<DPosition<2>::CoordinateType>::max())
}
END_SECTION

START_SECTION((CoordinateType getMinRT() const))
{
  PeakMap tmp;
  TEST_REAL_SIMILAR(tmp.getMinRT(),numeric_limits<DPosition<2>::CoordinateType>::max())
}
END_SECTION

START_SECTION((CoordinateType getMaxRT() const))
{
  PeakMap tmp;
  TEST_REAL_SIMILAR(tmp.getMaxRT(),-numeric_limits<DPosition<2>::CoordinateType>::max())
}
END_SECTION

START_SECTION((const std::vector<UInt>& getMSLevels() const))
{
  PeakMap tmp;
  TEST_EQUAL(tmp.getMSLevels().size(),0)
}
END_SECTION

START_SECTION((UInt64 getSize() const ))
{
  PeakMap tmp;
  TEST_EQUAL(tmp.getSize(),0)
}
END_SECTION

START_SECTION((const MSExperiment::RangeManagerType& MSExperiment::getRange() const))
{
  PeakMap tmp;
  TEST_EQUAL(tmp.getRange().hasRange() == HasRangeType::NONE, true)
}
END_SECTION

START_SECTION((virtual void updateRanges()))
{
  PeakMap tmp;
  MSSpectrum s;
  Peak1D p;

  s.setMSLevel(1);
  s.setRT(30.0);
  p.getPosition()[0] = 5.0;
  p.setIntensity(-5.0f);
  s.push_back(p);
  s.setDriftTime(99);
  tmp.addSpectrum(s);

  s.clear(true);
  s.setMSLevel(1);
  s.setRT(40.0);
  p.getPosition()[0] = 7.0;
  p.setIntensity(-7.0f);
  s.push_back(p);
  s.setDriftTime(99);
  tmp.addSpectrum(s);

  s.clear(true);
  s.setMSLevel(3);
  s.setRT(45.0);
  p.getPosition()[0] = 9.0;
  p.setIntensity(-10.0f);
  s.push_back(p);
  s.setDriftTime(199);
  tmp.addSpectrum(s);

  s.clear(true);
  s.setMSLevel(3);
  s.setRT(50.0);
  p.getPosition()[0] = 10.0;
  p.setIntensity(-9.0f);
  s.push_back(p);
  s.setDriftTime(66);
  tmp.addSpectrum(s);

  tmp.updateRanges();
  tmp.updateRanges(); //second time to check the initialization

  TEST_REAL_SIMILAR(tmp.getMinMZ(),5.0)
  TEST_REAL_SIMILAR(tmp.getMaxMZ(),10.0)
  TEST_REAL_SIMILAR(tmp.getMinIntensity(), -10.0)
  TEST_REAL_SIMILAR(tmp.getMaxIntensity(), -5.0)
  TEST_REAL_SIMILAR(tmp.getMinRT(),30.0)
  TEST_REAL_SIMILAR(tmp.getMaxRT(),50.0)
  TEST_EQUAL(tmp.getMSLevels().size(),2)
  TEST_EQUAL(tmp.getMSLevels()[0],1)
  TEST_EQUAL(tmp.getMSLevels()[1],3)
  TEST_EQUAL(tmp.getSize(),4)
  tmp.updateRanges();
  TEST_REAL_SIMILAR(tmp.getMinMZ(),5.0)
  TEST_REAL_SIMILAR(tmp.getMaxMZ(),10.0)
  TEST_REAL_SIMILAR(tmp.getMinIntensity(), -10.0)
  TEST_REAL_SIMILAR(tmp.getMaxIntensity(), -5.0)
  TEST_REAL_SIMILAR(tmp.getMinRT(),30.0)
  TEST_REAL_SIMILAR(tmp.getMaxRT(),50.0)

  TEST_REAL_SIMILAR(tmp.getRange().getMinMZ(), 5.0)
  TEST_REAL_SIMILAR(tmp.getRange().getMaxMZ(), 10.0)
  TEST_REAL_SIMILAR(tmp.getRange().getMinRT(), 30.0)
  TEST_REAL_SIMILAR(tmp.getRange().getMaxRT(), 50.0)
  TEST_REAL_SIMILAR(tmp.getRange().getMinMobility(), 66)
  TEST_REAL_SIMILAR(tmp.getRange().getMaxMobility(), 199)

  TEST_EQUAL(tmp.getMSLevels().size(),2)
  TEST_EQUAL(tmp.getMSLevels()[0],1)
  TEST_EQUAL(tmp.getMSLevels()[1],3)

  TEST_EQUAL(tmp.getSize(),4)

  //Update for MS level 1

  // Store initial MS levels
  std::vector<UInt> initial_ms_levels = tmp.getMSLevels();

  tmp.updateRanges(1);
  tmp.updateRanges(1); // Call twice to verify consistent behavior
  for (int l = 0; l < 2; ++l)
  {
    TEST_REAL_SIMILAR(tmp.getMinMZ(),5.0)
    TEST_REAL_SIMILAR(tmp.getMaxMZ(),7.0)
    TEST_REAL_SIMILAR(tmp.getMinIntensity(), -7.0)
    TEST_REAL_SIMILAR(tmp.getMaxIntensity(), -5.0)
    TEST_REAL_SIMILAR(tmp.getMinRT(),30.0)
    TEST_REAL_SIMILAR(tmp.getMaxRT(),40.0)
    TEST_REAL_SIMILAR(tmp.getRange().getMinMobility(), 99)
    TEST_REAL_SIMILAR(tmp.getRange().getMaxMobility(), 99)
    // Verify MS levels remain unchanged
    TEST_EQUAL(tmp.getMSLevels() == initial_ms_levels, true)
    TEST_EQUAL(tmp.getSize(),4)
    tmp.updateRanges(1);
  }

  // test with only one peak
  PeakMap tmp2;
  MSSpectrum s2;
  Peak1D p2;

  s2.setRT(30.0);
  p2.getPosition()[0] = 5.0;
  p2.setIntensity(-5.0f);
  s2.push_back(p2);
  s2.setDriftTime(99);
  tmp2.addSpectrum(s2);

  tmp2.updateRanges();
  TEST_REAL_SIMILAR(tmp2.getMinMZ(),5.0)
  TEST_REAL_SIMILAR(tmp2.getMaxMZ(),5.0)
  TEST_REAL_SIMILAR(tmp2.getMinIntensity(), -5.0)
  TEST_REAL_SIMILAR(tmp2.getMaxIntensity(), -5.0)
  TEST_REAL_SIMILAR(tmp2.getMinRT(),30.0)
  TEST_REAL_SIMILAR(tmp2.getMaxRT(),30.0)
  TEST_REAL_SIMILAR(tmp.getRange().getMinMobility(), 99)
  TEST_REAL_SIMILAR(tmp.getRange().getMaxMobility(), 99)

  tmp2.updateRanges(1);
  TEST_REAL_SIMILAR(tmp2.getMinMZ(),5.0)
  TEST_REAL_SIMILAR(tmp2.getMaxMZ(),5.0)
  TEST_REAL_SIMILAR(tmp2.getMinIntensity(), -5.0)
  TEST_REAL_SIMILAR(tmp2.getMaxIntensity(), -5.0)
  TEST_REAL_SIMILAR(tmp2.getMinRT(),30.0)
  TEST_REAL_SIMILAR(tmp2.getMaxRT(),30.0)
  TEST_REAL_SIMILAR(tmp.getRange().getMinMobility(), 99)
  TEST_REAL_SIMILAR(tmp.getRange().getMaxMobility(), 99)

  // test ranges with a chromatogram
  MSChromatogram chrom1, chrom2;
  ChromatogramPeak cp1, cp2, cp3;
  cp1.setRT(0.3);
  cp1.setIntensity(10.0f);
  cp2.setRT(0.2);
  cp2.setIntensity(10.2f);
  cp3.setRT(0.1);
  cp3.setIntensity(10.4f);

  Product prod1;
  prod1.setMZ(100.0);
  chrom1.setProduct(prod1);
  chrom1.push_back(cp1);
  chrom1.push_back(cp2);

  Product prod2;
  prod2.setMZ(80.0);
  chrom2.setProduct(prod2);
  chrom2.push_back(cp2);
  chrom2.push_back(cp3);

  vector<MSChromatogram> chroms;
  chroms.push_back(chrom1);
  chroms.push_back(chrom2);
  tmp2.setChromatograms(chroms);
  
  tmp2.updateRanges();
  TEST_REAL_SIMILAR(tmp2.getMinMZ(), 5.0)
  TEST_REAL_SIMILAR(tmp2.getMaxMZ(), 100.0)
  TEST_REAL_SIMILAR(tmp2.getMinIntensity(), -5.0)
  TEST_REAL_SIMILAR(tmp2.getMaxIntensity(), 10.4)
  TEST_REAL_SIMILAR(tmp2.getMinRT(), 0.1)
  TEST_REAL_SIMILAR(tmp2.getMaxRT(), 30.0)
}
END_SECTION

START_SECTION((void updateRanges(Int ms_level)))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((ConstAreaIterator areaEndConst() const))
{
  NOT_TESTABLE
}
END_SECTION

PeakMap exp_area;
{ // create some data for testing area iterator
  std::vector< Peak2D> plist;
  plist.push_back(Peak2D({-1.0, 2.0}, 0));
  plist.push_back(Peak2D({1.0, 2.0}, 0));
  plist.push_back(Peak2D({1.0, 3.0}, 0));
  plist.push_back(Peak2D({2.0, 10.0}, 0));
  plist.push_back(Peak2D({2.0, 11.0}, 0));
  plist.push_back(Peak2D({2.0, 12.0}, 0));
  exp_area.set2DData(plist);
  int dt{0};
  for (auto& spec : exp_area.getSpectra())
  {
    spec.setDriftTime(dt++);
  }
}

START_SECTION((ConstAreaIterator areaBeginConst(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz) const))
{
  const auto& exp = exp_area;
  PeakMap::ConstAreaIterator it = exp.areaBeginConst(0,15,0,11);

  TEST_EQUAL(it->getPosition()[0],2.0)
  it++;
  TEST_EQUAL(it->getPosition()[0],3.0)
  it++;
  TEST_EQUAL(it->getPosition()[0],10.0)
  it++;
  TEST_EQUAL(it->getPosition()[0],11.0)
  it++;
  TEST_EQUAL(it==exp.areaEndConst(),true)

  TEST_PRECONDITION_VIOLATED(exp.areaBeginConst(15,0,0,15))
  TEST_PRECONDITION_VIOLATED(exp.areaBeginConst(0,15,15,0))
  TEST_PRECONDITION_VIOLATED(exp.areaBeginConst(15,0,15,0))
}
END_SECTION

START_SECTION(MSExperiment::ConstAreaIterator MSExperiment::areaBeginConst(const RangeManagerType& range) const)
{
  MSExperiment::RangeManagerType rm;
  (RangeRT&)rm = RangeBase(0, 2);
  (RangeMZ&)rm = RangeBase(3, 11);

  const auto& exp = exp_area;
  PeakMap::ConstAreaIterator it = exp.areaBeginConst(rm);
  TEST_EQUAL(it->getPosition()[0], 3.0)
  ++it;
  TEST_EQUAL(it->getPosition()[0], 10.0)
  ++it;
  TEST_EQUAL(it->getPosition()[0], 11.0)
  ++it;
  TEST_EQUAL(it == exp.areaEndConst(), true)

  // add mobility as constraint
  (RangeMobility&)rm = RangeBase(2, 2);
  it = exp.areaBeginConst(rm);
  TEST_EQUAL(it->getPosition()[0], 10.0)
  ++it;
  TEST_EQUAL(it->getPosition()[0], 11.0)
  ++it;
  TEST_EQUAL(it == exp.areaEndConst(), true)
}
END_SECTION

START_SECTION((AreaIterator areaEnd()))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION((AreaIterator areaBegin(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz)))
{
  auto exp = exp_area;
  PeakMap::AreaIterator it = exp.areaBegin(0,15,0,11);

  TEST_EQUAL(it->getPosition()[0],2.0)
  it->getPosition()[0] = 4711.0;
  TEST_EQUAL(it->getPosition()[0],4711.0)
  ++it;
  TEST_EQUAL(it->getPosition()[0],3.0)
  ++it;
  TEST_EQUAL(it->getPosition()[0],10.0)
  ++it;
  TEST_EQUAL(it->getPosition()[0],11.0)
  ++it;
  TEST_EQUAL(it==exp.areaEnd(),true)

  TEST_PRECONDITION_VIOLATED(exp.areaBegin(15,0,0,15))
  TEST_PRECONDITION_VIOLATED(exp.areaBegin(0,15,15,0))
  TEST_PRECONDITION_VIOLATED(exp.areaBegin(15,0,15,0))
}
END_SECTION

START_SECTION(MSExperiment::AreaIterator MSExperiment::areaBegin(const RangeManagerType& range))
{
  MSExperiment::RangeManagerType rm;
  (RangeRT&)rm = RangeBase(0, 2);
  (RangeMZ&)rm = RangeBase(3, 11);

  auto exp = exp_area;
  PeakMap::AreaIterator it = exp.areaBegin(rm);
  TEST_EQUAL(it->getPosition()[0], 3.0)
  ++it;
  TEST_EQUAL(it->getPosition()[0], 10.0)
  ++it;
  TEST_EQUAL(it->getPosition()[0], 11.0)
  ++it;
  TEST_EQUAL(it == exp.areaEnd(), true)

  // add mobility as constraint
  (RangeMobility&)rm = RangeBase(2, 2);
  it = exp.areaBegin(rm);
  TEST_EQUAL(it->getPosition()[0], 10.0)
  ++it;
  TEST_EQUAL(it->getPosition()[0], 11.0)
  ++it;
  TEST_EQUAL(it == exp.areaEnd(), true)
}
END_SECTION

START_SECTION((Iterator RTBegin(CoordinateType rt)))
{
  PeakMap tmp = createPeakMapWithRTs({30, 40, 45, 50});
  PeakMap::Iterator it;

  it = tmp.RTBegin(20.0);
  TEST_REAL_SIMILAR(it->getRT(),30.0)
  it = tmp.RTBegin(30.0);
  TEST_REAL_SIMILAR(it->getRT(),30.0)
  it = tmp.RTBegin(31.0);
  TEST_REAL_SIMILAR(it->getRT(),40.0)
  TEST_TRUE(tmp.RTBegin(55.0) == tmp.end())
}
END_SECTION

START_SECTION((Iterator RTEnd(CoordinateType rt)))
{
  PeakMap tmp = createPeakMapWithRTs({30, 40, 45, 50});
  PeakMap::Iterator it;

  it = tmp.RTEnd(20.0);
  TEST_REAL_SIMILAR(it->getRT(),30.0)
  it = tmp.RTEnd(30.0);
  TEST_REAL_SIMILAR(it->getRT(),40.0)
  it = tmp.RTEnd(31.0);
  TEST_REAL_SIMILAR(it->getRT(),40.0)
  TEST_TRUE(tmp.RTEnd(55.0) == tmp.end())
}
END_SECTION

START_SECTION(ConstIterator getClosestSpectrumInRT(const double RT) const)
{
  const PeakMap tmp = createPeakMapWithRTs({30, 40, 45, 50});
  PeakMap::ConstIterator it;

  it = tmp.getClosestSpectrumInRT(-200.0);
  TEST_TRUE(it == tmp.cbegin());
  it = tmp.getClosestSpectrumInRT(20.0);
  TEST_TRUE(it == tmp.cbegin());
  it = tmp.getClosestSpectrumInRT(31.0);
  TEST_TRUE(it == tmp.cbegin());
  it = tmp.getClosestSpectrumInRT(34.9);
  TEST_TRUE(it == tmp.cbegin());

  it = tmp.getClosestSpectrumInRT(39);
  TEST_TRUE(it == tmp.cbegin() + 1);
  it = tmp.getClosestSpectrumInRT(41);
  TEST_TRUE(it == tmp.cbegin() + 1);
  it = tmp.getClosestSpectrumInRT(42.4);
  TEST_TRUE(it == tmp.cbegin() + 1);

  it = tmp.getClosestSpectrumInRT(44);
  TEST_TRUE(it == tmp.cbegin() + 2);
  it = tmp.getClosestSpectrumInRT(47);
  TEST_TRUE(it == tmp.cbegin() + 2);

  it = tmp.getClosestSpectrumInRT(47.6);
  TEST_TRUE(it == tmp.cbegin() + 3);
  it = tmp.getClosestSpectrumInRT(51);
  TEST_TRUE(it == tmp.cbegin() + 3);
  it = tmp.getClosestSpectrumInRT(5100000);
  TEST_TRUE(it == tmp.cbegin() + 3);


  const PeakMap tmp_empty;
  it = tmp_empty.getClosestSpectrumInRT(47.6);
  TEST_TRUE(it == tmp_empty.cend());
}
END_SECTION

START_SECTION(Iterator getClosestSpectrumInRT(const double RT))
{
  // minimal version of the above, just to see if iterator types are correct
  PeakMap tmp = createPeakMapWithRTs({30, 40, 45, 50});
  PeakMap::Iterator it;

  it = tmp.getClosestSpectrumInRT(-200.0);
  TEST_TRUE(it == tmp.begin());

  PeakMap tmp_empty;
  it = tmp_empty.getClosestSpectrumInRT(47.6);
  TEST_TRUE(it == tmp_empty.end());
}
END_SECTION


START_SECTION(ConstIterator getClosestSpectrumInRT(const double RT, UInt ms_level) const)
{
  const PeakMap tmp = setMSLevel(createPeakMapWithRTs({30, 31, 32,       40, 41,       50,       60, 61}), {1, 2, 2,     1, 2,     1,      1, 2});
  PeakMap::ConstIterator it;

  it = tmp.getClosestSpectrumInRT(-200.0, 0); // MS-level 0 does not exist --> cend()
  TEST_TRUE(it == tmp.cend());
  it = tmp.getClosestSpectrumInRT(-200.0, 1);
  TEST_TRUE(it == tmp.cbegin());
  it = tmp.getClosestSpectrumInRT(-200.0, 2);
  TEST_TRUE(it == tmp.cbegin() + 1);

  it = tmp.getClosestSpectrumInRT(20.0, 1);
  TEST_TRUE(it == tmp.cbegin());
  it = tmp.getClosestSpectrumInRT(31.0, 1);
  TEST_TRUE(it == tmp.cbegin());
  it = tmp.getClosestSpectrumInRT(34.9, 1);
  TEST_TRUE(it == tmp.cbegin());

  it = tmp.getClosestSpectrumInRT(20.0, 2);
  TEST_TRUE(it == tmp.cbegin() + 1);
  it = tmp.getClosestSpectrumInRT(31.0, 2);
  TEST_TRUE(it == tmp.cbegin() + 1);
  it = tmp.getClosestSpectrumInRT(31.4, 2);
  TEST_TRUE(it == tmp.cbegin() + 1);


  it = tmp.getClosestSpectrumInRT(39, 1);
  TEST_TRUE(it == tmp.cbegin() + 3);
  it = tmp.getClosestSpectrumInRT(41, 1);
  TEST_TRUE(it == tmp.cbegin() + 3);
  it = tmp.getClosestSpectrumInRT(42.4, 1);
  TEST_TRUE(it == tmp.cbegin() + 3);

  it = tmp.getClosestSpectrumInRT(45.5, 1);
  TEST_TRUE(it == tmp.cbegin() + 5);
  it = tmp.getClosestSpectrumInRT(49, 1);
  TEST_TRUE(it == tmp.cbegin() + 5);
  it = tmp.getClosestSpectrumInRT(54.5, 1);
  TEST_TRUE(it == tmp.cbegin() + 5);

  it = tmp.getClosestSpectrumInRT(55.1, 1);
  TEST_TRUE(it == tmp.cbegin() + 6);
  it = tmp.getClosestSpectrumInRT(59.1, 1);
  TEST_TRUE(it == tmp.cbegin() + 6);
  it = tmp.getClosestSpectrumInRT(5100000, 1);
  TEST_TRUE(it == tmp.cbegin() + 6);

  it = tmp.getClosestSpectrumInRT(58, 2);
  TEST_TRUE(it == tmp.cbegin() + 7);
  it = tmp.getClosestSpectrumInRT(63, 2);
  TEST_TRUE(it == tmp.cbegin() + 7);
  it = tmp.getClosestSpectrumInRT(5100000, 2);
  TEST_TRUE(it == tmp.cbegin() + 7);


  const PeakMap tmp_empty;
  it = tmp_empty.getClosestSpectrumInRT(47.6, 1);
  TEST_TRUE(it == tmp_empty.cend());
}
END_SECTION

START_SECTION(Iterator getClosestSpectrumInRT(const double RT, UInt ms_level))
{
  // minimal version of the above, just to see if iterator types are correct
  PeakMap tmp = setMSLevel(createPeakMapWithRTs({30, 31, 32, 40, 41, 50, 60, 61}), {1, 2, 2, 1, 2, 1, 1, 2});
  PeakMap::Iterator it;

  it = tmp.getClosestSpectrumInRT(-200.0, 0); // MS-level 0 does not exist --> cend()
  TEST_TRUE(it == tmp.cend());
  it = tmp.getClosestSpectrumInRT(-200.0, 1);
  TEST_TRUE(it == tmp.cbegin());
  it = tmp.getClosestSpectrumInRT(-200.0, 2);
  TEST_TRUE(it == tmp.cbegin() + 1);
  
  PeakMap tmp_empty;
  it = tmp_empty.getClosestSpectrumInRT(47.6, 1);
  TEST_TRUE(it == tmp_empty.cend());
}
END_SECTION

START_SECTION((Iterator IMBegin(CoordinateType im)))
{
  PeakMap tmp;
  MSSpectrum s;

  s.setDriftTime(30.0);
  tmp.addSpectrum(s);
  s.setDriftTime(40.0);
  tmp.addSpectrum(s);
  s.setDriftTime(45.0);
  tmp.addSpectrum(s);
  s.setDriftTime(50.0);
  tmp.addSpectrum(s);

  PeakMap::ConstIterator it;

  it = tmp.IMBegin(20.0);
  TEST_REAL_SIMILAR(it->getDriftTime(), 30.0)
  it = tmp.IMBegin(30.0);
  TEST_REAL_SIMILAR(it->getDriftTime(), 30.0)
  it = tmp.IMBegin(31.0);
  TEST_REAL_SIMILAR(it->getDriftTime(), 40.0)
  TEST_EQUAL(tmp.IMBegin(55.0) == tmp.end(), true)
}
END_SECTION

START_SECTION((Iterator IMEnd(CoordinateType rt)))
{
  PeakMap tmp;
  MSSpectrum s;

  s.setDriftTime(30.0);
  tmp.addSpectrum(s);
  s.setDriftTime(40.0);
  tmp.addSpectrum(s);
  s.setDriftTime(45.0);
  tmp.addSpectrum(s);
  s.setDriftTime(50.0);
  tmp.addSpectrum(s);

  PeakMap::ConstIterator it;

  it = tmp.IMEnd(20.0);
  TEST_REAL_SIMILAR(it->getDriftTime(), 30.0)
  it = tmp.IMEnd(30.0);
  TEST_REAL_SIMILAR(it->getDriftTime(), 40.0)
  it = tmp.IMEnd(31.0);
  TEST_REAL_SIMILAR(it->getDriftTime(), 40.0)
  TEST_EQUAL(tmp.IMEnd(55.0) == tmp.end(), true)
}
END_SECTION


START_SECTION((ConstIterator RTBegin(CoordinateType rt) const))
{
  PeakMap tmp;
  MSSpectrum s;

  s.setRT(30.0);
  tmp.addSpectrum(s);
  s.setRT(40.0);
  tmp.addSpectrum(s);
  s.setRT(45.0);
  tmp.addSpectrum(s);
  s.setRT(50.0);
  tmp.addSpectrum(s);

  PeakMap::Iterator it;

  it = tmp.RTBegin(20.0);
  TEST_REAL_SIMILAR(it->getRT(),30.0)
  it = tmp.RTBegin(30.0);
  TEST_REAL_SIMILAR(it->getRT(),30.0)
  it = tmp.RTBegin(31.0);
  TEST_REAL_SIMILAR(it->getRT(),40.0)
  TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true)
}
END_SECTION

START_SECTION((ConstIterator RTEnd(CoordinateType rt) const))
{
  PeakMap tmp;
  MSSpectrum s;

  s.setRT(30.0);
  tmp.addSpectrum(s);
  s.setRT(40.0);
  tmp.addSpectrum(s);
  s.setRT(45.0);
  tmp.addSpectrum(s);
  s.setRT(50.0);
  tmp.addSpectrum(s);

  PeakMap::Iterator it;

  it = tmp.RTEnd(20.0);
  TEST_REAL_SIMILAR(it->getRT(),30.0)
  it = tmp.RTEnd(30.0);
  TEST_REAL_SIMILAR(it->getRT(),40.0)
  it = tmp.RTEnd(31.0);
  TEST_REAL_SIMILAR(it->getRT(),40.0)
  TEST_EQUAL(tmp.RTBegin(55.0) == tmp.end(), true)
}
END_SECTION

START_SECTION((void sortSpectra(bool sort_mz = true)))
{
  std::vector< Peak2D> plist;

  Peak2D p1;
  p1.getPosition()[0] = 1.0;
  p1.getPosition()[1] = 5.0;
  plist.push_back(p1);

  Peak2D p2;
  p2.getPosition()[0] = 1.0;
  p2.getPosition()[1] = 3.0;
  plist.push_back(p2);

  Peak2D p3;
  p3.getPosition()[0] = 2.0;
  p3.getPosition()[1] = 14.0;
  plist.push_back(p3);

  Peak2D p4;
  p4.getPosition()[0] = 2.0;
  p4.getPosition()[1] = 11.0;
  plist.push_back(p4);

  PeakMap exp;
  exp.set2DData(plist);

  exp.sortSpectra(true);

  TEST_REAL_SIMILAR(exp[0][0].getMZ(),3.0);
  TEST_REAL_SIMILAR(exp[0][1].getMZ(),5.0);
  TEST_REAL_SIMILAR(exp[1][0].getMZ(),11.0);
  TEST_REAL_SIMILAR(exp[1][1].getMZ(),14.0);
}
END_SECTION

START_SECTION(bool isSorted(bool check_mz = true ) const)
{
  //make test dataset
  PeakMap exp;
  exp.resize(2);
  exp[0].setRT(1.0);
  exp[1].setRT(2.0);

  Peak1D p;
  p.setIntensity(1.0);
  p.setMZ(1000.0);
  exp[0].push_back(p);
  exp[1].push_back(p);

  p.setIntensity(1.0);
  p.setMZ(1001.0);
  exp[0].push_back(p);
  exp[1].push_back(p);

  p.setIntensity(1.0);
  p.setMZ(1002.0);
  exp[0].push_back(p);
  exp[1].push_back(p);

  //test with identical RTs
  TEST_EQUAL(exp.isSorted(false),true)
    TEST_EQUAL(exp.isSorted(),true)

    //test with acending RTs
    exp[0].setRT(1.0);
  exp[1].setRT(2.0);
  TEST_EQUAL(exp.isSorted(false),true)
    TEST_EQUAL(exp.isSorted(),true)

    //test with a reversed spectrum
    reverse(exp[0].begin(),exp[0].end());
  TEST_EQUAL(exp.isSorted(false),true)
    TEST_EQUAL(exp.isSorted(),false)

    //test with reversed RTs
    reverse(exp.begin(),exp.end());
  TEST_EQUAL(exp.isSorted(false),false)
    TEST_EQUAL(exp.isSorted(),false)
}
END_SECTION

START_SECTION((void reset()))
{
  std::vector< Peak2D> plist;

  Peak2D p;
  p.getPosition()[0] = 1.0;
  p.getPosition()[1] = 5.0;
  plist.push_back(p);
  p.getPosition()[0] = 2.0;
  p.getPosition()[1] = 3.0;
  plist.push_back(p);

  PeakMap exp;
  exp.set2DData(plist);
  exp.updateRanges();

  exp.reset();

  TEST_EQUAL(exp.empty(),true);
}
END_SECTION

START_SECTION((const ExperimentalSettings& getExperimentalSettings() const))
{
  PeakMap exp;
  exp.setComment("test");
  TEST_EQUAL(exp.getExperimentalSettings().getComment(),"test");
}
END_SECTION

START_SECTION((ExperimentalSettings& getExperimentalSettings()))
{
  PeakMap exp;
  exp.getExperimentalSettings().setComment("test");
  TEST_EQUAL(exp.getExperimentalSettings().getComment(),"test");
}
END_SECTION

START_SECTION((MSExperiment& operator=(const ExperimentalSettings &source)))
{
  PeakMap exp,exp2;
  exp.getExperimentalSettings().setComment("test");
  exp2 = exp.getExperimentalSettings();
  TEST_EQUAL(exp2.getExperimentalSettings().getComment(),"test");
}
END_SECTION

START_SECTION((ConstIterator getPrecursorSpectrum(ConstIterator iterator) const))
{
  PeakMap exp;
  exp.resize(10);
  exp[0].setMSLevel(1);
  exp[1].setMSLevel(2);
  exp[2].setMSLevel(1);
  exp[3].setMSLevel(2);
  exp[4].setMSLevel(2);

  TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin())==exp.end(),true)
  TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+1)==exp.begin(),true)
  TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+2)==exp.end(),true)
  TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+3)==exp.begin()+2,true)
  TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+4)==exp.begin()+2,true)
  TEST_EQUAL(exp.getPrecursorSpectrum(exp.end())==exp.end(),true)

  exp[0].setMSLevel(2);
  exp[1].setMSLevel(1);
  exp[2].setMSLevel(1);
  exp[3].setMSLevel(1);
  exp[4].setMSLevel(1);

  TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin())==exp.end(),true)
  TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+1)==exp.end(),true)
  TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+2)==exp.end(),true)
  TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+3)==exp.end(),true)
  TEST_EQUAL(exp.getPrecursorSpectrum(exp.begin()+4)==exp.end(),true)
  TEST_EQUAL(exp.getPrecursorSpectrum(exp.end())==exp.end(),true)
}
END_SECTION

START_SECTION((int getPrecursorSpectrum(int zero_based_index) const))
{
  PeakMap exp;
  exp.resize(10);
  exp[0].setMSLevel(1);
  exp[1].setMSLevel(2);
  exp[2].setMSLevel(1);
  exp[3].setMSLevel(2);
  exp[4].setMSLevel(2);

  TEST_EQUAL(exp.getPrecursorSpectrum(0) == -1,true)
  TEST_EQUAL(exp.getPrecursorSpectrum(1) == 0,true)
  TEST_EQUAL(exp.getPrecursorSpectrum(2) == -1,true)
  TEST_EQUAL(exp.getPrecursorSpectrum(3) == 2,true)
  TEST_EQUAL(exp.getPrecursorSpectrum(4) == 2,true)

  exp[0].setMSLevel(2);
  exp[1].setMSLevel(1);
  exp[2].setMSLevel(1);
  exp[3].setMSLevel(1);
  exp[4].setMSLevel(1);

  TEST_EQUAL(exp.getPrecursorSpectrum(0) == -1,true)
  TEST_EQUAL(exp.getPrecursorSpectrum(1) == -1,true)
  TEST_EQUAL(exp.getPrecursorSpectrum(2) == -1,true)
  TEST_EQUAL(exp.getPrecursorSpectrum(3) == -1,true)
  TEST_EQUAL(exp.getPrecursorSpectrum(4) == -1,true)
}
END_SECTION

START_SECTION((ConstIterator getFirstProductSpectrum(ConstIterator iterator) const))
{
  PeakMap exp;
  exp.resize(5);

  // Set MS levels
  exp[0].setMSLevel(1);
  exp[1].setMSLevel(2);
  exp[2].setMSLevel(1);
  exp[3].setMSLevel(2);
  exp[4].setMSLevel(2);

  // Set NativeIDs
  exp[0].setNativeID("scan=1");
  exp[1].setNativeID("scan=2");
  exp[2].setNativeID("scan=3");
  exp[3].setNativeID("scan=4");
  exp[4].setNativeID("scan=5");

  // Set up the 'spectrum_ref' in the precursors of the product spectra
  Precursor precursor1;
  precursor1.setMetaValue("spectrum_ref", "scan=1"); // Reference to exp[0]
  exp[1].getPrecursors().push_back(precursor1);

  Precursor precursor2;
  precursor2.setMetaValue("spectrum_ref", "scan=3"); // Reference to exp[2]
  exp[3].getPrecursors().push_back(precursor2);

  Precursor precursor3;
  precursor3.setMetaValue("spectrum_ref", "scan=3"); // Another reference to exp[2]
  exp[4].getPrecursors().push_back(precursor3);

  // Test getFirstProductSpectrum

  // From exp[0], expect to get exp[1] as the first product spectrum
  TEST_EQUAL(exp.getFirstProductSpectrum(exp.begin()) == exp.begin() + 1, true)

  // From exp[1], expect to get spectra_.end() since there is no higher MS level
  TEST_EQUAL(exp.getFirstProductSpectrum(exp.begin() + 1) == exp.end(), true)

  // From exp[2], expect to get exp[3] as the first product spectrum
  TEST_EQUAL(exp.getFirstProductSpectrum(exp.begin() + 2) == exp.begin() + 3, true)

  // From exp[3], expect to get spectra_.end() since there is no higher MS level
  TEST_EQUAL(exp.getFirstProductSpectrum(exp.begin() + 3) == exp.end(), true)

  // From exp[4], expect to get spectra_.end()
  TEST_EQUAL(exp.getFirstProductSpectrum(exp.begin() + 4) == exp.end(), true)

  // Test when iterator is spectra_.end()
  TEST_EQUAL(exp.getFirstProductSpectrum(exp.end()) == exp.end(), true)

  // Now change the MS levels to only MS1 spectra
  for (Size i = 0; i < exp.size(); ++i)
  {
    exp[i].setMSLevel(1);
  }

  // Test again with only MS1 spectra
  TEST_EQUAL(exp.getFirstProductSpectrum(exp.begin()) == exp.end(), true)
  TEST_EQUAL(exp.getFirstProductSpectrum(exp.begin() + 1) == exp.end(), true)
  TEST_EQUAL(exp.getFirstProductSpectrum(exp.begin() + 2) == exp.end(), true)
  TEST_EQUAL(exp.getFirstProductSpectrum(exp.begin() + 3) == exp.end(), true)
  TEST_EQUAL(exp.getFirstProductSpectrum(exp.begin() + 4) == exp.end(), true)
}
END_SECTION

START_SECTION((int getFirstProductSpectrum(int zero_based_index) const))
{
  PeakMap exp;
  exp.resize(5);

  // Set MS levels
  exp[0].setMSLevel(1);
  exp[1].setMSLevel(2);
  exp[2].setMSLevel(1);
  exp[3].setMSLevel(2);
  exp[4].setMSLevel(2);

  // Set NativeIDs
  exp[0].setNativeID("scan=1");
  exp[1].setNativeID("scan=2");
  exp[2].setNativeID("scan=3");
  exp[3].setNativeID("scan=4");
  exp[4].setNativeID("scan=5");

  // Set up the 'spectrum_ref' in the precursors of the product spectra
  Precursor precursor1;
  precursor1.setMetaValue("spectrum_ref", "scan=1"); // Reference to exp[0]
  exp[1].getPrecursors().push_back(precursor1);

  Precursor precursor2;
  precursor2.setMetaValue("spectrum_ref", "scan=3"); // Reference to exp[2]
  exp[3].getPrecursors().push_back(precursor2);

  Precursor precursor3;
  precursor3.setMetaValue("spectrum_ref", "scan=3"); // Another reference to exp[2]
  exp[4].getPrecursors().push_back(precursor3);

  // Test getFirstProductSpectrum

  // From index 0, expect to get index 1
  TEST_EQUAL(exp.getFirstProductSpectrum(0) == 1, true)

  // From index 1, expect to get -1 (no higher MS level with reference)
  TEST_EQUAL(exp.getFirstProductSpectrum(1) == -1, true)

  // From index 2, expect to get index 3
  TEST_EQUAL(exp.getFirstProductSpectrum(2) == 3, true)

  // From index 3, expect to get -1
  TEST_EQUAL(exp.getFirstProductSpectrum(3) == -1, true)

  // From index 4, expect to get -1
  TEST_EQUAL(exp.getFirstProductSpectrum(4) == -1, true)

  // Now change the MS levels to only MS1 spectra
  for (Size i = 0; i < exp.size(); ++i)
  {
    exp[i].setMSLevel(1);
  }

  // Test again with only MS1 spectra
  TEST_EQUAL(exp.getFirstProductSpectrum(0) == -1, true)
  TEST_EQUAL(exp.getFirstProductSpectrum(1) == -1, true)
  TEST_EQUAL(exp.getFirstProductSpectrum(2) == -1, true)
  TEST_EQUAL(exp.getFirstProductSpectrum(3) == -1, true)
  TEST_EQUAL(exp.getFirstProductSpectrum(4) == -1, true)
}
END_SECTION

START_SECTION((bool clearMetaDataArrays()))
{
  PeakMap exp;
  exp.resize(5);
  exp[0].getFloatDataArrays().resize(5);
  exp[0].getIntegerDataArrays().resize(5);
  exp[0].getStringDataArrays().resize(5);
  exp.clearMetaDataArrays();
  TEST_EQUAL(exp[0].getFloatDataArrays().size(),0)
    TEST_EQUAL(exp[0].getIntegerDataArrays().size(),0)
    TEST_EQUAL(exp[0].getStringDataArrays().size(),0)
}
END_SECTION

START_SECTION((void swap(MSExperiment &from)))
{
  PeakMap exp1, exp2;
  exp1.setComment("stupid comment");
  exp1.resize(1);
  exp1[0].setMSLevel(2);
  exp1[0].resize(2);
  exp1[0][0].setIntensity(0.5f);
  exp1[0][1].setIntensity(1.7f);
  exp1.updateRanges();

  exp1.swap(exp2);

  TEST_EQUAL(exp1.getComment(),"")
  TEST_EQUAL(exp1.size(),0)
  TEST_EQUAL(exp1.getRange().hasRange() == HasRangeType::NONE, true)
  TEST_EQUAL(exp1.getMSLevels().size(),0)
  TEST_EQUAL(exp1.getSize(),0);

  TEST_EQUAL(exp2.getComment(),"stupid comment")
  TEST_EQUAL(exp2.size(),1)
  TEST_REAL_SIMILAR(exp2.getMinIntensity(), 0.5)
  TEST_EQUAL(exp2.getMSLevels().size(),1)
  TEST_EQUAL(exp2.getSize(),2);
}
END_SECTION

START_SECTION(void clear(bool clear_meta_data))
{
  PeakMap edit;
  edit.getSample().setName("bla");
  edit.resize(5);
  edit.updateRanges();
  edit.setMetaValue("label",String("bla"));
  vector<MSChromatogram > tmp;
  tmp.resize(5);
  edit.setChromatograms(tmp);

  edit.clear(false);
  TEST_EQUAL(edit.size(),0)
  TEST_EQUAL(edit == MSExperiment(),false)

  edit.clear(true);
  TEST_EQUAL(edit.empty(),true)
  TEST_EQUAL(edit == MSExperiment(),true)
}
END_SECTION

START_SECTION(bool MSExperiment::isIMFrame() const)
{
  PeakMap tmp;
  MSSpectrum s;
  constexpr double rt_all = 1;
  s.setRT(rt_all);
  s.setDriftTime(30.0);
  tmp.addSpectrum(s);
  s.setDriftTime(40.0);
  tmp.addSpectrum(s);
  s.setDriftTime(45.0);
  tmp.addSpectrum(s);
  s.setDriftTime(50.0);
  tmp.addSpectrum(s);
  TEST_TRUE(tmp.isIMFrame())

  tmp[3].setRT(2);  // changing RT ...
  TEST_FALSE(tmp.isIMFrame());
  tmp[3].setRT(rt_all); // undo

  tmp[2].setDriftTime(tmp[1].getDriftTime()); // duplicate drift time
  TEST_FALSE(tmp.isIMFrame());

  tmp[3].setRT(2); // changing RT ...
  tmp[2].setDriftTime(tmp[1].getDriftTime()); // duplicate drift time
  TEST_FALSE(tmp.isIMFrame());
}
END_SECTION

START_SECTION((void sortChromatograms(bool sort_rt=true)))
{
  PeakMap exp;
  MSChromatogram chrom1, chrom2;
  ChromatogramPeak p1, p2, p3;
  p1.setRT(0.3);
  p1.setIntensity(10.0f);
  p2.setRT(0.2);
  p2.setIntensity(10.2f);
  p3.setRT(0.1);
  p3.setIntensity(10.4f);

  Product prod1;
  prod1.setMZ(100.0);
  chrom1.setProduct(prod1);
  chrom1.push_back(p1);
  chrom1.push_back(p2);

  Product prod2;
  prod2.setMZ(80.0);
  chrom2.setProduct(prod2);
  chrom2.push_back(p2);
  chrom2.push_back(p3);

  vector<MSChromatogram > chroms;
  chroms.push_back(chrom1);
  chroms.push_back(chrom2);
  exp.setChromatograms(chroms);
  TEST_EQUAL(exp.getChromatograms().size(), 2)
    TEST_REAL_SIMILAR(exp.getChromatograms()[0].getMZ(), 100.0)
    TEST_REAL_SIMILAR(exp.getChromatograms()[1].getMZ(), 80.0)

    // first sort without rt
    exp.sortChromatograms(false);
  TEST_REAL_SIMILAR(exp.getChromatograms()[0].getMZ(), 80.0)
    TEST_REAL_SIMILAR(exp.getChromatograms()[1].getMZ(), 100.0)

    TEST_REAL_SIMILAR(exp.getChromatograms()[1][0].getRT(), 0.3)
    TEST_REAL_SIMILAR(exp.getChromatograms()[1][1].getRT(), 0.2)

    // now also sort rt
    exp.sortChromatograms();

  TEST_REAL_SIMILAR(exp.getChromatograms()[0].getMZ(), 80.0)
    TEST_REAL_SIMILAR(exp.getChromatograms()[1].getMZ(), 100.0)

    TEST_REAL_SIMILAR(exp.getChromatograms()[1][0].getRT(), 0.2)
    TEST_REAL_SIMILAR(exp.getChromatograms()[1][1].getRT(), 0.3)

}
END_SECTION

START_SECTION((void setChromatograms(const std::vector< MSChromatogram > &chromatograms)))
{
  PeakMap exp;
  MSChromatogram chrom1, chrom2;
  ChromatogramPeak p1, p2, p3;
  p1.setRT(0.1);
  p1.setIntensity(10.0f);
  p2.setRT(0.2);
  p2.setIntensity(10.2f);
  p3.setRT(0.3);
  p3.setIntensity(10.4f);
  chrom1.push_back(p1);
  chrom1.push_back(p2);
  chrom2.push_back(p2);
  chrom2.push_back(p3);
  vector<MSChromatogram > chroms;
  chroms.push_back(chrom1);
  chroms.push_back(chrom2);
  exp.setChromatograms(chroms);
  TEST_EQUAL(exp.getChromatograms().size(), 2)
    TEST_EQUAL(exp.getChromatograms()[0] == chrom1, true)
    TEST_EQUAL(exp.getChromatograms()[1] == chrom2, true)
}
END_SECTION

START_SECTION((void addChromatogram(const MSChromatogram &chromatogram)))
{
  PeakMap exp;
  MSChromatogram chrom1, chrom2;
  ChromatogramPeak p1, p2, p3;
  p1.setRT(0.1);
  p1.setIntensity(10.0f);
  p2.setRT(0.2);
  p2.setIntensity(10.2f);
  p3.setRT(0.3);
  p3.setIntensity(10.4f);
  chrom1.push_back(p1);
  chrom1.push_back(p2);
  chrom2.push_back(p2);
  chrom2.push_back(p3);

  TEST_EQUAL(exp.getChromatograms().size(), 0)
    exp.addChromatogram(chrom1);
  TEST_EQUAL(exp.getChromatograms().size(), 1)
    TEST_EQUAL(exp.getChromatograms()[0] == chrom1, true)
    exp.addChromatogram(chrom2);
  TEST_EQUAL(exp.getChromatograms().size(), 2)
    TEST_EQUAL(exp.getChromatograms()[0] == chrom1, true)  
    TEST_EQUAL(exp.getChromatograms()[1] == chrom2, true)  
}
END_SECTION

START_SECTION((const std::vector<MSChromatogram >& getChromatograms() const))
{
  NOT_TESTABLE // tested above
}
END_SECTION

START_SECTION((std::vector<MSChromatogram >& getChromatograms()))
{
  PeakMap exp;
  vector<MSChromatogram > chromatograms(2);
  exp.getChromatograms().swap(chromatograms);
  TEST_EQUAL(exp.getChromatograms().size(), 2);
  TEST_EQUAL(chromatograms.size(), 0);
  exp.getChromatograms().swap(chromatograms);
  TEST_EQUAL(exp.getChromatograms().size(), 0);
  TEST_EQUAL(chromatograms.size(), 2);
}
END_SECTION

START_SECTION((const MSChromatogram calculateTIC(float rt_bin_size=0) const))
{
	MSChromatogram chrom;
  // Dummy peakmap
	PeakMap exp;
	exp.resize(4);
	Peak1D p;

	// MS spectrum at RT = 0
	p.setMZ(5.0);
	p.setIntensity(3);
	exp[0].push_back(p);
	p.setMZ(10.0);
	p.setIntensity(5);
	exp[0].push_back(p);
	exp[0].setMSLevel(1);
	exp[0].setRT(0);

	// MS spectrum at RT = 2
	p.setMZ(5.0);
	p.setIntensity(2);
	exp[1].push_back(p);
	exp[1].setMSLevel(1);
	exp[1].setRT(2);

	// MSMS spectrum at RT = 2
	p.setMZ(5.0);
	p.setIntensity(0.5);
	exp[2].push_back(p);
	exp[2].setMSLevel(2);
	exp[2].setRT(2);

	// MS spectrum at RT = 5
	p.setMZ(5.0);
	p.setIntensity(2.0);
	exp[3].push_back(p);
	p.setMZ(10.0);
	p.setIntensity(3.0);
	exp[3].push_back(p);
	p.setMZ(15.0);
	p.setIntensity(4.0);
	exp[3].push_back(p);
	exp[3].setMSLevel(1);
	exp[3].setRT(5);

	exp.updateRanges();

	// empty MSExperiment
	MSExperiment exp2;
	chrom = exp2.calculateTIC();
	TEST_EQUAL(chrom.empty(),true);

	// no binning
	chrom = exp.calculateTIC();
	ABORT_IF(chrom.size() != 3);
	TEST_EQUAL(chrom[0].getIntensity(),8);
	TEST_EQUAL(chrom[1].getIntensity(),2);
	TEST_EQUAL(chrom[2].getIntensity(),9);

	// bin size smaller than highest RT
	chrom = exp.calculateTIC(2.0);
	ABORT_IF(chrom.size() != 4);
	TEST_EQUAL(chrom[0].getIntensity(),8);
	TEST_EQUAL(chrom[1].getIntensity(),2);
	// Intensity at RT = 5 in between new data points at 4.0 and 6.0
	TEST_EQUAL(chrom[2].getIntensity(),4.5);
	TEST_EQUAL(chrom[3].getIntensity(),4.5);

	// bin size bigger than highest RT
	chrom = exp.calculateTIC(6.0);
	ABORT_IF(chrom.size() != 2);
	// Intensities at RT = 2 and RT = 5 in between new data points at 0.0 and 6.0
	TEST_REAL_SIMILAR(chrom[0].getIntensity(),8.0 + 2.0* 4.0/6.0 + 9 * 1.0/6.0);
	TEST_REAL_SIMILAR(chrom[1].getIntensity(),2.0* 2.0/6.0 + 9 * 5.0/6.0);

	// negative bin size
	chrom = exp.calculateTIC(-1.0);
	// should be like no bin size was given
	ABORT_IF(chrom.size() != 3);
	TEST_EQUAL(chrom[0].getIntensity(),8);
	TEST_EQUAL(chrom[1].getIntensity(),2);
	TEST_EQUAL(chrom[2].getIntensity(),9);
}
END_SECTION

START_SECTION( std::ostream& operator<<(std::ostream& os, const MSExperiment& chrom)) 
{
  PeakMap tmp;
  tmp.getContacts().resize(1);
  tmp.getContacts()[0].setFirstName("Name");
  tmp.resize(1);
  Peak1D p;
  p.setMZ(5.0);
  tmp[0].push_back(p);
  p.setMZ(10.77);
  tmp[0].push_back(p);

  MSChromatogram a;
  MSChromatogram::PeakType peak;
  peak.getPosition()[0] = 47.11;
  a.push_back(peak);
  tmp.addChromatogram(a);

  std::ostringstream os;
  os << tmp;

  TEST_EQUAL(String(os.str()).hasSubstring("MSEXPERIMENT BEGIN"), true);
  TEST_EQUAL(String(os.str()).hasSubstring("MSSPECTRUM BEGIN"), true);
  TEST_EQUAL(String(os.str()).hasSubstring("MSCHROMATOGRAM BEGIN"), true);
  TEST_EQUAL(String(os.str()).hasSubstring("47.11"), true);
  TEST_EQUAL(String(os.str()).hasSubstring("10.77"), true);
}
END_SECTION

START_SECTION((template<class MzReductionFunctionType> std::vector<std::vector<MSExperiment::CoordinateType>> aggregate(const std::vector<std::pair<RangeMZ, RangeRT>>& mz_rt_ranges, unsigned int ms_level, MzReductionFunctionType func_mz_reduction) const))
{
    // Create test experiment with known data
    PeakMap exp;
    exp.resize(4);

    // First spectrum (MS1) at RT=1.0
    exp[0] = MSSpectrum{
        {100.0, 1000.0},
        {200.0, 2000.0},
        {300.0, 3000.0}
    };    
    exp[0].setRT(1.0);
    exp[0].setMSLevel(1);

    // Second spectrum (MS2) at RT=2.0
    exp[1] = MSSpectrum{
      {150.0, 1500.0},
      {250.0, 2500.0}
    };
    exp[1].setRT(2.0);
    exp[1].setMSLevel(2);

    // Third spectrum (MS1) at RT=3.0
    exp[2] = MSSpectrum{
        {100.0, 1100.0},
        {200.0, 2100.0},
        {300.0, 3100.0}
    };
    exp[2].setRT(3.0);
    exp[2].setMSLevel(1);

    // Fourth spectrum (MS1) at RT=4.0
    exp[3] = MSSpectrum{
        {100.0, 1200.0},
        {200.0, 2200.0},
        {300.0, 3200.0}
    };
    exp[3].setRT(4.0);
    exp[3].setMSLevel(1);

    exp.updateRanges();

    // Test: Normal case - MS1 spectra
    {
        std::vector<std::pair<RangeMZ, RangeRT>> ranges;
        // Range 1: covers first peak of first and third spectrum
        ranges.push_back(std::make_pair(
            RangeMZ(90.0, 110.0),
            RangeRT(0.0, 3.5)
        ));
        // Range 2: covers second peak of all MS1 spectra
        ranges.push_back(std::make_pair(
            RangeMZ(190.0, 210.0),
            RangeRT(0.0, 5.0)
        ));

        // Simple intensity reduction function
        auto result = exp.aggregate(ranges, 1, 
          [](MSSpectrum::ConstIterator begin_it, MSSpectrum::ConstIterator /*end_it*/)->double // return first intensity of peaks in m/z range 
          { 
            return begin_it->getIntensity();
          });

        // Check results
        TEST_EQUAL(result.size(), 2);
        
        // Check Range 1 results
        TEST_EQUAL(result[0].size(), 2);  // Should cover 2 spectra
        TEST_EQUAL(result[0][0], 1000.0); // First spectrum intensity
        TEST_EQUAL(result[0][1], 1100.0); // Third spectrum intensity
        
        // Check Range 2 results
        TEST_EQUAL(result[1].size(), 3);   // Should cover 3 spectra
        TEST_EQUAL(result[1][0], 2000.0);  // First spectrum intensity
        TEST_EQUAL(result[1][1], 2100.0);  // Third spectrum intensity
        TEST_EQUAL(result[1][2], 2200.0);  // Fourth spectrum intensity
    }

    // Test 4: MS2 spectra
    {
        std::vector<std::pair<RangeMZ, RangeRT>> ranges;
        ranges.push_back(std::make_pair(
            RangeMZ(140.0, 160.0),
            RangeRT(1.5, 2.5)
        ));

        auto result = exp.aggregate(ranges, 2, 
          [](const MSSpectrum::ConstIterator begin_it, const MSSpectrum::ConstIterator /*end_it*/)->double // return first intensity of peaks in m/z range 
          { 
            return begin_it->getIntensity();
          });

        TEST_EQUAL(result.size(), 1);
        TEST_EQUAL(result[0].size(), 1);
        TEST_EQUAL(result[0][0], 1500.0);
    }

    // Test 5: Complex reduction function (average intensity)
    {
        std::vector<std::pair<RangeMZ, RangeRT>> ranges;
        ranges.push_back(std::make_pair(
            RangeMZ(90.0, 310.0),  // Covers all peaks
            RangeRT(0.0, 5.0)      // Covers all spectra
        ));

        // mean intensity in m/z range
        auto result = exp.aggregate(ranges, 1, 
            [](const MSSpectrum::ConstIterator begin_it, const MSSpectrum::ConstIterator end_it) 
            { 
                if (begin_it == end_it) return 0.0; // Check for empty range before accumulation

                double acc = std::accumulate(begin_it, end_it, 0.0, 
                    [](double a, const Peak1D& b) { return a + b.getIntensity(); });

                return acc / static_cast<double>(std::distance(begin_it, end_it));
            });
        TEST_EQUAL(result.size(), 1);
        TEST_EQUAL(result[0].size(), 3);
        TEST_REAL_SIMILAR(result[0][0], 2000.0);  // Average of first spectrum
        TEST_REAL_SIMILAR(result[0][1], 2100.0);  // Average of third spectrum
        TEST_REAL_SIMILAR(result[0][2], 2200.0);  // Average of fourth spectrum
    }
}
END_SECTION

START_SECTION((void get2DPeakDataPerSpectrum(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, Size ms_level, std::vector<float>& rt, std::vector<std::vector<float>>& mz, std::vector<std::vector<float>>& intensity) const))
{
  MSExperiment exp;
  
  // Create test spectra using initializer lists
  MSSpectrum s1{
    {100.0, 1000.0},
    {200.0, 2000.0}
  };
  s1.setRT(1.0);
  s1.setMSLevel(1);
  
  MSSpectrum s2{
    {150.0, 1500.0},
    {250.0, 2500.0}
  };
  s2.setRT(2.0);
  s2.setMSLevel(1);
  
  MSSpectrum s3{
    {175.0, 1750.0}
  };
  s3.setRT(3.0);
  s3.setMSLevel(2);
  
  exp.addSpectrum(s1);
  exp.addSpectrum(s2);
  exp.addSpectrum(s3);
  
  // Test 1: Full range, MS level 1
  {
    std::vector<float> rt;
    std::vector<std::vector<float>> mz, intensity;
    exp.get2DPeakDataPerSpectrum(0.0, 4.0, 0.0, 300.0, 1, rt, mz, intensity);
    
    TEST_EQUAL(rt.size(), 2)
    TEST_EQUAL(mz.size(), 2)
    TEST_EQUAL(intensity.size(), 2)
    
    // Check first spectrum
    TEST_REAL_SIMILAR(rt[0], 1.0)
    TEST_EQUAL(mz[0].size(), 2)
    TEST_REAL_SIMILAR(mz[0][0], 100.0)
    TEST_REAL_SIMILAR(mz[0][1], 200.0)
    TEST_REAL_SIMILAR(intensity[0][0], 1000.0)
    TEST_REAL_SIMILAR(intensity[0][1], 2000.0)
    
    // Check second spectrum
    TEST_REAL_SIMILAR(rt[1], 2.0)
    TEST_EQUAL(mz[1].size(), 2)
    TEST_REAL_SIMILAR(mz[1][0], 150.0)
    TEST_REAL_SIMILAR(mz[1][1], 250.0)
    TEST_REAL_SIMILAR(intensity[1][0], 1500.0)
    TEST_REAL_SIMILAR(intensity[1][1], 2500.0)
  }
  
  // Test 2: Limited RT range
  {
    std::vector<float> rt;
    std::vector<std::vector<float>> mz, intensity;
    exp.get2DPeakDataPerSpectrum(1.5, 2.5, 0.0, 300.0, 1, rt, mz, intensity);
    
    TEST_EQUAL(rt.size(), 1)
    TEST_EQUAL(mz.size(), 1)
    TEST_EQUAL(intensity.size(), 1)
    TEST_REAL_SIMILAR(rt[0], 2.0)
  }
  
  // Test 3: Limited MZ range
  {
    std::vector<float> rt;
    std::vector<std::vector<float>> mz, intensity;
    exp.get2DPeakDataPerSpectrum(0.0, 4.0, 120.0, 180.0, 1, rt, mz, intensity);
    
    TEST_EQUAL(rt.size(), 1)
    TEST_EQUAL(mz.size(), 1)
    TEST_EQUAL(intensity.size(), 1)
    TEST_REAL_SIMILAR(rt[0], 2.0)
    TEST_EQUAL(mz[0].size(), 1)
    TEST_REAL_SIMILAR(mz[0][0], 150.0)
  }
  
  // Test 4: MS level 2
  {
    std::vector<float> rt;
    std::vector<std::vector<float>> mz, intensity;
    exp.get2DPeakDataPerSpectrum(0.0, 4.0, 0.0, 300.0, 2, rt, mz, intensity);
    
    TEST_EQUAL(rt.size(), 1)
    TEST_EQUAL(mz.size(), 1)
    TEST_EQUAL(intensity.size(), 1)
    TEST_REAL_SIMILAR(rt[0], 3.0)
    TEST_EQUAL(mz[0].size(), 1)
    TEST_REAL_SIMILAR(mz[0][0], 175.0)
    TEST_REAL_SIMILAR(intensity[0][0], 1750.0)
  }
  
  // Test 5: Empty range
  {
    std::vector<float> rt;
    std::vector<std::vector<float>> mz, intensity;
    exp.get2DPeakDataPerSpectrum(5.0, 6.0, 0.0, 300.0, 1, rt, mz, intensity);
    
    TEST_EQUAL(rt.empty(), true)
    TEST_EQUAL(mz.empty(), true)
    TEST_EQUAL(intensity.empty(), true)
  }
}
END_SECTION

START_SECTION((void get2DPeakDataIMPerSpectrum(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, Size ms_level, std::vector<float>& rt, std::vector<std::vector<float>>& mz, std::vector<std::vector<float>>& intensity, std::vector<std::vector<float>>& ion_mobility) const))
{
  MSExperiment exp;
  
  // Create test spectra with ion mobility data
  MSSpectrum s1{
    {100.0, 1000.0},
    {200.0, 2000.0}
  };
  s1.setRT(1.0);
  s1.setMSLevel(1);
  DataArrays::FloatDataArray im1;
  im1.setName("Ion Mobility");
  im1.push_back(0.8f);  // IM value for first peak
  im1.push_back(1.2f);  // IM value for second peak
  im1.setMetaValue("unit", "millisecond");
  s1.getFloatDataArrays().push_back(im1);
  
  MSSpectrum s2{
    {150.0, 1500.0},
    {250.0, 2500.0}
  };
  s2.setRT(2.0);
  s2.setMSLevel(1);
  DataArrays::FloatDataArray im2;
  im2.setName("Ion Mobility");
  im2.push_back(1.5f);  // IM value for first peak
  im2.push_back(2.0f);  // IM value for second peak
  im2.setMetaValue("unit", "millisecond");
  s2.getFloatDataArrays().push_back(im2);
  
  MSSpectrum s3{
    {175.0, 1750.0}
  };
  s3.setRT(3.0);
  s3.setMSLevel(2);
  DataArrays::FloatDataArray im3;
  im3.setName("Ion Mobility");
  im3.push_back(1.8f);  // IM value for single peak
  im3.setMetaValue("unit", "millisecond");
  s3.getFloatDataArrays().push_back(im3);
  
  exp.addSpectrum(s1);
  exp.addSpectrum(s2);
  exp.addSpectrum(s3);
  
  // Test 1: Full range, MS level 1
  {
    std::vector<float> rt;
    std::vector<std::vector<float>> mz, intensity, ion_mobility;
    exp.get2DPeakDataIMPerSpectrum(0.0, 4.0, 0.0, 300.0, 1, rt, mz, intensity, ion_mobility);
    
    TEST_EQUAL(rt.size(), 2)
    TEST_EQUAL(mz.size(), 2)
    TEST_EQUAL(intensity.size(), 2)
    TEST_EQUAL(ion_mobility.size(), 2)
    
    // Check first spectrum
    TEST_REAL_SIMILAR(rt[0], 1.0)
    TEST_EQUAL(mz[0].size(), 2)
    TEST_REAL_SIMILAR(mz[0][0], 100.0)
    TEST_REAL_SIMILAR(mz[0][1], 200.0)
    TEST_REAL_SIMILAR(intensity[0][0], 1000.0)
    TEST_REAL_SIMILAR(intensity[0][1], 2000.0)
    TEST_REAL_SIMILAR(ion_mobility[0][0], 0.8)
    TEST_REAL_SIMILAR(ion_mobility[0][1], 1.2)
    
    // Check second spectrum
    TEST_REAL_SIMILAR(rt[1], 2.0)
    TEST_EQUAL(mz[1].size(), 2)
    TEST_REAL_SIMILAR(mz[1][0], 150.0)
    TEST_REAL_SIMILAR(mz[1][1], 250.0)
    TEST_REAL_SIMILAR(intensity[1][0], 1500.0)
    TEST_REAL_SIMILAR(intensity[1][1], 2500.0)
    TEST_REAL_SIMILAR(ion_mobility[1][0], 1.5)
    TEST_REAL_SIMILAR(ion_mobility[1][1], 2.0)
  }
  
  // Test 2: Limited RT range
  {
    std::vector<float> rt;
    std::vector<std::vector<float>> mz, intensity, ion_mobility;
    exp.get2DPeakDataIMPerSpectrum(1.5, 2.5, 0.0, 300.0, 1, rt, mz, intensity, ion_mobility);
    
    TEST_EQUAL(rt.size(), 1)
    TEST_EQUAL(mz.size(), 1)
    TEST_EQUAL(intensity.size(), 1)
    TEST_EQUAL(ion_mobility.size(), 1)
    TEST_REAL_SIMILAR(rt[0], 2.0)
    TEST_EQUAL(ion_mobility[0].size(), 2)
    TEST_REAL_SIMILAR(ion_mobility[0][0], 1.5)
    TEST_REAL_SIMILAR(ion_mobility[0][1], 2.0)
  }
  
  // Test 3: Limited MZ range
  {
    std::vector<float> rt;
    std::vector<std::vector<float>> mz, intensity, ion_mobility;
    exp.get2DPeakDataIMPerSpectrum(0.0, 4.0, 120.0, 180.0, 1, rt, mz, intensity, ion_mobility);
    
    TEST_EQUAL(rt.size(), 1)
    TEST_EQUAL(mz.size(), 1)
    TEST_EQUAL(intensity.size(), 1)
    TEST_EQUAL(ion_mobility.size(), 1)
    TEST_REAL_SIMILAR(rt[0], 2.0)
    TEST_EQUAL(mz[0].size(), 1)
    TEST_REAL_SIMILAR(mz[0][0], 150.0)
    TEST_REAL_SIMILAR(ion_mobility[0][0], 1.5)
  }
  
  // Test 4: MS level 2
  {
    std::vector<float> rt;
    std::vector<std::vector<float>> mz, intensity, ion_mobility;
    exp.get2DPeakDataIMPerSpectrum(0.0, 4.0, 0.0, 300.0, 2, rt, mz, intensity, ion_mobility);
    
    TEST_EQUAL(rt.size(), 1)
    TEST_EQUAL(mz.size(), 1)
    TEST_EQUAL(intensity.size(), 1)
    TEST_EQUAL(ion_mobility.size(), 1)
    TEST_REAL_SIMILAR(rt[0], 3.0)
    TEST_EQUAL(mz[0].size(), 1)
    TEST_REAL_SIMILAR(mz[0][0], 175.0)
    TEST_REAL_SIMILAR(intensity[0][0], 1750.0)
    TEST_REAL_SIMILAR(ion_mobility[0][0], 1.8)
  }
  
  // Test 5: Empty range
  {
    std::vector<float> rt;
    std::vector<std::vector<float>> mz, intensity, ion_mobility;
    exp.get2DPeakDataIMPerSpectrum(5.0, 6.0, 0.0, 300.0, 1, rt, mz, intensity, ion_mobility);
    
    TEST_EQUAL(rt.empty(), true)
    TEST_EQUAL(mz.empty(), true)
    TEST_EQUAL(intensity.empty(), true)
    TEST_EQUAL(ion_mobility.empty(), true)
  }

  // Test 6: Spectrum without ion mobility data
  {
    MSExperiment exp_no_im;
    MSSpectrum s_no_im{
      {100.0, 1000.0},
      {200.0, 2000.0}
    };
    s_no_im.setRT(1.0);
    s_no_im.setMSLevel(1);
    exp_no_im.addSpectrum(s_no_im);
    
    std::vector<float> rt;
    std::vector<std::vector<float>> mz, intensity, ion_mobility;
    exp_no_im.get2DPeakDataIMPerSpectrum(0.0, 4.0, 0.0, 300.0, 1, rt, mz, intensity, ion_mobility);
    
    TEST_EQUAL(rt.size(), 1)
    TEST_EQUAL(mz.size(), 1)
    TEST_EQUAL(intensity.size(), 1)
    TEST_EQUAL(ion_mobility.size(), 1)
    TEST_EQUAL(ion_mobility[0].size(), 2)
    TEST_REAL_SIMILAR(ion_mobility[0][0], -1.0)  // Should return -1.0 for missing IM data
    TEST_REAL_SIMILAR(ion_mobility[0][1], -1.0)
  }
}
END_SECTION

START_SECTION((void get2DPeakData(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, std::vector<float>& rt, std::vector<float>& mz, std::vector<float>& intensity) const))
{
  MSExperiment exp;
  
  // Create test spectra
  MSSpectrum s1{
    {100.0, 1000.0},
    {200.0, 2000.0}
  };
  s1.setRT(1.0);
  s1.setMSLevel(1);
  
  MSSpectrum s2{
    {150.0, 1500.0},
    {250.0, 2500.0}
  };
  s2.setRT(2.0);
  s2.setMSLevel(1);
  
  exp.addSpectrum(s1);
  exp.addSpectrum(s2);
  
  // Test 1: Full range
  {
    std::vector<float> rt, mz, intensity;
    exp.get2DPeakData(0.0, 4.0, 0.0, 300.0, 1, rt, mz, intensity);
    
    TEST_EQUAL(rt.size(), 4)
    TEST_EQUAL(mz.size(), 4)
    TEST_EQUAL(intensity.size(), 4)
    
    // Check all peaks in order
    TEST_REAL_SIMILAR(rt[0], 1.0)
    TEST_REAL_SIMILAR(mz[0], 100.0)
    TEST_REAL_SIMILAR(intensity[0], 1000.0)
    
    TEST_REAL_SIMILAR(rt[1], 1.0)
    TEST_REAL_SIMILAR(mz[1], 200.0)
    TEST_REAL_SIMILAR(intensity[1], 2000.0)
    
    TEST_REAL_SIMILAR(rt[2], 2.0)
    TEST_REAL_SIMILAR(mz[2], 150.0)
    TEST_REAL_SIMILAR(intensity[2], 1500.0)
    
    TEST_REAL_SIMILAR(rt[3], 2.0)
    TEST_REAL_SIMILAR(mz[3], 250.0)
    TEST_REAL_SIMILAR(intensity[3], 2500.0)
  }
  
  // Test 2: Limited RT range
  {
    std::vector<float> rt, mz, intensity;
    exp.get2DPeakData(1.5, 2.5, 0.0, 300.0, 1, rt, mz, intensity);
    
    TEST_EQUAL(rt.size(), 2)
    TEST_EQUAL(mz.size(), 2)
    TEST_EQUAL(intensity.size(), 2)
    
    TEST_REAL_SIMILAR(rt[0], 2.0)
    TEST_REAL_SIMILAR(mz[0], 150.0)
    TEST_REAL_SIMILAR(intensity[0], 1500.0)
  }
  
  // Test 3: Limited MZ range
  {
    std::vector<float> rt, mz, intensity;
    exp.get2DPeakData(0.0, 4.0, 120.0, 180.0, 1, rt, mz, intensity);
    
    TEST_EQUAL(rt.size(), 1)
    TEST_EQUAL(mz.size(), 1)
    TEST_EQUAL(intensity.size(), 1)
    
    TEST_REAL_SIMILAR(rt[0], 2.0)
    TEST_REAL_SIMILAR(mz[0], 150.0)
    TEST_REAL_SIMILAR(intensity[0], 1500.0)
  }
  
  // Test 4: Empty range
  {
    std::vector<float> rt, mz, intensity;
    exp.get2DPeakData(5.0, 6.0, 0.0, 300.0, 1, rt, mz, intensity);
    
    TEST_EQUAL(rt.empty(), true)
    TEST_EQUAL(mz.empty(), true)
    TEST_EQUAL(intensity.empty(), true)
  }
}
END_SECTION

START_SECTION((std::vector<std::vector<MSExperiment::CoordinateType>> aggregateFromMatrix(const Matrix<double>& ranges, unsigned int ms_level, const std::string& mz_agg) const))
{
    // Create test experiment with known data
    MSExperiment exp;
    exp.resize(4);

    // First spectrum (MS1) at RT=1.0
    exp[0] = MSSpectrum{
        {100.0, 1000.0},
        {200.0, 2000.0},
        {300.0, 3000.0}
    };    
    exp[0].setRT(1.0);
    exp[0].setMSLevel(1);

    // Second spectrum (MS2) at RT=2.0
    exp[1] = MSSpectrum{
      {150.0, 1500.0},
      {250.0, 2500.0}
    };
    exp[1].setRT(2.0);
    exp[1].setMSLevel(2);

    // Third spectrum (MS1) at RT=3.0
    exp[2] = MSSpectrum{
        {100.0, 1100.0},
        {200.0, 2100.0},
        {300.0, 3100.0}
    };
    exp[2].setRT(3.0);
    exp[2].setMSLevel(1);

    // Fourth spectrum (MS1) at RT=4.0
    exp[3] = MSSpectrum{
        {100.0, 1200.0},
        {200.0, 2200.0},
        {300.0, 3200.0}
    };
    exp[3].setRT(4.0);
    exp[3].setMSLevel(1);

    exp.updateRanges();

    // Test 1: Sum aggregation for MS1 spectra
    {
        Matrix<double> ranges(2, 4); // two rt-mz ranges, four columns with min_mz, max_mz, min_rt, max_rt
        // Range 1: m/z 90-110, RT 0-3.5 (covers first and third spectra)
        ranges(0, 0) = 90.0;  ranges(0, 1) = 110.0;
        ranges(0, 2) = 0.0;   ranges(0, 3) = 3.5;

        // Range 2: m/z 190-210, RT 0-5.0 (covers first, third, and fourth spectra)
        ranges(1, 0) = 190.0; ranges(1, 1) = 210.0;
        ranges(1, 2) = 0.0;   ranges(1, 3) = 5.0;

        auto result = exp.aggregateFromMatrix(ranges, 1, "sum");

        // Check results
        TEST_EQUAL(result.size(), 2);

        // Check Range 1 results
        TEST_EQUAL(result[0].size(), 2);  // Should cover 2 spectra (1 and 3)
        TEST_EQUAL(result[0][0], 1000.0); // First spectrum intensity RT = 1.0
        TEST_EQUAL(result[0][1], 1100.0); // Third spectrum intensity RT = 3.0

        // Check Range 2 results
        TEST_EQUAL(result[1].size(), 3);   // Should cover 3 spectra (1, 3, and 4)
        TEST_EQUAL(result[1][0], 2000.0);  // First spectrum intensity
        TEST_EQUAL(result[1][1], 2100.0);  // Third spectrum intensity
        TEST_EQUAL(result[1][2], 2200.0);  // Fourth spectrum intensity
    }

    // Test 2: Max aggregation for MS1 spectra
    {
        Matrix<double> ranges(1, 4);
        // Range: m/z 100-300, RT 1-4 (covers first, third, and fourth spectra)
        ranges(0,0) = 100.0; ranges(0,1) = 300.0;
        ranges(0,2) = 1.0;   ranges(0,3) = 4.0;

        auto result = exp.aggregateFromMatrix(ranges, 1, "max");

        // Check results
        TEST_EQUAL(result.size(), 1);
        TEST_EQUAL(result[0].size(), 3);
        TEST_EQUAL(result[0][0], 3000.0); // First spectrum max intensity in range
        TEST_EQUAL(result[0][1], 3100.0); // Third spectrum max intensity in range
        TEST_EQUAL(result[0][2], 3200.0); // Fourth spectrum max intensity in range
    }

    // Test 3: Min aggregation for MS1 spectra
    {
        Matrix<double> ranges(1, 4);
        // Range: m/z 100-300, RT 1-4 (covers first, third, and fourth spectra)
        ranges(0,0) = 100.0; ranges(0,1) = 300.0;
        ranges(0,2) = 1.0;   ranges(0,3) = 4.0;
        auto result = exp.aggregateFromMatrix(ranges, 1, "min");

        // Check results
        TEST_EQUAL(result.size(), 1);
        TEST_EQUAL(result[0].size(), 3);
        TEST_EQUAL(result[0][0], 1000.0); // First spectrum min intensity in range
        TEST_EQUAL(result[0][1], 1100.0); // Third spectrum min intensity in range
        TEST_EQUAL(result[0][2], 1200.0); // Fourth spectrum min intensity in range
    }

    // Test 4: Mean aggregation for MS1 spectra
    {
        Matrix<double> ranges(1, 4);
        // Range: m/z 100-300, RT 0-5 (covers all MS1 spectra)
        ranges(0,0) = 100.0; ranges(0,1) = 300.0;
        ranges(0,2) = 0.0;   ranges(0,3) = 5.0;
        auto result = exp.aggregateFromMatrix(ranges, 1, "mean");

        // Check results
        TEST_EQUAL(result.size(), 1);
        TEST_EQUAL(result[0].size(), 3);
        TEST_REAL_SIMILAR(result[0][0], 2000.0); // First spectrum mean intensity
        TEST_REAL_SIMILAR(result[0][1], 2100.0); // Third spectrum mean intensity
        TEST_REAL_SIMILAR(result[0][2], 2200.0); // Fourth spectrum mean intensity
    }

    // Test 5: Invalid aggregation function
    {
        Matrix<double> ranges(1, 4);
        ranges(0,0) = 100.0; ranges(0,1) = 200.0;
        ranges(0,2) = 1.0;   ranges(0,3) = 2.0;
        TEST_EXCEPTION(Exception::InvalidValue, exp.aggregateFromMatrix(ranges, 1, "invalid_agg"));
    }

    // Test 6: Invalid matrix dimensions (not 4 columns)
    {
        Matrix<double> ranges(1, 3);
        ranges(0,0) = 100.0; ranges(0,1) = 200.0; ranges(0,2) = 1.0;
        TEST_EXCEPTION(Exception::InvalidParameter, exp.aggregateFromMatrix(ranges, 1, "sum"));
    }

    // Test 7: Empty ranges matrix
    {
        Matrix<double> ranges(0, 4); // 0 rows, 4 columns
        auto result = exp.aggregateFromMatrix(ranges, 1, "sum");
        TEST_EQUAL(result.size(), 0);         // Expect an empty result
    }

    // Test 8: No matching spectra for given ms_level
    {
        Matrix<double> ranges(1, 4);
        ranges(0,0) = 100.0; ranges(0,1) = 300.0;
        ranges(0,2) = 0.0;   ranges(0,3) = 5.0;
        auto result = exp.aggregateFromMatrix(ranges, 3, "sum"); // No spectra with MS level 3

        // Failed to extract anything -> empty vector
        TEST_EQUAL(result.size(), 0);
    }

}
END_SECTION

START_SECTION((std::vector<MSChromatogram> extractXICsFromMatrix(const Matrix<double>& ranges, unsigned int ms_level, const std::string& mz_agg) const))
{
    // Create test experiment with known data
    MSExperiment exp;
    exp.resize(4);

     // First spectrum (MS1) at RT=1.0
    exp[0] = MSSpectrum{
        {100.0, 1000.0},
        {200.0, 2000.0},
        {300.0, 3000.0}
    };    
    exp[0].setRT(1.0);
    exp[0].setMSLevel(1);

    // Second spectrum (MS2) at RT=2.0
    exp[1] = MSSpectrum{
      {150.0, 1500.0},
      {250.0, 2500.0}
    };
    exp[1].setRT(2.0);
    exp[1].setMSLevel(2);

    // Third spectrum (MS1) at RT=3.0
    exp[2] = MSSpectrum{
        {100.0, 1100.0},
        {200.0, 2100.0},
        {300.0, 3100.0}
    };
    exp[2].setRT(3.0);
    exp[2].setMSLevel(1);

    // Fourth spectrum (MS1) at RT=4.0
    exp[3] = MSSpectrum{
        {100.0, 1200.0},
        {200.0, 2200.0},
        {300.0, 3200.0}
    };
    exp[3].setRT(4.0);
    exp[3].setMSLevel(1);

    exp.updateRanges();

    // Test 1: Sum aggregation for MS1 spectra
    {
        Matrix<double> ranges(2, 4);
        // Range 1: m/z 90-110, RT 0-3.5 (covers first and third spectra)
        ranges(0,0) = 90.0;  ranges(0,1) = 110.0;
        ranges(0,2) = 0.0;   ranges(0,3) = 3.5;
        // Range 2: m/z 190-210, RT 0-5.0 (covers first, third, and fourth spectra)
        ranges(1,0) = 190.0; ranges(1,1) = 210.0;
        ranges(1,2) = 0.0;   ranges(1,3) = 5.0;

        unsigned int ms_level = 1;
        std::string mz_agg = "sum";

        auto result = exp.extractXICsFromMatrix(ranges, ms_level, mz_agg);

        // Check results
        TEST_EQUAL(result.size(), 2);

        // Check Range 1 results
        TEST_EQUAL(result[0].size(), 2);  // Should cover 2 spectra
        TEST_REAL_SIMILAR(result[0][0].getIntensity(), 1000.0); // First spectrum intensity
        TEST_REAL_SIMILAR(result[0][1].getIntensity(), 1100.0); // Third spectrum intensity

        // Check Range 2 results
        TEST_EQUAL(result[1].size(), 3);   // Should cover 3 spectra
        TEST_REAL_SIMILAR(result[1][0].getIntensity(), 2000.0);  // First spectrum intensity
        TEST_REAL_SIMILAR(result[1][1].getIntensity(), 2100.0);  // Third spectrum intensity
        TEST_REAL_SIMILAR(result[1][2].getIntensity(), 2200.0);  // Fourth spectrum intensity
    }

}
END_SECTION

START_SECTION((void get2DPeakDataIM(CoordinateType min_rt, CoordinateType max_rt, CoordinateType min_mz, CoordinateType max_mz, Size ms_level, std::vector<float>& rt, std::vector<float>& mz, std::vector<float>& intensity, std::vector<float>& ion_mobility) const))
{
  MSExperiment exp;
  
  // Create test spectra with ion mobility data
  MSSpectrum s1{
    {100.0, 1000.0},
    {200.0, 2000.0}
  };
  s1.setRT(1.0);
  s1.setMSLevel(1);
  DataArrays::FloatDataArray im1;
  im1.setName("Ion Mobility");
  im1.push_back(0.8f);
  im1.push_back(1.2f);
  im1.setMetaValue("unit", "millisecond");
  s1.getFloatDataArrays().push_back(im1);
  
  MSSpectrum s2{
    {150.0, 1500.0},
    {250.0, 2500.0}
  };
  s2.setRT(2.0);
  s2.setMSLevel(1);
  DataArrays::FloatDataArray im2;
  im2.setName("Ion Mobility");
  im2.push_back(1.5f);
  im2.push_back(2.0f);
  im2.setMetaValue("unit", "millisecond");
  s2.getFloatDataArrays().push_back(im2);
  
  exp.addSpectrum(s1);
  exp.addSpectrum(s2);
  
  // Test 1: Full range, MS level 1
  {
    std::vector<float> rt, mz, intensity, ion_mobility;
    exp.get2DPeakDataIM(0.0, 4.0, 0.0, 300.0, 1, rt, mz, intensity, ion_mobility);
    
    TEST_EQUAL(rt.size(), 4)
    TEST_EQUAL(mz.size(), 4)
    TEST_EQUAL(intensity.size(), 4)
    TEST_EQUAL(ion_mobility.size(), 4)
    
    // Check all peaks in order
    TEST_REAL_SIMILAR(rt[0], 1.0)
    TEST_REAL_SIMILAR(mz[0], 100.0)
    TEST_REAL_SIMILAR(intensity[0], 1000.0)
    TEST_REAL_SIMILAR(ion_mobility[0], 0.8)
    
    TEST_REAL_SIMILAR(rt[1], 1.0)
    TEST_REAL_SIMILAR(mz[1], 200.0)
    TEST_REAL_SIMILAR(intensity[1], 2000.0)
    TEST_REAL_SIMILAR(ion_mobility[1], 1.2)
    
    TEST_REAL_SIMILAR(rt[2], 2.0)
    TEST_REAL_SIMILAR(mz[2], 150.0)
    TEST_REAL_SIMILAR(intensity[2], 1500.0)
    TEST_REAL_SIMILAR(ion_mobility[2], 1.5)
    
    TEST_REAL_SIMILAR(rt[3], 2.0)
    TEST_REAL_SIMILAR(mz[3], 250.0)
    TEST_REAL_SIMILAR(intensity[3], 2500.0)
    TEST_REAL_SIMILAR(ion_mobility[3], 2.0)
  }
  
  // Test 2: Limited RT range
  {
    std::vector<float> rt, mz, intensity, ion_mobility;
    exp.get2DPeakDataIM(1.5, 2.5, 0.0, 300.0, 1, rt, mz, intensity, ion_mobility);
    
    TEST_EQUAL(rt.size(), 2)
    TEST_EQUAL(mz.size(), 2)
    TEST_EQUAL(intensity.size(), 2)
    TEST_EQUAL(ion_mobility.size(), 2)
    
    TEST_REAL_SIMILAR(rt[0], 2.0)
    TEST_REAL_SIMILAR(mz[0], 150.0)
    TEST_REAL_SIMILAR(intensity[0], 1500.0)
    TEST_REAL_SIMILAR(ion_mobility[0], 1.5)
  }
  
  // Test 3: Spectrum without ion mobility data
  {
    MSExperiment exp_no_im;
    MSSpectrum s_no_im{
      {100.0, 1000.0},
      {200.0, 2000.0}
    };
    s_no_im.setRT(1.0);
    s_no_im.setMSLevel(1);
    exp_no_im.addSpectrum(s_no_im);
    
    std::vector<float> rt, mz, intensity, ion_mobility;
    exp_no_im.get2DPeakDataIM(0.0, 4.0, 0.0, 300.0, 1, rt, mz, intensity, ion_mobility);
    
    TEST_EQUAL(rt.size(), 2)
    TEST_EQUAL(mz.size(), 2)
    TEST_EQUAL(intensity.size(), 2)
    TEST_EQUAL(ion_mobility.size(), 2)
    TEST_REAL_SIMILAR(ion_mobility[0], -1.0)  // Should return -1.0 for missing IM data
    TEST_REAL_SIMILAR(ion_mobility[1], -1.0)
  }
}
END_SECTION

START_SECTION((template<class MzReductionFunctionType> std::vector<MSChromatogram> extractXICs(const std::vector<std::pair<RangeMZ, RangeRT>>& mz_rt_ranges, unsigned int ms_level, MzReductionFunctionType func_mz_reduction = SumIntensityReduction()) const))
{
    // Create test experiment with known data
    PeakMap exp;
    exp.resize(4);

    // First spectrum (MS1) at RT=1.0
    exp[0] = MSSpectrum{
        {100.0, 1000.0},
        {200.0, 2000.0},
        {300.0, 3000.0}
    };
    exp[0].setRT(1.0);
    exp[0].setMSLevel(1);

    // Second spectrum (MS2) at RT=2.0
    exp[1] = MSSpectrum{
        {150.0, 1500.0},
        {250.0, 2500.0}
    };
    exp[1].setRT(2.0);
    exp[1].setMSLevel(2);

    // Third spectrum (MS1) at RT=3.0
    exp[2] = MSSpectrum{
        {100.0, 1100.0},
        {200.0, 2100.0},
        {300.0, 3100.0}
    };
    exp[2].setRT(3.0);
    exp[2].setMSLevel(1);

    // Fourth spectrum (MS1) at RT=4.0
    exp[3] = MSSpectrum{
        {100.0, 1200.0},
        {200.0, 2200.0},
        {300.0, 3200.0}
    };
    exp[3].setRT(4.0);
    exp[3].setMSLevel(1);

    // Update the ranges of the experiment
    exp.updateRanges();

    // Test 1: Normal case - MS1 spectra using default reduction function
    {
        std::vector<std::pair<RangeMZ, RangeRT>> ranges;
        // Range 1: covers first peak of all MS1 spectra
        ranges.push_back(std::make_pair(
            RangeMZ(90.0, 110.0),
            RangeRT(0.0, 5.0)
        ));
        // Range 2: covers second peak of all MS1 spectra
        ranges.push_back(std::make_pair(
            RangeMZ(190.0, 210.0),
            RangeRT(0.0, 5.0)
        ));

        // Use default reduction function (SumIntensityReduction)
        auto chromatograms = exp.extractXICs(ranges, 1);

        // Check results
        TEST_EQUAL(chromatograms.size(), 2);

        // Check Range 1 chromatogram
        TEST_EQUAL(chromatograms[0].size(), 3); // Should cover 3 spectra
        TEST_REAL_SIMILAR(chromatograms[0][0].getRT(), 1.0);
        TEST_REAL_SIMILAR(chromatograms[0][0].getIntensity(), 1000.0);

        TEST_REAL_SIMILAR(chromatograms[0][1].getRT(), 3.0);
        TEST_REAL_SIMILAR(chromatograms[0][1].getIntensity(), 1100.0);

        TEST_REAL_SIMILAR(chromatograms[0][2].getRT(), 4.0);
        TEST_REAL_SIMILAR(chromatograms[0][2].getIntensity(), 1200.0);

        // Check m/z value of the chromatogram
        TEST_REAL_SIMILAR(chromatograms[0].getProduct().getMZ(), (90.0 + 110.0) / 2.0);

        // Check Range 2 chromatogram
        TEST_EQUAL(chromatograms[1].size(), 3); // Should cover 3 spectra
        TEST_REAL_SIMILAR(chromatograms[1][0].getRT(), 1.0);
        TEST_REAL_SIMILAR(chromatograms[1][0].getIntensity(), 2000.0);

        TEST_REAL_SIMILAR(chromatograms[1][1].getRT(), 3.0);
        TEST_REAL_SIMILAR(chromatograms[1][1].getIntensity(), 2100.0);

        TEST_REAL_SIMILAR(chromatograms[1][2].getRT(), 4.0);
        TEST_REAL_SIMILAR(chromatograms[1][2].getIntensity(), 2200.0);

        // Check m/z value of the chromatogram
        TEST_REAL_SIMILAR(chromatograms[1].getProduct().getMZ(), (190.0 + 210.0) / 2.0);
    }

    // Test 2: MS2 spectra
    {
        std::vector<std::pair<RangeMZ, RangeRT>> ranges;
        ranges.push_back(std::make_pair(
            RangeMZ(140.0, 160.0),
            RangeRT(1.5, 2.5)
        ));

        auto chromatograms = exp.extractXICs(ranges, 2);

        TEST_EQUAL(chromatograms.size(), 1);
        TEST_EQUAL(chromatograms[0].size(), 1);
        TEST_REAL_SIMILAR(chromatograms[0][0].getRT(), 2.0);
        TEST_REAL_SIMILAR(chromatograms[0][0].getIntensity(), 1500.0);

        // Check m/z value of the chromatogram
        TEST_REAL_SIMILAR(chromatograms[0].getProduct().getMZ(), (140.0 + 160.0) / 2.0);
    }

    // Test 3: Custom reduction function (average intensity)
    {
        std::vector<std::pair<RangeMZ, RangeRT>> ranges;
        ranges.push_back(std::make_pair(
            RangeMZ(90.0, 310.0),  // Covers all peaks
            RangeRT(0.0, 5.0)      // Covers all spectra
        ));

        // Mean intensity in m/z range
        auto chromatograms = exp.extractXICs(ranges, 1,
            [](MSSpectrum::ConstIterator begin_it, MSSpectrum::ConstIterator end_it) -> double
            {
                if (begin_it == end_it) return 0.0;

                double acc = std::accumulate(begin_it, end_it, 0.0,
                    [](double a, const Peak1D& b) { return a + b.getIntensity(); });
                return acc / static_cast<double>(std::distance(begin_it, end_it));
            });

        TEST_EQUAL(chromatograms.size(), 1);
        TEST_EQUAL(chromatograms[0].size(), 3);

        TEST_REAL_SIMILAR(chromatograms[0][0].getRT(), 1.0);
        TEST_REAL_SIMILAR(chromatograms[0][0].getIntensity(), 2000.0); // Average of [1000, 2000, 3000]

        TEST_REAL_SIMILAR(chromatograms[0][1].getRT(), 3.0);
        TEST_REAL_SIMILAR(chromatograms[0][1].getIntensity(), 2100.0); // Average of [1100, 2100, 3100]

        TEST_REAL_SIMILAR(chromatograms[0][2].getRT(), 4.0);
        TEST_REAL_SIMILAR(chromatograms[0][2].getIntensity(), 2200.0); // Average of [1200, 2200, 3200]

        // Check m/z value of the chromatogram
        TEST_REAL_SIMILAR(chromatograms[0].getProduct().getMZ(), (90.0 + 310.0) / 2.0);
    }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

