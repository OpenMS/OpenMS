// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/ConsensusFeature.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Element.h>
#include <OpenMS/KERNEL/Feature.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/METADATA/PeptideIdentification.h>

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusFeature, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusFeature* ptr = nullptr;
ConsensusFeature* nullPointer = nullptr;
START_SECTION((ConsensusFeature()))
  ptr = new ConsensusFeature();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~ConsensusFeature()))
  delete ptr;
END_SECTION

Feature tmp_feature;
tmp_feature.setRT(1);
tmp_feature.setMZ(2);
tmp_feature.setIntensity(200.0f);
tmp_feature.setUniqueId(3);

Feature tmp_feature2;
tmp_feature2.setRT(2);
tmp_feature2.setMZ(3);
tmp_feature2.setIntensity(300.0f);
tmp_feature2.setUniqueId(5);

Feature tmp_feature3;
tmp_feature3.setRT(3);
tmp_feature3.setMZ(4);
tmp_feature3.setIntensity(400.0f);
tmp_feature3.setUniqueId(7);


START_SECTION(([ConsensusFeature::SizeLess] bool operator () ( ConsensusFeature const & left, ConsensusFeature const & right ) const))
  ConsensusFeature c1(tmp_feature);
  c1.insert(1, tmp_feature);
  c1.insert(2, tmp_feature3);

  ConsensusFeature c2(tmp_feature2);
  c2.insert(1,tmp_feature2);

  ConsensusFeature::SizeLess sl;

  TEST_EQUAL(sl(c1,c2), 0);
  TEST_EQUAL(sl(c2,c1), 1);
END_SECTION

START_SECTION(([ConsensusFeature::SizeLess] bool operator () ( ConsensusFeature const & left, UInt64 const & right ) const))
  ConsensusFeature c1(tmp_feature);
  c1.insert(1, tmp_feature);
  c1.insert(2, tmp_feature3);

  ConsensusFeature c2(tmp_feature);
  c2.insert(1,tmp_feature);
  c2.insert(2,tmp_feature2);
  c2.insert(3,tmp_feature3);

  UInt64 rhs_size = c2.size();

  ConsensusFeature::SizeLess sl;

  TEST_EQUAL(sl(c1,rhs_size), 1);
  TEST_EQUAL(sl(c2,rhs_size), 0);
END_SECTION

START_SECTION(([ConsensusFeature::SizeLess] bool operator () ( UInt64 const & left, ConsensusFeature const & right ) const))
  ConsensusFeature c1(tmp_feature);
  c1.insert(1, tmp_feature);
  c1.insert(2, tmp_feature3);

  ConsensusFeature c2(tmp_feature);
  c2.insert(1,tmp_feature);
  c2.insert(2,tmp_feature2);
  c2.insert(3,tmp_feature3);

  UInt64 lhs_size = c1.size();

  ConsensusFeature::SizeLess sl;

  TEST_EQUAL(sl(lhs_size, c1), 0);
  TEST_EQUAL(sl(lhs_size, c2), 1);
END_SECTION

START_SECTION(([ConsensusFeature::SizeLess] bool operator () ( const UInt64 & left, const UInt64 & right ) const))
  ConsensusFeature c1(tmp_feature);
  c1.insert(1, tmp_feature);
  c1.insert(2, tmp_feature3);

  ConsensusFeature c2(tmp_feature);
  c2.insert(1,tmp_feature);
  c2.insert(2,tmp_feature2);
  c2.insert(3,tmp_feature3);

  UInt64 lhs_size = c1.size(), rhs_size = c2.size();

  ConsensusFeature::SizeLess sl;

  TEST_EQUAL(sl(lhs_size, rhs_size), 1);
  TEST_EQUAL(sl(rhs_size, lhs_size), 0);
END_SECTION

START_SECTION(([ConsensusFeature::MapsLess] bool operator () ( ConsensusFeature const & left, ConsensusFeature const & right ) const))
  ConsensusFeature c1(tmp_feature);
  c1.insert(1, tmp_feature);
  c1.insert(2, tmp_feature3);

  ConsensusFeature c2(tmp_feature);
  c2.insert(3,tmp_feature);
  c2.insert(4,tmp_feature2);
  c2.insert(5,tmp_feature3);

  ConsensusFeature::MapsLess ml;

  TEST_EQUAL(ml(c1,c1), 0);
  TEST_EQUAL(ml(c1,c2), 1);
  TEST_EQUAL(ml(c2,c1), 0);
  TEST_EQUAL(ml(c2,c2), 0);
END_SECTION

START_SECTION((ConsensusFeature& operator=(const ConsensusFeature &rhs)))
  ConsensusFeature cons(tmp_feature);
  cons.insert(1,tmp_feature);

  ConsensusFeature cons_copy;
  cons_copy = cons;

  TEST_REAL_SIMILAR(cons_copy.getRT(),1.0)
  TEST_REAL_SIMILAR(cons_copy.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons_copy.getIntensity(),200.0)
  TEST_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_EQUAL((cons_copy.begin())->getUniqueId(),3)
  TEST_EQUAL((cons_copy.begin())->getIntensity(),200)
END_SECTION

START_SECTION((ConsensusFeature(const ConsensusFeature &rhs)))
{
  ConsensusFeature cons(tmp_feature);
  cons.insert(1,tmp_feature);
  ConsensusFeature cons_copy(cons);

  TEST_REAL_SIMILAR(cons_copy.getRT(),1.0)
  TEST_REAL_SIMILAR(cons_copy.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons_copy.getIntensity(),200.0)
  TEST_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_EQUAL((cons_copy.begin())->getUniqueId(),3)
  TEST_EQUAL((cons_copy.begin())->getIntensity(),200)
}
END_SECTION

START_SECTION((ConsensusFeature(ConsensusFeature &&rhs)))
{
#ifndef OPENMS_COMPILER_MSVC
  // Ensure that ConsensusFeature has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  // Note that MSVS does not support noexcept move constructors for STL
  // constructs such as std::map.
  TEST_EQUAL(noexcept(ConsensusFeature(std::declval<ConsensusFeature&&>())), true)
#endif

  ConsensusFeature cons(tmp_feature);
  cons.insert(1,tmp_feature);

  ConsensusFeature cons_copy(std::move(cons));

  TEST_REAL_SIMILAR(cons_copy.getRT(),1.0)
  TEST_REAL_SIMILAR(cons_copy.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons_copy.getIntensity(),200.0)
  TEST_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_EQUAL((cons_copy.begin())->getUniqueId(),3)
  TEST_EQUAL((cons_copy.begin())->getIntensity(),200)
}
END_SECTION

START_SECTION((void insert(const HandleSetType &handle_set)))
  ConsensusFeature::HandleSetType hs;
  FeatureHandle fh;
  for ( UInt i = 0; i < 3; ++i )
  {
    fh.setRT(i*77.7);
    fh.setMapIndex(i+10);
    fh.setUniqueId(i+1000);
    hs.insert(fh);
  }
  ConsensusFeature cf;
  cf.insert(hs);

  TEST_EQUAL(cf.size(),3);
  TEST_EQUAL(cf.begin()->getMapIndex(),10);
  TEST_EQUAL(cf.rbegin()->getMapIndex(),12);

END_SECTION

START_SECTION((void insert(UInt64 map_index, const Peak2D &element, UInt64 element_index)))
  ConsensusFeature cf;
  Peak2D el;
  for ( UInt i = 0; i < 3; ++i )
  {
    el.setRT(i*77.7);
    cf.insert(10-i,el,i+1000);
    TEST_EQUAL(cf.size(),i+1);
    TEST_REAL_SIMILAR(cf.begin()->getRT(),i*77.7);
    TEST_EQUAL(cf.begin()->getMapIndex(),10-i);
    TEST_EQUAL(cf.begin()->getUniqueId(),i+1000);
  }

END_SECTION

START_SECTION((void insert(UInt64 map_index, const BaseFeature &element)))
  ConsensusFeature cf;
  BaseFeature el;
  for ( UInt i = 0; i < 3; ++i )
  {
    el.setRT(i*77.7);
    el.setCharge(2*i);
    el.setUniqueId(i+1000);
    cf.insert(10-i,el);
    TEST_EQUAL(cf.size(),i+1);
    TEST_REAL_SIMILAR(cf.begin()->getRT(),i*77.7);
    TEST_EQUAL(cf.begin()->getCharge(),2*i);
    TEST_EQUAL(cf.begin()->getMapIndex(),10-i);
    TEST_EQUAL(cf.begin()->getUniqueId(),i+1000);
  }
END_SECTION

START_SECTION((ConsensusFeature(const BaseFeature &feature)))
  BaseFeature f;
  f.setCharge(-17);
  f.setRT(44324.6);
  f.setMZ(867.4);
  f.setUniqueId(23);
  f.getPeptideIdentifications().resize(1);
  const BaseFeature& f_cref = f;
  ConsensusFeature cf(f_cref);

  TEST_EQUAL(cf.getRT(),44324.6);
  TEST_EQUAL(cf.getMZ(),867.4);
  TEST_EQUAL(cf.getCharge(),-17);
  TEST_EQUAL(cf.getPeptideIdentifications().size(), 1);
  TEST_EQUAL(cf.empty(), true);
END_SECTION

START_SECTION((ConsensusFeature(UInt64 map_index, const BaseFeature& element)))
{
  BaseFeature f;
  f.setCharge(-17);
  f.setRT(44324.6);
  f.setMZ(867.4);
  f.setIntensity(1000);
  f.setUniqueId(23);
  f.getPeptideIdentifications().resize(1);
  ConsensusFeature cf(99, f);

  TEST_EQUAL(cf.getRT(),44324.6);
  TEST_EQUAL(cf.getMZ(),867.4);
  TEST_EQUAL(cf.getCharge(),-17);
  TEST_EQUAL(cf.getPeptideIdentifications().size(), 1);
  ConsensusFeature::HandleSetType::const_iterator it = cf.begin();
  TEST_EQUAL(it->getMapIndex(), 99)
  TEST_EQUAL(it->getUniqueId(), 23)
  TEST_EQUAL(it->getIntensity(), 1000)
}
END_SECTION

START_SECTION([EXTRA](ConsensusFeature(UInt64 map_index, const Feature &element)))
   ConsensusFeature cons(1,tmp_feature);
  cons.setUniqueId(3);

  TEST_REAL_SIMILAR(cons.getRT(),1.0)
  TEST_REAL_SIMILAR(cons.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons.getIntensity(),200.0)
  ConsensusFeature::HandleSetType::const_iterator it = cons.begin();
  TEST_EQUAL(it->getMapIndex(),1)
  TEST_EQUAL(it->getUniqueId(),3)
  TEST_EQUAL(it->getIntensity(),200)
END_SECTION

START_SECTION((ConsensusFeature(UInt64 map_index, const Peak2D &element, UInt64 element_index)))
  Peak2D f;
  f.setIntensity(-17);
  const Peak2D& f_cref = f;
  ConsensusFeature cf(99,f_cref,23);

  ConsensusFeature::HandleSetType::const_iterator it = cf.begin();
  TEST_EQUAL(it->getMapIndex(),99);
  TEST_EQUAL(it->getUniqueId(),23);
  TEST_EQUAL(it->getIntensity(),-17);
END_SECTION

START_SECTION([EXTRA](ConsensusFeature(UInt64 map_index, const ConsensusFeature &element)))
  ConsensusFeature f;
  f.setUniqueId(23);
  f.setIntensity(-17);
  const ConsensusFeature& f_cref = f;
  ConsensusFeature cf(99,f_cref);

  ConsensusFeature::HandleSetType::const_iterator it = cf.begin();
  TEST_EQUAL(it->getMapIndex(),99);
  TEST_EQUAL(it->getUniqueId(),23);
  TEST_EQUAL(it->getIntensity(),-17);
END_SECTION

START_SECTION((DRange<1> getIntensityRange() const))
  ConsensusFeature cons;
  Feature f;
  f.setIntensity(0.0f);
  f.setUniqueId(0);
  cons.insert(0,f);
  f.setUniqueId(1);
  f.setIntensity(200.0f);
  cons.insert(0,f);

  TEST_REAL_SIMILAR(cons.getIntensityRange().minX(),0.0)
  TEST_REAL_SIMILAR(cons.getIntensityRange().maxX(),200.0)
END_SECTION

START_SECTION((DRange<2> getPositionRange() const))
  ConsensusFeature cons;
  Feature f;
  f.setRT(1.0);
  f.setMZ(500.0);
  f.setUniqueId(0);
  cons.insert(0,f);
  f.setRT(1000.0);
  f.setMZ(1500.0);
  f.setUniqueId(1);
  cons.insert(0,f);

  TEST_REAL_SIMILAR(cons.getPositionRange().minX(),1.0)
  TEST_REAL_SIMILAR(cons.getPositionRange().maxX(),1000.0)
  TEST_REAL_SIMILAR(cons.getPositionRange().minY(),500.0)
  TEST_REAL_SIMILAR(cons.getPositionRange().maxY(),1500.0)
END_SECTION

START_SECTION((const HandleSetType& getFeatures() const))
  ConsensusFeature cons;
  cons.insert(2,tmp_feature);
  const ConsensusFeature cons_copy(cons);

  ConsensusFeature::HandleSetType group = cons_copy.getFeatures();

  ConsensusFeature::HandleSetType::const_iterator it = group.begin();
  TEST_EQUAL(it->getMapIndex(),2)
  TEST_EQUAL(it->getUniqueId(),3)
  TEST_EQUAL(it->getIntensity(),200)
END_SECTION

START_SECTION(( std::vector<FeatureHandle> getFeatureList() const))
  ConsensusFeature cons;
  cons.insert(2,tmp_feature);
  const ConsensusFeature cons_copy(cons);

  std::vector<FeatureHandle> group = cons_copy.getFeatureList();

  TEST_EQUAL(group.size(),1)

  TEST_EQUAL(group[0].getMapIndex(),2)
  TEST_EQUAL(group[0].getUniqueId(),3)
  TEST_EQUAL(group[0].getIntensity(),200)
END_SECTION


START_SECTION((void insert(const ConsensusFeature& cf)))
  ConsensusFeature cons, cons_t;
  FeatureHandle h1(2,tmp_feature);
  h1.setUniqueId(3);
  FeatureHandle h2(4,tmp_feature);
  h2.setUniqueId(5);
  cons.insert(h1);
  cons.insert(h2);

  cons_t.insert(cons);

  ConsensusFeature::HandleSetType::const_iterator it = cons_t.begin();
  TEST_EQUAL(it->getMapIndex(),2)
  TEST_EQUAL(it->getUniqueId(),3)
  TEST_EQUAL(it->getIntensity(),200)
  ++it;
  TEST_EQUAL(it->getMapIndex(),4)
  TEST_EQUAL(it->getUniqueId(),5)
  TEST_EQUAL(it->getIntensity(),200)
  ++it;
  TEST_EQUAL(it==cons_t.end(), true)
END_SECTION

START_SECTION((void insert(const FeatureHandle &handle)))
  ConsensusFeature cons;
  FeatureHandle h1(2,tmp_feature);
  h1.setUniqueId(3);
  FeatureHandle h2(4,tmp_feature);
  h2.setUniqueId(5);
  cons.insert(h1);
  cons.insert(h2);

  ConsensusFeature::HandleSetType::const_iterator it = cons.begin();
  TEST_EQUAL(it->getMapIndex(),2)
  TEST_EQUAL(it->getUniqueId(),3)
  TEST_EQUAL(it->getIntensity(),200)
  ++it;
  TEST_EQUAL(it->getMapIndex(),4)
  TEST_EQUAL(it->getUniqueId(),5)
  TEST_EQUAL(it->getIntensity(),200)
  ++it;
  TEST_EQUAL(it==cons.end(), true)
END_SECTION

START_SECTION((void insert(UInt64 map_index, const BaseFeature &element)))
  ConsensusFeature cons;
  cons.insert(2, tmp_feature);

  ConsensusFeature::HandleSetType::const_iterator it = cons.begin();
  TEST_EQUAL(it->getMapIndex(),2)
  TEST_EQUAL(it->getUniqueId(),3)
  TEST_EQUAL(it->getIntensity(),200)
  ++it;
  TEST_EQUAL(it==cons.end(),true)
END_SECTION


START_SECTION((void computeConsensus()))
  ConsensusFeature cons;
  //one point
  cons.insert(2,tmp_feature);
  cons.computeConsensus();
  TEST_REAL_SIMILAR(cons.getIntensity(),200.0)
  TEST_REAL_SIMILAR(cons.getRT(),1.0)
  TEST_REAL_SIMILAR(cons.getMZ(),2.0)
  //two points
  cons.insert(4,tmp_feature2);
  cons.computeConsensus();
  TEST_REAL_SIMILAR(cons.getIntensity(),250.0)
  TEST_REAL_SIMILAR(cons.getRT(),1.5)
  TEST_REAL_SIMILAR(cons.getMZ(),2.5)
  //three points
  cons.insert(6,tmp_feature3);
  cons.computeConsensus();
  TEST_REAL_SIMILAR(cons.getIntensity(),300.0)
  TEST_REAL_SIMILAR(cons.getRT(),2.0)
  TEST_REAL_SIMILAR(cons.getMZ(),3.0)
END_SECTION

START_SECTION((void computeMonoisotopicConsensus()))
  ConsensusFeature cons;
  //one point
  cons.insert(2,tmp_feature);
  cons.computeMonoisotopicConsensus();
  TEST_REAL_SIMILAR(cons.getIntensity(),200.0)
  TEST_REAL_SIMILAR(cons.getRT(),1.0)
  TEST_REAL_SIMILAR(cons.getMZ(),2.0)
  //two points
  cons.insert(4,tmp_feature2);
  cons.computeMonoisotopicConsensus();
  TEST_REAL_SIMILAR(cons.getIntensity(),250.0)
  TEST_REAL_SIMILAR(cons.getRT(),1.5)
  TEST_REAL_SIMILAR(cons.getMZ(),2.0)
  //three points
  cons.insert(6,tmp_feature3);
  cons.computeMonoisotopicConsensus();
  TEST_REAL_SIMILAR(cons.getIntensity(),300.0)
  TEST_REAL_SIMILAR(cons.getRT(),2.0)
  TEST_REAL_SIMILAR(cons.getMZ(),2.0)
END_SECTION

START_SECTION((void computeDechargeConsensus(const FeatureMap& fm, bool intensity_weighted_averaging=false)))

  double proton_mass = ElementDB::getInstance()->getElement("H")->getMonoWeight();
  double natrium_mass = ElementDB::getInstance()->getElement("Na")->getMonoWeight();

  double m = 1000;
  double m1_add = 0.5;
  double mz1 = (m+m1_add+3*proton_mass) / 3;
  double m2_add = 1;
  double mz2 = (m+m2_add+1*proton_mass + 2*natrium_mass) / 3;
  double m3_add = -0.5;
  double mz3 = (m+m3_add+4*proton_mass + natrium_mass) / 5;

  FeatureMap fm;

  //one point
  ConsensusFeature cons;
  Feature tmp_feature;
  tmp_feature.setRT(100);
  tmp_feature.setMZ(mz1);
  tmp_feature.setIntensity(200.0f);
  tmp_feature.setCharge(3);
  tmp_feature.ensureUniqueId();
  fm.push_back(tmp_feature);
  cons.insert(2,tmp_feature);
  cons.computeDechargeConsensus(fm);
  TEST_REAL_SIMILAR(cons.getIntensity(),200.0)
  TEST_REAL_SIMILAR(cons.getRT(),100)
  TEST_REAL_SIMILAR(cons.getMZ(), m+m1_add);

  //two points
  Feature tmp_feature2;
  tmp_feature2.setRT(102);
  tmp_feature2.setMZ(mz2);
  tmp_feature2.setIntensity(400.0f);
  tmp_feature2.setCharge(3);
  tmp_feature2.ensureUniqueId();
  tmp_feature2.setMetaValue("dc_charge_adduct_mass", 2*natrium_mass + proton_mass);
  fm.push_back(tmp_feature2);
  cons.insert(4,tmp_feature2);
  cons.computeDechargeConsensus(fm, true);
  TEST_REAL_SIMILAR(cons.getIntensity(),600.0)
  TEST_REAL_SIMILAR(cons.getRT(),(100.0/3 + 102.0*2/3))
  TEST_REAL_SIMILAR(cons.getMZ(),((m+m1_add)/3 + (m+m2_add)*2/3))

  cons.computeDechargeConsensus(fm, false);
  TEST_REAL_SIMILAR(cons.getIntensity(),600.0)
  TEST_REAL_SIMILAR(cons.getRT(),(100.0/2 + 102.0/2))
  TEST_REAL_SIMILAR(cons.getMZ(),((m+m1_add)/2 + (m+m2_add)/2))

  //three points
  Feature tmp_feature3;
  tmp_feature3.setRT(101);
  tmp_feature3.setMZ(mz3);
  tmp_feature3.setIntensity(600.0f);
  tmp_feature3.setCharge(5);
  tmp_feature3.ensureUniqueId();
  tmp_feature3.setMetaValue("dc_charge_adduct_mass", 1*natrium_mass + 4*proton_mass);
  fm.push_back(tmp_feature3);
  cons.insert(4,tmp_feature3);
  cons.computeDechargeConsensus(fm, true);
  TEST_REAL_SIMILAR(cons.getIntensity(),1200.0)
  TEST_REAL_SIMILAR(cons.getRT(),(100.0/6 + 102.0/3 + 101.0/2))
  TEST_REAL_SIMILAR(cons.getMZ(),((m+m1_add)/6 + (m+m2_add)/3 + (m+m3_add)/2))

  cons.computeDechargeConsensus(fm, false);
  TEST_REAL_SIMILAR(cons.getIntensity(),1200.0)
  TEST_REAL_SIMILAR(cons.getRT(),(100.0/3 + 102.0/3 + 101.0/3))
  TEST_REAL_SIMILAR(cons.getMZ(),((m+m1_add)/3 + (m+m2_add)/3 + (m+m3_add)/3))

END_SECTION




START_SECTION((Size size() const))
{
  ConsensusFeature c1(tmp_feature);
  c1.insert(1, tmp_feature);
  c1.insert(2, tmp_feature3);
  TEST_EQUAL(c1.size(), 2)

  ConsensusFeature c2;
  TEST_EQUAL(c2.size(), 0)
  c2.insert(1,tmp_feature2);
  TEST_EQUAL(c2.size(), 1)
}
END_SECTION

START_SECTION((const_iterator begin() const))
{

  ConsensusFeature c;
  const ConsensusFeature& c2 = c;
  TEST_EQUAL(c2.begin()==c2.end(), true)
  c.insert(1,tmp_feature2);
  const ConsensusFeature& c3 = c;
  TEST_EQUAL(c3.begin()->getUniqueId(), 5)
}
END_SECTION

START_SECTION((iterator begin()))
{
  ConsensusFeature c;
  TEST_EQUAL(c.begin()==c.end(), true)
  c.insert(1,tmp_feature2);
  TEST_EQUAL(c.begin()->getUniqueId(), 5)
}
END_SECTION

START_SECTION((const_iterator end() const))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION((iterator end()))
  NOT_TESTABLE // tested above
END_SECTION


START_SECTION((const_reverse_iterator rbegin() const))
{

  ConsensusFeature c;
  const ConsensusFeature& c2 = c;
  TEST_EQUAL(c2.rbegin()==c2.rend(), true)
  c.insert(1,tmp_feature2);
  const ConsensusFeature& c3 = c;
  TEST_EQUAL(c3.rbegin()->getUniqueId(), 5)
}
END_SECTION

START_SECTION((reverse_iterator rbegin()))
{
  ConsensusFeature c;
  TEST_EQUAL(c.rbegin()==c.rend(), true)
  c.insert(1,tmp_feature2);
  TEST_EQUAL(c.rbegin()->getUniqueId(), 5)
}
END_SECTION


START_SECTION((const_reverse_iterator rend() const))
  NOT_TESTABLE // tested above
END_SECTION

START_SECTION((reverse_iterator rend()))
  NOT_TESTABLE // tested above
END_SECTION
START_SECTION((void clear()))
{
  ConsensusFeature c1(tmp_feature);
  c1.insert(1, tmp_feature);
  c1.insert(2, tmp_feature3);
  c1.clear();
  TEST_EQUAL(c1.size(), 0)

  ConsensusFeature c2;
  TEST_EQUAL(c2.size(), 0)
  c2.insert(1,tmp_feature2);
  c2.clear();
  TEST_EQUAL(c2.size(), 0)
}
END_SECTION

START_SECTION((bool empty() const))
{
  ConsensusFeature c1(tmp_feature);
  c1.insert(1, tmp_feature);
  c1.insert(2, tmp_feature3);
  TEST_EQUAL(c1.empty(), false)
  c1.clear();
  TEST_EQUAL(c1.empty(), true)

  ConsensusFeature c2;
  TEST_EQUAL(c2.size(), 0)
  c2.insert(1,tmp_feature2);
  TEST_EQUAL(c2.empty(), false)
  c2.clear();
  TEST_EQUAL(c2.empty(), true)
}
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
