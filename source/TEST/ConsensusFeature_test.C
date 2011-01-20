// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Erhan Kenar $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/KERNEL/ConsensusFeature.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/ElementDB.h>
#include <OpenMS/CHEMISTRY/Element.h>

using namespace OpenMS;
using namespace std;

START_TEST(ConsensusFeature, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConsensusFeature* ptr = 0;
START_SECTION((ConsensusFeature()))
	ptr = new ConsensusFeature();
	TEST_NOT_EQUAL(ptr, 0)
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

  ConsensusFeature cons(tmp_feature);
  cons.insert(1,tmp_feature);
  ConsensusFeature cons_copy(cons);

  TEST_REAL_SIMILAR(cons_copy.getRT(),1.0)
  TEST_REAL_SIMILAR(cons_copy.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons_copy.getIntensity(),200.0)
  TEST_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_EQUAL((cons_copy.begin())->getUniqueId(),3)
  TEST_EQUAL((cons_copy.begin())->getIntensity(),200)
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


START_SECTION((ConsensusFeature(const Peak2D &point)))
  ConsensusFeature cons(static_cast<Peak2D>(tmp_feature));
  TEST_REAL_SIMILAR(cons.getRT(),1.0)
  TEST_REAL_SIMILAR(cons.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons.getIntensity(),200.0)
  TEST_EQUAL(cons.empty(), true)
END_SECTION

START_SECTION((ConsensusFeature(const RichPeak2D &point)))

  ConsensusFeature cons(static_cast<RichPeak2D>(tmp_feature));
  TEST_REAL_SIMILAR(cons.getRT(),1.0)
  TEST_REAL_SIMILAR(cons.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons.getIntensity(),200.0)
  TEST_EQUAL(cons.empty(), true)
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

START_SECTION((void computeDechargeConsensus(const FeatureMap<>& fm, bool intensity_weighted_averaging=false)))
  
  DoubleReal proton_mass = ElementDB::getInstance()->getElement("H")->getMonoWeight();
  DoubleReal natrium_mass = ElementDB::getInstance()->getElement("Na")->getMonoWeight();
  
  DoubleReal m = 1000;
  DoubleReal m1_add = 0.5;
  DoubleReal mz1 = (m+m1_add+3*proton_mass) / 3;
  DoubleReal m2_add = 1;
  DoubleReal mz2 = (m+m2_add+1*proton_mass + 2*natrium_mass) / 3;
  DoubleReal m3_add = -0.5;
  DoubleReal mz3 = (m+m3_add+4*proton_mass + natrium_mass) / 5;
  
  FeatureMap<> fm;
  
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
