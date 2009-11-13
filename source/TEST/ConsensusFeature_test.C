// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Clemens Groepl $
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

Feature tmp_feature2;
tmp_feature2.setRT(2);
tmp_feature2.setMZ(3);
tmp_feature2.setIntensity(300.0f);

Feature tmp_feature3;
tmp_feature3.setRT(3);
tmp_feature3.setMZ(4);
tmp_feature3.setIntensity(400.0f);

START_SECTION((ConsensusFeature& operator=(const ConsensusFeature &rhs)))
  ConsensusFeature cons(tmp_feature);
  cons.insert(1,3,tmp_feature);

  ConsensusFeature cons_copy;
  cons_copy = cons;

  TEST_REAL_SIMILAR(cons_copy.getRT(),1.0)
  TEST_REAL_SIMILAR(cons_copy.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons_copy.getIntensity(),200.0)
  TEST_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_EQUAL((cons_copy.begin())->getElementIndex(),3)
  TEST_EQUAL((cons_copy.begin())->getIntensity(),200)
END_SECTION

START_SECTION((ConsensusFeature(const ConsensusFeature &rhs)))

  ConsensusFeature cons(tmp_feature);
  cons.insert(1,3,tmp_feature);
  ConsensusFeature cons_copy(cons);

  TEST_REAL_SIMILAR(cons_copy.getRT(),1.0)
  TEST_REAL_SIMILAR(cons_copy.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons_copy.getIntensity(),200.0)
  TEST_EQUAL((cons_copy.begin())->getMapIndex(),1)
  TEST_EQUAL((cons_copy.begin())->getElementIndex(),3)
  TEST_EQUAL((cons_copy.begin())->getIntensity(),200)
END_SECTION

START_SECTION((void insert(const HandleSetType &handle_set)))
  ConsensusFeature::HandleSetType hs;
  FeatureHandle fh;
  for ( UInt i = 0; i < 3; ++i )
  {
    fh.setRT(i*77.7);
    fh.setMapIndex(i+10);
    fh.setElementIndex(i+1000);
    hs.insert(fh);
  }
  ConsensusFeature cf;
  cf.insert(hs);

  TEST_EQUAL(cf.size(),3);
  TEST_EQUAL(cf.begin()->getMapIndex(),10);
  TEST_EQUAL(cf.rbegin()->getMapIndex(),12);

END_SECTION

START_SECTION((void insert(UInt64 map_index, UInt64 element_index, const Peak2D &element)))
  ConsensusFeature cf;
  Peak2D el;
  for ( UInt i = 0; i < 3; ++i )
  {
    el.setRT(i*77.7);
    cf.insert(10-i,i+1000,el);
    TEST_EQUAL(cf.size(),i+1);
    TEST_REAL_SIMILAR(cf.begin()->getRT(),i*77.7);
    TEST_EQUAL(cf.begin()->getMapIndex(),10-i);
    TEST_EQUAL(cf.begin()->getElementIndex(),i+1000);
  }

END_SECTION

START_SECTION((void insert(UInt64 map_index, UInt64 element_index, const ConsensusFeature &element)))
  ConsensusFeature cf;
  ConsensusFeature el;
  for ( UInt i = 0; i < 3; ++i )
  {
    el.setRT(i*77.7);
    el.setCharge(2*i);
    cf.insert(10-i,i+1000,el);
    TEST_EQUAL(cf.size(),i+1);
    TEST_REAL_SIMILAR(cf.begin()->getRT(),i*77.7);
    TEST_EQUAL(cf.begin()->getCharge(),2*i);
    TEST_EQUAL(cf.begin()->getMapIndex(),10-i);
    TEST_EQUAL(cf.begin()->getElementIndex(),i+1000);
  }

END_SECTION


START_SECTION((ConsensusFeature(const Peak2D &point)))

  ConsensusFeature cons(tmp_feature);
  TEST_REAL_SIMILAR(cons.getRT(),1.0)
  TEST_REAL_SIMILAR(cons.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons.getIntensity(),200.0)
  TEST_EQUAL(cons.empty(), true)
END_SECTION

START_SECTION((ConsensusFeature(const RichPeak2D &point)))

  ConsensusFeature cons(tmp_feature);
  TEST_REAL_SIMILAR(cons.getRT(),1.0)
  TEST_REAL_SIMILAR(cons.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons.getIntensity(),200.0)
  TEST_EQUAL(cons.empty(), true)
END_SECTION

START_SECTION((ConsensusFeature(const Feature &feature)))
  Feature f;
  f.setCharge(-17);
  f.setRT(44324.6);
  f.setMZ(867.4);
  const Feature& f_cref = f;
  ConsensusFeature cf(99,23,f_cref);

  TEST_EQUAL(cf.getRT(),44324.6);
  TEST_EQUAL(cf.getMZ(),867.4);
  TEST_EQUAL(cf.getCharge(),-17);
END_SECTION

START_SECTION((ConsensusFeature(UInt64 map_index, UInt64 element_index, const Feature &element)))
 	ConsensusFeature cons(1,3,tmp_feature);

  TEST_REAL_SIMILAR(cons.getRT(),1.0)
  TEST_REAL_SIMILAR(cons.getMZ(),2.0)
  TEST_REAL_SIMILAR(cons.getIntensity(),200.0)
  ConsensusFeature::HandleSetType::const_iterator it = cons.begin();
  TEST_EQUAL(it->getMapIndex(),1)
  TEST_EQUAL(it->getElementIndex(),3)
  TEST_EQUAL(it->getIntensity(),200)
END_SECTION

START_SECTION((ConsensusFeature(UInt64 map_index, UInt64 element_index, const Peak2D &element)))
  Peak2D f;
  f.setIntensity(-17);
  const Peak2D& f_cref = f;
  ConsensusFeature cf(99,23,f_cref);

  ConsensusFeature::HandleSetType::const_iterator it = cf.begin();
  TEST_EQUAL(it->getMapIndex(),99);
  TEST_EQUAL(it->getElementIndex(),23);
  TEST_EQUAL(it->getIntensity(),-17);
END_SECTION

START_SECTION((ConsensusFeature(UInt64 map_index, UInt64 element_index, const ConsensusFeature &element)))
  ConsensusFeature f;
  f.setIntensity(-17);
  const ConsensusFeature& f_cref = f;
  ConsensusFeature cf(99,23,f_cref);

  ConsensusFeature::HandleSetType::const_iterator it = cf.begin();
  TEST_EQUAL(it->getMapIndex(),99);
  TEST_EQUAL(it->getElementIndex(),23);
  TEST_EQUAL(it->getIntensity(),-17);
END_SECTION

START_SECTION((DRange<1> getIntensityRange() const))
  ConsensusFeature cons;
  Feature f;
  f.setIntensity(0.0f);
  cons.insert(0,0,f);
  f.setIntensity(200.0f);
  cons.insert(0,1,f);

  TEST_REAL_SIMILAR(cons.getIntensityRange().minX(),0.0)
  TEST_REAL_SIMILAR(cons.getIntensityRange().maxX(),200.0)
END_SECTION

START_SECTION((DRange<2> getPositionRange() const))
  ConsensusFeature cons;
  Feature f;
  f.setRT(1.0);
  f.setMZ(500.0);
  cons.insert(0,0,f);
  f.setRT(1000.0);
  f.setMZ(1500.0);
  cons.insert(0,1,f);

  TEST_REAL_SIMILAR(cons.getPositionRange().minX(),1.0)
  TEST_REAL_SIMILAR(cons.getPositionRange().maxX(),1000.0)
  TEST_REAL_SIMILAR(cons.getPositionRange().minY(),500.0)
  TEST_REAL_SIMILAR(cons.getPositionRange().maxY(),1500.0)
END_SECTION

START_SECTION((const HandleSetType& getFeatures() const))
  ConsensusFeature cons;
  cons.insert(2,3,tmp_feature);
  const ConsensusFeature cons_copy(cons);

  ConsensusFeature::HandleSetType group = cons_copy.getFeatures();

  ConsensusFeature::HandleSetType::const_iterator it = group.begin();
  TEST_EQUAL(it->getMapIndex(),2)
  TEST_EQUAL(it->getElementIndex(),3)
  TEST_EQUAL(it->getIntensity(),200)
END_SECTION


START_SECTION((void insert(const FeatureHandle &handle)))
  ConsensusFeature cons;
  FeatureHandle h1(2,3,tmp_feature);
  FeatureHandle h2(4,5,tmp_feature);
  cons.insert(h1);
  cons.insert(h2);

  ConsensusFeature::HandleSetType::const_iterator it = cons.begin();
  TEST_EQUAL(it->getMapIndex(),2)
  TEST_EQUAL(it->getElementIndex(),3)
  TEST_EQUAL(it->getIntensity(),200)
  ++it;
  TEST_EQUAL(it->getMapIndex(),4)
  TEST_EQUAL(it->getElementIndex(),5)
  TEST_EQUAL(it->getIntensity(),200)
  ++it;
  TEST_EQUAL(it==cons.end(), true)
END_SECTION

START_SECTION((void insert(UInt64 map_index, UInt64 element_index, const Feature &element)))
  ConsensusFeature cons;
  cons.insert(2,3,tmp_feature);

  ConsensusFeature::HandleSetType::const_iterator it = cons.begin();
  TEST_EQUAL(it->getMapIndex(),2)
  TEST_EQUAL(it->getElementIndex(),3)
  TEST_EQUAL(it->getIntensity(),200)
  ++it;
  TEST_EQUAL(it==cons.end(),true)
END_SECTION

START_SECTION((QualityType getQuality() const))
	ConsensusFeature cons;
	TEST_EQUAL(cons.getQuality(),0.0)
END_SECTION

START_SECTION((void setQuality(QualityType quality)))
	ConsensusFeature cons;
	cons.setQuality(4.5);
	TEST_REAL_SIMILAR(cons.getQuality(),4.5);
END_SECTION

START_SECTION((Int getCharge() const))
  ConsensusFeature const cons;
  TEST_EQUAL(cons.getCharge(),0);
END_SECTION

START_SECTION((void setCharge(Int charge)))
  ConsensusFeature cons;
  cons.setCharge(-567);
  TEST_EQUAL(cons.getCharge(),-567);
END_SECTION

START_SECTION((const std::vector<PeptideIdentification>& getPeptideIdentifications() const))
  ConsensusFeature const cons;
  TEST_EQUAL(cons.getPeptideIdentifications().size(),0);
END_SECTION

START_SECTION((std::vector<PeptideIdentification>& getPeptideIdentifications()))
  ConsensusFeature cons;
  cons.getPeptideIdentifications().resize(9);
  TEST_EQUAL(cons.getPeptideIdentifications().size(),9);
END_SECTION

START_SECTION((void setPeptideIdentifications(const std::vector< PeptideIdentification > &peptide_identifications)))
  ConsensusFeature cons;
  std::vector<PeptideIdentification> pi(99);
  cons.setPeptideIdentifications(pi);
  TEST_EQUAL(cons.getPeptideIdentifications().size(),99);
END_SECTION


START_SECTION((void computeConsensus()))
  ConsensusFeature cons;
  //one point
  cons.insert(2,3,tmp_feature);
	cons.computeConsensus();
	TEST_REAL_SIMILAR(cons.getIntensity(),200.0)
	TEST_REAL_SIMILAR(cons.getRT(),1.0)
	TEST_REAL_SIMILAR(cons.getMZ(),2.0)
	//two points
  cons.insert(4,5,tmp_feature2);
	cons.computeConsensus();
	TEST_REAL_SIMILAR(cons.getIntensity(),250.0)
	TEST_REAL_SIMILAR(cons.getRT(),1.5)
	TEST_REAL_SIMILAR(cons.getMZ(),2.5)
	//three points
  cons.insert(6,7,tmp_feature3);
	cons.computeConsensus();
	TEST_REAL_SIMILAR(cons.getIntensity(),300.0)
	TEST_REAL_SIMILAR(cons.getRT(),2.0)
	TEST_REAL_SIMILAR(cons.getMZ(),3.0)
END_SECTION

START_SECTION((void computeDechargeConsensus(const FeatureMap<>& fm)))
  
  DoubleReal proton_mass = ElementDB::getInstance()->getElement("H")->getMonoWeight();
  DoubleReal natrium_mass = ElementDB::getInstance()->getElement("Na")->getMonoWeight();
  
  DoubleReal m = 1000;
  DoubleReal mz1 = (m+3*proton_mass) / 3;
  DoubleReal mz2 = (m+1*proton_mass + 2*natrium_mass) / 3;
  DoubleReal mz3 = (m+4*proton_mass + natrium_mass) / 5;
  
  FeatureMap<> fm;
  
  //one point  
  ConsensusFeature cons;
  Feature tmp_feature;
	tmp_feature.setRT(100);
	tmp_feature.setMZ(mz1);
	tmp_feature.setIntensity(200.0f);
	tmp_feature.setCharge(3);
	fm.push_back(tmp_feature);
  cons.insert(2,0,tmp_feature);
	cons.computeDechargeConsensus(fm);
	TEST_REAL_SIMILAR(cons.getIntensity(),200.0)
	TEST_REAL_SIMILAR(cons.getRT(),100)
	TEST_REAL_SIMILAR(cons.getMZ(), m);
	
	//two points
  Feature tmp_feature2;
	tmp_feature2.setRT(102);
	tmp_feature2.setMZ(mz2);
	tmp_feature2.setIntensity(250.0f);
	tmp_feature2.setCharge(3);
	tmp_feature2.setMetaValue("dc_charge_adduct_mass", 2*natrium_mass + proton_mass);
	fm.push_back(tmp_feature2);
	cons.insert(4,1,tmp_feature2);
	cons.computeDechargeConsensus(fm);
	TEST_REAL_SIMILAR(cons.getIntensity(),450.0)
	TEST_REAL_SIMILAR(cons.getRT(),101)
	TEST_REAL_SIMILAR(cons.getMZ(), m)
	
	//three points
  Feature tmp_feature3;
	tmp_feature3.setRT(101);
	tmp_feature3.setMZ(mz3);
	tmp_feature3.setIntensity(350.0f);
	tmp_feature3.setCharge(5);
	tmp_feature3.setMetaValue("dc_charge_adduct_mass", 1*natrium_mass + 4*proton_mass);
	fm.push_back(tmp_feature3);
	cons.insert(4,2,tmp_feature3);
	cons.computeDechargeConsensus(fm);
	TEST_REAL_SIMILAR(cons.getIntensity(),800.0)
	TEST_REAL_SIMILAR(cons.getRT(),101)
	TEST_REAL_SIMILAR(cons.getMZ(), m)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



