// -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>
#include <OpenMS/KERNEL/ConsensusMap.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/DelaunayPairFinder.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DPosition<2> PositionType;

START_TEST(DelaunayPairFinder, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

DelaunayPairFinder* ptr = 0;
CHECK((DelaunayPairFinder()))
	ptr = new DelaunayPairFinder();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~DelaunayPairFinder()))
	delete ptr;
RESULT

DelaunayPairFinder::Point* ptr2 = 0;
CHECK(([EXTRA]Point()))
	ptr2 = new DelaunayPairFinder::Point();
	TEST_NOT_EQUAL(ptr2, 0)
RESULT

CHECK(([EXTRA]~Point()))
	delete ptr2;
RESULT

CHECK(([EXTRA]Point& operator = (const Point& source)))
  DelaunayPairFinder::Point p(1,2);
  DelaunayPairFinder::Point p_copy;
  p_copy = p;
  	
  TEST_REAL_EQUAL(p_copy.hx(),1)
  TEST_REAL_EQUAL(p_copy.hy(),2)
RESULT

CHECK(([EXTRA]Point(Base::RT hx, Base::RT hy)))
  DelaunayPairFinder::Point p(1,2);
  DelaunayPairFinder::Point p_copy;
  p_copy = p;
  	
  TEST_REAL_EQUAL(p_copy.hx(),1)
  TEST_REAL_EQUAL(p_copy.hy(),2)
RESULT

CHECK(([EXTRA]Point(Base::RT hx, Base::RT hy, const PointType& f, UInt k=0)))
  DelaunayPairFinder::Point p(1,2);
  DelaunayPairFinder::Point p_copy;
  p_copy = p;
  	
  TEST_REAL_EQUAL(p_copy.hx(),1)
  TEST_REAL_EQUAL(p_copy.hy(),2)
RESULT

CHECK(([EXTRA]Point(const Base& cgal_point)))
  CGAL::Point_2< CGAL::Cartesian<double> > cp(1,2);
  DelaunayPairFinder::Point p(cp);
    	
  TEST_REAL_EQUAL(p.hx(),1)
  TEST_REAL_EQUAL(p.hy(),2)
RESULT

CHECK(([EXTRA]Point(const Point& source)))
  DelaunayPairFinder::Point p(1,2);
  DelaunayPairFinder::Point p_copy(p);
  	
  TEST_REAL_EQUAL(p_copy.hx(),1)
  TEST_REAL_EQUAL(p_copy.hy(),2)
RESULT

CHECK(([EXTRA]Point_2 operator()(const Circle_2& c) const))
  //
RESULT

CHECK((static BasePairFinder* create()))
	BasePairFinder* base_ptr = 0;
	base_ptr = DelaunayPairFinder::create();
	TEST_NOT_EQUAL(base_ptr, 0)
RESULT

CHECK((static const String getProductName()))
	DelaunayPairFinder dpf;
	
  TEST_EQUAL(dpf.getName() == "delaunay",true)
RESULT

#if 0

CHECK((void findElementPairs()))
  FeatureMap<>scene;
  Feature feat1;
  Feature feat2;
  Feature feat3;
  PositionType pos1(0,0);
  PositionType pos2(200,300);
  PositionType pos3(400,500);
  feat1.setPosition(pos1);
  feat1.setIntensity(100);
  feat2.setPosition(pos2);
  feat2.setIntensity(300);
  feat3.setPosition(pos3);
  feat3.setIntensity(400);
  scene.push_back(feat1);
  scene.push_back(feat2);
  scene.push_back(feat3);
  
  FeatureMap<>model;
  Feature feat4;
  Feature feat5;
  Feature feat6;
  PositionType pos4(4,4);
  PositionType pos5(204,304);
  PositionType pos6(404,504);
  feat4.setPosition(pos4);
  feat4.setIntensity(100);
  feat5.setPosition(pos5);
  feat5.setIntensity(300);
  feat6.setPosition(pos6);
  feat6.setIntensity(400);
  model.push_back(feat4);
  model.push_back(feat5);
  model.push_back(feat6);
  
  DelaunayPairFinder dpf;
  // dpf.setDiffIntercept(0,1.0); // TODO use internal_mz_scaling_ instead
  // dpf.setDiffIntercept(1,1.0); // TODO use internal_mz_scaling_ instead
  dpf.precision_[0]=5.0;
  dpf.precision_[1]=5.0;
	ConsensusMap model2;
	DelaunayPairFinder::convert(0,model,model2);
  dpf.setModelMap(model2);
	ConsensusMap scene2;
	DelaunayPairFinder::convert(1,scene,scene2);
  dpf.setSceneMap(scene2);
	DelaunayPairFinder::ElementPairVectorType pairs;
  dpf.setElementPairs(pairs);
  dpf.findElementPairs();
  
  TEST_EQUAL((pairs.begin())->first == feat1, true)
  TEST_EQUAL((pairs.begin())->second == feat4, true)
  TEST_EQUAL((pairs.begin()+1)->first == feat2, true)
  TEST_EQUAL((pairs.begin()+1)->second == feat5, true)
  TEST_EQUAL((pairs.begin()+2)->first == feat3,true)
  TEST_EQUAL((pairs.begin()+2)->second == feat6,true)
RESULT

#endif

CHECK(void run(ConsensusMap& result))
{
  ConsensusMap scene;
  Feature feat1;
  Feature feat2;
  Feature feat3;
  PositionType pos1(0,0);
  PositionType pos2(200,300);
  PositionType pos3(400,500);
  feat1.setPosition(pos1);
  feat1.setIntensity(100);
  feat2.setPosition(pos2);
  feat2.setIntensity(300);
  feat3.setPosition(pos3);
  feat3.setIntensity(400); 
  ConsensusFeature cons1(0,0,feat1);
  ConsensusFeature cons2(0,1,feat2);
  ConsensusFeature cons3(0,2,feat3);
  scene.push_back(cons1);
  scene.push_back(cons2);
  scene.push_back(cons3);
  
  ConsensusMap model;
  Feature feat4;
  Feature feat5;
  Feature feat6;
  PositionType pos4(4,4);
  PositionType pos5(204,304);
  PositionType pos6(404,504);
  feat4.setPosition(pos4);
  feat4.setIntensity(100);
  feat5.setPosition(pos5);
  feat5.setIntensity(300);
  feat6.setPosition(pos6);
  feat6.setIntensity(400);
  ConsensusFeature cons4(1,0,feat4);
  ConsensusFeature cons5(1,1,feat5);
  ConsensusFeature cons6(1,2,feat6);
  model.push_back(cons4);
  model.push_back(cons5);
  model.push_back(cons6);
  
  DelaunayPairFinder dpf;
	Param param = dpf.getDefaults();
	param.setValue("similarity:internal_mz_scaling",1.0);
	param.setValue("similarity:max_pair_distance:RT",5.0);
	param.setValue("similarity:max_pair_distance:MZ",5.0);
	param.setValue("similarity:precision:RT",5.0);
	param.setValue("similarity:precision:MZ",5.0);
	dpf.setParameters(param);
	ConsensusMap const& model_cref(model);
	ConsensusMap const& scene_cref(scene);
	dpf.setModelMap(0,model_cref);
	dpf.setSceneMap(1,scene_cref);
	ConsensusMap result;
	dpf.run(result);
	TEST_EQUAL(result.size(),3);
	ABORT_IF(result.size()!=3);

  ConsensusFeature::HandleSetType group1 = result[0].getFeatures();
  ConsensusFeature::HandleSetType group2 = result[1].getFeatures();
  ConsensusFeature::HandleSetType group3 = result[2].getFeatures();
  
  FeatureHandle ind1(0,0,feat1);
  FeatureHandle ind2(0,1,feat2);
  FeatureHandle ind3(0,2,feat3);
  FeatureHandle ind4(1,0,feat4);
  FeatureHandle ind5(1,1,feat5);
  FeatureHandle ind6(1,2,feat6);

  ConsensusFeature::HandleSetType::const_iterator it;
	it = group1.begin();
  STATUS(*it);
	STATUS(ind1);
	TEST_EQUAL(*(it) == ind1, true)
	++it;
  STATUS(*it);
	STATUS(ind4);
  TEST_EQUAL(*(it) == ind4, true)
	it = group2.begin();
  STATUS(*it);
	STATUS(ind2);
  TEST_EQUAL(*(it) == ind2, true)
	++it;
  STATUS(*it);
	STATUS(ind5);
  TEST_EQUAL(*(it) == ind5, true)
  it = group3.begin();
  STATUS(*it);
	STATUS(ind3);
  TEST_EQUAL(*(it) == ind3, true)
	++it;
  STATUS(*it);
	STATUS(ind6);
  TEST_EQUAL(*(it) == ind6, true)
}
RESULT


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



