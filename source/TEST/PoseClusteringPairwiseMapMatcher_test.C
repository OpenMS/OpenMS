// // -*- Mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2007 -- Oliver Kohlbacher, Knut Reinert
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

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringPairwiseMapMatcher.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DFeature<2, KernelTraits> ElementType;
typedef DFeatureMap<2, ElementType> ElementMapType;
typedef DFeaturePair < 2, ElementType > ElementPairType;
typedef DFeaturePairVector < 2, ElementType > ElementPairVectorType;
typedef DPosition < 2, KernelTraits > PositionType;

START_TEST(PoseClusteringPairwiseMapMatcher<ElementMapType>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PoseClusteringPairwiseMapMatcher<ElementMapType>* ptr = 0;
CHECK((PoseClusteringPairwiseMapMatcher()))
	ptr = new PoseClusteringPairwiseMapMatcher<ElementMapType>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~PoseClusteringPairwiseMapMatcher()))
	delete ptr;
RESULT

CHECK((PoseClusteringPairwiseMapMatcher& operator= (const PoseClusteringPairwiseMapMatcher& source)))
  Param param;
  param.setValue("bla",3);
  ElementMapType first;
  ElementMapType second;
  
  PoseClusteringPairwiseMapMatcher<ElementMapType> pcpmm;
  pcpmm.setParameters(param);
  pcpmm.setElementMap(0,first);
  pcpmm.setElementMap(1,second);
  
  PoseClusteringPairwiseMapMatcher<ElementMapType> pcpmm_copy;
  pcpmm_copy = pcpmm;
  
  TEST_EQUAL(pcpmm.getParameters() == pcpmm_copy.getParameters(),true)
  TEST_EQUAL(&(pcpmm.getElementMap(0)) == &(pcpmm_copy.getElementMap(0)),true)
  TEST_EQUAL(&(pcpmm.getElementMap(1)) == &(pcpmm_copy.getElementMap(1)),true)
RESULT

CHECK((PoseClusteringPairwiseMapMatcher(const PoseClusteringPairwiseMapMatcher& source)))
  Param param;
  param.setValue("bla",3);
  ElementMapType first;
  ElementMapType second;
  
  PoseClusteringPairwiseMapMatcher<ElementMapType> pcpmm;
  pcpmm.setParameters(param);
  pcpmm.setElementMap(0,first);
  pcpmm.setElementMap(1,second);
  
  PoseClusteringPairwiseMapMatcher<ElementMapType> pcpmm_copy(pcpmm);
  
  TEST_EQUAL(pcpmm.getParameters() == pcpmm_copy.getParameters(),true)
  TEST_EQUAL(&(pcpmm.getElementMap(0)) == &(pcpmm_copy.getElementMap(0)),true)
  TEST_EQUAL(&(pcpmm.getElementMap(1)) == &(pcpmm_copy.getElementMap(1)),true)
RESULT

CHECK((static BasePairwiseMapMatcher<MapT>* create()))
  
RESULT

CHECK((static const String getName()))
  PoseClusteringPairwiseMapMatcher<ElementMapType> pcpmm;
  
  TEST_EQUAL(pcpmm.getName() == "poseclustering_pairwise",true)
RESULT

CHECK((void run()))
  Param param;
  param.setValue("superimposer","poseclustering_shift");
  param.setValue("pair_finder","simple");
  ElementMapType scene;
  ElementType feat1;
  ElementType feat2;
  PositionType pos1(0,0);
  PositionType pos2(200,300);
  feat1.setPosition(pos1);
  feat1.setIntensity(100);
  feat2.setPosition(pos2);
  feat2.setIntensity(300);
  scene.push_back(feat1);
  scene.push_back(feat2);
  ElementMapType modell = scene;
  
  ElementType feat3;
  ElementType feat4;
  PositionType pos3(2,5);
  PositionType pos4(20,30);
  feat3.setPosition(pos3);
  feat3.setIntensity(100);
  feat4.setPosition(pos4);
  feat4.setIntensity(300);
  scene.push_back(feat3);
  modell.push_back(feat4);
  
  PoseClusteringPairwiseMapMatcher<ElementMapType> pcpmm;
  pcpmm.setParameters(param);
  pcpmm.setElementMap(0,modell);
  pcpmm.setElementMap(1,scene);
  pcpmm.initGridTransformation(scene);
  pcpmm.run();
  
  TEST_EQUAL((pcpmm.getElementPairs().begin())->first == (pcpmm.getElementPairs().begin())->second,true)
  TEST_EQUAL((pcpmm.getElementPairs().begin()+1)->first == (pcpmm.getElementPairs().begin()+1)->second,true)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



