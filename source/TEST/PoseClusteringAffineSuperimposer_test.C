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

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DPosition <2> PositionType;

START_TEST(PoseClusteringAffineSuperimposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PoseClusteringAffineSuperimposer<FeatureMap<> >* ptr = 0;
CHECK((PoseClusteringAffineSuperimposer()))
	ptr = new PoseClusteringAffineSuperimposer<FeatureMap<> >();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~PoseClusteringAffineSuperimposer()))
	delete ptr;
RESULT

CHECK((static BaseSuperimposer<PointMapType>* create()))
  
RESULT

CHECK((static const String getProductName()))
  PoseClusteringAffineSuperimposer<FeatureMap<> > pcat;
  
  TEST_EQUAL(pcat.getName() == "poseclustering_affine",true)
RESULT

CHECK((virtual void run()))
  FeatureMap<> scene;
  Feature feat1;
  Feature feat2;
  PositionType pos1(1,1);
  PositionType pos2(5,5);
  feat1.setPosition(pos1);
  feat1.setIntensity(100);
  feat2.setPosition(pos2);
  feat2.setIntensity(100);
  scene.push_back(feat1);
  scene.push_back(feat2);
  
  FeatureMap<> modell;
  Feature feat3;
  Feature feat4;
  PositionType pos3(1.4,1.02);
  PositionType pos4(5.4,5.02);
  feat3.setPosition(pos3);
  feat3.setIntensity(100);
  feat4.setPosition(pos4);
  feat4.setIntensity(100);
  modell.push_back(feat3);
  modell.push_back(feat4);

  Param parameters;
	parameters.setValue(String("transformation_space:scaling_bucket_size:RT"), 0.01);
	parameters.setValue(String("transformation_space:scaling_bucket_size:MZ"), 0.01);
	parameters.setValue(String("transformation_space:shift_bucket_size:RT"), 0.01);
	parameters.setValue(String("transformation_space:shift_bucket_size:MZ"), 0.01);
	parameters.setValue(String("transformation_space:bucket_window_scaling:RT"), 0);
	parameters.setValue(String("transformation_space:bucket_window_scaling:MZ"), 0);
  
  PoseClusteringAffineSuperimposer<FeatureMap<> > pcat;
  pcat.setModelMap(modell);
  pcat.setSceneMap(scene);
  pcat.setParameters(parameters);
  
  LinearMapping rt_mapping;
  pcat.run(rt_mapping);
    
  TEST_REAL_EQUAL(rt_mapping.getSlope(),1.0)
  TEST_REAL_EQUAL(rt_mapping.getIntercept(),0.4)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



