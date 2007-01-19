// -*- Mode: C++; tab-width: 2; -*-
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
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringAffineSuperimposer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DPosition < 2, KernelTraits > PositionType;

START_TEST(PoseClusteringAffineSuperimposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PoseClusteringAffineSuperimposer<FeatureMap>* ptr = 0;
CHECK((PoseClusteringAffineSuperimposer()))
	ptr = new PoseClusteringAffineSuperimposer<FeatureMap>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((~PoseClusteringAffineSuperimposer()))
	delete ptr;
RESULT

CHECK((PoseClusteringAffineSuperimposer& operator = (const PoseClusteringAffineSuperimposer& source)))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat;
  
  pcat.setMzBucketSize(0.3);
  pcat.setShiftBucketSize(0,0.3);
  pcat.setShiftBucketSize(1,1.4);
  pcat.setScalingBucketSize(0,2.4);
  pcat.setScalingBucketSize(1,5.4);
  pcat.setBucketWindowShift(0,2);
  pcat.setBucketWindowShift(1,4);
  pcat.setBucketWindowScaling(0,6);
  pcat.setBucketWindowScaling(1,6);
  
  PoseClusteringAffineSuperimposer<FeatureMap> pcat_copy;
  pcat_copy = pcat;
  
  TEST_REAL_EQUAL(pcat_copy.getMzBucketSize(),0.3)
  TEST_REAL_EQUAL(pcat_copy.getShiftBucketSize(0),0.3)
  TEST_REAL_EQUAL(pcat_copy.getShiftBucketSize(1),1.4)
  TEST_REAL_EQUAL(pcat_copy.getScalingBucketSize(0),2.4)
  TEST_REAL_EQUAL(pcat_copy.getScalingBucketSize(1),5.4)
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowShift(0),2)
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowShift(1),4)
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowScaling(0),6)
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowScaling(1),6)
RESULT

CHECK((PoseClusteringAffineSuperimposer(const PoseClusteringAffineSuperimposer& source)))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat;
  
  pcat.setMzBucketSize(0.3);
  pcat.setShiftBucketSize(0,0.3);
  pcat.setShiftBucketSize(1,1.4);
  pcat.setScalingBucketSize(0,2.4);
  pcat.setScalingBucketSize(1,5.4);
  pcat.setBucketWindowShift(0,2);
  pcat.setBucketWindowShift(1,4);
  pcat.setBucketWindowScaling(0,6);
  pcat.setBucketWindowScaling(1,6);
  
  PoseClusteringAffineSuperimposer<FeatureMap> pcat_copy(pcat);
  
  TEST_REAL_EQUAL(pcat_copy.getMzBucketSize(),0.3)
  TEST_REAL_EQUAL(pcat_copy.getShiftBucketSize(0),0.3)
  TEST_REAL_EQUAL(pcat_copy.getShiftBucketSize(1),1.4)
  TEST_REAL_EQUAL(pcat_copy.getScalingBucketSize(0),2.4)
  TEST_REAL_EQUAL(pcat_copy.getScalingBucketSize(1),5.4)
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowShift(0),2)
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowShift(1),4)
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowScaling(0),6)
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowScaling(1),6)
RESULT

CHECK((static BaseSuperimposer<PointMapType>* create()))
  
RESULT

CHECK((static const String getName()))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat;
  
  TEST_EQUAL(pcat.getName() == "poseclustering_affine",true)
RESULT

CHECK((void run()))
  FeatureMap scene;
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
  
  FeatureMap modell;
  Feature feat3;
  Feature feat4;
  PositionType pos3(2.4,1.02);
  PositionType pos4(10.4,5.02);
  feat3.setPosition(pos3);
  feat3.setIntensity(100);
  feat4.setPosition(pos4);
  feat4.setIntensity(100);
  modell.push_back(feat3);
  modell.push_back(feat4);
  
  PoseClusteringAffineSuperimposer<FeatureMap> pcat;
  pcat.setElementMap(0,modell);
  pcat.setElementMap(1,scene);
  pcat.setScalingBucketSize(0,1);
  pcat.setScalingBucketSize(1,0.01);
  pcat.setShiftBucketSize(0,1);
  pcat.setShiftBucketSize(1,0.01);
  pcat.run();
    
  DLinearMapping<1> rt_mapping = pcat.getTransformation(0);
  DLinearMapping<1> mz_mapping = pcat.getTransformation(1);
  
  TEST_REAL_EQUAL(rt_mapping.getSlope(),2)
  TEST_REAL_EQUAL(rt_mapping.getIntercept(),0.4)
  TEST_REAL_EQUAL(mz_mapping.getSlope(),1)
  TEST_REAL_EQUAL(mz_mapping.getIntercept(),0.02)
RESULT

CHECK((UnsignedInt getBucketWindowScaling(UnsignedInt dim) const))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat_copy;
  
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowScaling(0),1)
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowScaling(1),1)
RESULT

CHECK((UnsignedInt getBucketWindowShift(UnsignedInt dim) const))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat_copy;
  
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowShift(0),1)
  TEST_REAL_EQUAL(pcat_copy.getBucketWindowShift(1),1)
RESULT

CHECK((double getMzBucketSize() const))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat_copy;
  
  TEST_REAL_EQUAL(pcat_copy.getMzBucketSize(),1)
RESULT

CHECK((double getScalingBucketSize(UnsignedInt dim) const))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat_copy;
  
  TEST_REAL_EQUAL(pcat_copy.getScalingBucketSize(0),0.5)
  TEST_REAL_EQUAL(pcat_copy.getScalingBucketSize(1),0.1)
RESULT

CHECK((double getShiftBucketSize(UnsignedInt dim) const))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat_copy;
  
  TEST_REAL_EQUAL(pcat_copy.getShiftBucketSize(0),1)
  TEST_REAL_EQUAL(pcat_copy.getShiftBucketSize(1),0.1)
RESULT

CHECK((void setBucketWindowScaling(UnsignedInt dim, UnsignedInt bucket_window_scaling)))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat;
  
  pcat.setBucketWindowScaling(0,6);
  pcat.setBucketWindowScaling(1,6);
  
  TEST_REAL_EQUAL(pcat.getBucketWindowScaling(0),6)
  TEST_REAL_EQUAL(pcat.getBucketWindowScaling(1),6)
RESULT

CHECK((void setBucketWindowShift(UnsignedInt dim, UnsignedInt bucket_window_shift)))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat;
  
  pcat.setBucketWindowShift(0,2);
  pcat.setBucketWindowShift(1,4);
  
  TEST_REAL_EQUAL(pcat.getBucketWindowShift(0),2)
  TEST_REAL_EQUAL(pcat.getBucketWindowShift(1),4)
RESULT

CHECK((void setMzBucketSize(double mz_bucket_size)))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat;
  
  pcat.setMzBucketSize(0.3);
  
  TEST_REAL_EQUAL(pcat.getMzBucketSize(),0.3)
RESULT

CHECK((void setScalingBucketSize(UnsignedInt dim, double scaling_bucket_size)))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat;
  
  pcat.setScalingBucketSize(0,2.4);
  pcat.setScalingBucketSize(1,5.4);
  
  TEST_REAL_EQUAL(pcat.getScalingBucketSize(0),2.4)
  TEST_REAL_EQUAL(pcat.getScalingBucketSize(1),5.4)
RESULT

CHECK((void setShiftBucketSize(UnsignedInt dim, double shift_bucket_size)))
  PoseClusteringAffineSuperimposer<FeatureMap> pcat;
  
  pcat.setShiftBucketSize(0,0.3);
  pcat.setShiftBucketSize(1,1.4);
  
  TEST_REAL_EQUAL(pcat.getShiftBucketSize(0),0.3)
  TEST_REAL_EQUAL(pcat.getShiftBucketSize(1),1.4)
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



