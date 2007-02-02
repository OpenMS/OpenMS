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
// $Maintainer: Clemens Groepl, Eva Lange $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringShiftSuperimposer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DPosition < 2, KernelTraits > PositionType;

START_TEST(PoseClusteringShiftSuperimposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PoseClusteringShiftSuperimposer<FeatureMap>* ptr = 0;
CHECK(PoseClusteringShiftSuperimposer())
	ptr = new PoseClusteringShiftSuperimposer<FeatureMap>();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK(~PoseClusteringShiftSuperimposer())
	delete ptr;
RESULT

PoseClusteringShiftSuperimposer<FeatureMap>::Shift* ptr2;
CHECK(Shift())
  ptr2 = new PoseClusteringShiftSuperimposer<FeatureMap>::Shift();
  TEST_NOT_EQUAL(ptr2, 0)
RESULT

CHECK(~Shift())
  delete ptr2;
RESULT

CHECK(Shift& operator= (Shift const & source))
  PoseClusteringShiftSuperimposer<FeatureMap>::Shift shift;
  shift.getQuality() = 0.1;
  shift.getPosition()[0] = 0.2;
  shift.getPosition()[1] = 0.4;
  
  PoseClusteringShiftSuperimposer<FeatureMap>::Shift shift_copy;
  shift_copy = shift;
  
  TEST_REAL_EQUAL(shift_copy.getQuality(), 0.1)
  TEST_REAL_EQUAL(shift_copy.getPosition()[0], 0.2)
  TEST_REAL_EQUAL(shift_copy.getPosition()[1], 0.4)
RESULT

CHECK(Shift(Shift const & source))
  PoseClusteringShiftSuperimposer<FeatureMap>::Shift shift;
  shift.getQuality() = 0.1;
  shift.getPosition()[0] = 0.2;
  shift.getPosition()[1] = 0.4;
  
  PoseClusteringShiftSuperimposer<FeatureMap>::Shift shift_copy(shift);
  
  TEST_REAL_EQUAL(shift_copy.getQuality(), 0.1)
  TEST_REAL_EQUAL(shift_copy.getPosition()[0], 0.2)
  TEST_REAL_EQUAL(shift_copy.getPosition()[1], 0.4)
RESULT

CHECK(QualityType& getQuality())
  PoseClusteringShiftSuperimposer<FeatureMap>::Shift shift;
  shift.getQuality() = 0.1;
  
  TEST_REAL_EQUAL(shift.getQuality(), 0.1)
RESULT


CHECK(PositionType& getPosition())
  PoseClusteringShiftSuperimposer<FeatureMap>::Shift shift;
  PositionType pos;
  pos[0] = 0.2;
  pos[1] = 0.4;
  shift.getPosition() = pos;
  
  TEST_REAL_EQUAL(shift.getPosition()[0], 0.2)
  TEST_REAL_EQUAL(shift.getPosition()[1], 0.4)
RESULT

CHECK(const PositionType& getPosition() const)
  PoseClusteringShiftSuperimposer<FeatureMap>::Shift shift;
  
  TEST_REAL_EQUAL(shift.getPosition()[0], 0)
  TEST_REAL_EQUAL(shift.getPosition()[1], 0)
RESULT

CHECK(const QualityType& getQuality() const)
  PoseClusteringShiftSuperimposer<FeatureMap>::Shift shift;
  
  TEST_REAL_EQUAL(shift.getQuality(), 0)
RESULT

CHECK(void setPosition(const PositionType& position))
  PoseClusteringShiftSuperimposer<FeatureMap>::Shift shift;
  PositionType pos;
  pos[0] = 0.2;
  pos[1] = 0.4;
  shift.setPosition(pos);
  
  TEST_REAL_EQUAL(shift.getPosition()[0], 0.2)
  TEST_REAL_EQUAL(shift.getPosition()[1], 0.4)
RESULT

CHECK(void setQuality(const QualityType& quality))
  PoseClusteringShiftSuperimposer<FeatureMap>::Shift shift;
  shift.setQuality(0.1);
  
  TEST_REAL_EQUAL(shift.getQuality(), 0.1)
RESULT

CHECK(PoseClusteringShiftSuperimposer(const PoseClusteringShiftSuperimposer& source))
  PoseClusteringShiftSuperimposer<FeatureMap> pcsi;
  pcsi.setShiftBucketSize(0,1.9);
  pcsi.setShiftBucketSize(1,2.9);
  pcsi.setElementBucketWindow(0,3);
  pcsi.setElementBucketWindow(1,4);
  pcsi.setShiftBucketWindow(0,4);
  pcsi.setShiftBucketWindow(1,5);
  
  cout << pcsi.getParameters() << endl;
  
  PoseClusteringShiftSuperimposer<FeatureMap> pcsi_copy(pcsi);
  
   cout << pcsi_copy.getParameters() << endl;
    
  TEST_REAL_EQUAL(pcsi_copy.getShiftBucketSize(0),1.9)
  TEST_REAL_EQUAL(pcsi_copy.getShiftBucketSize(1),2.9)
  TEST_REAL_EQUAL(pcsi_copy.getElementBucketWindow(0),3)
  TEST_REAL_EQUAL(pcsi_copy.getElementBucketWindow(1),4)
  TEST_REAL_EQUAL(pcsi_copy.getShiftBucketWindow(0),4)
  TEST_REAL_EQUAL(pcsi_copy.getShiftBucketWindow(1),5)
RESULT

CHECK(PoseClusteringAffineSuperimposer& operator = (const PoseClusteringAffineSuperimposer& source))
  PoseClusteringShiftSuperimposer<FeatureMap> pcsi;
  pcsi.setShiftBucketSize(0,1.9);
  pcsi.setShiftBucketSize(1,2.9);
  pcsi.setElementBucketWindow(0,3);
  pcsi.setElementBucketWindow(1,4);
  pcsi.setShiftBucketWindow(0,4);
  pcsi.setShiftBucketWindow(1,5);

  cout << pcsi.getParameters() << endl;
  
  PoseClusteringShiftSuperimposer<FeatureMap> pcsi_copy;
  pcsi_copy = pcsi;
  
  cout << pcsi_copy.getParameters() << endl;
  
  TEST_REAL_EQUAL(pcsi_copy.getShiftBucketSize(0),1.9)
  TEST_REAL_EQUAL(pcsi_copy.getShiftBucketSize(1),2.9)
  TEST_REAL_EQUAL(pcsi_copy.getElementBucketWindow(0),3)
  TEST_REAL_EQUAL(pcsi_copy.getElementBucketWindow(1),4)
  TEST_REAL_EQUAL(pcsi_copy.getShiftBucketWindow(0),4)
  TEST_REAL_EQUAL(pcsi_copy.getShiftBucketWindow(1),5)  
RESULT

CHECK(Size getElementBucketWindow(UnsignedInt dim) const)
  PoseClusteringShiftSuperimposer<FeatureMap> pcsi;
    
  TEST_REAL_EQUAL(pcsi.getElementBucketWindow(0),2)
  TEST_REAL_EQUAL(pcsi.getElementBucketWindow(1),1)
RESULT

CHECK(Size getShiftBucketWindow(UnsignedInt dim) const)
  PoseClusteringShiftSuperimposer<FeatureMap> pcsi;
      
  TEST_REAL_EQUAL(pcsi.getShiftBucketWindow(0),2)
  TEST_REAL_EQUAL(pcsi.getShiftBucketWindow(1),1)  
RESULT

CHECK(double getShiftBucketSize(UnsignedInt dim) const)
  PoseClusteringShiftSuperimposer<FeatureMap> pcsi;
    
  TEST_REAL_EQUAL(pcsi.getShiftBucketSize(0),5.0)
  TEST_REAL_EQUAL(pcsi.getShiftBucketSize(1),0.1)
RESULT

CHECK(static BaseSuperimposer<PointMapType>* create())
  
RESULT

CHECK(static const String getName())
  PoseClusteringShiftSuperimposer<FeatureMap> pcsi;
  
  TEST_EQUAL(pcsi.getName() == "poseclustering_shift",true)
RESULT

CHECK(void run())
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
  PositionType pos3(21.4,1.02);
  PositionType pos4(25.4,5.02);
  feat3.setPosition(pos3);
  feat3.setIntensity(100);
  feat4.setPosition(pos4);
  feat4.setIntensity(100);
  modell.push_back(feat3);
  modell.push_back(feat4);
  
  PoseClusteringShiftSuperimposer<FeatureMap> pcat;
  pcat.setElementMap(0,modell);
  pcat.setElementMap(1,scene);
  pcat.run();
    
  DLinearMapping<1> rt_mapping = pcat.getTransformation(0);
  DLinearMapping<1> mz_mapping = pcat.getTransformation(1);
  
  TEST_REAL_EQUAL(rt_mapping.getSlope(),1)
  TEST_REAL_EQUAL(rt_mapping.getIntercept(),20.4)
  TEST_REAL_EQUAL(mz_mapping.getSlope(),1)
  TEST_REAL_EQUAL(mz_mapping.getIntercept(),0.02)
RESULT

CHECK((void setElementBucketWindow(UnsignedInt dim, Size element_bucket_window)))
  PoseClusteringShiftSuperimposer<FeatureMap> pcsi;
  pcsi.setElementBucketWindow(0,3);
  pcsi.setElementBucketWindow(1,4);
    
  TEST_REAL_EQUAL(pcsi.getElementBucketWindow(0),3)
  TEST_REAL_EQUAL(pcsi.getElementBucketWindow(1),4)
RESULT

CHECK((void setShiftBucketSize(UnsignedInt dim, double shift_bucket_size)))
  PoseClusteringShiftSuperimposer<FeatureMap> pcsi;
  pcsi.setShiftBucketSize(0,1.9);
  pcsi.setShiftBucketSize(1,2.9);
    
  TEST_REAL_EQUAL(pcsi.getShiftBucketSize(0),1.9)
  TEST_REAL_EQUAL(pcsi.getShiftBucketSize(1),2.9)
RESULT

CHECK((void setShiftBucketWindow(UnsignedInt dim, Size shift_bucket_window)))
  PoseClusteringShiftSuperimposer<FeatureMap> pcsi;
  pcsi.setShiftBucketWindow(0,4);
  pcsi.setShiftBucketWindow(1,5);
  
  TEST_REAL_EQUAL(pcsi.getShiftBucketWindow(0),4)
  TEST_REAL_EQUAL(pcsi.getShiftBucketWindow(1),5)  
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



