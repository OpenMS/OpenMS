// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2012 -- Oliver Kohlbacher, Knut Reinert
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
// $Authors: Eva Lange, Clemens Groepl $
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

PoseClusteringAffineSuperimposer* ptr = 0;
PoseClusteringAffineSuperimposer* nullPointer = 0;
BaseSuperimposer* base_nullPointer = 0;

START_SECTION((PoseClusteringAffineSuperimposer()))
	ptr = new PoseClusteringAffineSuperimposer();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~PoseClusteringAffineSuperimposer()))
	delete ptr;
END_SECTION

START_SECTION((static BaseSuperimposer* create()))
  BaseSuperimposer* base_ptr = 0;
	base_ptr = PoseClusteringAffineSuperimposer::create();
  TEST_NOT_EQUAL(base_ptr, base_nullPointer)
END_SECTION

START_SECTION((static const String getProductName()))
  PoseClusteringAffineSuperimposer pcat;
  
  TEST_EQUAL(pcat.getName() == "poseclustering_affine",true)
END_SECTION

START_SECTION((virtual void run(const std::vector< ConsensusMap > &maps, std::vector< TransformationDescription > &transformations)))
  std::vector<ConsensusMap> input(2);
  Feature feat1;
  Feature feat2;
  PositionType pos1(1,1);
  PositionType pos2(5,5);
  feat1.setPosition(pos1);
  feat1.setIntensity(100.0f);
  feat2.setPosition(pos2);
  feat2.setIntensity(100.0f);
  input[0].push_back(feat1);
  input[0].push_back(feat2);
  
  Feature feat3;
  Feature feat4;
  PositionType pos3(1.4,1.02);
  PositionType pos4(5.4,5.02);
  feat3.setPosition(pos3);
  feat3.setIntensity(100.0f);
  feat4.setPosition(pos4);
  feat4.setIntensity(100.0f);
  input[1].push_back(feat3);
  input[1].push_back(feat4);

  Param parameters;
  parameters.setValue(String("scaling_bucket_size"), 0.01);
  parameters.setValue(String("shift_bucket_size"), 0.1);

  // If hashing goes wrong, get debug output with the following:
  //  parameters.setValue(String("dump_buckets"),"pcast_buckets");
  //  parameters.setValue(String("dump_pairs"),"pcast_pairs");

  std::vector<TransformationDescription> transformations;  
  PoseClusteringAffineSuperimposer pcat;
  pcat.setParameters(parameters);

  // That's a precondition for run()!  Now even documented :-)
  input[0].updateRanges();
  input[1].updateRanges();

  pcat.run(input, transformations);

	TEST_EQUAL(transformations.size(), 1)
  TEST_STRING_EQUAL(transformations[0].getModelType(), "linear")
	transformations[0].getModelParameters(parameters);
	TEST_EQUAL(parameters.size(), 2)    
  TEST_REAL_SIMILAR(parameters.getValue("slope"), 1.0)
  TEST_REAL_SIMILAR(parameters.getValue("intercept"), -0.4)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



