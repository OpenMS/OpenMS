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
#include <OpenMS/KERNEL/StandardTypes.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/PoseClusteringShiftSuperimposer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

typedef DPosition <2> PositionType;

START_TEST(PoseClusteringShiftSuperimposer, "$Id: PoseClusteringShiftSuperimposer_test.C 5994 2009-09-21 22:17:43Z groepl $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PoseClusteringShiftSuperimposer* ptr = 0;
START_SECTION((PoseClusteringShiftSuperimposer()))
	ptr = new PoseClusteringShiftSuperimposer();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~PoseClusteringShiftSuperimposer()))
	delete ptr;
END_SECTION

START_SECTION((static BaseSuperimposer* create()))
  BaseSuperimposer* base_ptr = 0;
	base_ptr = PoseClusteringShiftSuperimposer::create();
	TEST_NOT_EQUAL(base_ptr, 0)
  delete (base_ptr);
END_SECTION

START_SECTION((static const String getProductName()))
  PoseClusteringShiftSuperimposer pcsi;
  
  TEST_EQUAL(pcsi.getName() == "poseclustering_shift",true)
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
  PositionType pos3(21.4,1.02);
  PositionType pos4(25.4,5.02);
  feat3.setPosition(pos3);
  feat3.setIntensity(100.0f);
  feat4.setPosition(pos4);
  feat4.setIntensity(100.0f);
  input[1].push_back(feat3);
  input[1].push_back(feat4);

  std::vector<TransformationDescription> transformations;
  PoseClusteringShiftSuperimposer pcat;
#if 0 // switch this on for debugging
	Param params;
  params.setValue("dump_buckets","tmp_PoseClusteringShiftSuperimposer_buckets");
  params.setValue("dump_pairs","tmp_PoseClusteringShiftSuperimposer_pairs");
  pcat.setParameters(params);
#endif
	pcat.run(input,transformations);
  
	TEST_EQUAL(transformations.size(),1)
  TEST_STRING_EQUAL(transformations[0].getName(),"linear")
	TEST_EQUAL(transformations[0].getParameters().size(),2)
  TEST_REAL_SIMILAR(transformations[0].getParameters().getValue("slope"),1.0)
  TEST_REAL_SIMILAR(transformations[0].getParameters().getValue("intercept"),-20.4)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



