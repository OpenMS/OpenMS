// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(MapAlignmentEvaluationAlgorithmRecall, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MapAlignmentEvaluationAlgorithmRecall* ptr = nullptr;
MapAlignmentEvaluationAlgorithmRecall* nullPointer = nullptr;

START_SECTION((MapAlignmentEvaluationAlgorithmRecall()))
	ptr = new MapAlignmentEvaluationAlgorithmRecall();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~MapAlignmentEvaluationAlgorithmRecall()))
	delete ptr;
END_SECTION

MapAlignmentEvaluationAlgorithm* base_nullPointer = nullptr;
START_SECTION((static MapAlignmentEvaluationAlgorithm* create()))
	MapAlignmentEvaluationAlgorithm* ptr2 = nullptr;
	ptr2 = MapAlignmentEvaluationAlgorithmRecall::create();
  TEST_NOT_EQUAL(ptr2, base_nullPointer)
  delete ptr2;
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(MapAlignmentEvaluationAlgorithmRecall::getProductName(),"recall")
END_SECTION

START_SECTION((virtual void evaluate(const ConsensusMap &consensus_map_in, const ConsensusMap &consensus_map_gt, const double &rt_dev, const double &mz_dev, const Peak2D::IntensityType &int_dev, const bool use_charge, double &out)))
	MapAlignmentEvaluationAlgorithmRecall maea;
	ConsensusMap in;
	ConsensusMap gt;
	double out;

	ConsensusXMLFile consensus_xml_file_in;
	consensus_xml_file_in.load( OPENMS_GET_TEST_DATA_PATH("MapAlignmentEvaluationAlgorithm_in.consensusXML"), in );

	ConsensusXMLFile consensus_xml_file_gt;
	consensus_xml_file_gt.load( OPENMS_GET_TEST_DATA_PATH("MapAlignmentEvaluationAlgorithm_gt.consensusXML"), gt );

	maea.evaluate(in, gt, 0.1, 0.1, 100, true, out);

	TEST_REAL_SIMILAR(out, 0.5)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
