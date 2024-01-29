// Copyright (c) 2002-present, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmPrecision.h>

#include <OpenMS/FORMAT/ConsensusXMLFile.h>

///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(MapAlignmentEvaluationAlgorithmPrecision, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MapAlignmentEvaluationAlgorithmPrecision* ptr = nullptr;
MapAlignmentEvaluationAlgorithmPrecision* nullPointer = nullptr;

START_SECTION((MapAlignmentEvaluationAlgorithmPrecision()))
	ptr = new MapAlignmentEvaluationAlgorithmPrecision();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~MapAlignmentEvaluationAlgorithmPrecision()))
	delete ptr;
END_SECTION

MapAlignmentEvaluationAlgorithm* base_nullPointer = nullptr;
START_SECTION((static MapAlignmentEvaluationAlgorithm* create()))
	MapAlignmentEvaluationAlgorithm* ptr2 = nullptr;
	ptr2 = MapAlignmentEvaluationAlgorithmPrecision::create();
  TEST_NOT_EQUAL(ptr2, base_nullPointer)
  delete ptr2;
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(MapAlignmentEvaluationAlgorithmPrecision::getProductName(),"precision")
END_SECTION

START_SECTION((virtual void evaluate(const ConsensusMap &consensus_map_in, const ConsensusMap &consensus_map_gt, const double &rt_dev, const double &mz_dev, const Peak2D::IntensityType &int_dev, const bool use_charge, double &out)))
	MapAlignmentEvaluationAlgorithmPrecision maea;
	ConsensusMap in;
	ConsensusMap gt;
	double out;

	ConsensusXMLFile consensus_xml_file_in;
	consensus_xml_file_in.load( OPENMS_GET_TEST_DATA_PATH("MapAlignmentEvaluationAlgorithm_in.consensusXML"), in );

	ConsensusXMLFile consensus_xml_file_gt;
	consensus_xml_file_gt.load( OPENMS_GET_TEST_DATA_PATH("MapAlignmentEvaluationAlgorithm_gt.consensusXML"), gt );

	maea.evaluate(in, gt, 0.1, 0.1, 100, true, out);

	TEST_REAL_SIMILAR(out, 0.757143)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

