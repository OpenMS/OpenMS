// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
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

	// degenerate inputs (previously a division by zero): empty tool map -> recall 0; empty ground truth -> exception
	ConsensusMap empty_map;
	double out_empty = -1.0;
	maea.evaluate(empty_map, gt, 0.1, 0.1, 100, true, out_empty);
	TEST_REAL_SIMILAR(out_empty, 0.0)
	TEST_EXCEPTION(Exception::InvalidParameter, maea.evaluate(in, empty_map, 0.1, 0.1, 100, true, out_empty))
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
