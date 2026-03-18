// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Eva Lange, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelLinear.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace std;
using namespace OpenMS;

/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithmPoseClustering, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithmPoseClustering* ptr = nullptr;
MapAlignmentAlgorithmPoseClustering* nullPointer = nullptr;
START_SECTION((MapAlignmentAlgorithmPoseClustering()))
	ptr = new MapAlignmentAlgorithmPoseClustering();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~MapAlignmentAlgorithmPoseClustering()))
	delete ptr;
END_SECTION

START_SECTION((template <typename MapType> void setReference(const MapType& map)))
{
  NOT_TESTABLE // tested together with "align"
}
END_SECTION

START_SECTION((void align(const PeakMap& map, TransformationDescription& trafo)))
{
  MzMLFile f;
  std::vector<PeakMap > maps(2);
  f.load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmPoseClustering_in1.mzML.gz"), maps[0]);
  f.load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmPoseClustering_in2.mzML.gz"), maps[1]);

  MapAlignmentAlgorithmPoseClustering aligner;
  aligner.setReference(maps[0]);

  TransformationDescription trafo;
  aligner.align(maps[1], trafo);

  TEST_EQUAL(trafo.getModelType(), "linear");
  TEST_EQUAL(trafo.getDataPoints().size(), 307);
  
  // @TODO: can we get the slope/intercept without fitting a model again?
  TransformationModelLinear lm(trafo.getDataPoints(),
                               trafo.getModelParameters());
  double slope, intercept;
  String x_weight, y_weight;
  double x_datum_min, x_datum_max, y_datum_min, y_datum_max;
  lm.getParameters(slope, intercept, x_weight, y_weight, x_datum_min, x_datum_max, y_datum_min, y_datum_max);
  TEST_REAL_SIMILAR(slope, 1.01164);
  TEST_REAL_SIMILAR(intercept, -32.0912);
}
END_SECTION

START_SECTION((void align(const FeatureMap& map, TransformationDescription& trafo)))
{
  // Tested extensively in TEST/TOPP
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((void align(const ConsensusMap& map, TransformationDescription& trafo)))
{
  // Tested extensively in TEST/TOPP
  NOT_TESTABLE;
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
