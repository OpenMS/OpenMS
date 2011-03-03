// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
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

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/FORMAT/MzMLFile.h>

using namespace std;
using namespace OpenMS;


/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithmPoseClustering, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithmPoseClustering* ptr = 0;
START_SECTION((MapAlignmentAlgorithmPoseClustering()))
	ptr = new MapAlignmentAlgorithmPoseClustering();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~MapAlignmentAlgorithmPoseClustering()))
	delete ptr;
END_SECTION

START_SECTION((static MapAlignmentAlgorithm* create()))
	TEST_NOT_EQUAL(MapAlignmentAlgorithmPoseClustering::create(),0)
END_SECTION

START_SECTION((static String getProductName()))
	TEST_EQUAL(MapAlignmentAlgorithmPoseClustering::getProductName(), "pose_clustering")
END_SECTION

START_SECTION((virtual void setReference(Size reference_index=0, const String& reference_file="")))
{
	NOT_TESTABLE; // only some internal variables are set
}
END_SECTION

START_SECTION((virtual void alignPeakMaps(std::vector< MSExperiment<> > &, std::vector< TransformationDescription > &)))
{
  MzMLFile f;
  std::vector< MSExperiment<> > peak_maps(2);
  f.load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmPoseClustering_in1.mzML.gz"), peak_maps[0]);
  f.load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmPoseClustering_in2.mzML.gz"), peak_maps[1]);
  
  MapAlignmentAlgorithm* alignment = Factory<MapAlignmentAlgorithm>::create("pose_clustering");
  std::vector<TransformationDescription> transformations;
  // Trafo cannot be computed, due to too few datapoints
  // -- the ideal solution would be to fix the trafo estimation
  TEST_EXCEPTION(Exception::InvalidValue, alignment->alignPeakMaps(peak_maps,transformations));
}
END_SECTION

START_SECTION((virtual void alignFeatureMaps(std::vector< FeatureMap<> > &, std::vector< TransformationDescription > &)))
{
  // Tested extensively in TEST/TOPP.  See MapAligner_test.
  NOT_TESTABLE;
}
END_SECTION

START_SECTION((virtual void getDefaultModel(String& model_type, Param& params)))
{
	String model_type;
	Param params;
	MapAlignmentAlgorithmPoseClustering aligner;
	aligner.getDefaultModel(model_type, params);
	TEST_EQUAL(model_type, "linear");
	TEST_EQUAL(params.getValue("symmetric_regression"), "true");
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
