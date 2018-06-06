// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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

#include <OpenMS/CONCEPT/Factory.h>

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
