// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>
#include <OpenMS/FORMAT/IdXMLFile.h>

#include <iostream>

using namespace std;
using namespace OpenMS;

/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithmIdentification, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithmIdentification* ptr = nullptr;
MapAlignmentAlgorithmIdentification* nullPointer = nullptr;
START_SECTION((MapAlignmentAlgorithmIdentification()))
	ptr = new MapAlignmentAlgorithmIdentification();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION


START_SECTION((virtual ~MapAlignmentAlgorithmIdentification()))
	delete ptr;
END_SECTION

vector<PeptideIdentificationList > peptides(2);
vector<ProteinIdentification> proteins;
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmIdentification_test_1.idXML"),	proteins, peptides[0]);
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmIdentification_test_2.idXML"),	proteins, peptides[1]);

MapAlignmentAlgorithmIdentification aligner;
aligner.setLogType(ProgressLogger::CMD);
Param params = aligner.getParameters();
params.setValue("peptide_score_threshold", 0.0);
aligner.setParameters(params);
vector<double> reference_rts; // needed later

START_SECTION((template <typename DataType> void align(std::vector<DataType>& data, std::vector<TransformationDescription>& transformations, Int reference_index = -1)))
{
  // alignment without reference:
	vector<TransformationDescription> transforms;
	aligner.align(peptides, transforms);

  TEST_EQUAL(transforms.size(), 2);
  TEST_EQUAL(transforms[0].getDataPoints().size(), 10);
  TEST_EQUAL(transforms[1].getDataPoints().size(), 10);

  reference_rts.reserve(10);
  for (Size i = 0; i < transforms[0].getDataPoints().size(); ++i)
  {
    // both RT transforms should map to a common RT scale:
    TEST_REAL_SIMILAR(transforms[0].getDataPoints()[i].second,
                      transforms[1].getDataPoints()[i].second);
    reference_rts.push_back(transforms[0].getDataPoints()[i].first);
  }

  // alignment with internal reference:
  transforms.clear();
  aligner.align(peptides, transforms, 0);

  TEST_EQUAL(transforms.size(), 2);
  TEST_EQUAL(transforms[0].getModelType(), "identity");
  TEST_EQUAL(transforms[1].getDataPoints().size(), 10);

  for (Size i = 0; i < transforms[1].getDataPoints().size(); ++i)
  {
    // RT transform should map to RT scale of the reference:
    TEST_REAL_SIMILAR(transforms[1].getDataPoints()[i].second,
                      reference_rts[i]);
  }

  // algorithm works the same way for other input data types -> no extra tests
}
END_SECTION


START_SECTION([EXTRA] repeated align() with internal reference does not leak stale reference state)
{
  // Regression: checkParameters_ used to inspect reference_ left over from a
  // previous align() call when the user asked for an internal reference, so
  // the effective run count was inflated by 1 on every subsequent call. With
  // min_run_occur larger than the actual number of input maps, the cap that
  // normally clamps it to the run count no longer fired, and the per-sequence
  // filter in computeTransformations_ dropped every peptide -> identity
  // transforms. Two consecutive calls with the same inputs must therefore
  // produce identical, non-empty alignments.
  MapAlignmentAlgorithmIdentification repeat_aligner;
  repeat_aligner.setLogType(ProgressLogger::CMD);
  Param repeat_params = repeat_aligner.getParameters();
  // Force the cap path: min_run_occur (3) > data.size() (2).
  repeat_params.setValue("min_run_occur", 3);
  repeat_aligner.setParameters(repeat_params);

  vector<TransformationDescription> first_transforms;
  repeat_aligner.align(peptides, first_transforms, 0);
  vector<TransformationDescription> second_transforms;
  repeat_aligner.align(peptides, second_transforms, 0);

  TEST_EQUAL(first_transforms.size(), 2);
  TEST_EQUAL(second_transforms.size(), 2);
  TEST_EQUAL(first_transforms[0].getModelType(), "identity"); // reference map
  TEST_EQUAL(second_transforms[0].getModelType(), "identity");
  const auto& first_points = first_transforms[1].getDataPoints();
  const auto& second_points = second_transforms[1].getDataPoints();
  TEST_NOT_EQUAL(first_points.size(), 0);
  TEST_EQUAL(second_points.size(), first_points.size());
  for (Size i = 0; i < first_points.size(); ++i)
  {
    TEST_REAL_SIMILAR(second_points[i].first, first_points[i].first);
    TEST_REAL_SIMILAR(second_points[i].second, first_points[i].second);
  }
  for (Size i = 0; i < first_transforms[1].getDataPoints().size(); ++i)
  {
    TEST_REAL_SIMILAR(second_transforms[1].getDataPoints()[i].first,
                      first_transforms[1].getDataPoints()[i].first);
    TEST_REAL_SIMILAR(second_transforms[1].getDataPoints()[i].second,
                      first_transforms[1].getDataPoints()[i].second);
  }
}
END_SECTION


START_SECTION((template <typename DataType> void setReference(DataType& data)))
{
  // alignment with external reference:
  aligner.setReference(peptides[0]);
  peptides.erase(peptides.begin());

  vector<TransformationDescription> transforms;
  aligner.align(peptides, transforms);

  TEST_EQUAL(transforms.size(), 1);
  TEST_EQUAL(transforms[0].getDataPoints().size(), 10);

  for (Size i = 0; i < transforms[0].getDataPoints().size(); ++i)
  {
    // RT transform should map to RT scale of the reference:
    TEST_REAL_SIMILAR(transforms[0].getDataPoints()[i].second,
                      reference_rts[i]);
  }
}
END_SECTION


// can't test protected methods...

// START_SECTION((void computeMedians_(SeqToList&, SeqToValue&, bool)))
// {
// 	map<String, DoubleList> seq_to_list;
// 	map<String, double> seq_to_value;
// 	seq_to_list["ABC"] << -1.0 << 2.5 << 0.5 << -3.5;
// 	seq_to_list["DEF"] << 1.5 << -2.5 << -1;
// 	computeMedians_(seq_to_list, seq_to_value, false);
// 	TEST_EQUAL(seq_to_value.size(), 2);
// 	TEST_EQUAL(seq_to_value["ABC"], -0.25);
// 	TEST_EQUAL(seq_to_value["DEF"], -1);
// 	seq_to_value.clear();
// 	computeMedians_(seq_to_list, seq_to_value, true); // should be sorted now
// 	TEST_EQUAL(seq_to_value.size(), 2);
// 	TEST_EQUAL(seq_to_value["ABC"], -0.25);
// 	TEST_EQUAL(seq_to_value["DEF"], -1);
// }
// END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
