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
  // alignment without reference, to a consensus of the input maps:
  Param consensus_params = params;
  consensus_params.setValue("auto_reference", "consensus");
  aligner.setParameters(consensus_params);
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
  aligner.setParameters(params); // back to the default ("best_run")

  // alignment with internal reference:
  transforms.clear();
  aligner.align(peptides, transforms, 0);

  TEST_EQUAL(transforms.size(), 2);
  TEST_EQUAL(transforms[0].getModelType(), "identity");
  TEST_EQUAL(transforms[1].getDataPoints().size(), 10);

  map<std::string, double> rts_second; // RTs in the second map, per sequence
  for (Size i = 0; i < transforms[1].getDataPoints().size(); ++i)
  {
    // RT transform should map to RT scale of the reference:
    TEST_REAL_SIMILAR(transforms[1].getDataPoints()[i].second,
                      reference_rts[i]);
    rts_second[transforms[1].getDataPoints()[i].note] =
      transforms[1].getDataPoints()[i].first;
  }

  // alignment without reference, to the map that shares the most sequences
  // with every other map - with two maps a tie, and both have ten sequences,
  // so the first map is used:
  transforms.clear();
  aligner.align(peptides, transforms);

  TEST_EQUAL(transforms.size(), 2);
  TEST_EQUAL(transforms[0].getModelType(), "identity");
  TEST_EQUAL(transforms[1].getDataPoints().size(), 10);
  for (Size i = 0; i < transforms[1].getDataPoints().size(); ++i)
  {
    TEST_REAL_SIMILAR(transforms[1].getDataPoints()[i].second,
                      reference_rts[i]);
  }

  // with one ID less in the first map, the second one (more IDs) wins the tie:
  vector<PeptideIdentificationList> fewer_in_first = peptides;
  fewer_in_first[0].erase(fewer_in_first[0].begin());
  transforms.clear();
  aligner.align(fewer_in_first, transforms);

  TEST_EQUAL(transforms.size(), 2);
  TEST_EQUAL(transforms[0].getDataPoints().size(), 9);
  TEST_EQUAL(transforms[1].getModelType(), "identity");
  for (const auto& point : transforms[0].getDataPoints())
  {
    // RT transform should map to RT scale of the second map:
    TEST_REAL_SIMILAR(point.second, rts_second[point.note]);
  }

  // the reference picked in one call must not carry over into the next one
  // (it would be used as an external reference, so no map would get the
  // identity transformation):
  transforms.clear();
  aligner.align(peptides, transforms);

  TEST_EQUAL(transforms.size(), 2);
  TEST_EQUAL(transforms[0].getModelType(), "identity");
  TEST_EQUAL(transforms[1].getDataPoints().size(), 10);

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


START_SECTION([EXTRA] the automatic reference shares IDs with every other input)
{
  const std::string residues = "ACDEFGHIKLMNPQRSTVWY";
  auto sequence = [&residues](Size i)
  {
    return "PEPTIDE" + std::string(1, residues[i / 20]) + residues[i % 20];
  };
  auto add_id = [](PeptideIdentificationList& run, const std::string& seq, double rt)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(seq));
    hit.setScore(1.0);
    PeptideIdentification pep;
    pep.setRT(rt);
    pep.setHits({hit});
    run.push_back(pep);
  };

  // A has the most IDs (20 shared with B, 40 only its own), but shares none
  // with C; B shares 20 with each of the others, so it becomes the reference:
  vector<PeptideIdentificationList> runs(3);
  for (Size i = 0; i < 20; ++i)
  {
    double rt = 100.0 + 30.0 * i;
    add_id(runs[0], sequence(i), rt);
    add_id(runs[1], sequence(i), rt + 10.0);
    add_id(runs[1], sequence(20 + i), rt + 10.0);
    add_id(runs[2], sequence(20 + i), rt + 50.0);
  }
  for (Size i = 0; i < 40; ++i)
  {
    add_id(runs[0], sequence(40 + i), 100.0 + 15.0 * i);
  }

  MapAlignmentAlgorithmIdentification auto_aligner;
  vector<TransformationDescription> transforms;
  auto_aligner.align(runs, transforms);

  TEST_EQUAL(transforms.size(), 3);
  TEST_EQUAL(transforms[0].getDataPoints().size(), 20);
  TEST_EQUAL(transforms[1].getModelType(), "identity");
  TEST_EQUAL(transforms[2].getDataPoints().size(), 20);
  for (TransformationDescription& trafo : transforms)
  {
    trafo.fitModel("b_spline"); // throws if a run has too few data points
  }

  // each pair of runs shares only one sequence, so no run can be the
  // reference for both others - the runs are aligned to a consensus instead:
  vector<PeptideIdentificationList> sparse(3);
  add_id(sparse[0], sequence(0), 100.0);
  add_id(sparse[1], sequence(0), 110.0);
  add_id(sparse[1], sequence(1), 200.0);
  add_id(sparse[2], sequence(1), 250.0);
  add_id(sparse[2], sequence(2), 300.0);
  add_id(sparse[0], sequence(2), 290.0);

  transforms.clear();
  auto_aligner.align(sparse, transforms);

  TEST_EQUAL(transforms.size(), 3);
  for (const TransformationDescription& trafo : transforms)
  {
    TEST_EQUAL(trafo.getDataPoints().size(), 2);
  }
}
END_SECTION


START_SECTION([EXTRA] fallback to a consensus if the automatic reference provides too few alignment points)
{
  const std::string residues = "ACDEFGHIKLMNPQRSTVWY";
  auto sequence = [&residues](Size i)
  {
    return "PEPTIDE" + std::string(1, residues[i / 20]) + residues[i % 20];
  };
  auto add_id = [](PeptideIdentificationList& run, const std::string& seq, double rt)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(seq));
    hit.setScore(1.0);
    PeptideIdentification pep;
    pep.setRT(rt);
    pep.setHits({hit});
    run.push_back(pep);
  };

  // four runs, each pair shares four sequences that no other run has: any run
  // as reference gives the others only four points each, while a consensus
  // gives every run all of its twelve sequences:
  vector<PeptideIdentificationList> runs(4);
  Size pair_index = 0;
  for (Size i = 0; i < 4; ++i)
  {
    for (Size j = i + 1; j < 4; ++j, ++pair_index)
    {
      for (Size k = 0; k < 4; ++k)
      {
        Size s = 4 * pair_index + k;
        add_id(runs[i], sequence(s), 100.0 + 30.0 * s + 10.0 * i);
        add_id(runs[j], sequence(s), 100.0 + 30.0 * s + 10.0 * j);
      }
    }
  }

  MapAlignmentAlgorithmIdentification fallback_aligner;
  vector<TransformationDescription> transforms;
  fallback_aligner.align(runs, transforms);

  TEST_EQUAL(transforms.size(), 4);
  for (const TransformationDescription& trafo : transforms)
  {
    TEST_EQUAL(trafo.getDataPoints().size(), 12);
  }

  // four points are enough if the threshold says so - then the first run
  // (tie) stays the reference:
  Param fallback_params = fallback_aligner.getParameters();
  fallback_params.setValue("auto_reference_min_points", 4);
  fallback_aligner.setParameters(fallback_params);
  transforms.clear();
  fallback_aligner.align(runs, transforms);

  TEST_EQUAL(transforms.size(), 4);
  TEST_EQUAL(transforms[0].getModelType(), "identity");
  for (Size i = 1; i < 4; ++i)
  {
    TEST_EQUAL(transforms[i].getDataPoints().size(), 4);
  }

  // a run with few IDs overall gets no more points from a consensus, so it
  // does not pull the other runs away from the reference:
  vector<PeptideIdentificationList> with_sparse(3);
  for (Size s = 0; s < 33; ++s)
  {
    add_id(with_sparse[0], sequence(s), 100.0 + 30.0 * s);
    add_id(with_sparse[1], sequence(s), 110.0 + 30.0 * s);
    if (s >= 30) add_id(with_sparse[2], sequence(s), 150.0 + 30.0 * s);
  }
  MapAlignmentAlgorithmIdentification sparse_aligner;
  transforms.clear();
  sparse_aligner.align(with_sparse, transforms);

  TEST_EQUAL(transforms.size(), 3);
  TEST_EQUAL(transforms[0].getModelType(), "identity");
  TEST_EQUAL(transforms[1].getDataPoints().size(), 33);
  TEST_EQUAL(transforms[2].getDataPoints().size(), 3);
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
// 	map<std::string, DoubleList> seq_to_list;
// 	map<std::string, double> seq_to_value;
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
