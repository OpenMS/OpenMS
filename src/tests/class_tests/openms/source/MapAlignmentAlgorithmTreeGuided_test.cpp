// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julia Thueringer $
// $Authors: Julia Thueringer $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmTreeGuided.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>

#include <iostream>

using namespace std;
using namespace OpenMS;

/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithmTreeGuided, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithmTreeGuided* ptr = nullptr;
MapAlignmentAlgorithmTreeGuided* nullPointer = nullptr;
START_SECTION((MapAlignmentAlgorithmTreeGuided()))
  ptr = new MapAlignmentAlgorithmTreeGuided();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION


START_SECTION((virtual ~MapAlignmentAlgorithmTreeGuided()))
  delete ptr;
END_SECTION

vector<FeatureMap> maps(3);
FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmTreeGuided_test_in0.featureXML"), maps[0]);
FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmTreeGuided_test_in1.featureXML"), maps[1]);
FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("MapAlignmentAlgorithmTreeGuided_test_in2.featureXML"), maps[2]);
// copy maps for computeTrafosByOriginalRT and computeTransformedFeatureMaps
vector<FeatureMap> maps_orig = maps;

MapAlignmentAlgorithmTreeGuided aligner;
aligner.setLogType(ProgressLogger::CMD);
Param params = aligner.getParameters();
aligner.setParameters(params);

vector<BinaryTreeNode> result_tree;
vector<vector<double>> maps_ranges(3);

FeatureMap map_transformed;
vector<Size> trafo_order;

vector<TransformationDescription> trafos(3);

START_SECTION((static void buildTree(std::vector<FeatureMap>& feature_maps, std::vector<BinaryTreeNode>& tree, std::vector<std::vector<double>>& maps_ranges)))
{
  // test of protected nested class PeptideIdentificationsPearsonDistance_ that functions as comparator for ClusterHierarchical with AverageLinkage:
  // input map in0 and in2 are nearly identical with in2 having larger RT range, in1 has largest rt range and differs in identifications
  vector<BinaryTreeNode> test_tree;
  test_tree.emplace_back(0, 2, 1.84834e-04);
  test_tree.emplace_back(0, 1, 0.505752);
  OpenMS::MapAlignmentAlgorithmTreeGuided::buildTree(maps, result_tree, maps_ranges);

  TEST_EQUAL(result_tree.size(), test_tree.size());
  for (Size i = 0; i < result_tree.size(); ++i)
  {
    TEST_EQUAL(test_tree[i].left_child, result_tree[i].left_child);
    TEST_EQUAL(test_tree[i].right_child, result_tree[i].right_child);
    TEST_REAL_SIMILAR(test_tree[i].distance, result_tree[i].distance);
  }

  TEST_EQUAL(maps_ranges.size(), 3);
  // peptide identification counts for indirect test of protected methods extractSeqAndRt_ and addPeptideSequences_
  TEST_EQUAL(maps_ranges[0].size(), 6);
  TEST_EQUAL(maps_ranges[1].size(), 5);
  TEST_EQUAL(maps_ranges[2].size(), 5);
}
END_SECTION

START_SECTION((void treeGuidedAlignment(const std::vector<BinaryTreeNode>& tree, std::vector<FeatureMap>& feature_maps_transformed,
        std::vector<std::vector<double>>& maps_ranges, FeatureMap& map_transformed, std::vector<Size>& trafo_order)))
{
  aligner.treeGuidedAlignment(result_tree, maps, maps_ranges, map_transformed, trafo_order);

  TEST_EQUAL(map_transformed.size(), 15); // contains 3*5 feature from input maps
  // map_transformed contains all input map feature in order of trafo_order
  // trafo order should be: (1, (2, 0)), because cluster with larger rt is reference in alignment and other cluster is attached to it
  TEST_EQUAL(trafo_order[0], 1);
  TEST_EQUAL(map_transformed[0].getUniqueId(), 20);
  TEST_EQUAL(trafo_order[2], 0);
  TEST_EQUAL(map_transformed.back().getUniqueId(), 14);

  // order of aligned features should correspond to trafo_order
  // check indirectly with the existence of meta value "original_RT"
  // RTs of in1 (first 5 features) should be unchanged (no meta value) because map is last cluster and a reference
  for (Size i = 0; i < 5; ++i)
  {
    TEST_EQUAL(map_transformed[i].metaValueExists("original_RT"), false);
  }
  // feature RTs of maps 0 and 2 should be corrected -> meta value exists
  for (Size i = 5; i < 15; ++i)
  {
    TEST_EQUAL(map_transformed[i].metaValueExists("original_RT"), true);
  }
}
END_SECTION

START_SECTION((void computeTrafosByOriginalRT(std::vector<FeatureMap>& feature_maps, FeatureMap& map_transformed,
        std::vector<TransformationDescription>& transformations, const std::vector<Size>& trafo_order)))
{
  aligner.computeTrafosByOriginalRT(maps_orig, map_transformed, trafos, trafo_order);

  TEST_EQUAL(trafos.size(), 3);

  for (Size i = 0; i < maps.size(); ++i)
  {
    // first rt in trafo should be the same as in original map
    Size j = 0;
    for (auto feature_it = maps_orig[i].begin(); feature_it < maps_orig[i].end(); ++feature_it)
    {
      TEST_REAL_SIMILAR(trafos[i].getDataPoints()[j].first, feature_it->getRT());
      ++j;
    }
  }
}
END_SECTION

START_SECTION((static void computeTransformedFeatureMaps(std::vector<FeatureMap>& feature_maps, const std::vector<TransformationDescription>& transformations)))
{
  OpenMS::MapAlignmentAlgorithmTreeGuided::computeTransformedFeatureMaps(maps_orig, trafos);

  // check storing of original RTs:
  for (auto& map : maps_orig)
  {
    for (auto feat_it = map.begin(); feat_it < map.end(); ++feat_it)
    {
      TEST_EQUAL(feat_it->metaValueExists("original_RT"), true);
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
