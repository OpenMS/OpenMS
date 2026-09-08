// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/MAPMATCHING/TransformationModelDefaults.h>
#include <OpenMS/APPLICATIONS/MapAlignerBase.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmTreeGuided.h>
#include <OpenMS/CONCEPT/ClassTest.h>

using namespace OpenMS;

START_TEST(TransformationModelDefaults, "$Id$")

START_SECTION((model selection and compatibility parameter trees))
  const std::vector<std::string> models = {"linear", "b_spline", "lowess", "interpolated"};
  for (const std::string selected : {"linear", "b_spline", "lowess", "interpolated", "none", "custom"})
  {
    const Param params = TransformationModelDefaults::getDefaults(selected);
    TEST_TRUE(params == MapAlignerBase::getModelDefaults(selected))
    TEST_EQUAL(params.getValue("type"), selected)
    TEST_EQUAL(params.getDescription("type"), "Type of model")
    auto expected = models;
    if (selected == "none" || selected == "custom")
    {
      expected.insert(expected.begin(), selected);
    }
    TEST_TRUE(params.getValidStrings("type") == expected)
    for (const auto& model : models)
    {
      TEST_EQUAL(params.getSectionDescription(model), "Parameters for '" + model + "' model")
    }
    TEST_EQUAL(params.getValue("b_spline:num_nodes"), 5)
    TEST_EQUAL(params.getValue("b_spline:boundary_condition"), 2)
    TEST_REAL_SIMILAR(static_cast<double>(params.getValue("lowess:span")), 2.0 / 3.0)
    TEST_EQUAL(params.getMinFloat("lowess:span"), 0.0)
    TEST_EQUAL(params.getMaxFloat("lowess:span"), 1.0)
    TEST_EQUAL(params.getValue("interpolated:interpolation_type"), "cspline")
  }
END_SECTION

START_SECTION((tree-guided alignment retains its model parameter subtree))
  MapAlignmentAlgorithmTreeGuided algorithm;
  const Param model = algorithm.getDefaults().copy("model:", true);
  TEST_TRUE(model == TransformationModelDefaults::getDefaults("b_spline"))
END_SECTION

END_TEST
