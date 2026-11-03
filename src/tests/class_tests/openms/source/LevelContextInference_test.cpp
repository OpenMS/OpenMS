// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/OPENSWATH/LevelContextInference.h>

#include <algorithm>
#include <cmath>

using namespace OpenMS;

START_TEST(LevelContextInference, "$Id$")

TOLERANCE_ABSOLUTE(1e-10)
TOLERANCE_RELATIVE(1.0 + 1e-10)

START_SECTION(static std::vector<LevelContextResultRow> infer(const std::vector<LevelContextInputRow>&, const LevelContextInferenceConfig&))
{
  LevelContextInferenceConfig global_config;
  global_config.level = InferenceLevel::Peptide;
  global_config.context = InferenceContext::Global;

  {
    const std::vector<LevelContextResultRow> empty = LevelContextInference::infer({}, global_config);
    TEST_EQUAL(empty.empty(), true)
  }

  {
    const std::vector<LevelContextInputRow> global_rows =
    {
      {std::nullopt, "101", 101, false, 5.0, InferenceContext::Global},
      {std::nullopt, "102", 102, false, 4.0, InferenceContext::Global},
      {std::nullopt, "201", 201, true, 2.0, InferenceContext::Global},
      {std::nullopt, "202", 202, true, 1.0, InferenceContext::Global}
    };
    const std::vector<LevelContextResultRow> results = LevelContextInference::infer(global_rows, global_config);
    TEST_EQUAL(results.size(), global_rows.size())
    for (const auto& row : results)
    {
      TEST_EQUAL(row.run_id.has_value(), false)
      TEST_EQUAL(row.context, InferenceContext::Global)
      TEST_EQUAL(std::isfinite(row.pvalue), true)
      TEST_EQUAL(std::isfinite(row.qvalue), true)
      TEST_EQUAL(std::isfinite(row.pep), true)
    }
    auto best = std::find_if(results.begin(), results.end(), [](const auto& row) { return row.entity_id == 101; });
    auto second = std::find_if(results.begin(), results.end(), [](const auto& row) { return row.entity_id == 102; });
    TEST_EQUAL(best != results.end(), true)
    TEST_EQUAL(second != results.end(), true)
    if (best != results.end() && second != results.end())
    {
      TEST_EQUAL(best->qvalue <= second->qvalue, true)
    }
  }

  {
    LevelContextInferenceConfig run_specific_config;
    run_specific_config.level = InferenceLevel::Peptide;
    run_specific_config.context = InferenceContext::RunSpecific;

    const std::vector<LevelContextInputRow> run_rows =
    {
      {1, "1_101", 101, false, 5.0, InferenceContext::RunSpecific},
      {1, "1_102", 102, false, 4.0, InferenceContext::RunSpecific},
      {1, "1_201", 201, true, 2.0, InferenceContext::RunSpecific},
      {1, "1_202", 202, true, 1.0, InferenceContext::RunSpecific},
      {2, "2_301", 301, false, 5.0, InferenceContext::RunSpecific},
      {2, "2_302", 302, false, 4.0, InferenceContext::RunSpecific},
      {2, "2_401", 401, true, 2.0, InferenceContext::RunSpecific},
      {2, "2_402", 402, true, 1.0, InferenceContext::RunSpecific}
    };
    const std::vector<LevelContextResultRow> run_results = LevelContextInference::infer(run_rows, run_specific_config);
    TEST_EQUAL(run_results.size(), run_rows.size())

    auto run1_best = std::find_if(run_results.begin(), run_results.end(), [](const auto& row) { return row.run_id == 1 && row.entity_id == 101; });
    auto run2_best = std::find_if(run_results.begin(), run_results.end(), [](const auto& row) { return row.run_id == 2 && row.entity_id == 301; });
    auto run1_second = std::find_if(run_results.begin(), run_results.end(), [](const auto& row) { return row.run_id == 1 && row.entity_id == 102; });
    auto run2_second = std::find_if(run_results.begin(), run_results.end(), [](const auto& row) { return row.run_id == 2 && row.entity_id == 302; });
    TEST_EQUAL(run1_best != run_results.end(), true)
    TEST_EQUAL(run2_best != run_results.end(), true)
    TEST_EQUAL(run1_second != run_results.end(), true)
    TEST_EQUAL(run2_second != run_results.end(), true)
    if (run1_best != run_results.end() && run2_best != run_results.end())
    {
      TEST_REAL_SIMILAR(run1_best->pvalue, run2_best->pvalue)
      TEST_REAL_SIMILAR(run1_best->qvalue, run2_best->qvalue)
      TEST_REAL_SIMILAR(run1_best->pep, run2_best->pep)
    }
    if (run1_second != run_results.end() && run2_second != run_results.end())
    {
      TEST_REAL_SIMILAR(run1_second->qvalue, run2_second->qvalue)
    }
  }

  {
    LevelContextInferenceConfig experiment_config;
    experiment_config.level = InferenceLevel::Peptide;
    experiment_config.context = InferenceContext::ExperimentWide;

    const std::vector<LevelContextInputRow> experiment_rows =
    {
      {1, "1_101", 101, false, 5.0, InferenceContext::ExperimentWide},
      {1, "1_102", 102, false, 4.0, InferenceContext::ExperimentWide},
      {1, "1_201", 201, true, 2.0, InferenceContext::ExperimentWide},
      {1, "1_202", 202, true, 1.0, InferenceContext::ExperimentWide},
      {2, "2_301", 301, false, 5.0, InferenceContext::ExperimentWide},
      {2, "2_302", 302, false, 4.0, InferenceContext::ExperimentWide},
      {2, "2_401", 401, true, 2.0, InferenceContext::ExperimentWide},
      {2, "2_402", 402, true, 1.0, InferenceContext::ExperimentWide}
    };
    const std::vector<LevelContextResultRow> results = LevelContextInference::infer(experiment_rows, experiment_config);
    TEST_EQUAL(results.size(), experiment_rows.size())
    auto run1 = std::find_if(results.begin(), results.end(), [](const auto& row) { return row.run_id == 1 && row.entity_id == 101; });
    auto run2 = std::find_if(results.begin(), results.end(), [](const auto& row) { return row.run_id == 2 && row.entity_id == 301; });
    auto run1_second = std::find_if(results.begin(), results.end(), [](const auto& row) { return row.run_id == 1 && row.entity_id == 102; });
    auto run2_second = std::find_if(results.begin(), results.end(), [](const auto& row) { return row.run_id == 2 && row.entity_id == 302; });
    TEST_EQUAL(run1 != results.end(), true)
    TEST_EQUAL(run2 != results.end(), true)
    TEST_EQUAL(run1_second != results.end(), true)
    TEST_EQUAL(run2_second != results.end(), true)
    if (run1 != results.end() && run2 != results.end())
    {
      TEST_REAL_SIMILAR(run1->pvalue, run2->pvalue)
      TEST_REAL_SIMILAR(run1->qvalue, run2->qvalue)
      TEST_REAL_SIMILAR(run1->pep, run2->pep)
    }
    if (run1_second != results.end() && run2_second != results.end())
    {
      TEST_REAL_SIMILAR(run1_second->qvalue, run2_second->qvalue)
    }
  }

  {
    LevelContextInferenceConfig invalid_run_config;
    invalid_run_config.level = InferenceLevel::Peptide;
    invalid_run_config.context = InferenceContext::RunSpecific;
    TEST_EXCEPTION(Exception::Precondition,
      LevelContextInference::infer({{std::nullopt, "101", 101, false, 5.0, InferenceContext::RunSpecific}}, invalid_run_config))
  }

}
END_SECTION

END_TEST
