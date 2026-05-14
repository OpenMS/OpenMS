// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathGeneInference.h>
#include <OpenMS/ANALYSIS/OPENSWATH/LevelContextInference.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathProteinInference.h>
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPeptideInference.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cmath>

using namespace OpenMS;

namespace
{
  void copySharedInferenceFixture_(const String& filename)
  {
    File::remove(filename);
    TEST_TRUE(File::copy(OPENMS_GET_TEST_DATA_PATH("PyProphet_inference_test.osw"), filename))
  }
} // namespace

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

  {
    String tmp_osw;
    NEW_TMP_FILE(tmp_osw);
    copySharedInferenceFixture_(tmp_osw);

    OSWFile osw(tmp_osw);
    OpenSwathPeptideInference peptide_inference;
    OpenSwathProteinInference protein_inference;
    OpenSwathGeneInference gene_inference;

    const auto peptide_input = osw.readLevelContextData(InferenceLevel::Peptide, InferenceContext::Global);
    const std::vector<LevelContextResultRow> peptide_results = peptide_inference.infer(peptide_input, global_config);
    TEST_EQUAL(peptide_results.size(), 682)
    auto peptide_row = std::find_if(peptide_results.begin(), peptide_results.end(), [](const auto& row) { return row.entity_id == 2; });
    TEST_EQUAL(peptide_row != peptide_results.end(), true)
    if (peptide_row != peptide_results.end())
    {
      TEST_EQUAL(std::fabs(peptide_row->score - 5.53895807266235) < 1e-8, true)
      TEST_EQUAL(std::fabs(peptide_row->pvalue - 0.0029) < 5e-4, true)
      // q-values and local-FDR/PEP estimates depend on spline/KDE details and
      // vary slightly across platforms, so keep these checks coarse.
      TEST_EQUAL(std::isfinite(peptide_row->qvalue), true)
      TEST_EQUAL(std::isfinite(peptide_row->pep), true)
      TEST_EQUAL(peptide_row->qvalue >= 0.0 && peptide_row->qvalue <= 1.0, true)
      TEST_EQUAL(peptide_row->pep >= 0.0 && peptide_row->pep <= 1.0, true)
      TEST_EQUAL(peptide_row->qvalue <= 0.05, true)
      TEST_EQUAL(peptide_row->pep <= 0.05, true)
    }

    LevelContextInferenceConfig peptide_run_config = global_config;
    peptide_run_config.context = InferenceContext::RunSpecific;
    const auto peptide_run_input = osw.readLevelContextData(InferenceLevel::Peptide, InferenceContext::RunSpecific);
    const std::vector<LevelContextResultRow> peptide_run_results = peptide_inference.infer(peptide_run_input, peptide_run_config);
    TEST_EQUAL(peptide_run_results.size(), 682)
    TEST_EQUAL(std::count_if(peptide_run_results.begin(), peptide_run_results.end(),
                             [](const auto& row) { return row.run_id.has_value(); }),
               682)

    LevelContextInferenceConfig protein_global_config;
    protein_global_config.level = InferenceLevel::Protein;
    protein_global_config.context = InferenceContext::Global;
    const auto protein_input = osw.readLevelContextData(InferenceLevel::Protein, InferenceContext::Global);
    const std::vector<LevelContextResultRow> protein_results = protein_inference.infer(protein_input, protein_global_config);
    TEST_EQUAL(protein_results.size(), 32)
    auto protein_row = std::find_if(protein_results.begin(), protein_results.end(), [](const auto& row) { return row.entity_id == 0; });
    TEST_EQUAL(protein_row != protein_results.end(), true)
    if (protein_row != protein_results.end())
    {
      TEST_EQUAL(std::fabs(protein_row->score - 5.84014081954956) < 1e-8, true)
      TEST_EQUAL(std::fabs(protein_row->pvalue - 0.0625) < 5e-4, true)
      TEST_EQUAL(std::isfinite(protein_row->qvalue), true)
      TEST_EQUAL(std::isfinite(protein_row->pep), true)
      // Protein-level local FDR is estimated from only 32 rows in this fixture,
      // so the exact PEP can shift noticeably across platforms. Keep only basic
      // sanity checks here and leave tighter numeric checks to the score and
      // p-value, which are stable for this fixture.
      TEST_EQUAL(protein_row->qvalue >= 0.0 && protein_row->qvalue <= 1.0, true)
      TEST_EQUAL(protein_row->pep >= 0.0 && protein_row->pep <= 1.0, true)
      TEST_EQUAL(protein_row->qvalue <= 0.5, true)
    }

    LevelContextInferenceConfig gene_global_config;
    gene_global_config.level = InferenceLevel::Gene;
    gene_global_config.context = InferenceContext::Global;
    const auto gene_input = osw.readLevelContextData(InferenceLevel::Gene, InferenceContext::Global);
    const std::vector<LevelContextResultRow> gene_results = gene_inference.infer(gene_input, gene_global_config);
    TEST_EQUAL(gene_results.size(), 32)
    auto gene_row = std::find_if(gene_results.begin(), gene_results.end(), [](const auto& row) { return row.entity_id == 0; });
    TEST_EQUAL(gene_row != gene_results.end(), true)
    if (gene_row != gene_results.end())
    {
      TEST_EQUAL(std::fabs(gene_row->score - 5.84014081954956) < 1e-8, true)
    }
  }
}
END_SECTION

START_SECTION("OpenSwathPeptideInference / OpenSwathProteinInference / OpenSwathGeneInference wrappers")
{
  OpenSwathPeptideInference peptide_inference;
  OpenSwathProteinInference protein_inference;
  OpenSwathGeneInference gene_inference;

  LevelContextInferenceConfig peptide_config;
  peptide_config.level = InferenceLevel::Peptide;
  peptide_config.context = InferenceContext::Global;
  TEST_EQUAL(peptide_inference.infer({}, peptide_config).empty(), true)

  LevelContextInferenceConfig protein_config;
  protein_config.level = InferenceLevel::Protein;
  protein_config.context = InferenceContext::Global;
  TEST_EQUAL(protein_inference.infer({}, protein_config).empty(), true)

  LevelContextInferenceConfig gene_config;
  gene_config.level = InferenceLevel::Gene;
  gene_config.context = InferenceContext::Global;
  TEST_EQUAL(gene_inference.infer({}, gene_config).empty(), true)

  TEST_EXCEPTION(Exception::Precondition, peptide_inference.infer({}, protein_config))
  TEST_EXCEPTION(Exception::Precondition, protein_inference.infer({}, gene_config))
  TEST_EXCEPTION(Exception::Precondition, gene_inference.infer({}, peptide_config))
}
END_SECTION

END_TEST
