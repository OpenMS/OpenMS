// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include "OpenSwathGeneInference.h"
#include <OpenMS/ANALYSIS/OPENSWATH/LevelContextInference.h>
#include "OpenSwathProteinInference.h"
#include "OpenSwathPeptideInference.h"
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>
#include <cmath>

using namespace OpenMS;

namespace
{
  void copySharedInferenceFixture_(const std::string& filename)
  {
    File::remove(filename);
    TEST_TRUE(File::copy(OPENMS_GET_TEST_DATA_PATH("PyProphet_inference_test.osw"), filename))
  }
} // namespace

START_TEST(OpenSwathInference, "$Id$")

TOLERANCE_ABSOLUTE(1e-10)
TOLERANCE_RELATIVE(1.0 + 1e-10)

START_SECTION("Inference wrappers on the shared OSW fixture")
{
  LevelContextInferenceConfig global_config;
  global_config.level = InferenceLevel::Peptide;
  global_config.context = InferenceContext::Global;
  {
    std::string tmp_osw;
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
