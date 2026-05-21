// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathPeptidoformInference.h>
#include <OpenMS/FORMAT/OSWFile.h>
#include <OpenMS/SYSTEM/File.h>

#include <algorithm>

using namespace OpenMS;

namespace
{
  void copySharedInferenceFixture_(const String& filename)
  {
    File::remove(filename);
    TEST_TRUE(File::copy(OPENMS_GET_TEST_DATA_PATH("PyProphet_inference_test.osw"), filename))
  }
} // namespace

START_TEST(OpenSwathPeptidoformInference, "$Id$")

TOLERANCE_ABSOLUTE(1e-10)
TOLERANCE_RELATIVE(1.0 + 1e-10)

START_SECTION(static std::vector<double> computeModelFDR(const std::vector<double>& pep_values))
{
  const std::vector<double> pep_values = {0.2, 0.1, 0.1, 0.4};
  const std::vector<double> qvalues = OpenSwathPeptidoformInference::computeModelFDR(pep_values);
  TEST_EQUAL(qvalues.size(), pep_values.size())
  TEST_REAL_SIMILAR(qvalues[0], 0.4 / 3.0)
  TEST_REAL_SIMILAR(qvalues[1], 0.1)
  TEST_REAL_SIMILAR(qvalues[2], 0.1)
  TEST_REAL_SIMILAR(qvalues[3], 0.2)
}
END_SECTION

START_SECTION(static std::vector<PosteriorRow> applyBM(const std::vector<BayesianModelRow>& rows))
{
  const std::vector<OpenSwathPeptidoformInference::BayesianModelRow> rows =
  {
    {1, 0, 0.3, 0.5},
    {1, 1, 0.7, 0.5},
    {1, 1, 0.8, 0.5}
  };
  const std::vector<OpenSwathPeptidoformInference::PosteriorRow> posteriors = OpenSwathPeptidoformInference::applyBM(rows);
  TEST_EQUAL(posteriors.size(), 2)
  TEST_EQUAL(posteriors[0].feature_id, 1)
  TEST_EQUAL(posteriors[1].feature_id, 1)
  TEST_REAL_SIMILAR(posteriors[0].posterior + posteriors[1].posterior, 1.0)

  auto hyp0 = std::find_if(posteriors.begin(), posteriors.end(), [](const auto& row) { return row.hypothesis == 0; });
  auto hyp1 = std::find_if(posteriors.begin(), posteriors.end(), [](const auto& row) { return row.hypothesis == 1; });
  TEST_EQUAL(hyp0 != posteriors.end(), true)
  TEST_EQUAL(hyp1 != posteriors.end(), true)
  if (hyp0 != posteriors.end())
  {
    TEST_REAL_SIMILAR(hyp0->posterior, 0.15 / 0.325)
  }
  if (hyp1 != posteriors.end())
  {
    TEST_REAL_SIMILAR(hyp1->posterior, 0.175 / 0.325)
  }

  const std::vector<OpenSwathPeptidoformInference::BayesianModelRow> zero_rows =
  {
    {2, 0, 0.5, 0.0},
    {2, 1, 0.5, 0.0}
  };
  const std::vector<OpenSwathPeptidoformInference::PosteriorRow> zero_posteriors = OpenSwathPeptidoformInference::applyBM(zero_rows);
  TEST_EQUAL(zero_posteriors.size(), 2)
  TEST_REAL_SIMILAR(zero_posteriors[0].posterior, 0.0)
  TEST_REAL_SIMILAR(zero_posteriors[1].posterior, 0.0)
}
END_SECTION

START_SECTION(std::vector<IPFPrecursorProbabilityRow> precursorInference(const std::vector<IPFPrecursorRow>&, const PeptidoformInferenceConfig&) const)
{
  OpenSwathPeptidoformInference inference;

  {
    PeptidoformInferenceConfig config;
    config.ipf_ms1_scoring = true;
    config.ipf_ms2_scoring = false;
    config.ipf_max_precursor_peakgroup_pep = 1.0;
    const std::vector<IPFPrecursorProbabilityRow> results = inference.precursorInference({{1, 0.2, 0.1, std::nullopt}}, config);
    TEST_EQUAL(results.size(), 1)
    if (!results.empty())
    {
      TEST_REAL_SIMILAR(results[0].precursor_peakgroup_pep, 1.0 - (0.72 / 0.74))
    }
  }

  {
    PeptidoformInferenceConfig config;
    config.ipf_ms1_scoring = false;
    config.ipf_ms2_scoring = true;
    config.ipf_max_precursor_peakgroup_pep = 1.0;
    const std::vector<IPFPrecursorProbabilityRow> results = inference.precursorInference({{1, 0.2, std::nullopt, 0.25}}, config);
    TEST_EQUAL(results.size(), 1)
    if (!results.empty())
    {
      TEST_REAL_SIMILAR(results[0].precursor_peakgroup_pep, 1.0 - (0.6 / 0.65))
    }
  }

  {
    PeptidoformInferenceConfig config;
    config.ipf_ms1_scoring = true;
    config.ipf_ms2_scoring = true;
    config.ipf_max_precursor_peakgroup_pep = 1.0;
    const std::vector<IPFPrecursorProbabilityRow> results = inference.precursorInference({{1, 0.3, 0.1, 0.2}}, config);
    TEST_EQUAL(results.size(), 1)
    if (!results.empty())
    {
      TEST_REAL_SIMILAR(results[0].precursor_peakgroup_pep, 1.0 - (0.504 / 0.51))
    }
  }

  {
    PeptidoformInferenceConfig config;
    config.ipf_ms1_scoring = false;
    config.ipf_ms2_scoring = false;
    config.ipf_max_precursor_peakgroup_pep = 0.3;
    const std::vector<IPFPrecursorProbabilityRow> results = inference.precursorInference(
      {{1, 0.25, std::nullopt, std::nullopt}, {2, 0.35, std::nullopt, std::nullopt}},
      config
    );
    TEST_EQUAL(results.size(), 1)
    if (!results.empty())
    {
      TEST_EQUAL(results[0].feature_id, 1)
      TEST_REAL_SIMILAR(results[0].precursor_peakgroup_pep, 0.25)
    }
  }

  {
    PeptidoformInferenceConfig config;
    config.ipf_ms1_scoring = true;
    config.ipf_ms2_scoring = false;
    config.ipf_max_precursor_pep = 0.05;
    config.ipf_max_precursor_peakgroup_pep = 1.0;
    const std::vector<IPFPrecursorProbabilityRow> results = inference.precursorInference({{1, 0.2, 0.1, std::nullopt}}, config);
    // MS1 evidence above the configured precursor threshold is treated as absent.
    // PyProphet then falls back to the synthetic true/false precursor rows, which
    // yield a precursor PEP of 1.0 here and are filtered by the strict "<" cutoff.
    TEST_EQUAL(results.empty(), true)
  }
}
END_SECTION

START_SECTION(std::vector<IPFResultRow> infer(const std::vector<IPFPrecursorRow>&, const std::vector<IPFTransitionRow>&, const std::vector<IPFAlignmentRow>&, const PeptidoformInferenceConfig&) const)
{
  OpenSwathPeptidoformInference inference;

  {
    PeptidoformInferenceConfig config;
    config.ipf_ms1_scoring = false;
    config.ipf_ms2_scoring = false;
    config.ipf_max_precursor_peakgroup_pep = 1.0;

    const std::vector<IPFResultRow> results = inference.infer(
      {{1, 0.2, std::nullopt, std::nullopt}},
      {
        {1, 11, 101, 0.1, 1, 2, std::nullopt},
        {1, 11, 102, 0.1, 0, 2, std::nullopt},
        {1, 11, -1, 0.1, 0, 2, std::nullopt}
      },
      {},
      config
    );

    TEST_EQUAL(results.size(), 2)
    auto pep101 = std::find_if(results.begin(), results.end(), [](const auto& row) { return row.peptide_id == 101; });
    auto pep102 = std::find_if(results.begin(), results.end(), [](const auto& row) { return row.peptide_id == 102; });
    TEST_EQUAL(pep101 != results.end(), true)
    TEST_EQUAL(pep102 != results.end(), true)
    if (pep101 != results.end())
    {
      TEST_REAL_SIMILAR(pep101->precursor_peakgroup_pep, 0.2)
      TEST_REAL_SIMILAR(pep101->pep, 1.0 - (0.36 / 0.42))
      TEST_REAL_SIMILAR(pep101->qvalue, 0.14285714285714285)
    }
    if (pep102 != results.end())
    {
      TEST_REAL_SIMILAR(pep102->pep, 1.0 - (0.04 / 0.42))
      TEST_REAL_SIMILAR(pep102->qvalue, 0.5238095238095238)
    }
  }

  {
    PeptidoformInferenceConfig grouped_config;
    grouped_config.ipf_ms1_scoring = false;
    grouped_config.ipf_ms2_scoring = false;
    grouped_config.ipf_max_precursor_peakgroup_pep = 1.0;
    grouped_config.ipf_grouped_fdr = true;

    const std::vector<IPFTransitionRow> transition_rows =
    {
      {1, 11, 101, 0.1, 1, 2, std::nullopt},
      {1, 11, -1, 0.1, 0, 2, std::nullopt},
      {2, 21, 201, 0.1, 1, 3, std::nullopt},
      {2, 21, -1, 0.1, 0, 3, std::nullopt}
    };
    const std::vector<IPFPrecursorRow> precursor_rows =
    {
      {1, 0.2, std::nullopt, std::nullopt},
      {2, 0.2, std::nullopt, std::nullopt}
    };

    const std::vector<IPFResultRow> grouped = inference.infer(precursor_rows, transition_rows, {}, grouped_config);
    TEST_EQUAL(grouped.size(), 2)
    if (grouped.size() == 2)
    {
      TEST_REAL_SIMILAR(grouped[0].qvalue, grouped[0].pep)
      TEST_REAL_SIMILAR(grouped[1].qvalue, grouped[1].pep)
    }

    PeptidoformInferenceConfig nongrouped_config = grouped_config;
    nongrouped_config.ipf_grouped_fdr = false;
    const std::vector<IPFResultRow> nongrouped = inference.infer(precursor_rows, transition_rows, {}, nongrouped_config);
    TEST_EQUAL(nongrouped.size(), 2)
    if (nongrouped.size() == 2 && grouped.size() == 2)
    {
      TEST_REAL_SIMILAR(nongrouped[0].qvalue, nongrouped[0].pep)
      TEST_REAL_SIMILAR(nongrouped[1].qvalue, (grouped[0].pep + grouped[1].pep) / 2.0)
      TEST_EQUAL(nongrouped[1].qvalue < grouped[1].qvalue, true)
    }
  }

  {
    PeptidoformInferenceConfig base_config;
    base_config.ipf_ms1_scoring = false;
    base_config.ipf_ms2_scoring = false;
    base_config.ipf_max_precursor_peakgroup_pep = 1.0;

    const std::vector<IPFPrecursorRow> precursor_rows =
    {
      {1, 0.2, std::nullopt, std::nullopt},
      {2, 0.2, std::nullopt, std::nullopt}
    };
    const std::vector<IPFTransitionRow> transition_rows =
    {
      {1, 31, 301, 0.01, 1, 1, std::nullopt},
      {1, 31, -1, 0.01, 0, 1, std::nullopt},
      {2, 31, 301, 0.8, 1, 1, std::nullopt},
      {2, 31, -1, 0.8, 0, 1, std::nullopt}
    };
    const std::vector<IPFAlignmentRow> alignment_rows = {{7, 1}, {7, 2}};

    const std::vector<IPFResultRow> without_propagation = inference.infer(precursor_rows, transition_rows, alignment_rows, base_config);
    auto without_feature2 = std::find_if(without_propagation.begin(), without_propagation.end(),
                                         [](const auto& row) { return row.feature_id == 2 && row.peptide_id == 301; });
    TEST_EQUAL(without_feature2 != without_propagation.end(), true)
    if (without_feature2 != without_propagation.end())
    {
      TEST_REAL_SIMILAR(without_feature2->pep, 0.5)
    }

    PeptidoformInferenceConfig propagation_config = base_config;
    propagation_config.propagate_signal_across_runs = true;
    propagation_config.across_run_confidence_threshold = 0.05;
    const std::vector<IPFResultRow> with_propagation = inference.infer(precursor_rows, transition_rows, alignment_rows, propagation_config);
    auto with_feature2 = std::find_if(with_propagation.begin(), with_propagation.end(),
                                      [](const auto& row) { return row.feature_id == 2 && row.peptide_id == 301; });
    TEST_EQUAL(with_feature2 != with_propagation.end(), true)
    if (with_feature2 != with_propagation.end() && without_feature2 != without_propagation.end())
    {
      TEST_EQUAL(with_feature2->pep <= 0.01, true)
      TEST_EQUAL(with_feature2->pep < without_feature2->pep, true)
    }
  }

  {
    String tmp_osw;
    NEW_TMP_FILE(tmp_osw);
    copySharedInferenceFixture_(tmp_osw);

    OSWFile osw(tmp_osw);
    PeptidoformInferenceConfig config;
    config.ipf_ms1_scoring = false;
    config.ipf_ms2_scoring = false;

    const std::vector<IPFPrecursorRow> precursor_rows = osw.readIPFPrecursorData(config);
    const std::vector<IPFTransitionRow> transition_rows = osw.readIPFTransitionData(config);
    const std::vector<IPFResultRow> results = inference.infer(precursor_rows, transition_rows, {}, config);

    TEST_EQUAL(results.size(), 333)
    auto pyprophet_row = std::find_if(results.begin(), results.end(),
                                      [](const auto& row) { return row.feature_id == -9078977811506172301LL && row.peptide_id == 305; });
    TEST_EQUAL(pyprophet_row != results.end(), true)
    if (pyprophet_row != results.end())
    {
      TEST_REAL_SIMILAR(pyprophet_row->precursor_peakgroup_pep, 0.00314231908146708)
      TEST_EQUAL(pyprophet_row->qvalue <= 0.01, true)
      TEST_EQUAL(pyprophet_row->pep <= 0.01, true)
    }
  }
}
END_SECTION

END_TEST
