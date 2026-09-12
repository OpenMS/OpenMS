// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/ID/PeptDeepRescoring.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/SYSTEM/File.h>

#include <cmath>

///////////////////////////

using namespace OpenMS;
using namespace std;

// The models the build downloads for the ONNX tests; the class-test CMake copies them
// next to this test's data. Same files PeptDeepInference_test loads.
const std::string ms2_model_path = "data/peptdeep_ms2_dynamic.onnx";
const std::string rt_model_path = "data/peptdeep_rt_dynamic.onnx";

namespace
{
  /// A PSM whose peak annotations are the b/y ions of @p seq, so the observed side is a
  /// plausible spectrum for it rather than noise.
  PeptideHit makeHit_(const std::string& seq, int charge, double score,
                      const std::vector<std::pair<std::string, double>>& peaks)
  {
    PeptideHit hit;
    hit.setSequence(AASequence::fromString(seq));
    hit.setCharge(charge);
    hit.setScore(score);
    std::vector<PeptideHit::PeakAnnotation> anns;
    for (const auto& [name, intensity] : peaks)
    {
      PeptideHit::PeakAnnotation pa;
      pa.annotation = name;
      pa.charge = 1;
      pa.mz = 100.0 + anns.size();
      pa.intensity = intensity;
      anns.push_back(pa);
    }
    hit.setPeakAnnotations(anns);
    return hit;
  }

  PeptideIdentification makeId_(const std::string& run_id, double rt, const PeptideHit& hit)
  {
    PeptideIdentification pi;
    pi.setIdentifier(run_id);
    pi.setRT(rt);
    pi.setHigherScoreBetter(true);
    pi.setHits({hit});
    return pi;
  }
}

START_TEST(PeptDeepRescoring, "$Id$")

/////////////////////////////////////////////////////////////

PeptDeepRescoring* ptr = nullptr;
PeptDeepRescoring* null_ptr = nullptr;

START_SECTION((PeptDeepRescoring()))
  ptr = new PeptDeepRescoring();
  TEST_NOT_EQUAL(ptr, null_ptr)
END_SECTION

START_SECTION((~PeptDeepRescoring()))
  delete ptr;
END_SECTION

START_SECTION((static StringList featureNames()))
{
  StringList names = PeptDeepRescoring::featureNames();
  TEST_EQUAL(names.size(), 5)
  // The order is part of the contract: it is what gets appended to extra_features.
  TEST_STRING_EQUAL(names[0], Constants::UserParam::MS2_COSINE)
  TEST_STRING_EQUAL(names[1], Constants::UserParam::MS2_SPECTRAL_ANGLE)
  TEST_STRING_EQUAL(names[2], Constants::UserParam::MS2_PEARSON)
  TEST_STRING_EQUAL(names[3], Constants::UserParam::MS2_FRAC_PRED_FOUND)
  TEST_STRING_EQUAL(names[4], Constants::UserParam::RT_ABS_ERROR)
}
END_SECTION

START_SECTION([EXTRA] default parameters)
{
  PeptDeepRescoring r;
  Param p = r.getParameters();
  // Models have no sensible default, so they start empty and annotate() refuses to run.
  TEST_STRING_EQUAL(p.getValue("ms2_model").toString(), "")
  TEST_STRING_EQUAL(p.getValue("rt_model").toString(), "")
  // A negative NCE means "select it from the data".
  TEST_REAL_SIMILAR((double)p.getValue("nce"), -1.0)
  TEST_STRING_EQUAL(p.getValue("rt_model_type").toString(), "b_spline")
  TEST_STRING_EQUAL(p.getValue("instrument").toString(), "QE")
  // Nothing has been run yet.
  TEST_REAL_SIMILAR(r.getUsedNCE(), -1.0)
}
END_SECTION

START_SECTION([EXTRA] annotate() refuses to run without models)
{
  PeptDeepRescoring r;
  PeakMap exp;
  vector<ProteinIdentification> prot(1);
  PeptideIdentificationList peps;
  PeptideIdentification pi;
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTIDEK"));
  hit.setCharge(2);
  pi.setHits({hit});
  peps.push_back(pi);

  // Both model paths empty: the features cannot be computed at all.
  TEST_EXCEPTION(Exception::MissingInformation, r.annotate(exp, prot, peps))

  // A configured but absent model must fail loudly rather than silently emit zeros.
  Param p = r.getParameters();
  p.setValue("ms2_model", "does_not_exist_ms2.onnx");
  p.setValue("rt_model", "does_not_exist_rt.onnx");
  r.setParameters(p);
  TEST_EXCEPTION(Exception::FileNotFound, r.annotate(exp, prot, peps))
}
END_SECTION

START_SECTION([EXTRA] annotate() is a no-op on empty input)
{
  PeptDeepRescoring r;
  Param p = r.getParameters();
  p.setValue("ms2_model", OPENMS_GET_TEST_DATA_PATH("PeptDeepRescoring_test.txt"));
  p.setValue("rt_model", OPENMS_GET_TEST_DATA_PATH("PeptDeepRescoring_test.txt"));
  r.setParameters(p);

  PeakMap exp;
  vector<ProteinIdentification> prot(1);
  PeptideIdentificationList peps;
  // No PSMs: returns before touching the ONNX runtime, so the placeholder files above
  // are never actually loaded.
  r.annotate(exp, prot, peps);
  TEST_EQUAL(peps.size(), 0)
}
END_SECTION

START_SECTION([EXTRA] annotate() writes every feature and registers them per run)
{
  // Skip rather than fail if the downloaded models are not beside the test binary: the
  // build normally provides them, but a hand-run binary may not see them.
  if (!File::exists(ms2_model_path) || !File::exists(rt_model_path))
  {
    STATUS("PeptDeep ONNX models not found next to the test; skipping the end-to-end section.")
  }
  else
  {
    // Two runs, so the per-run grouping, calibration and feature registration are all
    // exercised rather than assumed.
    const std::vector<std::pair<std::string, double>> peaks =
      {{"b2", 500.0}, {"y3", 900.0}, {"y4", 750.0}, {"b3", 300.0}, {"y5", 620.0}, {"y6", 410.0}};

    PeptideIdentificationList peps;
    // Enough PSMs that each run's confident half (calibration_quantile 0.5) still clears
    // the four points the retention-time fit needs -- otherwise calibration is skipped
    // and rt_abs_error is a constant zero that asserts nothing.
    const std::vector<std::string> seqs =
      {"PEPTIDEK", "TESTPEPTIDER", "ELVISLIVESK", "SAMPLERPEPTIDEK",
       "LNGGKPVDEK", "VATVSLPRK", "AGGDLSTVEK", "YLDGTSLSPK",
       "IADPEHLVK", "GFTVSGNLTK", "SLYPQEDVK", "NVLDTGAPIK",
       "TFSDLPHLDK", "QALEGSVTLK", "MDSTEPPYSQK", "HVLTSIGEK",
       "WGADLSPVEK", "CFAENGTLVK", "RPLDSNVTEK", "DYSAQLTPGK"};
    for (Size i = 0; i < seqs.size(); ++i)
    {
      const std::string run = (i % 2 == 0) ? "run_A" : "run_B";
      // Retention times increase with index so the calibration has a real trend to fit.
      peps.push_back(makeId_(run, 100.0 + 40.0 * i, makeHit_(seqs[i], 2, 20.0 - i, peaks)));
    }

    std::vector<ProteinIdentification> prot(2);
    prot[0].setIdentifier("run_A");
    prot[1].setIdentifier("run_B");

    PeptDeepRescoring r;
    Param p = r.getParameters();
    p.setValue("ms2_model", ms2_model_path);
    p.setValue("rt_model", rt_model_path);
    // Fixing the NCE keeps the test off the grid scan, which would otherwise dominate
    // its runtime for no extra coverage of the write-back path.
    p.setValue("nce", 30.0);
    p.setValue("rt_model_type", "linear");
    r.setParameters(p);

    PeakMap exp;
    r.annotate(exp, prot, peps);

    // Every hit carries every feature. A feature present on only some PSMs of a run is
    // dropped for all of them by Percolator, so partial annotation is a real failure.
    for (const PeptideIdentification& pi : peps)
    {
      for (const PeptideHit& hit : pi.getHits())
      {
        for (const std::string& f : PeptDeepRescoring::featureNames())
        {
          TEST_EQUAL(hit.metaValueExists(f), true)
        }
      }
    }

    // Both runs must be told about the features, not just the first ProteinIdentification.
    for (const ProteinIdentification& prot_run : prot)
    {
      const ProteinIdentification::SearchParameters& sp = prot_run.getSearchParameters();
      TEST_EQUAL(sp.metaValueExists("extra_features"), true)
      StringList registered =
        ListUtils::create<std::string>(sp.getMetaValue("extra_features").toString(), ',');
      for (const std::string& f : PeptDeepRescoring::featureNames())
      {
        TEST_EQUAL(std::find(registered.begin(), registered.end(), f) != registered.end(), true)
      }
    }

    // The fixed NCE must be reported back rather than a stale or unset value.
    TEST_REAL_SIMILAR(r.getUsedNCE(), 30.0)

    // A negative calibration error means the fit was skipped, which would leave
    // rt_abs_error a constant zero and make the checks below vacuous.
    TEST_EQUAL(r.getRTCalibrationError() >= 0.0, true)

    // Cosine is a similarity, so it has to stay in range whatever the model emits.
    for (const PeptideIdentification& pi : peps)
    {
      const double cos = pi.getHits()[0].getMetaValue(Constants::UserParam::MS2_COSINE);
      TEST_EQUAL(cos >= -1.0 && cos <= 1.0, true)
      const double frac = pi.getHits()[0].getMetaValue(Constants::UserParam::MS2_FRAC_PRED_FOUND);
      TEST_EQUAL(frac >= 0.0 && frac <= 1.0, true)
      const double rt_err = pi.getHits()[0].getMetaValue(Constants::UserParam::RT_ABS_ERROR);
      TEST_EQUAL(rt_err >= 0.0, true)
    }
  }
}
END_SECTION

START_SECTION([EXTRA] a run without peak annotations is reported rather than silently zeroed)
{
  if (!File::exists(ms2_model_path) || !File::exists(rt_model_path))
  {
    STATUS("PeptDeep ONNX models not found next to the test; skipping the mixed-annotation section.")
  }
  else
  {
    // run_A carries annotations, run_B does not. The global "no PSM has annotations"
    // guard does not fire here, so without a per-run check run_B would receive four
    // MS2 features that are exactly zero for every PSM -- a constant block Percolator
    // reads as a run separator rather than as evidence.
    const std::vector<std::pair<std::string, double>> peaks =
      {{"b2", 400.0}, {"y3", 900.0}, {"y4", 650.0}, {"b3", 280.0}, {"y5", 520.0}};
    const std::vector<std::pair<std::string, double>> none;

    PeptideIdentificationList peps;
    const std::vector<std::string> seqs =
      {"PEPTIDEK", "TESTPEPTIDER", "ELVISLIVESK", "LNGGKPVDEK",
       "VATVSLPRK", "AGGDLSTVEK", "YLDGTSLSPK", "IADPEHLVK"};
    for (Size i = 0; i < seqs.size(); ++i)
    {
      const bool run_a = (i % 2 == 0);
      peps.push_back(makeId_(run_a ? "run_A" : "run_B", 150.0 + 35.0 * i,
                             makeHit_(seqs[i], 2, 20.0 - i, run_a ? peaks : none)));
    }

    std::vector<ProteinIdentification> prot(2);
    prot[0].setIdentifier("run_A");
    prot[1].setIdentifier("run_B");

    PeptDeepRescoring r;
    Param p = r.getParameters();
    p.setValue("ms2_model", ms2_model_path);
    p.setValue("rt_model", rt_model_path);
    p.setValue("nce", 30.0);
    p.setValue("rt_model_type", "linear");
    r.setParameters(p);

    PeakMap exp;
    // Must not throw: one run is annotated, so the input is usable -- the unannotated
    // run is a warning, not a fatal error.
    r.annotate(exp, prot, peps);

    // The features are still written everywhere (dropping them for one run would leave
    // Percolator with a ragged feature set), but run_B's MS2 values are all zero and
    // run_A's are not. That asymmetry is exactly what the warning exists to announce.
    double run_a_cosine_sum = 0.0, run_b_cosine_sum = 0.0;
    for (const PeptideIdentification& pi : peps)
    {
      for (const PeptideHit& hit : pi.getHits())
      {
        TEST_EQUAL(hit.metaValueExists(Constants::UserParam::MS2_COSINE), true)
        const double c = hit.getMetaValue(Constants::UserParam::MS2_COSINE);
        if (pi.getIdentifier() == "run_A") { run_a_cosine_sum += c; } else { run_b_cosine_sum += c; }
      }
    }
    TEST_REAL_SIMILAR(run_b_cosine_sum, 0.0)
    TEST_EQUAL(run_a_cosine_sum > 0.0, true)
  }
}
END_SECTION

START_SECTION([EXTRA] every run lacking annotations is still a hard error)
{
  PeptDeepRescoring r;
  Param p = r.getParameters();
  p.setValue("ms2_model", ms2_model_path);
  p.setValue("rt_model", rt_model_path);
  r.setParameters(p);

  // Two runs, neither annotated: nothing can be computed, so this must still throw
  // rather than degrade to a file of zeros.
  PeptideIdentificationList peps;
  peps.push_back(makeId_("run_A", 100.0, makeHit_("PEPTIDEK", 2, 20.0, {})));
  peps.push_back(makeId_("run_B", 200.0, makeHit_("TESTPEPTIDER", 2, 19.0, {})));
  std::vector<ProteinIdentification> prot(2);
  prot[0].setIdentifier("run_A");
  prot[1].setIdentifier("run_B");

  PeakMap exp;
  TEST_EXCEPTION(Exception::MissingInformation, r.annotate(exp, prot, peps))
}
END_SECTION

START_SECTION([EXTRA] automatic NCE selection scans the grid and reports its choice)
{
  if (!File::exists(ms2_model_path) || !File::exists(rt_model_path))
  {
    STATUS("PeptDeep ONNX models not found next to the test; skipping the NCE section.")
  }
  else
  {
    // No collision energy anywhere (the PeakMap is empty), so the grid centres on
    // nce_fallback and the scan runs over centre +/- halfwidth. This is the only section
    // that reaches the candidate loop; a fixed NCE skips it entirely.
    const std::vector<std::pair<std::string, double>> peaks =
      {{"b2", 400.0}, {"y3", 900.0}, {"y4", 650.0}, {"b3", 280.0}, {"y5", 520.0}};

    PeptideIdentificationList peps;
    const std::vector<std::string> seqs =
      {"PEPTIDEK", "TESTPEPTIDER", "ELVISLIVESK", "LNGGKPVDEK",
       "VATVSLPRK", "AGGDLSTVEK", "YLDGTSLSPK", "IADPEHLVK"};
    for (Size i = 0; i < seqs.size(); ++i)
    {
      peps.push_back(makeId_("run_A", 200.0 + 30.0 * i, makeHit_(seqs[i], 2, 20.0 - i, peaks)));
    }

    std::vector<ProteinIdentification> prot(1);
    prot[0].setIdentifier("run_A");

    PeptDeepRescoring r;
    Param p = r.getParameters();
    p.setValue("ms2_model", ms2_model_path);
    p.setValue("rt_model", rt_model_path);
    p.setValue("nce", -1.0);            // negative: select from the data
    p.setValue("nce_fallback", 30.0);
    p.setValue("nce_grid_halfwidth", 6.0);
    p.setValue("nce_grid_step", 2.0);
    p.setValue("rt_model_type", "linear");
    r.setParameters(p);

    PeakMap exp;
    r.annotate(exp, prot, peps);

    // Whatever it picked has to come from the grid it was given.
    const double used = r.getUsedNCE();
    TEST_EQUAL(used >= 24.0 && used <= 36.0, true)
    // And it must be a grid point, not an interpolation.
    TEST_REAL_SIMILAR(std::fmod(used - 24.0, 2.0), 0.0)
  }
}
END_SECTION

START_SECTION([EXTRA] a matching peptide scores a higher cosine than a mismatched one)
{
  if (!File::exists(ms2_model_path) || !File::exists(rt_model_path))
  {
    STATUS("PeptDeep ONNX models not found next to the test; skipping the discrimination section.")
  }
  else
  {
    // Give both hits the annotations of the SAME observed spectrum. Only the first hit's
    // sequence produced it, so its prediction should agree with it better than the other's.
    // This is what the feature exists to do; if it fails, the ion mapping is wrong.
    const std::string observed_seq = "SAMPLERPEPTIDEK";
    const std::vector<std::pair<std::string, double>> peaks =
      {{"y2", 800.0}, {"y3", 950.0}, {"y4", 700.0}, {"y5", 500.0},
       {"b2", 300.0}, {"b3", 450.0}, {"b4", 260.0}};

    PeptideIdentificationList peps;
    peps.push_back(makeId_("run_A", 300.0, makeHit_(observed_seq, 2, 20.0, peaks)));
    peps.push_back(makeId_("run_A", 300.0, makeHit_("KEDITPEPRELPMAS", 2, 19.0, peaks)));
    // A few more so the run has enough points to calibrate at all.
    const std::vector<std::string> filler = {"PEPTIDEK", "TESTPEPTIDER", "ELVISLIVESK", "AGGDLSTVEK"};
    for (Size i = 0; i < filler.size(); ++i)
    {
      peps.push_back(makeId_("run_A", 400.0 + 50.0 * i, makeHit_(filler[i], 2, 18.0 - i, peaks)));
    }

    std::vector<ProteinIdentification> prot(1);
    prot[0].setIdentifier("run_A");

    PeptDeepRescoring r;
    Param p = r.getParameters();
    p.setValue("ms2_model", ms2_model_path);
    p.setValue("rt_model", rt_model_path);
    p.setValue("nce", 30.0);
    p.setValue("rt_model_type", "linear");
    r.setParameters(p);

    PeakMap exp;
    r.annotate(exp, prot, peps);

    const double matched = peps[0].getHits()[0].getMetaValue(Constants::UserParam::MS2_COSINE);
    const double mismatched = peps[1].getHits()[0].getMetaValue(Constants::UserParam::MS2_COSINE);
    TEST_EQUAL(matched > mismatched, true)
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
END_TEST
