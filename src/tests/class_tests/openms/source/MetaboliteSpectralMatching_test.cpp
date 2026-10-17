// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/MetaboliteSpectralMatching.h>
///////////////////////////

#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/Precursor.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <cmath>
#include <limits>

using namespace OpenMS;
using namespace std;

// helper to build a small spectral library with two compounds at almost identical precursor m/z
// but different CCS values (compound_A: 150, compound_B: 200), stored as "CCS" meta value.
static PeakMap makeLibrary()
{
  PeakMap spec_db;

  MSSpectrum lib1;
  lib1.setName("compound_A");
  lib1.setMetaValue("Metabolite_Name", "compound_A");
  lib1.setMetaValue("CCS", "150.0");
  // run() reads these string meta values for every scored candidate and converts them to String,
  // which throws for non-string/absent values - so populate them as a real MSP/GNPS library would.
  lib1.setMetaValue("HMDB_ID", "HMDB0000001");
  lib1.setMetaValue(Constants::UserParam::MSM_SUM_FORMULA, "C10H15N");
  lib1.setMetaValue(Constants::UserParam::MSM_INCHI_STRING, "InChI=1S/compoundA");
  lib1.setMetaValue(Constants::UserParam::MSM_SMILES_STRING, "CCO");
  lib1.setMetaValue(Constants::UserParam::MSM_PRECURSOR_ADDUCT, "[M+H]+");
  Precursor p1;
  p1.setMZ(300.0);
  p1.setCharge(1);
  lib1.setPrecursors({p1});
  lib1.push_back(Peak1D(100.0, 999.0));
  lib1.push_back(Peak1D(150.0, 500.0));
  lib1.push_back(Peak1D(200.0, 300.0));
  lib1.push_back(Peak1D(250.0, 100.0));
  lib1.setRT(0);
  spec_db.addSpectrum(lib1);

  MSSpectrum lib2;
  lib2.setName("compound_B");
  lib2.setMetaValue("Metabolite_Name", "compound_B");
  lib2.setMetaValue("CCS", "200.0");
  lib2.setMetaValue("HMDB_ID", "HMDB0000002");
  lib2.setMetaValue(Constants::UserParam::MSM_SUM_FORMULA, "C12H17N");
  lib2.setMetaValue(Constants::UserParam::MSM_INCHI_STRING, "InChI=1S/compoundB");
  lib2.setMetaValue(Constants::UserParam::MSM_SMILES_STRING, "CCCO");
  lib2.setMetaValue(Constants::UserParam::MSM_PRECURSOR_ADDUCT, "[M+H]+");
  Precursor p2;
  p2.setMZ(300.1);
  p2.setCharge(1);
  lib2.setPrecursors({p2});
  lib2.push_back(Peak1D(100.0, 999.0));
  lib2.push_back(Peak1D(150.0, 500.0));
  lib2.push_back(Peak1D(200.0, 300.0));
  lib2.push_back(Peak1D(250.0, 100.0));
  lib2.setRT(1);
  spec_db.addSpectrum(lib2);

  return spec_db;
}

static Param makeParams(double ccs_error_percent)
{
  // start from the algorithm defaults so the fixture inherits future parameter additions,
  // then override only the knobs these tests rely on.
  Param params = MetaboliteSpectralMatching().getDefaults();
  params.setValue("prec_mass_error_value", 1.0);
  params.setValue("frag_mass_error_value", 0.5);
  params.setValue("mass_error_unit", "Da");
  params.setValue("ionization_mode", "positive");
  params.setValue("report_mode", "top3");
  params.setValue("merge_spectra", "false");
  params.setValue("ccs_error_percent", ccs_error_percent);
  return params;
}

START_TEST(MetaboliteSpectralMatching, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MetaboliteSpectralMatching* ptr = nullptr;
MetaboliteSpectralMatching* null_ptr = nullptr;
START_SECTION(MetaboliteSpectralMatching())
{
	ptr = new MetaboliteSpectralMatching();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~MetaboliteSpectralMatching())
{
	delete ptr;
}
END_SECTION

START_SECTION(SpectralMatch CCS getters and setters)
{
  SpectralMatch sm;
  // default CCS values should be "not set"
  TEST_REAL_SIMILAR(sm.getObservedCCS(), IMTypes::DRIFTTIME_NOT_SET)
  TEST_REAL_SIMILAR(sm.getFoundCCS(), IMTypes::DRIFTTIME_NOT_SET)

  sm.setObservedCCS(167.3);
  sm.setFoundCCS(170.1);
  TEST_REAL_SIMILAR(sm.getObservedCCS(), 167.3)
  TEST_REAL_SIMILAR(sm.getFoundCCS(), 170.1)

  // copy constructor
  SpectralMatch sm2(sm);
  TEST_REAL_SIMILAR(sm2.getObservedCCS(), 167.3)
  TEST_REAL_SIMILAR(sm2.getFoundCCS(), 170.1)

  // assignment operator
  SpectralMatch sm3;
  sm3 = sm;
  TEST_REAL_SIMILAR(sm3.getObservedCCS(), 167.3)
  TEST_REAL_SIMILAR(sm3.getFoundCCS(), 170.1)
}
END_SECTION

START_SECTION(MetaboliteSpectralMatching default parameters)
{
  MetaboliteSpectralMatching msm;
  Param p = msm.getDefaults();
  // CCS filtering is disabled by default
  TEST_REAL_SIMILAR((double)p.getValue("ccs_error_percent"), 0.0)
}
END_SECTION

START_SECTION(CCS filtering disabled reports all matches with CCS columns)
{
  PeakMap spec_db = makeLibrary();

  // experimental spectrum with CCS = 152.0 (as "CCS" meta value)
  PeakMap msexp;
  MSSpectrum exp_spec;
  exp_spec.setMSLevel(2);
  exp_spec.setMetaValue("CCS", "152.0");
  Precursor exp_prec;
  exp_prec.setMZ(300.05);
  exp_prec.setCharge(1);
  exp_spec.setPrecursors({exp_prec});
  exp_spec.push_back(Peak1D(100.0, 800.0));
  exp_spec.push_back(Peak1D(150.0, 400.0));
  exp_spec.push_back(Peak1D(200.0, 250.0));
  exp_spec.push_back(Peak1D(250.0, 80.0));
  exp_spec.setRT(10.0);
  msexp.addSpectrum(exp_spec);

  MetaboliteSpectralMatching msm;
  msm.setParameters(makeParams(0.0)); // filtering disabled

  MzTab mztab;
  std::string out_spectra;
  msm.run(msexp, spec_db, mztab, out_spectra);

  const auto& rows = mztab.getSmallMoleculeSectionRows();
  TEST_EQUAL(rows.size(), 2) // both compounds match

  for (const auto& row : rows)
  {
    bool found_obs_ccs = false, found_lib_ccs = false, found_ccs_err = false;
    for (const auto& opt : row.opt_)
    {
      if (opt.first == "opt_observed_ccs")
      {
        TEST_EQUAL(opt.second.toCellString(), "152.0")
        found_obs_ccs = true;
      }
      if (opt.first == "opt_library_ccs")
      {
        TEST_EQUAL(opt.second.toCellString() == "150.0" || opt.second.toCellString() == "200.0", true)
        found_lib_ccs = true;
      }
      if (opt.first == "opt_ccs_error_percent")
      {
        TEST_NOT_EQUAL(opt.second.toCellString(), "null")
        found_ccs_err = true;
      }
    }
    TEST_EQUAL(found_obs_ccs, true)
    TEST_EQUAL(found_lib_ccs, true)
    TEST_EQUAL(found_ccs_err, true)
  }
}
END_SECTION

START_SECTION(CCS filtering enabled keeps only within-tolerance matches)
{
  PeakMap spec_db = makeLibrary();

  PeakMap msexp;
  MSSpectrum exp_spec;
  exp_spec.setMSLevel(2);
  exp_spec.setMetaValue("CCS", "152.0"); // within 5% of compound_A (150), far from compound_B (200)
  Precursor exp_prec;
  exp_prec.setMZ(300.05);
  exp_prec.setCharge(1);
  exp_spec.setPrecursors({exp_prec});
  exp_spec.push_back(Peak1D(100.0, 800.0));
  exp_spec.push_back(Peak1D(150.0, 400.0));
  exp_spec.push_back(Peak1D(200.0, 250.0));
  exp_spec.push_back(Peak1D(250.0, 80.0));
  exp_spec.setRT(10.0);
  msexp.addSpectrum(exp_spec);

  MetaboliteSpectralMatching msm;
  msm.setParameters(makeParams(5.0)); // 5% CCS tolerance

  MzTab mztab;
  std::string out_spectra;
  msm.run(msexp, spec_db, mztab, out_spectra);

  const auto& rows = mztab.getSmallMoleculeSectionRows();
  TEST_EQUAL(rows.size(), 1) // only compound_A matches

  if (rows.size() == 1)
  {
    for (const auto& opt : rows[0].opt_)
    {
      if (opt.first == "opt_library_ccs") { TEST_EQUAL(opt.second.toCellString(), "150.0") }
      if (opt.first == "opt_observed_ccs") { TEST_EQUAL(opt.second.toCellString(), "152.0") }
    }
  }
}
END_SECTION

START_SECTION(CCS filtering converts 1/K0 (VSSC) drift time to CCS)
{
  PeakMap spec_db = makeLibrary();

  // experimental spectrum carries a 1/K0 (VSSC) drift time instead of a CCS meta value.
  // For m/z 300.05, z=1, 1/K0 = 0.726 converts to ~152 Angstrom^2 (within 5% of compound_A only).
  PeakMap msexp;
  MSSpectrum exp_spec;
  exp_spec.setMSLevel(2);
  exp_spec.setDriftTime(0.726);
  exp_spec.setDriftTimeUnit(DriftTimeUnit::VSSC);
  Precursor exp_prec;
  exp_prec.setMZ(300.05);
  exp_prec.setCharge(1);
  exp_spec.setPrecursors({exp_prec});
  exp_spec.push_back(Peak1D(100.0, 800.0));
  exp_spec.push_back(Peak1D(150.0, 400.0));
  exp_spec.push_back(Peak1D(200.0, 250.0));
  exp_spec.push_back(Peak1D(250.0, 80.0));
  exp_spec.setRT(10.0);
  msexp.addSpectrum(exp_spec);

  MetaboliteSpectralMatching msm;
  msm.setParameters(makeParams(5.0));

  MzTab mztab;
  std::string out_spectra;
  msm.run(msexp, spec_db, mztab, out_spectra);

  // Without the 1/K0 -> CCS conversion the observed CCS would be unknown and both compounds
  // would pass; getting a single row proves the conversion + filtering worked.
  const auto& rows = mztab.getSmallMoleculeSectionRows();
  TEST_EQUAL(rows.size(), 1)

  if (rows.size() == 1)
  {
    for (const auto& opt : rows[0].opt_)
    {
      if (opt.first == "opt_library_ccs") { TEST_EQUAL(opt.second.toCellString(), "150.0") }
      if (opt.first == "opt_observed_ccs") { TEST_NOT_EQUAL(opt.second.toCellString(), "null") }
    }
  }
}
END_SECTION

START_SECTION(CCS filtering keeps matches when CCS is missing)
{
  PeakMap spec_db = makeLibrary();

  // experimental spectrum without any CCS information
  PeakMap msexp;
  MSSpectrum exp_spec;
  exp_spec.setMSLevel(2);
  Precursor exp_prec;
  exp_prec.setMZ(300.05);
  exp_prec.setCharge(1);
  exp_spec.setPrecursors({exp_prec});
  exp_spec.push_back(Peak1D(100.0, 800.0));
  exp_spec.push_back(Peak1D(150.0, 400.0));
  exp_spec.push_back(Peak1D(200.0, 250.0));
  exp_spec.push_back(Peak1D(250.0, 80.0));
  exp_spec.setRT(10.0);
  msexp.addSpectrum(exp_spec);

  MetaboliteSpectralMatching msm;
  msm.setParameters(makeParams(5.0)); // enabled, but no experimental CCS -> not filtered

  MzTab mztab;
  std::string out_spectra;
  msm.run(msexp, spec_db, mztab, out_spectra);

  const auto& rows = mztab.getSmallMoleculeSectionRows();
  TEST_EQUAL(rows.size(), 2) // both kept (CCS does not penalise non-IM data)

  for (const auto& row : rows)
  {
    for (const auto& opt : row.opt_)
    {
      if (opt.first == "opt_observed_ccs") { TEST_EQUAL(opt.second.toCellString(), "null") }
      if (opt.first == "opt_ccs_error_percent") { TEST_EQUAL(opt.second.toCellString(), "null") }
    }
  }
}
END_SECTION

START_SECTION((static double computeHyperScore(double, bool, const MSSpectrum&, const MSSpectrum&, double, bool)))
{
  MSSpectrum exp_spec;
  exp_spec.push_back(Peak1D(100.0, 100.0));
  exp_spec.push_back(Peak1D(200.0, 200.0));
  exp_spec.push_back(Peak1D(300.0, 300.0));
  exp_spec.push_back(Peak1D(400.0, 50.0));

  MSSpectrum db_spec;
  db_spec.push_back(Peak1D(100.0, 999.0));
  db_spec.push_back(Peak1D(200.0, 999.0));
  db_spec.push_back(Peak1D(300.0, 999.0));
  db_spec.push_back(Peak1D(400.0, 999.0));

  double score = MetaboliteSpectralMatching::computeHyperScore(0.5, false, exp_spec, db_spec, 0.0);
  TEST_EQUAL(score > 0, true)

  // an empty experimental spectrum scores 0
  MSSpectrum empty_spec;
  score = MetaboliteSpectralMatching::computeHyperScore(0.5, false, empty_spec, db_spec, 0.0);
  TEST_REAL_SIMILAR(score, 0.0)
}
END_SECTION

START_SECTION(([EXTRA] computeHyperScore rejects invalid tolerances and accepts exact matching))
{
  MSSpectrum spectrum;
  for (double mz : {100.0, 200.0, 300.0})
  {
    spectrum.push_back(Peak1D(mz, 10.0));
  }
  const MSSpectrum empty;
  const MSSpectrum* experimental_spectra[] = {&spectrum, &empty};
  for (bool ppm : {false, true})
  {
    for (const MSSpectrum* experimental : experimental_spectra)
    {
      for (double tolerance : {-1.0, std::numeric_limits<double>::quiet_NaN(),
                               std::numeric_limits<double>::infinity(),
                               -std::numeric_limits<double>::infinity()})
      {
        std::vector<PeptideHit::PeakAnnotation> annotations(1);
        TEST_EXCEPTION(Exception::InvalidParameter,
          MetaboliteSpectralMatching::computeHyperScore(tolerance, ppm, *experimental, spectrum))
        TEST_EXCEPTION(Exception::InvalidParameter,
          MetaboliteSpectralMatching::computeHyperScore(tolerance, ppm, *experimental, spectrum, annotations))
        TEST_EQUAL(annotations.size(), 1)
      }
    }
    for (double tolerance : {0.0, 0.5})
    {
      TEST_REAL_SIMILAR(MetaboliteSpectralMatching::computeHyperScore(tolerance, ppm, spectrum, spectrum),
        std::log(300.0) + std::log(6.0))
      TEST_REAL_SIMILAR(MetaboliteSpectralMatching::computeHyperScore(tolerance, ppm, empty, spectrum), 0.0)
    }
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION((static double computeMVHScore(double, bool, const MSSpectrum&, const MSSpectrum&, Size, double, double, bool, bool)))
{
  // Intensity classes must follow intensity, not the spectrum's m/z order.
  MSSpectrum exp_spec;
  for (Size i = 0; i < 7; ++i)
  {
    exp_spec.push_back(Peak1D(100.0 + i, i == 5 ? 1000.0 : 10.0 + i));
  }
  MSSpectrum strong_match;
  strong_match.push_back(Peak1D(105.0, 1.0));
  MSSpectrum weak_match;
  weak_match.push_back(Peak1D(100.0, 1.0));
  const double strong_score = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, strong_match, 3, 1.0, 0.0, false);
  const double weak_score = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, weak_match, 3, 1.0, 0.0, false);
  TEST_TRUE(strong_score > weak_score)

  // All retained peaks must be divided in the 1:2:4 ratio. With 70 peaks,
  // ranks 5 and 15 belong to different classes (10/20/40), not one remainder.
  MSSpectrum dense_spec;
  for (Size i = 0; i < 70; ++i)
  {
    dense_spec.push_back(Peak1D(100.0 + i, 70.0 - i));
  }
  MSSpectrum top_class_match;
  top_class_match.push_back(Peak1D(104.0, 1.0));
  MSSpectrum second_class_match;
  second_class_match.push_back(Peak1D(114.0, 1.0));
  const double top_score = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, dense_spec, top_class_match, 3, 1.0, 0.0, false);
  const double second_score = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, dense_spec, second_class_match, 3, 1.0, 0.0, false);
  TEST_TRUE(top_score > second_score)

  // Two ions matching the same peak must receive the same credit as one
  // matched ion and one unmatched ion, with the theoretical count held fixed.
  // The strongest class contains only one peak, so double credit previously
  // made its binomial coefficient invalid and returned a zero score.
  MSSpectrum shared_peak;
  shared_peak.push_back(Peak1D(104.98, 1.0));
  shared_peak.push_back(Peak1D(105.005, 1.0));
  MSSpectrum single_match;
  single_match.push_back(Peak1D(105.005, 1.0));
  single_match.push_back(Peak1D(105.5, 1.0));
  for (bool ppm : {false, true})
  {
    const double tolerance = ppm ? 1000.0 : 0.1;
    for (bool accuracy : {false, true})
    {
      const double shared_score = MetaboliteSpectralMatching::computeMVHScore(
        tolerance, ppm, exp_spec, shared_peak, 3, 1.0, 0.0, accuracy);
      const double single_score = MetaboliteSpectralMatching::computeMVHScore(
        tolerance, ppm, exp_spec, single_match, 3, 1.0, 0.0, accuracy);
      TEST_TRUE(shared_score > 0.0)
      TEST_REAL_SIMILAR(shared_score, single_score)
    }
  }

  // Even when duplicate credits fit within a class, they must not inflate
  // either the hypergeometric term or the mass-accuracy bonus.
  const double shared_score = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, shared_peak, 1, 1.0, 0.0, true);
  const double single_score = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, single_match, 1, 1.0, 0.0, true);
  TEST_REAL_SIMILAR(shared_score, single_score)

  // The disabled mode reproduces the historical per-ion score, including
  // its zero result when duplicate credits exceed an intensity class size.
  const double legacy_invalid = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, shared_peak, 3, 1.0, 0.0, true, false);
  TEST_REAL_SIMILAR(legacy_invalid, 0.0)
  const double legacy_shared = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, shared_peak, 1, 1.0, 0.0, true, false);
  const Size locations = static_cast<Size>(((106.0 + 0.1) - (100.0 - 0.1)) / 0.1);
  const double accuracy = (exp(-2.0 * pow((105.0 - 104.98) / 0.1, 2)) +
                           exp(-2.0 * pow((105.0 - 105.005) / 0.1, 2))) / 2.0;
  const double legacy_expected = log(double(locations * (locations - 1)) / 42.0) + accuracy * log(3.0);
  TEST_REAL_SIMILAR(legacy_shared, legacy_expected)
  TEST_TRUE(legacy_shared > shared_score)
  const double explicit_single = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, shared_peak, 1, 1.0, 0.0, true, true);
  TEST_REAL_SIMILAR(explicit_single, shared_score)

  // A matched peak excluded by the TIC cutoff supplies no accuracy credit.
  MSSpectrum outside_tic;
  outside_tic.push_back(Peak1D(100.0, 1.0));
  outside_tic.push_back(Peak1D(105.0, 1.0));
  MSSpectrum unmatched_low_peak;
  unmatched_low_peak.push_back(Peak1D(100.5, 1.0));
  unmatched_low_peak.push_back(Peak1D(105.0, 1.0));
  const double outside_score = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, outside_tic, 1, 0.5, 0.0, true);
  const double unmatched_score = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, unmatched_low_peak, 1, 0.5, 0.0, true);
  TEST_REAL_SIMILAR(outside_score, unmatched_score)
}
END_SECTION

START_SECTION((static double computeMVHScore(double, bool, const MSSpectrum&, const MSSpectrum&, std::vector<PeptideHit::PeakAnnotation>&, Size, double, double, bool, bool)))
{
  MSSpectrum exp_spec;
  for (Size i = 0; i < 7; ++i)
  {
    exp_spec.push_back(Peak1D(100.0 + i, i == 5 ? 1000.0 : 10.0 + i));
  }
  MSSpectrum db_spec;
  db_spec.push_back(Peak1D(104.98, 1.0));
  db_spec.push_back(Peak1D(105.005, 1.0));
  MSSpectrum::StringDataArray names;
  names.setName("IonNames");
  names.push_back("c1");
  names.push_back("c2");
  db_spec.getStringDataArrays().push_back(names);
  MSSpectrum::IntegerDataArray charges;
  charges.setName("Charges");
  charges.push_back(-1);
  charges.push_back(-2);
  db_spec.getIntegerDataArrays().push_back(charges);

  vector<PeptideHit::PeakAnnotation> annotations;
  const double annotated_score = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, db_spec, annotations, 3, 1.0, 0.0, true);
  const double plain_score = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, db_spec, 3, 1.0, 0.0, true);
  TEST_TRUE(annotated_score > 0.0)
  TEST_REAL_SIMILAR(annotated_score, plain_score)
  TEST_EQUAL(annotations.size(), 2)
  if (annotations.size() == 2)
  {
    TEST_STRING_EQUAL(annotations[0].annotation, "c1")
    TEST_EQUAL(annotations[0].charge, -1)
    TEST_REAL_SIMILAR(annotations[0].mz, 105.0)
    TEST_STRING_EQUAL(annotations[1].annotation, "c2")
    TEST_EQUAL(annotations[1].charge, -2)
    TEST_REAL_SIMILAR(annotations[1].mz, 105.0)
  }
  vector<PeptideHit::PeakAnnotation> legacy_annotations;
  const double legacy_annotated = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, db_spec, legacy_annotations, 1, 1.0, 0.0, true, false);
  const double legacy_plain = MetaboliteSpectralMatching::computeMVHScore(
    0.1, false, exp_spec, db_spec, 1, 1.0, 0.0, true, false);
  TEST_REAL_SIMILAR(legacy_annotated, legacy_plain)
  TEST_TRUE(legacy_annotations == annotations)

  // Theoretical ions in the mass range, but with no fragment matches, must
  // not produce a positive identification score or matching annotations.
  MSSpectrum no_match = db_spec;
  no_match[0].setMZ(100.5);
  no_match[1].setMZ(101.5);
  for (bool ppm : {false, true})
  {
    const double tolerance = ppm ? 100.0 : 0.1;
    for (bool accuracy : {false, true})
    {
      for (bool single_credit : {false, true})
      {
        vector<PeptideHit::PeakAnnotation> no_annotations;
        const double no_match_score = MetaboliteSpectralMatching::computeMVHScore(
          tolerance, ppm, exp_spec, no_match, no_annotations, 3, 1.0, 0.0, accuracy, single_credit);
        TEST_EQUAL(no_match_score, 0.0)
        TEST_TRUE(no_annotations.empty())
        TEST_EQUAL(MetaboliteSpectralMatching::computeMVHScore(
          tolerance, ppm, exp_spec, no_match, 3, 1.0, 0.0, accuracy, single_credit), 0.0)
      }
    }
  }
}
END_SECTION

END_TEST
