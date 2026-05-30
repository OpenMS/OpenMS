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
  String out_spectra;
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
  String out_spectra;
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
  String out_spectra;
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
  String out_spectra;
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

START_SECTION((static double computeHyperScore(double, bool, const MSSpectrum&, const MSSpectrum&, double)))
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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
