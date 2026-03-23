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
#include <OpenMS/FORMAT/MzTab.h>

using namespace OpenMS;
using namespace std;

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
  // default CCS values should be -1.0
  TEST_REAL_SIMILAR(sm.getObservedCCS(), -1.0)
  TEST_REAL_SIMILAR(sm.getFoundCCS(), -1.0)

  sm.setObservedCCS(167.3);
  sm.setFoundCCS(170.1);
  TEST_REAL_SIMILAR(sm.getObservedCCS(), 167.3)
  TEST_REAL_SIMILAR(sm.getFoundCCS(), 170.1)

  // test copy constructor
  SpectralMatch sm2(sm);
  TEST_REAL_SIMILAR(sm2.getObservedCCS(), 167.3)
  TEST_REAL_SIMILAR(sm2.getFoundCCS(), 170.1)

  // test assignment operator
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
  TEST_REAL_SIMILAR((double)p.getValue("ccs_error_percent"), 0.0)
}
END_SECTION

START_SECTION(CCS filtering in spectral matching)
{
  // Create a simple spectral library at similar precursor m/z but different CCS values
  PeakMap spec_db;

  MSSpectrum lib1;
  lib1.setName("compound_A");
  lib1.setMetaValue("Metabolite_Name", "compound_A");
  lib1.setMetaValue("CCS", "150.0");
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

  // Create experimental spectrum with CCS = 152.0
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

  {
    MetaboliteSpectralMatching msm;
    Param params;
    params.setValue("prec_mass_error_value", 1.0, "");
    params.setValue("frag_mass_error_value", 0.5, "");
    params.setValue("mass_error_unit", "Da", "");
    params.setValue("ionization_mode", "positive", "");
    params.setValue("report_mode", "top3", "");
    params.setValue("merge_spectra", "false", "");
    params.setValue("ccs_error_percent", 0.0, "");
    msm.setParameters(params);

    PeakMap msexp_copy(msexp);
    PeakMap spec_db_copy(spec_db);
    MzTab mztab;
    String out_spectra;
    msm.run(msexp_copy, spec_db_copy, mztab, out_spectra);

    const auto& rows = mztab.getSmallMoleculeSectionRows();
    TEST_EQUAL(rows.size(), 2) // both compounds match

    // verify CCS values in MzTab optional columns
    for (const auto& row : rows)
    {
      bool found_obs_ccs = false;
      bool found_lib_ccs = false;
      bool found_ccs_err = false;
      for (const auto& opt : row.opt_)
      {
        if (opt.first == "opt_observed_ccs")
        {
          TEST_EQUAL(opt.second.toCellString(), "152")
          found_obs_ccs = true;
        }
        if (opt.first == "opt_library_ccs")
        {
          // should be either 150 or 200
          TEST_EQUAL(opt.second.toCellString() == "150" || opt.second.toCellString() == "200", true)
          found_lib_ccs = true;
        }
        if (opt.first == "opt_ccs_error_percent")
        {
          // should be a numeric value, not "null"
          TEST_NOT_EQUAL(opt.second.toCellString(), "null")
          found_ccs_err = true;
        }
      }
      TEST_EQUAL(found_obs_ccs, true)
      TEST_EQUAL(found_lib_ccs, true)
      TEST_EQUAL(found_ccs_err, true)
    }
  }

  // Test 2: With CCS filtering (5%), only compound_A should match
  {
    MetaboliteSpectralMatching msm;
    Param params;
    params.setValue("prec_mass_error_value", 1.0, "");
    params.setValue("frag_mass_error_value", 0.5, "");
    params.setValue("mass_error_unit", "Da", "");
    params.setValue("ionization_mode", "positive", "");
    params.setValue("report_mode", "top3", "");
    params.setValue("merge_spectra", "false", "");
    params.setValue("ccs_error_percent", 5.0, "");
    msm.setParameters(params);

    PeakMap msexp_copy(msexp);
    PeakMap spec_db_copy(spec_db);
    MzTab mztab;
    String out_spectra;
    msm.run(msexp_copy, spec_db_copy, mztab, out_spectra);

    const auto& rows = mztab.getSmallMoleculeSectionRows();
    TEST_EQUAL(rows.size(), 1) // only compound_A matches

    // verify the matched library CCS is 150 (compound_A)
    for (const auto& opt : rows[0].opt_)
    {
      if (opt.first == "opt_library_ccs")
      {
        TEST_EQUAL(opt.second.toCellString(), "150")
      }
      if (opt.first == "opt_observed_ccs")
      {
        TEST_EQUAL(opt.second.toCellString(), "152")
      }
    }
  }

  // Test 3: When experimental spectrum has no CCS
  {
    PeakMap msexp_no_ccs;
    MSSpectrum exp_no_ccs;
    exp_no_ccs.setMSLevel(2);

    Precursor exp_prec_no_ccs;
    exp_prec_no_ccs.setMZ(300.05);
    exp_prec_no_ccs.setCharge(1);
    exp_no_ccs.setPrecursors({exp_prec_no_ccs});
    exp_no_ccs.push_back(Peak1D(100.0, 800.0));
    exp_no_ccs.push_back(Peak1D(150.0, 400.0));
    exp_no_ccs.push_back(Peak1D(200.0, 250.0));
    exp_no_ccs.push_back(Peak1D(250.0, 80.0));
    exp_no_ccs.setRT(10.0);
    msexp_no_ccs.addSpectrum(exp_no_ccs);

    MetaboliteSpectralMatching msm;
    Param params;
    params.setValue("prec_mass_error_value", 1.0, "");
    params.setValue("frag_mass_error_value", 0.5, "");
    params.setValue("mass_error_unit", "Da", "");
    params.setValue("ionization_mode", "positive", "");
    params.setValue("report_mode", "top3", "");
    params.setValue("merge_spectra", "false", "");
    params.setValue("ccs_error_percent", 5.0, "");
    msm.setParameters(params);

    PeakMap spec_db_copy(spec_db);
    MzTab mztab;
    String out_spectra;
    msm.run(msexp_no_ccs, spec_db_copy, mztab, out_spectra);

    const auto& rows = mztab.getSmallMoleculeSectionRows();
    TEST_EQUAL(rows.size(), 2)

    // verify observed CCS is "null" when not available
    for (const auto& opt : rows[0].opt_)
    {
      if (opt.first == "opt_observed_ccs")
      {
        TEST_EQUAL(opt.second.toCellString(), "null")
      }
      if (opt.first == "opt_ccs_error_percent")
      {
        TEST_EQUAL(opt.second.toCellString(), "null")
      }
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

  MSSpectrum empty_spec;
  score = MetaboliteSpectralMatching::computeHyperScore(0.5, false, empty_spec, db_spec, 0.0);
  TEST_REAL_SIMILAR(score, 0.0)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



