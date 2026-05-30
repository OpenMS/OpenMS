// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <map>

///////////////////////////
#include <OpenMS/ANALYSIS/DECHARGING/MetaboliteFeatureDeconvolution.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/ConsensusXMLFile.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>

///////////////////////////


namespace OpenMS
{
  class MetaboliteFeatureDeconvolutionTest
    : public MetaboliteFeatureDeconvolution
  {
  public:
      /// List of adducts used to explain mass differences
      MassExplainer::AdductsType getPotentialAdducts()
      { return potential_adducts_;}
      /// labeling table
      std::map<Size, String> getMapLabels()
      { return map_label_;}

      /// labeling table inverse
      std::map<String, Size> getMapLabelInverse()
      { return map_label_inverse_;}

			/// status of intensity filter for edges
			bool isIntensityFilterEnabled()
      { return enable_intensity_filter_;}

			/// status of charge discovery
			CHARGEMODE getChargeMode()
      { return q_try_;}

  };

}

using namespace OpenMS;
using namespace std;

START_TEST(MetaboliteFeatureDeconvolution, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MetaboliteFeatureDeconvolution* ptr = nullptr;
MetaboliteFeatureDeconvolution* nullPointer = nullptr;
START_SECTION(MetaboliteFeatureDeconvolution())
	ptr = new MetaboliteFeatureDeconvolution();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~MetaboliteFeatureDeconvolution())
	delete ptr;
END_SECTION

START_SECTION([EXTRA](void updateMembers_()))
  MetaboliteFeatureDeconvolutionTest fdt;

  Param p;
  p.setValue("charge_min", 11, "minimal possible charge");
  p.setValue("charge_max", 13, "maximal possible charge");
  p.setValue("retention_max_diff", 1.0, "maximum allowed RT difference between any two features if their relation shall be determined");
  p.setValue("retention_max_diff_local", 2.0, "maxi");
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:0.7","Na:+:0.3","(2)H4H-4:0:0.2:-2:heavy"}, "Ad");
  fdt.setParameters(p);

  {
  MassExplainer::AdductsType adducts = fdt.getPotentialAdducts();
  std::map<Size, String> map = fdt.getMapLabels();
  std::map<String, Size> map_i = fdt.getMapLabelInverse();
  bool b_filter = fdt.isIntensityFilterEnabled();
  MetaboliteFeatureDeconvolution::CHARGEMODE cm = fdt.getChargeMode();

  TEST_EQUAL(adducts.size(), 3)
  TEST_EQUAL(adducts[0].getFormula(), "H1");
  TEST_EQUAL(adducts[0].getRTShift(), 0);
  TEST_EQUAL(adducts[0].getCharge(), 1);
  TEST_REAL_SIMILAR(adducts[0].getLogProb(), log(0.7));
  TEST_EQUAL(adducts[1].getFormula(), "Na1");
  TEST_EQUAL(adducts[1].getRTShift(), 0);
  TEST_EQUAL(adducts[1].getCharge(), 1);
  TEST_REAL_SIMILAR(adducts[1].getLogProb(), log(0.3));
  TEST_EQUAL(adducts[2].getFormula(), "(2)H4H-4");
  TEST_EQUAL(adducts[2].getRTShift(), -2);
  TEST_EQUAL(adducts[2].getCharge(), 0);
  TEST_REAL_SIMILAR(adducts[2].getLogProb(), log(0.2));
  TEST_EQUAL(cm, MetaboliteFeatureDeconvolution::CHARGEMODE::QFROMFEATURE)
  TEST_EQUAL(map.size(), 2)
  TEST_EQUAL(map_i.size(), 2)
  TEST_EQUAL(map[0], "decharged features");
  TEST_EQUAL(map_i["decharged features"], 0);
  TEST_EQUAL(map[1], "heavy");
  TEST_EQUAL(map_i["heavy"], 1);
  TEST_EQUAL(b_filter, false)
  Param p_internal = fdt.getParameters();
  TEST_REAL_SIMILAR((double) p_internal.getValue("retention_max_diff"), 1.0);
  TEST_REAL_SIMILAR((double) p_internal.getValue("retention_max_diff_local"), 1.0);
  }

  // second param set
  p.setValue("charge_min", 11, "minimal possible charge");
  p.setValue("charge_max", 13, "maximal possible charge");
  p.setValue("q_try", "heuristic", "Try dif");
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:0.9","Na:++:0.1"});
  p.setValue("retention_max_diff", 1.0, "maximum ");
  p.setValue("retention_max_diff_local", 1.0, "maxim");
  p.setValue("intensity_filter", "true", "Enable");
  p.setValue("default_map_label", "mylabel", "Label");
  p.setValue("retention_max_diff", 2.0, "maximum allowed RT difference between any two features if their relation shall be determined");
  p.setValue("retention_max_diff_local", 5.0, "maxi");

  fdt.setParameters(p);
  {
  MassExplainer::AdductsType adducts = fdt.getPotentialAdducts();
  std::map<Size, String> map = fdt.getMapLabels();
  std::map<String, Size> map_i = fdt.getMapLabelInverse();
  bool b_filter = fdt.isIntensityFilterEnabled();
  MetaboliteFeatureDeconvolution::CHARGEMODE cm = fdt.getChargeMode();

  TEST_EQUAL(adducts.size(), 2)
  TEST_EQUAL(adducts[0].getFormula(), "H1");
  TEST_EQUAL(adducts[0].getRTShift(), 0);
  TEST_EQUAL(adducts[0].getCharge(), 1);
  TEST_REAL_SIMILAR(adducts[0].getLogProb(), log(0.9));
  TEST_EQUAL(adducts[1].getFormula(), "Na1");
  TEST_EQUAL(adducts[1].getRTShift(), 0);
  TEST_EQUAL(adducts[1].getCharge(), 2);
  TEST_REAL_SIMILAR(adducts[1].getLogProb(), log(0.1));

  TEST_EQUAL(cm, MetaboliteFeatureDeconvolution::CHARGEMODE::QHEURISTIC)
  TEST_EQUAL(map.size(), 1)
  TEST_EQUAL(map_i.size(), 1)
  TEST_EQUAL(map[0], "mylabel");
  TEST_EQUAL(map_i["mylabel"], 0);
  TEST_EQUAL(b_filter, true)
  Param p_internal = fdt.getParameters();
  TEST_REAL_SIMILAR((double) p_internal.getValue("retention_max_diff"), 2.0);
  TEST_REAL_SIMILAR((double) p_internal.getValue("retention_max_diff_local"), 2.0);

  }

END_SECTION

START_SECTION(MetaboliteFeatureDeconvolution(const MetaboliteFeatureDeconvolution &source))
  MetaboliteFeatureDeconvolution fd;
  Param p;
  p.setValue("charge_min", 11, "minimal possible charge");
  p.setValue("charge_max", 13, "maximal possible charge");
  fd.setParameters(p);
  MetaboliteFeatureDeconvolution fd2(fd);
  MetaboliteFeatureDeconvolution fd_untouched;

  TEST_EQUAL(fd2.getParameters(), fd.getParameters())
  TEST_NOT_EQUAL(fd2.getParameters(), fd_untouched.getParameters())

END_SECTION

START_SECTION(MetaboliteFeatureDeconvolution& operator=(const MetaboliteFeatureDeconvolution &source))
  MetaboliteFeatureDeconvolution fd;
  Param p;
  p.setValue("charge_min", 11, "minimal possible charge");
  p.setValue("charge_max", 13, "maximal possible charge");
  fd.setParameters(p);
  MetaboliteFeatureDeconvolution fd2 = fd;
  MetaboliteFeatureDeconvolution fd_untouched;

  TEST_EQUAL(fd2.getParameters(), fd.getParameters())
  TEST_NOT_EQUAL(fd2.getParameters(), fd_untouched.getParameters())
END_SECTION


START_SECTION(void compute(const FeatureMapType &fm_in, FeatureMapType &fm_out, ConsensusMap &cons_map, ConsensusMap &cons_map_p))
//_CrtSetDbgFlag(_CrtSetDbgFlag(0)|_CRTDBG_CHECK_ALWAYS_DF);

  MetaboliteFeatureDeconvolution fd;
  Param p;
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:0.7","Na:+:0.3","(2)H4H-4:0:0.2:-2:heavy"}, "Ad");
  p.setValue("mass_max_diff", 0.1);
  p.setValue("use_minority_bound","true","enable bound");
  fd.setParameters(p);

  FeatureMap fm_in, fm_out;
  ConsensusMap cm, cm2;
  FeatureXMLFile fl;
  fl.load(OPENMS_GET_TEST_DATA_PATH("FeatureDeconvolution_easy_input.featureXML"), fm_in);
  fd.compute(fm_in, fm_out, cm, cm2);

  String out_file;
  NEW_TMP_FILE(out_file)
  ConsensusXMLFile c1;
  c1.store(out_file,cm);

  WHITELIST("xml-stylesheet,consensusXML version=,consensusElement id=,<UserParam type=");
  // WHITELIST("xml-stylesheet,map id,consensusElement id=");
  TEST_FILE_SIMILAR(out_file, OPENMS_GET_TEST_DATA_PATH("MetaboliteFeatureDeconvolution_easy_output.consensusXML"));

  //small pos test file with specific ions
  Param p_pos;
  p_pos.setValue("potential_adducts", std::vector<std::string>{"H:+:0.6","Na:+:0.2","NH4:+:0.1","K:+:0.1","C2H3N:0:0.05","H-2O-1:0:0.05","H-1Na:0:0.05"}, "Ad_p");
  p_pos.setValue("charge_min", 1, "minimal possible charge");
  p_pos.setValue("charge_max", 3, "maximal possible charge");
  p_pos.setValue("charge_span_max", 3);
  p_pos.setValue("max_neutrals", 1);
  p_pos.setValue("q_try", "feature");
  p_pos.setValue("mass_max_diff", 0.05);
  p_pos.setValue("retention_max_diff", 1.0);
  p_pos.setValue("retention_max_diff_local", 1.0);
  p_pos.setValue("intensity_filter", "false");
  p_pos.setValue("use_minority_bound", "false");

  fd.setParameters(p_pos);

  FeatureMap fm_p_in, fm_p_out;
  ConsensusMap cm_p, cm_p2;
  FeatureXMLFile fl_p;
  fl_p.load(OPENMS_GET_TEST_DATA_PATH("MetaboliteFeatureDeconvolution_test.featureXML"), fm_p_in);
  fd.compute(fm_p_in, fm_p_out, cm_p, cm_p2);

  String out_file_p;
  NEW_TMP_FILE(out_file_p)
  ConsensusXMLFile c_p;
  c_p.store(out_file_p,cm_p);

  WHITELIST("xml-stylesheet,consensusXML version=,consensusElement id=,<UserParam type=");
  TEST_FILE_SIMILAR(out_file_p, OPENMS_GET_TEST_DATA_PATH("MetaboliteFeatureDeconvolution_pos_output.consensusXML"));


  //small neg test file with specific ions
  Param p_neg;
  p_neg.setValue("potential_adducts", std::vector<std::string>{"H-1:-:0.6","Cl:-:0.2","Br:-:0.2","CH2O2:0:0.05","H-2O-1:0:0.05","H-1Na:0:0.05","H-1K:0:0.05"}, "Ad_n");
  p_neg.setValue("charge_min", -3, "minimal possible charge");
  p_neg.setValue("charge_max", -1, "maximal possible charge");
  p_neg.setValue("charge_span_max", 3);
  p_neg.setValue("max_neutrals", 1);
  p_neg.setValue("q_try", "feature");
  p_neg.setValue("mass_max_diff", 0.05);
  p_neg.setValue("retention_max_diff", 1.0);
  p_neg.setValue("retention_max_diff_local", 1.0);
  p_neg.setValue("intensity_filter", "false");
  p_neg.setValue("use_minority_bound", "false");
  p_neg.setValue("negative_mode", "true");

  fd.setParameters(p_neg);

  FeatureMap fm_n_in, fm_n_out;
  ConsensusMap cm_n, cm_n2;
  FeatureXMLFile fl_n;
  fl_n.load(OPENMS_GET_TEST_DATA_PATH("MetaboliteFeatureDeconvolution_test.featureXML"), fm_n_in);
  fd.compute(fm_n_in, fm_n_out, cm_n, cm_n2);

  String out_file_n;
  NEW_TMP_FILE(out_file_n)
  ConsensusXMLFile c_n;
  c_n.store(out_file_n,cm_n);

  WHITELIST("xml-stylesheet,consensusXML version=,consensusElement id=,<UserParam type=");
  TEST_FILE_SIMILAR(out_file_n, OPENMS_GET_TEST_DATA_PATH("MetaboliteFeatureDeconvolution_neg_output.consensusXML"));


  //small pos test file with specific ions and ppm error
  Param p_pos_ppm;
  p_pos_ppm.setValue("potential_adducts", std::vector<std::string>{"H:+:0.6","Na:+:0.4"}, "Ad_p");
  p_pos_ppm.setValue("charge_min", 1, "minimal possible charge");
  p_pos_ppm.setValue("charge_max", 3, "maximal possible charge");
  p_pos_ppm.setValue("charge_span_max", 3);
  p_pos_ppm.setValue("max_neutrals", 1);
  p_pos_ppm.setValue("q_try", "feature");
  p_pos_ppm.setValue("mass_max_diff", 50.0);
  p_pos_ppm.setValue("unit", "ppm");
  p_pos_ppm.setValue("retention_max_diff", 1.0);
  p_pos_ppm.setValue("retention_max_diff_local", 1.0);
  p_pos_ppm.setValue("intensity_filter", "false");
  p_pos_ppm.setValue("use_minority_bound", "false");

  fd.setParameters(p_pos_ppm);

  FeatureMap fm_ppm_in, fm_ppm_out;
  ConsensusMap cm_ppm, cm_ppm2;
  FeatureXMLFile fl_ppm;
  fl_ppm.load(OPENMS_GET_TEST_DATA_PATH("MetaboliteFeatureDeconvolution_test_ppm.featureXML"), fm_ppm_in);
  fd.compute(fm_ppm_in, fm_ppm_out, cm_ppm, cm_ppm2);

  String out_file_ppm;
  NEW_TMP_FILE(out_file_ppm)
  ConsensusXMLFile f_ppm;
  f_ppm.store(out_file_ppm,cm_ppm);

  WHITELIST("xml-stylesheet,consensusXML version=,consensusElement id=,<UserParam type=");
  TEST_FILE_SIMILAR(out_file_ppm, OPENMS_GET_TEST_DATA_PATH("MetaboliteFeatureDeconvolution_ppm_output.consensusXML"));

END_SECTION

START_SECTION(([EXTRA] Multimer detection with max_multimer=2))
{
  // Cross-adduct multimer detection test.
  // MassExplainer needs at least 2 charged adducts to produce valid net_charge=0
  // compomers (e.g. Na+LEFT / H+RIGHT). These compomers are required for
  // queryMultimer() to match feature pairs with different multipliers.
  //
  // Setup: M = 500.0 (neutral mass)
  // Feature 1: monomer [M+Na]+  -> m/z ~ M + Na_adduct_mass ~ 522.989
  // Feature 2: dimer   [2M+H]+  -> m/z ~ 2*M + H_adduct_mass ~ 1001.007
  //
  // The multimer equation: n2*m1 - n1*m2 = n2*left_mass - n1*right_mass
  // With n1=1 (monomer), n2=2 (dimer):
  //   2*(M+Na) - 1*(2M+H) = 2*Na - H  (M cancels exactly)
  // This matches the Na+LEFT/H+RIGHT compomer.

  // Approximate adduct masses (H+ ~ proton mass, Na+ ~ Na - electron)
  double H_adduct = 1.00728;
  double Na_adduct = 22.98922;
  double M = 500.0;

  FeatureMap fm;
  Feature f1;
  f1.setMZ(M + Na_adduct);   // monomer with Na+ adduct
  f1.setRT(100.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(1);
  fm.push_back(f1);

  Feature f2;
  f2.setMZ(2.0 * M + H_adduct);  // dimer with H+ adduct
  f2.setRT(100.0);
  f2.setIntensity(500.0f);
  f2.setCharge(1);
  fm.push_back(f2);

  // Run decharger with multimer detection enabled
  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 3);       // >= 2 needed for valid compomers
  p.setValue("charge_span_max", 3);
  p.setValue("max_multimer", 2);
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:0.7","Na:+:0.3"});
  p.setValue("mass_max_diff", 0.5);  // generous tolerance for approximate masses
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // Check that exactly one consensus feature groups the two features
  Size group_count = 0;
  Size group_size = 0;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 2)
    {
      ++group_count;
      group_size = cm_out[i].size();
    }
  }
  TEST_EQUAL(group_count, 1);
  TEST_EQUAL(group_size, 2);

  // Check annotation: one feature should have mol_multiplier=2
  bool found_dimer_annotation = false;
  for (Size i = 0; i < fm_out.size(); ++i)
  {
    if (fm_out[i].metaValueExists("mol_multiplier"))
    {
      TEST_EQUAL((Int)fm_out[i].getMetaValue("mol_multiplier"), 2);
      found_dimer_annotation = true;
    }
  }
  TEST_EQUAL(found_dimer_annotation, true);

  // Regression: max_multimer=1 should NOT find multimer relationships
  MetaboliteFeatureDeconvolution mfd2;
  Param p2 = mfd2.getDefaults();
  p2.setValue("charge_min", 1);
  p2.setValue("charge_max", 3);
  p2.setValue("charge_span_max", 3);
  p2.setValue("max_multimer", 1);
  p2.setValue("potential_adducts", std::vector<std::string>{"H:+:0.7","Na:+:0.3"});
  p2.setValue("mass_max_diff", 0.5);
  p2.setValue("retention_max_diff", 2.0);
  p2.setValue("retention_max_diff_local", 2.0);
  mfd2.setParameters(p2);

  FeatureMap fm_out2;
  ConsensusMap cm_out2_a, cm_out2_b;
  mfd2.compute(fm, fm_out2, cm_out2_a, cm_out2_b);

  // With max_multimer=1, cross-multiplier edges are never generated
  bool found_group_no_multimer = false;
  for (Size i = 0; i < cm_out2_a.size(); ++i)
  {
    if (cm_out2_a[i].size() >= 2)
    {
      found_group_no_multimer = true;
    }
  }
  TEST_EQUAL(found_group_no_multimer, false);
}
END_SECTION

START_SECTION(([EXTRA] Same-adduct multimer detection))
{
  // Same-adduct multimer: [M+H]+ (monomer) and [2M+H]+ (dimer)
  // This requires identity compomers (H+LEFT / H+RIGHT) in MassExplainer.
  // With identity, the equation gives: 2*H - H = H, which matches.
  double H_adduct = 1.00728;   // approximate proton mass
  double M = 500.0;

  FeatureMap fm;
  Feature f1;
  f1.setMZ(M + H_adduct);           // monomer [M+H]+
  f1.setRT(100.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(1);
  fm.push_back(f1);

  Feature f2;
  f2.setMZ(2.0 * M + H_adduct);     // dimer [2M+H]+
  f2.setRT(100.0);
  f2.setIntensity(500.0f);
  f2.setCharge(1);
  fm.push_back(f2);

  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 1);
  p.setValue("charge_span_max", 1);
  p.setValue("max_multimer", 2);
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:1.0"});
  p.setValue("mass_max_diff", 0.5);
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // The two features should be grouped (exactly one group of size 2)
  Size group_count = 0;
  Size group_size = 0;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 2)
    {
      ++group_count;
      group_size = cm_out[i].size();
    }
  }
  TEST_EQUAL(group_count, 1);
  TEST_EQUAL(group_size, 2);

  // The dimer feature should have mol_multiplier=2
  bool found_dimer = false;
  for (Size i = 0; i < fm_out.size(); ++i)
  {
    if (fm_out[i].metaValueExists("mol_multiplier"))
    {
      TEST_EQUAL((Int)fm_out[i].getMetaValue("mol_multiplier"), 2);
      found_dimer = true;
    }
  }
  TEST_EQUAL(found_dimer, true);
}
END_SECTION

START_SECTION(([EXTRA] Multimer annotation strings and trimer detection))
{
  // Test that annotation strings are correct and trimers (n=3) work.
  // Three features: monomer, dimer, trimer — all with H+ adduct.
  double H_adduct = 1.00728;
  double M = 300.0;

  FeatureMap fm;
  Feature f1;
  f1.setMZ(M + H_adduct);
  f1.setRT(100.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(1);
  fm.push_back(f1);

  Feature f2;
  f2.setMZ(2.0 * M + H_adduct);
  f2.setRT(100.0);
  f2.setIntensity(500.0f);
  f2.setCharge(1);
  fm.push_back(f2);

  Feature f3;
  f3.setMZ(3.0 * M + H_adduct);
  f3.setRT(100.0);
  f3.setIntensity(300.0f);
  f3.setCharge(1);
  fm.push_back(f3);

  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 1);
  p.setValue("charge_span_max", 1);
  p.setValue("max_multimer", 3);
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:1.0"});
  p.setValue("mass_max_diff", 0.5);
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // All three features should be in one group
  Size group_count = 0;
  Size group_size = 0;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 2)
    {
      ++group_count;
      group_size = cm_out[i].size();
    }
  }
  TEST_EQUAL(group_count, 1);
  TEST_EQUAL(group_size, 3);

  // Check annotation content
  bool found_dimer_annotation = false;
  bool found_trimer_annotation = false;
  for (Size i = 0; i < fm_out.size(); ++i)
  {
    if (fm_out[i].metaValueExists("mol_multiplier"))
    {
      Int mult = (Int)fm_out[i].getMetaValue("mol_multiplier");
      if (mult == 2)
      {
        found_dimer_annotation = true;
        TEST_EQUAL(fm_out[i].metaValueExists("adducts"), true);
        StringList adducts = fm_out[i].getMetaValue("adducts");
        TEST_EQUAL(adducts[0].hasSubstring("2M"), true);
      }
      if (mult == 3)
      {
        found_trimer_annotation = true;
        TEST_EQUAL(fm_out[i].metaValueExists("adducts"), true);
        StringList adducts = fm_out[i].getMetaValue("adducts");
        TEST_EQUAL(adducts[0].hasSubstring("3M"), true);
      }
    }
  }
  TEST_EQUAL(found_dimer_annotation, true);
  TEST_EQUAL(found_trimer_annotation, true);
}
END_SECTION

START_SECTION(([EXTRA] Multimer penalty suppresses dimer when penalty is extreme))
{
  // Create features where multimer detection is possible.
  // With extreme penalty (e.g., -100), the ILP should not select multimer edges.
  double H_adduct = 1.00728;
  double M = 500.0;

  FeatureMap fm;
  Feature f1;
  f1.setMZ(M + H_adduct);
  f1.setRT(100.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(1);
  fm.push_back(f1);

  Feature f2;
  f2.setMZ(2.0 * M + H_adduct);
  f2.setRT(100.0);
  f2.setIntensity(500.0f);
  f2.setCharge(1);
  fm.push_back(f2);

  // Extreme penalty: exp(-100) ~ 3.7e-44, effectively zero edge score
  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 1);
  p.setValue("charge_span_max", 1);
  p.setValue("max_multimer", 2);
  p.setValue("multimer_log_penalty", -100.0);
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:1.0"});
  p.setValue("mass_max_diff", 0.5);
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // With extreme penalty, no grouping should occur
  bool found_group = false;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 2) found_group = true;
  }
  TEST_EQUAL(found_group, false);

  // No mol_multiplier annotation should exist
  bool found_multiplier = false;
  for (Size i = 0; i < fm_out.size(); ++i)
  {
    if (fm_out[i].metaValueExists("mol_multiplier")) found_multiplier = true;
  }
  TEST_EQUAL(found_multiplier, false);
}
END_SECTION

START_SECTION(([EXTRA] ILP prevents contradictory multiplier assignments))
{
  // Regression test for ILP variant key bug: without mol_multiplier in the
  // variant key, a feature can be assigned contradictory multipliers from
  // different edges simultaneously.
  //
  // Setup: three features where the middle one (B) participates in two
  // multimer edges with different multipliers:
  //   Edge (A,B): A=monomer(mult=1), B=dimer(mult=2) of A
  //   Edge (B,C): B=monomer(mult=1), C=dimer(mult=2) of B
  //
  // Feature B cannot be BOTH a dimer (mult=2) and a monomer (mult=1).
  // The ILP must choose one interpretation, so at most one edge is active.
  //
  // Feature A: m/z = 100 + H  (M_A = 100)
  // Feature B: m/z = 200 + H  (= 2*M_A, or M_B = 200)
  // Feature C: m/z = 400 + H  (= 2*M_B)

  double H_adduct = 1.00728;

  FeatureMap fm;
  Feature fA;
  fA.setMZ(100.0 + H_adduct);
  fA.setRT(100.0);
  fA.setIntensity(1000.0f);
  fA.setCharge(1);
  fm.push_back(fA);

  Feature fB;
  fB.setMZ(200.0 + H_adduct);
  fB.setRT(100.0);
  fB.setIntensity(800.0f);
  fB.setCharge(1);
  fm.push_back(fB);

  Feature fC;
  fC.setMZ(400.0 + H_adduct);
  fC.setRT(100.0);
  fC.setIntensity(500.0f);
  fC.setCharge(1);
  fm.push_back(fC);

  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 1);
  p.setValue("charge_span_max", 1);
  p.setValue("max_multimer", 2);
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:1.0"});
  p.setValue("mass_max_diff", 0.5);
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // With correct ILP modeling, Feature B must be assigned a single multiplier.
  // The ILP can activate at most one of the two edges, so no group of size 3
  // should exist. A group of 3 indicates the bug: both edges active, B has
  // contradictory multiplier assignments.
  bool found_group_of_3 = false;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 3)
    {
      found_group_of_3 = true;
    }
  }
  TEST_EQUAL(found_group_of_3, false);

  // Exactly one group of size 2 should exist (either A-B or B-C, not both)
  Size groups_of_2 = 0;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() == 2) ++groups_of_2;
  }
  TEST_EQUAL(groups_of_2, 1);
}
END_SECTION

START_SECTION(([EXTRA] Negative mode multimer detection))
{
  // Negative mode: [M-H]- (monomer) and [2M-H]- (dimer)
  // H-1 adduct has charge -1 and mass = -PROTON_MASS_U
  double H_adduct = 1.00728;  // proton mass (will be subtracted in negative mode)
  double M = 500.0;

  FeatureMap fm;
  Feature f1;
  f1.setMZ(M - H_adduct);           // monomer [M-H]- (deprotonated)
  f1.setRT(100.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(1);                   // FFM reports absolute charge
  fm.push_back(f1);

  Feature f2;
  f2.setMZ(2.0 * M - H_adduct);     // dimer [2M-H]-
  f2.setRT(100.0);
  f2.setIntensity(500.0f);
  f2.setCharge(1);
  fm.push_back(f2);

  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", -1);
  p.setValue("charge_max", -1);
  p.setValue("charge_span_max", 1);
  p.setValue("max_multimer", 2);
  p.setValue("negative_mode", "true");
  p.setValue("potential_adducts", std::vector<std::string>{"H-1:-:1.0"});
  p.setValue("mass_max_diff", 0.5);
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // The two features should be grouped as monomer/dimer (one group of size 2)
  Size group_count = 0;
  Size group_size = 0;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 2)
    {
      ++group_count;
      group_size = cm_out[i].size();
    }
  }
  TEST_EQUAL(group_count, 1);
  TEST_EQUAL(group_size, 2);

  // Dimer annotation should exist
  bool found_dimer = false;
  for (Size i = 0; i < fm_out.size(); ++i)
  {
    if (fm_out[i].metaValueExists("mol_multiplier"))
    {
      TEST_EQUAL((Int)fm_out[i].getMetaValue("mol_multiplier"), 2);
      found_dimer = true;
    }
  }
  TEST_EQUAL(found_dimer, true);
}
END_SECTION

START_SECTION(([EXTRA] Mixed-charge multimer detection))
{
  // Mixed charge: [M+H]+ (charge 1) and [2M+2H]2+ (charge 2)
  // Both use H+ adduct but at different charges.
  // mz1 = M + H = 501.007, q1 = 1, m1 = 501.007
  // mz2 = (2M + 2H) / 2 = M + H = 501.007, q2 = 2, m2 = 1002.014
  // For n1=1, n2=2: observed = 2*501.007 - 1*1002.014 = 0.0
  // Requires H1_LEFT/H2_RIGHT compomer: expected = 2*H - 1*(2*H) = 0. Match!
  // This compomer is generated by the same-adduct enumeration in compute()
  // (the i!=j case). The standard enumeration never puts the same adduct on
  // both LEFT and RIGHT sides, so this relies on the multimer compomer loop.
  //
  // NOTE: mass_max_diff must be tight (0.01 Da) so that the wrong compomer
  // H1_RIGHT (expected = -H ≈ -1.007, residual ≈ 1.007) does NOT match.
  // Only the correct H1_LEFT/H2_RIGHT compomer (residual = 0) should pass.

  double H_adduct = 1.00728;
  double M = 500.0;

  FeatureMap fm;
  Feature f1;
  f1.setMZ(M + H_adduct);                     // [M+H]+, charge 1
  f1.setRT(100.0);
  f1.setIntensity(1000.0f);
  f1.setCharge(1);
  fm.push_back(f1);

  Feature f2;
  f2.setMZ((2.0 * M + 2.0 * H_adduct) / 2.0);  // [2M+2H]2+, charge 2
  f2.setRT(100.0);
  f2.setIntensity(500.0f);
  f2.setCharge(2);
  fm.push_back(f2);

  MetaboliteFeatureDeconvolution mfd;
  Param p = mfd.getDefaults();
  p.setValue("charge_min", 1);
  p.setValue("charge_max", 2);
  p.setValue("charge_span_max", 2);
  p.setValue("max_multimer", 2);
  p.setValue("potential_adducts", std::vector<std::string>{"H:+:1.0"});
  p.setValue("mass_max_diff", 0.01);
  p.setValue("retention_max_diff", 2.0);
  p.setValue("retention_max_diff_local", 2.0);
  mfd.setParameters(p);

  FeatureMap fm_out;
  ConsensusMap cm_out, cm_out2;
  mfd.compute(fm, fm_out, cm_out, cm_out2);

  // Features should be grouped (one group of size 2)
  Size group_count = 0;
  Size group_size = 0;
  for (Size i = 0; i < cm_out.size(); ++i)
  {
    if (cm_out[i].size() >= 2)
    {
      ++group_count;
      group_size = cm_out[i].size();
    }
  }
  TEST_EQUAL(group_count, 1);
  TEST_EQUAL(group_size, 2);

  // Dimer annotation should exist
  bool found_dimer = false;
  for (Size i = 0; i < fm_out.size(); ++i)
  {
    if (fm_out[i].metaValueExists("mol_multiplier"))
    {
      TEST_EQUAL((Int)fm_out[i].getMetaValue("mol_multiplier"), 2);
      found_dimer = true;
    }
  }
  TEST_EQUAL(found_dimer, true);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
