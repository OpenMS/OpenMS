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
  TEST_EQUAL(cm, MetaboliteFeatureDeconvolution::QFROMFEATURE)
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

  TEST_EQUAL(cm, MetaboliteFeatureDeconvolution::QHEURISTIC)
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
