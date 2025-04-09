// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/ANALYSIS/OPENSWATH/ConfidenceScoring.h>

///////////////////////////

using namespace OpenMS;

std::vector<TargetedExperiment::RetentionTime> get_rts_(double rt_val)
{
  // add retention time for the peptide
  std::vector<TargetedExperiment::RetentionTime> retention_times;
  TargetedExperiment::RetentionTime retention_time;
  retention_time.setRT(rt_val);
  retention_time.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::NORMALIZED;
  retention_times.push_back(retention_time);
  return retention_times;
}

START_TEST(ConfidenceScoring<D>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ConfidenceScoring* confidence_scoring_ptr = nullptr;
ConfidenceScoring* confidence_scoring_nullPointer = nullptr;

START_SECTION((explicit ConfidenceScoring(bool test_mode_=false)))
  confidence_scoring_ptr = new ConfidenceScoring;
  TEST_NOT_EQUAL(confidence_scoring_ptr, confidence_scoring_nullPointer)
END_SECTION

START_SECTION((virtual ~ConfidenceScoring()))
    delete confidence_scoring_ptr;
END_SECTION


START_SECTION((void initialize(TargetedExperiment library, Size n_decoys, Size n_transitions, TransformationDescription rt_trafo)))
  ConfidenceScoring scoring;
  TargetedExperiment library;
  TransformationDescription rt_trafo;
  scoring.initialize(library, 0, 0, rt_trafo);
  TEST_NOT_EQUAL(&scoring, confidence_scoring_nullPointer)
END_SECTION

START_SECTION((void initializeGlm(double intercept, double rt_coef, double int_coef)))
  ConfidenceScoring scoring;
  scoring.initializeGlm(0.0, -1.0, -1.0);
  TEST_NOT_EQUAL(&scoring, confidence_scoring_nullPointer)
END_SECTION

START_SECTION((void scoreMap(FeatureMap & features)))
{
  ConfidenceScoring scoring(true); // initialize with test mode
  TargetedExperiment library;
  TransformationDescription rt_trafo;
  scoring.initialize(library, 0, 0, rt_trafo);
  scoring.initializeGlm(0.0, -1.0, -1.0);
  FeatureMap features;
  TEST_EXCEPTION(Exception::IllegalArgument, scoring.scoreMap(features))

  // The input to the program is 
  // - a transition library which contains peptides with corresponding assays
  // - a feature map where each feature corresponds to an assay (mapped with
  //   MetaValue "PeptideRef") and each feature has as many subordinates as the
  //   assay has transitions (mapped with MetaValue "native_id")

  // In this case we have 2 assays (pep_1 and pep_2) with 1 transition each
  // (tr_10 for pep_1 and tr_20 for pep_2).
  {
    TargetedExperiment::Peptide p;

    p.id = "pep_1";
    p.rts = get_rts_(50.0);
    library.addPeptide(p);

    ReactionMonitoringTransition rm_trans;
    rm_trans.setNativeID("tr_10");
    rm_trans.setPrecursorMZ(400.0);
    rm_trans.setProductMZ(500.0);
    rm_trans.setPeptideRef(p.id);
    rm_trans.setLibraryIntensity(500.0);
    library.addTransition(rm_trans);
  }
  {
    TargetedExperiment::Peptide p;
    p.id = "pep_2";
    p.rts = get_rts_(60.0);
    library.addPeptide(p);


    ReactionMonitoringTransition rm_trans;
    rm_trans.setNativeID("tr_20");
    rm_trans.setPrecursorMZ(400.0);
    rm_trans.setProductMZ(500.0);
    rm_trans.setPeptideRef(p.id);
    rm_trans.setLibraryIntensity(500.0);
    library.addTransition(rm_trans);

  }

  {
    Feature f;
    f.setRT(60.0);
    f.setMetaValue("PeptideRef", "pep_1");
    f.setOverallQuality(-1);

    std::vector<Feature> subordinates;
    Feature sub;
    sub.setIntensity(1);
    sub.setMZ(500);
    sub.setMetaValue("native_id", "tr_10");
    subordinates.push_back(sub);
    f.setSubordinates(subordinates);

    features.push_back(f);
  }
  {
    Feature f;
    f.setRT(60.0);
    f.setMetaValue("PeptideRef", "pep_2");
    f.setOverallQuality(-1);

    std::vector<Feature> subordinates;
    Feature sub;
    sub.setIntensity(1);
    sub.setMZ(500);
    sub.setMetaValue("native_id", "tr_20");
    subordinates.push_back(sub);
    f.setSubordinates(subordinates);

    features.push_back(f);
  }

  scoring.initialize(library, 0, 0, rt_trafo);
  scoring.scoreMap(features);

  TEST_REAL_SIMILAR(features[0].getOverallQuality(), 0.0);
  TEST_REAL_SIMILAR(features[1].getOverallQuality(), 1.0);

  // the absolute computed score for each feature
  TEST_REAL_SIMILAR(features[0].getMetaValue("GLM_score"), 0.0);
  TEST_REAL_SIMILAR(features[1].getMetaValue("GLM_score"), 0.5);
  // the local fdr score (1-quality)
  TEST_REAL_SIMILAR(features[0].getMetaValue("local_FDR"), 1.0);
  TEST_REAL_SIMILAR(features[1].getMetaValue("local_FDR"), 0.0);
}
END_SECTION


START_SECTION(([EXTRA] test exceptions))
{
  ConfidenceScoring scoring(true); // initialize with test mode
  TargetedExperiment library;
  TransformationDescription rt_trafo;
  scoring.initialize(library, 0, 0, rt_trafo);
  scoring.initializeGlm(0.0, -1.0, -1.0);
  FeatureMap features;

  {
    TargetedExperiment::Peptide p;

    p.id = "pep_1";
    p.rts = get_rts_(50.0);
    library.addPeptide(p);

    ReactionMonitoringTransition rm_trans;
    rm_trans.setNativeID("tr_10");
    rm_trans.setPrecursorMZ(400.0);
    rm_trans.setProductMZ(500.0);
    rm_trans.setPeptideRef(p.id);
    rm_trans.setLibraryIntensity(500.0);
    library.addTransition(rm_trans);
  }
  {
    TargetedExperiment::Peptide p;
    p.id = "pep_2";
    p.rts = get_rts_(60.0);
    library.addPeptide(p);

    ReactionMonitoringTransition rm_trans;
    rm_trans.setNativeID("tr_20");
    rm_trans.setPrecursorMZ(400.0);
    rm_trans.setProductMZ(500.0);
    rm_trans.setPeptideRef(p.id);
    rm_trans.setLibraryIntensity(500.0);
    library.addTransition(rm_trans);

  }

  // If no meta value is present for the featuere, we cannot map it to the assay
  {
    Feature f;
    f.setRT(60.0);
    f.setOverallQuality(-1);
    //f.setMetaValue("PeptideRef", "pep_1");
    features.push_back(f);
  }
  {
    Feature f;
    f.setRT(60.0);
    f.setOverallQuality(-1);
    //f.setMetaValue("PeptideRef", "pep_2");
    features.push_back(f);
  }

  scoring.initialize(library, 0, 0, rt_trafo);
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, scoring.scoreMap(features), "Feature does not contain meta value 'PeptideRef' (reference to assay)")

  // After we add the meta value, we still should get an exception
  features[0].setMetaValue("PeptideRef", "pep_1");
  features[1].setMetaValue("PeptideRef", "pep_2");
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, scoring.scoreMap(features), "Feature intensities were empty - please provide feature subordinate with intensities")

  // An exception should be thrown if the sub-features cannot be mapped to the
  // transitions (e.g. the metavalue "native_id" is missing)
  {
    std::vector<Feature> subordinates;
    Feature sub;
    sub.setIntensity(1);
    sub.setMZ(500);
    // sub.setMetaValue("native_id", "tr_10");
    subordinates.push_back(sub);
    features[0].setSubordinates(subordinates);
  }
  {
    std::vector<Feature> subordinates;
    Feature sub;
    sub.setIntensity(1);
    sub.setMZ(500);
    //sub.setMetaValue("native_id", "tr_20");
    subordinates.push_back(sub);
    features[1].setSubordinates(subordinates);
  }
  TEST_EXCEPTION_WITH_MESSAGE(Exception::IllegalArgument, scoring.scoreMap(features), "Did not find a feature for each assay provided - each feature needs to have n subordinates with the meta-value 'native_id' set to the corresponding transition.")

  {
    std::vector<Feature> subordinates;
    Feature sub;
    sub.setIntensity(1);
    sub.setMZ(500);
    sub.setMetaValue("native_id", "tr_10");
    subordinates.push_back(sub);
    features[0].setSubordinates(subordinates);
  }
  {
    std::vector<Feature> subordinates;
    Feature sub;
    sub.setIntensity(1);
    sub.setMZ(500);
    sub.setMetaValue("native_id", "tr_20");
    subordinates.push_back(sub);
    features[1].setSubordinates(subordinates);
  }
  scoring.scoreMap(features);
  TEST_REAL_SIMILAR(features[0].getOverallQuality(), 0.0);
  TEST_REAL_SIMILAR(features[1].getOverallQuality(), 1.0);


}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
