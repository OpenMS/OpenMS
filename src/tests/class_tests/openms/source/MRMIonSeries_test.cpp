// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/MRMIonSeries.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MRMIonSeries, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MRMIonSeries * ptr = nullptr;
MRMIonSeries* nullPointer = nullptr;

START_SECTION(MRMIonSeries())
{
  ptr = new MRMIonSeries();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MRMIonSeries())
{
  delete ptr;
}

END_SECTION

START_SECTION((boost::unordered_map<String, double> MRMIonSeries::getIonSeries(AASequence sequence, size_t precursor_charge, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses)))
{
  MRMIonSeries mrmis;
  std::vector<String> fragment_types;
  fragment_types.push_back(String("b"));
  fragment_types.push_back(String("y"));

  std::vector<size_t> fragment_charges;
  fragment_charges.push_back(3);
  fragment_charges.push_back(2);
  fragment_charges.push_back(1);

  // Standard peptide
  MRMIonSeries::IonSeries ionseries1 = mrmis.getIonSeries(AASequence::fromString(String("PEPTIDEK")), 3, fragment_types, fragment_charges, false, false);

  TEST_EQUAL(ionseries1.size(), 42)
  TEST_REAL_SIMILAR(ionseries1["b2^1"], 227.10263491)
  TEST_REAL_SIMILAR(ionseries1["b2^2"], 114.05495569)
  TEST_REAL_SIMILAR(ionseries1["b2^3"], 76.37239595)
  TEST_REAL_SIMILAR(ionseries1["b5^1"], 538.28714271)
  TEST_REAL_SIMILAR(ionseries1["b5^2"], 269.64720959)
  TEST_REAL_SIMILAR(ionseries1["b5^3"], 180.10056521)
  TEST_REAL_SIMILAR(ionseries1["b7^1"], 782.35668109)
  TEST_REAL_SIMILAR(ionseries1["b7^2"], 391.68197878)
  TEST_REAL_SIMILAR(ionseries1["b7^3"], 261.45707801)
  TEST_REAL_SIMILAR(ionseries1["y1^1"], 147.11280491)
  TEST_REAL_SIMILAR(ionseries1["y1^2"], 74.06004069)
  TEST_REAL_SIMILAR(ionseries1["y1^3"], 49.70911928)
  TEST_REAL_SIMILAR(ionseries1["y5^1"], 605.31408687)
  TEST_REAL_SIMILAR(ionseries1["y5^2"], 303.16068167)
  TEST_REAL_SIMILAR(ionseries1["y5^3"], 202.44287993)
  TEST_REAL_SIMILAR(ionseries1["y7^1"], 831.40944532)
  TEST_REAL_SIMILAR(ionseries1["y7^2"], 416.20836089)
  TEST_REAL_SIMILAR(ionseries1["y7^3"], 277.80799942)

  // Enable neutral losses
  MRMIonSeries::IonSeries ionseries2 = mrmis.getIonSeries(AASequence::fromString(String("PEPTIDEK")), 3, fragment_types, fragment_charges, true, true);

  TEST_EQUAL(ionseries2.size(), 99)
  TEST_REAL_SIMILAR(ionseries2["b5^2"], 269.64720959)
  TEST_REAL_SIMILAR(ionseries2["b5-H2O1^2"], 260.64192709)
  TEST_REAL_SIMILAR(ionseries2["b5-C1H4O1S1^2"], 0) // no oxidized methionine in peptide
  TEST_REAL_SIMILAR(ionseries2["b5-H1O3P1^2"], 0) // no phosphorylation in peptide
  TEST_REAL_SIMILAR(ionseries2["b5-H3O4P1^2"], 0) // no phosphorylation in peptide
  TEST_REAL_SIMILAR(ionseries2["b5-C1O2^2"], 0)


  MRMIonSeries::IonSeries ionseries3 = mrmis.getIonSeries(AASequence::fromString(String("ES(Phospho)")), 3, fragment_types, fragment_charges, true, true);

  TEST_EQUAL(ionseries3.size(), 12)

  TEST_REAL_SIMILAR(ionseries3["y1^1"], 186.0162)
  TEST_REAL_SIMILAR(ionseries3["y1-H3O4P1^2"], 44.5233)
  TEST_REAL_SIMILAR(ionseries3["y1^2"], 93.5117)
  TEST_REAL_SIMILAR(ionseries3["y1^3"], 62.6769)
  TEST_REAL_SIMILAR(ionseries3["b1-H2O1^1"], 112.0393)
  TEST_REAL_SIMILAR(ionseries3["b1-H2O1^2"], 56.5233)
  TEST_REAL_SIMILAR(ionseries3["y1-H3O4P1^1"], 88.0393)
  TEST_REAL_SIMILAR(ionseries3["y1-H3O4P1^3"], 30.018)
  TEST_REAL_SIMILAR(ionseries3["b1^2"], 65.5285)
  TEST_REAL_SIMILAR(ionseries3["b1-H2O1^3"], 38.018)
  TEST_REAL_SIMILAR(ionseries3["b1^1"], 130.0498)
  TEST_REAL_SIMILAR(ionseries3["b1^3"], 44.0214)

  MRMIonSeries::IonSeries ionseries4 = mrmis.getIonSeries(AASequence::fromString(String("ES")), 3, fragment_types, fragment_charges, true, true);

  TEST_REAL_SIMILAR(ionseries4["y1-H2O1^1"], 88.0393)
  TEST_REAL_SIMILAR(ionseries4["y1^1"], 106.0499)
  TEST_REAL_SIMILAR(ionseries4["y1-H2O1^2"], 44.5233)
  TEST_REAL_SIMILAR(ionseries4["y1^2"], 53.5286)
  TEST_REAL_SIMILAR(ionseries4["y1-H2O1^3"], 30.0180)
  TEST_REAL_SIMILAR(ionseries4["y1^3"], 36.0215)
  TEST_REAL_SIMILAR(ionseries4["b1-H2O1^1"], 112.0393)
  TEST_REAL_SIMILAR(ionseries4["b1-H2O1^2"], 56.5233)
  TEST_REAL_SIMILAR(ionseries4["b1^2"], 65.5286)
  TEST_REAL_SIMILAR(ionseries4["b1-H2O1^3"], 38.0180)
  TEST_REAL_SIMILAR(ionseries4["b1^1"], 130.0499)
  TEST_REAL_SIMILAR(ionseries4["b1^3"], 44.0215)

  MRMIonSeries::IonSeries ionseries5 = mrmis.getIonSeries(AASequence::fromString(String("ES(Phospho)")), 3, fragment_types, fragment_charges, true, false);

  TEST_EQUAL(ionseries5.size(), 9)

  TEST_REAL_SIMILAR(ionseries5["y1^1"], 186.0162)
  TEST_REAL_SIMILAR(ionseries5["y1-H3O4P1^2"], 44.5233)
  TEST_REAL_SIMILAR(ionseries5["y1^2"], 93.5117)
  TEST_REAL_SIMILAR(ionseries5["y1^3"], 62.6769)
  TEST_REAL_SIMILAR(ionseries5["y1-H3O4P1^1"], 88.0393)
  TEST_REAL_SIMILAR(ionseries5["y1-H3O4P1^3"], 30.018)
  TEST_REAL_SIMILAR(ionseries5["b1^2"], 65.5285)
  TEST_REAL_SIMILAR(ionseries5["b1^1"], 130.0498)
  TEST_REAL_SIMILAR(ionseries5["b1^3"], 44.0214)

  MRMIonSeries::IonSeries ionseries6 = mrmis.getIonSeries(AASequence::fromString(String("ES(Phospho)")), 3, fragment_types, fragment_charges, false, true);

  TEST_EQUAL(ionseries6.size(), 9)

  TEST_REAL_SIMILAR(ionseries6["y1^1"], 186.0162)
  TEST_REAL_SIMILAR(ionseries6["y1^2"], 93.5117)
  TEST_REAL_SIMILAR(ionseries6["y1^3"], 62.6769)
  TEST_REAL_SIMILAR(ionseries6["b1-H2O1^1"], 112.0393)
  TEST_REAL_SIMILAR(ionseries6["b1-H2O1^2"], 56.5233)
  TEST_REAL_SIMILAR(ionseries6["b1^2"], 65.5285)
  TEST_REAL_SIMILAR(ionseries6["b1-H2O1^3"], 38.018)
  TEST_REAL_SIMILAR(ionseries6["b1^1"], 130.0498)
  TEST_REAL_SIMILAR(ionseries6["b1^3"], 44.0214)

}

END_SECTION

START_SECTION((std::pair<String, double> MRMIonSeries::annotateIon(IonSeries ionseries, double ProductMZ, double mz_threshold)))
{
  MRMIonSeries mrmis;
  std::vector<String> fragment_types;
  fragment_types.push_back(String("b"));
  fragment_types.push_back(String("y"));

  std::vector<size_t> fragment_charges;
  fragment_charges.push_back(3);
  fragment_charges.push_back(2);
  fragment_charges.push_back(1);

  // Standard peptide
  MRMIonSeries::IonSeries ionseries1 = mrmis.getIonSeries(AASequence::fromString(String("PEPTIDEK")), 3, fragment_types, fragment_charges, false, false);

  std::pair<String, double> ion1 = mrmis.annotateIon(ionseries1, 202.44287993, 0.05);
  TEST_EQUAL(ion1.first, "y5^3")
  TEST_REAL_SIMILAR(ion1.second, 202.44287993)

  std::pair<String, double> ion2 = mrmis.annotateIon(ionseries1, 202.44287993, 0);
  TEST_EQUAL(ion2.first, "unannotated")
  TEST_REAL_SIMILAR(ion2.second, -1)

  std::pair<String, double> ion3 = mrmis.annotateIon(ionseries1, 202.4, 0.05);
  TEST_EQUAL(ion3.first, "y5^3")
  TEST_REAL_SIMILAR(ion3.second, 202.44287993)
}

END_SECTION

START_SECTION((std::pair<String, double> MRMIonSeries::getIon(IonSeries ionseries, String ionid)))
{
  MRMIonSeries mrmis;
  std::vector<String> fragment_types;
  fragment_types.push_back(String("b"));
  fragment_types.push_back(String("y"));

  std::vector<size_t> fragment_charges;
  fragment_charges.push_back(3);
  fragment_charges.push_back(2);
  fragment_charges.push_back(1);

  // Standard peptide
  MRMIonSeries::IonSeries ionseries1 = mrmis.getIonSeries(AASequence::fromString(String("PEPTIDEK")), 3, fragment_types, fragment_charges, false, false);

  std::pair<String, double> ion1 = mrmis.getIon(ionseries1, "y5^3");
  TEST_EQUAL(ion1.first, "y5^3")
  TEST_REAL_SIMILAR(ion1.second, 202.44287993)
}

END_SECTION

START_SECTION((void MRMIonSeries::annotateTransitionCV(ReactionMonitoringTransition & tr, String annotation)))
{
  MRMIonSeries mrmis;
  ReactionMonitoringTransition tr, tr2, tr3;

  mrmis.annotateTransitionCV(tr, "y5^3");
  mrmis.annotateTransitionCV(tr2, "y5-H2O1^3");
  mrmis.annotateTransitionCV(tr3, "y5-18^3");

  TEST_EQUAL(tr.getProduct().getChargeState(), 3)
  TEST_EQUAL(tr.getProduct().getInterpretationList()[0].iontype, TargetedExperiment::IonType::YIon);
  TEST_EQUAL(tr.getProduct().getInterpretationList()[0].ordinal, 5);
  TEST_EQUAL(tr.getProduct().getInterpretationList()[0].hasCVTerm("MS:1001524"), false)  // no neutral loss

  TEST_EQUAL(tr2.getProduct().getChargeState(), 3)
  TEST_EQUAL(tr2.getProduct().getInterpretationList()[0].iontype, TargetedExperiment::IonType::YIon);
  TEST_EQUAL(tr2.getProduct().getInterpretationList()[0].ordinal, 5)
  TEST_EQUAL(tr2.getProduct().getInterpretationList()[0].hasCVTerm("MS:1001524"), true)
  TEST_REAL_SIMILAR(tr2.getProduct().getInterpretationList()[0].getCVTerms().at("MS:1001524")[0].getValue(), -18.0105650638)

  TEST_EQUAL(tr3.getProduct().getChargeState(), 3)
  TEST_EQUAL(tr3.getProduct().getInterpretationList()[0].iontype, TargetedExperiment::IonType::YIon);
  TEST_EQUAL(tr3.getProduct().getInterpretationList()[0].ordinal, 5);
  TEST_EQUAL(tr3.getProduct().getInterpretationList()[0].hasCVTerm("MS:1001524"), true)
  TEST_REAL_SIMILAR(tr3.getProduct().getInterpretationList()[0].getCVTerms().at("MS:1001524")[0].getValue(), -18)
}

END_SECTION

START_SECTION((void MRMIonSeries::annotateTransition(ReactionMonitoringTransition & tr, const TargetedExperiment::Peptide peptide, const double precursor_mz_threshold, const double product_mz_threshold, bool enable_reannotation, std::vector<String> fragment_types, std::vector<size_t> fragment_charges, bool enable_specific_losses, bool enable_unspecific_losses)))
{
  MRMIonSeries mrmis;
  ReactionMonitoringTransition tr, tr2, tr3;
  TargetedExperiment::Peptide peptide;
  peptide.sequence = "PEPTIDEK";
  peptide.setChargeState(3);

  std::vector<String> fragment_types;
  fragment_types.push_back(String("b"));
  fragment_types.push_back(String("y"));

  std::vector<size_t> fragment_charges;
  fragment_charges.push_back(3);
  fragment_charges.push_back(2);
  fragment_charges.push_back(1);

  tr.setProductMZ(202.44);
  mrmis.annotateTransition(tr, peptide, 0.05, 0.05, true, fragment_types, fragment_charges, false, false);

  TEST_REAL_SIMILAR(tr.getProductMZ(), 202.442879934638)
  TEST_EQUAL(tr.getProduct().getChargeState(), 3)
  TEST_EQUAL(tr.getProduct().getInterpretationList()[0].iontype, TargetedExperiment::IonType::YIon);
  TEST_EQUAL(tr.getProduct().getInterpretationList()[0].ordinal, 5)
  TEST_EQUAL(tr.getProduct().getInterpretationList()[0].hasCVTerm("MS:1001524"), false) // no neutral loss

  tr2.setProductMZ(196.44287993);
  mrmis.annotateTransition(tr2, peptide, 0.05, 0.05, true, fragment_types, fragment_charges, true, true);

  TEST_EQUAL(tr2.getProduct().getChargeState(), 3)
  TEST_EQUAL(tr2.getProduct().getInterpretationList()[0].iontype, TargetedExperiment::IonType::YIon);
  TEST_EQUAL(tr2.getProduct().getInterpretationList()[0].ordinal, 5)
  TEST_EQUAL(tr2.getProduct().getInterpretationList()[0].hasCVTerm("MS:1001524"), true)
  TEST_EQUAL((int)tr2.getProduct().getInterpretationList()[0].getCVTerms().at("MS:1001524")[0].getValue().toString().toDouble(), -18)

  tr3.setProductMZ(202.44);
  mrmis.annotateTransition(tr3, peptide, 0.05, 0.05, false, fragment_types, fragment_charges, false, false);

  TEST_REAL_SIMILAR(tr3.getProductMZ(), 202.44)
  TEST_EQUAL(tr3.getProduct().hasCharge(), false)
  TEST_EQUAL(tr3.getProduct().getInterpretationList()[0].iontype, TargetedExperiment::IonType::NonIdentified);
  TEST_EQUAL(tr3.getProduct().getInterpretationList()[0].hasCVTerm("MS:1001524"), false)
}

END_SECTION

END_TEST
