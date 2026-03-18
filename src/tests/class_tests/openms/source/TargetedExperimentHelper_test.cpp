// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperimentHelper.h>
///////////////////////////

#include <unordered_set>

using namespace OpenMS;
using namespace std;
using namespace OpenMS::TargetedExperimentHelper;

START_TEST(TargetedExperimentHelper, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(TargetedExperimentHelper::Configuration())
{
  // Ensure that Configuration has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Configuration(std::declval<TargetedExperimentHelper::Configuration&&>())), true)

  auto ptr = new TargetedExperimentHelper::Configuration();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::CV())
{
  // Ensure that CV has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::CV(std::declval<TargetedExperimentHelper::CV&&>())), true)

  auto ptr = new TargetedExperimentHelper::CV("", "", "", "");
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Protein())
{
  // Ensure that Protein has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Protein(std::declval<TargetedExperimentHelper::Protein&&>())), true)

  auto ptr = new TargetedExperimentHelper::Protein();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::RetentionTime())
{
  // Ensure that RetentionTime has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::RetentionTime(std::declval<TargetedExperimentHelper::RetentionTime&&>())), true)

  auto ptr = new TargetedExperimentHelper::RetentionTime();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::PeptideCompound())
{
  // Ensure that PeptideCompound has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::PeptideCompound(std::declval<TargetedExperimentHelper::PeptideCompound&&>())), true)

  auto ptr = new TargetedExperimentHelper::PeptideCompound();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Peptide())
{
  // Ensure that Peptide has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Peptide(std::declval<TargetedExperimentHelper::Peptide&&>())), true)

  auto ptr = new TargetedExperimentHelper::Peptide();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Compound())
{
  // Ensure that Compound has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Compound(std::declval<TargetedExperimentHelper::Compound&&>())), true)

  auto ptr = new TargetedExperimentHelper::Compound();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Contact())
{
  // Ensure that Contact has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Contact(std::declval<TargetedExperimentHelper::Contact&&>())), true)

  auto ptr = new TargetedExperimentHelper::Contact();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Publication())
{
  // Ensure that Publication has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Publication(std::declval<TargetedExperimentHelper::Publication&&>())), true)

  auto ptr = new TargetedExperimentHelper::Publication();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Instrument())
{
  // Ensure that Instrument has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Instrument(std::declval<TargetedExperimentHelper::Instrument&&>())), true)

  auto ptr = new TargetedExperimentHelper::Instrument();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Prediction())
{
  // Ensure that Prediction has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Prediction(std::declval<TargetedExperimentHelper::Prediction&&>())), true)

  auto ptr = new TargetedExperimentHelper::Prediction();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::Interpretation())
{
  // Ensure that Interpretation has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::Interpretation(std::declval<TargetedExperimentHelper::Interpretation&&>())), true)

  auto ptr = new TargetedExperimentHelper::Interpretation();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

START_SECTION(TargetedExperimentHelper::TraMLProduct())
{
  // Ensure that TraMLProduct has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(TargetedExperimentHelper::TraMLProduct(std::declval<TargetedExperimentHelper::TraMLProduct&&>())), true)

  auto ptr = new TargetedExperimentHelper::TraMLProduct();
  TEST_FALSE(ptr == nullptr)
  delete ptr;
}
END_SECTION

TargetedExperimentHelper::Peptide* ptr = nullptr;
TargetedExperimentHelper::Peptide* null_ptr = nullptr;

START_SECTION(TargetedExperimentHelper::Peptide())
{
  ptr = new TargetedExperimentHelper::Peptide();
  TEST_NOT_EQUAL(ptr, null_ptr)

  // Ensure that Peptide has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(Peptide(std::declval<Peptide&&>())), true)
}
END_SECTION

START_SECTION(~TargetedExperimentHelper::Peptide())
{
  delete ptr;
}
END_SECTION

START_SECTION((TargetedExperiment::RetentionTime))
{
  TargetedExperimentHelper::RetentionTime rt;
    
  TEST_EQUAL(rt.isRTset(), false)

  rt.setRT(5.0);
  TEST_EQUAL(rt.isRTset(), true)
  TEST_REAL_SIMILAR (rt.getRT(), 5.0)

}
END_SECTION

START_SECTION((TargetedExperiment::Peptide))
{
  TargetedExperimentHelper::Peptide p;
    
  TEST_EQUAL(p.hasRetentionTime(), false)
  TEST_EQUAL(p.rts.size(), 0)

  // add a RT
  TargetedExperimentHelper::RetentionTime rt;
  rt.setRT(5.0);
  rt.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
  rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::PREDICTED;
  p.rts.push_back(rt);

  // test the RT methods
  TEST_EQUAL(p.rts.size(), 1)
  TEST_EQUAL(p.rts[0] == rt, true)
  TEST_EQUAL(p.rts[0].retention_time_unit == TargetedExperimentHelper::RetentionTime::RTUnit::SECOND, true)
  TEST_EQUAL(p.rts[0].retention_time_type == TargetedExperimentHelper::RetentionTime::RTType::PREDICTED, true)
  TEST_REAL_SIMILAR(p.rts[0].getRT(), 5.0)

  // test the Peptide methods
  TEST_EQUAL(p.hasRetentionTime(), true)
  TEST_REAL_SIMILAR(p.getRetentionTime(), 5.0)
  TEST_EQUAL(p.getRetentionTimeUnit() == TargetedExperimentHelper::RetentionTime::RTUnit::SECOND, true)
  TEST_EQUAL(p.getRetentionTimeType() == TargetedExperimentHelper::RetentionTime::RTType::PREDICTED, true)

  TEST_EQUAL(p.getPeptideGroupLabel(), "")
  p.setPeptideGroupLabel("test1");
  TEST_EQUAL(p.getPeptideGroupLabel(), "test1")

  TEST_EQUAL(p.hasCharge(), false)
  p.setChargeState(-1);
  TEST_EQUAL(p.getChargeState(), -1)
  p.setChargeState(2);
  TEST_EQUAL(p.getChargeState(), 2)
}
END_SECTION

START_SECTION((TargetedExperiment::Compound))
{
  TargetedExperimentHelper::Compound p;
    
  TEST_EQUAL(p.hasRetentionTime(), false)
  TEST_EQUAL(p.rts.size(), 0)

  // add a RT
  TargetedExperimentHelper::RetentionTime rt;
  rt.setRT(5.0);
  rt.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
  rt.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::PREDICTED;
  p.rts.push_back(rt);

  // test the RT methods
  TEST_EQUAL(p.rts.size(), 1)
  TEST_EQUAL(p.rts[0] == rt, true)
  TEST_EQUAL(p.rts[0].retention_time_unit == TargetedExperimentHelper::RetentionTime::RTUnit::SECOND, true)
  TEST_EQUAL(p.rts[0].retention_time_type == TargetedExperimentHelper::RetentionTime::RTType::PREDICTED, true)
  TEST_REAL_SIMILAR(p.rts[0].getRT(), 5.0)

  // test the Compound methods
  TEST_EQUAL(p.hasRetentionTime(), true)
  TEST_REAL_SIMILAR(p.getRetentionTime(), 5.0)
  TEST_EQUAL(p.getRetentionTimeUnit() == TargetedExperimentHelper::RetentionTime::RTUnit::SECOND, true)
  TEST_EQUAL(p.getRetentionTimeType() == TargetedExperimentHelper::RetentionTime::RTType::PREDICTED, true)

  // TEST_EQUAL(p.getPeptideGroupLabel(), "")
  // p.setPeptideGroupLabel("test1");
  // TEST_EQUAL(p.getPeptideGroupLabel(), "test1")

  TEST_EQUAL(p.hasCharge(), false)
  p.setChargeState(-1);
  TEST_EQUAL(p.getChargeState(), -1)
  p.setChargeState(2);
  TEST_EQUAL(p.getChargeState(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Hash function tests
/////////////////////////////////////////////////////////////

START_SECTION((std::hash<TargetedExperimentHelper::CV>))
{
  TargetedExperimentHelper::CV cv1("id1", "fullname1", "version1", "URI1");
  TargetedExperimentHelper::CV cv2("id1", "fullname1", "version1", "URI1");
  TargetedExperimentHelper::CV cv3("id2", "fullname2", "version2", "URI2");

  std::hash<TargetedExperimentHelper::CV> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(cv1), hasher(cv2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(cv1), hasher(cv3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::CV> cv_set;
  cv_set.insert(cv1);
  cv_set.insert(cv2);
  cv_set.insert(cv3);
  TEST_EQUAL(cv_set.size(), 2) // cv1 and cv2 are equal
}
END_SECTION

START_SECTION((std::hash<TargetedExperimentHelper::Protein>))
{
  TargetedExperimentHelper::Protein p1;
  p1.id = "protein1";
  p1.sequence = "ACDEFGH";

  TargetedExperimentHelper::Protein p2;
  p2.id = "protein1";
  p2.sequence = "ACDEFGH";

  TargetedExperimentHelper::Protein p3;
  p3.id = "protein2";
  p3.sequence = "IKJLMNPQRST";

  std::hash<TargetedExperimentHelper::Protein> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(p1), hasher(p2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(p1), hasher(p3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::Protein> protein_set;
  protein_set.insert(p1);
  protein_set.insert(p2);
  protein_set.insert(p3);
  TEST_EQUAL(protein_set.size(), 2)
}
END_SECTION

START_SECTION((std::hash<TargetedExperimentHelper::RetentionTime>))
{
  TargetedExperimentHelper::RetentionTime rt1;
  rt1.setRT(5.0);
  rt1.software_ref = "software1";
  rt1.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
  rt1.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::PREDICTED;

  TargetedExperimentHelper::RetentionTime rt2;
  rt2.setRT(5.0);
  rt2.software_ref = "software1";
  rt2.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::SECOND;
  rt2.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::PREDICTED;

  TargetedExperimentHelper::RetentionTime rt3;
  rt3.setRT(10.0);
  rt3.software_ref = "software2";
  rt3.retention_time_unit = TargetedExperimentHelper::RetentionTime::RTUnit::MINUTE;
  rt3.retention_time_type = TargetedExperimentHelper::RetentionTime::RTType::LOCAL;

  std::hash<TargetedExperimentHelper::RetentionTime> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(rt1), hasher(rt2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(rt1), hasher(rt3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::RetentionTime> rt_set;
  rt_set.insert(rt1);
  rt_set.insert(rt2);
  rt_set.insert(rt3);
  TEST_EQUAL(rt_set.size(), 2)
}
END_SECTION

START_SECTION((std::hash<TargetedExperimentHelper::PeptideCompound>))
{
  TargetedExperimentHelper::PeptideCompound pc1;
  pc1.id = "compound1";
  pc1.setChargeState(2);

  TargetedExperimentHelper::PeptideCompound pc2;
  pc2.id = "compound1";
  pc2.setChargeState(2);

  TargetedExperimentHelper::PeptideCompound pc3;
  pc3.id = "compound2";
  pc3.setChargeState(3);

  std::hash<TargetedExperimentHelper::PeptideCompound> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(pc1), hasher(pc2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(pc1), hasher(pc3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::PeptideCompound> pc_set;
  pc_set.insert(pc1);
  pc_set.insert(pc2);
  pc_set.insert(pc3);
  TEST_EQUAL(pc_set.size(), 2)
}
END_SECTION

START_SECTION((std::hash<TargetedExperimentHelper::Compound>))
{
  TargetedExperimentHelper::Compound c1;
  c1.id = "compound1";
  c1.molecular_formula = "C6H12O6";
  c1.smiles_string = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O";
  c1.theoretical_mass = 180.063;

  TargetedExperimentHelper::Compound c2;
  c2.id = "compound1";
  c2.molecular_formula = "C6H12O6";
  c2.smiles_string = "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O";
  c2.theoretical_mass = 180.063;

  TargetedExperimentHelper::Compound c3;
  c3.id = "compound2";
  c3.molecular_formula = "H2O";
  c3.smiles_string = "O";
  c3.theoretical_mass = 18.015;

  std::hash<TargetedExperimentHelper::Compound> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(c1), hasher(c2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(c1), hasher(c3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::Compound> compound_set;
  compound_set.insert(c1);
  compound_set.insert(c2);
  compound_set.insert(c3);
  TEST_EQUAL(compound_set.size(), 2)
}
END_SECTION

START_SECTION((std::hash<TargetedExperimentHelper::Peptide>))
{
  TargetedExperimentHelper::Peptide pep1;
  pep1.id = "peptide1";
  pep1.sequence = "ACDEFGHIK";
  pep1.setPeptideGroupLabel("group1");

  TargetedExperimentHelper::Peptide pep2;
  pep2.id = "peptide1";
  pep2.sequence = "ACDEFGHIK";
  pep2.setPeptideGroupLabel("group1");

  TargetedExperimentHelper::Peptide pep3;
  pep3.id = "peptide2";
  pep3.sequence = "LMNPQRSTVWY";
  pep3.setPeptideGroupLabel("group2");

  std::hash<TargetedExperimentHelper::Peptide> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(pep1), hasher(pep2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(pep1), hasher(pep3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::Peptide> peptide_set;
  peptide_set.insert(pep1);
  peptide_set.insert(pep2);
  peptide_set.insert(pep3);
  TEST_EQUAL(peptide_set.size(), 2)
}
END_SECTION

START_SECTION((std::hash<TargetedExperimentHelper::Contact>))
{
  TargetedExperimentHelper::Contact c1;
  c1.id = "contact1";

  TargetedExperimentHelper::Contact c2;
  c2.id = "contact1";

  TargetedExperimentHelper::Contact c3;
  c3.id = "contact2";

  std::hash<TargetedExperimentHelper::Contact> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(c1), hasher(c2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(c1), hasher(c3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::Contact> contact_set;
  contact_set.insert(c1);
  contact_set.insert(c2);
  contact_set.insert(c3);
  TEST_EQUAL(contact_set.size(), 2)
}
END_SECTION

START_SECTION((std::hash<TargetedExperimentHelper::Publication>))
{
  TargetedExperimentHelper::Publication p1;
  p1.id = "pub1";

  TargetedExperimentHelper::Publication p2;
  p2.id = "pub1";

  TargetedExperimentHelper::Publication p3;
  p3.id = "pub2";

  std::hash<TargetedExperimentHelper::Publication> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(p1), hasher(p2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(p1), hasher(p3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::Publication> pub_set;
  pub_set.insert(p1);
  pub_set.insert(p2);
  pub_set.insert(p3);
  TEST_EQUAL(pub_set.size(), 2)
}
END_SECTION

START_SECTION((std::hash<TargetedExperimentHelper::Instrument>))
{
  TargetedExperimentHelper::Instrument i1;
  i1.id = "inst1";

  TargetedExperimentHelper::Instrument i2;
  i2.id = "inst1";

  TargetedExperimentHelper::Instrument i3;
  i3.id = "inst2";

  std::hash<TargetedExperimentHelper::Instrument> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(i1), hasher(i2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(i1), hasher(i3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::Instrument> inst_set;
  inst_set.insert(i1);
  inst_set.insert(i2);
  inst_set.insert(i3);
  TEST_EQUAL(inst_set.size(), 2)
}
END_SECTION

START_SECTION((std::hash<TargetedExperimentHelper::Prediction>))
{
  TargetedExperimentHelper::Prediction pred1;
  pred1.software_ref = "software1";
  pred1.contact_ref = "contact1";

  TargetedExperimentHelper::Prediction pred2;
  pred2.software_ref = "software1";
  pred2.contact_ref = "contact1";

  TargetedExperimentHelper::Prediction pred3;
  pred3.software_ref = "software2";
  pred3.contact_ref = "contact2";

  std::hash<TargetedExperimentHelper::Prediction> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(pred1), hasher(pred2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(pred1), hasher(pred3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::Prediction> pred_set;
  pred_set.insert(pred1);
  pred_set.insert(pred2);
  pred_set.insert(pred3);
  TEST_EQUAL(pred_set.size(), 2)
}
END_SECTION

START_SECTION((std::hash<TargetedExperimentHelper::Interpretation>))
{
  TargetedExperimentHelper::Interpretation interp1;
  interp1.ordinal = 1;
  interp1.rank = 1;
  interp1.iontype = Residue::BIon;

  TargetedExperimentHelper::Interpretation interp2;
  interp2.ordinal = 1;
  interp2.rank = 1;
  interp2.iontype = Residue::BIon;

  TargetedExperimentHelper::Interpretation interp3;
  interp3.ordinal = 2;
  interp3.rank = 2;
  interp3.iontype = Residue::YIon;

  std::hash<TargetedExperimentHelper::Interpretation> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(interp1), hasher(interp2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(interp1), hasher(interp3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::Interpretation> interp_set;
  interp_set.insert(interp1);
  interp_set.insert(interp2);
  interp_set.insert(interp3);
  TEST_EQUAL(interp_set.size(), 2)
}
END_SECTION

START_SECTION((std::hash<TargetedExperimentHelper::TraMLProduct>))
{
  TargetedExperimentHelper::TraMLProduct prod1;
  prod1.setMZ(500.5);
  prod1.setChargeState(2);

  TargetedExperimentHelper::TraMLProduct prod2;
  prod2.setMZ(500.5);
  prod2.setChargeState(2);

  TargetedExperimentHelper::TraMLProduct prod3;
  prod3.setMZ(600.6);
  prod3.setChargeState(3);

  std::hash<TargetedExperimentHelper::TraMLProduct> hasher;

  // Equal objects should have equal hashes
  TEST_EQUAL(hasher(prod1), hasher(prod2))

  // Different objects should (very likely) have different hashes
  TEST_NOT_EQUAL(hasher(prod1), hasher(prod3))

  // Test usability in unordered_set
  std::unordered_set<TargetedExperimentHelper::TraMLProduct> prod_set;
  prod_set.insert(prod1);
  prod_set.insert(prod2);
  prod_set.insert(prod3);
  TEST_EQUAL(prod_set.size(), 2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



