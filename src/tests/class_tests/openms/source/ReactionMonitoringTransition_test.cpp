// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
///////////////////////////

#include <unordered_set>
#include <unordered_map>

using namespace OpenMS;
using namespace std;

START_TEST(ReactionMonitoringTransition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ReactionMonitoringTransition* ptr = nullptr;
ReactionMonitoringTransition* nullPointer = nullptr;

START_SECTION(ReactionMonitoringTransition())
{
  ptr = new ReactionMonitoringTransition();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~ReactionMonitoringTransition())
{
  delete ptr;
}
END_SECTION

OpenMS::ReactionMonitoringTransition transition = ReactionMonitoringTransition();

CVTerm charge_cv;
String charge_cv_acc = "MS:1000041";
charge_cv.setCVIdentifierRef("MS");
charge_cv.setAccession(charge_cv_acc);
charge_cv.setName("charge state");
charge_cv.setValue(3);

/////////////////////////////////////////////////////////////
// Copy constructor, move constructor, assignment operator, move assignment operator, equality

START_SECTION((ReactionMonitoringTransition(const ReactionMonitoringTransition &rhs)))
{
  ReactionMonitoringTransition tr1, tr2, tr3;

  tr1.addPrecursorCVTerm(charge_cv);
  tr1.setPrecursorMZ(42.0);
	tr2 = ReactionMonitoringTransition(tr1);
  TEST_TRUE(tr1 == tr2)
  ReactionMonitoringTransition::Prediction p;
  p.contact_ref = "dummy";
  tr1.setPrediction(p);
  tr1.setIdentifyingTransition(false);
  tr1.setDetectingTransition(false);
  tr1.setQuantifyingTransition(false);
	tr3 = ReactionMonitoringTransition(tr1);
  TEST_TRUE(tr1 == tr3)
  TEST_EQUAL(tr1 == tr2, false)


}
END_SECTION

START_SECTION((ReactionMonitoringTransition(ReactionMonitoringTransition &&rhs)))
{
  // Ensure that ReactionMonitoringTransition has a no-except move constructor (otherwise
  // std::vector is inefficient and will copy instead of move).
  TEST_EQUAL(noexcept(ReactionMonitoringTransition(std::declval<ReactionMonitoringTransition&&>())), true)

  ReactionMonitoringTransition tr1, tr2, tr3;

  ReactionMonitoringTransition::Prediction pred;
  pred.contact_ref = "dummy";
  tr1.setPrediction(pred);
  tr1.addPrecursorCVTerm(charge_cv);
  tr1.setPrecursorMZ(42.0);
  tr1.setCompoundRef("test_ref");

  auto orig = tr1;
	tr2 = ReactionMonitoringTransition(std::move(tr1));
  TEST_TRUE(orig == tr2);

  TEST_EQUAL(tr2.hasPrecursorCVTerms(), true);
  TEST_EQUAL(tr2.hasPrediction(), true);
  TEST_EQUAL(tr2.getPrediction().contact_ref, "dummy");
  TEST_EQUAL(tr2.getPrecursorCVTermList().hasCVTerm(charge_cv_acc), true)
  TEST_EQUAL(tr2.getCompoundRef(), "test_ref")

  TEST_EQUAL(tr1.hasPrecursorCVTerms(), false); // its gone
  TEST_EQUAL(tr1.hasPrediction(), false); // its gone
  TEST_EQUAL(tr1.getCompoundRef(), "") // its gone

  ReactionMonitoringTransition::Prediction p;
  p.contact_ref = "dummy";
  orig.setPrediction(p);
  orig.setIdentifyingTransition(false);
  orig.setDetectingTransition(false);
  orig.setQuantifyingTransition(false);
  tr1 = orig;
	tr3 = ReactionMonitoringTransition(std::move(tr1));
  TEST_TRUE(orig == tr3)
  TEST_EQUAL(orig == tr2, false)


}
END_SECTION

START_SECTION((ReactionMonitoringTransition& operator=(const ReactionMonitoringTransition &rhs)))
{
  ReactionMonitoringTransition tr1, tr2, tr3;

  tr1.addPrecursorCVTerm(charge_cv);
  tr1.setPrecursorMZ(42.0);
	tr2 = tr1;
  TEST_TRUE(tr1 == tr2)
  ReactionMonitoringTransition::Prediction p;
  p.contact_ref = "dummy";
  tr1.setPrediction(p);
  tr3 = tr1;
  TEST_TRUE(tr1 == tr3)
  TEST_EQUAL(tr1 == tr2, false)

  tr1.setDetectingTransition(false);
  TEST_EQUAL(tr1 == tr3, false)
  tr1.setIdentifyingTransition(true);
  tr1.setQuantifyingTransition(false);
  TEST_EQUAL(tr1 == tr3, false)
  tr3 = tr1;
  TEST_TRUE(tr1 == tr3)
}
END_SECTION


START_SECTION((void setName(const String &name)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setName("test_tr");

  TEST_EQUAL(tr.getName(), "test_tr")
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  TEST_EQUAL(transition.getName(), "")
}
END_SECTION

START_SECTION((void setPeptideRef(const String &peptide_ref)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setPeptideRef("test_ref");

  TEST_EQUAL(tr.getPeptideRef(), "test_ref")
}
END_SECTION

START_SECTION((const String& getPeptideRef() const ))
{
  TEST_EQUAL(transition.getPeptideRef(), "")
}
END_SECTION

START_SECTION((void setCompoundRef(const String &compound_ref)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setCompoundRef("test_ref");

  TEST_EQUAL(tr.getCompoundRef(), "test_ref")
}
END_SECTION

START_SECTION((const String& getCompoundRef() const ))
{
  TEST_EQUAL(transition.getCompoundRef(), "")
}
END_SECTION

START_SECTION((void setPrecursorMZ(double mz)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setPrecursorMZ(42.0);

  TEST_REAL_SIMILAR(tr.getPrecursorMZ(), 42.0)
}
END_SECTION

START_SECTION((double getPrecursorMZ() const ))
{
  TEST_REAL_SIMILAR(transition.getPrecursorMZ(), 0.0)
}
END_SECTION

START_SECTION((void setPrecursorCVTermList(const CVTermList &list)))
{
  CVTermList list;
  list.addCVTerm(charge_cv);
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setPrecursorCVTermList(list);

  TEST_EQUAL(tr.getPrecursorCVTermList().hasCVTerm(charge_cv_acc ), true)
}
END_SECTION

START_SECTION((bool hasPrecursorCVTerms() const))
{
  CVTermList list;
  list.addCVTerm(charge_cv);
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  TEST_EQUAL(tr.hasPrecursorCVTerms(), false)
  tr.setPrecursorCVTermList(list);
  TEST_EQUAL(tr.hasPrecursorCVTerms(), true)
}
END_SECTION

START_SECTION((void addPrecursorCVTerm(const CVTerm &cv_term)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  TEST_EQUAL(tr.hasPrecursorCVTerms(), false)
  tr.addPrecursorCVTerm(charge_cv);
  TEST_EQUAL(tr.hasPrecursorCVTerms(), true)

  TEST_EQUAL(tr.getPrecursorCVTermList().hasCVTerm(charge_cv_acc ), true)
}
END_SECTION

START_SECTION((const CVTermList& getPrecursorCVTermList() const ))
{
  CVTermList list;
  list.addCVTerm(charge_cv);
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setPrecursorCVTermList(list);

  TEST_EQUAL(tr.getPrecursorCVTermList() == list, true)
}
END_SECTION

START_SECTION((void setProductMZ(double mz)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setProductMZ(42.0);

  TEST_REAL_SIMILAR(tr.getProductMZ(), 42.0)
}
END_SECTION

START_SECTION((double getProductMZ() const ))
{
  TEST_REAL_SIMILAR(transition.getProductMZ(), 0.0)
}
END_SECTION

START_SECTION((void setProduct(Product product)))
{
  auto product = OpenMS::ReactionMonitoringTransition::Product();
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setProduct(product);

  TEST_EQUAL(tr.getProduct() == product, true)
}
END_SECTION

START_SECTION((const Product & getProduct() const))
{
  TEST_EQUAL(transition.getProduct() == OpenMS::ReactionMonitoringTransition::Product(), true)
}
END_SECTION

START_SECTION((void addProductCVTerm(const CVTerm &cv_term)))
{
  // TODO
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::isDetectingTransition() const))
{
  TEST_EQUAL(transition.isDetectingTransition(), true)
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::setDetectingTransition(bool val)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setDetectingTransition(false);
  TEST_EQUAL(tr.isDetectingTransition(), false)
  tr.setDetectingTransition(true);
  TEST_EQUAL(tr.isDetectingTransition(), true)
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::isIdentifyingTransition() const))
{
  TEST_EQUAL(transition.isIdentifyingTransition(), false)
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::setIdentifyingTransition(bool val)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setIdentifyingTransition(true);
  TEST_EQUAL(tr.isIdentifyingTransition(), true)
  tr.setIdentifyingTransition(false);
  TEST_EQUAL(tr.isIdentifyingTransition(), false)
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::isQuantifyingTransition() const))
{
  TEST_EQUAL(transition.isQuantifyingTransition(), true)
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::setQuantifyingTransition(bool val)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setQuantifyingTransition(false);
  TEST_EQUAL(tr.isQuantifyingTransition(), false)
  tr.setQuantifyingTransition(true);
  TEST_EQUAL(tr.isQuantifyingTransition(), true)
}
END_SECTION


START_SECTION((bool operator==(const ReactionMonitoringTransition &rhs) const ))
{
  ReactionMonitoringTransition tr1, tr2;
  TEST_TRUE(tr1 == tr2)

  tr1.addPrecursorCVTerm(charge_cv);
  TEST_EQUAL(tr1 == tr2, false)
  tr2.addPrecursorCVTerm(charge_cv);
  TEST_TRUE(tr1 == tr2)

  tr1.setDetectingTransition(false);
  TEST_EQUAL(tr1 == tr2, false)
  tr2.setDetectingTransition(false);
  TEST_TRUE(tr1 == tr2)

}
END_SECTION

START_SECTION((bool operator!=(const ReactionMonitoringTransition &rhs) const ))
{
  ReactionMonitoringTransition tr1, tr2;
  TEST_EQUAL(tr1 != tr2, false)

  tr1.addPrecursorCVTerm(charge_cv);
  TEST_FALSE(tr1 == tr2)
}
END_SECTION

/////////////////////////////////////////////////////////////
// Hash function tests
/////////////////////////////////////////////////////////////

START_SECTION(([EXTRA] std::hash<ReactionMonitoringTransition>))
{
  // Test that equal objects have equal hashes
  ReactionMonitoringTransition tr1, tr2;
  std::hash<ReactionMonitoringTransition> hasher;

  // Default constructed objects should have equal hashes
  TEST_EQUAL(hasher(tr1), hasher(tr2))

  // Modify both identically and verify hashes still match
  tr1.setName("test_transition");
  tr2.setName("test_transition");
  TEST_EQUAL(hasher(tr1), hasher(tr2))

  tr1.setPeptideRef("peptide_ref_1");
  tr2.setPeptideRef("peptide_ref_1");
  TEST_EQUAL(hasher(tr1), hasher(tr2))

  tr1.setPrecursorMZ(500.0);
  tr2.setPrecursorMZ(500.0);
  TEST_EQUAL(hasher(tr1), hasher(tr2))

  tr1.setProductMZ(250.0);
  tr2.setProductMZ(250.0);
  TEST_EQUAL(hasher(tr1), hasher(tr2))

  // Test with CV terms
  CVTerm cv_term;
  cv_term.setCVIdentifierRef("MS");
  cv_term.setAccession("MS:1000041");
  cv_term.setName("charge state");
  cv_term.setValue(2);
  tr1.addPrecursorCVTerm(cv_term);
  tr2.addPrecursorCVTerm(cv_term);
  TEST_EQUAL(hasher(tr1), hasher(tr2))

  // Test with prediction
  ReactionMonitoringTransition::Prediction pred;
  pred.contact_ref = "contact_1";
  pred.software_ref = "software_1";
  tr1.setPrediction(pred);
  tr2.setPrediction(pred);
  TEST_EQUAL(hasher(tr1), hasher(tr2))

  // Test with transition flags
  tr1.setDetectingTransition(false);
  tr2.setDetectingTransition(false);
  TEST_EQUAL(hasher(tr1), hasher(tr2))

  tr1.setIdentifyingTransition(true);
  tr2.setIdentifyingTransition(true);
  TEST_EQUAL(hasher(tr1), hasher(tr2))

  // Test that different objects have different hashes (not guaranteed but expected)
  ReactionMonitoringTransition tr3;
  tr3.setName("different_name");
  TEST_NOT_EQUAL(hasher(tr1), hasher(tr3))

  ReactionMonitoringTransition tr4;
  tr4.setPrecursorMZ(999.0);
  TEST_NOT_EQUAL(hasher(tr1), hasher(tr4))

  // Test use in unordered_set
  std::unordered_set<ReactionMonitoringTransition> transition_set;
  transition_set.insert(tr1);
  transition_set.insert(tr2); // Should not insert (equal to tr1)
  transition_set.insert(tr3);
  transition_set.insert(tr4);
  TEST_EQUAL(transition_set.size(), 3) // tr1 and tr2 are equal, so only 3 unique

  // Test use in unordered_map
  std::unordered_map<ReactionMonitoringTransition, int> transition_map;
  transition_map[tr1] = 1;
  transition_map[tr2] = 2; // Should overwrite value for tr1
  transition_map[tr3] = 3;
  TEST_EQUAL(transition_map.size(), 2)
  TEST_EQUAL(transition_map[tr1], 2) // Value should be 2 (overwritten)
  TEST_EQUAL(transition_map[tr3], 3)

  // Test lookup
  TEST_EQUAL(transition_set.count(tr1), 1)
  TEST_EQUAL(transition_set.count(tr2), 1) // Equal to tr1
  TEST_EQUAL(transition_set.count(tr3), 1)

  ReactionMonitoringTransition tr5;
  tr5.setName("not_in_set");
  TEST_EQUAL(transition_set.count(tr5), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



