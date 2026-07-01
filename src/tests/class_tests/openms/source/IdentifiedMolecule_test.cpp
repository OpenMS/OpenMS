// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/ID/IdentifiedMolecule.h>
#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/NASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>

///////////////////////////

START_TEST(IdentifiedMolecule, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

using ID = IdentificationData;
namespace IDI = OpenMS::IdentificationDataInternal;

// Shared fixture: the IdentificationData object owns the containers the
// references point into, so it must outlive the references below.
ID data;
ID::IdentifiedPeptideRef pref = data.registerIdentifiedPeptide(ID::IdentifiedPeptide(AASequence::fromString("PEPTIDE")));
ID::IdentifiedCompoundRef cref = data.registerIdentifiedCompound(ID::IdentifiedCompound("comp_id", EmpiricalFormula("C6H12O6"), "glucose"));
ID::IdentifiedOligoRef oref = data.registerIdentifiedOligo(ID::IdentifiedOligo(NASequence::fromString("ACGU")));

ID::IdentifiedMolecule* ptr = nullptr;
ID::IdentifiedMolecule* null_ptr = nullptr;

START_SECTION((IdentifiedMolecule()))
{
  ptr = new ID::IdentifiedMolecule();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((~IdentifiedMolecule()))
{
  delete ptr;
}
END_SECTION

START_SECTION((IdentifiedMolecule(IdentifiedPeptideRef ref)))
{
  ID::IdentifiedMolecule m(pref);
  TEST_EQUAL(m.getMoleculeType() == IDI::MoleculeType::PROTEIN, true)
  TEST_EQUAL(m.getIdentifiedPeptideRef() == pref, true)
}
END_SECTION

START_SECTION((IdentifiedMolecule(IdentifiedCompoundRef ref)))
{
  ID::IdentifiedMolecule m(cref);
  TEST_EQUAL(m.getMoleculeType() == IDI::MoleculeType::COMPOUND, true)
  TEST_EQUAL(m.getIdentifiedCompoundRef() == cref, true)
}
END_SECTION

START_SECTION((IdentifiedMolecule(IdentifiedOligoRef ref)))
{
  ID::IdentifiedMolecule m(oref);
  TEST_EQUAL(m.getMoleculeType() == IDI::MoleculeType::RNA, true)
  TEST_EQUAL(m.getIdentifiedOligoRef() == oref, true)
}
END_SECTION

START_SECTION((IdentifiedMolecule(const IdentifiedMolecule&)))
{
  ID::IdentifiedMolecule m(pref);
  ID::IdentifiedMolecule copy(m);
  TEST_EQUAL(copy.getMoleculeType() == IDI::MoleculeType::PROTEIN, true)
  TEST_EQUAL(copy == m, true)
}
END_SECTION

START_SECTION((MoleculeType getMoleculeType() const))
{
  TEST_EQUAL(ID::IdentifiedMolecule(pref).getMoleculeType() == IDI::MoleculeType::PROTEIN, true)
  TEST_EQUAL(ID::IdentifiedMolecule(cref).getMoleculeType() == IDI::MoleculeType::COMPOUND, true)
  TEST_EQUAL(ID::IdentifiedMolecule(oref).getMoleculeType() == IDI::MoleculeType::RNA, true)
  // default-constructed variant holds index 0 -> a (default) peptide reference
  TEST_EQUAL(ID::IdentifiedMolecule().getMoleculeType() == IDI::MoleculeType::PROTEIN, true)
}
END_SECTION

START_SECTION((IdentifiedPeptideRef getIdentifiedPeptideRef() const))
{
  ID::IdentifiedMolecule m(pref);
  TEST_EQUAL(m.getIdentifiedPeptideRef() == pref, true)
  // wrong molecule type throws
  TEST_EXCEPTION(Exception::IllegalArgument, ID::IdentifiedMolecule(cref).getIdentifiedPeptideRef())
  TEST_EXCEPTION(Exception::IllegalArgument, ID::IdentifiedMolecule(oref).getIdentifiedPeptideRef())
}
END_SECTION

START_SECTION((IdentifiedCompoundRef getIdentifiedCompoundRef() const))
{
  ID::IdentifiedMolecule m(cref);
  TEST_EQUAL(m.getIdentifiedCompoundRef() == cref, true)
  TEST_EXCEPTION(Exception::IllegalArgument, ID::IdentifiedMolecule(pref).getIdentifiedCompoundRef())
  TEST_EXCEPTION(Exception::IllegalArgument, ID::IdentifiedMolecule(oref).getIdentifiedCompoundRef())
}
END_SECTION

START_SECTION((IdentifiedOligoRef getIdentifiedOligoRef() const))
{
  ID::IdentifiedMolecule m(oref);
  TEST_EQUAL(m.getIdentifiedOligoRef() == oref, true)
  TEST_EXCEPTION(Exception::IllegalArgument, ID::IdentifiedMolecule(pref).getIdentifiedOligoRef())
  TEST_EXCEPTION(Exception::IllegalArgument, ID::IdentifiedMolecule(cref).getIdentifiedOligoRef())
}
END_SECTION

START_SECTION((std::string toString() const))
{
  TEST_STRING_EQUAL(ID::IdentifiedMolecule(pref).toString(), "PEPTIDE")
  // compounds are represented by their identifier
  TEST_STRING_EQUAL(ID::IdentifiedMolecule(cref).toString(), "comp_id")
  TEST_STRING_EQUAL(ID::IdentifiedMolecule(oref).toString(), "ACGU")
}
END_SECTION

START_SECTION((EmpiricalFormula getFormula(Size fragment_type = 0, Int charge = 0) const))
{
  // peptide: full (default fragment_type/charge) -> sequence formula
  TEST_EQUAL(ID::IdentifiedMolecule(pref).getFormula() == AASequence::fromString("PEPTIDE").getFormula(), true)
  // compound: stored empirical formula
  TEST_EQUAL(ID::IdentifiedMolecule(cref).getFormula() == EmpiricalFormula("C6H12O6"), true)
  // oligo: full sequence formula
  TEST_EQUAL(ID::IdentifiedMolecule(oref).getFormula() == NASequence::fromString("ACGU").getFormula(), true)
}
END_SECTION

START_SECTION((bool operator==(const IdentifiedMolecule& a, const IdentifiedMolecule& b)))
{
  ID::IdentifiedMolecule a(pref), b(pref), c(cref);
  TEST_EQUAL(a == b, true)
  TEST_EQUAL(a == c, false)
  // different molecule types compare unequal
  TEST_EQUAL(ID::IdentifiedMolecule(cref) == ID::IdentifiedMolecule(oref), false)
}
END_SECTION

START_SECTION((bool operator!=(const IdentifiedMolecule& a, const IdentifiedMolecule& b)))
{
  ID::IdentifiedMolecule a(pref), b(pref), c(cref);
  TEST_EQUAL(a != b, false)
  TEST_EQUAL(a != c, true)
}
END_SECTION

START_SECTION((bool operator<(const IdentifiedMolecule& a, const IdentifiedMolecule& b)))
{
  ID::IdentifiedMolecule p(pref), c(cref);
  // irreflexive: equal molecules are not ordered either way
  TEST_EQUAL(p < ID::IdentifiedMolecule(pref), false)
  TEST_EQUAL(ID::IdentifiedMolecule(pref) < p, false)
  // antisymmetric/total for distinct molecules: exactly one direction holds
  bool lt = p < c;
  bool gt = c < p;
  TEST_EQUAL(lt != gt, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
