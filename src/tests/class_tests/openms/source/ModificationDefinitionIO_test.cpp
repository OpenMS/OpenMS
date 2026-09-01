// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/ModificationDefinitionIO.h>
///////////////////////////

#include <OpenMS/CHEMISTRY/AASequence.h>
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/ResidueModification.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

using namespace OpenMS;
using namespace std;

namespace
{
  // A defined modification with a formula, registered under a test-only namespace. ModificationsDB is
  // a process-wide singleton, so every section uses its own names.
  const ResidueModification* define(const std::string& id, char origin, const std::string& formula)
  {
    ResidueModification d;
    d.setId(id);
    d.setOrigin(origin);
    d.setTermSpecificity(ResidueModification::ANYWHERE);
    d.setFullId();
    d.setDiffFormula(EmpiricalFormula(formula));
    d.setDiffMonoMass(EmpiricalFormula(formula).getMonoWeight());
    return ModificationsDB::getInstance()->registerDefinition(d);
  }
}

START_TEST(ModificationDefinitionIO, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION([EXTRA] collect - only the DEFINED named modifications a run references)
{
  const ResidueModification* adduct = define("TestIO:Adduct", 'K', "C9H11N2O8P");
  const ResidueModification* var_only = define("TestIO:VarOnly", 'S', "HPO3");
  TEST_TRUE(adduct != nullptr)
  TEST_TRUE(var_only != nullptr)

  ProteinIdentification prot;
  prot.setIdentifier("run1");
  ProteinIdentification::SearchParameters sp;
  sp.variable_modifications.push_back("Oxidation (M)");       // CV: must NOT be collected
  sp.variable_modifications.push_back("TestIO:VarOnly (S)");  // defined, but on no hit: MUST be collected
  prot.setSearchParameters(sp);

  PeptideIdentification pep;
  pep.setIdentifier("run1");
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("PEPTK(TestIO:Adduct)M(Oxidation)IDE"));
  pep.insertHit(hit);

  // a second run that references nothing must not appear at all
  ProteinIdentification other;
  other.setIdentifier("run2");

  std::vector<ProteinIdentification> prots;
  prots.push_back(prot);
  prots.push_back(other);
  PeptideIdentificationList peps;
  peps.push_back(pep);

  auto defs = ModificationDefinitionIO::collect(prots, peps);
  TEST_EQUAL(defs.size(), 1)
  TEST_EQUAL(defs.count("run1"), 1)
  TEST_EQUAL(defs.count("run2"), 0)

  const std::set<const ResidueModification*>& run_defs = defs["run1"];
  TEST_EQUAL(run_defs.size(), 2)
  TEST_EQUAL(run_defs.count(adduct), 1)
  TEST_EQUAL(run_defs.count(var_only), 1)
  for (const ResidueModification* m : run_defs)
  {
    TEST_FALSE(m->getId() == "Oxidation")
  }
}
END_SECTION

START_SECTION([EXTRA] encode - deterministic bytes for the same definitions)
{
  const ResidueModification* a = define("TestIO:EncA", 'K', "C2H2O");
  const ResidueModification* b = define("TestIO:EncB", 'K', "CH2");
  std::set<const ResidueModification*> defs;
  defs.insert(a);
  defs.insert(b);

  const std::string once = ModificationDefinitionIO::encode(defs);
  const std::string twice = ModificationDefinitionIO::encode(defs);
  TEST_EQUAL(once, twice)
  TEST_TRUE(once.find("TestIO:EncA") != std::string::npos)
  TEST_TRUE(once.find("TestIO:EncB") != std::string::npos)

  std::vector<std::string> records = ResidueModification::splitDefinitionRecords(once);
  TEST_EQUAL(records.size(), 2)
  TEST_TRUE(records[0] < records[1]) // sorted
  TEST_EQUAL(ModificationDefinitionIO::encode(std::set<const ResidueModification*>()), "")
}
END_SECTION

START_SECTION([EXTRA] attach - unions with what is already stored and never duplicates)
{
  const ResidueModification* a = define("TestIO:AttA", 'K', "C2H2O");
  const ResidueModification* b = define("TestIO:AttB", 'R', "C2H2O");
  const std::string& key = Constants::UserParam::MODIFICATION_DEFINITIONS;

  ProteinIdentification::SearchParameters sp;
  ModificationDefinitionIO::attach(sp, std::set<const ResidueModification*>());
  TEST_FALSE(sp.metaValueExists(key)) // nothing to say, nothing written

  std::set<const ResidueModification*> only_a;
  only_a.insert(a);
  ModificationDefinitionIO::attach(sp, only_a);
  TEST_TRUE(sp.metaValueExists(key))
  const std::string v1 = sp.getMetaValue(key).toString();
  TEST_EQUAL(ResidueModification::splitDefinitionRecords(v1).size(), 1)

  // attaching the same thing again changes nothing
  ModificationDefinitionIO::attach(sp, only_a);
  TEST_EQUAL(sp.getMetaValue(key).toString(), v1)

  // attaching something new keeps what was there
  std::set<const ResidueModification*> only_b;
  only_b.insert(b);
  ModificationDefinitionIO::attach(sp, only_b);
  const std::string v2 = sp.getMetaValue(key).toString();
  TEST_EQUAL(ResidueModification::splitDefinitionRecords(v2).size(), 2)
  TEST_TRUE(v2.find("TestIO:AttA") != std::string::npos)
  TEST_TRUE(v2.find("TestIO:AttB") != std::string::npos)
}
END_SECTION

START_SECTION([EXTRA] registerFrom - registers the good records and skips a bad one loudly)
{
  const ModificationsDB* db = ModificationsDB::getInstance();
  // Build the records by hand so this section does not depend on define() having registered them.
  ResidueModification d;
  d.setId("TestIO:FromBlob");
  d.setOrigin('K');
  d.setTermSpecificity(ResidueModification::ANYWHERE);
  d.setFullId();
  d.setDiffFormula(EmpiricalFormula("C9H11N2O8P"));
  d.setDiffMonoMass(EmpiricalFormula("C9H11N2O8P").getMonoWeight());
  TEST_FALSE(db->hasDefinedModification("TestIO:FromBlob"))

  const std::string blob = d.toDefinitionString() + ";this is not a record;" + d.toDefinitionString();
  Size n = ModificationDefinitionIO::registerFrom(blob);
  TEST_EQUAL(n, 2) // the duplicate is registered idempotently, the garbage is skipped
  TEST_TRUE(db->hasDefinedModification("TestIO:FromBlob"))

  // ... and a peptide naming it now resolves, with the right chemistry
  AASequence seq = AASequence::fromString("AEADNLDDK(TestIO:FromBlob)K");
  TEST_REAL_SIMILAR(seq.getMonoWeight(), 1423.5504442334)
  TEST_EQUAL(seq.getFormula(Residue::Full, 0).toString(), "C54H86N15O28P1")

  // the SearchParameters overload
  ProteinIdentification::SearchParameters sp;
  TEST_EQUAL(ModificationDefinitionIO::registerFrom(sp), 0)
  sp.setMetaValue(Constants::UserParam::MODIFICATION_DEFINITIONS, d.toDefinitionString());
  TEST_EQUAL(ModificationDefinitionIO::registerFrom(sp), 1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
