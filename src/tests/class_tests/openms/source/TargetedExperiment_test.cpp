// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TargetedExperiment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TargetedExperiment* ptr = nullptr;
TargetedExperiment* nullPointer = nullptr;
START_SECTION(TargetedExperiment())
{
	ptr = new TargetedExperiment();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~TargetedExperiment())
{
	delete ptr;
}
END_SECTION

TargetedExperiment te;

START_SECTION((TargetedExperiment(const TargetedExperiment &rhs)))
{
  TargetedExperiment t; 
  ReactionMonitoringTransition tr;
  t.addTransition(tr);

  {
    TargetedExperiment::Peptide p;
    p.id = "myPep";
    t.addPeptide(p);
  }

  {
    TargetedExperiment::Protein p;
    p.id = "myProtein";
    t.addProtein(p);
  }

  {
    TargetedExperiment::Compound c;
    c.id = "myCompound";
    t.addCompound(c);
  }

  TargetedExperiment t2(t); 
    
  TEST_TRUE(t2 == t)
}
END_SECTION

START_SECTION((bool operator==(const TargetedExperiment &rhs) const ))
{
  TargetedExperiment t; 
  ReactionMonitoringTransition tr;
  t.addTransition(tr);

  {
    TargetedExperiment::Peptide p;
    p.id = "myPep";
    t.addPeptide(p);
  }

  {
    TargetedExperiment::Protein p;
    p.id = "myProtein";
    t.addProtein(p);
  }

  {
    TargetedExperiment::Compound c;
    c.id = "myCompound";
    t.addCompound(c);
  }

  TargetedExperiment t2; 
  t2 = t;
    
  TEST_TRUE(t2 == t)
}
END_SECTION

START_SECTION((void setCVs(const std::vector< CV > &cvs)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<CV>& getCVs() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addCV(const CV &cv)))
{
  // TODO
}
END_SECTION

START_SECTION(( void setContacts(const std::vector< CVTermList > &contacts)))
{
  // TODO
}
END_SECTION

START_SECTION(const std::vector<CVTermList>& getContacts() const)
{
  // TODO
}
END_SECTION

START_SECTION( void addContact(const CVTermList &contact))
{
  // TODO
}
END_SECTION

START_SECTION( void setPublications(const std::vector< CVTermList > &publications))
{
  // TODO
}
END_SECTION

START_SECTION( const std::vector<CVTermList>& getPublications() const)
{
  // TODO
}
END_SECTION

START_SECTION(void addPublication(const CVTermList &publication))
{
  // TODO
}
END_SECTION

START_SECTION( void setInstruments(const std::vector< CVTermList > &instruments))
{
  // TODO
}
END_SECTION

START_SECTION(const std::vector<CVTermList>& getInstruments() const)
{
  // TODO
}
END_SECTION

START_SECTION(void addInstrument(const CVTermList &instrument))
{
  // TODO
}
END_SECTION

START_SECTION((void setSoftware(const std::vector< Software > &software)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Software>& getSoftware() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addSoftware(const Software &software)))
{
  // TODO
}
END_SECTION

START_SECTION((void setProteins(const std::vector< Protein > &proteins)))
{
  TargetedExperiment t; 
  TargetedExperiment::Protein p;

  std::vector<TargetedExperiment::Protein> proteins;
  proteins.push_back(p);
  t.setProteins(proteins);
    
  TEST_EQUAL(t.getProteins().size(), 1)
}
END_SECTION

START_SECTION((const std::vector<Protein>& getProteins() const ))
{
  TEST_EQUAL(te.getProteins().size(), 0)
}
END_SECTION

START_SECTION((bool hasProtein(const String & ref) const))
{
  TargetedExperiment t; 
  TargetedExperiment::Protein p;
  p.id = "myProtein";
  t.addProtein(p);
  TEST_EQUAL(t.hasProtein("myProtein"), true)
}
END_SECTION

START_SECTION((void addProtein(const Protein &protein)))
{
  TargetedExperiment t; 
  TargetedExperiment::Protein p;
  t.addProtein(p);
  TEST_EQUAL(t.getProteins().size(), 1)
}
END_SECTION

START_SECTION((void setCompounds(const std::vector< Compound > &rhs)))
{
  TargetedExperiment t; 
  TargetedExperiment::Compound c;

  std::vector<TargetedExperiment::Compound> compounds;
  compounds.push_back(c);
  t.setCompounds(compounds);
    
  TEST_EQUAL(t.getCompounds().size(), 1)
}
END_SECTION

START_SECTION((const std::vector<Compound>& getCompounds() const ))
{
  TEST_EQUAL(te.getCompounds().size(), 0)
}
END_SECTION

START_SECTION((bool hasCompound(const String & ref) const))
{
  TargetedExperiment t; 
  TargetedExperiment::Compound c;
  c.id = "myCompound";
  t.addCompound(c);
  TEST_EQUAL(t.hasCompound("myCompound"), true)
}
END_SECTION

START_SECTION((void addCompound(const Compound &rhs)))
{
  TargetedExperiment t; 
  TargetedExperiment::Compound c;
  t.addCompound(c);
  TEST_EQUAL(t.getCompounds().size(), 1)
}
END_SECTION

START_SECTION((void setPeptides(const std::vector< Peptide > &rhs)))
{
  TargetedExperiment t; 
  TargetedExperiment::Peptide p;

  std::vector<TargetedExperiment::Peptide> peptides;
  peptides.push_back(p);
  t.setPeptides(peptides);
    
  TEST_EQUAL(t.getPeptides().size(), 1)
}
END_SECTION

START_SECTION((const std::vector<Peptide>& getPeptides() const ))
{
  TEST_EQUAL(te.getPeptides().size(), 0)
}
END_SECTION

START_SECTION((void addPeptide(const Peptide &rhs)))
{
  TargetedExperiment t; 
  TargetedExperiment::Peptide p;
  t.addPeptide(p);
    
  TEST_EQUAL(t.getPeptides().size(), 1)
}
END_SECTION

START_SECTION((bool hasPeptide(const String & ref) const))
{
  TargetedExperiment t; 
  TargetedExperiment::Peptide p;
  p.id = "myPep";
  t.addPeptide(p);
  TEST_EQUAL(t.hasPeptide("myPep"), true)
}
END_SECTION

START_SECTION((void setTransitions(const std::vector< ReactionMonitoringTransition > &transitions)))
{
  TargetedExperiment t; 
  ReactionMonitoringTransition tr;

  std::vector<ReactionMonitoringTransition> transitions;
  transitions.push_back(tr);
  t.setTransitions(transitions);
    
  TEST_EQUAL(t.getTransitions().size(), 1)
}
END_SECTION

START_SECTION((const std::vector<ReactionMonitoringTransition>& getTransitions() const ))
{
  TEST_EQUAL(te.getTransitions().size(), 0)
}
END_SECTION

START_SECTION((void addTransition(const ReactionMonitoringTransition &transition)))
{
  TargetedExperiment t; 
  ReactionMonitoringTransition tr;

  t.addTransition(tr);
    
  TEST_EQUAL(t.getTransitions().size(), 1)
}
END_SECTION

START_SECTION((TargetedExperiment& operator=(const TargetedExperiment &rhs)))
{
  TargetedExperiment t; 
  ReactionMonitoringTransition tr;
  t.addTransition(tr);

  {
    TargetedExperiment::Peptide p;
    p.id = "myPep";
    t.addPeptide(p);
  }

  {
    TargetedExperiment::Protein p;
    p.id = "myProtein";
    t.addProtein(p);
  }

  {
    TargetedExperiment::Compound c;
    c.id = "myCompound";
    t.addCompound(c);
  }

  TargetedExperiment t2; 
  t2 = t;
    
  TEST_TRUE(t2 == t)
}
END_SECTION

START_SECTION( bool TargetedExperiment::containsInvalidReferences() )
{
  // wrong references
  {
    TargetedExperiment tr;
    TEST_EQUAL(tr.containsInvalidReferences(), false)

    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.id = "peptide1";
    peptide.protein_refs.push_back("protein1");
    tr.addPeptide(peptide);

    // now should be invalid due to a missing protein
    TEST_EQUAL(tr.containsInvalidReferences(), true)

    OpenMS::TargetedExperiment::Protein protein;
    protein.id = "protein1";
    tr.addProtein(protein);

    // now should be valid again
    TEST_EQUAL(tr.containsInvalidReferences(), false)

    OpenMS::ReactionMonitoringTransition t;
    t.setNativeID("tr1");
    t.setPeptideRef("peptide1");
    tr.addTransition(t);
    TEST_EQUAL(tr.containsInvalidReferences(), false)
    OpenMS::ReactionMonitoringTransition t2;
    t2.setNativeID("tr2");
    t2.setCompoundRef("compound1");
    tr.addTransition(t2);
    TEST_EQUAL(tr.containsInvalidReferences(), true)

    OpenMS::TargetedExperiment::Compound comp;
    comp.id = "compound1";
    tr.addCompound(comp);

    // now should be valid again
    TEST_EQUAL(tr.containsInvalidReferences(), false)
  }

  // duplications
  {
    TargetedExperiment tr;
    TEST_EQUAL(tr.containsInvalidReferences(), false)

    OpenMS::TargetedExperiment::Peptide peptide;
    peptide.id = "peptide1";
    tr.addPeptide(peptide);
    TEST_EQUAL(tr.containsInvalidReferences(), false)
    tr.addPeptide(peptide);
    TEST_EQUAL(tr.containsInvalidReferences(), true)
  }

  {
    TargetedExperiment tr;
    TEST_EQUAL(tr.containsInvalidReferences(), false)

    OpenMS::ReactionMonitoringTransition t;
    t.setNativeID("tr1");
    tr.addTransition(t);
    TEST_EQUAL(tr.containsInvalidReferences(), false)
    tr.addTransition(t);
    TEST_EQUAL(tr.containsInvalidReferences(), true)
  }

  {
    TargetedExperiment tr;
    TEST_EQUAL(tr.containsInvalidReferences(), false)

    OpenMS::TargetedExperiment::Protein protein;
    protein.id = "protein1";
    tr.addProtein(protein);
    TEST_EQUAL(tr.containsInvalidReferences(), false)
    tr.addProtein(protein);
    TEST_EQUAL(tr.containsInvalidReferences(), true)
  }
}
END_SECTION


START_SECTION(OpenMS::AASequence getAASequence(const OpenMS::TargetedExperiment::Peptide &peptide))
{
  OpenMS::TargetedExperiment::Peptide peptide;
  peptide.sequence = "TESTPEPTIDE";
  OpenMS::TargetedExperiment::Peptide::Modification modification;
  modification.avg_mass_delta = 79.9799;
  modification.location = 2;
  modification.mono_mass_delta = 79.966331;
  peptide.mods.push_back(modification);

  OpenMS::AASequence aas = TargetedExperimentHelper::getAASequence(peptide);
  OpenMS::String modified_sequence = "TES(Phospho)TPEPTIDE";
  TEST_EQUAL(aas.toUnmodifiedString(),peptide.sequence)
  //TEST_EQUAL(aas.toString(),modified_sequence)

  OpenMS::TargetedExperiment::Peptide peptide2;
  peptide2.sequence = "TESTPEPTIDER";
  OpenMS::TargetedExperiment::Peptide::Modification modification2;
  modification2.avg_mass_delta = 9.9296;
  modification2.location = 11;
  modification2.mono_mass_delta = 10.008269;
  peptide2.mods.push_back(modification2);

  OpenMS::AASequence aas2 = TargetedExperimentHelper::getAASequence(peptide2);
  OpenMS::String modified_sequence2 = "TESTPEPTIDER(Label:13C(6)15N(4))";
  TEST_EQUAL(aas2.toUnmodifiedString(),peptide2.sequence)

  OpenMS::TargetedExperiment::Peptide peptide3;
  peptide3.sequence = "TESTMPEPTIDE";
  OpenMS::TargetedExperiment::Peptide::Modification modification3;
  modification3.avg_mass_delta = 15.9994;
  modification3.location = 4;
  modification3.mono_mass_delta = 15.994915;
  peptide3.mods.push_back(modification3);

  OpenMS::AASequence aas3 = TargetedExperimentHelper::getAASequence(peptide3);
  OpenMS::String modified_sequence3 = "TESTM(Oxidation)PEPTIDE";
  TEST_EQUAL(aas3.toUnmodifiedString(),peptide3.sequence)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


