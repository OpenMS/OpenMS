// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: George Rosenberger $
// $Authors: George Rosenberger, Hannes Roest $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/TraMLFile.h>
#include <OpenMS/FORMAT/SqliteConnector.h>
#include <OpenMS/KERNEL/MRMTransitionGroup.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/SYSTEM/File.h>

#include <boost/assign/std/vector.hpp>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/TransitionPQPFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TransitionPQPFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TransitionPQPFile* ptr = nullptr;
TransitionPQPFile* nullPointer = nullptr;

START_SECTION(TransitionPQPFile())
{
  ptr = new TransitionPQPFile();
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~TransitionPQPFile())
{
  delete ptr;
}
END_SECTION

START_SECTION( void convertTargetedExperimentToPQP(const char * filename, OpenMS::TargetedExperiment & targeted_exp))
{
  // see TOPP tool test
  NOT_TESTABLE
}
END_SECTION

START_SECTION( void convertPQPToTargetedExperiment(const char * filename, OpenMS::TargetedExperiment & targeted_exp))
{
  // see TOPP tool test
  NOT_TESTABLE
}
END_SECTION

START_SECTION( void validateTargetedExperiment(OpenMS::TargetedExperiment & targeted_exp))
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION([EXTRA] test reading PQP with empty GENE table (issue #8687))
{
  // Create a minimal PQP file with an empty GENE table
  // This tests the fix for issue #8687 where an empty GENE table caused
  // INNER JOIN to return zero rows
  String temp_file;
  NEW_TMP_FILE(temp_file);

  // Create the PQP database structure with empty GENE table
  {
    SqliteConnector conn(temp_file);

    // Create all required tables (based on PQP schema)
    const char* create_sql =
      "CREATE TABLE VERSION(ID INT NOT NULL);"
      "CREATE TABLE GENE(ID INT PRIMARY KEY NOT NULL, GENE_NAME TEXT NOT NULL, DECOY INT NOT NULL);"
      "CREATE TABLE PEPTIDE_GENE_MAPPING(PEPTIDE_ID INT NOT NULL, GENE_ID INT NOT NULL);"
      "CREATE TABLE PROTEIN(ID INT PRIMARY KEY NOT NULL, PROTEIN_ACCESSION TEXT NOT NULL, DECOY INT NOT NULL);"
      "CREATE TABLE PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID INT NOT NULL, PROTEIN_ID INT NOT NULL);"
      "CREATE TABLE PEPTIDE(ID INT PRIMARY KEY NOT NULL, UNMODIFIED_SEQUENCE TEXT NOT NULL, MODIFIED_SEQUENCE TEXT NOT NULL, DECOY INT NOT NULL);"
      "CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID INT NOT NULL, PEPTIDE_ID INT NOT NULL);"
      "CREATE TABLE COMPOUND(ID INT PRIMARY KEY NOT NULL, COMPOUND_NAME TEXT NOT NULL, SUM_FORMULA TEXT NOT NULL, SMILES TEXT NOT NULL, ADDUCTS TEXT NOT NULL, DECOY INT NOT NULL);"
      "CREATE TABLE PRECURSOR_COMPOUND_MAPPING(PRECURSOR_ID INT NOT NULL, COMPOUND_ID INT NOT NULL);"
      "CREATE TABLE PRECURSOR(ID INT PRIMARY KEY NOT NULL, TRAML_ID TEXT NULL, GROUP_LABEL TEXT NULL, PRECURSOR_MZ REAL NOT NULL, CHARGE INT NULL, LIBRARY_INTENSITY REAL NULL, LIBRARY_RT REAL NULL, LIBRARY_DRIFT_TIME REAL NULL, DECOY INT NOT NULL);"
      "CREATE TABLE TRANSITION_PRECURSOR_MAPPING(TRANSITION_ID INT NOT NULL, PRECURSOR_ID INT NOT NULL);"
      "CREATE TABLE TRANSITION_PEPTIDE_MAPPING(TRANSITION_ID INT NOT NULL, PEPTIDE_ID INT NOT NULL);"
      "CREATE TABLE TRANSITION(ID INT PRIMARY KEY NOT NULL, TRAML_ID TEXT NULL, PRODUCT_MZ REAL NOT NULL, CHARGE INT NULL, TYPE CHAR(1) NULL, ANNOTATION TEXT NULL, ORDINAL INT NULL, DETECTING INT NOT NULL, IDENTIFYING INT NOT NULL, QUANTIFYING INT NOT NULL, LIBRARY_INTENSITY REAL NULL, DECOY INT NOT NULL);";

    conn.executeStatement(create_sql);

    // Insert version
    conn.executeStatement("INSERT INTO VERSION (ID) VALUES (3);");

    // Insert a protein
    conn.executeStatement("INSERT INTO PROTEIN (ID, PROTEIN_ACCESSION, DECOY) VALUES (0, 'TestProtein', 0);");

    // Insert a peptide
    conn.executeStatement("INSERT INTO PEPTIDE (ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, DECOY) VALUES (0, 'PEPTIDE', 'PEPTIDE', 0);");

    // Insert peptide-protein mapping
    conn.executeStatement("INSERT INTO PEPTIDE_PROTEIN_MAPPING (PEPTIDE_ID, PROTEIN_ID) VALUES (0, 0);");

    // Insert a precursor
    conn.executeStatement("INSERT INTO PRECURSOR (ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_RT, LIBRARY_DRIFT_TIME, DECOY) VALUES (0, 'precursor_1', 'group1', 500.0, 2, 100.0, -1, 0);");

    // Insert precursor-peptide mapping
    conn.executeStatement("INSERT INTO PRECURSOR_PEPTIDE_MAPPING (PRECURSOR_ID, PEPTIDE_ID) VALUES (0, 0);");

    // Insert a transition
    conn.executeStatement("INSERT INTO TRANSITION (ID, TRAML_ID, PRODUCT_MZ, CHARGE, TYPE, ORDINAL, DETECTING, IDENTIFYING, QUANTIFYING, LIBRARY_INTENSITY, DECOY) VALUES (0, 'transition_1', 300.0, 1, 'y', 3, 1, 0, 1, 1000.0, 0);");

    // Insert transition-precursor mapping
    conn.executeStatement("INSERT INTO TRANSITION_PRECURSOR_MAPPING (TRANSITION_ID, PRECURSOR_ID) VALUES (0, 0);");

    // NOTE: GENE table and PEPTIDE_GENE_MAPPING are intentionally left empty!
  }

  // Now try to read the PQP file - this should NOT fail even with empty GENE table
  TransitionPQPFile pqp_reader;
  TargetedExperiment targeted_exp;

  pqp_reader.convertPQPToTargetedExperiment(temp_file.c_str(), targeted_exp, true);

  // Verify we got the transition data
  TEST_EQUAL(targeted_exp.getTransitions().size(), 1)
  TEST_EQUAL(targeted_exp.getPeptides().size(), 1)
  TEST_EQUAL(targeted_exp.getProteins().size(), 1)

  // Test with LightTargetedExperiment as well
  OpenSwath::LightTargetedExperiment light_exp;
  pqp_reader.convertPQPToTargetedExperiment(temp_file.c_str(), light_exp, true);

  TEST_EQUAL(light_exp.transitions.size(), 1)
  TEST_EQUAL(light_exp.compounds.size(), 1)
  TEST_EQUAL(light_exp.proteins.size(), 1)

}
END_SECTION

START_SECTION([EXTRA] test reading PQP with multiple gene mappings does not duplicate transitions)
{
  String temp_file;
  NEW_TMP_FILE(temp_file);

  {
    SqliteConnector conn(temp_file);

    const char* create_sql =
      "CREATE TABLE VERSION(ID INT NOT NULL);"
      "CREATE TABLE GENE(ID INT PRIMARY KEY NOT NULL, GENE_NAME TEXT NOT NULL, DECOY INT NOT NULL);"
      "CREATE TABLE PEPTIDE_GENE_MAPPING(PEPTIDE_ID INT NOT NULL, GENE_ID INT NOT NULL);"
      "CREATE TABLE PROTEIN(ID INT PRIMARY KEY NOT NULL, PROTEIN_ACCESSION TEXT NOT NULL, DECOY INT NOT NULL);"
      "CREATE TABLE PEPTIDE_PROTEIN_MAPPING(PEPTIDE_ID INT NOT NULL, PROTEIN_ID INT NOT NULL);"
      "CREATE TABLE PEPTIDE(ID INT PRIMARY KEY NOT NULL, UNMODIFIED_SEQUENCE TEXT NOT NULL, MODIFIED_SEQUENCE TEXT NOT NULL, DECOY INT NOT NULL);"
      "CREATE TABLE PRECURSOR_PEPTIDE_MAPPING(PRECURSOR_ID INT NOT NULL, PEPTIDE_ID INT NOT NULL);"
      "CREATE TABLE COMPOUND(ID INT PRIMARY KEY NOT NULL, COMPOUND_NAME TEXT NOT NULL, SUM_FORMULA TEXT NOT NULL, SMILES TEXT NOT NULL, ADDUCTS TEXT NOT NULL, DECOY INT NOT NULL);"
      "CREATE TABLE PRECURSOR_COMPOUND_MAPPING(PRECURSOR_ID INT NOT NULL, COMPOUND_ID INT NOT NULL);"
      "CREATE TABLE PRECURSOR(ID INT PRIMARY KEY NOT NULL, TRAML_ID TEXT NULL, GROUP_LABEL TEXT NULL, PRECURSOR_MZ REAL NOT NULL, CHARGE INT NULL, LIBRARY_INTENSITY REAL NULL, LIBRARY_RT REAL NULL, LIBRARY_DRIFT_TIME REAL NULL, DECOY INT NOT NULL);"
      "CREATE TABLE TRANSITION_PRECURSOR_MAPPING(TRANSITION_ID INT NOT NULL, PRECURSOR_ID INT NOT NULL);"
      "CREATE TABLE TRANSITION_PEPTIDE_MAPPING(TRANSITION_ID INT NOT NULL, PEPTIDE_ID INT NOT NULL);"
      "CREATE TABLE TRANSITION(ID INT PRIMARY KEY NOT NULL, TRAML_ID TEXT NULL, PRODUCT_MZ REAL NOT NULL, CHARGE INT NULL, TYPE CHAR(1) NULL, ANNOTATION TEXT NULL, ORDINAL INT NULL, DETECTING INT NOT NULL, IDENTIFYING INT NOT NULL, QUANTIFYING INT NOT NULL, LIBRARY_INTENSITY REAL NULL, DECOY INT NOT NULL);";

    conn.executeStatement(create_sql);
    conn.executeStatement("INSERT INTO VERSION (ID) VALUES (3);");
    conn.executeStatement("INSERT INTO GENE (ID, GENE_NAME, DECOY) VALUES (0, 'GeneA', 0);");
    conn.executeStatement("INSERT INTO GENE (ID, GENE_NAME, DECOY) VALUES (1, 'GeneB', 0);");
    conn.executeStatement("INSERT INTO PROTEIN (ID, PROTEIN_ACCESSION, DECOY) VALUES (0, 'TestProtein', 0);");
    conn.executeStatement("INSERT INTO PEPTIDE (ID, UNMODIFIED_SEQUENCE, MODIFIED_SEQUENCE, DECOY) VALUES (0, 'PEPTIDE', 'PEPTIDE', 0);");
    conn.executeStatement("INSERT INTO PEPTIDE_GENE_MAPPING (PEPTIDE_ID, GENE_ID) VALUES (0, 0);");
    conn.executeStatement("INSERT INTO PEPTIDE_GENE_MAPPING (PEPTIDE_ID, GENE_ID) VALUES (0, 1);");
    conn.executeStatement("INSERT INTO PEPTIDE_PROTEIN_MAPPING (PEPTIDE_ID, PROTEIN_ID) VALUES (0, 0);");
    conn.executeStatement("INSERT INTO PRECURSOR (ID, TRAML_ID, GROUP_LABEL, PRECURSOR_MZ, CHARGE, LIBRARY_RT, LIBRARY_DRIFT_TIME, DECOY) VALUES (0, 'precursor_1', 'group1', 500.0, 2, 100.0, -1, 0);");
    conn.executeStatement("INSERT INTO PRECURSOR_PEPTIDE_MAPPING (PRECURSOR_ID, PEPTIDE_ID) VALUES (0, 0);");
    conn.executeStatement("INSERT INTO TRANSITION (ID, TRAML_ID, PRODUCT_MZ, CHARGE, TYPE, ORDINAL, DETECTING, IDENTIFYING, QUANTIFYING, LIBRARY_INTENSITY, DECOY) VALUES (0, 'transition_1', 300.0, 1, 'y', 3, 1, 0, 1, 1000.0, 0);");
    conn.executeStatement("INSERT INTO TRANSITION_PRECURSOR_MAPPING (TRANSITION_ID, PRECURSOR_ID) VALUES (0, 0);");
  }

  TransitionPQPFile pqp_reader;
  TargetedExperiment targeted_exp;
  pqp_reader.convertPQPToTargetedExperiment(temp_file.c_str(), targeted_exp, true);

  TEST_EQUAL(targeted_exp.getTransitions().size(), 1)
  TEST_EQUAL(targeted_exp.getPeptides().size(), 1)

  OpenSwath::LightTargetedExperiment light_exp;
  pqp_reader.convertPQPToTargetedExperiment(temp_file.c_str(), light_exp, true);

  TEST_EQUAL(light_exp.transitions.size(), 1)
  TEST_EQUAL(light_exp.compounds.size(), 1)
  TEST_EQUAL(light_exp.proteins.size(), 1)

  MRMTransitionGroup<MSChromatogram, OpenSwath::LightTransition> transition_group;
  for (const auto& transition : light_exp.transitions)
  {
    transition_group.addTransition(transition, transition.getNativeID());
  }
  TEST_EQUAL(transition_group.getTransitions().size(), 1)
  TEST_EQUAL(transition_group.hasTransition("transition_1"), true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
