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
#include <OpenMS/FORMAT/TraMLFile.h>
///////////////////////////

#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>

using namespace OpenMS;
using namespace std;

START_TEST(TraMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


TraMLFile * ptr = nullptr;
TraMLFile * nullPointer = nullptr;

START_SECTION((TraMLFile()))
{
  ptr = new TraMLFile;
  TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((~TraMLFile()))
{
  delete ptr;
}
END_SECTION

START_SECTION((void load(const String &filename, TargetedExperiment & id)))
{
  NOT_TESTABLE // tested below
}
END_SECTION

START_SECTION((void store(const String &filename, const TargetedExperiment &id) const))
{
  TraMLFile file;

  //load map
  TargetedExperiment exp_original;
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.traML"), exp_original);

  //store map
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename, exp_original);
  //load written map
  TargetedExperiment exp;
  file.load(tmp_filename, exp);

  //test if everything worked
  TEST_TRUE(exp == exp_original)

  // Test storing a minimal example
  {
    TargetedExperiment minimal_exp;
    NEW_TMP_FILE(tmp_filename);
    file.store(tmp_filename, minimal_exp);

    TargetedExperiment newexp;
    file.load(tmp_filename, newexp);

    // Test if everything worked
    //
    // The two objects are not exactly identical, while storing some CVs are
    // added that are not present in the newly instantiated object but get
    // added to the object when loaded.
    minimal_exp.setCVs(newexp.getCVs());
    TEST_TRUE(newexp == minimal_exp)
  }

  // Test storing a minimal example (with one protein/peptide/transition)
  {
    TargetedExperiment minimal_exp;
    TargetedExperimentHelper::Protein protein; 
    TargetedExperimentHelper::Peptide peptide; 
    ReactionMonitoringTransition tr; 
    minimal_exp.addProtein(protein);
    minimal_exp.addPeptide(peptide);
    minimal_exp.addTransition(tr);
    NEW_TMP_FILE(tmp_filename);
    file.store(tmp_filename, minimal_exp);

    TargetedExperiment newexp;
    file.load(tmp_filename, newexp);

    // Test if everything worked
    //
    // The two objects are not exactly identical, while storing some CVs are
    // added that are not present in the newly instantiated object but get
    // added to the object when loaded.
    minimal_exp.setCVs(newexp.getCVs()); 
    TEST_TRUE(newexp == minimal_exp)
  }
}
END_SECTION

START_SECTION((void equal()))
{
  TraMLFile file;

  TargetedExperiment exp_original;
  TargetedExperiment exp_second;
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.traML"), exp_original);
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.traML"), exp_second);

  TEST_TRUE(exp_second == exp_original)
}
END_SECTION

START_SECTION((void assign()))
{
  TraMLFile file;

  //load map
  TargetedExperiment exp_original;
  TargetedExperiment exp_added;
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.traML"), exp_original);

  //store map
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  exp_added = exp_original;
  TEST_EQUAL(exp_original.getTargetCVTerms().getCVTerms().size(), 1)
  TEST_EQUAL(exp_added.getTargetCVTerms().getCVTerms().size(), 1)

  TEST_TRUE(exp_added == exp_original)
}
END_SECTION

START_SECTION((void add()))
{
  TraMLFile file;

  //load map
  TargetedExperiment exp_original;
  TargetedExperiment exp_added;
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.traML"), exp_original);

  //store map
  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  exp_added += exp_original;

  TEST_TRUE(exp_added == exp_original)
}
END_SECTION

START_SECTION([EXTRA] bool isValid(const String & filename, std::ostream & os = std::cerr))
{
  std::string tmp_filename;
  TraMLFile file;
  TargetedExperiment e;

//written empty file
  NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename, e);
  TEST_EQUAL(file.isValid(tmp_filename, std::cerr), true);

//written filled file
  NEW_TMP_FILE(tmp_filename);
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.traML"), e);
  file.store(tmp_filename, e);
  TEST_EQUAL(file.isValid(tmp_filename, std::cerr), true);
}
END_SECTION

START_SECTION(bool isSemanticallyValid(const String & filename, StringList & errors, StringList & warnings))
{
  std::string tmp_filename;
  TraMLFile file;
  StringList errors, warnings;
  TargetedExperiment e;

//written empty file
  NEW_TMP_FILE(tmp_filename);
  file.store(tmp_filename, e);
  TEST_EQUAL(file.isSemanticallyValid(tmp_filename, errors, warnings), true);
  TEST_EQUAL(errors.size(), 0)
  TEST_EQUAL(warnings.size(), 0)

//written filled file
  NEW_TMP_FILE(tmp_filename);
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.traML"), e);
  file.store(tmp_filename, e);
//TEST_EQUAL(file.isSemanticallyValid(tmp_filename, errors, warnings),true);
  TEST_EQUAL(errors.size(), 0)
  TEST_EQUAL(warnings.size(), 0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
