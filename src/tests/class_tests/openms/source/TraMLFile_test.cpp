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
#include <OpenMS/CONCEPT/LogStream.h>

#include <sstream>
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

START_SECTION((void load(const std::string &filename, TargetedExperiment & id)))
{
  NOT_TESTABLE // tested below
}
END_SECTION

START_SECTION((void store(const std::string &filename, const TargetedExperiment &id) const))
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

START_SECTION([EXTRA] bool isValid(const std::string & filename, std::ostream & os = std::cerr))
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

START_SECTION(bool isSemanticallyValid(const std::string & filename, StringList & errors, StringList & warnings))
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
START_SECTION([EXTRA] store escapes XML special characters in CV term values)
{
  // A CV value with '&', '<' and '"' must be written escaped, or the stored file cannot be read.
  TargetedExperiment exp_original;
  TargetedExperiment::Contact contact;
  contact.id = "contact_1";
  CVTerm name_term("MS:1000586", "contact name", "MS", "A & B <lab> \"C\"", CVTerm::Unit());
  contact.addCVTerm(name_term);
  exp_original.addContact(contact);

  std::string tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  TraMLFile file;
  file.store(tmp_filename, exp_original);

  TargetedExperiment exp;
  file.load(tmp_filename, exp);
  ABORT_IF(exp.getContacts().size() != 1)
  const auto& terms = exp.getContacts()[0].getCVTerms();
  ABORT_IF(terms.count("MS:1000586") != 1)
  TEST_EQUAL(terms.at("MS:1000586")[0].getValue().toString(), "A & B <lab> \"C\"")
}
END_SECTION

START_SECTION([EXTRA] void load(const std::string& filename, TargetedExperiment& id) - CV term value checks)
{
  // The value checks run for every CV term and only warn: terms with an invalid value are kept,
  // obsolete or not. Values that OpenMS writes itself must not warn.
  std::ostringstream log; // XMLHandler::warning() logs to the warning stream (Debug) or the debug stream (Release)
  getGlobalLogWarn().insert(log);
  getGlobalLogDebug().insert(log);
  TargetedExperiment exp;
  TraMLFile().load(OPENMS_GET_TEST_DATA_PATH("TraMLFile_test_cv_values.TraML"), exp);
  getThreadLocalLogWarn().flush();
  getThreadLocalLogDebug().flush();
  getGlobalLogWarn().remove(log);
  getGlobalLogDebug().remove(log);
  const std::string msg = log.str();

  TEST_EQUAL(exp.getTransitions().size(), 1)
  ABORT_IF(exp.getTransitions().size() != 1)
  const ReactionMonitoringTransition& tr = exp.getTransitions()[0];
  TEST_REAL_SIMILAR(tr.getLibraryIntensity(), 100.0)
  TEST_EQUAL(tr.getCVTerms().count("MS:1000045"), 1) // current term, invalid value
  TEST_EQUAL(tr.getCVTerms().count("MS:1000023"), 1) // obsolete term, invalid value: kept
  TEST_EQUAL(exp.getPeptides().size(), 1)
  ABORT_IF(exp.getPeptides().size() != 1)
  TEST_REAL_SIMILAR(exp.getPeptides()[0].getRetentionTime(), 44.5)

  TEST_TRUE(msg.find("'MS:1000045 - collision energy' used in tag 'Transition' must have a floating-point value") != std::string::npos)
  TEST_TRUE(msg.find("'MS:1000023 - isolation width' used in tag 'Transition' must have a floating-point value") != std::string::npos)
  TEST_TRUE(msg.find("Obsolete CV term 'MS:1000023 - isolation width' used in tag 'Transition'") != std::string::npos)
  TEST_EQUAL(msg.find("must not have a value"), std::string::npos) // product ion intensity, iRT
  TEST_EQUAL(msg.find("MS:1002217"), std::string::npos) // boolean '1'
  TEST_EQUAL(msg.find("MS:1002031"), std::string::npos) // URI
}
END_SECTION

END_TEST
