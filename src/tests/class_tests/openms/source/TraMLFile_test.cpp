// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
  TraMLFile file;
  TargetedExperiment exp;
  file.load(OPENMS_GET_TEST_DATA_PATH("ToyExample1.traML"), exp);
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
  //file.store("include.traML", TargetedExperiment());
  //load written map
  TargetedExperiment exp;
  file.load(tmp_filename, exp);

  //test if everything worked
  TEST_EQUAL(exp == exp_original, true)

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
    TEST_EQUAL(newexp == minimal_exp, true)
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
    TEST_EQUAL(newexp == minimal_exp, true)
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

  TEST_EQUAL(exp_second == exp_original, true)
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

  TEST_EQUAL(exp_added == exp_original, true)
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

  TEST_EQUAL(exp_added == exp_original, true)
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
