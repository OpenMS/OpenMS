// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/SearchEngineBase.h>
///////////////////////////

#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/MSExperiment.h>

#include <string>

using namespace OpenMS;
using namespace std;

//test class with optional parameters
class SearchEngineBaseTest
  : public SearchEngineBase
{
  public:
    SearchEngineBaseTest()
      : SearchEngineBase("SearchEngineBaseTest", "A test class", false, {}, false)
    {
      char* var = (char*)("OPENMS_DISABLE_UPDATE_CHECK=ON");
#ifdef OPENMS_WINDOWSPLATFORM
      _putenv(var);
#else
      putenv(var);
#endif
      //main(0,nullptr);
    }

    void registerOptionsAndFlags_() override
    {
      // these flags are currently used in public functions of SEB
      
      // used in: getRawfileName()
      registerInputFile_("in", "<file>", "", "Input file");
      setValidFormats_("in", { "mzML" } );

      // used in: getDBFilename()
      registerInputFile_("database", "<file>", "", "FASTA file", true, false, {"skipexists"});
      setValidFormats_("database", { "FASTA" } );
    }


    ExitCodes main_(int /*argc*/ , const char** /*argv*/) override
    {
      // check raw file (must contain centroided MS2 spectra)
      String in = getRawfileName();

      // check if DB is found (no PATH lookup possible here, since we do not control the OpenMS.ini; so usefulness is limited)
      String db = getDBFilename();

      return EXECUTION_OK;
    }
};

std::vector<const char*> toArgV(const StringList& args)
{
  std::vector<const char*> args_p;
  for (auto& a : args) args_p.push_back(a.c_str());
  return args_p;
}

/////////////////////////////////////////////////////////////

START_TEST(SearchEngineBase, "$Id$");

/////////////////////////////////////////////////////////////



SearchEngineBaseTest* ptr = nullptr;
SearchEngineBaseTest* nullPointer = nullptr;
START_SECTION(SearchEngineBase(const String& name, const String& description, bool official = true, const std::vector<Citation>& citations = {}, bool toolhandler_test = true))
	ptr = new SearchEngineBaseTest();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~SearchEngineBase()))
	delete ptr;
END_SECTION

START_SECTION(String getRawfileName(int ms_level = 2) const;)
  // collect arguments to TOPP/SEB
  String db(OPENMS_GET_TEST_DATA_PATH("degenerate_cases/empty.fasta"));
  String infile_empty(OPENMS_GET_TEST_DATA_PATH("degenerate_cases/empty_spec.mzML"));
  String infile_profile(OPENMS_GET_TEST_DATA_PATH("Single_MS2_profileMode.mzML"));
  StringList args;

  // need local scopes; calling TOPPBase::main() twice leads to error

  //
  // TEST UNKNOWN-only DATA (via an empty spectrum):
  //
  {
    SearchEngineBaseTest instance;
    args = { { "main.exe", "-in", infile_empty, "-database", db } };
    // run it
    // --> fails, since the experiment does not contain centroided spectra (only one 'UNKNOWN')
    TEST_EQUAL(TOPPBase::ExitCodes::UNKNOWN_ERROR == instance.main(args.size(), toArgV(args).data()), true);
  }
  {
    SearchEngineBaseTest instance;
    args = { { "main.exe", "-in", infile_empty, "-database", db, "-force" } };
    // run it
    // --> ok, due to -force
    TEST_EQUAL(TOPPBase::ExitCodes::EXECUTION_OK == instance.main(args.size(), toArgV(args).data()), true);
  }

  //
  // TEST PROFILE-only DATA:
  //
  {
    SearchEngineBaseTest instance;
    args = { { "main.exe", "-in", infile_profile, "-database", db } };
    // run it
    // --> fails, since the experiment contains a spectrum of type 'PROFILE'
    TEST_EQUAL(TOPPBase::ExitCodes::UNKNOWN_ERROR == instance.main(args.size(), toArgV(args).data()), true);
  }
  {
    SearchEngineBaseTest instance;
    args = { { "main.exe", "-in", infile_profile, "-database", db, "-force" } };
    // run it
    // --> ok, due to -force
    TEST_EQUAL(TOPPBase::ExitCodes::EXECUTION_OK == instance.main(args.size(), toArgV(args).data()), true);
  }

END_SECTION


START_SECTION(String getDBFilename(String db = "") const)
  NOT_TESTABLE // tested above
END_SECTION



/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



