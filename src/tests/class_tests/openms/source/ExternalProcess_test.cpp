// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
#include <OpenMS/SYSTEM/ExternalProcess.h>
///////////////////////////

#include <OpenMS/config.h>

using namespace OpenMS;
using namespace std;

// we just need ANY commandline tool available on (hopefully) all boxes.
// note that commands like "dir" or "type" are only known within cmd.exe and are not actual executables (unlike on Linux)
#ifdef OPENMS_WINDOWSPLATFORM
  const QString exe = "cmd";
  const QStringList args = QStringList() << "/C" << "echo hi";
  const QStringList args_broken = QStringList() << "/C" << "doesnotexist";
#else
  const QString exe = "ls";
  const QStringList args("-l");
  const QStringList args_broken = QStringList() << "-0";
#endif // 

START_TEST(ExternalProcess, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
START_SECTION(ExternalProcess())
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(ExternalProcess(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr))
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(~ExternalProcess())
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(void setCallbacks(std::function<void(const String&)> callbackStdOut, std::function<void(const String&)> callbackStdErr))
  NOT_TESTABLE; // tested below
END_SECTION

START_SECTION(RETURNSTATE run(const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose, String& error_msg))
{
  String error_msg;
  { // without callbacks
    ExternalProcess ep;
    String error_msg;
    auto r = ep.run(exe, args, "", true, error_msg);
    TEST_EQUAL(r == ExternalProcess::RETURNSTATE::SUCCESS, true)
    TEST_EQUAL(error_msg.size(), 0)

    r = ep.run("this_exe_does_not_exist", args, "", true, error_msg);
    TEST_EQUAL(r == ExternalProcess::RETURNSTATE::FAILED_TO_START, true)
    TEST_NOT_EQUAL(error_msg.size(), 0);

    r = ep.run(exe, args_broken, "", true, error_msg);
    TEST_EQUAL(r == ExternalProcess::RETURNSTATE::NONZERO_EXIT, true)
    TEST_NOT_EQUAL(error_msg.size(), 0);
  }
  { // with callbacks
    String all_out, all_err;
    auto l_out = [&](const String& out) {all_out += out;};
    auto l_err = [&](const String& out) {all_err += out;};
    ExternalProcess ep(l_out, l_err);
    auto r = ep.run(exe, args, "", true, error_msg);
    TEST_EQUAL(r == ExternalProcess::RETURNSTATE::SUCCESS, true)
    TEST_EQUAL(error_msg.size(), 0);
    TEST_NOT_EQUAL(all_out.size(), 0)
    TEST_EQUAL(all_err.size(), 0)
    all_out.clear();
    all_err.clear();

    r = ep.run(exe, args_broken, "", false, error_msg);
    TEST_EQUAL(r == ExternalProcess::RETURNSTATE::NONZERO_EXIT, true)
    TEST_NOT_EQUAL(error_msg.size(), 0);
    TEST_EQUAL(all_out.size(), 0)
    std::cout << all_out << "\n\n";
    TEST_NOT_EQUAL(all_err.size(), 0)
    all_out.clear();
    all_err.clear();

    ep.setCallbacks(l_err, l_out); // swap callbacks
    r = ep.run(exe, args_broken, "", false, error_msg);
    TEST_EQUAL(r == ExternalProcess::RETURNSTATE::NONZERO_EXIT, true)
    TEST_NOT_EQUAL(error_msg.size(), 0);
    TEST_NOT_EQUAL(all_out.size(), 0)
    TEST_EQUAL(all_err.size(), 0)
    all_out.clear();
    all_err.clear();
  }
}
END_SECTION

START_SECTION(ExternalProcess::RETURNSTATE run(QWidget* parent, const QString& exe, const QStringList& args, const QString& working_dir, const bool verbose = false))
 NOT_TESTABLE // tested above..
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST


