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

/////////////////////////////////////////////////////////////

#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/SYSTEM/PythonInfo.h>

#include <QDir>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(TextFile, "$Id$")

/////////////////////////////////////////////////////////////

START_SECTION((static bool canRun(String& python_executable, String& error_msg)))
  // test for missing python executable
  String py = "does_not_exist_@@";
  String error_msg;
  TEST_EQUAL(PythonInfo::canRun(py, error_msg), false)
  TEST_EQUAL(error_msg.hasSubstring("Python not found at"), true)

  auto tmp_file = File::getTemporaryFile();
  ofstream f(tmp_file); // create the file
  f.close(); 
  TEST_EQUAL(PythonInfo::canRun(tmp_file, error_msg), false)
  TEST_EQUAL(error_msg.hasSubstring("failed to run"), true)  

  py = "python";
  if (PythonInfo::canRun(py, error_msg))
  { 
    TEST_EQUAL(File::exists(py), true)
    TEST_EQUAL(QDir::isRelativePath(py.toQString()), false)
  }

END_SECTION

START_SECTION(bool PythonInfo::isPackageInstalled(const String& python_executable, const String& package_name))
  String error_msg;
  String py = "python";
  if (PythonInfo::canRun(py, error_msg))
  {
    TEST_EQUAL(PythonInfo::isPackageInstalled(py, "veryWeirdPackage___@@__@"), false)
    TEST_EQUAL(PythonInfo::isPackageInstalled(py, "math"), true)
  }
END_SECTION

START_SECTION(static String getVersion(const String& python_executable))
  
  String py = "python";
  String error_msg;
  if (PythonInfo::canRun(py, error_msg))
  {
    String version = PythonInfo::getVersion(py);
    TEST_EQUAL(version.empty(), false)
  }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
