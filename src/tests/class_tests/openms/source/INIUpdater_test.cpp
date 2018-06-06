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
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/INIUpdater.h>
///////////////////////////

#include <OpenMS/DATASTRUCTURES/ListUtils.h>

using namespace OpenMS;
using namespace std;

START_TEST(INIUpdater, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

INIUpdater* ptr = nullptr;
INIUpdater* null_ptr = nullptr;
START_SECTION(INIUpdater())
{
	ptr = new INIUpdater();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~INIUpdater())
{
	delete ptr;
}
END_SECTION

START_SECTION((StringList getToolNamesFromINI(const Param &ini) const))
{
  Param p;
  INIUpdater i;
  StringList names = i.getToolNamesFromINI(p);

  TEST_EQUAL(names.size(), 0)

  p.setValue("FeatureFinder:version","1.9");
  p.setValue("SomeTool:version","whatever");
  names = i.getToolNamesFromINI(p);

  TEST_EQUAL(names.size(), 2)

  p.setValue("BrokenTool:version2","1.9");
  names = i.getToolNamesFromINI(p);

  TEST_EQUAL(names.size(), 2)

}
END_SECTION

START_SECTION((const ToolMapping& getNameMapping()))
{
  INIUpdater i;
  ToolMapping m = i.getNameMapping();

  TEST_NOT_EQUAL(m.size(), 0)
  TEST_EQUAL(m[Internal::ToolDescriptionInternal("FeatureFinder",ListUtils::create<String>("centroided"))]
             == Internal::ToolDescriptionInternal("FeatureFinderCentroided",ListUtils::create<String>("")), true)

}
END_SECTION

START_SECTION((bool getNewToolName(const String &old_name, const String &tools_type, String &new_name)))
{
  INIUpdater i;
  String new_name;
  i.getNewToolName("FeatureFinder", "centroided", new_name);
  TEST_EQUAL(new_name, "FeatureFinderCentroided");

  i.getNewToolName("PeakPicker", "wavelet", new_name);
  TEST_EQUAL(new_name, "PeakPickerWavelet");

  i.getNewToolName("FileInfo", "", new_name);
  TEST_EQUAL(new_name, "FileInfo");

  i.getNewToolName("FileInfo", "bogus type", new_name); // type will be ignored - ok
  TEST_EQUAL(new_name, "FileInfo");

  TEST_EQUAL(i.getNewToolName("UNKNOWNTOOL", "bogus type", new_name), false);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



