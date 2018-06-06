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
#include <OpenMS/APPLICATIONS/ToolHandler.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ToolHandler, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ToolHandler* ptr = nullptr;
ToolHandler* null_ptr = nullptr;
START_SECTION(ToolHandler())
{
	ptr = new ToolHandler();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ToolHandler())
{
	delete ptr;
}
END_SECTION

START_SECTION((static ToolListType getTOPPToolList(const bool includeGenericWrapper=false)))
{
  ToolListType list = ToolHandler::getTOPPToolList();
  TEST_EQUAL(list.has("FeatureFinderMRM"), true)
  TEST_EQUAL(list.has("GenericWrapper"), false)
  TEST_EQUAL(list.size() > 30, true)  // assume we have over 30 tools in there
  list = ToolHandler::getTOPPToolList(true);
  TEST_EQUAL(list.has("FeatureFinderMRM"), true)
  TEST_EQUAL(list.has("GenericWrapper"), true)
  TEST_EQUAL(list.size() > 30, true)  // assume we have over 30 tools in there
}
END_SECTION

START_SECTION((static ToolListType getUtilList()))
{
  ToolListType list = ToolHandler::getUtilList();
  TEST_EQUAL(list.has("SemanticValidator"), true)
  TEST_EQUAL(list.has("FFEval"), true)
  TEST_EQUAL(list.size() > 10, true)  // assume we have over 10 tools in there
}
END_SECTION

START_SECTION((static StringList getTypes(const String &toolname)))
{
  TEST_EQUAL(ToolHandler::getTypes("IsobaricAnalyzer") == StringList(), true);
  TEST_EQUAL(ToolHandler::getTypes("IDMapper") == StringList(), true);
}
END_SECTION

START_SECTION((static String getExternalToolsPath()))
{
  TEST_NOT_EQUAL(ToolHandler::getExternalToolsPath(), String())
}
END_SECTION

START_SECTION((static String getInternalToolsPath()))
{
  TEST_NOT_EQUAL(ToolHandler::getExternalToolsPath(), String())
}
END_SECTION

START_SECTION((static String getCategory(const String &toolname)))
{
  TEST_EQUAL(ToolHandler::getCategory("PepNovoAdapter"), "Identification")
  TEST_EQUAL(ToolHandler::getCategory("DOESNOTEXIST"), "")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



