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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/METADATA/MetaInfoRegistry.h>

///////////////////////////

START_TEST(MetaInfoRegistry, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

MetaInfoRegistry* test = nullptr;
MetaInfoRegistry* nullPointer = nullptr;
START_SECTION((MetaInfoRegistry()))
	test = new MetaInfoRegistry;
  TEST_NOT_EQUAL(test, nullPointer)
END_SECTION


START_SECTION((~MetaInfoRegistry()))
	delete test;
END_SECTION

MetaInfoRegistry mir;

START_SECTION((UInt registerName(const String& name, const String& description = "", const String& unit = "")))
{
	UInt testname = mir.registerName("testname", "this is just a test");
  TEST_EQUAL(testname, 1024);
	UInt retention_time = mir.registerName("retention time", "this is just another test", "sec");
  TEST_EQUAL(retention_time, 1025);
	UInt another_testname = mir.registerName("another testname", "i will be set later", "me too");
  TEST_EQUAL(another_testname, 1026);
}
END_SECTION

START_SECTION((void setDescription(UInt index, const String& description)))
	mir.setDescription(1026, "foo");
	TEST_STRING_EQUAL(mir.getDescription(1026), "foo")
END_SECTION

START_SECTION((void setDescription(const String& name, const String& description)))
	mir.setDescription("another testname", "bar");
	TEST_STRING_EQUAL(mir.getDescription(1026), "bar")
END_SECTION

START_SECTION((void setUnit(UInt index, const String& unit)))
	mir.setUnit(1026, "foo");
	TEST_STRING_EQUAL(mir.getUnit(1026), "foo")
END_SECTION

START_SECTION((void setUnit(const String& name, const String& unit)))
	mir.setUnit("another testname", "bar");
	TEST_STRING_EQUAL(mir.getUnit(1026), "bar")
END_SECTION

START_SECTION((UInt getIndex(const String& name) const))
	TEST_EQUAL(mir.getIndex("testname"), 1024)
	TEST_EQUAL(mir.getIndex("retention time"), 1025)
	TEST_EQUAL(mir.getIndex("isotopic_range"), 1)
	TEST_EQUAL(mir.getIndex("cluster_id"), 2)
  TEST_EQUAL(mir.getIndex("unregistered name"), UInt(-1))
END_SECTION

START_SECTION((String getName(UInt index) const))
	TEST_STRING_EQUAL(mir.getName(1), "isotopic_range")
	TEST_STRING_EQUAL(mir.getName(2), "cluster_id")
	TEST_STRING_EQUAL(mir.getName(3), "label")
	TEST_STRING_EQUAL(mir.getName(4), "icon")
	TEST_STRING_EQUAL(mir.getName(1024), "testname")
	TEST_STRING_EQUAL(mir.getName(1025), "retention time")
END_SECTION

START_SECTION((String getDescription(UInt index) const))
	TEST_STRING_EQUAL(mir.getDescription(1024), "this is just a test")
	TEST_STRING_EQUAL(mir.getDescription(1025), "this is just another test")
	TEST_STRING_EQUAL(mir.getDescription(1), "consecutive numbering of the peaks in an isotope pattern. 0 is the monoisotopic peak")
	TEST_STRING_EQUAL(mir.getDescription(2), "consecutive numbering of isotope clusters in a spectrum")
END_SECTION

START_SECTION((String getDescription(const String& name) const))
	TEST_STRING_EQUAL(mir.getDescription("testname"), "this is just a test")
	TEST_STRING_EQUAL(mir.getDescription("retention time"), "this is just another test")
	TEST_STRING_EQUAL(mir.getDescription("isotopic_range"), "consecutive numbering of the peaks in an isotope pattern. 0 is the monoisotopic peak")
	TEST_STRING_EQUAL(mir.getDescription("cluster_id"), "consecutive numbering of isotope clusters in a spectrum")
END_SECTION

START_SECTION((String getUnit(UInt index) const))
	TEST_STRING_EQUAL(mir.getUnit(1024), "")
	TEST_STRING_EQUAL(mir.getUnit(1025), "sec")
	TEST_STRING_EQUAL(mir.getUnit(1), "")
	TEST_STRING_EQUAL(mir.getUnit(2), "")
END_SECTION

START_SECTION((String getUnit(const String& name) const))
	TEST_STRING_EQUAL(mir.getUnit("testname"), "")
  TEST_STRING_EQUAL(mir.getUnit("retention time"), "sec")
	TEST_STRING_EQUAL(mir.getUnit("isotopic_range"), "")
	TEST_STRING_EQUAL(mir.getUnit("cluster_id"), "")
END_SECTION

START_SECTION((MetaInfoRegistry(const MetaInfoRegistry& rhs)))
	MetaInfoRegistry mir2(mir);
  TEST_EQUAL(mir2.getIndex("testname"), 1024)
  TEST_EQUAL(mir2.getIndex("retention time"), 1025)
	TEST_STRING_EQUAL(mir2.getName(1), "isotopic_range")
	TEST_STRING_EQUAL(mir2.getName(1024), "testname")
	TEST_STRING_EQUAL(mir2.getName(1025), "retention time")
	TEST_STRING_EQUAL(mir2.getDescription(1024), "this is just a test")
	TEST_STRING_EQUAL(mir2.getDescription(1025), "this is just another test")
	TEST_STRING_EQUAL(mir2.getDescription("testname"), "this is just a test")
	TEST_STRING_EQUAL(mir2.getDescription("retention time"), "this is just another test")
	TEST_STRING_EQUAL(mir2.getUnit(1024), "")
	TEST_STRING_EQUAL(mir2.getUnit(1025), "sec")
	TEST_STRING_EQUAL(mir2.getUnit("testname"), "")
	TEST_STRING_EQUAL(mir2.getUnit("retention time"), "sec")
END_SECTION

START_SECTION((MetaInfoRegistry& operator=(const MetaInfoRegistry& rhs)))
	MetaInfoRegistry mir2;
	mir2 = mir;
  TEST_EQUAL(mir2.getIndex("testname"), 1024)
  TEST_EQUAL(mir2.getIndex("retention time"), 1025)
	TEST_STRING_EQUAL(mir2.getName(1), "isotopic_range")
	TEST_STRING_EQUAL(mir2.getName(1024), "testname")
	TEST_STRING_EQUAL(mir2.getName(1025), "retention time")
	TEST_STRING_EQUAL(mir2.getDescription(1024), "this is just a test")
	TEST_STRING_EQUAL(mir2.getDescription(1025), "this is just another test")
	TEST_STRING_EQUAL(mir2.getDescription("testname"), "this is just a test")
	TEST_STRING_EQUAL(mir2.getDescription("retention time"), "this is just another test")
	TEST_STRING_EQUAL(mir2.getUnit(1024), "")
	TEST_STRING_EQUAL(mir2.getUnit(1025), "sec")
	TEST_STRING_EQUAL(mir2.getUnit("testname"), "")
	TEST_STRING_EQUAL(mir2.getUnit("retention time"), "sec")
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
