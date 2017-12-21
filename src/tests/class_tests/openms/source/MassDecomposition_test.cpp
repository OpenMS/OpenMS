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
#include <OpenMS/CHEMISTRY/MASSDECOMPOSITION/MassDecomposition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MassDecomposition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MassDecomposition* ptr = nullptr;
MassDecomposition* nullPointer = nullptr;
START_SECTION(MassDecomposition())
{
	ptr = new MassDecomposition();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~MassDecomposition())
{
	delete ptr;
}
END_SECTION

START_SECTION((MassDecomposition(const MassDecomposition &deco)))
{
  MassDecomposition md("C3 M4 S200");
	TEST_EQUAL(md.getNumberOfMaxAA(), 200)
	TEST_STRING_EQUAL(md.toString(), "C3 M4 S200")

	MassDecomposition md2(md);
	TEST_EQUAL(md2.getNumberOfMaxAA(), 200)
	TEST_STRING_EQUAL(md2.toString(), "C3 M4 S200")
}
END_SECTION

START_SECTION((MassDecomposition(const String &deco)))
{
  MassDecomposition md("C3 M4 S200");
	TEST_EQUAL(md.getNumberOfMaxAA(), 200)
	TEST_STRING_EQUAL(md.toString(), "C3 M4 S200")
}
END_SECTION

START_SECTION((MassDecomposition& operator=(const MassDecomposition &rhs)))
{
  MassDecomposition md("C3 M4 S200");
	MassDecomposition md2;

	md2 = md;
	TEST_EQUAL(md2.getNumberOfMaxAA(), 200)
	TEST_STRING_EQUAL(md2.toString(), "C3 M4 S200")
}
END_SECTION

START_SECTION((MassDecomposition& operator+=(const MassDecomposition &d)))
{
  MassDecomposition md;
	MassDecomposition md1("C3");
	MassDecomposition md2("M4");
	MassDecomposition md3("S200");
	md += md1;
	TEST_EQUAL(md.getNumberOfMaxAA(), 3)
	TEST_STRING_EQUAL(md.toString(), "C3")
	md += md2;
	TEST_EQUAL(md.getNumberOfMaxAA(), 4)
	TEST_STRING_EQUAL(md.toString(), "C3 M4")
	md += md3;
	TEST_EQUAL(md.getNumberOfMaxAA(), 200)
	TEST_STRING_EQUAL(md.toString(), "C3 M4 S200")
}
END_SECTION

START_SECTION((String toString() const ))
{
  MassDecomposition md1("C3");
  MassDecomposition md2("C3 M4");
  MassDecomposition md3("C3 M4 S200");
  
	TEST_EQUAL(md1.getNumberOfMaxAA(), 3)
  TEST_STRING_EQUAL(md1.toString(), "C3")
  
	TEST_EQUAL(md2.getNumberOfMaxAA(), 4)
  TEST_STRING_EQUAL(md2.toString(), "C3 M4")
 
 	TEST_EQUAL(md3.getNumberOfMaxAA(), 200)
  TEST_STRING_EQUAL(md3.toString(), "C3 M4 S200")
}
END_SECTION

START_SECTION((String toExpandedString() const ))
{
  MassDecomposition md1("C3");
  MassDecomposition md2("C3 M4");

  TEST_STRING_EQUAL(md1.toExpandedString(), "CCC")
  TEST_STRING_EQUAL(md2.toExpandedString(), "CCCMMMM")
}
END_SECTION

START_SECTION((MassDecomposition operator+(const MassDecomposition &rhs) const ))
{
  MassDecomposition md;
  MassDecomposition md1("C3");
  MassDecomposition md2("M4");
  MassDecomposition md3("S200");
  MassDecomposition md5 = md + md1;
  TEST_EQUAL(md5.getNumberOfMaxAA(), 3)
  TEST_STRING_EQUAL(md5.toString(), "C3")
  
	MassDecomposition md6 = md1 + md2;
  TEST_EQUAL(md6.getNumberOfMaxAA(), 4)
  TEST_STRING_EQUAL(md6.toString(), "C3 M4")

  MassDecomposition md7 = md1 + md2 + md3;
  TEST_EQUAL(md7.getNumberOfMaxAA(), 200)
  TEST_STRING_EQUAL(md7.toString(), "C3 M4 S200")
}
END_SECTION

START_SECTION((Size getNumberOfMaxAA() const ))
{
  MassDecomposition md;
  MassDecomposition md1("C3");
  MassDecomposition md2("M4");
  MassDecomposition md3("S200");
  md += md1;
  TEST_EQUAL(md.getNumberOfMaxAA(), 3)
  TEST_STRING_EQUAL(md.toString(), "C3")
  md += md2;
  TEST_EQUAL(md.getNumberOfMaxAA(), 4)
  TEST_STRING_EQUAL(md.toString(), "C3 M4")
  md += md3;
  TEST_EQUAL(md.getNumberOfMaxAA(), 200)
  TEST_STRING_EQUAL(md.toString(), "C3 M4 S200")
}
END_SECTION

START_SECTION((bool operator<(const MassDecomposition &rhs) const ))
{
  MassDecomposition md;
  MassDecomposition md1("C3");
  MassDecomposition md2("M4");
  md += md1;
	TEST_EQUAL(md2 < md1, false)
	TEST_EQUAL(md < md2, true)
  md += md2;
	TEST_EQUAL(md < md2, true)
}
END_SECTION

START_SECTION((bool operator==(const String &deco) const ))
{
	MassDecomposition md;
  MassDecomposition md1("C3");
	TEST_EQUAL(md == md1.toString(), false)
  MassDecomposition md2("M4");
	md = md2;
	TEST_STRING_EQUAL(md2.toString(), "M4")
	TEST_EQUAL(md == md2.toString(), true)
  MassDecomposition md3("S200");
	md = md2 + md3;
	TEST_EQUAL(md == md3.toString(), false)
}
END_SECTION

START_SECTION((bool containsTag(const String &tag) const ))
{
  MassDecomposition md;
  MassDecomposition md1("C3");
  MassDecomposition md2("C3 M4");
  MassDecomposition md3("C3 M4 S200");
	TEST_EQUAL(md.containsTag("C"), false)
	TEST_EQUAL(md.containsTag("CCC"), false)
	TEST_EQUAL(md1.containsTag("CCC"), true)
	TEST_EQUAL(md1.containsTag("CCCC"), false)
	TEST_EQUAL(md2.containsTag("CMC"), true)
	TEST_EQUAL(md3.containsTag("CCCSSMSSSSSSSSSSSSSSM"), true)
}
END_SECTION

START_SECTION((bool compatible(const MassDecomposition &deco) const ))
{
  MassDecomposition md;
  MassDecomposition md1("C3");
  MassDecomposition md2("C3 M4");
  MassDecomposition md3("C3 M4 S200");
  MassDecomposition md4("M4 S200");
  MassDecomposition md5("C3 S200");
  MassDecomposition md6("S2");
	TEST_EQUAL(md.compatible(MassDecomposition("")), true)
	TEST_EQUAL(md.compatible(MassDecomposition("C1")), false)
	TEST_EQUAL(md1.compatible(MassDecomposition("C1")), true)
	TEST_EQUAL(md2.compatible(MassDecomposition("C2 M4")), true)
	TEST_EQUAL(md2.compatible(MassDecomposition("C2 M5")), false)
	TEST_EQUAL(md3.compatible(md5), true)
	TEST_EQUAL(md3.compatible(md2), true)
	TEST_EQUAL(md3.compatible(md6), true)
	TEST_EQUAL(md3.compatible(md4), true)
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



