// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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
// $Authors:  Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/ConsoleUtils.h>
#include <OpenMS/CONCEPT/Colorizer.h>
///////////////////////////
#include <fstream>
#include <iostream>
#include <vector>
#include <string>
#include <sstream>

using namespace OpenMS;
using namespace std;

START_TEST(ConsoleUtils, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(const int getConsoleSize())
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(ConsoleUtils getInstance())
{
	NOT_TESTABLE
}
END_SECTION

#ifdef OPENMS_WINDOWSPLATTFORM

START_SECTION(void resetConsoleColor())
{
	ConsoleUtils c(ConsoleUtils const);
	c.ConsoleUtils::setCoutColor(16);
	c.ConsoleUtils::setCerrColor(9);
	int def_cout = c.ConsoleUtils::getCoutColor();
	int def_cerr = c.ConsoleUtils::getCerrColor();
	TEST_EQUAL(def_cout,16)
	TEST_EQUAL(def_cerr,10)
}
END_SECTION

START_SECTION(void setCoutColor())
{
	ConsoleUtils c(ConsoleUtils const);
	c.ConsoleUtils::setCoutColor(16);
	int def_ = c.ConsoleUtils::getCoutColor();
	TEST_EQUAL(def,16);
}
END_SECTION

START_SECTION(void setCerrColor())
{
	ConsoleUtils c(ConsoleUtils const);
	c.ConsoleUtils::setCerrColor(11);
	int def_ = c.ConsoleUtils::getCoutColor();
	TEST_EQUAL(def,11);
}
END_SECTION
#endif

START_SECTION((static OpenMS::String breakString(const String& input,
										const Size indentation, 
										const Size max_lines,
										const Size curser_pos = 0)))
{
	// we cannot predict which shape the broken string will have, so testing is rather limited
	// we cannot predict which shape the broken string will have, so testing is rather limited
	String test_string = "This is a test string which should be broken up into multiple lines.";
	String broken_string = ConsoleUtils::breakString(test_string, 0, 10);

	TEST_EQUAL(test_string.length() <= broken_string.length(), true)

}
END_SECTION

START_SECTION(static OpenMS::StringList breakStringList(const String& input,
										const Size indentation, 
										const Size max_lines,
										const Size curser_pos = 0))
{
	String test_string = "This is a test string which should be broken up into multiple lines.";
	OpenMS::StringList broken_string = ConsoleUtils::breakStringList(test_string, 0, 10);

	String broken_string_string = broken_string.StringList::at(0);
	TEST_EQUAL(test_string.length() <= broken_string_string.length(), true)
}
END_SECTION

END_TEST



