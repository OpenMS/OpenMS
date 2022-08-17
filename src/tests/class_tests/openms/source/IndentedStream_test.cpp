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
#include <OpenMS/FORMAT/IndentedStream.h>
///////////////////////////

#include <OpenMS/CONCEPT/Colorizer.h>

#include <sstream>

using namespace OpenMS;
using namespace std;


const int TEST_CONSOLE_WIDTH = 9;
namespace OpenMS
{
  struct ConsoleWidthTest {
    ConsoleWidthTest()
    {
      auto& t = ConsoleUtils::getInstance(); // make sure the singleton is initialized
      const_cast<ConsoleUtils&>(t).console_width_ = TEST_CONSOLE_WIDTH;
    }
  };

  ConsoleWidthTest instance; // set console width!
} // namespace OpenMS

START_TEST(IndentedStream, "$Id$")


// test this first, because all the other tests rely on it
START_SECTION([EXTRA] int getConsoleWidth() const)
{
  auto& t = ConsoleUtils::getInstance();
  TEST_EQUAL(t.getConsoleWidth(), TEST_CONSOLE_WIDTH)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(IndentedStream(std::ostream& stream, const UInt indentation, const UInt max_lines))
{
	NOT_TESTABLE // tested below
}
END_SECTION

START_SECTION(IndentedStream& operator<<(Colorizer& colorizer))
{
  { // inline color
    stringstream ss; // will contain ANSI color codes... but they are not counted as characters on the current line
    IndentedStream is(ss, 3, 10);
    is << "12" << red("red") << "6789ab";
    TEST_EQUAL(ss.str(), "12\033[91mred\033[39m6789\n   ab") // the first line has more than #TEST_CONSOLE_WIDTH chars, but the ANSI codes do not count
  }
  
  { // color until revoked
    stringstream ss;
    IndentedStream is(ss, 3, 10);
    is << "12" << red() << "red" << red.undo() << "6789ab";
    TEST_EQUAL(ss.str(), "12\033[91mred\033[39m6789\n   ab") // the first line has more than #TEST_CONSOLE_WIDTH chars, but the ANSI codes do not count
  }
}
END_SECTION

START_SECTION(IndentedStream& operator<<(IndentedStream& self))
{
  NOT_TESTABLE // tested below
}
END_SECTION


const String x20(TEST_CONSOLE_WIDTH * 2 + 1, 'x'); // test string (2 full lines plus one 'x')
const String xC(TEST_CONSOLE_WIDTH, 'x');          // full console width of 'x'

START_SECTION((template<typename T> IndentedStream & operator<<(const T& data)))
{
  stringstream ss;
  int indent = 3;
  String s_indent(indent, ' ');
  IndentedStream is(ss, indent, 10);
  is << 3 << " " << '\n' << 'c';
  TEST_EQUAL(ss.str(), "3 \n" + s_indent + "c")
}
END_SECTION

START_SECTION(IndentedStream& operator<<(StreamManipulator manip))
{
  stringstream ss;
  IndentedStream is(ss, 3, 10);
  is << "xx" << std::endl << 'y';
  TEST_EQUAL(ss.str(), "xx\ny")
}
END_SECTION


START_SECTION(IndentedStream& indent(const UInt new_indent))
{
  stringstream ss;
  int indent = 3;
  int indent_new = 5;
  String s_indent(indent, ' ');
  String s_indent_new(indent_new, ' ');
  IndentedStream is(ss, indent, 10);
  is << xC  // a full line
     << 'y' // indented 'y'
     << is.indent(indent_new) // 'announce' that we want a new indentation for the next line
     << xC; // a full line which does not fit and triggers a linebreak
  TEST_EQUAL(ss.str(), xC  + "\n" 
                       + s_indent + "y" + (xC.data() + indent + 1) + "\n" // +1 to skip one 'x'
                         + s_indent_new + xC.suffix(indent + 1))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
