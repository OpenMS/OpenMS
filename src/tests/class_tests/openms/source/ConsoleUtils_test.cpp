// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors:  Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/APPLICATIONS/ConsoleUtils.h>
///////////////////////////


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

START_TEST(ConsoleUtils, "$Id$")


// test this first, because all the other tests rely on it
START_SECTION(int getConsoleWidth() const)
{
  auto& t = ConsoleUtils::getInstance();
  TEST_EQUAL(t.getConsoleWidth(), TEST_CONSOLE_WIDTH)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(ConsoleUtils())
{
	NOT_TESTABLE // private
}
END_SECTION

START_SECTION(~ConsoleUtils())
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(static const ConsoleUtils& getInstance())
{
  NOT_TESTABLE // tested below
}
END_SECTION


const String x20(TEST_CONSOLE_WIDTH * 2 + 1, 'x'); // test string (2 full lines plus one 'x')
const String xC(TEST_CONSOLE_WIDTH, 'x');          // full console width of 'x'

START_SECTION((static StringList breakStringList(const String& input, const Size indentation, const Size max_lines, const Size first_line_prefill = 0)))
{
  // we actually test the concatenation using breakString() since its easier to write
  String broken_string;
  { // test with indent = 0
    broken_string = ConsoleUtils::breakString(x20, 0, 10);
    TEST_EQUAL(broken_string, xC + '\n' + xC + '\n' + "x")
  }
  { // try again ...
    broken_string = ConsoleUtils::breakString(x20, 0, 10);
    TEST_EQUAL(broken_string, xC + '\n' + xC + '\n' + "x")
  }
  { // test with indent = 3
    int indent = 3;
    String shortX(TEST_CONSOLE_WIDTH - indent, 'x');
    String s_indent(indent, ' ');
    broken_string = ConsoleUtils::breakString(x20, indent, 10);
    TEST_EQUAL(broken_string, xC + '\n' + s_indent + shortX + '\n' + s_indent + "xxxx")
  }
  { // test with prefilled first line
    int indent = 3;
    int prefill = 5;
    String firstX(TEST_CONSOLE_WIDTH - prefill, 'x');
    String shortX(TEST_CONSOLE_WIDTH - indent, 'x');
    String s_indent(indent, ' ');
    broken_string = ConsoleUtils::breakString(x20, indent, 10, prefill);
    TEST_EQUAL(broken_string, firstX + '\n'                 // 4x
                              + s_indent + shortX + '\n'    //+6x
                              + s_indent + shortX + '\n'    //+6x
                              + s_indent + "xxx")           //+3x
  }
  { // test with manual linebreaks in between
    { // just a linebreak
      int indent = 0;
      int prefill = 0;
      broken_string = ConsoleUtils::breakString("\n", indent, 10, prefill);
      TEST_EQUAL(broken_string, '\n')
    }
    { // just a linebreak with indent
      int indent = 3;
      int prefill = 0;
      broken_string = ConsoleUtils::breakString("\n", indent, 10, prefill);
      TEST_EQUAL(broken_string, '\n' + String(indent, ' '))
    }
    { // prefilled linebreak with indent (should not make a difference)
      int indent = 3;
      int prefill = 5;
      broken_string = ConsoleUtils::breakString("\n", indent, 10, prefill);
      TEST_EQUAL(broken_string, '\n' + String(indent, ' '))
    }
    { // text with a linebreak with indent and prefill
      int indent = 3;
      int prefill = TEST_CONSOLE_WIDTH - 1; // one char left on first line
      broken_string = ConsoleUtils::breakString("xxx\n", indent, 10, prefill);
      TEST_EQUAL(broken_string, "x\n" + String(indent, ' ') + "xx\n" + String(indent, ' '))
    }
    
    { // some corner cases (only one char per line)
      int indent = TEST_CONSOLE_WIDTH - 1;
      int prefill = indent; // one char left on first line
      broken_string = ConsoleUtils::breakString("xxx\n", indent, 10, prefill);
      TEST_EQUAL(broken_string, "x\n" 
                                + String(indent, ' ') + "x\n" 
                                + String(indent, ' ') + "x\n" 
                                + String(indent, ' ') + "\n"  // manual linebreak right after breakString()'s linebreak -- you get two lines ...
                                + String(indent, ' '))
    }

    
  }
  { // test max_lines
    int indent = TEST_CONSOLE_WIDTH - 2;
    int prefill = indent;  // two chars per EVERY line
    broken_string = ConsoleUtils::breakString(String(99, 'x'), indent, 3, prefill);
    TEST_EQUAL(broken_string, "xx\n" + String(indent, ' ') + "...\n" + String(indent, ' ') + 'x')
  }

}
END_SECTION

START_SECTION(static String breakString(const String& input, const Size indentation, const Size max_lines, const Size first_line_prefill = 0))
{
  NOT_TESTABLE // tested above
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



