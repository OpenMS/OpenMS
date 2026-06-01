// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/Colorizer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;


START_TEST(Colorizer, "$Id$")

// Test variables
int test_int = 15;
//string test_string = " !#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

// ANSI codes
string const redANSI = "\033[91m";
string const greenANSI = "\033[92m";
string const yellowANSI = "\033[93m";
string const blueANSI = "\033[94m";
string const magentaANSI = "\033[95m";
string const cyanANSI = "\033[96m";
string const resetColorANSI = "\033[39m";


START_SECTION(Colorizer(const ConsoleColor color))
{
  Colorizer test(ConsoleColor::BLUE);
  stringstream s;
  s << "-" << test("test") << "-";
  TEST_EQUAL(s.str(), String("-") + blueANSI + "test" + resetColorANSI + "-")
}
END_SECTION

START_SECTION(~Colorizer())
{
  NOT_TESTABLE
}
END_SECTION

START_SECTION(Colorizer& operator()())
{
  Colorizer test(ConsoleColor::BLUE);
  stringstream s;
  s << "-" << test() << "test-" << test_int;
  TEST_EQUAL(s.str(), String("-") + blueANSI + "test-" + test_int)
}
END_SECTION

START_SECTION(template<typename T> Colorizer& operator()(T s))
{
  Colorizer test(ConsoleColor::MAGENTA);
  stringstream s;
  s << "-" << test(test_int) << "test-" << test_int;
  TEST_EQUAL(s.str(), String("-") + magentaANSI + test_int + resetColorANSI + "test-" + test_int)
}
END_SECTION

START_SECTION(Colorizer& undo())
{
  Colorizer test(ConsoleColor::CYAN);
  {
    stringstream s;
    s << "-" << test() << "test" << test_int << test.undo() << "nocol";
    TEST_EQUAL(s.str(), String("-") + cyanANSI + "test" + test_int + resetColorANSI + "nocol")
  }
  
  // test double coloring + reset (using any Colorizer)
  Colorizer yellow_test(ConsoleColor::YELLOW);
  stringstream s;
  s << "-" << test() << "test" << yellow_test() << test_int << yellow_test.undo() << "nocol";
  TEST_EQUAL(s.str(), String("-") + cyanANSI + "test" + yellowANSI + test_int + resetColorANSI + "nocol")
}
END_SECTION

START_SECTION([EXTRA] visual inspection(only works when not redirecting cout/cerr to file))
{
  std::cerr << "\n\n --- for COUT ---\n" << std::flush;
  std::cout << red.undoAll();
  std::cout << "\n-" << red("red inline text") << " -" << red() << " red to infinity " << red.undo();
  std::cout << "\n-" << green("green inline  text") << " -" << green() << " green to infinity " << green.undo();
  std::cout << "\n-" << yellow("yellow inline  text") << " -" << yellow() << " yellow to infinity " << yellow.undo();
  std::cout << "\n-" << blue("blue inline  text") << " -" << blue() << " blue to infinity " << blue.undo();
  std::cout << "\n-" << magenta("magenta inline  text") << " -" << magenta() << " magenta to infinity " << magenta.undo();
  std::cout << "\n-" << cyan("cyan inline  text") << " -" << cyan() << " cyan to infinity " << cyan.undo();
  std::cout << "\n-" << invert("inverted inline  text") << " -" << invert() << " invert to infinity " << invert.undo();
  std::cout << "\n-" << bright("bright inline  text") << " -" << bright() << " bright to infinity " << bright.undo();
  std::cout << "\n-" << underline("underline inline  text") << " -" << underline() << " underline to infinity " << underline.undo();
  std::cout << "\n-" << underline() << bright() << green() <<" underlined, bright, green to infinity " << underline.undoAll();
  std::cout << std::flush; // make sure the ANSI code made it to the stream

  std::cerr << "\n\n --- for CERR ---\n";
  std::cerr << red.undoAll();
  std::cerr << "\n-" << red("red inline text") << " -" << red() << " red to infinity " << red.undo();
  std::cerr << "\n-" << green("green inline  text") << " -" << green() << " green to infinity " << green.undo();
  std::cerr << "\n-" << yellow("yellow inline  text") << " -" << yellow() << " yellow to infinity " << yellow.undo();
  std::cerr << "\n-" << blue("blue inline  text") << " -" << blue() << " blue to infinity " << blue.undo();
  std::cerr << "\n-" << magenta("magenta inline  text") << " -" << magenta() << " magenta to infinity " << magenta.undo();
  std::cerr << "\n-" << cyan("cyan inline  text") << " -" << cyan() << " cyan to infinity " << cyan.undo();
  std::cerr << "\n-" << invert("inverted inline  text") << " -" << invert() << " invert to infinity " << invert.undo();
  std::cerr << "\n-" << bright("bright inline  text") << " -" << bright() << " bright to infinity " << bright.undo();
  std::cerr << "\n-" << underline("underline inline  text") << " -" << underline() << " underline to infinity " << underline.undo();
  std::cerr << "\n-" << underline() << bright() << green() << " underlined, bright, green to infinity " << underline.undoAll();
  std::cerr << std::flush; // make sure the ANSI code made it to the stream
}
END_SECTION

END_TEST
