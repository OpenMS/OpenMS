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
