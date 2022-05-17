// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2021.
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
// $Maintainer: $
// $Authors: Moritz Berger, Tetana Krymovska$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/CONCEPT/Colorizer.h>
///////////////////////////

#include <iostream>
#include <vector>
#include <string>
#include <sstream>
#include <fstream>

using namespace OpenMS;
using namespace std;

/* TO-DO
+ Check if string inputs are unchanged, and ANSI codes inside/not inside
- Windows test - define, copy, TEST if works
+/- Test every color function for all variable sets
- Test every other public function
- Test background color change test case - DEFINE functionality first
- ConsoleUtils Erweiterung
- mention "NOT TESTABLE" methods in sections
*/


START_TEST(Colorizer(),"$Id$")

 //Test variables
char tchar = 'a';
unsigned char unsignedChar = 'i';
signed char signedChar = 's';
short int shortInt = 32766;
unsigned short int unsShortInt = 600;
signed short int signShortInt = -32560;
unsigned long int unsLongInt = 40000000000;
long long int longLongInt = 981278728478274;
unsigned long long int unsLongLongInt = -78273829375;
float flt = 2094.5892;
double dbl = -253575634.98925;
long double longDbl= 135315.2929849375;
wchar_t wideChar = L'S';

/*
//////Functions to replace manual creation of stringstreams
//////and comparison strings for TEST_EQUAL for coloured outputs for Linux
//////
//////UNFINISHED!!!

template<typename T>
stringstream colorStream(T const& testVariable, Colorizer)
//saves a colored output of an instance of Colorizer into a stream 
{
    stringstream coloredStream;
    stringstream outputStream;

    coloredStream << Colorizer(testVariable);
    coloredStream >> outputStream;

    return outputStream; 
}

string createComparisonANSIString(string testVariable,Colorizer){
    //creates a model string with according ANSI codes to be
    //compared to stream created by colorStream
    
    string comparisonString;   
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(tchar)); 
    //to_string is causing formatting issues with char variables,
    // -> improve
    comparisonString.append(resetColorANSI);

  return conparisonString;
}

*/

START_SECTION(Colorizer::colorStream(ostream& stream) const)
{
    NOT_TESTABLE //is testable?
}
END_SECTION

START_SECTION(outputToStream(ostream& o_stream))
{
    NOT_TESTABLE //is testable?
}
END_SECTION

START_SECTION(resetColor(ostream& stream))
{
    NOT_TESTABLE //is testable?
}
END_SECTION



//    string asciiString = " !#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

#if defined(_WIN32) || defined(_WIN32) || defined(__WIN32__) || defined(_WIN64)

//Check that the colorized input contains the original text and DOES NOT contain ASCI codes
START_SECTION(Colorizer)
{
    string asciiString = "!#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

    stringstream blackStream;
    stringstream redStream;
    stringstream greenStream;
    stringstream yellowStream;
    stringstream blueStream;
    stringstream magentaStream;
    stringstream cyanStream;
    stringstream whiteStream;

    string outStream;

    blackStream << black(asciiString);
    blackStream >> outStream;
    TEST_EQUAL(outStream, asciiString)

    redStream << red(asciiString);
    redStream >> outStream;
    TEST_EQUAL(outStream, asciiString)

    greenStream << green(asciiString);
    greenStream >> outStream;
    TEST_EQUAL(outStream, asciiString)

    yellowStream << yellow(asciiString);
    yellowStream >> outStream;
    TEST_EQUAL(outStream, asciiString)

    blueStream << blue(asciiString);
    blueStream >> outStream;
    TEST_EQUAL(outStream, asciiString)

    magentaStream << magenta(asciiString);
    magentaStream >> outStream;
    TEST_EQUAL(outStream, asciiString)

    cyanStream << cyan(asciiString);
    cyanStream >> outStream;
    TEST_EQUAL(outStream, asciiString)

    whiteStream << white(asciiString);
    whiteStream >> outStream;
    TEST_EQUAL(outStream, asciiString)
}
END_SECTION

#elif defined(__linux__) || defined(__OSX__)

//ANSI codes
string blackANSI        = "\e[30m";
string redANSI          = "\e[31m";
string greenANSI        = "\e[32m";
string yellowANSI       = "\e[33m";
string blueANSI         = "\e[34m";
string magentaANSI      = "\e[35m";
string cyanANSI         = "\e[36m";;
string whiteANSI        = "\e[37m";
string resetColorANSI   = "\e[0m";

 //Check that the colorized input contains the original text and according ASCI codes
START_SECTION(Colorizer)
{
    string asciiString = "!#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";

    stringstream blackStream;
    stringstream redStream;
    stringstream greenStream;
    stringstream yellowStream;
    stringstream blueStream;
    stringstream magentaStream;
    stringstream cyanStream;
    stringstream whiteStream;

    string outStream;

    blackStream << black(asciiString);
    blackStream >> outStream;
    TEST_EQUAL(outStream, "\e[30m"+asciiString+"\e[0m")

    redStream << red(asciiString);
    redStream >> outStream;
    TEST_EQUAL(outStream, "\e[31m"+asciiString+"\e[0m")

    greenStream << green(asciiString);
    greenStream >> outStream;
    TEST_EQUAL(outStream, "\e[32m"+asciiString+"\e[0m")

    yellowStream << yellow(asciiString);
    yellowStream >> outStream;
    TEST_EQUAL(outStream, "\e[33m"+asciiString+"\e[0m")

    blueStream << blue(asciiString);
    blueStream >> outStream;
    TEST_EQUAL(outStream, "\e[34m"+asciiString+"\e[0m")

    magentaStream << magenta(asciiString);
    magentaStream >> outStream;
    TEST_EQUAL(outStream, "\e[35m"+asciiString+"\e[0m")

    cyanStream << cyan(asciiString);
    cyanStream >> outStream;
    TEST_EQUAL(outStream, "\e[36m"+asciiString+"\e[0m")

    whiteStream << white(asciiString);
    whiteStream >> outStream;
    TEST_EQUAL(outStream, "\e[37m"+asciiString+"\e[0m")

}
END_SECTION

//testing various inputs for colorizing
START_SECTION(Colorizer::black())
{

    string coloredStream;
    string comparisonString;
    stringstream blackStream;

    //char////////////////////////////////////
    blackStream << black(tchar);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(tchar));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //unsignedChar////////////////////////////
    blackStream << black(unsignedChar);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(unsignedChar));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //signed char///////////////////////////////
    blackStream << black(signedChar);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(signedChar));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //wide char///////////////////////////////
    blackStream << black(wideChar);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(wideChar));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //short int///////////////////////////////
    blackStream << black(shortInt);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(shortInt));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //unsigned short int///////////////////////////////
    blackStream << black(unsShortInt);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(unsShortInt));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //signed short int///////////////////////////////
    blackStream << black(signShortInt);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(signShortInt));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //unsigned long int///////////////////////////////
    blackStream << black(unsLongInt);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(unsLongInt));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //long long int///////////////////////////////
    blackStream << black(longLongInt);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(longLongInt));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //unsigned long long int///////////////////////////////
    blackStream << black(unsLongLongInt);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(unsLongLongInt));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //float///////////////////////////////
    blackStream << black(flt);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(flt));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //double///////////////////////////////
    blackStream << black(dbl);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(dbl));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();

    //long double///////////////////////////////
    blackStream << black(longDbl);
    blackStream >> coloredStream;
    comparisonString.append(blackANSI);
    comparisonString.append(to_string(longDbl));
    comparisonString.append(resetColorANSI);
    TEST_EQUAL(coloredStream, comparisonString)

    //clearing streams
    coloredStream.clear();
    comparisonString = "";
    blackStream.clear();



}
END_SECTION

//test case - resetting color only for given stream - implement?
//Or is testing of only class methods sufficient?


#endif


END_TEST

