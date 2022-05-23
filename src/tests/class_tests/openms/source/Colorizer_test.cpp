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

/* TO-DO 1
+ Check if string inputs are unchanged, and ANSI codes inside/not inside
- Windows test - define, copy, TEST if works
+ Test every color function for all variable sets
- Test every other public function
- mention "NOT TESTABLE" methods in sections
-!!!ALL changes requested after last push
*/

/* TO-DO 2
- UPDATE OpenMS Version (consult) and Colorizer version
- Test background color change test case - DEFINE functionality first
- ConsoleUtils Erweiterung
*/

template <typename T>
//convert any input to string
string convertToString ( T var_input )
{
    ostringstream ss;
    ss << var_input;
    return ss.str();
}

START_TEST(Colorizer(),"$Id$")

 //Test variables
char test_char = 'a';
int test_int = 15;
float test_float = 2094.5892;
string test_string = " !#$%&\'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\\]^_`abcdefghijklmnopqrstuvwxyz{|}~";
// string test_string = "ABCDE";

//ANSI codes
string const blackANSI        = "\e[30m";
string const redANSI          = "\e[31m";
string const greenANSI        = "\e[32m";
string const yellowANSI       = "\e[33m";
string const blueANSI         = "\e[34m";
string const magentaANSI      = "\e[35m";
string const cyanANSI         = "\e[36m";;
string const whiteANSI        = "\e[37m";
string const resetColorANSI   = "\e[0m";



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

  return comparisonString;
}

*/

START_SECTION(Colorizer::colorStream(ostream& stream) const)
{
    NOT_TESTABLE //is testable?
}
END_SECTION

START_SECTION(Colorizer::outputToStream(ostream& o_stream))
{
    // Colorizer test_colorizer;
    // stringstream test_stream;
    // stringstream test_stream << "ABCDEFGHIJKLM";
    // test_colorizer.outputToStream(test_stream);
    // TEST_EQUAL(test_colortestizer,"ABCDEFGHIJKLM");


}
END_SECTION

START_SECTION(Colorizer::resetColor(ostream& stream))
{
    // stringstream test_stream;
    // Colorizer test_colorizer;

    // test_colorizer.resetColor();
    // test_colorizer >> test_stream;
    // TEST_EQUAL(test_stream, resetColorANSI)

    #if defined(__linux__) || defined(__OSX__)

    #endif 

}
END_SECTION

START_SECTION("Testing Colorizer instances")
{
    #ifdef OPENMS_WINDOWSPLATFORM 
    //Check that the colorized input contains the original text and DOES NOT contain ASCI codes
    //TEST ON Windows!S
        stringstream colored_stream;

        colored_stream << black(test_string);
        TEST_EQUAL(colored_stream.str(), test_string)
        colored_stream.clear();

        colored_stream << red(test_string);
        TEST_EQUAL(colored_stream.str(), test_string)
        colored_stream.clear();

        colored_stream << green(test_string);
        TEST_EQUAL(colored_stream.str(), test_string)
        colored_stream.clear();

        colored_stream << yellow(test_string);
        TEST_EQUAL(colored_stream.str(), test_string)
        colored_stream.clear();

        colored_stream << blue(test_string);
        TEST_EQUAL(colored_stream.str(), test_string)
        colored_stream.clear();

        colored_stream << magenta(test_string);
        TEST_EQUAL(colored_stream.str(), test_string)
        colored_stream.clear();

        colored_stream << cyan(test_string);
        TEST_EQUAL(colored_stream.str(), test_string)
        colored_stream.clear();

        colored_stream << white(test_string);
        TEST_EQUAL(colored_stream.str(), test_string)
        colored_stream.clear();

    #elif defined(__linux__) || defined(__OSX__)
        //LINUX//////LINUX////LINUX////LINUX////LINUX////LINUX////LINUX////
        //Check that the colorized input contains the original text and according ASCI codes

        stringstream colored_stream;

        colored_stream << black(test_string);
        TEST_EQUAL(colored_stream.str(), blackANSI+test_string+resetColorANSI)
        colored_stream.str(string());
        colored_stream.clear();

        colored_stream << red(test_string);
        TEST_EQUAL(colored_stream.str(), redANSI+test_string+resetColorANSI)
        colored_stream.str(string());
        colored_stream.clear();

        colored_stream << green(test_string);
        TEST_EQUAL(colored_stream.str(), greenANSI+test_string+resetColorANSI)
        colored_stream.str(string());
        colored_stream.clear();

        colored_stream << yellow(test_string);
        TEST_EQUAL(colored_stream.str(), yellowANSI+test_string+resetColorANSI)
        colored_stream.str(string());
        colored_stream.clear();

        colored_stream << blue(test_string);
        TEST_EQUAL(colored_stream.str(), blueANSI+test_string+resetColorANSI)
        colored_stream.str(string());
        colored_stream.clear();

        colored_stream << magenta(test_string);
        TEST_EQUAL(colored_stream.str(), magentaANSI+test_string+resetColorANSI)
        colored_stream.str(string());
        colored_stream.clear();

        colored_stream << cyan(test_string);
        TEST_EQUAL(colored_stream.str(), cyanANSI+test_string+resetColorANSI)
        colored_stream.str(string());
        colored_stream.clear();

        colored_stream << white(test_string);
        TEST_EQUAL(colored_stream.str(), whiteANSI+test_string+resetColorANSI)
        colored_stream.str(string());
        colored_stream.clear();


    #endif

}
END_SECTION

//testing various inputs for colorizing
START_SECTION("Testing Colorizer inputs")
{
    #ifdef OPENMS_WINDOWSPLATTFORM

        stringstream colored_stream;
        string comparison_string;

        //char////////////////////////////////////
        colored_stream << black(test_char);
        comparison_string.append(convertToString(test_char));
        TEST_EQUAL(colored_stream.str(), comparison_string)
        
        //clearing streams
        colored_stream.str(string());
        colored_stream.clear();
        comparison_string = "";

        //int/////////////////////////////////////
        colored_stream << cyan(test_int);
        comparison_string.append(convertToString(test_int));
        TEST_EQUAL(colored_stream.str(), comparison_string)
        
        //clearing streams
        colored_stream.str(string());
        colored_stream.clear();
        comparison_string = "";

        //floag///////////////////////////////////
        colored_stream << magenta(test_float);
        comparison_string.append(convertToString(test_float));
        TEST_EQUAL(colored_stream.str(), comparison_string)
        
        //clearing streams
        colored_stream.str(string());
        colored_stream.clear();
        comparison_string = "";

    #elif defined(__linux__) || defined(__OSX__)
        stringstream colored_stream;
        string comparison_string;

        //char////////////////////////////////////
        colored_stream << yellow(test_char);
        comparison_string.append(yellowANSI);
        comparison_string.append(convertToString(test_char));
        comparison_string.append(resetColorANSI);
        TEST_EQUAL(colored_stream.str(), comparison_string)

        //clearing streams
        colored_stream.str(string());
        colored_stream.clear();
        comparison_string = "";

        //int///////////////////////////////
        colored_stream << green(test_int);
        comparison_string.append(greenANSI);
        comparison_string.append(convertToString(test_int));
        comparison_string.append(resetColorANSI);
        TEST_EQUAL(colored_stream.str(), comparison_string)

        //clearing streams
        colored_stream.str(string());
        colored_stream.clear();
        comparison_string = "";

        //float///////////////////////////////
        colored_stream << red(test_float);
        comparison_string.append(redANSI);
        comparison_string.append(convertToString(test_float));
        comparison_string.append(resetColorANSI);
        TEST_EQUAL(colored_stream.str(), comparison_string)

        //clearing streams
        colored_stream.str(string());
        colored_stream.clear();
        comparison_string = "";
    #endif

}
END_SECTION

START_SECTION(Colorizer& operator()())
{
    // stringstream test_stream;
    // test_stream << blue();
    // test_stream << "ABDCE";
    // TEST_EQUAL(test_stream.str(),blueANSI+"ABCDE")

    // test_stream << reset_color();
    // TEST_EQUAL(test_stream.str(),blueANSI+"ABCDE"+resetColorANSI)
    
}
END_SECTION

END_TEST

