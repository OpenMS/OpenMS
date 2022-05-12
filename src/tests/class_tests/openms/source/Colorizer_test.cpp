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


START_TEST(Colorizer(),"$Id$")

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
START_SECTION(Colorizer black)
{
    //Test variables
    char tchar = 'a';
    // unsigned char unsignedChar = 'i';
    // signed char signedChar = 's';
    // short int shortInt = 32766;
    // unsigned short int unsShortInt = 600;
    // signed short int signShortInt = -32560;
    // unsigned long int unsLongInt = 40000000000;
    // long long int longLongInt = 981278728478274;
    // unsigned long long int unsLongLongInt = -78273829375;
    // float flt = 2094.5892;
    // double dbl = -253575634.98925;
    // long double longDbl= 135315.2929849375;
    // wchar_t wideChar = L'S';

    string outStream;
    string outString;
    stringstream blackStream;

    blackStream << black(tchar);
    blackStream >> outStream;
    outString = "\e[30m" + tchar;
    outString.append("\e[0m");
    TEST_EQUAL(outStream, outString)


}
END_SECTION

#endif



//Testing other variables
/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

