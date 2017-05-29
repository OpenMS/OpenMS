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
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/KERNEL/MSSpectrum.h>
#include <OpenMS/KERNEL/MSChromatogram.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(SpectrumHelper, "$Id$")

START_SECTION((MSSpectrum<>::FloatDataArrays::iterator getDataArrayByName(MSSpectrum::FloatDataArrays& a, const String& name)))
	MSSpectrum<> ds;
	MSSpectrum<>::FloatDataArray float_array;
        float_array.push_back(56);  
        float_array.push_back(201); 
        float_array.push_back(31);  
        float_array.push_back(201); 
        float_array.push_back(201);  
	
        ds.getFloatDataArrays() = std::vector<MSSpectrum<>::FloatDataArray>(3, float_array);
	ds.getFloatDataArrays()[0].setName("f1");
	ds.getFloatDataArrays()[1].setName("f2");
	ds.getFloatDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f2") == ds.getFloatDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "NOT_THERE") == ds.getFloatDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f1") - ds.getFloatDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f2") - ds.getFloatDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f3") - ds.getFloatDataArrays().begin(), 2);

        // test const version
        const MSSpectrum<>::FloatDataArrays cfa(ds.getFloatDataArrays());
        TEST_EQUAL(getDataArrayByName(cfa, "f2") == cfa.end(), false);
        TEST_EQUAL(getDataArrayByName(cfa, "f1") - cfa.begin(), 0);
        TEST_EQUAL(getDataArrayByName(cfa, "NOT_THERE") == cfa.end(), true);

END_SECTION

START_SECTION((MSSpectrum<>::StringDataArrays::iterator getDataArrayByName(MSSpectrum::StringDataArrays& a, const String& name)))
	MSSpectrum<> ds;
	MSSpectrum<>::StringDataArray string_array;
        string_array.push_back("56");  
        string_array.push_back("201"); 
        string_array.push_back("31");  
        string_array.push_back("201"); 
        string_array.push_back("201");  
	ds.getStringDataArrays() = std::vector<MSSpectrum<>::StringDataArray>(3, string_array);
	ds.getStringDataArrays()[0].setName("f1");
	ds.getStringDataArrays()[1].setName("f2");
	ds.getStringDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f2") == ds.getStringDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "NOT_THERE") == ds.getStringDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f1") - ds.getStringDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f2") - ds.getStringDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f3") - ds.getStringDataArrays().begin(), 2);
END_SECTION

START_SECTION((MSSpectrum<>::IntegerDataArrays::iterator getDataArrayByName(MSSpectrum::IntegerDataArrays& a, const String& name)))
	MSSpectrum<> ds;
	MSSpectrum<>::IntegerDataArray int_array;
        int_array.push_back(56);  
        int_array.push_back(201); 
        int_array.push_back(31);  
        int_array.push_back(201); 
        int_array.push_back(201);  

	ds.getIntegerDataArrays() = std::vector<MSSpectrum<>::IntegerDataArray>(3, int_array);
	ds.getIntegerDataArrays()[0].setName("f1");
	ds.getIntegerDataArrays()[1].setName("f2");
	ds.getIntegerDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f2") == ds.getIntegerDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "NOT_THERE") == ds.getIntegerDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f1") - ds.getIntegerDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f2") - ds.getIntegerDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f3") - ds.getIntegerDataArrays().begin(), 2);
END_SECTION

START_SECTION((MSChromatogram<>::FloatDataArrays::iterator getDataArrayByName(const MSChromatogram::FloatDataArrays& a, const String& name)))
	MSChromatogram<> ds;
	MSChromatogram<>::FloatDataArray float_array;
        float_array.push_back(56);  
        float_array.push_back(201); 
        float_array.push_back(31);  
        float_array.push_back(201); 
        float_array.push_back(201);  
	
        ds.getFloatDataArrays() = std::vector<MSChromatogram<>::FloatDataArray>(3, float_array);
	ds.getFloatDataArrays()[0].setName("f1");
	ds.getFloatDataArrays()[1].setName("f2");
	ds.getFloatDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f2") == ds.getFloatDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "NOT_THERE") == ds.getFloatDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f1") - ds.getFloatDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f2") - ds.getFloatDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getFloatDataArrays(), "f3") - ds.getFloatDataArrays().begin(), 2);
END_SECTION

START_SECTION((MSChromatogram<>::StringDataArrays::iterator getDataArrayByName(MSChromatogram::StringDataArrays& a, const String& name)))
	MSChromatogram<> ds;
	MSChromatogram<>::StringDataArray string_array;
        string_array.push_back("56");  
        string_array.push_back("201"); 
        string_array.push_back("31");  
        string_array.push_back("201"); 
        string_array.push_back("201");  
	ds.getStringDataArrays() = std::vector<MSChromatogram<>::StringDataArray>(3, string_array);
	ds.getStringDataArrays()[0].setName("f1");
	ds.getStringDataArrays()[1].setName("f2");
	ds.getStringDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f2") == ds.getStringDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "NOT_THERE") == ds.getStringDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f1") - ds.getStringDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f2") - ds.getStringDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getStringDataArrays(), "f3") - ds.getStringDataArrays().begin(), 2);
END_SECTION

START_SECTION((MSChromatogram<>::IntegerDataArrays::iterator getDataArrayByName(MSChromatogram::IntegerDataArrays& a, const String& name)))
	MSChromatogram<> ds;
	MSChromatogram<>::IntegerDataArray int_array;
        int_array.push_back(56);  
        int_array.push_back(201); 
        int_array.push_back(31);  
        int_array.push_back(201); 
        int_array.push_back(201);  

	ds.getIntegerDataArrays() = std::vector<MSChromatogram<>::IntegerDataArray>(3, int_array);
	ds.getIntegerDataArrays()[0].setName("f1");
	ds.getIntegerDataArrays()[1].setName("f2");
	ds.getIntegerDataArrays()[2].setName("f3");

        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f2") == ds.getIntegerDataArrays().end(), false);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "NOT_THERE") == ds.getIntegerDataArrays().end(), true);

        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f1") - ds.getIntegerDataArrays().begin(), 0);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f2") - ds.getIntegerDataArrays().begin(), 1);
        TEST_EQUAL(getDataArrayByName(ds.getIntegerDataArrays(), "f3") - ds.getIntegerDataArrays().begin(), 2);
END_SECTION


END_TEST

