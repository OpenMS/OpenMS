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
// $Authors: Erhan Kenar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/FuzzyStringComparator.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>
#include <OpenMS/KERNEL/MSExperiment.h>

///////////////////////////
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFindingMetabo, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFindingMetabo* ptr = nullptr;
FeatureFindingMetabo* null_ptr = nullptr;
START_SECTION(FeatureFindingMetabo())
{
    ptr = new FeatureFindingMetabo();
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~FeatureFindingMetabo())
{
    delete ptr;
}
END_SECTION

// load a mzML file for testing the algorithm
PeakMap input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_input1.mzML"), input);

FeatureMap test_fm;

std::vector<MassTrace> output_mt, splitted_mt, filtered_mt;

std::vector<std::vector< OpenMS::MSChromatogram > > chromatograms;

MassTraceDetection test_mtd;
test_mtd.run(input, output_mt);

ElutionPeakDetection test_epd;
test_epd.detectPeaks(output_mt, splitted_mt);

FuzzyStringComparator fsc;
fsc.setAcceptableRelative(1.001);
fsc.setAcceptableAbsolute(1);
StringList sl;
sl.push_back("xml-stylesheet");
sl.push_back("<featureMap");
sl.push_back("<feature id");
fsc.setWhitelist(sl);

//std::cout << "\n\n" << fsc.compareStrings("529090", "529091") << "\n\n\n";

START_SECTION((void run(std::vector< MassTrace > &, FeatureMap &, chromatograms &)))
{
    FeatureFindingMetabo test_ffm;
    // run with non-default setting (C13 isotope distance)
    Param p = test_ffm.getParameters();
    p.setValue("mz_scoring_13C", "true");
    test_ffm.setParameters(p);
    test_ffm.run(splitted_mt, test_fm, chromatograms);
    TEST_EQUAL(test_fm.size(), 93);

    // run with default settings (from paper using charge+isotope# dependent distances)
    p.setValue("report_convex_hulls", "true");
    p.setValue("mz_scoring_13C", "false");
    test_ffm.setParameters(p);
    test_ffm.run(splitted_mt, test_fm, chromatograms);
    TEST_EQUAL(test_fm.size(), 91);
    // --> this gives less features, i.e. more isotope clusters (but the input data is simulated and highly weird -- should be replaced at some point)

    // test annotation of input
    String tmp_file;
    NEW_TMP_FILE(tmp_file);
    FeatureXMLFile().store(tmp_file, test_fm);
    TEST_EQUAL(fsc.compareFiles(tmp_file, OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_output1.featureXML")), true);

}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



