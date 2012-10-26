// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Erhan Kenar$
// $Authors: Erhan Kenar$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FILTERING/DATAREDUCTION/MassTraceDetection.h>
#include <OpenMS/FILTERING/DATAREDUCTION/ElutionPeakDetection.h>

///////////////////////////
#include <OpenMS/FILTERING/DATAREDUCTION/FeatureFindingMetabo.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FeatureFindingMetabo, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FeatureFindingMetabo* ptr = 0;
FeatureFindingMetabo* null_ptr = 0;
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
MSExperiment<Peak1D> input;
MzMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_input1.mzML"), input);

FeatureMap<> exp_fm, test_fm;
FeatureXMLFile().load(OPENMS_GET_TEST_DATA_PATH("FeatureFindingMetabo_output1.featureXML"), exp_fm);
// exp_fm.sortByMZ();

std::vector<MassTrace> output_mt, splitted_mt, filtered_mt;

MassTraceDetection test_mtd;
test_mtd.run(input, output_mt);

ElutionPeakDetection test_epd;
test_epd.detectPeaks(output_mt, splitted_mt);
// test_epd.filterByPeakWidth(splitted_mt, filtered_mt);

std::cout << "!!!!" << splitted_mt.size() << std::endl;

START_SECTION((void run(std::vector< MassTrace > &, FeatureMap<> &)))
{
    FeatureFindingMetabo test_ffm;
    test_ffm.run(splitted_mt, test_fm);
    // test_fm.sortByMZ();

    TEST_EQUAL(exp_fm.size(), test_fm.size());

    for (Size i = 0; i < exp_fm.size(); ++i)
    {
        TEST_EQUAL(exp_fm[i].getMetaValue(3), test_fm[i].getMetaValue(3));
        TEST_REAL_SIMILAR(exp_fm[i].getRT(), test_fm[i].getRT());
        TEST_REAL_SIMILAR(exp_fm[i].getMZ(), test_fm[i].getMZ());
        TEST_REAL_SIMILAR(exp_fm[i].getIntensity(), test_fm[i].getIntensity());
    }
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



