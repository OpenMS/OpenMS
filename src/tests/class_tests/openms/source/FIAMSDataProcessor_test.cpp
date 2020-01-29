// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Authors: Erhan Kenar, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/ID/FIAMSDataProcessor.h>

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/PeakWidthEstimator.h>
#include <OpenMS/TRANSFORMATIONS/RAW2PEAK/PeakPickerCWT.h>
#include <OpenMS/ANALYSIS/OPENSWATH/TargetedSpectraExtractor.h>
#include <OpenMS/ANALYSIS/ID/AccurateMassSearchEngine.h>
#include <OpenMS/FORMAT/MzTabFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/FILTERING/NOISEESTIMATION/SignalToNoiseEstimatorMedianRapid.h>
#include <OpenMS/FILTERING/BASELINE/MorphologicalFilter.h>
#include <OpenMS/KERNEL/SpectrumHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/SpectrumAddition.h>
#include <OpenMS/FILTERING/SMOOTHING/SavitzkyGolayFilter.h>


///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(FIAMSDataProcessor, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FIAMSDataProcessor* ptr = nullptr;
FIAMSDataProcessor* null_ptr = nullptr;
START_SECTION(FIAMSDataProcessor())
{
    ptr = new FIAMSDataProcessor(120000.0);
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~FIAMSDataProcessor())
{
    delete ptr;
}
END_SECTION

FIAMSDataProcessor fia_processor(120000.0);
PeakMap input;
Peak1D p;
std::vector<float> ints {100, 120, 130, 140, 150, 100, 60, 50, 30};
std::vector<float> rts {10, 20, 30, 40};
for (Size i = 0; i < rts.size(); ++i) {
    MSSpectrum s;
    for (Size j = 0; j < ints.size(); ++j) {
        p.setIntensity(ints[j]); p.setMZ(100 + j*2);
        s.push_back(p);
    }
    s.setRT(rts[i]);
    input.addSpectrum(s);
}

START_SECTION((void cutForTime(const MSExperiment & experiment, vector<MSSpectrum> & output, float n_seconds)))
{
    vector<MSSpectrum> output1;
    fia_processor.cutForTime(input, output1, 0);
    TEST_EQUAL(output1.size(), 0);
    vector<MSSpectrum> output2;
    fia_processor.cutForTime(input, output2, 25);
    TEST_EQUAL(output2.size(), 2);
    vector<MSSpectrum> output3;
    fia_processor.cutForTime(input, output3, 100);
    TEST_EQUAL(output3.size(), 4);
    PeakMap empty_input;
    vector<MSSpectrum> output4;
    fia_processor.cutForTime(empty_input, output4, 100);
    TEST_EQUAL(output4.size(), 0);

}
END_SECTION

START_SECTION((mergeAlongTime(
  const std::vector<OpenMS::MSSpectrum> & input, 
  )))
{
    vector<MSSpectrum> output_cut;
    fia_processor.cutForTime(input, output_cut, 100);
    MSSpectrum output = fia_processor.mergeAlongTime(output_cut);
    TEST_EQUAL(output.size(), 0);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST