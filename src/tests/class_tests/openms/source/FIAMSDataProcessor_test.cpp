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
    ptr = new FIAMSDataProcessor();
    TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(virtual ~FIAMSDataProcessor())
{
    delete ptr;
}
END_SECTION

// START_SECTION([EXTRA]AdductInfo)
// {
//   EmpiricalFormula ef_empty;
//   // make sure an empty formula has no weight (we rely on that in AdductInfo's getMZ() and getNeutralMass()
//   TEST_EQUAL(ef_empty.getMonoWeight(), 0)

//   // now we test if converting from neutral mass to m/z and back recovers the input value using different adducts
//   {
//   // testing M;-2  // intrinsic doubly negative charge
//     AdductInfo ai("TEST_INTRINSIC", ef_empty, -2, 1);
//     double neutral_mass=1000; // some mass...
//     double mz = ai.getMZ(neutral_mass);
//     double neutral_mass_recon = ai.getNeutralMass(mz);
//     TEST_REAL_SIMILAR(neutral_mass, neutral_mass_recon);
//   }
//   { // testing M+Na+H;+2
//     EmpiricalFormula simpleAdduct("HNa");
//     AdductInfo ai("TEST_WITHADDUCT", simpleAdduct, 2, 1);
//     double neutral_mass=1000; // some mass...
//     double mz = ai.getMZ(neutral_mass);
//     double neutral_mass_recon = ai.getNeutralMass(mz);
//     TEST_REAL_SIMILAR(neutral_mass, neutral_mass_recon);
//   }

// }
// END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST