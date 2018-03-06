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
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Chris Bielow, Sandro Andreotti $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/KERNEL/RangeUtils.h>

///////////////////////////
#include <OpenMS/SIMULATION/RawTandemMSSignalSimulation.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(RawTandemMSSignalSimulation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

RawTandemMSSignalSimulation* ptr = nullptr;
RawTandemMSSignalSimulation* null_ptr = nullptr;
SimTypes::MutableSimRandomNumberGeneratorPtr rng (new SimTypes::SimRandomNumberGenerator);
rng->initialize(false, false);

START_SECTION((RawTandemMSSignalSimulation(SimTypes::SimRandomNumberGeneratorPtr rng)))
{
	ptr = new RawTandemMSSignalSimulation(rng);
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~RawTandemMSSignalSimulation())
{
	delete ptr;
}
END_SECTION

START_SECTION((RawTandemMSSignalSimulation(const RawTandemMSSignalSimulation &source)))
{
    ptr = new RawTandemMSSignalSimulation(rng);
    Param tmp_par = ptr->getParameters();
    tmp_par.setValue("status", "precursor");
    ptr->setParameters(tmp_par);

    RawTandemMSSignalSimulation copy(*ptr);
    TEST_EQUAL(copy.getParameters(), ptr->getParameters())
}
END_SECTION

START_SECTION((RawTandemMSSignalSimulation& operator=(const RawTandemMSSignalSimulation &source)))
{
    RawTandemMSSignalSimulation copy(rng);
    copy = *ptr;
    TEST_EQUAL(copy.getParameters(), ptr->getParameters())
}
END_SECTION

START_SECTION((void generateRawTandemSignals(const SimTypes::FeatureMapSim &, SimTypes::MSSimExperiment &, SimTypes::MSSimExperiment &)))
{
    rng->initialize(false, false);

    //Load featureXML and MSExperiment from MSSimulator run without MS2 simulation
    String feature_filename, exp_no_ms2_file, exp_with_ms2_file;
    feature_filename = OPENMS_GET_TEST_DATA_PATH("RawTandemMSSignalSimulation_no_ms2.featureXML");
    exp_no_ms2_file = OPENMS_GET_TEST_DATA_PATH("RawTandemMSSignalSimulation_no_ms2.mzML");
    exp_with_ms2_file = OPENMS_GET_TEST_DATA_PATH("RawTandemMSSignalSimulation_with_ms2.mzML");;
    SimTypes::FeatureMapSim features;
    SimTypes::MSSimExperiment exp_no_ms2, exp_with_ms2, peak_map;
    FeatureXMLFile().load(feature_filename, features);
    MzMLFile().load(exp_no_ms2_file, exp_no_ms2);
    MzMLFile().load(exp_with_ms2_file, exp_with_ms2);

    RawTandemMSSignalSimulation sim(rng);
    Param p;
    p.setValue("status", "precursor");
    p.setValue("tandem_mode", 2);
    p.setValue("TandemSim:SVM:hide_losses", "true");
    p.setValue("Precursor:Exclusion:use_dynamic_exclusion", "true");
    p.setValue("Precursor:Exclusion:exclusion_time", 50.0);
    p.setValue("Precursor:min_mz_peak_distance", 3.0);
    sim.setParameters(p);

    sim.generateRawTandemSignals(features, exp_no_ms2, peak_map);
    IntList levels;
    levels.push_back(1);
    exp_no_ms2.getSpectra().erase(remove_if(exp_no_ms2.begin(), exp_no_ms2.end(), InMSLevelRange<SimTypes::MSSimExperiment::SpectrumType>(levels)), exp_no_ms2.end());
    exp_with_ms2.getSpectra().erase(remove_if(exp_with_ms2.begin(), exp_with_ms2.end(), InMSLevelRange<SimTypes::MSSimExperiment::SpectrumType>(levels)), exp_with_ms2.end());
//    MzMLFile().store(OPENMS_GET_TEST_DATA_PATH("RawTandemMSSignalSimulation_with_ms2.mzML"), exp_no_ms2);

    TEST_EQUAL(exp_with_ms2.size(), exp_no_ms2.size());
    TEST_EQUAL(exp_with_ms2[0].size(), exp_no_ms2[0].size());
#if OPENMS_BOOST_VERSION_MINOR < 56
    TEST_EQUAL(exp_with_ms2[1].size(), exp_no_ms2[1].size());
#else
    TEST_EQUAL(exp_with_ms2[1].size(), exp_no_ms2[1].size()-1);
#endif
    TEST_REAL_SIMILAR(exp_with_ms2[0].getPrecursors()[0].getMZ(), exp_no_ms2[0].getPrecursors()[0].getMZ());
    TEST_REAL_SIMILAR(exp_with_ms2[1].getPrecursors()[0].getMZ(), exp_no_ms2[1].getPrecursors()[0].getMZ());
    TEST_REAL_SIMILAR(exp_with_ms2[0][0].getIntensity(), exp_no_ms2[0][0].getIntensity());
    TEST_REAL_SIMILAR(exp_with_ms2[0][1].getIntensity(), exp_no_ms2[0][1].getIntensity());
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



