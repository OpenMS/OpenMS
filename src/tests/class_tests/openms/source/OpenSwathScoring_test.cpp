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
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------
#include <boost/shared_ptr.hpp>

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathScoring.h>
///////////////////////////

// we dont want any inclusion of OpenMS Kernel classes here ...
#ifdef OPENMS_KERNEL_MSSPECTRUM_H
ThisShouldFailAtCompileTime = 0
#endif
#ifdef OPENMS_KERNEL_MRMFEATURE_H
ThisShouldFailAtCompileTime = 0
#endif

#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/DataAccessHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/DATAACCESS/SimpleOpenMSSpectraAccessFactory.h>

using namespace OpenMS;
using namespace std;

START_TEST(OpenSwathScoring, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

OpenSwathScoring* ptr = nullptr;
OpenSwathScoring* nullPointer = nullptr;

START_SECTION(OpenSwathScoring())
{
	ptr = new OpenSwathScoring();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(~OpenSwathScoring())
{
  delete ptr;
}
END_SECTION

START_SECTION((void initialize(double rt_normalization_factor_, int add_up_spectra_, double spacing_for_spectra_resampling_, OpenSwath_Scores_Usage & su_)))
{
	ptr = new OpenSwathScoring();
  OpenSwath_Scores_Usage su;
	TEST_NOT_EQUAL(ptr, nullPointer)
  ptr->initialize(100.0, 1, 0.01, su);
  delete ptr;
}
END_SECTION

START_SECTION((void calculateChromatographicScores( OpenSwath::IMRMFeature* imrmfeature, const std::vector<std::string>& native_ids, const std::vector<double>& normalized_library_intensity, std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators, OpenSwath_Scores & scores) ))
{
  NOT_TESTABLE // see MRMFeatureFinderScoring_test.cpp
  // - the OpenSwathScoring is a facade and thus does not need testing on its own
}
END_SECTION

START_SECTION((void calculateChromatographicIdScores( OpenSwath::IMRMFeature* imrmfeature, const std::vector<std::string>& native_ids_identification,, const std::vector<std::string>& native_ids_detection, std::vector<OpenSwath::ISignalToNoisePtr>& signal_noise_estimators, OpenSwath_Scores & idscores) ))
{
  NOT_TESTABLE // see MRMFeatureFinderScoring_test.cpp
  // - the OpenSwathScoring is a facade and thus does not need testing on its own
}
END_SECTION

START_SECTION((void calculateLibraryScores( OpenSwath::IMRMFeature* imrmfeature, const std::vector<TransitionType> & transitions, const PeptideType& pep, const double normalized_feature_rt, OpenSwath_Scores & scores)))
{
  NOT_TESTABLE // see MRMFeatureFinderScoring_test.cpp
  // - the OpenSwathScoring is a facade and thus does not need testing on its own
}
END_SECTION

START_SECTION((void calculateDIAScores(OpenSwath::IMRMFeature* imrmfeature, const std::vector<TransitionType> & transitions, OpenSwath::SpectrumAccessPtr swath_map, OpenMS::DIAScoring & diascoring, const PeptideType& pep, OpenSwath_Scores & scores)))
{
  NOT_TESTABLE // see MRMFeatureFinderScoring_test.cpp
  // - the OpenSwathScoring is a facade and thus does not need testing on its own
}
END_SECTION

START_SECTION((void getNormalized_library_intensities_(const std::vector<TransitionType> & transitions, std::vector<double>& normalized_library_intensity)))
{
  NOT_TESTABLE // see MRMFeatureFinderScoring_test.cpp
  // - the OpenSwathScoring is a facade and thus does not need testing on its own
}
END_SECTION

START_SECTION((OpenSwath::SpectrumPtr getAddedSpectra_(OpenSwath::SpectrumAccessPtr swath_map, double RT, int nr_spectra_to_add)))
{

  // test result for empty map
  {
    boost::shared_ptr<PeakMap > swath_map (new PeakMap);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    OpenSwathScoring sc;
    OpenSwath::SpectrumPtr sp = sc.getAddedSpectra_(swath_ptr, 20.0, 1);

    TEST_EQUAL(sp->getMZArray()->data.empty(), true);
  }

  // test result for map with single spectrum
  {
    PeakMap* eptr = new PeakMap;
    MSSpectrum s;
    Peak1D p;
    p.setMZ(20.0);
    p.setIntensity(200.0);
    s.push_back(p);
    s.setRT(20.0);
    eptr->addSpectrum(s);
    boost::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    TEST_EQUAL(swath_ptr->getNrSpectra(), 1)
    OpenSwathScoring sc;
    OpenSwath::SpectrumPtr sp = sc.getAddedSpectra_(swath_ptr, 20.0, 1);

    TEST_EQUAL(sp->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp->getIntensityArray()->data.size(), 1);

    TEST_REAL_SIMILAR(sp->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp->getIntensityArray()->data[0], 200.0);
  }

  // test result for map with three spectra
  {
    PeakMap* eptr = new PeakMap;
    MSSpectrum s;
    Peak1D p;
    p.setMZ(20.0);
    p.setIntensity(200.0);
    s.push_back(p);
    s.setRT(10.0);
    eptr->addSpectrum(s);
    s.setRT(20.0);
    eptr->addSpectrum(s);
    s.setRT(30.0);
    eptr->addSpectrum(s);
    boost::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    TEST_EQUAL(swath_ptr->getNrSpectra(), 3)
    OpenSwathScoring sc;
    OpenSwath::SpectrumPtr sp = sc.getAddedSpectra_(swath_ptr, 20.0, 3);

    TEST_EQUAL(sp->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp->getIntensityArray()->data.size(), 1);

    TEST_REAL_SIMILAR(sp->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp->getIntensityArray()->data[0], 600.0);
  }

  // test result for map with uneven number of spectra
  {
    PeakMap* eptr = new PeakMap;
    MSSpectrum s;
    Peak1D p;
    p.setMZ(20.0);
    p.setIntensity(200.0);
    s.push_back(p);
    s.setRT(10.0);
    eptr->addSpectrum(s);
    s.setRT(20.0);
    eptr->addSpectrum(s);
    /*
    s.setRT(30.0);
    eptr->addSpectrum(s);
    */
    boost::shared_ptr<PeakMap > swath_map (eptr);
    OpenSwath::SpectrumAccessPtr swath_ptr = SimpleOpenMSSpectraFactory::getSpectrumAccessOpenMSPtr(swath_map);

    TEST_EQUAL(swath_ptr->getNrSpectra(), 2)
    OpenSwathScoring sc;
    OpenSwath::SpectrumPtr sp = sc.getAddedSpectra_(swath_ptr, 20.0, 3);

    TEST_EQUAL(sp->getMZArray()->data.size(), 1);
    TEST_EQUAL(sp->getIntensityArray()->data.size(), 1);

    TEST_REAL_SIMILAR(sp->getMZArray()->data[0], 20.0);
    TEST_REAL_SIMILAR(sp->getIntensityArray()->data[0], 400.0);
  }
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

