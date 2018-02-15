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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/FORMAT/MzDataFile.h>

#include <cmath>
#include <fstream>

using namespace OpenMS;
using namespace std;

START_TEST(IsotopeWaveletTransform, "$Id$")

PeakMap map;
MzDataFile file; file.load (OPENMS_GET_TEST_DATA_PATH("IsotopeWaveletTransform_test.mzData"), map);
map.updateRanges();
IsotopeWaveletTransform<Peak1D>* iw = nullptr;
IsotopeWaveletTransform<Peak1D>* nullIw = nullptr;
IsotopeWaveletTransform<Peak1D>::TransSpectrum* test2 = nullptr;
MSSpectrum* spec = new MSSpectrum (map[0]);

START_SECTION([IsotopeWaveletTransform::TransSpectrum] TransSpectrum())
	IsotopeWaveletTransform<Peak1D>::TransSpectrum test;
	NOT_TESTABLE
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] TransSpectrum(const MSSpectrum* reference))
	test2 = new IsotopeWaveletTransform<Peak1D>::TransSpectrum (spec);
	const MSSpectrum* ref = test2->getRefSpectrum();
  const MSSpectrum* nullPtr = nullptr;
	TEST_NOT_EQUAL (ref, nullPtr)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] double getRT () const )
	TEST_EQUAL (test2->getRT(), 100.00)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] double getMZ (const UInt i) const)
	TEST_EQUAL ((int)(test2->getMZ(0)*10), 14200)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] double getRefIntensity (const UInt i) const)
	TEST_EQUAL ((int)(test2->getRefIntensity(0)*100), 39)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] double getTransIntensity (const UInt i) const)
	TEST_EQUAL (test2->getTransIntensity(0), 0)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] void setTransIntensity (const UInt i, const double intens))
	test2->setTransIntensity(0,-1);
	TEST_EQUAL (test2->getTransIntensity(0), -1)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] Size size() const)
	TEST_EQUAL (test2->size(), spec->size())
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] const MSSpectrum* getRefSpectrum ())
	const MSSpectrum* ref = test2->getRefSpectrum();
	TEST_EQUAL (ref, spec)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] const MSSpectrum* getRefSpectrum () const)
	const IsotopeWaveletTransform<Peak1D>::TransSpectrum* test3 = new IsotopeWaveletTransform<Peak1D>::TransSpectrum (spec);
	const MSSpectrum* ref = test3->getRefSpectrum();
	TEST_EQUAL (ref, spec)
	delete (test3);
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] MSSpectrum::const_iterator MZBegin (const double mz) const)
	TEST_EQUAL((int)(test2->MZBegin(1420)->getMZ()*10), 14200)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] MSSpectrum::const_iterator MZEnd (const double mz) const)
	TEST_EQUAL((int)(test2->MZEnd(1420.01)->getMZ()*100), 142001)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] MSSpectrum::const_iterator begin () const)
	TEST_EQUAL((int)(test2->begin()->getMZ()*10), 14200)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] MSSpectrum::const_iterator end () const)
	TEST_EQUAL((int)((--test2->end())->getMZ()*10), 14349)
END_SECTION

START_SECTION((IsotopeWaveletTransform(const double min_mz, const double max_mz, const UInt max_charge, const Size max_scan_size=0, const bool hr_data=false, const String intenstype="ref")))
	iw = new IsotopeWaveletTransform<Peak1D> (map[0].begin()->getMZ(), (map[0].end()-1)->getMZ(), 1);
	TEST_NOT_EQUAL (iw, nullIw)
END_SECTION

START_SECTION(void initializeScan(const MSSpectrum &c_ref, const UInt c=0))
	iw->initializeScan (map[0]);
	NOT_TESTABLE
END_SECTION

START_SECTION(Size getMaxScanSize () const)
	TEST_EQUAL (iw->getMaxScanSize(), 0);
	NOT_TESTABLE
END_SECTION

START_SECTION(void computeMinSpacing (const MSSpectrum& c_ref))
	iw->computeMinSpacing (map[0]);
	NOT_TESTABLE
END_SECTION

START_SECTION(double getMinSpacing () const)
	TEST_EQUAL ((int)(iw->getMinSpacing()*100), 1);
END_SECTION

START_SECTION(void getTransformHighRes(MSSpectrum &c_trans, const MSSpectrum &c_ref, const UInt c))
	iw->getTransformHighRes (*spec, map[0], 0);
	TEST_EQUAL (*spec!= map[0], true)
END_SECTION

START_SECTION(void getTransform(MSSpectrum &c_trans, const MSSpectrum &c_ref, const UInt c))
	iw->getTransform (*spec, map[0], 0);
	TEST_EQUAL (*spec!= map[0], true)
END_SECTION

START_SECTION(void setSigma (const double sigma))
	iw->setSigma (1);
	NOT_TESTABLE
END_SECTION

START_SECTION(double getSigma () const)
	TEST_EQUAL (iw->getSigma(), 1)
END_SECTION

START_SECTION(void identifyCharge(const MSSpectrum &candidates, const MSSpectrum &ref, const UInt scan_index, const UInt c, const double ampl_cutoff, const bool check_PPMs))
	iw->identifyCharge (*spec, map[0], 0, 0, 0, false);
	NOT_TESTABLE
END_SECTION

START_SECTION(void updateBoxStates(const MSExperiment< PeakType > &map, const Size scan_index, const UInt RT_interleave, const UInt RT_votes_cutoff, const Int front_bound=-1, const Int end_bound=-1))
	iw->updateBoxStates(map, INT_MAX, 0, 0);
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual std::multimap<double, Box> getClosedBoxes ()))
	TEST_EQUAL (iw->getClosedBoxes().size(), 1)
END_SECTION

START_SECTION(FeatureMap< Feature > mapSeeds2Features(const MSExperiment< PeakType > &map, const UInt RT_votes_cutoff))
	FeatureMap f = iw->mapSeeds2Features(map, 0);
	TEST_EQUAL (f.size(), 1)
END_SECTION

START_SECTION(void mergeFeatures(IsotopeWaveletTransform< PeakType > *later_iwt, const UInt RT_interleave, const UInt RT_votes_cutoff))
	NOT_TESTABLE //only via CUDA
END_SECTION

START_SECTION(double getLinearInterpolation(const typename MSSpectrum::const_iterator &left_iter, const double mz_pos, const typename MSSpectrum::const_iterator &right_iter))
	TEST_EQUAL((int)(iw->getLinearInterpolation(map[0].begin(), 1420.02, (map[0].begin()+1))*10),5)
END_SECTION

START_SECTION(double getLinearInterpolation(const double mz_a, const double intens_a, const double mz_pos, const double mz_b, const double intens_b))
	TEST_EQUAL(iw->getLinearInterpolation(1,1, 1.5, 2, 2), 1.5)
END_SECTION


START_SECTION(~IsotopeWaveletTransform())
	delete (iw);
END_SECTION


START_SECTION([IsotopeWaveletTransform::TransSpectrum] void destroy ())
	test2->destroy();
	NOT_TESTABLE
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] virtual ~TransSpectrum())
	delete(test2);
END_SECTION

END_TEST
