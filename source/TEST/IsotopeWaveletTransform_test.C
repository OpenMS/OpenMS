// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Rene Hussong $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWaveletTransform.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeWavelet.h>
#include <OpenMS/KERNEL/Peak1D.h>
#include <OpenMS/FORMAT/MzDataFile.h>

#include <math.h>
#include <fstream>

using namespace OpenMS;
using namespace std;

START_TEST(IsotopeWaveletTransform, "$Id$")

MSExperiment<> map;
MzDataFile file; file.load (OPENMS_GET_TEST_DATA_PATH("IsotopeWaveletTransform_test.mzData"), map);
map.updateRanges();
IsotopeWaveletTransform<Peak1D>* iw = 0;
IsotopeWaveletTransform<Peak1D>::TransSpectrum* test2 = 0;
MSSpectrum<Peak1D>* spec = new MSSpectrum<Peak1D> (map[0]);

START_SECTION([IsotopeWaveletTransform::TransSpectrum] TransSpectrum())
	IsotopeWaveletTransform<Peak1D>::TransSpectrum test;
	NOT_TESTABLE
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] TransSpectrum(const MSSpectrum<PeakType>* reference))
	test2 = new IsotopeWaveletTransform<Peak1D>::TransSpectrum (spec);
	const MSSpectrum<Peak1D>* ref = test2->getRefSpectrum();
	TEST_NOT_EQUAL (ref, 0)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] DoubleReal getRT () const )
	TEST_EQUAL (test2->getRT(), 100.00)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] DoubleReal getMZ (const UInt i) const)
	TEST_EQUAL ((int)(test2->getMZ(0)*10), 14200)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] DoubleReal getRefIntensity (const UInt i) const)
	TEST_EQUAL ((int)(test2->getRefIntensity(0)*100), 39)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] DoubleReal getTransIntensity (const UInt i) const)
	TEST_EQUAL (test2->getTransIntensity(0), 0)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] void setTransIntensity (const UInt i, const DoubleReal intens))
	test2->setTransIntensity(0,-1);
	TEST_EQUAL (test2->getTransIntensity(0), -1)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] Size size() const)
	TEST_EQUAL (test2->size(), spec->size())
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] const MSSpectrum<PeakType>* getRefSpectrum ())
	const MSSpectrum<Peak1D>* ref = test2->getRefSpectrum();
	TEST_EQUAL (ref, spec)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] const MSSpectrum<PeakType>* getRefSpectrum () const)
	const IsotopeWaveletTransform<Peak1D>::TransSpectrum* test3 = new IsotopeWaveletTransform<Peak1D>::TransSpectrum (spec);
	const MSSpectrum<Peak1D>* ref = test3->getRefSpectrum();
	TEST_EQUAL (ref, spec)
	delete (test3);
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] MSSpectrum<PeakType>::const_iterator MZBegin (const DoubleReal mz) const)
	TEST_EQUAL((int)(test2->MZBegin(1420)->getMZ()*10), 14200)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] MSSpectrum<PeakType>::const_iterator MZEnd (const DoubleReal mz) const)
	TEST_EQUAL((int)(test2->MZEnd(1420.01)->getMZ()*100), 142001)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] MSSpectrum<PeakType>::const_iterator begin () const)
	TEST_EQUAL((int)(test2->begin()->getMZ()*10), 14200)
END_SECTION

START_SECTION([IsotopeWaveletTransform::TransSpectrum] MSSpectrum<PeakType>::const_iterator end () const)
	TEST_EQUAL((int)((--test2->end())->getMZ()*10), 14349)
END_SECTION

START_SECTION((IsotopeWaveletTransform(const DoubleReal min_mz, const DoubleReal max_mz, const UInt max_charge, const Size max_scan_size=0, const bool use_cuda=false, const bool hr_data=false, const String intenstype="ref")))
	iw = new IsotopeWaveletTransform<Peak1D> (map[0].begin()->getMZ(), (map[0].end()-1)->getMZ(), 1);
	TEST_NOT_EQUAL (iw, 0)
END_SECTION

START_SECTION(void initializeScan(const MSSpectrum< PeakType > &c_ref, const UInt c=0))
	iw->initializeScan (map[0]);
	NOT_TESTABLE
END_SECTION

START_SECTION(Size getMaxScanSize () const)
	TEST_EQUAL (iw->getMaxScanSize(), 0);
	NOT_TESTABLE
END_SECTION

START_SECTION(void computeMinSpacing (const MSSpectrum<PeakType>& c_ref))
	iw->computeMinSpacing (map[0]);
	NOT_TESTABLE
END_SECTION

START_SECTION(DoubleReal getMinSpacing () const)
	TEST_EQUAL ((int)(iw->getMinSpacing()*100), 1);
END_SECTION

START_SECTION(void getTransformHighRes(MSSpectrum<PeakType> &c_trans, const MSSpectrum<PeakType> &c_ref, const UInt c))
	iw->getTransformHighRes (*spec, map[0], 0);
	TEST_EQUAL (*spec!= map[0], true)
END_SECTION

START_SECTION(void getTransform(MSSpectrum<PeakType> &c_trans, const MSSpectrum<PeakType> &c_ref, const UInt c))
	iw->getTransform (*spec, map[0], 0);
	TEST_EQUAL (*spec!= map[0], true)
END_SECTION

START_SECTION(void setSigma (const DoubleReal sigma))
	iw->setSigma (1);
	NOT_TESTABLE
END_SECTION

START_SECTION(DoubleReal getSigma () const)
	TEST_EQUAL (iw->getSigma(), 1)
END_SECTION

START_SECTION(void identifyCharge(const MSSpectrum< PeakType > &candidates, const MSSpectrum< PeakType > &ref, const UInt scan_index, const UInt c, const DoubleReal ampl_cutoff, const bool check_PPMs))
	iw->identifyCharge (*spec, map[0], 0, 0, 0, false);
	NOT_TESTABLE
END_SECTION

START_SECTION(void updateBoxStates(const MSExperiment< PeakType > &map, const Size scan_index, const UInt RT_interleave, const UInt RT_votes_cutoff, const Int front_bound=-1, const Int end_bound=-1))
	iw->updateBoxStates(map, INT_MAX, 0, 0);
	NOT_TESTABLE
END_SECTION

START_SECTION((virtual std::multimap<DoubleReal, Box> getClosedBoxes ()))
	TEST_EQUAL (iw->getClosedBoxes().size(), 1)
END_SECTION

START_SECTION(FeatureMap< Feature > mapSeeds2Features(const MSExperiment< PeakType > &map, const UInt RT_votes_cutoff))
	FeatureMap<Feature> f = iw->mapSeeds2Features(map, 0);
	TEST_EQUAL (f.size(), 1)
END_SECTION

START_SECTION(void mergeFeatures(IsotopeWaveletTransform< PeakType > *later_iwt, const UInt RT_interleave, const UInt RT_votes_cutoff))
	NOT_TESTABLE //only via CUDA
END_SECTION

START_SECTION(DoubleReal getLinearInterpolation(const typename MSSpectrum< PeakType >::const_iterator &left_iter, const DoubleReal mz_pos, const typename MSSpectrum< PeakType >::const_iterator &right_iter))
	TEST_EQUAL((int)(iw->getLinearInterpolation(map[0].begin(), 1420.02, (map[0].begin()+1))*10),5)
END_SECTION

START_SECTION(DoubleReal getLinearInterpolation(const DoubleReal mz_a, const DoubleReal intens_a, const DoubleReal mz_pos, const DoubleReal mz_b, const DoubleReal intens_b))
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
