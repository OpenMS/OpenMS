// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2009 -- Oliver Kohlbacher, Knut Reinert
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

#include <fstream>

using namespace OpenMS;
using namespace std;

START_TEST(IsotopeWaveletTransform, "$Id$")

MSExperiment<> map;
MzDataFile file; file.load (OPENMS_GET_TEST_DATA_PATH("IsotopeWaveletTransform_test.mzData"), map);
map.updateRanges();
IsotopeWaveletTransform<Peak1D>* iw = NULL;
MSSpectrum<Peak1D>* spec = new MSSpectrum<Peak1D> (map[0]);

START_SECTION(IsotopeWaveletTransform(const DoubleReal min_mz, const DoubleReal max_mz, const UInt max_charge, const UInt max_scan_size=0))
	iw = new IsotopeWaveletTransform<Peak1D> (map[0].begin()->getMZ(), (map[0].end()-1)->getMZ(), 1);
	TEST_NOT_EQUAL (iw, NULL)
END_SECTION

START_SECTION(void initializeScan(const MSSpectrum< PeakType > &c_ref))
	iw->initializeScan (map[0]);
	NOT_TESTABLE
END_SECTION

START_SECTION(void getTransform(MSSpectrum<PeakType> &c_trans, const MSSpectrum<PeakType> &c_ref, const UInt c))
	iw->getTransform (*spec, map[0], 0);
	TEST_EQUAL (*spec!= map[0], true)
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

START_SECTION(~IsotopeWaveletTransform())
	delete (spec);
	delete (iw);
END_SECTION

END_TEST
