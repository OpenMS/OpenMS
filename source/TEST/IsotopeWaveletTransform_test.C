// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2008 -- Oliver Kohlbacher, Knut Reinert
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
MzDataFile file; file.load ("data/IsotopeWaveletTransform_test_2.mzData", map);
map.updateRanges();
IsotopeWaveletTransform<Peak1D>* trans = 0;
START_SECTION(IsotopeWaveletTransform (const DoubleReal min_mz, const DoubleReal max_mz, const UInt max_charge))
	trans = new IsotopeWaveletTransform<Peak1D> (map.getMinMZ(), map.getMaxMZ(), 1);
	TEST_NOT_EQUAL(trans, 0)
END_SECTION


START_SECTION((virtual std::multimap<DoubleReal, Box> getClosedBoxes()))
	TEST_EQUAL(trans->getClosedBoxes().size(), 0)
END_SECTION

std::vector<MSSpectrum<> > pwts (1, map[0]);
START_SECTION((virtual void getTransforms (const MSSpectrum<PeakType>& scan, std::vector<MSSpectrum<PeakType> > &transforms, const UInt max_charge, const Int mode)))
	trans->getTransforms (map[0], pwts, 1, 1);
	TEST_NOT_EQUAL (trans, 0)
END_SECTION


START_SECTION((virtual void identifyCharges (const std::vector<MSSpectrum<PeakType> >& candidates, const MSSpectrum<PeakType>& ref, const UInt scan_index, const DoubleReal ampl_cutoff=0)))
	trans->identifyCharges (pwts, map[0], 0, 5);
	TEST_NOT_EQUAL (trans, 0)
END_SECTION
	

START_SECTION((void updateBoxStates (const MSExperiment<PeakType>& map, const UInt scan_index, const UInt RT_interleave, const UInt RT_votes_cutoff)))
	trans->updateBoxStates(map, 0, 0, 0);
	trans->updateBoxStates(map, INT_MAX, 0, 0);
	TEST_NOT_EQUAL (trans, 0)
END_SECTION


START_SECTION((FeatureMap<Feature> mapSeeds2Features (const MSExperiment<PeakType>& map, const UInt max_charge, const UInt RT_votes_cutoff)))
	FeatureMap<Feature> features = trans->mapSeeds2Features (map, 1, 0);
	FeatureMap<Feature>::iterator iter;
	std::ifstream ifile ("data/IsotopeWaveletTransform.out");
	DoubleReal tmp;
	TOLERANCE_ABSOLUTE (1e-1);
	for (iter=features.begin(); iter!=features.end(); ++iter)
	{
		ifile >> tmp;
		TEST_REAL_SIMILAR (iter->getMZ(), tmp);
	}
	ifile.close();
END_SECTION


START_SECTION(UInt getPeakCutOff(const DoubleReal mass, const UInt z))
	TEST_EQUAL(trans->getPeakCutOff(2000, 1), 4)
END_SECTION


START_SECTION((virtual ~IsotopeWaveletTransform ()))
	delete(trans);
END_SECTION


END_TEST
