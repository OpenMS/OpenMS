// -*- Mode: C++; tab-width: 2; -*-
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



IsotopeWaveletTransform<Peak2D>* trans2 = 0;
CHECK(IsotopeWaveletTransform)
	trans2 = new IsotopeWaveletTransform<Peak2D> (800, 4000, 2);
	TEST_NOT_EQUAL(trans2, 0)
	TEST_EQUAL(IsotopeWavelet::getMaxCharge(), 2)
RESULT


CHECK(getClosedBoxes)
	TEST_EQUAL(trans2->getClosedBoxes().size(), 0)
RESULT


CHECK(~IsotopeWaveletTransform)
	delete(trans2);
RESULT


MSExperiment<Peak1D> map;
MzDataFile file; file.load ("data/IsotopeWaveletTransform_test.mzData", map);
map.updateRanges();
std::cout << map.getMinMZ() << "\t" << map.getMaxMZ() << std::endl;
IsotopeWaveletTransform<Peak1D>* trans = 0;
CHECK(IsotopeWaveletTransform)
	trans = new IsotopeWaveletTransform<Peak1D> (map.getMinMZ(), map.getMaxMZ(), 1);
	TEST_NOT_EQUAL(trans, 0)
RESULT


std::vector<MSSpectrum<Peak1D> > pwts (1, map[0]);
CHECK(getTransforms)
	trans->getTransforms (map[0], pwts, 1, 1);
RESULT


CHECK(identifyCharges)
	trans->identifyCharges (pwts, map[0], 0, 5);
RESULT
	

CHECK(updateBoxStates)
	trans->updateBoxStates(0, 0, 0);
RESULT

CHECK(updateBoxStates)
	trans->updateBoxStates(INT_MAX, 0, 0);
RESULT

CHECK(mapSeeds2Features)
	FeatureMap<Feature> features = trans->mapSeeds2Features (map, 1, 0);
	FeatureMap<Feature>::iterator iter;
	std::ifstream ifile ("data/IsotopeWaveletTransform.out");
	DoubleReal tmp;
	PRECISION (1e-1);
	for (iter=features.begin(); iter!=features.end(); ++iter)
	{
		ifile >> tmp;
		TEST_REAL_EQUAL (iter->getMZ(), tmp);
		ifile >> tmp;
		TEST_REAL_EQUAL (iter->getIntensity(), tmp);
	}
	ifile.close();
RESULT


CHECK(~IsotopeWaveletTransform)
	delete(trans);
RESULT


END_TEST
