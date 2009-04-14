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

START_SECTION(IsotopeWaveletTransform(const DoubleReal min_mz, const DoubleReal max_mz, const UInt max_charge, const UInt max_scan_size=0))
	IsotopeWaveletTransform<Peak1D> iw (800, 1000, 2);
	TEST_EQUAL (iw.getClosedBoxes().size(), 0)
END_SECTION

START_SECTION(~IsotopeWaveletTransform())
	IsotopeWaveletTransform<Peak1D>* iw2 = new IsotopeWaveletTransform<Peak1D> (800, 1000, 2);
	delete (iw2);
	iw2 = NULL;
END_SECTION

START_SECTION(void mergeFeatures(IsotopeWaveletTransform< PeakType > *later_iwt, const UInt RT_votes_cutoff))
//unable to test here something, since this is only testable with CUDA and TBB
	TEST_EQUAL (0, 0)
END_SECTION

END_TEST
