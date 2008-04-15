
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
// $Maintainer: Eva Lange $
// --------------------------------------------------------------------------

#include<OpenMS/CONCEPT/ClassTest.h>

#include<OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>

using namespace std;
using namespace OpenMS;


/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithmSpectrumAlignment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithmSpectrumAlignment* ptr = 0;
CHECK((MapAlignmentAlgorithmSpectrumAlignment()))
	ptr = new MapAlignmentAlgorithmSpectrumAlignment();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~MapAlignmentAlgorithmSpectrumAlignment()))
	delete ptr;
RESULT

CHECK(static MapAlignmentAlgorithm* create())
	TEST_NOT_EQUAL(MapAlignmentAlgorithmSpectrumAlignment::create(),0)
RESULT

CHECK(static String getProductName())
	TEST_EQUAL(MapAlignmentAlgorithmSpectrumAlignment::getProductName(), "spectrum_alignment")
RESULT

CHECK(void  alignPeakMaps(std::vector< MSExperiment<> >&))
  MapAlignmentAlgorithmSpectrumAlignment ma;
  std::vector< MSExperiment<> > maps;
  ma.alignPeakMaps(maps);
RESULT

CHECK(void alignFeatureMaps(std::vector< FeatureMap<> >&))
  MapAlignmentAlgorithmSpectrumAlignment ma;
  std::vector< FeatureMap<> > maps;
  TEST_EXCEPTION(Exception::NotImplemented, ma.alignFeatureMaps(maps));
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
