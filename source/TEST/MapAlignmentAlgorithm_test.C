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
// $Maintainer: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>

using namespace std;
using namespace OpenMS;


/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithm* ptr = 0;
CHECK((MapAlignmentAlgorithm()))
	ptr = new MapAlignmentAlgorithm();
	TEST_NOT_EQUAL(ptr, 0)
RESULT

CHECK((virtual ~MapAlignmentAlgorithm()))
	delete ptr;
RESULT

CHECK((virtual void alignPeakMaps(std::vector< MSExperiment<> > &, std::vector< TransformationDescription > &)))
  MapAlignmentAlgorithm ma;
  std::vector< MSExperiment<> > maps;
  std::vector<TransformationDescription> transformations;
  TEST_EXCEPTION(Exception::NotImplemented, ma.alignPeakMaps(maps,transformations));
RESULT

CHECK((virtual void alignFeatureMaps(std::vector< FeatureMap<> > &, std::vector< TransformationDescription > &)))
  MapAlignmentAlgorithm ma;
  std::vector< FeatureMap<> > maps;
  std::vector<TransformationDescription> transformations;
  TEST_EXCEPTION(Exception::NotImplemented, ma.alignFeatureMaps(maps,transformations));
RESULT

CHECK((static void registerChildren()))
{
  // I do not know why the classes show up in this particular order.
  TEST_STRING_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts()[0],MapAlignmentAlgorithmPoseClustering::getProductName());
  TEST_STRING_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts()[1],MapAlignmentAlgorithmSpectrumAlignment::getProductName());
  TEST_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts().size(),2)
}
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
