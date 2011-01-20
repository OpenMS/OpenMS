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
// $Maintainer: Clemens Groepl $
// $Authors: Marc Sturm, Clemens Groepl $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithm.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmApplyGivenTrafo.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmPoseClustering.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmIdentification.h>

using namespace std;
using namespace OpenMS;


/////////////////////////////////////////////////////////////

START_TEST(MapAlignmentAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


MapAlignmentAlgorithm* ptr = 0;
START_SECTION((MapAlignmentAlgorithm()))
	ptr = new MapAlignmentAlgorithm();
	TEST_NOT_EQUAL(ptr, 0)
END_SECTION

START_SECTION((virtual ~MapAlignmentAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION((virtual void alignPeakMaps(std::vector< MSExperiment<> > &, std::vector< TransformationDescription > &)))
  MapAlignmentAlgorithm ma;
  std::vector< MSExperiment<> > maps;
  std::vector<TransformationDescription> transformations;
  TEST_EXCEPTION(Exception::NotImplemented, ma.alignPeakMaps(maps,transformations));
END_SECTION

START_SECTION((virtual void alignFeatureMaps(std::vector< FeatureMap<> > &, std::vector< TransformationDescription > &)))
  MapAlignmentAlgorithm ma;
  std::vector< FeatureMap<> > maps;
  std::vector<TransformationDescription> transformations;
  TEST_EXCEPTION(Exception::NotImplemented, ma.alignFeatureMaps(maps,transformations));
END_SECTION

START_SECTION((virtual void alignPeptideIdentifications(std::vector< std::vector< PeptideIdentification > >&, std::vector<TransformationDescription>&)))
  MapAlignmentAlgorithm ma;
  std::vector< std::vector< PeptideIdentification > > maps;
  std::vector<TransformationDescription> transformations;
  TEST_EXCEPTION(Exception::NotImplemented, ma.alignPeptideIdentifications(maps,transformations));
END_SECTION

START_SECTION((static void registerChildren()))
{
  // I do not know why the classes show up in this particular order (sorted by name?).
	TEST_STRING_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts()[0],MapAlignmentAlgorithmApplyGivenTrafo::getProductName());
	TEST_STRING_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts()[1],MapAlignmentAlgorithmIdentification::getProductName());	
	TEST_STRING_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts()[2],MapAlignmentAlgorithmPoseClustering::getProductName());
	TEST_STRING_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts()[3],MapAlignmentAlgorithmSpectrumAlignment::getProductName());
  TEST_EQUAL(Factory<MapAlignmentAlgorithm>::registeredProducts().size(),4)
}
END_SECTION

START_SECTION((virtual void setReference(Size, const String&)))
{
	MapAlignmentAlgorithm ma;
	ma.setReference(); // no exception, nothing happens
	TEST_EXCEPTION(Exception::InvalidParameter, ma.setReference(1));
	TEST_EXCEPTION(Exception::InvalidParameter, ma.setReference(0, "test"));
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
