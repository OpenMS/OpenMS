
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
// $Maintainer: Vipul Patel $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentAlgorithmSpectrumAlignment.h>

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

CHECK((static MapAlignmentAlgorithm* create()))
	TEST_NOT_EQUAL(MapAlignmentAlgorithmSpectrumAlignment::create(),0)
RESULT

CHECK((static String getProductName()))
	TEST_EQUAL(MapAlignmentAlgorithmSpectrumAlignment::getProductName(), "spectrum_alignment")
RESULT

CHECK((virtual void alignPeakMaps(std::vector< MSExperiment<> > &, std::vector< TransformationDescription > &)))
  MapAlignmentAlgorithmSpectrumAlignment ma;
  std::vector< MSExperiment<> > maps;
	PeakMap map1;
	PeakMap map2;	
	for(UInt i= 0; i< 15; ++i)
	{
		for(UInt j =1 ; j < 5; ++j)
		{
			PeakSpectrum spectrum;
			spectrum.setRT(i);
			spectrum.setMSLevel(j);
		
			for (Real mz=500.0; mz<=900; mz+=100.0)
		    { 
				Peak1D peak;
				peak.setMZ(mz+i);
				peak.setIntensity(mz+i);
				spectrum.push_back(peak);  
		    }
		    map1.push_back(spectrum);
		}
	}
	for(UInt i= 0; i< 15; ++i)
		{
			for(UInt j =1 ; j < 5; ++j)
			{
				PeakSpectrum spectrum;
				spectrum.setRT(i*1.2+200);
				spectrum.setMSLevel(j);
			
				for (Real mz=500.0; mz<=900; mz+=100.0)
			    { 
					Peak1D peak;
					peak.setMZ(mz+i);
					peak.setIntensity(mz+i);
					spectrum.push_back(peak);  
			    }
			    map2.push_back(spectrum);
			}
		}
	
	maps.push_back(map1);
	maps.push_back(map2);
	std::vector<TransformationDescription> transformations;
  ma.alignPeakMaps(maps,transformations);
  Int counter =0;
	maps[0].updateRanges(-1);
	maps[1].updateRanges(-1);
  for(UInt i=0; i< maps[0].size(); ++i)
  {
		if((maps[0])[i].getMSLevel() <2)
		{
	  	if((maps[0])[i].getRT() != (maps[1])[i].getRT())
	  	{
	  		++counter;
	  	}
		}
  }
	TEST_REAL_EQUAL(counter, 0)
RESULT

CHECK([EXTRA] void alignFeatureMaps(std::vector< FeatureMap<> >&))
  MapAlignmentAlgorithmSpectrumAlignment ma;
  std::vector< FeatureMap<> > maps;
	std::vector<TransformationDescription> transformations;
  TEST_EXCEPTION(Exception::NotImplemented, ma.alignFeatureMaps(maps,transformations));
RESULT

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

