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
// $Maintainer: Clemens Groepl$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm_impl.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
	template <class PeakType, class FeatureType>
	class FFA
		:public FeatureFinderAlgorithm<PeakType,FeatureType>
	{
		public:
			FFA()
				: FeatureFinderAlgorithm<PeakType,FeatureType>()
			{
			}
			
			~FFA()
			{
			}
			
			virtual void run()
			{
				
			}
			
			virtual Param getDefaultParameters() const
			{
				Param tmp;
				tmp.setValue("bla","bluff");
				return tmp;
			}
			
			const MSExperiment<PeakType>* getMap()
			{
				return this->map_;
			}
		
			const FeatureMap<Feature>* getFeatures()
			{
				return this->features_;
			}
		
			const FeatureFinder* getFF()
			{
				return this->ff_;
			}
	};
}

START_TEST(FeatureFinderAlgorithm, "$Id FeatureFinder_test.C 139 2006-07-14 10:08:39Z ole_st $")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FFA<Peak1D,Feature>* ptr = 0;
FFA<Peak1D,Feature>* nullPointer = 0;

MSExperiment<Peak1D>* map_nullPointer = 0;
FeatureMap<Feature>*  featureMap_nullPointer = 0;
FeatureFinder*        ff_nullPointer = 0;

START_SECTION((FeatureFinderAlgorithm()))
	ptr = new FFA<Peak1D,Feature>();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureFinderAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION([EXTRA] FeatureFinderAlgorithmPicked() - with RichPeak1D)
	FeatureFinderAlgorithmPicked<RichPeak1D,Feature> ffa;
END_SECTION

START_SECTION((virtual void run()=0))
	FFA<Peak1D,Feature> ffa;
	ffa.run();
END_SECTION

START_SECTION((virtual Param getDefaultParameters() const))
	FFA<Peak1D,Feature> ffa;
	TEST_EQUAL(String(ffa.getDefaultParameters().getValue("bla")),"bluff")
END_SECTION

START_SECTION((void setData(const MapType& map, FeatureMapType& features, FeatureFinder& ff)))
	FFA<Peak1D,Feature> ffa;
  TEST_EQUAL(ffa.getMap(),map_nullPointer)
  TEST_EQUAL(ffa.getFeatures(),featureMap_nullPointer)
  TEST_EQUAL(ffa.getFF(),ff_nullPointer)
	
	MSExperiment<Peak1D> map;
	FeatureMap<Feature> features;
	FeatureFinder ff;
	ffa.setData(map, features, ff);

  TEST_NOT_EQUAL(ffa.getMap(),map_nullPointer)
  TEST_NOT_EQUAL(ffa.getFeatures(),featureMap_nullPointer)
  TEST_NOT_EQUAL(ffa.getFF(),ff_nullPointer)
END_SECTION

START_SECTION((virtual void setSeeds(const FeatureMapType& seeds)))
	FFA<Peak1D,Feature> ffa;
	FeatureMap<Feature> seeds;
	seeds.resize(4);
	TEST_EXCEPTION(Exception::IllegalArgument,ffa.setSeeds(seeds))	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



