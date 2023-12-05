// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/FeatureFinderAlgorithm.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
	class FFA
    :public FeatureFinderAlgorithm
	{
		public:
			FFA()
        : FeatureFinderAlgorithm()
			{
			}

			~FFA() override
			{
			}

			void run() override
			{

			}

			Param getDefaultParameters() const override
			{
				Param tmp;
				tmp.setValue("bla","bluff");
				return tmp;
			}

      const PeakMap* getMap()
			{
				return this->map_;
			}

			const FeatureMap* getFeatures()
			{
				return this->features_;
			}

			const FeatureFinder* getFF()
			{
				return this->ff_;
			}
	};
}

START_TEST(FeatureFinderAlgorithm, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

FFA* ptr = nullptr;
FFA* nullPointer = nullptr;

PeakMap* map_nullPointer = nullptr;
FeatureMap*  featureMap_nullPointer = nullptr;
FeatureFinder*        ff_nullPointer = nullptr;

START_SECTION((FeatureFinderAlgorithm()))
  ptr = new FFA();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~FeatureFinderAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION((virtual void run()=0))
  FFA ffa;
	ffa.run();
END_SECTION

START_SECTION((virtual Param getDefaultParameters() const))
  FFA ffa;
	TEST_EQUAL(String(ffa.getDefaultParameters().getValue("bla").toString()),"bluff")
END_SECTION

START_SECTION((void setData(const MapType& map, FeatureMap features, FeatureFinder& ff)))
  FFA ffa;
  TEST_EQUAL(ffa.getMap(),map_nullPointer)
  TEST_EQUAL(ffa.getFeatures(),featureMap_nullPointer)
  TEST_EQUAL(ffa.getFF(),ff_nullPointer)

  PeakMap map;
	FeatureMap features;
	FeatureFinder ff;
	ffa.setData(map, features, ff);

  TEST_NOT_EQUAL(ffa.getMap(),map_nullPointer)
  TEST_NOT_EQUAL(ffa.getFeatures(),featureMap_nullPointer)
  TEST_NOT_EQUAL(ffa.getFF(),ff_nullPointer)
END_SECTION

START_SECTION((virtual void setSeeds(const FeatureMap& seeds)))
  FFA ffa;
	FeatureMap seeds;
	seeds.resize(4);
	TEST_EXCEPTION(Exception::IllegalArgument,ffa.setSeeds(seeds))
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



