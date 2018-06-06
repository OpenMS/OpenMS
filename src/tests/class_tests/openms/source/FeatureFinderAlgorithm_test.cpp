// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
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
	TEST_EQUAL(String(ffa.getDefaultParameters().getValue("bla")),"bluff")
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



