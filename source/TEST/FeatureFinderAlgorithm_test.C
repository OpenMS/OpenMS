// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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



