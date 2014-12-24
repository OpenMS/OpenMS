// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Clemens Groepl $
// $Authors: Katharina Albers $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithm.h>
///////////////////////////

#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmPrecision.h>
#include <OpenMS/ANALYSIS/MAPMATCHING/MapAlignmentEvaluationAlgorithmRecall.h>

#include <OpenMS/CONCEPT/Factory.h>
#include <OpenMS/KERNEL/Feature.h>

using namespace OpenMS;
using namespace std;

namespace OpenMS
{
	class MAEA
	 : public MapAlignmentEvaluationAlgorithm
	{
		public:
			void evaluate(const ConsensusMap&, const ConsensusMap&, const double&, const double&, const Peak2D::IntensityType&, const bool use_charge, double& real)
			{
				bool x = use_charge;
				x=!x;
				real = 1.5;
			}
	};
}

START_TEST(MapAlignmentEvaluation, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MAEA* ptr = 0;
MAEA* nullPointer = 0;
START_SECTION((MapAlignmentEvaluationAlgorithm()))
	ptr = new MAEA();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~MapAlignmentEvaluationAlgorithm()))
	delete ptr;
END_SECTION

START_SECTION((virtual void evaluate(const ConsensusMap &conensus_map_in, const ConsensusMap &consensus_map_gt, const double &rt_dev, const double &mz_dev, const Peak2D::IntensityType &int_dev, const bool use_charge, double &out)=0))
	MAEA maea;
	ConsensusMap map1;
	ConsensusMap map2;
	double rt_dev, mz_dev;
	Peak2D::IntensityType int_dev;
	double real;
	maea.evaluate(map1, map2, rt_dev, mz_dev, int_dev, false, real);
	TEST_EQUAL(real, 1.5)
END_SECTION

START_SECTION((bool isSameHandle(const FeatureHandle &lhs, const FeatureHandle &rhs, const double &rt_dev, const double &mz_dev, const Peak2D::IntensityType &int_dev, const bool use_charge)))
{
	Feature tmp_feature;
	tmp_feature.setRT(100);
	tmp_feature.setMZ(555);
	tmp_feature.setIntensity(200.0f);
	tmp_feature.setCharge(3);
	tmp_feature.setUniqueId(1);

  Feature tmp_feature2;
	tmp_feature2.setRT(101);
	tmp_feature2.setMZ(556);
	tmp_feature2.setIntensity(1199.0f);
	tmp_feature2.setCharge(4);
  tmp_feature2.setUniqueId(2);

	FeatureHandle a(0,tmp_feature);
	FeatureHandle b(0,tmp_feature2);

	MAEA maea;

	TEST_EQUAL(maea.isSameHandle(a, b, 2, 1.5, 1000, false), true);
	TEST_EQUAL(maea.isSameHandle(a, b, 2, 1.5, 1000, true), false);

	tmp_feature2.setCharge(3); // now charge is equal
	FeatureHandle b2(0,tmp_feature2);

	TEST_EQUAL(maea.isSameHandle(a, b2, 2, 1.5, 1000, false), true);
	TEST_EQUAL(maea.isSameHandle(a, b2, 2, 1.5, 1000, true), true);

}
END_SECTION

START_SECTION((static void registerChildren()))
{
	TEST_STRING_EQUAL(Factory<MapAlignmentEvaluationAlgorithm>::registeredProducts()[0],MapAlignmentEvaluationAlgorithmPrecision::getProductName());
	TEST_STRING_EQUAL(Factory<MapAlignmentEvaluationAlgorithm>::registeredProducts()[1],MapAlignmentEvaluationAlgorithmRecall::getProductName());
	TEST_EQUAL(Factory<MapAlignmentEvaluationAlgorithm>::registeredProducts().size(),2)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
