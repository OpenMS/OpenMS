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
// $Maintainer: Timo Sachsenberg $
// $Authors: Clemens Groepl, Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MAPMATCHING/BaseSuperimposer.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

class TestSuperimposer 
	: public BaseSuperimposer
{
  public:
	TestSuperimposer() 
		: BaseSuperimposer()
	{ 
		check_defaults_ = false; 
	}

	void run(const ConsensusMap& , const ConsensusMap& , TransformationDescription& transformation) override
	{
		Param params;
		params.setValue("slope",1.1);
		params.setValue("intercept", 5.0);
		transformation.fitModel("linear", params);
	}
};

START_TEST(BaseSuperimposer, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TestSuperimposer* ptr = nullptr;
TestSuperimposer* nullPointer = nullptr;
START_SECTION((BaseSuperimposer()))
	ptr = new TestSuperimposer();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~BaseSuperimposer()))
	delete ptr;
END_SECTION

START_SECTION((virtual void run(const ConsensusMap& map_model, const ConsensusMap& map_scene, TransformationDescription& transformation)=0))
{
  TransformationDescription transformation;
  TestSuperimposer si;
	std::vector<ConsensusMap> maps;
	maps.resize(2);
  si.run(maps[0], maps[1], transformation);
  TEST_STRING_EQUAL(transformation.getModelType(), "linear");
  Param params = transformation.getModelParameters();
  TEST_REAL_SIMILAR(params.getValue("slope"), 1.1)
  TEST_REAL_SIMILAR(params.getValue("intercept"), 5.0)
}
END_SECTION

START_SECTION((static void registerChildren()))
  NOT_TESTABLE
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



