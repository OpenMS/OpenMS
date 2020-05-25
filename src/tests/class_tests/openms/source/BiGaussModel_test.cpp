// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/BiGaussModel.h>


///////////////////////////

START_TEST(BiGaussModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using std::stringstream;

// default ctor
BiGaussModel* ptr = nullptr;
BiGaussModel* nullPointer = nullptr;
START_SECTION((BiGaussModel()))
	ptr = new BiGaussModel();
        TEST_EQUAL(ptr->getName(), "BiGaussModel")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~BiGaussModel()))
	delete ptr;
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(BiGaussModel::getProductName(),"BiGaussModel")
	TEST_EQUAL(BiGaussModel().getName(),"BiGaussModel")
END_SECTION

START_SECTION( static BaseModel<1>* create() )
	BaseModel<1>* ptr = BiGaussModel::create();
	TEST_EQUAL(ptr->getName(), "BiGaussModel")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// assignment operator
START_SECTION((virtual BiGaussModel& operator=(const BiGaussModel &source)))
	BiGaussModel bgm1;
	bgm1.setScalingFactor(10.0);
	bgm1.setInterpolationStep(0.3);
	
	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm1.setParameters(tmp);

  BiGaussModel bgm2;
  bgm2 = bgm1;

  BiGaussModel bgm3;
	bgm3.setScalingFactor(10.0);
	bgm3.setInterpolationStep(0.3);
	bgm3.setParameters(tmp);

  bgm1 = BiGaussModel();
	TEST_EQUAL(bgm3.getParameters(), bgm2.getParameters())
END_SECTION

// copy ctor
START_SECTION((BiGaussModel(const BiGaussModel& source)))
	BiGaussModel bgm1;
	BasicStatistics<>  stat;
	bgm1.setScalingFactor(10.0);
	bgm1.setInterpolationStep(0.3);

	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm1.setParameters(tmp);

	BiGaussModel bgm2(bgm1);
  BiGaussModel bgm3;
	bgm3.setScalingFactor(10.0);
	bgm3.setInterpolationStep(0.3);
	bgm3.setParameters(tmp);
  bgm1 = BiGaussModel();
	TEST_EQUAL(bgm3.getParameters(), bgm2.getParameters())
END_SECTION

START_SECTION([EXTRA] DefaultParamHandler::setParameters(...))
	TOLERANCE_ABSOLUTE(0.001)
	BiGaussModel bgm1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm1.setParameters(tmp);
	bgm1.setOffset(680.0);

	BiGaussModel bgm2;
	bgm2.setParameters(bgm1.getParameters());
	TEST_REAL_SIMILAR(bgm1.getCenter(), 681.2)

	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
	bgm1.getSamples(dpa1);
	bgm2.getSamples(dpa2);

	TOLERANCE_ABSOLUTE(0.0001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
END_SECTION

START_SECTION((void setOffset(CoordinateType offset)))
	BiGaussModel bgm1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm1.setParameters(tmp);
	bgm1.setOffset(680.9);

	BiGaussModel bgm2;
	tmp.setValue("bounding_box:min", 680.9);
	tmp.setValue("bounding_box:max", 791.0);
	tmp.setValue("statistics:mean", 682.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm2.setParameters(tmp);

	TEST_EQUAL(bgm1.getParameters(), bgm2.getParameters())
	TEST_REAL_SIMILAR(bgm1.getCenter(), bgm2.getCenter())
	TEST_REAL_SIMILAR(bgm1.getCenter(), 682.1)

	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
	bgm1.getSamples(dpa1);
	bgm2.getSamples(dpa2);

	TOLERANCE_ABSOLUTE(0.001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}

	tmp.setValue("bounding_box:min", -4.0);
	tmp.setValue("bounding_box:max", 4.001);
	tmp.setValue("statistics:mean", 0.0 );
	tmp.setValue("statistics:variance1",  0.81);
	tmp.setValue("statistics:variance2", 0.81 );
	bgm1.setParameters(tmp);
	bgm1.setOffset(0.123);
	TEST_REAL_SIMILAR(bgm1.getCenter(), 4.123)

	TOLERANCE_ABSOLUTE(0.001)
	TEST_REAL_SIMILAR(bgm1.getIntensity(4.123), 0.4432692);
	TEST_REAL_SIMILAR(bgm1.getIntensity(4.223), bgm1.getIntensity(4.023));
	TEST_REAL_SIMILAR(bgm1.getIntensity(3.123), bgm1.getIntensity(5.123));

END_SECTION

START_SECTION( CoordinateType getCenter() const )
	
	// already test above, but just for the sake of it
	TOLERANCE_ABSOLUTE(0.001)
	BiGaussModel bgm1;
	
	Param tmp;
	tmp.setValue("bounding_box:min", 	678.9);
	tmp.setValue("bounding_box:max", 789.0);
	tmp.setValue("statistics:mean", 680.1 );
	tmp.setValue("statistics:variance1",  2.0);
	tmp.setValue("statistics:variance2", 5.0 );
	bgm1.setParameters(tmp);
	bgm1.setOffset(680.0);
	TEST_REAL_SIMILAR(bgm1.getCenter(), 681.2)

END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
