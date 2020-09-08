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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussModel.h>

///////////////////////////

START_TEST(GaussModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace OpenMS::Math;
using std::stringstream;

// default ctor
GaussModel* ptr = nullptr;
GaussModel* nullPointer = nullptr;
START_SECTION((GaussModel()))
	ptr = new GaussModel();
  TEST_EQUAL(ptr->getName(), "GaussModel")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~GaussModel()))
	delete ptr;
END_SECTION


START_SECTION((static const String getProductName()))
	TEST_EQUAL(GaussModel::getProductName(),"GaussModel")
	TEST_EQUAL(GaussModel().getProductName(),"GaussModel")
END_SECTION

START_SECTION(static BaseModel<1>* create())
	BaseModel<1>* ptr = GaussModel::create();
	TEST_EQUAL(ptr->getName(), "GaussModel")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// assignment operator
START_SECTION((virtual GaussModel& operator=(const GaussModel &source)))
	GaussModel gm1;
	BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	gm1.setScalingFactor(10.0);
	gm1.setInterpolationStep(0.3);
	
	Param tmp;
	tmp.setValue("bounding_box:min",678.9);
	tmp.setValue("bounding_box:max",789.0 );
	tmp.setValue("statistics:variance",stat.variance() );
	tmp.setValue("statistics:mean",stat.mean() );			
	gm1.setParameters(tmp);

  GaussModel gm2;
  gm2 = gm1;

  GaussModel gm3;
	gm3.setScalingFactor(10.0);
	gm3.setInterpolationStep(0.3);
	gm3.setParameters(tmp);

  gm1 = GaussModel();
	TEST_EQUAL(gm3.getParameters(), gm2.getParameters())
END_SECTION

// copy ctor
START_SECTION((GaussModel(const GaussModel& source)))
	GaussModel gm1;
	BasicStatistics<>  stat;
	stat.setMean(680.1);
	stat.setVariance(2.0);
	gm1.setScalingFactor(10.0);
	gm1.setInterpolationStep(0.3);
	
	Param tmp;
	tmp.setValue("bounding_box:min",678.9);
	tmp.setValue("bounding_box:max",789.0 );
	tmp.setValue("statistics:variance",stat.variance() );
	tmp.setValue("statistics:mean",stat.mean() );			
	gm1.setParameters(tmp);	

	GaussModel gm2(gm1);
  GaussModel gm3;
	gm3.setScalingFactor(10.0);
	gm3.setInterpolationStep(0.3);
	gm3.setParameters(tmp);

  	gm1 = GaussModel();
	TEST_EQUAL(gm3.getParameters(), gm2.getParameters())
END_SECTION

START_SECTION([EXTRA] DefaultParamHandler::setParameters(...))
	TOLERANCE_ABSOLUTE(0.001)
	GaussModel gm1;
	
	gm1.setScalingFactor(10.0);
	
	Param tmp;
	tmp.setValue("bounding_box:min",678.9);
	tmp.setValue("bounding_box:max",789.0 );
	tmp.setValue("statistics:variance",2.0);
	tmp.setValue("statistics:mean",679.1);			
	gm1.setParameters(tmp);
	gm1.setOffset(680.0);

	TEST_REAL_SIMILAR(gm1.getCenter(), 680.2)

	GaussModel gm2;
	gm2.setParameters( gm1.getParameters() );

	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
	gm1.getSamples(dpa1);
	gm2.getSamples(dpa2);

	TOLERANCE_ABSOLUTE(0.0000001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
END_SECTION

START_SECTION(([EXTRA]void setParam(const BasicStatistics&,CoordinateType,CoordinateType)))
	GaussModel gm1;
	BasicStatistics<>  stat;
	stat.setMean(0.0);
	stat.setVariance(1.0);
	gm1.setInterpolationStep(0.001);
	
	Param tmp;
	tmp.setValue("bounding_box:min",-4.0);
	tmp.setValue("bounding_box:max",4.0 );
	tmp.setValue("statistics:variance",stat.variance() );
	tmp.setValue("statistics:mean",stat.mean() );			
	gm1.setParameters(tmp);	

	TEST_REAL_SIMILAR(gm1.getCenter(), 0.0)

	TOLERANCE_ABSOLUTE(0.001)
	TEST_REAL_SIMILAR(gm1.getIntensity(-1.0), 0.24197072);
	TEST_REAL_SIMILAR(gm1.getIntensity(0.0), 0.39894228);
	TEST_REAL_SIMILAR(gm1.getIntensity(1.0), 0.24197072);
	TEST_REAL_SIMILAR(gm1.getIntensity(2.0), 0.05399097);

	gm1.setInterpolationStep(0.2);
	gm1.setSamples();

	TEST_REAL_SIMILAR(gm1.getIntensity(-1.0), 0.24197072);
	TEST_REAL_SIMILAR(gm1.getIntensity(0.0), 0.39894228);
	TEST_REAL_SIMILAR(gm1.getIntensity(1.0), 0.24197072);
	TEST_REAL_SIMILAR(gm1.getIntensity(2.0), 0.05399097);

	gm1.setScalingFactor(10.0);
	gm1.setSamples();

	TEST_REAL_SIMILAR(gm1.getIntensity(-1.0), 2.4197072);
	TEST_REAL_SIMILAR(gm1.getIntensity(0.0), 3.9894228);
	TEST_REAL_SIMILAR(gm1.getIntensity(1.0), 2.4197072);
	TEST_REAL_SIMILAR(gm1.getIntensity(2.0), 0.5399097);
END_SECTION


START_SECTION((void setOffset(CoordinateType offset)))

	TOLERANCE_ABSOLUTE(0.001)
	GaussModel gm1;
	
	gm1.setScalingFactor(10.0);
	
	Param tmp;
	tmp.setValue("bounding_box:min",678.9);
	tmp.setValue("bounding_box:max",789.0 );
	tmp.setValue("statistics:variance",2.0);
	tmp.setValue("statistics:mean",679.1);			
	gm1.setParameters(tmp);
	gm1.setOffset(680.0);

	TEST_REAL_SIMILAR(gm1.getCenter(), 680.2)

END_SECTION

START_SECTION( CoordinateType getCenter() const )
	TOLERANCE_ABSOLUTE(0.001)
	GaussModel gm1;
	
	gm1.setScalingFactor(10.0);
	
	Param tmp;
	tmp.setValue("bounding_box:min",650.0);
	tmp.setValue("bounding_box:max",750.0 );
	tmp.setValue("statistics:variance",2.0);
	tmp.setValue("statistics:mean",679.1);			
	gm1.setParameters(tmp);

	TEST_REAL_SIMILAR(gm1.getCenter(), 679.1)	
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
