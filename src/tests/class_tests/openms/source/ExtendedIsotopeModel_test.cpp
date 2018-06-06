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
// $Maintainer: Timo Sachsenberg  $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeModel.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <sstream>


///////////////////////////

START_TEST(ExtendedIsotopeModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

// default ctor
ExtendedIsotopeModel* ptr = nullptr;
ExtendedIsotopeModel* nullPointer = nullptr;
START_SECTION((ExtendedIsotopeModel()))
	ptr = new ExtendedIsotopeModel();
  TEST_EQUAL(ptr->getName(), "ExtendedIsotopeModel")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

// destructor
START_SECTION((virtual ~ExtendedIsotopeModel()))
	delete ptr;
END_SECTION

START_SECTION(static BaseModel<1>* create())
	BaseModel<1>* ptr = ExtendedIsotopeModel::create();
	TEST_EQUAL(ptr->getName(), "ExtendedIsotopeModel")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((static const String getProductName()))
	TEST_EQUAL(ExtendedIsotopeModel::getProductName(),"ExtendedIsotopeModel")
	TEST_EQUAL(ExtendedIsotopeModel().getName(),"ExtendedIsotopeModel")
END_SECTION

// assignment operator
START_SECTION((virtual ExtendedIsotopeModel& operator=(const ExtendedIsotopeModel &source)))
	ExtendedIsotopeModel im1;
	
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:stdev",0.8);
  tmp.setValue("isotope:monoisotopic_mz", 670.5);
	im1.setParameters(tmp);

  ExtendedIsotopeModel im2;
  im2 = im1;

  ExtendedIsotopeModel im3;
	im3.setParameters(tmp);

  im1 = ExtendedIsotopeModel();
	TEST_EQUAL(im3.getParameters(), im2.getParameters())
END_SECTION

// copy ctor
START_SECTION((ExtendedIsotopeModel(const ExtendedIsotopeModel& source)))
	ExtendedIsotopeModel im1;
	
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:stdev",0.8);
  tmp.setValue("isotope:monoisotopic_mz", 670.5);
	im1.setParameters(tmp);

	ExtendedIsotopeModel im2(im1);
  ExtendedIsotopeModel im3;
	im3.setParameters(tmp);

  im1 = ExtendedIsotopeModel();
	TEST_EQUAL(im3.getParameters(), im2.getParameters())
END_SECTION

START_SECTION([EXTRA] DefaultParamHandler::setParameters(...))
	TOLERANCE_ABSOLUTE(0.001)
	ExtendedIsotopeModel im1;
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:stdev",0.8);
  tmp.setValue("isotope:monoisotopic_mz", 670.5);
	im1.setParameters(tmp);

	ExtendedIsotopeModel im2;
	im2.setParameters(im1.getParameters());

	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
	im1.getSamples(dpa1);
	im2.getSamples(dpa2);

	TOLERANCE_ABSOLUTE(0.00001)
	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
END_SECTION

START_SECTION(UInt getCharge() )
	// can only reliably be tested after fitting, only sanity check here
	ExtendedIsotopeModel im1;
	TEST_EQUAL(im1.getCharge() == 1, true)		// default charge is 1
END_SECTION

START_SECTION( CoordinateType getCenter() const )
	// can only reliably be tested after fitting, only sanity check here
	ExtendedIsotopeModel im1;
  TEST_EQUAL(im1.getCenter() == 1, true) // default charge is 1 and hence center mus be 1
END_SECTION

START_SECTION( void setOffset(CoordinateType offset) )
	TOLERANCE_ABSOLUTE(0.1)
	ExtendedIsotopeModel im1;
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:stdev",0.8);
  tmp.setValue("isotope:monoisotopic_mz", 670.5);
	im1.setParameters(tmp);
	im1.setOffset( 673.5 );
	
	ExtendedIsotopeModel im2;
	im2.setParameters(im1.getParameters());
	im2.setOffset( 673.5 );
	
	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
	im1.getSamples(dpa1);
	im2.getSamples(dpa2);

	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}
END_SECTION

START_SECTION( CoordinateType getOffset() )
	TOLERANCE_ABSOLUTE(0.1)
	ExtendedIsotopeModel im1;
	Param tmp;
	tmp.setValue("charge", 3);
	tmp.setValue("isotope:stdev",0.8);
  tmp.setValue("isotope:monoisotopic_mz", 670.5);
	im1.setParameters(tmp);
	im1.setOffset( 673.5 );
	
	ExtendedIsotopeModel im2;
	im2.setParameters(im1.getParameters());
	im2.setOffset( im1.getOffset() );
	
	std::vector<Peak1D> dpa1;
	std::vector<Peak1D> dpa2;
	im1.getSamples(dpa1);
	im2.getSamples(dpa2);

	TEST_EQUAL(dpa1.size(),dpa2.size())
	ABORT_IF(dpa1.size()!=dpa2.size());
	for (Size i=0; i<dpa1.size(); ++i)
	{
		TEST_REAL_SIMILAR(dpa1[i].getPosition()[0],dpa2[i].getPosition()[0])
		TEST_REAL_SIMILAR(dpa1[i].getIntensity(),dpa2[i].getIntensity())
	}	
END_SECTION

START_SECTION((void setSamples()))
{
  // dummy subtest
	TEST_EQUAL(1,1)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
