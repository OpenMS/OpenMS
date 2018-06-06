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
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/GaussFitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GaussFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

GaussFitter1D* ptr = nullptr;
GaussFitter1D* nullPointer = nullptr;
START_SECTION(GaussFitter1D())
{
	ptr = new GaussFitter1D();
  TEST_EQUAL(ptr->getName(), "GaussFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((GaussFitter1D(const  GaussFitter1D &source)))
	GaussFitter1D gf1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
	gf1.setParameters(param);

	GaussFitter1D gf2(gf1);
  GaussFitter1D gf3;
	gf3.setParameters(param);
  gf1 = GaussFitter1D();
	TEST_EQUAL(gf3.getParameters(), gf2.getParameters())
END_SECTION

START_SECTION((virtual ~GaussFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual GaussFitter1D& operator=(const  GaussFitter1D &source)))
 	GaussFitter1D gf1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
	gf1.setParameters(param);

  GaussFitter1D gf2;
  gf2 = gf1;

  GaussFitter1D gf3;
	gf3.setParameters(param);

  gf1 = GaussFitter1D();
	TEST_EQUAL(gf3.getParameters(), gf2.getParameters())
END_SECTION

START_SECTION((QualityType fit1d(const  RawDataArrayType &range, InterpolationModel *&model)))
	// dummy subtest
	TEST_EQUAL(1,1)
END_SECTION

START_SECTION((Fitter1D* create()))
{
  Fitter1D* ptr = GaussFitter1D::create();
  TEST_EQUAL(ptr->getName(), "GaussFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((const String getProductName()))
{
  TEST_EQUAL(GaussFitter1D::getProductName(),"GaussFitter1D")
  TEST_EQUAL(GaussFitter1D().getName(),"GaussFitter1D")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



