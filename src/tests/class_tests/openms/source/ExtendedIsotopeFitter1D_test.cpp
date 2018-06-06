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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ExtendedIsotopeFitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ExtendedIsotopeFitter1D, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ExtendedIsotopeFitter1D* ptr = nullptr;
ExtendedIsotopeFitter1D* nullPointer = nullptr;
START_SECTION(ExtendedIsotopeFitter1D())
{
	ptr = new ExtendedIsotopeFitter1D();
        TEST_EQUAL(ptr->getName(), "ExtendedIsotopeFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((ExtendedIsotopeFitter1D(const  ExtendedIsotopeFitter1D &source)))
	ExtendedIsotopeFitter1D eisof1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "charge", 1 );
  param.setValue( "isotope:stdev", 0.04 );
  param.setValue( "isotope:maximum", 20 );                	
	eisof1.setParameters(param);

	ExtendedIsotopeFitter1D eisof2(eisof1);
  ExtendedIsotopeFitter1D eisof3;
	eisof3.setParameters(param);
  eisof1 = ExtendedIsotopeFitter1D();
	TEST_EQUAL(eisof3.getParameters(), eisof2.getParameters())

END_SECTION

START_SECTION((virtual ~ExtendedIsotopeFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual ExtendedIsotopeFitter1D& operator=(const  ExtendedIsotopeFitter1D &source)))
	ExtendedIsotopeFitter1D eisof1;
	
	Param param;
	param.setValue( "tolerance_stdev_bounding_box", 1.0);
  param.setValue( "statistics:mean", 680.1 );
  param.setValue( "statistics:variance", 2.0 );
  param.setValue( "interpolation_step", 1.0 );
  param.setValue( "charge", 1 );
  param.setValue( "isotope:stdev", 0.04 );
  param.setValue( "isotope:maximum", 20 );                	
	eisof1.setParameters(param);

  ExtendedIsotopeFitter1D eisof2;
  eisof2 = eisof1;

  ExtendedIsotopeFitter1D eisof3;
	eisof3.setParameters(param);

  eisof1 = ExtendedIsotopeFitter1D();
	TEST_EQUAL(eisof3.getParameters(), eisof3.getParameters())
END_SECTION

START_SECTION((QualityType fit1d(const  RawDataArrayType &range, InterpolationModel *&model)))
	// dummy subtest TODO
	TEST_EQUAL(1,1)
END_SECTION

START_SECTION((Fitter1D* create()))
  Fitter1D* ptr = ExtendedIsotopeFitter1D::create();
  TEST_EQUAL(ptr->getName(), "ExtendedIsotopeFitter1D")
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((const String getProductName()))
  TEST_EQUAL(ExtendedIsotopeFitter1D::getProductName(),"ExtendedIsotopeFitter1D")
  TEST_EQUAL(ExtendedIsotopeFitter1D().getName(),"ExtendedIsotopeFitter1D")
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



