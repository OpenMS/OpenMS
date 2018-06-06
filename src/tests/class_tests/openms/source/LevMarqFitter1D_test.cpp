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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/InterpolationModel.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/LevMarqFitter1D.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/Fitter1D.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

///////////////////////////


class TestModel : public LevMarqFitter1D
{
  public: TestModel() : LevMarqFitter1D()
  {
    setName(getProductName());
    check_defaults_ = false;
    defaultsToParam_();
  }


  TestModel(const TestModel& source) : LevMarqFitter1D(source)
  {
    updateMembers_();
  }

  ~TestModel() override
  {
  }

  virtual TestModel& operator = (const TestModel& source)
  {
    if (&source == this) return *this;

    LevMarqFitter1D::operator = (source);
    updateMembers_();

    return *this;
  }

  void updateMembers_() override
  {
     LevMarqFitter1D::updateMembers_();
  }

  QualityType fit1d(const RawDataArrayType& /*range*/, InterpolationModel*&  /*model*/) override
  {
//    double center = 0.0;
//    center = model->getCenter();

    return 1.0;
  }

  void optimize_()
  {
  }

  static const String getProductName()
  {
    return "TestModel";
  }

};

/////////////////////////////////////////////////////////////

START_TEST(LevMarqFitter1D, "$Id$")

///////////////////////////


/////////////////////////////////////////////////////////////


TestModel* ptr = nullptr;
TestModel* nullPointer = nullptr;
START_SECTION((LevMarqFitter1D()))
	ptr = new TestModel();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((LevMarqFitter1D(const  LevMarqFitter1D &source)))
	TestModel tm1;

  TestModel tm2(tm1);
	TEST_EQUAL(tm1.getProductName(),tm2.getProductName())
END_SECTION

START_SECTION((virtual ~LevMarqFitter1D()))
	delete ptr;
END_SECTION

START_SECTION((virtual LevMarqFitter1D& operator=(const  LevMarqFitter1D &source)))
	TestModel tm1;
  TestModel tm2;

  tm2 = tm1;
	TEST_EQUAL(tm1.getProductName(),tm2.getProductName())
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



