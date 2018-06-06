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
// $Authors: Stephan Aiche$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/SIMULATION/EGHModel.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(EGHModel, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

EGHModel* ptr = nullptr;
EGHModel* nullPointer = nullptr;
START_SECTION(EGHModel())
{
	ptr = new EGHModel();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((EGHModel(const EGHModel &source)))
{
  EGHModel egh1;
  egh1.setInterpolationStep(0.2);

  Param tmp;
  tmp.setValue("statistics:mean", 680.1 );
  tmp.setValue("statistics:variance",  2.0);
  tmp.setValue("egh:retention", 680.1);
  tmp.setValue("egh:height",100000.0);
  tmp.setValue("egh:A",150.0);
  tmp.setValue("egh:B",100.0);
  tmp.setValue("egh:alpha",0.4);
  egh1.setParameters(tmp);

  EGHModel egh2(egh1);
  EGHModel egh3;
  egh3.setInterpolationStep(0.2);
  egh3.setParameters(tmp);

  egh1 = EGHModel();
  TEST_EQUAL(egh3.getParameters(), egh2.getParameters())
  TEST_EQUAL(egh3==egh2,true)
}
END_SECTION

START_SECTION((virtual ~EGHModel()))
{
  delete ptr;
}
END_SECTION

START_SECTION((virtual EGHModel& operator=(const EGHModel &source)))
{
  EGHModel egh1;
  egh1.setInterpolationStep(0.2);

  Param tmp;
  tmp.setValue("statistics:mean", 680.1 );
  tmp.setValue("statistics:variance",  2.0);
  tmp.setValue("egh:retention", 680.1);
  tmp.setValue("egh:height",100000.0);
  tmp.setValue("egh:A",150.0);
  tmp.setValue("egh:B",100.0);
  tmp.setValue("egh:alpha",0.4);
  egh1.setParameters(tmp);

  EGHModel egh2;
  egh2 = egh1;

  EGHModel egh3;
  egh3.setInterpolationStep(0.2);
  egh3.setParameters(tmp);

  egh1 = EGHModel();
  TEST_EQUAL(egh3.getParameters(), egh2.getParameters())
  TEST_EQUAL(egh3==egh2,true)
}
END_SECTION

START_SECTION((void setOffset(CoordinateType offset)))
{
  //
  EGHModel egh1;
  egh1.setInterpolationStep(0.2);

  Param tmp;
  tmp.setValue("statistics:mean", 680.1 );
  tmp.setValue("statistics:variance",  2.0);
  tmp.setValue("egh:retention", 680.1);
  tmp.setValue("egh:height",100000.0);
  tmp.setValue("egh:A",150.0);
  tmp.setValue("egh:B",100.0);
  tmp.setValue("egh:alpha",0.4);
  egh1.setParameters(tmp);

  //
  EGHModel::CoordinateType current_offset = egh1.getInterpolation().getOffset();
  EGHModel::CoordinateType current_mean = egh1.getCenter();
  EGHModel::CoordinateType new_offset = current_offset + 10.0;
  egh1.setOffset(new_offset);

  TEST_REAL_SIMILAR(egh1.getInterpolation().getOffset(), new_offset)
  TEST_REAL_SIMILAR(egh1.getCenter(), current_mean + 10.0)
}
END_SECTION

START_SECTION((void setSamples()))
{
  EGHModel egh1;

  Param tmp;
  tmp.setValue("statistics:mean", 1000.0 );
  tmp.setValue("statistics:variance",  2.0);
  tmp.setValue("egh:retention", 1000.0);
  tmp.setValue("egh:height",100.0);
  tmp.setValue("egh:A",10.0);
  tmp.setValue("egh:B",20.0);
  tmp.setValue("egh:alpha",0.5);
  egh1.setInterpolationStep(0.2);
  egh1.setParameters(tmp); // setSamples() is called here

  TEST_REAL_SIMILAR(egh1.getInterpolation().value(1000.0), 100.0)
  TEST_REAL_SIMILAR(egh1.getInterpolation().value(990.0), 50.0) // corresponds to A_
  TEST_REAL_SIMILAR(egh1.getInterpolation().value(1020.0), 50.0) // corresponds to B_
}
END_SECTION

START_SECTION((CoordinateType getCenter() const ))
{
  //
  EGHModel egh1;
  egh1.setInterpolationStep(0.2);

  Param tmp;
  tmp.setValue("statistics:mean", 680.1 );
  tmp.setValue("statistics:variance",  2.0);
  tmp.setValue("egh:retention", 680.1);
  tmp.setValue("egh:height",100000.0);
  tmp.setValue("egh:A",150.0);
  tmp.setValue("egh:B",100.0);
  tmp.setValue("egh:alpha",0.4);
  egh1.setParameters(tmp);

  //
  EGHModel::CoordinateType current_offset = egh1.getInterpolation().getOffset();
  EGHModel::CoordinateType current_mean = egh1.getCenter();
  EGHModel::CoordinateType new_offset = current_offset + 10.0;
  egh1.setOffset(new_offset);

  TEST_REAL_SIMILAR(egh1.getInterpolation().getOffset(), new_offset)
  TEST_REAL_SIMILAR(egh1.getCenter(), current_mean + 10.0)
}
END_SECTION

START_SECTION((static BaseModel<1>* create()))
{
  BaseModel<1>* ptr = EGHModel::create();
  TEST_EQUAL(ptr->getName(), "EGHModel")
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION((static const String getProductName()))
{
  TEST_EQUAL(EGHModel::getProductName(),"EGHModel")
  TEST_EQUAL(EGHModel().getName(),"EGHModel")
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



