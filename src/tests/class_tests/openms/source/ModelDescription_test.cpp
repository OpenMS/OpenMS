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

#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/ModelDescription.h>
#include <OpenMS/TRANSFORMATIONS/FEATUREFINDER/IsotopeModel.h>
#include <OpenMS/CONCEPT/Exception.h>


///////////////////////////

START_TEST(ModelDescription<2>, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using std::stringstream;

ModelDescription<2>* ptr = nullptr;
ModelDescription<2>* nullPointer = nullptr;
START_SECTION((ModelDescription()))
	ptr = new ModelDescription<2>();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((virtual ~ModelDescription()))
	delete ptr;
END_SECTION

START_SECTION(BaseModel<D>* createModel())
	BaseModel<2>* ptr = ModelDescription<2>().createModel();
  BaseModel<2>* baseModel_nullPointer = nullptr;
  TEST_EQUAL(ptr, baseModel_nullPointer)	// no name is set, should be zero pointer
END_SECTION

START_SECTION( virtual bool operator==(const ModelDescription &rhs) const )
	ModelDescription<2> fp1,fp2;	
	TEST_EQUAL(fp1==fp2,true)
  
  fp1.setName("halligalli2000");
	TEST_EQUAL(fp1==fp2,false)
  
  fp2.setName("halligalli2000");
  TEST_EQUAL(fp1==fp2,true)
  
  Param param;
  param.setValue("bla","bluff");
	fp1.setParam(param);
	TEST_EQUAL(fp1==fp2,false)
  
	fp2.setParam(param);
  TEST_EQUAL(fp1==fp2,true)
END_SECTION

START_SECTION( virtual bool operator!=(const ModelDescription &rhs) const )
	ModelDescription<2> fp1, fp2;
	TEST_EQUAL(fp1!=fp2,false)
	
  fp1.setName("halligalli2000");
	TEST_EQUAL(fp1!=fp2,true)

  fp2.setName("halligalli2000");
	TEST_EQUAL(fp1!=fp2,false)
END_SECTION

START_SECTION((virtual ModelDescription& operator=(const ModelDescription &source)))
	ModelDescription<2> tm1;
  tm1.setName("halligalli");
	Param param;
	param.setValue("test","test");
	tm1.setParam(param);
	
  ModelDescription<2> tm2;
  tm2 = tm1;

	TEST_EQUAL(tm1==tm2,true)
END_SECTION

START_SECTION((ModelDescription(const ModelDescription &source)))
	ModelDescription<2> tm1;
  tm1.setName("halligalli");
	Param param;
	param.setValue("test","test");
	tm1.setParam(param);
	
  ModelDescription<2> tm2(tm1);

	TEST_EQUAL(tm1==tm2,true)
END_SECTION

START_SECTION( ModelDescription(const BaseModel< D > *model) )
	const BaseModel<1> * bm = new IsotopeModel();

  ModelDescription<1> md(bm);
	
	BaseModel<1>* ptr = md.createModel();
	TEST_EQUAL( *ptr == *bm, true)	
END_SECTION

START_SECTION((const String& getName() const ))
  const ModelDescription<2> m;
  TEST_EQUAL(m.getName(), "")
END_SECTION

START_SECTION((void setName(const String &name)))
  ModelDescription<2> m;
	m.setName("halligalli2006");
  TEST_EQUAL(m.getName(), "halligalli2006")
END_SECTION

START_SECTION((const Param& getParam() const ))
  ModelDescription<2> m;
	Param p;
	p.setValue("x1",1.0);
	p.setValue("x2",2.0);
	m.setParam(p);
  TEST_EQUAL(m.getParam(), p)
END_SECTION

START_SECTION( String& getName() )
  ModelDescription<2> m;
	m.setName("halligalli2006");
  TEST_EQUAL(m.getName(), "halligalli2006")
END_SECTION

START_SECTION( Param& getParam() )
	ModelDescription<2> m;
	Param p;
	p.setValue("x1",1.0);
	p.setValue("x2",2.0);
	m.setParam(p);
  TEST_EQUAL(m.getParam(), p)
END_SECTION

START_SECTION( void setParam(const Param &param) )
	ModelDescription<2> m;
	Param p;
	p.setValue("x1",1.0);
	p.setValue("x2",2.0);
	m.setParam(p);
  TEST_EQUAL(m.getParam(), p)
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
