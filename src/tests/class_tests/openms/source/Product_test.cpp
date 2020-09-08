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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/Product.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Product, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Product* ptr = nullptr;
Product* nullPointer = nullptr;
START_SECTION((Product()))
	ptr = new Product();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~Product()))
	delete ptr;
END_SECTION

START_SECTION((double getMZ() const ))
  Product tmp;
  TEST_EQUAL(tmp.getMZ(),0);
END_SECTION

START_SECTION((void setMZ(double mz)))
  Product tmp;
  tmp.setMZ(47.11);
  TEST_REAL_SIMILAR(tmp.getMZ(),47.11);
END_SECTION

START_SECTION((double getIsolationWindowUpperOffset() const ))
  Product tmp;
  TEST_REAL_SIMILAR(tmp.getIsolationWindowUpperOffset(), 0);
END_SECTION

START_SECTION((void setIsolationWindowUpperOffset(double bound)))
  Product tmp;
  tmp.setIsolationWindowUpperOffset(22.7);
  TEST_REAL_SIMILAR(tmp.getIsolationWindowUpperOffset(), 22.7);
END_SECTION

START_SECTION((double getIsolationWindowLowerOffset() const ))
  Product tmp;
  TEST_REAL_SIMILAR(tmp.getIsolationWindowLowerOffset(), 0);
END_SECTION

START_SECTION((void setIsolationWindowLowerOffset(double bound)))
  Product tmp;
  tmp.setIsolationWindowLowerOffset(22.8);
  TEST_REAL_SIMILAR(tmp.getIsolationWindowLowerOffset(), 22.8);
END_SECTION

START_SECTION((Product(const Product& source)))
	Product tmp;
	tmp.setMZ(47.11);
  tmp.setIsolationWindowUpperOffset(22.7);
  tmp.setIsolationWindowLowerOffset(22.8);
	tmp.setMetaValue("label",String("label"));
	
	Product tmp2(tmp);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 22.7);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 22.8);
	TEST_REAL_SIMILAR(tmp2.getMZ(),47.11);
END_SECTION

START_SECTION((Product& operator= (const Product& source)))
	Product tmp;
	tmp.setMZ(47.11);
  tmp.setIsolationWindowUpperOffset(22.7);
  tmp.setIsolationWindowLowerOffset(22.8);
	tmp.setMetaValue("label",String("label"));
	
	//normal assignment
	Product tmp2;
	tmp2 = tmp;
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 22.7);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 22.8);
	TEST_REAL_SIMILAR(tmp2.getMZ(),47.11);
		
	//assignment of empty object
	tmp2 = Product();
	TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowUpperOffset(), 0.0);
	TEST_REAL_SIMILAR(tmp2.getIsolationWindowLowerOffset(), 0.0);
	TEST_REAL_SIMILAR(tmp2.getMZ(),0.0);
END_SECTION

START_SECTION((bool operator== (const Product& rhs) const))
	Product tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp2.setMZ(47.11);
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
  tmp.setIsolationWindowUpperOffset(22.7);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
  tmp.setIsolationWindowLowerOffset(22.8);
	TEST_EQUAL(tmp==tmp2, false);

	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp==tmp2, false);
END_SECTION

START_SECTION((bool operator!= (const Product& rhs) const))
	Product tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp2.setMZ(47.11);
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;	tmp2 = tmp;
  tmp.setIsolationWindowUpperOffset(22.7);
	TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;	tmp2 = tmp;
  tmp.setIsolationWindowLowerOffset(22.8);
	TEST_EQUAL(tmp!=tmp2, true);

	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp!=tmp2, true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



