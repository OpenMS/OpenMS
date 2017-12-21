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
// $Authors: Marc Sturm $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/Acquisition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(Acquisition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

Acquisition* ptr = nullptr;
Acquisition* nullPointer = nullptr;
START_SECTION(Acquisition())
	ptr = new Acquisition();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~Acquisition())
	delete ptr;
END_SECTION

START_SECTION(const String& getIdentifier() const)
  Acquisition tmp;
  TEST_EQUAL(tmp.getIdentifier(), "");
END_SECTION

START_SECTION(void setIdentifier(const String& identifier))
	Acquisition tmp;
	tmp.setIdentifier("5");
  TEST_EQUAL(tmp.getIdentifier(), "5");
END_SECTION

START_SECTION(Acquisition(const Acquisition& source))
	Acquisition tmp;
	tmp.setIdentifier("5");
	tmp.setMetaValue("label",String("label"));
	Acquisition tmp2(tmp);
	TEST_EQUAL(tmp2.getIdentifier(), "5");
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
END_SECTION

START_SECTION(Acquisition& operator= (const Acquisition& source))
	Acquisition tmp,tmp2,tmp3;
	// assignment of a modified object
	tmp2.setIdentifier("5");
	tmp2.setMetaValue("label",String("label"));
	tmp = tmp2;
	TEST_EQUAL(tmp.getIdentifier(), "5");
	TEST_EQUAL((String)(tmp.getMetaValue("label")), String("label"));
	
	// assignment of a default-constructed object
	tmp = tmp3;
	TEST_EQUAL(tmp.getIdentifier(), "");
	TEST_EQUAL(tmp.isMetaEmpty(), true);	
END_SECTION

START_SECTION(bool operator== (const Acquisition& rhs) const)
	Acquisition tmp,tmp2;
	
	TEST_EQUAL(tmp==tmp2, true);
	
	tmp2.setIdentifier("5");
	TEST_EQUAL(tmp==tmp2, false);
	
	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp==tmp2, false);
END_SECTION

START_SECTION(bool operator!= (const Acquisition& rhs) const)
	Acquisition tmp,tmp2;
	
	TEST_EQUAL(tmp!=tmp2, false);
	
	tmp2.setIdentifier("5");
	TEST_EQUAL(tmp!=tmp2, true);
	
	tmp2 = tmp;
	tmp.setMetaValue("label",String("label"));
	TEST_EQUAL(tmp!=tmp2, true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



