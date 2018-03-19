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
#include <OpenMS/METADATA/MetaInfoDescription.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MetaInfoDescription, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

MetaInfoDescription* ptr = nullptr;
MetaInfoDescription* nullPointer = nullptr;
START_SECTION((MetaInfoDescription()))
	ptr = new MetaInfoDescription();
	TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION((~MetaInfoDescription()))
	delete ptr;
END_SECTION

START_SECTION((const String& getName() const))
  MetaInfoDescription tmp;
  TEST_EQUAL(tmp.getName(),"");
END_SECTION

START_SECTION((void setName(const String& name)))
  MetaInfoDescription tmp;
  tmp.setName("name");
  TEST_EQUAL(tmp.getName(),"name");
END_SECTION

START_SECTION((const std::vector<DataProcessing>& getDataProcessing() const))
  MetaInfoDescription tmp;
  TEST_EQUAL(tmp.getDataProcessing().size(),0);
END_SECTION

START_SECTION((void setDataProcessing(const std::vector< DataProcessing > &data_processing)))
  MetaInfoDescription tmp;
  std::vector<DataProcessingPtr> dummy;
  dummy.resize(1);
  tmp.setDataProcessing(dummy);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION((std::vector<DataProcessing>& getDataProcessing()))
  MetaInfoDescription tmp;
  tmp.getDataProcessing().resize(1);
  TEST_EQUAL(tmp.getDataProcessing().size(),1);
END_SECTION

START_SECTION((MetaInfoDescription(const MetaInfoDescription& source)))
  MetaInfoDescription tmp;
  tmp.setName("bla2");
  tmp.getDataProcessing().resize(1);
  tmp.setMetaValue("label",String("label"));
  
  MetaInfoDescription tmp2(tmp);
  TEST_EQUAL(tmp2.getName(),"bla2");
	TEST_EQUAL(tmp.getDataProcessing().size(),1);
  TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
END_SECTION

START_SECTION((MetaInfoDescription& operator= (const MetaInfoDescription& source)))
  MetaInfoDescription tmp;
  tmp.setName("bla2");
  tmp.getDataProcessing().resize(1);
  tmp.setMetaValue("label",String("label"));
  
  MetaInfoDescription tmp2;
  tmp2 = tmp;
  TEST_EQUAL(tmp2.getName(),"bla2");
	TEST_EQUAL(tmp.getDataProcessing().size(),1);
	TEST_EQUAL((String)(tmp2.getMetaValue("label")), "label");
	
  tmp2 = MetaInfoDescription();
  TEST_EQUAL(tmp2.getName(),"");
	TEST_EQUAL(tmp2.getDataProcessing().size(),0);
  TEST_EQUAL(tmp2.getMetaValue("label").isEmpty(), true);
END_SECTION

START_SECTION((bool operator== (const MetaInfoDescription& rhs) const))
  MetaInfoDescription edit, empty;
  
  TEST_EQUAL(edit==empty, true);
  
  edit = empty;
  edit.setName("bla2");
	TEST_EQUAL(edit==empty, false);

  edit = empty;
  edit.setMetaValue("label",String("label"));
  TEST_EQUAL(edit==empty, false);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



