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
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ReactionMonitoringTransition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ReactionMonitoringTransition* ptr = nullptr;
ReactionMonitoringTransition* nullPointer = nullptr;

START_SECTION(ReactionMonitoringTransition())
  ptr = new ReactionMonitoringTransition();
  TEST_NOT_EQUAL(ptr, nullPointer)
END_SECTION

START_SECTION(~ReactionMonitoringTransition())
  delete ptr;
END_SECTION

OpenMS::ReactionMonitoringTransition transition = ReactionMonitoringTransition();

CVTerm charge_cv;
String charge_cv_acc = "MS:1000041";
charge_cv.setCVIdentifierRef("MS");
charge_cv.setAccession(charge_cv_acc);
charge_cv.setName("charge state");
charge_cv.setValue(3);

START_SECTION((ReactionMonitoringTransition(const ReactionMonitoringTransition &rhs)))
{
  ReactionMonitoringTransition tr1, tr2, tr3;

  tr1.addPrecursorCVTerm(charge_cv);
  tr1.setPrecursorMZ(42.0);
	tr2 = ReactionMonitoringTransition(tr1);
  TEST_EQUAL(tr1 == tr2, true)
  ReactionMonitoringTransition::Prediction p;
  p.contact_ref = "dummy";
  tr1.setPrediction(p);
  tr1.setIdentifyingTransition(false);
  tr1.setDetectingTransition(false);
  tr1.setQuantifyingTransition(false);
	tr3 = ReactionMonitoringTransition(tr1);
  TEST_EQUAL(tr1 == tr3, true)
  TEST_EQUAL(tr1 == tr2, false)


}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setName("test_tr");

  TEST_EQUAL(tr.getName(), "test_tr")
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  TEST_EQUAL(transition.getName(), "")
}
END_SECTION

START_SECTION((void setPeptideRef(const String &peptide_ref)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setPeptideRef("test_ref");

  TEST_EQUAL(tr.getPeptideRef(), "test_ref")
}
END_SECTION

START_SECTION((const String& getPeptideRef() const ))
{
  TEST_EQUAL(transition.getPeptideRef(), "")
}
END_SECTION

START_SECTION((void setCompoundRef(const String &compound_ref)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setCompoundRef("test_ref");

  TEST_EQUAL(tr.getCompoundRef(), "test_ref")
}
END_SECTION

START_SECTION((const String& getCompoundRef() const ))
{
  TEST_EQUAL(transition.getCompoundRef(), "")
}
END_SECTION

START_SECTION((void setPrecursorMZ(double mz)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setPrecursorMZ(42.0);

  TEST_REAL_SIMILAR(tr.getPrecursorMZ(), 42.0)
}
END_SECTION

START_SECTION((double getPrecursorMZ() const ))
{
  TEST_REAL_SIMILAR(transition.getPrecursorMZ(), 0.0)
}
END_SECTION

START_SECTION((void setPrecursorCVTermList(const CVTermList &list)))
{
  CVTermList list;
  list.addCVTerm(charge_cv);
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setPrecursorCVTermList(list);

  TEST_EQUAL(tr.getPrecursorCVTermList().hasCVTerm(charge_cv_acc ), true)
}
END_SECTION

START_SECTION((bool hasPrecursorCVTerms() const))
{
  CVTermList list;
  list.addCVTerm(charge_cv);
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  TEST_EQUAL(tr.hasPrecursorCVTerms(), false)
  tr.setPrecursorCVTermList(list);
  TEST_EQUAL(tr.hasPrecursorCVTerms(), true)
}
END_SECTION

START_SECTION((void addPrecursorCVTerm(const CVTerm &cv_term)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  TEST_EQUAL(tr.hasPrecursorCVTerms(), false)
  tr.addPrecursorCVTerm(charge_cv);
  TEST_EQUAL(tr.hasPrecursorCVTerms(), true)

  TEST_EQUAL(tr.getPrecursorCVTermList().hasCVTerm(charge_cv_acc ), true)
}
END_SECTION

START_SECTION((const CVTermList& getPrecursorCVTermList() const ))
{
  CVTermList list;
  list.addCVTerm(charge_cv);
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setPrecursorCVTermList(list);

  TEST_EQUAL(tr.getPrecursorCVTermList() == list, true)
}
END_SECTION

START_SECTION((void setProductMZ(double mz)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setProductMZ(42.0);

  TEST_REAL_SIMILAR(tr.getProductMZ(), 42.0)
}
END_SECTION

START_SECTION((double getProductMZ() const ))
{
  TEST_REAL_SIMILAR(transition.getProductMZ(), 0.0)
}
END_SECTION

START_SECTION((void addProductCVTerm(const CVTerm &cv_term)))
{
  // TODO
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::isDetectingTransition() const))
{
  TEST_EQUAL(transition.isDetectingTransition(), true)
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::setDetectingTransition(bool val)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setDetectingTransition(false);
  TEST_EQUAL(tr.isDetectingTransition(), false)
  tr.setDetectingTransition(true);
  TEST_EQUAL(tr.isDetectingTransition(), true)
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::isIdentifyingTransition() const))
{
  TEST_EQUAL(transition.isIdentifyingTransition(), false)
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::setIdentifyingTransition(bool val)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setIdentifyingTransition(true);
  TEST_EQUAL(tr.isIdentifyingTransition(), true)
  tr.setIdentifyingTransition(false);
  TEST_EQUAL(tr.isIdentifyingTransition(), false)
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::isQuantifyingTransition() const))
{
  TEST_EQUAL(transition.isQuantifyingTransition(), true)
}
END_SECTION

START_SECTION((ReactionMonitoringTransition::setQuantifyingTransition(bool val)))
{
  OpenMS::ReactionMonitoringTransition tr = ReactionMonitoringTransition();
  tr.setQuantifyingTransition(false);
  TEST_EQUAL(tr.isQuantifyingTransition(), false)
  tr.setQuantifyingTransition(true);
  TEST_EQUAL(tr.isQuantifyingTransition(), true)
}
END_SECTION


START_SECTION((bool operator==(const ReactionMonitoringTransition &rhs) const ))
{
  ReactionMonitoringTransition tr1, tr2;
  TEST_EQUAL(tr1 == tr2, true)

  tr1.addPrecursorCVTerm(charge_cv);
  TEST_EQUAL(tr1 == tr2, false)
  tr2.addPrecursorCVTerm(charge_cv);
  TEST_EQUAL(tr1 == tr2, true)

  tr1.setDetectingTransition(false);
  TEST_EQUAL(tr1 == tr2, false)
  tr2.setDetectingTransition(false);
  TEST_EQUAL(tr1 == tr2, true)

}
END_SECTION

START_SECTION((bool operator!=(const ReactionMonitoringTransition &rhs) const ))
{
  ReactionMonitoringTransition tr1, tr2;
  TEST_EQUAL(tr1 != tr2, false)

  tr1.addPrecursorCVTerm(charge_cv);
  TEST_EQUAL(tr1 != tr2, true)
}
END_SECTION

START_SECTION((ReactionMonitoringTransition& operator=(const ReactionMonitoringTransition &rhs)))
{
  ReactionMonitoringTransition tr1, tr2, tr3;

  tr1.addPrecursorCVTerm(charge_cv);
  tr1.setPrecursorMZ(42.0);
	tr2 = tr1;
  TEST_EQUAL(tr1 == tr2, true)
  ReactionMonitoringTransition::Prediction p;
  p.contact_ref = "dummy";
  tr1.setPrediction(p);
  tr3 = tr1;
  TEST_EQUAL(tr1 == tr3, true)
  TEST_EQUAL(tr1 == tr2, false)

  tr1.setDetectingTransition(false);
  TEST_EQUAL(tr1 == tr3, false)
  tr1.setIdentifyingTransition(true);
  tr1.setQuantifyingTransition(false);
  TEST_EQUAL(tr1 == tr3, false)
  tr3 = tr1;
  TEST_EQUAL(tr1 == tr3, true)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



