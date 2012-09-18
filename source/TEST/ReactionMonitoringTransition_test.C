// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry               
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2012.
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
// $Maintainer: Andreas Bertsch $
// $Authors: Andreas Bertsch $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////
#include <OpenMS/ANALYSIS/MRM/ReactionMonitoringTransition.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(ReactionMonitoringTransition, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ReactionMonitoringTransition* ptr = 0;
ReactionMonitoringTransition* nullPointer = 0;
START_SECTION(ReactionMonitoringTransition())
{
	ptr = new ReactionMonitoringTransition();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~ReactionMonitoringTransition())
{
	delete ptr;
}
END_SECTION

START_SECTION((ReactionMonitoringTransition(const ReactionMonitoringTransition &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void setName(const String &name)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getName() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPeptideRef(const String &peptide_ref)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getPeptideRef() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setCompoundRef(const String &compound_ref)))
{
  // TODO
}
END_SECTION

START_SECTION((const String& getCompoundRef() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrecursorMZ(DoubleReal mz)))
{
  // TODO
}
END_SECTION

START_SECTION((DoubleReal getPrecursorMZ() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setPrecursorCVTermList(const CVTermList &list)))
{
  // TODO
}
END_SECTION

START_SECTION((void addPrecursorCVTerm(const CVTerm &cv_term)))
{
  // TODO
}
END_SECTION

START_SECTION((const CVTermList& getPrecursorCVTermList() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProductMZ(DoubleReal mz)))
{
  // TODO
}
END_SECTION

START_SECTION((DoubleReal getProductMZ() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setProductCVTermList(const CVTermList &list)))
{
  // TODO
}
END_SECTION

START_SECTION((void addProductCVTerm(const CVTerm &cv_term)))
{
  // TODO
}
END_SECTION

START_SECTION((const CVTermList& getProductCVTermList() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setInterpretations(const std::vector< TransitionInterpretation > &interpretations)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<TransitionInterpretation>& getInterpretations() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addInterpretation(const TransitionInterpretation &interpretation)))
{
  // TODO
}
END_SECTION

START_SECTION((void setConfigurations(const std::vector< Configuration > &configuration)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Configuration>& getConfigurations() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addConfiguration(const Configuration &configuration)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const ReactionMonitoringTransition &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator!=(const ReactionMonitoringTransition &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((ReactionMonitoringTransition& operator=(const ReactionMonitoringTransition &rhs)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



