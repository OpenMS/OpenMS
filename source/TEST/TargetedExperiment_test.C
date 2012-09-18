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
#include <OpenMS/ANALYSIS/TARGETED/TargetedExperiment.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(TargetedExperiment, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

TargetedExperiment* ptr = 0;
TargetedExperiment* nullPointer = 0;
START_SECTION(TargetedExperiment())
{
	ptr = new TargetedExperiment();
	TEST_NOT_EQUAL(ptr, nullPointer)
}
END_SECTION

START_SECTION(virtual ~TargetedExperiment())
{
	delete ptr;
}
END_SECTION

START_SECTION((TargetedExperiment(const TargetedExperiment &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((bool operator==(const TargetedExperiment &rhs) const ))
{
  // TODO
}
END_SECTION

START_SECTION((void setCVs(const std::vector< CV > &cvs)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<CV>& getCVs() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addCV(const CV &cv)))
{
  // TODO
}
END_SECTION

START_SECTION((void setContacts(const std::vector< CVTermList > &contacts)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<CVTermList>& getContacts() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addContact(const CVTermList &contact)))
{
  // TODO
}
END_SECTION

START_SECTION((void setPublications(const std::vector< CVTermList > &publications)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<CVTermList>& getPublications() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addPublication(const CVTermList &publication)))
{
  // TODO
}
END_SECTION

START_SECTION((void setInstruments(const std::vector< CVTermList > &instruments)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<CVTermList>& getInstruments() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addInstrument(const CVTermList &instrument)))
{
  // TODO
}
END_SECTION

START_SECTION((void setSoftware(const std::vector< Software > &software)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Software>& getSoftware() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addSoftware(const Software &software)))
{
  // TODO
}
END_SECTION

START_SECTION((void setProteins(const std::vector< Protein > &proteins)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Protein>& getProteins() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addProtein(const Protein &protein)))
{
  // TODO
}
END_SECTION

START_SECTION((void setCompounds(const std::vector< Compound > &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Compound>& getCompounds() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addCompound(const Compound &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void setPeptides(const std::vector< Peptide > &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<Peptide>& getPeptides() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addPeptide(const Peptide &rhs)))
{
  // TODO
}
END_SECTION

START_SECTION((void setTransitions(const std::vector< ReactionMonitoringTransition > &transitions)))
{
  // TODO
}
END_SECTION

START_SECTION((const std::vector<ReactionMonitoringTransition>& getTransitions() const ))
{
  // TODO
}
END_SECTION

START_SECTION((void addTransition(const ReactionMonitoringTransition &transition)))
{
  // TODO
}
END_SECTION

START_SECTION((TargetedExperiment& operator=(const TargetedExperiment &rhs)))
{
  // TODO
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



