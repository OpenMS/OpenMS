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
// $Authors: Alexandra Zerck$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/TARGETED/PSProteinInference.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/FORMAT/IdXMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PSProteinInference, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PSProteinInference* ptr = nullptr;
PSProteinInference* null_ptr = nullptr;
START_SECTION(PSProteinInference())
{
	ptr = new PSProteinInference();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION((virtual ~PSProteinInference()))
{
	delete ptr;
}
END_SECTION

vector<PeptideIdentification> pep_ids;
vector<ProteinIdentification> prot_ids;
String document_id;
IdXMLFile().load(OPENMS_GET_TEST_DATA_PATH("PSProteinInference_test_input.iDXML"),prot_ids,pep_ids,document_id);
ptr = new PSProteinInference();
START_SECTION((Size findMinimalProteinList(const std::vector< PeptideIdentification > &peptide_ids)))
{
  Size num = ptr->findMinimalProteinList(pep_ids);
  TEST_EQUAL(num,3)
}
END_SECTION

START_SECTION((void calculateProteinProbabilities(const std::vector< PeptideIdentification > &ids)))
{
  ptr->calculateProteinProbabilities(pep_ids);
  TEST_REAL_SIMILAR(ptr->getProteinProbability("A2RUR9"),1.)
  TEST_REAL_SIMILAR(ptr->getProteinProbability("O00762"),0.132644)
  TEST_REAL_SIMILAR(ptr->getProteinProbability("O14795"),0.99999)
  TEST_REAL_SIMILAR(ptr->getProteinProbability("O15085"),1.)
}
END_SECTION

START_SECTION((double getProteinProbability(const String &acc)))
{
  TEST_REAL_SIMILAR(ptr->getProteinProbability("A2RUR9"),1.)
}
END_SECTION

START_SECTION((bool isProteinInMinimalList(const String &acc)))
{
  TEST_EQUAL(ptr->isProteinInMinimalList("A2RUR9"),true)
  TEST_EQUAL(ptr->isProteinInMinimalList("O00762"),true)
  TEST_EQUAL(ptr->isProteinInMinimalList("O14795"),false)
  TEST_EQUAL(ptr->isProteinInMinimalList("O15085"),true)
}
END_SECTION

START_SECTION((Int getNumberOfProtIds(double protein_id_threshold)))
{
  TEST_EQUAL(ptr->getNumberOfProtIds(.95),2)
}
END_SECTION

START_SECTION((Int getNumberOfProtIdsPeptideRule(Int min_peptides, std::map< String, std::set< String > > &prot_id_counter)))
{
  std::map< String, std::set< String > > prot_id_counter;
  set< String > peps;
  peps.insert("PEPTIDEONE");
  prot_id_counter.insert(make_pair("A2RUR9",peps));
  peps.insert("PEPTIDETWO");
  prot_id_counter.insert(make_pair("O00762",peps));
  peps.insert("PEPTIDETHREE");
  prot_id_counter.insert(make_pair("O15085",peps));
  TEST_EQUAL(ptr->getNumberOfProtIdsPeptideRule(3,prot_id_counter),1)
}
END_SECTION

START_SECTION((void setSolver(LPWrapper::SOLVER solver)))
{
  ptr->setSolver(LPWrapper::SOLVER_GLPK);
  TEST_EQUAL(ptr->getSolver(),LPWrapper::SOLVER_GLPK)
#if COINOR_SOLVER == 1
  ptr->setSolver(LPWrapper::SOLVER_COINOR);
  TEST_EQUAL(ptr->getSolver(),LPWrapper::SOLVER_COINOR)
#endif
}
END_SECTION

START_SECTION((LPWrapper::SOLVER getSolver()))
{
  ptr->setSolver(LPWrapper::SOLVER_GLPK);
  TEST_EQUAL(ptr->getSolver(),LPWrapper::SOLVER_GLPK)
#if COINOR_SOLVER == 1
  ptr->setSolver(LPWrapper::SOLVER_COINOR);
  TEST_EQUAL(ptr->getSolver(),LPWrapper::SOLVER_COINOR)
#endif
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



