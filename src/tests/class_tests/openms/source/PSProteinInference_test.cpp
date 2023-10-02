// Copyright (c) 2002-2023, The OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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

START_SECTION((LPWrapper::SOLVER getSolver()))
{
#if COINOR_SOLVER == 1
  TEST_EQUAL(ptr->getSolver(),LPWrapper::SOLVER_GLPK)
#else
  TEST_EQUAL(ptr->getSolver(), LPWrapper::SOLVER_GLPK)
#endif
}
END_SECTION

delete ptr;

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



