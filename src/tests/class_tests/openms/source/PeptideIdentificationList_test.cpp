// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/PeptideIdentificationList.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(PeptideIdentificationList, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

PeptideIdentificationList* ptr = nullptr;
PeptideIdentificationList* null_ptr = nullptr;

START_SECTION(PeptideIdentificationList())
{
  ptr = new PeptideIdentificationList();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~PeptideIdentificationList())
{
  delete ptr;
}
END_SECTION

START_SECTION((PeptideIdentificationList(const PeptideIdentificationList& vec)))
{
  PeptideIdentificationList base_vector;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  base_vector.push_back(pep_id);
  
  PeptideIdentificationList pep_ids(base_vector);
  TEST_EQUAL(pep_ids.size(), 1)
  TEST_REAL_SIMILAR(pep_ids[0].getRT(), 1234.5)
  TEST_REAL_SIMILAR(pep_ids[0].getMZ(), 567.8)
}
END_SECTION

START_SECTION((PeptideIdentificationList(PeptideIdentificationList&& vec)))
{
  PeptideIdentificationList base_vector;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  base_vector.push_back(pep_id);
  
  PeptideIdentificationList pep_ids(std::move(base_vector));
  TEST_EQUAL(pep_ids.size(), 1)
  TEST_REAL_SIMILAR(pep_ids[0].getRT(), 1234.5)
  TEST_REAL_SIMILAR(pep_ids[0].getMZ(), 567.8)
}
END_SECTION

START_SECTION((PeptideIdentificationList(size_type count)))
{
  PeptideIdentificationList pep_ids(5);
  TEST_EQUAL(pep_ids.size(), 5)
}
END_SECTION

START_SECTION((PeptideIdentificationList(size_type count, const PeptideIdentification& value)))
{
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  
  PeptideIdentificationList pep_ids(3, pep_id);
  TEST_EQUAL(pep_ids.size(), 3)
  for (size_t i = 0; i < 3; ++i)
  {
    TEST_REAL_SIMILAR(pep_ids[i].getRT(), 1234.5)
    TEST_REAL_SIMILAR(pep_ids[i].getMZ(), 567.8)
  }
}
END_SECTION

START_SECTION((PeptideIdentificationList(InputIt first, InputIt last)))
{
  PeptideIdentificationList base_vector;
  for (int i = 0; i < 3; ++i)
  {
    PeptideIdentification pep_id;
    pep_id.setRT(i * 100.0);
    pep_id.setMZ(i * 200.0);
    base_vector.push_back(pep_id);
  }
  
  PeptideIdentificationList pep_ids(base_vector.begin(), base_vector.end());
  TEST_EQUAL(pep_ids.size(), 3)
  for (size_t i = 0; i < 3; ++i)
  {
    TEST_REAL_SIMILAR(pep_ids[i].getRT(), i * 100.0)
    TEST_REAL_SIMILAR(pep_ids[i].getMZ(), i * 200.0)
  }
}
END_SECTION

START_SECTION((PeptideIdentificationList(std::initializer_list<PeptideIdentification> init)))
{
  PeptideIdentification pep_id1, pep_id2;
  pep_id1.setRT(100.0);
  pep_id1.setMZ(200.0);
  pep_id2.setRT(300.0);
  pep_id2.setMZ(400.0);
  
  PeptideIdentificationList pep_ids{pep_id1, pep_id2};
  TEST_EQUAL(pep_ids.size(), 2)
  TEST_REAL_SIMILAR(pep_ids[0].getRT(), 100.0)
  TEST_REAL_SIMILAR(pep_ids[0].getMZ(), 200.0)
  TEST_REAL_SIMILAR(pep_ids[1].getRT(), 300.0)
  TEST_REAL_SIMILAR(pep_ids[1].getMZ(), 400.0)
}
END_SECTION

START_SECTION((PeptideIdentificationList(const PeptideIdentificationList&)))
{
  PeptideIdentificationList pep_ids1;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  pep_ids1.push_back(pep_id);
  
  PeptideIdentificationList pep_ids2(pep_ids1);
  TEST_EQUAL(pep_ids2.size(), 1)
  TEST_REAL_SIMILAR(pep_ids2[0].getRT(), 1234.5)
  TEST_REAL_SIMILAR(pep_ids2[0].getMZ(), 567.8)
}
END_SECTION

START_SECTION((PeptideIdentificationList(PeptideIdentificationList&&)))
{
  PeptideIdentificationList pep_ids1;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  pep_ids1.push_back(pep_id);
  
  PeptideIdentificationList pep_ids2(std::move(pep_ids1));
  TEST_EQUAL(pep_ids2.size(), 1)
  TEST_REAL_SIMILAR(pep_ids2[0].getRT(), 1234.5)
  TEST_REAL_SIMILAR(pep_ids2[0].getMZ(), 567.8)
}
END_SECTION

START_SECTION((PeptideIdentificationList& operator=(const PeptideIdentificationList&)))
{
  PeptideIdentificationList pep_ids1;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  pep_ids1.push_back(pep_id);
  
  PeptideIdentificationList pep_ids2;
  pep_ids2 = pep_ids1;
  TEST_EQUAL(pep_ids2.size(), 1)
  TEST_REAL_SIMILAR(pep_ids2[0].getRT(), 1234.5)
  TEST_REAL_SIMILAR(pep_ids2[0].getMZ(), 567.8)
}
END_SECTION

START_SECTION((PeptideIdentificationList& operator=(PeptideIdentificationList&&)))
{
  PeptideIdentificationList pep_ids1;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  pep_ids1.push_back(pep_id);
  
  PeptideIdentificationList pep_ids2;
  pep_ids2 = std::move(pep_ids1);
  TEST_EQUAL(pep_ids2.size(), 1)
  TEST_REAL_SIMILAR(pep_ids2[0].getRT(), 1234.5)
  TEST_REAL_SIMILAR(pep_ids2[0].getMZ(), 567.8)
}
END_SECTION

START_SECTION((PeptideIdentificationList& operator=(const PeptideIdentificationList& vec)))
{
  PeptideIdentificationList base_vector;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  base_vector.push_back(pep_id);
  
  PeptideIdentificationList pep_ids;
  pep_ids = base_vector;
  TEST_EQUAL(pep_ids.size(), 1)
  TEST_REAL_SIMILAR(pep_ids[0].getRT(), 1234.5)
  TEST_REAL_SIMILAR(pep_ids[0].getMZ(), 567.8)
}
END_SECTION

START_SECTION((PeptideIdentificationList& operator=(PeptideIdentificationList&& vec)))
{
  PeptideIdentificationList base_vector;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  base_vector.push_back(pep_id);
  
  PeptideIdentificationList pep_ids;
  pep_ids = std::move(base_vector);
  TEST_EQUAL(pep_ids.size(), 1)
  TEST_REAL_SIMILAR(pep_ids[0].getRT(), 1234.5)
  TEST_REAL_SIMILAR(pep_ids[0].getMZ(), 567.8)
}
END_SECTION

START_SECTION((PeptideIdentificationList& operator=(std::initializer_list<PeptideIdentification> init)))
{
  PeptideIdentification pep_id1, pep_id2;
  pep_id1.setRT(100.0);
  pep_id1.setMZ(200.0);
  pep_id2.setRT(300.0);
  pep_id2.setMZ(400.0);
  
  PeptideIdentificationList pep_ids;
  pep_ids = {pep_id1, pep_id2};
  TEST_EQUAL(pep_ids.size(), 2)
  TEST_REAL_SIMILAR(pep_ids[0].getRT(), 100.0)
  TEST_REAL_SIMILAR(pep_ids[0].getMZ(), 200.0)
  TEST_REAL_SIMILAR(pep_ids[1].getRT(), 300.0)
  TEST_REAL_SIMILAR(pep_ids[1].getMZ(), 400.0)
}
END_SECTION

START_SECTION((Test vector functionality))
{
  // Test that it behaves like a vector
  PeptideIdentificationList pep_ids;
  
  // Test push_back
  PeptideIdentification pep_id1;
  pep_id1.setRT(100.0);
  pep_ids.push_back(pep_id1);
  TEST_EQUAL(pep_ids.size(), 1)
  
  // Test clear
  pep_ids.clear();
  TEST_EQUAL(pep_ids.size(), 0)
  TEST_EQUAL(pep_ids.empty(), true)
  
  // Test resize
  pep_ids.resize(5);
  TEST_EQUAL(pep_ids.size(), 5)
}
END_SECTION

START_SECTION((Test comparison operators))
{
  PeptideIdentificationList pep_ids1, pep_ids2;
  
  // Test equality of empty lists
  TEST_EQUAL(pep_ids1 == pep_ids2, true)
  TEST_EQUAL(pep_ids1 != pep_ids2, false)
  
  // Add element to one list
  PeptideIdentification pep_id;
  pep_id.setRT(100.0);
  pep_ids1.push_back(pep_id);
  
  // Test inequality
  TEST_EQUAL(pep_ids1 == pep_ids2, false)
  TEST_EQUAL(pep_ids1 != pep_ids2, true)
  
  // Add same element to second list
  pep_ids2.push_back(pep_id);
  
  // Test equality again
  TEST_EQUAL(pep_ids1 == pep_ids2, true)
  TEST_EQUAL(pep_ids1 != pep_ids2, false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST