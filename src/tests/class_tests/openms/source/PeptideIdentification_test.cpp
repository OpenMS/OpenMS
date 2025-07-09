// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/PeptideIdentification.h>
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

START_SECTION((PeptideIdentificationList(const Base& base)))
{
  PeptideIdentificationList base_vector;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  base_vector.push_back(pep_id);
  
  PeptideIdentificationList spectra_ids(base_vector);
  TEST_EQUAL(spectra_ids.size(), 1)
  TEST_REAL_SIMILAR(spectra_ids[0].getRT(), 1234.5)
  TEST_REAL_SIMILAR(spectra_ids[0].getMZ(), 567.8)
}
END_SECTION

START_SECTION((PeptideIdentificationList(Base&& base)))
{
  PeptideIdentificationList base_vector;
  PeptideIdentification pep_id;
  pep_id.setRT(9876.5);
  pep_id.setMZ(432.1);
  base_vector.push_back(pep_id);
  
  PeptideIdentificationList spectra_ids(std::move(base_vector));
  TEST_EQUAL(spectra_ids.size(), 1)
  TEST_REAL_SIMILAR(spectra_ids[0].getRT(), 9876.5)
  TEST_REAL_SIMILAR(spectra_ids[0].getMZ(), 432.1)
}
END_SECTION

START_SECTION((PeptideIdentificationList(size_type count)))
{
  PeptideIdentificationList spectra_ids(3);
  TEST_EQUAL(spectra_ids.size(), 3)
}
END_SECTION

START_SECTION((PeptideIdentificationList(size_type count, const PeptideIdentification& value)))
{
  PeptideIdentification pep_id;
  pep_id.setRT(555.5);
  pep_id.setMZ(777.7);
  
  PeptideIdentificationList spectra_ids(2, pep_id);
  TEST_EQUAL(spectra_ids.size(), 2)
  TEST_REAL_SIMILAR(spectra_ids[0].getRT(), 555.5)
  TEST_REAL_SIMILAR(spectra_ids[1].getRT(), 555.5)
  TEST_REAL_SIMILAR(spectra_ids[0].getMZ(), 777.7)
  TEST_REAL_SIMILAR(spectra_ids[1].getMZ(), 777.7)
}
END_SECTION

START_SECTION((PeptideIdentificationList(std::initializer_list<PeptideIdentification> init)))
{
  PeptideIdentification pep_id1, pep_id2;
  pep_id1.setRT(111.1);
  pep_id1.setMZ(222.2);
  pep_id2.setRT(333.3);
  pep_id2.setMZ(444.4);
  
  PeptideIdentificationList spectra_ids{pep_id1, pep_id2};
  TEST_EQUAL(spectra_ids.size(), 2)
  TEST_REAL_SIMILAR(spectra_ids[0].getRT(), 111.1)
  TEST_REAL_SIMILAR(spectra_ids[1].getRT(), 333.3)
}
END_SECTION

START_SECTION((PeptideIdentificationList& operator=(const Base& base)))
{
  PeptideIdentificationList base_vector;
  PeptideIdentification pep_id;
  pep_id.setRT(654.3);
  pep_id.setMZ(987.6);
  base_vector.push_back(pep_id);
  
  PeptideIdentificationList spectra_ids;
  spectra_ids = base_vector;
  TEST_EQUAL(spectra_ids.size(), 1)
  TEST_REAL_SIMILAR(spectra_ids[0].getRT(), 654.3)
  TEST_REAL_SIMILAR(spectra_ids[0].getMZ(), 987.6)
}
END_SECTION

START_SECTION((PeptideIdentificationList& operator=(Base&& base)))
{
  PeptideIdentificationList base_vector;
  PeptideIdentification pep_id;
  pep_id.setRT(147.2);
  pep_id.setMZ(258.3);
  base_vector.push_back(pep_id);
  
  PeptideIdentificationList spectra_ids;
  spectra_ids = std::move(base_vector);
  TEST_EQUAL(spectra_ids.size(), 1)
  TEST_REAL_SIMILAR(spectra_ids[0].getRT(), 147.2)
  TEST_REAL_SIMILAR(spectra_ids[0].getMZ(), 258.3)
}
END_SECTION

START_SECTION((PeptideIdentificationList& operator=(std::initializer_list<PeptideIdentification> init)))
{
  PeptideIdentification pep_id1, pep_id2;
  pep_id1.setRT(369.7);
  pep_id1.setMZ(741.8);
  pep_id2.setRT(852.9);
  pep_id2.setMZ(963.1);
  
  PeptideIdentificationList spectra_ids;
  spectra_ids = {pep_id1, pep_id2};
  TEST_EQUAL(spectra_ids.size(), 2)
  TEST_REAL_SIMILAR(spectra_ids[0].getRT(), 369.7)
  TEST_REAL_SIMILAR(spectra_ids[1].getRT(), 852.9)
}
END_SECTION

START_SECTION((Vector inheritance test))
{
  PeptideIdentificationList spectra_ids;
  
  // Test that all vector methods are available
  PeptideIdentification pep_id;
  pep_id.setRT(100.0);
  pep_id.setMZ(200.0);
  
  // Test push_back
  spectra_ids.push_back(pep_id);
  TEST_EQUAL(spectra_ids.size(), 1)
  
  // Test clear
  spectra_ids.clear();
  TEST_EQUAL(spectra_ids.size(), 0)
  TEST_EQUAL(spectra_ids.empty(), true)
  
  // Test reserve and capacity
  spectra_ids.reserve(10);
  TEST_EQUAL(spectra_ids.capacity() >= 10, true)
  
  // Test iterators
  spectra_ids.push_back(pep_id);
  spectra_ids.push_back(pep_id);
  
  size_t count = 0;
  for (auto it = spectra_ids.begin(); it != spectra_ids.end(); ++it)
  {
    count++;
  }
  TEST_EQUAL(count, 2)
  
  // Test range-based for loop
  count = 0;
  for (const auto& id : spectra_ids)
  {
    count++;
  }
  TEST_EQUAL(count, 2)
}
END_SECTION

START_SECTION((Type compatibility test))
{
  // Test that PeptideIdentificationList can be used where PeptideIdentificationList is expected
  PeptideIdentificationList spectra_ids;
  PeptideIdentification pep_id;
  pep_id.setRT(500.0);
  pep_id.setMZ(600.0);
  spectra_ids.push_back(pep_id);
  
  // Test assignment to base type
  PeptideIdentificationList base_vector = spectra_ids;
  TEST_EQUAL(base_vector.size(), 1)
  TEST_REAL_SIMILAR(base_vector[0].getRT(), 500.0)
  TEST_REAL_SIMILAR(base_vector[0].getMZ(), 600.0)
  
  // Test const reference compatibility
  const PeptideIdentificationList& const_ref = spectra_ids;
  TEST_EQUAL(const_ref.size(), 1)
  TEST_REAL_SIMILAR(const_ref[0].getRT(), 500.0)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST