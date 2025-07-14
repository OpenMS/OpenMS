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
  PeptideIdentificationList id_list;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  id_list.push_back(pep_id);
  
  PeptideIdentificationList pep_ids(id_list);
  TEST_EQUAL(pep_ids.size(), 1)
  TEST_REAL_SIMILAR(pep_ids[0].getRT(), 1234.5)
  TEST_REAL_SIMILAR(pep_ids[0].getMZ(), 567.8)
}
END_SECTION

START_SECTION((PeptideIdentificationList(PeptideIdentificationList&& vec)))
{
  PeptideIdentificationList id_list;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  id_list.push_back(pep_id);
  
  PeptideIdentificationList pep_ids(std::move(id_list));
  TEST_EQUAL(pep_ids.size(), 1)
  TEST_EQUAL(id_list.capacity(), 0)  // this is a move constructor, so the original id_list should be empty and not occupy any memory on the heap
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
  PeptideIdentificationList id_list;
  for (int i = 0; i < 3; ++i)
  {
    PeptideIdentification pep_id;
    pep_id.setRT(i * 100.0);
    pep_id.setMZ(i * 200.0);
    id_list.push_back(pep_id);
  }
  
  PeptideIdentificationList pep_ids(id_list.begin(), id_list.end());
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
  TEST_EQUAL(pep_ids1.capacity(), 0)  // this is a move constructor, so the original pep_ids1 should be empty and not occupy any memory on the heap
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
  PeptideIdentificationList id_list;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  id_list.push_back(pep_id);
  
  PeptideIdentificationList pep_ids;
  pep_ids = id_list;
  TEST_EQUAL(pep_ids.size(), 1)
  TEST_REAL_SIMILAR(pep_ids[0].getRT(), 1234.5)
  TEST_REAL_SIMILAR(pep_ids[0].getMZ(), 567.8)
}
END_SECTION

START_SECTION((PeptideIdentificationList& operator=(PeptideIdentificationList&& vec)))
{
  PeptideIdentificationList id_list;
  PeptideIdentification pep_id;
  pep_id.setRT(1234.5);
  pep_id.setMZ(567.8);
  id_list.push_back(pep_id);
  
  PeptideIdentificationList pep_ids;
  pep_ids = std::move(id_list);
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

START_SECTION((Test list-level score metadata methods))
{
  PeptideIdentificationList pep_ids;
  
  // Test default state
  TEST_EQUAL(pep_ids.usesListLevelScores(), false)
  TEST_EQUAL(pep_ids.getListScoreType(), "")
  TEST_EQUAL(pep_ids.isListHigherScoreBetter(), false)
  
  // Test setting list-level metadata
  pep_ids.setListScoreType("Mascot");
  pep_ids.setListHigherScoreBetter(true);
  pep_ids.enableListLevelScores(true);
  
  TEST_EQUAL(pep_ids.usesListLevelScores(), true)
  TEST_EQUAL(pep_ids.getListScoreType(), "Mascot")
  TEST_EQUAL(pep_ids.isListHigherScoreBetter(), true)
  
  // Test disabling list-level scores
  pep_ids.enableListLevelScores(false);
  TEST_EQUAL(pep_ids.usesListLevelScores(), false)
}
END_SECTION

START_SECTION((Test score consistency validation))
{
  PeptideIdentificationList pep_ids;
  
  // Test empty list consistency
  TEST_EQUAL(pep_ids.hasConsistentScoreMetadata(), true)
  pep_ids.validateScoreConsistency();  // Should not throw
  
  // Test single element consistency
  PeptideIdentification pep_id1;
  pep_id1.setScoreType("Mascot");
  pep_id1.setHigherScoreBetter(true);
  pep_ids.push_back(pep_id1);
  
  TEST_EQUAL(pep_ids.hasConsistentScoreMetadata(), true)
  pep_ids.validateScoreConsistency();  // Should not throw
  
  // Test consistent multiple elements
  PeptideIdentification pep_id2;
  pep_id2.setScoreType("Mascot");
  pep_id2.setHigherScoreBetter(true);
  pep_ids.push_back(pep_id2);
  
  TEST_EQUAL(pep_ids.hasConsistentScoreMetadata(), true)
  pep_ids.validateScoreConsistency();  // Should not throw
  
  // Test inconsistent score types
  PeptideIdentification pep_id3;
  pep_id3.setScoreType("Sequest");
  pep_id3.setHigherScoreBetter(true);
  pep_ids.push_back(pep_id3);
  
  TEST_EQUAL(pep_ids.hasConsistentScoreMetadata(), false)
  TEST_EXCEPTION(Exception::InvalidValue, pep_ids.validateScoreConsistency())
  
  // Fix inconsistency by changing score type
  pep_ids[2].setScoreType("Mascot");
  TEST_EQUAL(pep_ids.hasConsistentScoreMetadata(), true)
  
  // Test inconsistent score orientations
  pep_ids[2].setHigherScoreBetter(false);
  TEST_EQUAL(pep_ids.hasConsistentScoreMetadata(), false)
  TEST_EXCEPTION(Exception::InvalidValue, pep_ids.validateScoreConsistency())
}
END_SECTION

START_SECTION((Test migration from individual scores))
{
  PeptideIdentificationList pep_ids;
  
  // Test migration with empty list
  pep_ids.migrateFromIndividualScores();
  TEST_EQUAL(pep_ids.usesListLevelScores(), true)
  TEST_EQUAL(pep_ids.getListScoreType(), "")
  TEST_EQUAL(pep_ids.isListHigherScoreBetter(), false)
  
  // Reset and add consistent elements
  pep_ids.clear();
  pep_ids.enableListLevelScores(false);
  PeptideIdentification pep_id1, pep_id2;
  pep_id1.setScoreType("Mascot");
  pep_id1.setHigherScoreBetter(true);
  pep_id2.setScoreType("Mascot");
  pep_id2.setHigherScoreBetter(true);
  pep_ids.push_back(pep_id1);
  pep_ids.push_back(pep_id2);
  
  // Test successful migration
  pep_ids.migrateFromIndividualScores();
  TEST_EQUAL(pep_ids.usesListLevelScores(), true)
  TEST_EQUAL(pep_ids.getListScoreType(), "Mascot")
  TEST_EQUAL(pep_ids.isListHigherScoreBetter(), true)
  
  // Test migration with inconsistent scores (should throw)
  pep_ids.enableListLevelScores(false);
  pep_ids[1].setScoreType("Sequest");
  TEST_EXCEPTION(Exception::InvalidValue, pep_ids.migrateFromIndividualScores())
}
END_SECTION

START_SECTION((Test sync to individual scores))
{
  PeptideIdentificationList pep_ids;
  
  // Setup list with different individual scores
  PeptideIdentification pep_id1, pep_id2;
  pep_id1.setScoreType("Mascot");
  pep_id1.setHigherScoreBetter(true);
  pep_id2.setScoreType("Sequest");
  pep_id2.setHigherScoreBetter(false);
  pep_ids.push_back(pep_id1);
  pep_ids.push_back(pep_id2);
  
  // Set list-level metadata
  pep_ids.setListScoreType("Percolator");
  pep_ids.setListHigherScoreBetter(false);
  pep_ids.enableListLevelScores(true);
  
  // Sync to individual scores
  pep_ids.syncToIndividualScores();
  
  // Check that individual scores were updated
  TEST_EQUAL(pep_ids[0].getScoreType(), "Percolator")
  TEST_EQUAL(pep_ids[0].isHigherScoreBetter(), false)
  TEST_EQUAL(pep_ids[1].getScoreType(), "Percolator")
  TEST_EQUAL(pep_ids[1].isHigherScoreBetter(), false)
  
  // Test sync when list-level scores are disabled (should be no-op)
  pep_ids.enableListLevelScores(false);
  pep_ids.setListScoreType("NewScore");
  pep_ids.syncToIndividualScores();
  
  // Should not have changed
  TEST_EQUAL(pep_ids[0].getScoreType(), "Percolator")
  TEST_EQUAL(pep_ids[1].getScoreType(), "Percolator")
}
END_SECTION

START_SECTION((Test file I/O helper methods))
{
  PeptideIdentificationList pep_ids;
  
  // Test tryMigrateAfterLoad with empty list
  pep_ids.tryMigrateAfterLoad();
  TEST_EQUAL(pep_ids.usesListLevelScores(), false)
  
  // Test tryMigrateAfterLoad with consistent scores
  PeptideIdentification pep_id1, pep_id2;
  pep_id1.setScoreType("Mascot");
  pep_id1.setHigherScoreBetter(true);
  pep_id2.setScoreType("Mascot");
  pep_id2.setHigherScoreBetter(true);
  pep_ids.push_back(pep_id1);
  pep_ids.push_back(pep_id2);
  
  pep_ids.tryMigrateAfterLoad();
  TEST_EQUAL(pep_ids.usesListLevelScores(), true)
  TEST_EQUAL(pep_ids.getListScoreType(), "Mascot")
  TEST_EQUAL(pep_ids.isListHigherScoreBetter(), true)
  
  // Test tryMigrateAfterLoad with inconsistent scores (should not migrate)
  pep_ids.enableListLevelScores(false);
  pep_ids[1].setScoreType("Sequest");
  pep_ids.tryMigrateAfterLoad();
  TEST_EQUAL(pep_ids.usesListLevelScores(), false)
  
  // Test force migration with inconsistent scores
  TEST_EXCEPTION(Exception::InvalidValue, pep_ids.tryMigrateAfterLoad(true))
  
  // Test prepareForSave
  pep_ids[1].setScoreType("Mascot");  // Fix inconsistency
  pep_ids.migrateFromIndividualScores();
  pep_ids.setListScoreType("NewScore");
  pep_ids.setListHigherScoreBetter(false);
  
  pep_ids.prepareForSave();
  
  // Check that individual scores were synced
  TEST_EQUAL(pep_ids[0].getScoreType(), "NewScore")
  TEST_EQUAL(pep_ids[0].isHigherScoreBetter(), false)
  TEST_EQUAL(pep_ids[1].getScoreType(), "NewScore")
  TEST_EQUAL(pep_ids[1].isHigherScoreBetter(), false)
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST