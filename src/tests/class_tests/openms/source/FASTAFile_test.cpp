// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Nico Pfeifer, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <string>

#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <vector>

///////////////////////////

START_TEST(FASTAFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

FASTAFile *ptr = nullptr;

START_SECTION((FASTAFile()))
  ptr = new FASTAFile();
  TEST_EQUAL(ptr == nullptr, false)
END_SECTION

START_SECTION((~FASTAFile()))
  delete (ptr);
END_SECTION

START_SECTION([FASTAFile::FASTAEntry] FASTAEntry())
  FASTAFile::FASTAEntry *ptr_e;
  ptr_e = new FASTAFile::FASTAEntry();
  TEST_EQUAL(ptr_e == nullptr, false)
  delete ptr_e;
END_SECTION


START_SECTION([FASTAFile::FASTAEntry] FASTAEntry(String id, String desc, String seq))
  FASTAFile::FASTAEntry entry("ID", "DESC", "DAVLDELNER");
  TEST_EQUAL(entry.identifier, "ID")
  TEST_EQUAL(entry.description, "DESC")
  TEST_EQUAL(entry.sequence, "DAVLDELNER")
END_SECTION

START_SECTION([FASTAFile::FASTAEntry] bool operator==(const FASTAEntry &rhs) const)
  FASTAFile::FASTAEntry entry1("ID", "DESC", "DAV*LDELNER");
  FASTAFile::FASTAEntry entry2("ID", "DESC", "DAV*LDELNER");
  FASTAFile::FASTAEntry entry3("ID2", "DESC", "DAV*LDELNER");
  TEST_TRUE(entry1 == entry2)
  TEST_EQUAL(entry1 == entry3, false)
END_SECTION


START_SECTION((void
                  load(
                  const String &filename, std::vector<FASTAEntry>
                  &data)))
  vector<FASTAFile::FASTAEntry> data;
  FASTAFile file;

  TEST_EXCEPTION(Exception::FileNotFound, file.load("FASTAFile_test_this_file_does_not_exist", data))

  file.load(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"), data);
  vector<FASTAFile::FASTAEntry>::const_iterator sequences_iterator = data.begin();
  TEST_EQUAL(data.size(), 5)
  TEST_EQUAL(sequences_iterator->identifier, String("P68509|1433F_BOVIN"))
  TEST_EQUAL(sequences_iterator->description, String("This is the description of the first protein"))
  TEST_EQUAL(sequences_iterator->sequence, String("GDREQLLQRARLAEQAERYDDMASAMKAVTEL") +
                                           String(
                                               "NEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVC") +
                                           String(
                                               "NDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEIS") +
                                           String(
                                               "KEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQL") +
                                           String("LRDNLTLWTSDQQDEEAGEGN"))
  sequences_iterator++;
  TEST_EQUAL(sequences_iterator->identifier, "Q9CQV8|1433B_MOUSE")
  TEST_EQUAL(sequences_iterator->sequence, String("TMDKSELVQKAKLAEQAERYDDMAAAMKAVTE") +
                                           String(
                                               "QGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICND") +
                                           String(
                                               "VLELLDKYLILNATQAESKVFYLKMKGDYFRYLSEVASGENKQTTVSNSQQAYQEAFEISKKEMQ") +
                                           String(
                                               "PTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLT") +
                                           String("LWTSENQGDEGDAGEGEN"))
  sequences_iterator++;
  TEST_EQUAL(sequences_iterator->identifier, "sp|P31946|1433B_HUMAN")
  TEST_EQUAL(sequences_iterator->description,
             String("14-3-3 protein beta/alpha OS=Homo sapiens GN=YWHAB PE=1 SV=3"))
  TEST_EQUAL(sequences_iterator->sequence, String("MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSS")
                                           +
                                           String("WRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFY")
                                           +
                                           String("LKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFY")
                                           +
                                           String("YEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGD")
                                           + String("AGEGEN"))

  sequences_iterator++;
  TEST_EQUAL(sequences_iterator->identifier, "sp|P00000|0000A_UNKNOWN")
  TEST_EQUAL(sequences_iterator->description, String("Artificially modified version of sp|P31946|1433B_HUMAN"))
  TEST_EQUAL(sequences_iterator->sequence,
             String("(ICPL:13C(6))MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSS")
             + String("WRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFY")
             + String("LKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFY")
             + String("YEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGD")

             + String("AGEGEN"))




  // test if the modifed sequence is convertable
  AASequence aa = AASequence::fromString(sequences_iterator->sequence);
  TEST_EQUAL(aa.toUnmodifiedString(), String("MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSS")
                                      + String("WRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFY")
                                      + String("LKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFY")
                                      + String("YEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGD")
                                      + String("AGEGEN"))

  TEST_EQUAL(aa.isModified(), true)
  String expectedModification = ModificationsDB::getInstance()->getModification("ICPL:13C(6)", "",
                                                                                ResidueModification::N_TERM)->getId();
  TEST_EQUAL(aa.getNTerminalModificationName(), expectedModification)

  sequences_iterator++;
  TEST_EQUAL(sequences_iterator->identifier, "test")
  TEST_EQUAL(sequences_iterator->description, String(" ##0"))
  TEST_EQUAL(sequences_iterator->sequence,
             String("GSMTVDMQEIGSTEMPYEVPTQPNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSFKEKVSELVS")
             + String("PAVYTFGLFVQNASESLTSDDPSDVPTQRTFKSDFQSVGSMTVDMQEIGSTEMPYEVPTQ")
             + String("PNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSFKEKVSELVSPAVYTFGL")
             + String("FVQNASESLTSDDPSDVPTQRTFKSDFQSVAXXSTFDFYQRRLVTLAESPRAPSPGSMTV")
             + String("DMQEIGSTEMPYEVPTQPNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSF")
             + String("KEKVSELVSPAVYTFGLFVQNASESLTSDDPSDVPTQRTFKSDFQSV"))
END_SECTION


START_SECTION((void
                  store(
                  const String &filename,
                  const std::vector<FASTAEntry> &data) const))
  vector<FASTAFile::FASTAEntry> data, data2;
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  FASTAFile file;

  file.load(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"), data);
  TEST_EXCEPTION(Exception::UnableToCreateFile, file.store("/bla/bluff/blblb/sdfhsdjf/test.txt", data))

  file.store(tmp_filename, data);

  file.load(tmp_filename, data2);

  file.load(tmp_filename, data2);

  TEST_TRUE(data == data2);
END_SECTION


START_SECTION([EXTRA] test_strange_symbols_in_sequence)
  // test if * is read correctly (not changed into something weird like 'X')
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  FASTAFile file;
  vector<FASTAFile::FASTAEntry> data, data2;
  FASTAFile::FASTAEntry temp_entry;
  temp_entry.identifier = String("P68509|1433F_BOVIN");
  temp_entry.description = String("This is the description of the first protein");
  temp_entry.sequence = String("GDREQLLQRAR*LAEQ*AERYDDMASAMKAVTEL");
  data.push_back(temp_entry);
  data.push_back(temp_entry); // twice

  file.store(tmp_filename, data);
  file.load(tmp_filename, data2);

  ABORT_IF(data2.size() != 2);
  TEST_EQUAL(data2[0] == temp_entry, true);
  TEST_EQUAL(data2[1] == temp_entry, true);
END_SECTION


START_SECTION([EXTRA] test_strange_symbols_in_sequence)
  // test if * is read correctly (not changed into something weird like 'X')
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  FASTAFile file;
  vector<FASTAFile::FASTAEntry> data, data2;
  FASTAFile::FASTAEntry temp_entry;
  temp_entry.identifier = String("P68509|1433F_BOVIN");
  temp_entry.description = String("This is the description of the first protein");
  temp_entry.sequence = String("GDREQLLQRAR*LAEQ*AERYDDMASAMKAVTEL");
  data.push_back(temp_entry);
  data.push_back(temp_entry); // twice


  file.store(tmp_filename, data);
  file.load(tmp_filename, data2);

  ABORT_IF(data2.size() != 2);
  TEST_EQUAL(data2[0] == temp_entry, true);
  TEST_EQUAL(data2[1] == temp_entry, true);
END_SECTION


START_SECTION(test_white_spaces)
  //test if spaces and tabulators are removed correctly
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  FASTAFile file;
  vector<FASTAFile::FASTAEntry> data, data2;

  FASTAFile::FASTAEntry temp_entry;
  temp_entry.identifier = String("P68509|1433F_BOVIN");
  temp_entry.description = String("This is the description of the first protein");
  temp_entry.sequence = String("GDREQLLQRAR LAEQ\tAERYDDMASAMKAVTEL");
  data.push_back(temp_entry);

  file.store(tmp_filename, data);
  file.load(tmp_filename, data2);

  ABORT_IF(data2.size() != 1);
  TEST_EQUAL(data2[0].sequence == string("GDREQLLQRARLAEQAERYDDMASAMKAVTEL"), true);
END_SECTION

START_SECTION([EXTRA] test_position)
  // test if setPosition() works correctly
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  FASTAFile file;

  vector<pair<streampos, FASTAFile::FASTAEntry>> data1;
  vector<pair<streampos, FASTAFile::FASTAEntry>> data2;
  FASTAFile::FASTAEntry temp_entry;
  file.readStart(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"));


  for (int i = 0; i < 4; i++) {
    file.readNext(temp_entry);
    data1.push_back(std::make_pair(file.position(), temp_entry));
  }
  file.setPosition(data1[0].first);
  for (int i = 0; i < 3; i++) {
    file.readNext(temp_entry);
    data2.push_back(std::make_pair(file.position(), temp_entry));
  }

  ABORT_IF(data1.size() != 4 || data2.size() != 3);

  for (Size i = 1; i < data1.size(); i++) {
    TEST_EQUAL(data1[i].second.identifier, data2[i - 1].second.identifier);
    TEST_EQUAL(data1[i].second.description, data2[i - 1].second.description);
    TEST_EQUAL(data1[i].second.sequence, data2[i - 1].second.sequence);
    TEST_EQUAL(data1[i].first, data2[i - 1].first);
  }
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
