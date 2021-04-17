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

FASTAFile* ptr = nullptr; //wird der konstruktor getestet
START_SECTION((FASTAFile()))
  ptr = new FASTAFile();
  TEST_EQUAL(ptr == nullptr, false)
END_SECTION

START_SECTION((~FASTAFile())) //wird der destruktor getestet
  delete(ptr);
END_SECTION

START_SECTION([FASTAFile::FASTAEntry] FASTAEntry())//fasta entry wird getestet
  FASTAFile::FASTAEntry * ptr_e;
  ptr_e = new FASTAFile::FASTAEntry();
  TEST_EQUAL(ptr_e == nullptr, false)
END_SECTION

START_SECTION([FASTAFile::FASTAEntry] FASTAEntry(String id, String desc, String seq))
  FASTAFile::FASTAEntry entry("ID", "DESC", "DAVLDELNER"); //fasta entry wird befüllt
  TEST_EQUAL(entry.identifier, "ID") //hat das befüllen richtig geklappt
  TEST_EQUAL(entry.description, "DESC")
  TEST_EQUAL(entry.sequence, "DAVLDELNER")
END_SECTION

START_SECTION([FASTAFile::FASTAEntry] bool operator==(const FASTAEntry &rhs) const)// == operator wird getestet
  FASTAFile::FASTAEntry entry1("ID", "DESC", "DAV*LDELNER");
  FASTAFile::FASTAEntry entry2("ID", "DESC", "DAV*LDELNER");
  FASTAFile::FASTAEntry entry3("ID2", "DESC", "DAV*LDELNER");
  TEST_EQUAL(entry1==entry2, true)
  TEST_EQUAL(entry1==entry3, false)
END_SECTION


START_SECTION((void load(const String& filename, std::vector< FASTAEntry > &data)))//load funktion wird getestet
  vector<FASTAFile::FASTAEntry> data;
  FASTAFile file;

  TEST_EXCEPTION(Exception::FileNotFound, file.load("FASTAFile_test_this_file_does_not_exist",data))

  file.load(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"),data); //in diese Datei auch falsche Zeichen einfügen und einen peff header?//in data ligt jetzt die test fasta datei
  vector< FASTAFile::FASTAEntry >::const_iterator sequences_iterator = data.begin();//zum ersten protein gehen
  TEST_EQUAL(data.size(), 5)//5 proteine
  TEST_EQUAL(sequences_iterator->identifier, String("P68509|1433F_BOVIN")) //erstes protein vergleichen
  TEST_EQUAL(sequences_iterator->description, String("This is the description of the first protein"))
  TEST_EQUAL(sequences_iterator->sequence, String("GDREQLLQRARLAEQAERYDDMASAMKAVTEL") +
    String("NEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVC") +
    String("NDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEIS") +
    String("KEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQL") +
    String("LRDNLTLWTSDQQDEEAGEGN"))
  sequences_iterator++; //zum nächsten protein
  TEST_EQUAL(sequences_iterator->identifier, "Q9CQV8|1433B_MOUSE")
  TEST_EQUAL(sequences_iterator->sequence, String("TMDKSELVQKAKLAEQAERYDDMAAAMKAVTE") +
    String("QGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICND") +
    String("VLELLDKYLILNATQAESKVFYLKMKGDYFRYLSEVASGENKQTTVSNSQQAYQEAFEISKKEMQ") +
    String("PTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLT") +
    String("LWTSENQGDEGDAGEGEN"))
  sequences_iterator++; //in diese sequenz wurden falsche zeichen eingefügt und es wird getestet ob sie entfernt wurden
  TEST_EQUAL(sequences_iterator->identifier, "sp|P31946|1433B_HUMAN")
  TEST_EQUAL(sequences_iterator->description, String("14-3-3 protein beta/alpha OS=Homo sapiens GN=YWHAB PE=1 SV=3"))
  TEST_EQUAL(sequences_iterator->sequence, String("MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSS")
    + String("WRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFY")
    + String("LKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFY")
    + String("YEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGD")
    + String("AGEGEN"))

  sequences_iterator++;//hier wurde ein PEFF header angefügt, der übersprungen werden soll
  TEST_EQUAL(sequences_iterator->identifier, "sp|P00000|0000A_UNKNOWN")
  TEST_EQUAL(sequences_iterator->description, String("Artificially modified version of sp|P31946|1433B_HUMAN"))
  TEST_EQUAL(sequences_iterator->sequence, String("(ICPL:13C(6))MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSS")
    + String("WRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFY")
    + String("LKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFY")
    + String("YEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGD")
    + String("AGEGEN"))

  // test if the modifed sequence is convertable
  AASequence aa = AASequence::fromString(sequences_iterator->sequence);//wieso war die modifiziert vorher? woran erkennt man das?-> im header!
  TEST_EQUAL(aa.toUnmodifiedString(), String("MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSS")
    + String("WRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFY")
    + String("LKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFY")
    + String("YEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGD")
    + String("AGEGEN"))

  TEST_EQUAL(aa.isModified(), true) //wieso soll das true sein wenn es vorher zu unmodified converted wurde
  String expectedModification = ModificationsDB::getInstance()->getModification("ICPL:13C(6)", "", ResidueModification::N_TERM)->getId();
  TEST_EQUAL(aa.getNTerminalModificationName(), expectedModification)

  sequences_iterator++; //last sequence, wenn fasta format keine zeilenumbrüche hat
  TEST_EQUAL(sequences_iterator->identifier, "test")
  TEST_EQUAL(sequences_iterator->description, String(" ##0"))
  TEST_EQUAL(sequences_iterator->sequence, String("GSMTVDMQEIGSTEMPYEVPTQPNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSFKEKVSELVS")
    + String("PAVYTFGLFVQNASESLTSDDPSDVPTQRTFKSDFQSVGSMTVDMQEIGSTEMPYEVPTQ")
    + String("PNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSFKEKVSELVSPAVYTFGL")
    + String("FVQNASESLTSDDPSDVPTQRTFKSDFQSVAXXSTFDFYQRRLVTLAESPRAPSPGSMTV")
    + String("DMQEIGSTEMPYEVPTQPNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSF")
    + String("KEKVSELVSPAVYTFGLFVQNASESLTSDDPSDVPTQRTFKSDFQSV"))

END_SECTION

START_SECTION((void store(const String& filename, const std::vector< FASTAEntry > &data) const))
  vector<FASTAFile::FASTAEntry> data, data2;
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  FASTAFile file;

  file.load(OPENMS_GET_TEST_DATA_PATH("FASTAFile_test.fasta"),data);
  TEST_EXCEPTION(Exception::UnableToCreateFile, file.store("/bla/bluff/blblb/sdfhsdjf/test.txt",data))

  file.store(tmp_filename,data);
  file.load(tmp_filename,data2);
  TEST_EQUAL(data==data2,true); // vectoren mit fasta entries auf gleichheit testen
END_SECTION

/*
START_SECTION([EXTRA] test_strange_symbols_in_sequence)
  // test if * is read correctly (not changed into something weird like 'X') //diesen Teil hier anpassen bzw löschen, weil das oben getestet wurde und * nicht erlaubt sein sollte und außerdem während des einlesens gecheckt wird und nicht erst danach beim pushbacken
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
  */

END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST

