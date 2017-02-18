// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2016.
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
// $Authors: Marie Hoffmann $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/FASTQFile.h>
#include <OpenMS/DATASTRUCTURES/FASTQEntry.h>
#include <OpenMS/DATASTRUCTURES/String.h>
#include <OpenMS/CHEMISTRY/ModificationsDB.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

#include <vector>
#include <string>
#include <iostream>
#include <ios>

///////////////////////////

START_TEST(FASTQFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;

FASTQFile* ptr = 0;
START_SECTION((FASTQFile()))
  ptr = new FASTQFile();
  TEST_EQUAL(ptr == 0, false)
END_SECTION


START_SECTION((~FASTQFile()))
  delete(ptr);
END_SECTION


FASTQFile file;
std::vector< FASTQEntry > sequences;
std::vector< FASTQEntry >::const_iterator sequences_iterator;
FASTQEntry temp_entry;


START_SECTION([FASTQEntry] FASTQEntry())
  FASTQEntry * ptr_e;
  ptr_e = new FASTQEntry();
  TEST_EQUAL(ptr_e == 0, false)
END_SECTION


START_SECTION([FASTQEntry] FASTQEntry(seqan::CharString id, seqan::CharString desc, seqan::CharString seq))
  FASTQEntry entry("ID", "DESC", "DAVLDELNER", "@@IIJJ!!((");
  TEST_EQUAL(entry.identifier, "ID")
  TEST_EQUAL(entry.description, "DESC")
  TEST_EQUAL(entry.sequence, "DAVLDELNER")
  TEST_EQUAL(entry.quality, "@@IIJJ!!((")
END_SECTION

START_SECTION([FASTQEntry] bool operator==(const FASTQEntry &rhs) const)
  FASTQEntry entry1("ID", "DESC", "DAVLDELNER", "@@IIJJ!!((");
  FASTQEntry entry2("ID", "DESC", "DAVLDELNER", "@@IIJJ!!((");
  FASTQEntry entry3("ID2", "DESC", "DAVLDELNER", "@@IIJJ!!((");
  TEST_EQUAL(entry1 == entry2, true)
  TEST_EQUAL(entry1 == entry3, false)
END_SECTION

START_SECTION([FASTQEntry] std::vector<int> qual2phred())
  FASTQEntry entry("ID", "DESC", "DAVLDELNERDAVLDELNERDAVLDELNERDAVLDELNERAA", "!\"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJ");
  std::vector<int> phreds = entry.qual2phred();
  for (int i = 0; i < length(entry.quality); ++i)
    TEST_EQUAL(i == phreds[i], true);
END_SECTION

temp_entry.identifier = seqan::CharString("P68509|1433F_BOVIN");
temp_entry.description = seqan::CharString("This is the description of the first protein");
temp_entry.sequence = seqan::CharString("GDREQLLQRARLAEQAERYDDMASAMKAVTELNEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVCNDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEISKEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDQQDEEAGEGN");
temp_entry.quality =  seqan::CharString("7?'@?J.6'%$99F6.%->I8@1.4-2JG;A(+8(F*<;=0#!5C6:J'>B:IAFI@BHI3-6;7DCA>&-94@'G)&!+-@2/@)+?985I?D6!I@GH7EII\"%#J8=J(ED)!@<<B83?3B1>6B-4,+/E7=!91!H,-4F>14''2?&B(A&49)/)(8HA13B-><<*E6A,C-I7\">$'H5$%IF0@B/J#8C%+=/3;6(\";6F0%?75,<;!(<933JA)A%7GC?F<H\"C25(>");

sequences.push_back(temp_entry);


START_SECTION((void load(const String& filename, std::vector< FASTQEntry > &data)))
  std::vector<FASTQEntry> data;

  TEST_EXCEPTION(Exception::FileNotFound, file.load("FASTQFile_test_this_file_does_not_exist", data))
  file.load(OPENMS_GET_TEST_DATA_PATH("FASTQFile_test.fastq"), data);
  sequences_iterator = data.begin();
  TEST_EQUAL(data.size(), 5)

  TEST_EQUAL(sequences_iterator->identifier, seqan::CharString("P68509|1433F_BOVIN"))
  TEST_EQUAL(sequences_iterator->description, seqan::CharString("protein 1, quality=one_letter"))
  TEST_EQUAL(sequences_iterator->sequence, seqan::CharString("GDREQLLQRARLAEQAERYDDMASAMKAVTELNEPLSNEDRNLLSVAYKNVVGARRSSWRVISSIEQKTMADGNEKKLEKVKAYREKIEKELETVCNDVLALLDKFLIKNCNDFQYESKVFYLKMKGDYYRYLAEVASGEKKNSVVEASEAAYKEAFEISKEHMQPTHPIRLGLALNFSVFYYEIQNAPEQACLLAKQAFDDAIAELDTLNEDSYKDSTLIMQLLRDNLTLWTSDQQDEEAGEGN"))
  TEST_EQUAL(sequences_iterator->quality, seqan::CharString("IIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIIII"))

  ++sequences_iterator;

  TEST_EQUAL(sequences_iterator->identifier, seqan::CharString("Q9CQV8|1433B_MOUSE"))
  TEST_EQUAL(sequences_iterator->description, seqan::CharString("protein 2, quality=five_letters"))
  TEST_EQUAL(sequences_iterator->sequence, seqan::CharString("TMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLILNATQAESKVFYLKMKGDYFRYLSEVASGENKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGDAGEGEN"))
  TEST_EQUAL(sequences_iterator->quality, seqan::CharString("ABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDEABCDE"))

  ++sequences_iterator;

  TEST_EQUAL(sequences_iterator->identifier, seqan::CharString("sp|P31946|1433B_HUMAN"))
  TEST_EQUAL(sequences_iterator->description, seqan::CharString("protein 3, quality='!' to 'J', no_quotes_no_at"))
  TEST_EQUAL(sequences_iterator->sequence, seqan::CharString("MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFYLKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGDAGEGEN"))
  TEST_EQUAL(sequences_iterator->quality, seqan::CharString("F:6C&+C2CA15=&36+##H+#DC(2:#>B5A:D&<9:,<9-;*GA9D,I#%568F-<%<>:A5F*=$<BC;9AD*6-#+7D7-$8,4D;I#)!4I-E$<B4B*B>4=8><+H*,5-B9E;7DG81C<%;+-$$CI$3F+!&C;04(30+&?BDB1=F-C**967F&H.C-5JJA4B(2I2!E&,0A0;*D-D2-6=3!G)-,4III:C:DE.6!DF(8D=442!,;DF<BB?4DB)%G13E$!3)"));

  ++sequences_iterator;

  TEST_EQUAL(sequences_iterator->identifier, seqan::CharString("sp|P00000|0000A_UNKNOWN"))
  TEST_EQUAL(sequences_iterator->description, seqan::CharString("protein 4, quality='!'-'J' no_quotes"))
  TEST_EQUAL(sequences_iterator->sequence, seqan::CharString("MTMDKSELVQKAKLAEQAERYDDMAAAMKAVTEQGHELSNEERNLLSVAYKNVVGARRSSWRVISSIEQKTERNEKKQQMGKEYREKIEAELQDICNDVLELLDKYLIPNATQPESKVFYLKMKGDYFRYLSEVASGDNKQTTVSNSQQAYQEAFEISKKEMQPTHPIRLGLALNFSVFYYEILNSPEKACSLAKTAFDEAIAELDTLNEESYKDSTLIMQLLRDNLTLWTSENQGDEGDAGEGEN"))
  TEST_EQUAL(sequences_iterator->quality, seqan::CharString("4%)0++4I.<7!%--#3@9B51&A-1CE(+JC%6*9-B16D<B%:CC?G$B;2$C;:2+..H=$48.4D+?6!F0<-;*B2?*B45EG2)>7=7%AGA%0DD$C)>CGA&?>4;FF>-CB-.GB=)?3E$<=DD@):H9=.4E$D6@1!+$C;I<<A*IAD<-D3B=H6H<>)2*51$7?D>+$9>3.2*H0@4%($B@H-;!!D3I0IHBBD8DAI#H?<(-#F+?G#@2)#5%1,06J2HA%D%"))

  ++sequences_iterator;
  TEST_EQUAL(sequences_iterator->identifier, seqan::CharString("test"))
  TEST_EQUAL(sequences_iterator->description, seqan::CharString("protein 5, quality='!' to 'J'"))
  TEST_EQUAL(sequences_iterator->sequence, seqan::CharString("GSMTVDMQEIGSTEMPYEVPTQPNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSFKEKVSELVSPAVYTFGLFVQNASESLTSDDPSDVPTQRTFKSDFQSVGSMTVDMQEIGSTEMPYEVPTQPNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSFKEKVSELVSPAVYTFGLFVQNASESLTSDDPSDVPTQRTFKSDFQSVAXXSTFDFYQRRLVTLAESPRAPSPGSMTVDMQEIGSTEMPYEVPTQPNATSASAGRGWFDGPSFKVPSVPTRPSGIFRRPSRIKPEFSFKEKVSELVSPAVYTFGLFVQNASESLTSDDPSDVPTQRTFKSDFQSV"))
  TEST_EQUAL(sequences_iterator->quality, seqan::CharString("->1F:%.#B+,3?CHH;@G1\"+3=2*CJ-J>1<:3\"F\"E\"2BG.06:E!I!<)>9/\")'53)7&7*#@J#-3!DH4I<-2F0D-5AI9#->-5?ED30:$=E8?5/5J4H!JD71>::0%7%(I8#9'5)D%?6/#3#9<!3J.BC8!5D:&2,G1$,$0<83>62*H.,-'5!6CI9(&FH6H,F7B;;(AB-'%$;*$-&1=/43.(J-0G9JB@A-=2D,D+4341B>,9H:/II\"J1(@%D1=A6%7HBB!9\")D60;D7,44-8?C>1B4G#%-,BI/78->CE?<.9.$&H;\">%CA*%J%9IH&.8#H5GH:<@&'>@8,6IHJ@;4'J;G@+A+(4&8*8+G6('4E:!E=5-"))

END_SECTION


START_SECTION((void store(const String& filename, const std::vector< FASTQEntry > &data) const))
  std::vector<FASTQEntry> data, data2;
  String tmp_filename;
  NEW_TMP_FILE(tmp_filename);
  FASTQFile file;
  file.load(OPENMS_GET_TEST_DATA_PATH("FASTQFile_test.fastq"), data);
  TEST_EXCEPTION(Exception::UnableToCreateFile, file.store("/bla/bluff/blblb/sdfhsdjf/test.txt", data))
  file.store(tmp_filename, data);
  file.load(tmp_filename, data2);
  TEST_EQUAL(data == data2, true);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
