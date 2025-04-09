// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: $
// $Authors: $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/QcMLFile.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(QcMLFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

QcMLFile* ptr = 0;
QcMLFile* null_ptr = 0;
START_SECTION(QcMLFile())
{
	ptr = new QcMLFile();
	TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~QcMLFile())
{
	delete ptr;
}
END_SECTION

START_SECTION((~QcMLFile()))
{
  // uh, twice?! No!
}
END_SECTION

QcMLFile qcmlfile;
qcmlfile.registerRun("123","testrun1");
qcmlfile.registerRun("456","testrun2");


START_SECTION((void registerRun(const String id, const String name)))
{
  QcMLFile qcmlfile1;
  qcmlfile1.registerRun("abc","somerun");
  TEST_EQUAL(qcmlfile1.existsRun("abc"), true)
  TEST_EQUAL(qcmlfile1.existsRun("somerun",true), true)
}
END_SECTION


START_SECTION((void registerSet(const String id, const String name, const std::set< String > &names)))
{
  QcMLFile qcmlfile1;
  std::set<String> n;
  n.insert("somerun1");
  n.insert("somerun2");
  qcmlfile1.registerSet("def","someset", n);
  TEST_EQUAL(qcmlfile1.existsSet("def"), true)
  TEST_EQUAL(qcmlfile1.existsSet("someset",true), true)
}
END_SECTION

START_SECTION((void addRunQualityParameter(String r, QualityParameter qp)))
NOT_TESTABLE
END_SECTION

START_SECTION((void addRunAttachment(String r, Attachment at)))
NOT_TESTABLE
END_SECTION

START_SECTION((void addSetQualityParameter(String r, QualityParameter qp)))
NOT_TESTABLE
END_SECTION

START_SECTION((void addSetAttachment(String r, Attachment at)))
NOT_TESTABLE
END_SECTION

START_SECTION((void removeAttachment(String r, std::vector< String > &ids, String at="")))
NOT_TESTABLE
END_SECTION

START_SECTION((void removeAttachment(String r, String at)))
NOT_TESTABLE
END_SECTION

START_SECTION((void removeAllAttachments(String at)))
NOT_TESTABLE
END_SECTION

START_SECTION((void removeQualityParameter(String r, std::vector< String > &ids)))
NOT_TESTABLE
END_SECTION

START_SECTION((void merge(const QcMLFile &addendum, String setname="")))
NOT_TESTABLE
END_SECTION

START_SECTION((void collectSetParameter(const String setname, const String qp, std::vector< String > &ret)))
NOT_TESTABLE
END_SECTION

START_SECTION((String exportAttachment(const String filename, const String qpname) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((String exportQP(const String filename, const String qpname) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((String exportQPs(const String filename, const StringList qpnames) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((String map2csv(const std::map< String, std::map< String, String > > &cvs_table, const String &separator) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((String exportIDstats(const String &filename) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((void getRunIDs(std::vector< String > &ids) const ))
{
  std::vector< String > ids, test;
  test.push_back("123");
  test.push_back("456");
  qcmlfile.getRunIDs(ids);
  TEST_EQUAL(StringList(ids),StringList(test))
}
END_SECTION

START_SECTION((void getRunNames(std::vector< String > &ids) const ))
{
  std::vector< String > ids, test;
  test.push_back("testrun1");
  test.push_back("testrun2");
  qcmlfile.getRunNames(ids);
  TEST_EQUAL(ids==test,true)
}
END_SECTION

START_SECTION((bool existsRun(const String filename, bool checkname=false) const ))
{
  QcMLFile qcmlfile1;
  qcmlfile.registerRun("abc","somerun");
  TEST_EQUAL(qcmlfile.existsRun("abc"), true)
  TEST_EQUAL(qcmlfile.existsRun("somerun",true), true)
}
END_SECTION

START_SECTION((bool existsSet(const String filename, bool checkname=false) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((void existsRunQualityParameter(const String filename, const String qpname, std::vector< String > &ids) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((void existsSetQualityParameter(const String filename, const String qpname, std::vector< String > &ids) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((void store(const String &filename) const ))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((void load(const String &filename)))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(([QcMLFile::Attachment] Attachment()))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(([QcMLFile::Attachment] Attachment(const Attachment &rhs)))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(([QcMLFile::Attachment] Attachment& operator=(const Attachment &rhs)))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(([QcMLFile::Attachment] bool operator==(const Attachment &rhs) const ))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(([QcMLFile::Attachment] bool operator<(const Attachment &rhs) const ))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(([QcMLFile::Attachment] bool operator>(const Attachment &rhs) const ))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(([QcMLFile::Attachment] String toXMLString(UInt indentation_level) const ))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(([QcMLFile::Attachment] String toCSVString(String separator) const ))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(([QcMLFile::QualityParameter] QualityParameter()))
	NOT_TESTABLE
END_SECTION

START_SECTION(([QcMLFile::QualityParameter] QualityParameter(const QualityParameter &rhs)))
{
  QcMLFile::QualityParameter qp1;
  qp1.name = "somename"; ///< Name
  qp1.id = "id"; ///< Identifier
  qp1.cvRef = "MS"; ///< cv reference
  qp1.cvAcc = "MS:1000577";
  qp1.value = "somevalue";
  
  QcMLFile::QualityParameter qp2 = QcMLFile::QualityParameter(qp1);
  TEST_EQUAL(qp1.name, qp2.name)
  TEST_EQUAL(qp1.id, qp2.id)
  TEST_EQUAL(qp1.value, qp2.value)
}
END_SECTION

START_SECTION(([QcMLFile::QualityParameter] QualityParameter& operator=(const QualityParameter &rhs)))
{
  QcMLFile::QualityParameter qp1;
  qp1.name = "somename"; ///< Name
  qp1.id = "id"; ///< Identifier
  qp1.cvRef = "MS"; ///< cv reference
  qp1.cvAcc = "MS:1000577";
  qp1.value = "somevalue";
  
  QcMLFile::QualityParameter qp2 = QcMLFile::QualityParameter(qp1);
  qp2.name = "someothername"; ///< Name
  qp2.id = "otherid"; ///< Identifier
  qp2.cvRef = "MS"; ///< cv reference
  qp2.cvAcc = "MS:1000577";
  qp2.value = "someothervalue";
  
  qp2 = qp1;
  
  TEST_EQUAL(qp1.name, qp2.name)
  TEST_EQUAL(qp1.id, qp2.id)
  TEST_EQUAL(qp1.value, qp2.value)
}
END_SECTION

START_SECTION(([QcMLFile::QualityParameter] bool operator==(const QualityParameter &rhs) const ))
{
  QcMLFile::QualityParameter qp1;
  qp1.name = "somename"; ///< Name
  
  QcMLFile::QualityParameter qp2 = QcMLFile::QualityParameter(qp1);
  TEST_TRUE(qp1 == qp2)
}
END_SECTION

START_SECTION(([QcMLFile::QualityParameter] bool operator<(const QualityParameter &rhs) const ))
{
	QcMLFile::QualityParameter qp1;
  qp1.name = "somename"; ///< Name
  
  QcMLFile::QualityParameter qp2;
  qp2.name = "tomename"; ///< Name

  TEST_EQUAL(qp1<qp2, true)
}
END_SECTION

START_SECTION(([QcMLFile::QualityParameter] bool operator>(const QualityParameter &rhs) const ))
{
	QcMLFile::QualityParameter qp1;
  qp1.name = "somename"; ///< Name
  
  QcMLFile::QualityParameter qp2;
  qp2.name = "romename"; ///< Name

  TEST_EQUAL(qp1>qp2, true)
}
END_SECTION

START_SECTION(([QcMLFile::QualityParameter] String toXMLString(UInt indentation_level) const ))
{
	NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



