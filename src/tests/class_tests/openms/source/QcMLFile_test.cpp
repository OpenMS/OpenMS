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

#include <algorithm>
#include <vector>

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


START_SECTION((void registerRun(const std::string id, const std::string name)))
{
  QcMLFile qcmlfile1;
  qcmlfile1.registerRun("abc","somerun");
  TEST_EQUAL(qcmlfile1.existsRun("abc"), true)
  TEST_EQUAL(qcmlfile1.existsRun("somerun",true), true)
}
END_SECTION


START_SECTION((void registerSet(const std::string id, const std::string name, const std::set<std::string> &names)))
{
  QcMLFile qcmlfile1;
  std::set<std::string> n;
  n.insert("somerun1");
  n.insert("somerun2");
  qcmlfile1.registerSet("def","someset", n);
  TEST_EQUAL(qcmlfile1.existsSet("def"), true)
  TEST_EQUAL(qcmlfile1.existsSet("someset",true), true)
}
END_SECTION

START_SECTION((void addRunQualityParameter(std::string r, QualityParameter qp)))
NOT_TESTABLE
END_SECTION

START_SECTION((void addRunAttachment(std::string r, Attachment at)))
NOT_TESTABLE
END_SECTION

START_SECTION((void addSetQualityParameter(std::string r, QualityParameter qp)))
NOT_TESTABLE
END_SECTION

START_SECTION((void addSetAttachment(std::string r, Attachment at)))
NOT_TESTABLE
END_SECTION

START_SECTION((void removeAttachment(std::string r, std::vector<std::string> &ids, std::string at="")))
NOT_TESTABLE
END_SECTION

START_SECTION((void removeAttachment(std::string r, std::string at)))
NOT_TESTABLE
END_SECTION

START_SECTION((void removeAllAttachments(std::string at)))
NOT_TESTABLE
END_SECTION

START_SECTION((void removeQualityParameter(std::string r, std::vector<std::string> &ids)))
NOT_TESTABLE
END_SECTION

START_SECTION((void merge(const QcMLFile &addendum, std::string setname="")))
NOT_TESTABLE
END_SECTION

START_SECTION((void collectSetParameter(const std::string setname, const std::string qp, std::vector<std::string> &ret)))
NOT_TESTABLE
END_SECTION

START_SECTION((std::string exportAttachment(const std::string filename, const std::string qpname) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((std::string exportQP(const std::string filename, const std::string qpname) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((std::string exportQPs(const std::string filename, const StringList qpnames) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((std::string map2csv(const std::map<std::string, std::map<std::string, std::string> > &cvs_table, const std::string &separator) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((std::string exportIDstats(const std::string &filename) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((void getRunIDs(std::vector<std::string> &ids) const ))
{
  std::vector<std::string> ids, test;
  test.push_back("123");
  test.push_back("456");
  qcmlfile.getRunIDs(ids);
  TEST_EQUAL(StringList(ids),StringList(test))
}
END_SECTION

START_SECTION((void getRunNames(std::vector<std::string> &ids) const ))
{
  std::vector<std::string> ids, test;
  test.push_back("testrun1");
  test.push_back("testrun2");
  qcmlfile.getRunNames(ids);
  TEST_EQUAL(ids==test,true)
}
END_SECTION

START_SECTION((bool existsRun(const std::string filename, bool checkname=false) const ))
{
  QcMLFile qcmlfile1;
  qcmlfile.registerRun("abc","somerun");
  TEST_EQUAL(qcmlfile.existsRun("abc"), true)
  TEST_EQUAL(qcmlfile.existsRun("somerun",true), true)
}
END_SECTION

START_SECTION((bool existsSet(const std::string filename, bool checkname=false) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((void existsRunQualityParameter(const std::string filename, const std::string qpname, std::vector<std::string> &ids) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((void existsSetQualityParameter(const std::string filename, const std::string qpname, std::vector<std::string> &ids) const ))
NOT_TESTABLE
END_SECTION

START_SECTION((void store(const std::string &filename) const ))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION((void load(const std::string &filename)))
{
  // Regression: loading a second file into the same object must NOT leave stale
  // run/set names from the first file behind (the derived name->ID maps have to
  // be cleared in load() alongside the data maps).
  QcMLFile qcml;

  // --- first file A ---
  qcml.load(OPENMS_GET_TEST_DATA_PATH("QcMLFile_reload_A.qcML"));

  std::vector<std::string> run_names_a;
  qcml.getRunNames(run_names_a);
  TEST_EQUAL(run_names_a.size(), 1)
  TEST_TRUE(std::find(run_names_a.begin(), run_names_a.end(), "runAlpha") != run_names_a.end())
  TEST_TRUE(qcml.existsRun("runAlpha", true))
  TEST_TRUE(qcml.existsSet("setAlpha", true))

  // --- second file B (disjoint run/set names) into the SAME object ---
  qcml.load(OPENMS_GET_TEST_DATA_PATH("QcMLFile_reload_B.qcML"));

  std::vector<std::string> run_names_b;
  qcml.getRunNames(run_names_b);
  // only B's run name must remain, none of A's
  TEST_EQUAL(run_names_b.size(), 1)
  TEST_TRUE(std::find(run_names_b.begin(), run_names_b.end(), "runBeta") != run_names_b.end())
  TEST_TRUE(std::find(run_names_b.begin(), run_names_b.end(), "runAlpha") == run_names_b.end())

  // stale name lookups from file A must no longer resolve
  TEST_FALSE(qcml.existsRun("runAlpha", true))
  TEST_FALSE(qcml.existsSet("setAlpha", true))
  // B's names resolve
  TEST_TRUE(qcml.existsRun("runBeta", true))
  TEST_TRUE(qcml.existsSet("setBeta", true))
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

START_SECTION(([QcMLFile::Attachment] std::string toXMLString(UInt indentation_level) const ))
{
	NOT_TESTABLE
}
END_SECTION

START_SECTION(([QcMLFile::Attachment] std::string toCSVString(std::string separator) const ))
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

START_SECTION(([QcMLFile::QualityParameter] std::string toXMLString(UInt indentation_level) const ))
{
	NOT_TESTABLE
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



