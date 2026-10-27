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
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <OpenMS/KERNEL/MSExperiment.h>
#include <OpenMS/METADATA/PeptideIdentificationList.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
///////////////////////////

#include <algorithm>
#include <fstream>
#include <sstream>
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
  // store() must tag a set member's quality parameter with the same CV
  // accession load() looks for inside a setQuality element ("MS:1000577",
  // the same one already used, consistently, for a run's own name) -- not
  // "QC:0000005", which load() never recognizes as a member at all. [EXTRA]
  QcMLFile qcml;
  qcml.registerRun("run1_id", "run1_name");

  QcMLFile::QualityParameter run_name_qp;
  run_name_qp.name = "mzML file";
  run_name_qp.id = "run1_run_name";
  run_name_qp.cvRef = "MS";
  run_name_qp.cvAcc = "MS:1000577";
  run_name_qp.value = "run1_name";
  qcml.addRunQualityParameter("run1_id", run_name_qp);

  std::set<std::string> members;
  members.insert("run1_id");
  qcml.registerSet("set1_id", "set1_name", members);

  std::string filename;
  NEW_TMP_FILE(filename);
  qcml.store(filename);

  ifstream is(filename.c_str());
  std::stringstream contents;
  contents << is.rdbuf();
  std::string written = contents.str();

  // The embedded qcML report stylesheet carries its own unrelated
  // "MS:1000577" occurrences, so check the member entry's exact
  // serialized form rather than a whole-file accession count.
  bool member_tagged_correctly =
      written.find("ID=\"run1_id\" cvRef=\"MS\" accession=\"MS:1000577\"") != std::string::npos;
  TEST_EQUAL(member_tagged_correctly, true)
  TEST_EQUAL(written.find("QC:0000005") == std::string::npos, true)
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
  // CPP-155: same unitRef/unitAcc vs unitCvRef/unitAccession mismatch as QualityParameter::toXMLString.
  QcMLFile::Attachment at;
  at.unitRef = "UO";
  at.unitAcc = "UO:0000169";
  at.binary = "somebinarydata"; // toXMLString returns "" for an attachment with neither binary nor a table

  std::string xml = at.toXMLString(0);
  TEST_EQUAL(xml.find("unitCvRef=\"UO\"") != std::string::npos, true)
  TEST_EQUAL(xml.find("unitAccession=\"UO:0000169\"") != std::string::npos, true)
  TEST_EQUAL(xml.find("unitRef=\"") != std::string::npos, false)
  TEST_EQUAL(xml.find("unitAcc=\"") != std::string::npos, false)

  // CPP-156: the table-row writer computed a space-substituted copy_row but then concatenated
  // the original, unmodified row (*it) instead, so a cell containing a space (e.g. "hello world")
  // is written with the space intact and becomes two tokens when the row is split on whitespace
  // while reloading.
  QcMLFile::Attachment table_at;
  table_at.colTypes.push_back("col1");
  table_at.colTypes.push_back("col2");
  std::vector<std::string> row;
  row.push_back("hello world");
  row.push_back("x");
  table_at.tableRows.push_back(row);

  std::string table_xml = table_at.toXMLString(0);
  TEST_EQUAL(table_xml.find("hello_world") != std::string::npos, true)
  TEST_EQUAL(table_xml.find("hello world") != std::string::npos, false)
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
  // CPP-155: the writer emitted unit information under unitRef/unitAcc, but the reader (and the
  // qcML schema) look for unitCvRef/unitAccession, so every unit was lost on a store/load round
  // trip.
  QcMLFile::QualityParameter qp;
  qp.unitRef = "UO";
  qp.unitAcc = "UO:0000169";

  std::string xml = qp.toXMLString(0);
  TEST_EQUAL(xml.find("unitCvRef=\"UO\"") != std::string::npos, true)
  TEST_EQUAL(xml.find("unitAccession=\"UO:0000169\"") != std::string::npos, true)
  TEST_EQUAL(xml.find("unitRef=\"") != std::string::npos, false)
  TEST_EQUAL(xml.find("unitAcc=\"") != std::string::npos, false)
}
END_SECTION

START_SECTION(([EXTRA] collectQCData reports TIC and RIC slumps as percentages))
{
  // 4 MS1 spectra (1 with a summed intensity below 10k) and 4 MS2 spectra, plus a TIC with 5 points (2 below 10k)
  MSExperiment exp;
  for (Size i = 0; i < 8; ++i)
  {
    MSSpectrum spec;
    spec.setRT(10.0 * i);
    spec.setMSLevel(i % 2 == 0 ? 1 : 2);
    if (spec.getMSLevel() == 2)
    {
      Precursor prec;
      prec.setMZ(500.0);
      prec.setCharge(2);
      spec.setPrecursors({prec});
    }
    const double intensity = (i == 0) ? 1000.0 : 50000.0;
    spec.push_back(Peak1D(400.0, intensity));
    spec.push_back(Peak1D(600.0, intensity));
    exp.addSpectrum(spec);
  }
  MSChromatogram tic;
  tic.setChromatogramType(ChromatogramSettings::ChromatogramType::TOTAL_ION_CURRENT_CHROMATOGRAM);
  for (Size i = 0; i < 5; ++i)
  {
    tic.push_back(ChromatogramPeak(10.0 * i, i < 2 ? 5000.0 : 50000.0));
  }
  exp.addChromatogram(tic);
  exp.updateRanges();

  std::vector<ProteinIdentification> prot_ids;
  PeptideIdentificationList pep_ids;
  QcMLFile qcml;
  qcml.collectQCData(prot_ids, pep_ids, FeatureMap(), ConsensusMap(), "slumps.mzML", false, exp);

  TEST_REAL_SIMILAR(std::stod(qcml.exportQP("slumps", "QC:0000057")), 25.0) // 1 of 4 MS1 spectra
  TEST_REAL_SIMILAR(std::stod(qcml.exportQP("slumps", "QC:0000023")), 40.0) // 2 of 5 TIC points
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



