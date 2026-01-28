// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Stephan Aiche, Andreas Bertsch, Marc Sturm, Chris Bielow $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/FileHandler.h>
#include <OpenMS/FORMAT/FileTypes.h>


///////////////////////////

START_TEST(FileHandler, "Id")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

START_SECTION((static String typeToName(Type type)))
{
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::UNKNOWN), "unknown");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::DTA), "dta");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::DTA2D), "dta2d");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::MZDATA), "mzData");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::MZXML), "mzXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::MZML), "mzML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::FEATUREXML), "featureXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::IDXML), "idXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::CONSENSUSXML), "consensusXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::TRANSFORMATIONXML), "trafoXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::INI), "ini");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::TOPPAS), "toppas");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::PNG), "png");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::TXT), "txt");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::CSV), "csv");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::Type::MZTAB), "mzTab");

  // try them all, just to make sure they are all there
  for (int i = 0; i < (int)FileTypes::Type::SIZE_OF_TYPE; ++i)
  {
    TEST_EQUAL(FileTypes::nameToType(FileTypes::typeToName(FileTypes::Type(i))), FileTypes::Type(i));
  }
}
END_SECTION

START_SECTION((static Type nameToType(const String& name)))
  TEST_EQUAL(FileTypes::typeToDescription(FileTypes::Type::DTA2D), "dta2d raw data file");
  TEST_EQUAL(FileTypes::typeToDescription(FileTypes::Type::UNKNOWN), "unknown file extension");
END_SECTION


START_SECTION((static Type nameToType(const String& name)))
{
  TEST_EQUAL(FileTypes::Type::UNKNOWN, FileTypes::nameToType("unknown"));
  TEST_EQUAL(FileTypes::Type::DTA, FileTypes::nameToType("dta"));
  TEST_EQUAL(FileTypes::Type::DTA2D, FileTypes::nameToType("dta2d"));
  TEST_EQUAL(FileTypes::Type::MZDATA, FileTypes::nameToType("mzData"));
  TEST_EQUAL(FileTypes::Type::MZXML, FileTypes::nameToType("mzXML"));
  TEST_EQUAL(FileTypes::Type::FEATUREXML, FileTypes::nameToType("featureXML"));
  TEST_EQUAL(FileTypes::Type::IDXML, FileTypes::nameToType("idXmL")); // case-insensitivity
  TEST_EQUAL(FileTypes::Type::CONSENSUSXML, FileTypes::nameToType("consensusXML"));
  TEST_EQUAL(FileTypes::Type::MGF, FileTypes::nameToType("mgf"));
  TEST_EQUAL(FileTypes::Type::INI, FileTypes::nameToType("ini"));
  TEST_EQUAL(FileTypes::Type::TOPPAS, FileTypes::nameToType("toppas"));
  TEST_EQUAL(FileTypes::Type::TRANSFORMATIONXML, FileTypes::nameToType("trafoXML"));
  TEST_EQUAL(FileTypes::Type::MZML, FileTypes::nameToType("mzML"));
  TEST_EQUAL(FileTypes::Type::MS2, FileTypes::nameToType("ms2"));
  TEST_EQUAL(FileTypes::Type::PEPXML, FileTypes::nameToType("pepXML"));
  TEST_EQUAL(FileTypes::Type::PROTXML, FileTypes::nameToType("protXML"));
  TEST_EQUAL(FileTypes::Type::MZIDENTML, FileTypes::nameToType("mzid"));
  TEST_EQUAL(FileTypes::Type::GELML, FileTypes::nameToType("gelML"));
  TEST_EQUAL(FileTypes::Type::TRAML, FileTypes::nameToType("traML"));
  TEST_EQUAL(FileTypes::Type::MSP, FileTypes::nameToType("msp"));
  TEST_EQUAL(FileTypes::Type::OMSSAXML, FileTypes::nameToType("omssaXML"));
  TEST_EQUAL(FileTypes::Type::PNG, FileTypes::nameToType("png"));
  TEST_EQUAL(FileTypes::Type::XMASS, FileTypes::nameToType("fid"));
  TEST_EQUAL(FileTypes::Type::TSV, FileTypes::nameToType("tsv"));
  TEST_EQUAL(FileTypes::Type::PEPLIST, FileTypes::nameToType("peplist"));
  TEST_EQUAL(FileTypes::Type::HARDKLOER, FileTypes::nameToType("hardkloer"));
  TEST_EQUAL(FileTypes::Type::KROENIK, FileTypes::nameToType("kroenik"));
  TEST_EQUAL(FileTypes::Type::FASTA, FileTypes::nameToType("fasta"));
  TEST_EQUAL(FileTypes::Type::EDTA, FileTypes::nameToType("edta"));
  TEST_EQUAL(FileTypes::Type::CSV, FileTypes::nameToType("csv"));
  TEST_EQUAL(FileTypes::Type::TXT, FileTypes::nameToType("txt"));
  TEST_EQUAL(FileTypes::Type::PARQUET, FileTypes::nameToType("parquet"));
  TEST_EQUAL(FileTypes::Type::PARQUET, FileTypes::nameToType("pqt")); // Test alternate extension

  TEST_EQUAL(FileTypes::Type::UNKNOWN, FileTypes::nameToType("somethingunknown"));
}
END_SECTION

START_SECTION([EXTRA] FileTypes::FileTypeList)
  FileTypeList list({ FileTypes::Type::MZML, FileTypes::Type::BZ2 });
  TEST_EQUAL(list.contains(FileTypes::Type::MZML), true);
  TEST_EQUAL(list.contains(FileTypes::Type::BZ2), true);
  TEST_EQUAL(list.contains(FileTypes::Type::MZDATA), false);

  TEST_EQUAL(list.toFileDialogFilter(FilterLayout::BOTH, true), "all readable files (*.mzML *.bz2);;mzML raw data file (*.mzML);;bzip2 compressed file (*.bz2);;all files (*)")
  TEST_EQUAL(list.toFileDialogFilter(FilterLayout::COMPACT, true), "all readable files (*.mzML *.bz2);;all files (*)")
  TEST_EQUAL(list.toFileDialogFilter(FilterLayout::ONE_BY_ONE, true), "mzML raw data file (*.mzML);;bzip2 compressed file (*.bz2);;all files (*)")
  TEST_EQUAL(list.toFileDialogFilter(FilterLayout::BOTH, false), "all readable files (*.mzML *.bz2);;mzML raw data file (*.mzML);;bzip2 compressed file (*.bz2)")

  // testing Type FileTypeList::fromFileDialogFilter(const String& filter, const Type fallback = Type::UNKNOWN) const
  TEST_EQUAL(list.fromFileDialogFilter("all readable files (*.mzML *.bz2)"), FileTypes::Type::UNKNOWN);
  TEST_EQUAL(list.fromFileDialogFilter("all files (*)"), FileTypes::Type::UNKNOWN);
  TEST_EQUAL(list.fromFileDialogFilter("mzML raw data file (*.mzML)"), FileTypes::Type::MZML);
  TEST_EQUAL(list.fromFileDialogFilter("bzip2 compressed file (*.bz2)"), FileTypes::Type::BZ2);
  TEST_EXCEPTION(Exception::ElementNotFound, list.fromFileDialogFilter("not a valid filter"));
  // with default
  TEST_EQUAL(list.fromFileDialogFilter("all readable files (*.mzML *.bz2)", FileTypes::Type::CONSENSUSXML), FileTypes::Type::CONSENSUSXML);
  TEST_EQUAL(list.fromFileDialogFilter("all files (*)", FileTypes::Type::CONSENSUSXML), FileTypes::Type::CONSENSUSXML);
  TEST_EQUAL(list.fromFileDialogFilter("mzML raw data file (*.mzML)", FileTypes::Type::CONSENSUSXML), FileTypes::Type::MZML);
  TEST_EQUAL(list.fromFileDialogFilter("bzip2 compressed file (*.bz2)", FileTypes::Type::CONSENSUSXML), FileTypes::Type::BZ2);
  TEST_EXCEPTION(Exception::ElementNotFound, list.fromFileDialogFilter("not a valid filter", FileTypes::Type::CONSENSUSXML));

  END_SECTION

  START_SECTION(static FileTypes::FileTypeList typesWithProperties(const std::vector<FileProperties>& features))
  {
    std::vector<FileTypes::FileProperties> f;
    f.push_back(FileTypes::FileProperties::READABLE);
    FileTypeList g = FileTypeList::typesWithProperties(f);
    TEST_EQUAL(g.getTypes().size(), 38);
    // Test that empty filter returns the full list
    TEST_EQUAL(FileTypeList::typesWithProperties({}).size(), 61);
    // Test that the full list is equal to the list of known file types
    TEST_EQUAL(FileTypeList::typesWithProperties({}).size(),static_cast<size_t>(FileTypes::Type::SIZE_OF_TYPE));
    // Check that we don't have duplicate Types in our type_with_annotation__
    vector<FileTypes::Type> vec = FileTypeList::typesWithProperties({});
    sort(vec.begin(),vec.end());
    auto it = std::unique(vec.begin(), vec.end());
    TEST_TRUE(it ==vec.end());
  }
  END_SECTION

END_TEST
