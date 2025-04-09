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
  TEST_EQUAL(FileTypes::typeToName(FileTypes::UNKNOWN), "unknown");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::DTA), "dta");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::DTA2D), "dta2d");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::MZDATA), "mzData");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::MZXML), "mzXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::MZML), "mzML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::FEATUREXML), "featureXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::IDXML), "idXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::CONSENSUSXML), "consensusXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::TRANSFORMATIONXML), "trafoXML");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::INI), "ini");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::TOPPAS), "toppas");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::PNG), "png");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::TXT), "txt");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::CSV), "csv");
  TEST_EQUAL(FileTypes::typeToName(FileTypes::MZTAB), "mzTab");

  // try them all, just to make sure they are all there
  for (int i = 0; i < (int)FileTypes::SIZE_OF_TYPE; ++i)
  {
    TEST_EQUAL(FileTypes::nameToType(FileTypes::typeToName(FileTypes::Type(i))), FileTypes::Type(i));
  }
}
END_SECTION

START_SECTION((static Type nameToType(const String& name)))
  TEST_EQUAL(FileTypes::typeToDescription(FileTypes::DTA2D), "dta2d raw data file");
  TEST_EQUAL(FileTypes::typeToDescription(FileTypes::UNKNOWN), "unknown file extension");
END_SECTION


START_SECTION((static Type nameToType(const String& name)))
{
  TEST_EQUAL(FileTypes::UNKNOWN, FileTypes::nameToType("unknown"));
  TEST_EQUAL(FileTypes::DTA, FileTypes::nameToType("dta"));
  TEST_EQUAL(FileTypes::DTA2D, FileTypes::nameToType("dta2d"));
  TEST_EQUAL(FileTypes::MZDATA, FileTypes::nameToType("mzData"));
  TEST_EQUAL(FileTypes::MZXML, FileTypes::nameToType("mzXML"));
  TEST_EQUAL(FileTypes::FEATUREXML, FileTypes::nameToType("featureXML"));
  TEST_EQUAL(FileTypes::IDXML, FileTypes::nameToType("idXmL")); // case-insensitivity
  TEST_EQUAL(FileTypes::CONSENSUSXML, FileTypes::nameToType("consensusXML"));
  TEST_EQUAL(FileTypes::MGF, FileTypes::nameToType("mgf"));
  TEST_EQUAL(FileTypes::INI, FileTypes::nameToType("ini"));
  TEST_EQUAL(FileTypes::TOPPAS, FileTypes::nameToType("toppas"));
  TEST_EQUAL(FileTypes::TRANSFORMATIONXML, FileTypes::nameToType("trafoXML"));
  TEST_EQUAL(FileTypes::MZML, FileTypes::nameToType("mzML"));
  TEST_EQUAL(FileTypes::MS2, FileTypes::nameToType("ms2"));
  TEST_EQUAL(FileTypes::PEPXML, FileTypes::nameToType("pepXML"));
  TEST_EQUAL(FileTypes::PROTXML, FileTypes::nameToType("protXML"));
  TEST_EQUAL(FileTypes::MZIDENTML, FileTypes::nameToType("mzid"));
  TEST_EQUAL(FileTypes::GELML, FileTypes::nameToType("gelML"));
  TEST_EQUAL(FileTypes::TRAML, FileTypes::nameToType("traML"));
  TEST_EQUAL(FileTypes::MSP, FileTypes::nameToType("msp"));
  TEST_EQUAL(FileTypes::OMSSAXML, FileTypes::nameToType("omssaXML"));
  TEST_EQUAL(FileTypes::PNG, FileTypes::nameToType("png"));
  TEST_EQUAL(FileTypes::XMASS, FileTypes::nameToType("fid"));
  TEST_EQUAL(FileTypes::TSV, FileTypes::nameToType("tsv"));
  TEST_EQUAL(FileTypes::PEPLIST, FileTypes::nameToType("peplist"));
  TEST_EQUAL(FileTypes::HARDKLOER, FileTypes::nameToType("hardkloer"));
  TEST_EQUAL(FileTypes::KROENIK, FileTypes::nameToType("kroenik"));
  TEST_EQUAL(FileTypes::FASTA, FileTypes::nameToType("fasta"));
  TEST_EQUAL(FileTypes::EDTA, FileTypes::nameToType("edta"));
  TEST_EQUAL(FileTypes::CSV, FileTypes::nameToType("csv"));
  TEST_EQUAL(FileTypes::TXT, FileTypes::nameToType("txt"));

  TEST_EQUAL(FileTypes::UNKNOWN, FileTypes::nameToType("somethingunknown"));
}
END_SECTION

START_SECTION([EXTRA] FileTypes::FileTypeList)
  FileTypeList list({ FileTypes::MZML, FileTypes::BZ2 });
  TEST_EQUAL(list.contains(FileTypes::MZML), true);
  TEST_EQUAL(list.contains(FileTypes::BZ2), true);
  TEST_EQUAL(list.contains(FileTypes::MZDATA), false);

  TEST_EQUAL(list.toFileDialogFilter(FilterLayout::BOTH, true), "all readable files (*.mzML *.bz2);;mzML raw data file (*.mzML);;bzip2 compressed file (*.bz2);;all files (*)")
  TEST_EQUAL(list.toFileDialogFilter(FilterLayout::COMPACT, true), "all readable files (*.mzML *.bz2);;all files (*)")
  TEST_EQUAL(list.toFileDialogFilter(FilterLayout::ONE_BY_ONE, true), "mzML raw data file (*.mzML);;bzip2 compressed file (*.bz2);;all files (*)")
  TEST_EQUAL(list.toFileDialogFilter(FilterLayout::BOTH, false), "all readable files (*.mzML *.bz2);;mzML raw data file (*.mzML);;bzip2 compressed file (*.bz2)")

  // testing Type FileTypeList::fromFileDialogFilter(const String& filter, const Type fallback = Type::UNKNOWN) const
  TEST_EQUAL(list.fromFileDialogFilter("all readable files (*.mzML *.bz2)"), FileTypes::UNKNOWN);
  TEST_EQUAL(list.fromFileDialogFilter("all files (*)"), FileTypes::UNKNOWN);
  TEST_EQUAL(list.fromFileDialogFilter("mzML raw data file (*.mzML)"), FileTypes::MZML);
  TEST_EQUAL(list.fromFileDialogFilter("bzip2 compressed file (*.bz2)"), FileTypes::BZ2);
  TEST_EXCEPTION(Exception::ElementNotFound, list.fromFileDialogFilter("not a valid filter"));
  // with default
  TEST_EQUAL(list.fromFileDialogFilter("all readable files (*.mzML *.bz2)", FileTypes::CONSENSUSXML), FileTypes::CONSENSUSXML);
  TEST_EQUAL(list.fromFileDialogFilter("all files (*)", FileTypes::CONSENSUSXML), FileTypes::CONSENSUSXML);
  TEST_EQUAL(list.fromFileDialogFilter("mzML raw data file (*.mzML)", FileTypes::CONSENSUSXML), FileTypes::MZML);
  TEST_EQUAL(list.fromFileDialogFilter("bzip2 compressed file (*.bz2)", FileTypes::CONSENSUSXML), FileTypes::BZ2);
  TEST_EXCEPTION(Exception::ElementNotFound, list.fromFileDialogFilter("not a valid filter", FileTypes::CONSENSUSXML));

  END_SECTION

  START_SECTION(static FileTypes::FileTypeList typesWithProperties(const std::vector<FileProperties>& features))
  {
    std::vector<FileTypes::FileProperties> f;
    f.push_back(FileTypes::FileProperties::READABLE);
    FileTypeList g = FileTypeList::typesWithProperties(f);
    TEST_EQUAL(g.getTypes().size(), 37);
    // Test that empty filter returns the full list
    TEST_EQUAL(FileTypeList::typesWithProperties({}).size(), 60);
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
