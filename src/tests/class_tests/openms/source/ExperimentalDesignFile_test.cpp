// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// 
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg$
// $Authors: Timo Sachsenberg$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/METADATA/ExperimentalDesign.h>
#include <OpenMS/FORMAT/ExperimentalDesignFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/SYSTEM/TempFiles.h>
#include <fstream>
///////////////////////////

using namespace OpenMS;
using namespace std;

// Writes @p content to a temporary design file and returns its path
std::string writeDesign(const std::string& content)
{
  const std::string path = TempFiles::getTemporaryFile();
  std::ofstream os(path.c_str());
  os << content;
  return path;
}

// Returns the message of the Exception::ParseError that loading @p design_file throws; empty if it loads
std::string parseErrorMessage(const std::string& design_file)
{
  try
  {
    ExperimentalDesignFile::load(design_file, false);
  }
  catch (const Exception::ParseError& e)
  {
    return e.what();
  }
  return "";
}

START_TEST(ExperimentalDesignFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

ExperimentalDesignFile* ptr = 0;
ExperimentalDesignFile* null_ptr = 0;
START_SECTION(ExperimentalDesignFile())
{
  ptr = new ExperimentalDesignFile();
  TEST_NOT_EQUAL(ptr, null_ptr)
}
END_SECTION

START_SECTION(~ExperimentalDesignFile())
{
  delete ptr;
}
END_SECTION

START_SECTION((static ExperimentalDesign load(const std::string &tsv_file, bool require_spectra_files)))
{
ExperimentalDesign design = ExperimentalDesignFile::load(
  OPENMS_GET_TEST_DATA_PATH("ExperimentalDesign_input_1.tsv"), false);
  // tested in ExperimentalDesign_test
}
END_SECTION

START_SECTION((static ExperimentalDesign load(const TextFile&, bool, String) rejects multiplex one-table design without Sample column))
{
  TextFile tf;
  tf.addLine("Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tMSstats_Condition");
  tf.addLine("1\t1\tmix_a.mzML\t1\tA");
  tf.addLine("1\t1\tmix_a.mzML\t2\tA");

  TEST_EXCEPTION(Exception::ParseError, ExperimentalDesignFile::load(tf, false, "inline_multiplex_no_sample.tsv"));
}
END_SECTION

START_SECTION(([EXTRA] load rejects a sample-section row of a two-table design that ends in an empty cell))
{
  // Lines are whitespace-trimmed before they are split into cells, so a row whose LAST cell is
  // empty comes out one cell short. Such a sample row used to be stored as it was, and the later
  // getFactorValue() for the missing column (MSstatsConverter asks for MSstats_BioReplicate) read
  // the row past its end: undefined behaviour, in practice a crash or a garbage string.
  const std::string file_section = "Fraction_Group\tFraction\tSpectra_Filepath\tLabel\tSample\n"
                                   "1\t1\ta.mzML\t1\tS1\n"
                                   "\n"
                                   "Sample\tMSstats_Condition\tMSstats_BioReplicate\n";

  // line 5: empty MSstats_BioReplicate, i.e. 2 cells for a 3-column header
  const std::string short_row = writeDesign(file_section + "S1\tA\t\n");
  TEST_EXCEPTION(Exception::ParseError, ExperimentalDesignFile::load(short_row, false))
  // ... and the error names the line and the expected versus the actual number of cells
  const std::string msg = parseErrorMessage(short_row);
  TEST_TRUE(msg.find("in line 5 (") != std::string::npos)
  TEST_TRUE(msg.find("expected 3 ") != std::string::npos)
  TEST_TRUE(msg.find("found 2.") != std::string::npos)

  // positive control: with the last cell filled the same design loads and answers as before
  const ExperimentalDesign design = ExperimentalDesignFile::load(writeDesign(file_section + "S1\tA\t1\n"), false);
  TEST_EQUAL(design.getMSFileSection().size(), 1)
  TEST_STRING_EQUAL(design.getMSFileSection()[0].sample_name, "S1")
  TEST_EQUAL(design.getMSFileSection()[0].sample, 0)
  TEST_EQUAL(design.getMSFileSection()[0].label, 1)
  TEST_STRING_EQUAL(design.getSampleSection().getFactorValue("S1", "MSstats_Condition"), "A")
  TEST_STRING_EQUAL(design.getSampleSection().getFactorValue("S1", "MSstats_BioReplicate"), "1")
  TEST_STRING_EQUAL(design.getSampleSection().getFactorValue(0u, "MSstats_BioReplicate"), "1")
}
END_SECTION

START_SECTION(([EXTRA] load rejects a row of a label-free one-table design that ends in an empty cell))
{
  // Same trimming effect in the one-table layout. Here the parser appends the implicit Label
  // column ("1") to every row and then used to read cells[Label] BEFORE checking the width of
  // the row: for this 5-column header the trimmed 4-cell row grows to 5 cells, and cells[5] is
  // past its end.
  const std::string header = "Fraction_Group\tFraction\tSpectra_Filepath\tSample\tCondition\n";

  // line 2: empty Condition, i.e. 4 cells for a 5-column header
  const std::string short_row = writeDesign(header + "1\t1\ta.mzML\tS1\t\n");
  TEST_EXCEPTION(Exception::ParseError, ExperimentalDesignFile::load(short_row, false))
  const std::string msg = parseErrorMessage(short_row);
  TEST_TRUE(msg.find("in line 2 (") != std::string::npos)
  TEST_TRUE(msg.find("expected 5 ") != std::string::npos)
  TEST_TRUE(msg.find("found 4.") != std::string::npos)

  // positive control: the well-formed row still yields the implicit label, the sample and its factor
  const ExperimentalDesign design = ExperimentalDesignFile::load(writeDesign(header + "1\t1\ta.mzML\tS1\tA\n"), false);
  TEST_EQUAL(design.getMSFileSection().size(), 1)
  TEST_EQUAL(design.getMSFileSection()[0].label, 1)
  TEST_EQUAL(design.getMSFileSection()[0].fraction_group, 1)
  TEST_EQUAL(design.getMSFileSection()[0].fraction, 1)
  TEST_EQUAL(design.getMSFileSection()[0].sample, 0)
  TEST_STRING_EQUAL(design.getMSFileSection()[0].sample_name, "S1")
  TEST_EQUAL(design.getSampleSection().getFactors().size(), 2) // Sample and Condition
  TEST_STRING_EQUAL(design.getSampleSection().getFactorValue("S1", "Condition"), "A")
  TEST_STRING_EQUAL(design.getSampleSection().getFactorValue("S1", "Sample"), "S1")

  // a row that is too LONG is rejected as before, now before it is indexed
  TEST_EXCEPTION(Exception::ParseError, ExperimentalDesignFile::load(writeDesign(header + "1\t1\ta.mzML\tS1\tA\textra\n"), false))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
