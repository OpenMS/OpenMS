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
///////////////////////////

using namespace OpenMS;
using namespace std;

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

START_SECTION((static ExperimentalDesign load(const String &tsv_file, bool require_spectra_files)))
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


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
