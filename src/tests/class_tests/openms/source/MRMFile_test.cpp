// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/MRMFile.h>
#include <OpenMS/METADATA/ExperimentalSettings.h>

using namespace OpenMS;
using namespace std;

START_TEST(MRMFile, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(static std::vector<::OpenSwath::SwathMap> loadMzML(...))
{
  std::shared_ptr<ExperimentalSettings> exp_meta;

  // Test with non-existent file - should throw exception
  TEST_EXCEPTION(Exception::FileNotFound, MRMFile::loadMzML("nonexistent.mzML", "", exp_meta))

  // Test with empty filename - should throw exception
  TEST_EXCEPTION(Exception::FileNotFound, MRMFile::loadMzML("", "", exp_meta))
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST