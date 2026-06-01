// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Oliver Alka$
// $Authors: Oliver Alka$
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/FORMAT/MzTabMFile.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/OMSFile.h>
#include <OpenMS/FORMAT/TextFile.h>
#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(MzTabMFile, "$Id$")
/////////////////////////////////////////////////////////////
MzTabMFile* ptr = nullptr;
MzTabMFile* null_ptr = nullptr;

START_SECTION(MzTabMFile())
    {
      ptr = new MzTabMFile();
      TEST_NOT_EQUAL(ptr, null_ptr)
    }
END_SECTION

START_SECTION(~MzTabFile())
    {
      delete ptr;
    }
END_SECTION

START_SECTION(void store(const String& filename, MzTabM& mztab_m))
    {
      FeatureMap feature_map;
      MzTabM mztabm;

      OMSFile().load(OPENMS_GET_TEST_DATA_PATH("MzTabMFile_input_1.oms"), feature_map);

      mztabm = MzTabM::exportFeatureMapToMzTabM(feature_map);

      String mztabm_tmpfile;
      NEW_TMP_FILE(mztabm_tmpfile);
      MzTabMFile().store(mztabm_tmpfile, mztabm);

      TEST_FILE_SIMILAR(mztabm_tmpfile.c_str(), OPENMS_GET_TEST_DATA_PATH("MzTabMFile_output_1.mztab"));
    }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST