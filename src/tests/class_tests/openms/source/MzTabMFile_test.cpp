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
#include <OpenMS/CHEMISTRY/EmpiricalFormula.h>
#include <OpenMS/FORMAT/FeatureXMLFile.h>
#include <OpenMS/FORMAT/MzTabMFile.h>
#include <OpenMS/FORMAT/FileHandler.h>
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

START_SECTION(void store(const std::string& filename, MzTabM& mztab_m))
    {
      FeatureMap feature_map;
      MzTabM mztabm;

      // AccurateMassSearch result (ID format) with owning identification data
      FileHandler().loadFeatures(OPENMS_GET_TEST_DATA_PATH("MzTabMFile_input_1.featureparquet"), feature_map, {FileTypes::FEATUREPARQUET});

      mztabm = MzTabM::exportFeatureMapToMzTabM(feature_map);

      // Evidence scores and ion m/z must survive the owning-model conversion.
      const auto& rows = mztabm.getMSmallMoleculeEvidenceSectionRows();
      TEST_EQUAL(rows.size(), 312);
      TEST_EQUAL(mztabm.getMSmallMoleculeFeatureSectionRows().size(), 83);
      TEST_EQUAL(mztabm.getMSmallMoleculeSectionRows().size(), 83);
      TEST_EQUAL(rows.front().id_confidence_measure.size(), 2);
      TEST_REAL_SIMILAR(rows.front().id_confidence_measure.at(1).get(), 0.225411404792002);
      TEST_REAL_SIMILAR(rows.front().id_confidence_measure.at(2).get(), 2.661798863812237e-5);
      TEST_EQUAL(rows.front().charge.get(), 1);
      const auto& ref = *feature_map.front().getIDMatches().begin();
      const auto& match = *feature_map.getIdentificationData().findRunByUuid(ref.run_uuid)->findMatch(ref.match);
      TEST_EQUAL(match.details->adduct.has_value(), true);
      TEST_REAL_SIMILAR(rows.front().calc_mass_to_charge.get(), match.details->adduct->getMZ(EmpiricalFormula(*match.details->formula).getMonoWeight()));

      std::string mztabm_tmpfile;
      NEW_TMP_FILE(mztabm_tmpfile);
      MzTabMFile().store(mztabm_tmpfile, mztabm);

      TEST_FILE_SIMILAR(mztabm_tmpfile.c_str(), OPENMS_GET_TEST_DATA_PATH("MzTabMFile_output_1.mztab"));
    }
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST