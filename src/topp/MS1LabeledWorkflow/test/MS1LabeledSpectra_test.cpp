// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Timo Sachsenberg $
// --------------------------------------------------------------------------

#include "MS1LabeledSpectra.h"

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/FORMAT/MzMLFile.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <limits>

using namespace OpenMS;

START_TEST(MS1LabeledSpectra, "$Id$")

START_SECTION(load mzML retains all scan metadata and only MS1 peaks)
{
  MSExperiment original;
  for (int i = 0; i < 4; ++i)
  {
    MSSpectrum spectrum;
    spectrum.setMSLevel(i % 2 + 1);
    spectrum.setRT(100.0 + i * 0.1);
    spectrum.setNativeID("controllerType=0 controllerNumber=1 scan=" + std::to_string(i + 1));
    spectrum.setType(SpectrumSettings::SpectrumType::CENTROID);
    spectrum.setDriftTime(-45.0 - 20.0 * (i / 2));
    spectrum.setDriftTimeUnit(DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE);
    spectrum.getAcquisitionInfo().push_back(Acquisition());
    spectrum.emplace_back(500.0, 100.0 + i);
    original.addSpectrum(spectrum);
  }
  std::string filename;
  NEW_TMP_FILE(filename)
  MzMLFile().store(filename, original);
  MSExperiment metadata, ms1;
  MS1LabeledSpectra::load(filename, metadata, ms1, ProgressLogger::NONE);
  TEST_EQUAL(metadata.size(), 4)
  TEST_EQUAL(ms1.size(), 2)
  for (Size i = 0; i < metadata.size(); ++i)
  {
    TEST_TRUE(metadata[i].empty())
    TEST_EQUAL(metadata[i].getNativeID(), original[i].getNativeID())
    TEST_REAL_SIMILAR(metadata[i].getRT(), original[i].getRT())
    TEST_REAL_SIMILAR(metadata[i].getDriftTime(), original[i].getDriftTime())
    TEST_EQUAL(metadata[i].getDriftTimeUnit(), DriftTimeUnit::FAIMS_COMPENSATION_VOLTAGE)
  }
  TEST_EQUAL(ms1[0].getMSLevel(), 1)
  TEST_EQUAL(ms1[1].getMSLevel(), 1)
  TEST_REAL_SIMILAR(ms1[1][0].getIntensity(), 102.0)
  // Reusing the outputs must replace rather than append spectra.
  MS1LabeledSpectra::load(filename, metadata, ms1, ProgressLogger::NONE);
  TEST_EQUAL(metadata.size(), 4)
  TEST_EQUAL(ms1.size(), 2)
}
END_SECTION

START_SECTION(repair missing references from MS2 metadata without reopening a file)
{
  MSExperiment metadata;
  MSSpectrum ms1, ms2;
  ms1.setRT(100.0);
  ms1.setMSLevel(1);
  ms1.setNativeID("controllerType=0 controllerNumber=1 scan=1");
  ms2.setRT(100.004);
  ms2.setMSLevel(2);
  ms2.setNativeID("controllerType=0 controllerNumber=1 scan=2");
  metadata.addSpectrum(ms1);
  metadata.addSpectrum(ms2);
  PeptideIdentificationList ids(3);
  for (auto& id : ids)
  {
    id.setRT(100.0);
  }
  ids[1].setSpectrumReference("");
  ids[2].setSpectrumReference("scan=existing");
  MS1LabeledSpectra::addMissingSpectrumReferences(metadata, ids);
  TEST_EQUAL(ids[0].getSpectrumReference(), ms2.getNativeID())
  TEST_EQUAL(ids[1].getSpectrumReference(), ms2.getNativeID())
  TEST_EQUAL(ids[2].getSpectrumReference(), "scan=existing")

  ids.resize(1);
  ids[0].setSpectrumReference("");
  ids[0].setRT(200.0);
  TEST_EXCEPTION(Exception::MissingInformation, MS1LabeledSpectra::addMissingSpectrumReferences(metadata, ids))
  ids[0].setRT(std::numeric_limits<double>::quiet_NaN());
  TEST_EXCEPTION(Exception::MissingInformation, MS1LabeledSpectra::addMissingSpectrumReferences(metadata, ids))
}
END_SECTION

#if defined(WITH_THERMO_RAW) && defined(THERMO_RAW_TEST_DATA)
START_SECTION(load native Thermo RAW and repair references using its MS2 metadata)
{
  // Shared opt-in PNNL Angiotensin fixture: 87 MS1 scans and over 1000 scans in total.
  MSExperiment metadata, ms1;
  MS1LabeledSpectra::load(THERMO_RAW_TEST_DATA, metadata, ms1, ProgressLogger::NONE);
  TEST_TRUE(metadata.size() > 1000)
  TEST_EQUAL(ms1.size(), 87)
  TEST_FALSE(ms1.getSourceFiles().empty())
  Size ms1_index = 0;
  bool repaired = false;
  for (const auto& spectrum : metadata)
  {
    TEST_TRUE(spectrum.empty())
    TEST_FALSE(spectrum.getNativeID().empty())
    if (spectrum.getMSLevel() == 1)
    {
      TEST_EQUAL(ms1[ms1_index].getMSLevel(), 1)
      TEST_EQUAL(ms1[ms1_index].getNativeID(), spectrum.getNativeID())
      TEST_REAL_SIMILAR(ms1[ms1_index].getRT(), spectrum.getRT())
      TEST_FALSE(ms1[ms1_index].empty())
      ++ms1_index;
    }
    else if (spectrum.getMSLevel() == 2 && ! repaired)
    {
      PeptideIdentificationList ids(1);
      ids[0].setRT(spectrum.getRT());
      MS1LabeledSpectra::addMissingSpectrumReferences(metadata, ids);
      TEST_EQUAL(ids[0].getSpectrumReference(), spectrum.getNativeID())
      repaired = true;
    }
  }
  TEST_TRUE(repaired)
}
END_SECTION
#endif

END_TEST
