// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/OpenSwathHelper.h>
#include <OpenMS/ANALYSIS/OPENSWATH/PeakMapExtractor.h>
#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/FORMAT/DATAACCESS/XIPMParquetConsumer.h>
#include <OpenMS/FORMAT/XIPMParquetFile.h>

using namespace OpenMS;
using namespace std;

namespace
{
  OpenSwath::LightTargetedExperiment makeExperiment_()
  {
    OpenSwath::LightTargetedExperiment exp;

    OpenSwath::LightCompound compound;
    compound.id = "pep1";
    compound.sequence = "PEPTIDE";
    compound.charge = 2;
    compound.rt = 100.0;
    compound.drift_time = 1.1;
    exp.compounds.push_back(compound);

    OpenSwath::LightTransition transition;
    transition.transition_name = "tr1";
    transition.peptide_ref = "pep1";
    transition.precursor_mz = 600.2;
    transition.product_mz = 500.2;
    transition.fragment_charge = 1;
    transition.fragment_nr = 7;
    transition.setFragmentType("y");
    transition.setDetectingTransition(true);
    exp.transitions.push_back(transition);

    return exp;
  }

  PeakMapExtractor::ExtractedPeakMap makeTransitionPeakMap_()
  {
    PeakMapExtractor::ExtractedPeakMap peak_map;
    peak_map.native_id = "tr1";
    peak_map.target_mz = 500.2;
    peak_map.target_rt = 100.0;
    peak_map.target_ion_mobility = 1.1;
    peak_map.rt_start = 95.0;
    peak_map.rt_end = 105.0;
    peak_map.mz = {500.19, 500.20};
    peak_map.rt = {100.0, 101.0};
    peak_map.ion_mobility = {1.05, 1.08};
    peak_map.intensity = {1000.0, 900.0};
    return peak_map;
  }

  PeakMapExtractor::ExtractedPeakMap makePrecursorPeakMap_()
  {
    PeakMapExtractor::ExtractedPeakMap peak_map;
    peak_map.native_id = OpenSwathHelper::computePrecursorId("pep1", 0);
    peak_map.target_mz = 600.2;
    peak_map.target_rt = 100.0;
    peak_map.target_ion_mobility = 1.1;
    peak_map.rt_start = 95.0;
    peak_map.rt_end = 105.0;
    peak_map.mz = {600.19};
    peak_map.rt = {100.0};
    peak_map.ion_mobility = {1.06};
    peak_map.intensity = {500.0};
    return peak_map;
  }

  String writeTestFile_(const UInt64 run_id, const String& source_file)
  {
    String tmp;
    NEW_TMP_FILE(tmp);
    const String out = tmp + ".xipm";

    const auto light_exp = makeExperiment_();
    XIPMParquetConsumer consumer(out, light_exp);
    consumer.consumePeakMap(makeTransitionPeakMap_(), run_id, source_file, 2);
    consumer.consumePeakMap(makePrecursorPeakMap_(), run_id, source_file, 1);
    consumer.finalize();
    return out;
  }
}

START_TEST(XIPMParquetFile, "$Id$")

START_SECTION(void load(std::vector<XIPMPeakMap>& output) const)
{
  const String file = writeTestFile_(7, "run1.mzML");
  XIPMParquetFile xipm(file);

  std::vector<XIPMParquetFile::XIPMPeakMap> peak_maps;
  xipm.load(peak_maps);
  TEST_EQUAL(peak_maps.size(), 2)
  TEST_EQUAL(peak_maps[0].mz.size(), peak_maps[0].rt.size())
  TEST_EQUAL(peak_maps[0].mz.size(), peak_maps[0].ion_mobility.size())
  TEST_EQUAL(peak_maps[0].mz.size(), peak_maps[0].intensity.size())
  TEST_EQUAL(peak_maps[0].has_target_rt, true)
  TEST_REAL_SIMILAR(peak_maps[0].target_rt, 100.0)
}
END_SECTION

START_SECTION(void getPeakMaps(...filters...) const)
{
  const String file = writeTestFile_(7, "run1.mzML");
  XIPMParquetFile xipm(file);

  std::vector<XIPMParquetFile::XIPMPeakMap> precursor_peak_maps;
  xipm.getPeakMaps(precursor_peak_maps, -1, -1, "", -1, -1, 1, 7, "precursor");
  TEST_EQUAL(precursor_peak_maps.size(), 1)
  TEST_STRING_EQUAL(precursor_peak_maps[0].peakmap_type, "precursor")

  std::vector<XIPMParquetFile::XIPMPeakMap> transition_peak_maps;
  xipm.getPeakMaps(transition_peak_maps, -1, 1, "", -1, -1, 2, 7, "transition");
  TEST_EQUAL(transition_peak_maps.size(), 1)
  TEST_STRING_EQUAL(transition_peak_maps[0].annotation, "y7^1")
  TEST_EQUAL(transition_peak_maps[0].has_target_rt, true)
  TEST_REAL_SIMILAR(transition_peak_maps[0].target_rt, 100.0)
}
END_SECTION

START_SECTION(void getPeakMaps_multi_file)
{
  const String file = writeTestFile_(7, "run1.mzML");
  std::vector<String> files = {file, file};
  XIPMParquetFile xipm(files);

  std::vector<XIPMParquetFile::XIPMPeakMap> peak_maps;
  xipm.getPeakMaps(peak_maps);
  TEST_EQUAL(peak_maps.size(), 4)
}
END_SECTION

START_SECTION(void getRuns(std::vector<XIPMRunInfo>& output) const)
{
  const String file = writeTestFile_(7, "run1.mzML");
  XIPMParquetFile xipm(file);

  std::vector<XIPMParquetFile::XIPMRunInfo> runs;
  xipm.getRuns(runs);
  TEST_EQUAL(runs.size(), 1)
  TEST_EQUAL(runs[0].run_id, 7)
}
END_SECTION

START_SECTION(void getColumns(std::vector<String>& output) const)
{
  const String file = writeTestFile_(7, "run1.mzML");
  XIPMParquetFile xipm(file);

  std::vector<String> columns;
  xipm.getColumns(columns);
  TEST_EQUAL(columns.empty(), false)

  bool has_mz = false;
  bool has_rt = false;
  bool has_im = false;
  bool has_target_rt = false;
  for (const auto& col : columns)
  {
    if (col == "MZ_DATA") has_mz = true;
    if (col == "RT_DATA") has_rt = true;
    if (col == "MOBILITY_DATA") has_im = true;
    if (col == "TARGET_RT") has_target_rt = true;
  }
  TEST_EQUAL(has_mz, true)
  TEST_EQUAL(has_rt, true)
  TEST_EQUAL(has_im, true)
  TEST_EQUAL(has_target_rt, true)
}
END_SECTION

START_SECTION(void load_invalid_path)
{
  TEST_EXCEPTION(Exception::FileNotFound, XIPMParquetFile("no_such_file.xipm"))
}
END_SECTION

END_TEST
