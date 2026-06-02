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
#include <OpenMS/FORMAT/DATAACCESS/XIPMParquetConsumer.h>
#include <OpenMS/FORMAT/XIPMParquetFile.h>
#include <OpenMS/SYSTEM/File.h>

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

  OpenSwath::LightTargetedExperiment makeMixedIPFExperiment_()
  {
    OpenSwath::LightTargetedExperiment exp = makeExperiment_();

    OpenSwath::LightTransition identifying_decoy;
    identifying_decoy.transition_name = "tr1_ident_decoy";
    identifying_decoy.peptide_ref = "pep1";
    identifying_decoy.precursor_mz = 600.2;
    identifying_decoy.product_mz = 501.2;
    identifying_decoy.fragment_charge = 1;
    identifying_decoy.fragment_nr = 3;
    identifying_decoy.setFragmentType("b");
    identifying_decoy.setDetectingTransition(false);
    identifying_decoy.setIdentifyingTransition(true);
    identifying_decoy.setDecoy(true);
    exp.transitions.push_back(identifying_decoy);

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

  PeakMapExtractor::ExtractedPeakMap makeTransitionPeakMap_(const String& native_id,
                                                            const double target_mz,
                                                            const std::vector<double>& mz_values)
  {
    PeakMapExtractor::ExtractedPeakMap peak_map;
    peak_map.native_id = native_id;
    peak_map.target_mz = target_mz;
    peak_map.target_rt = 100.0;
    peak_map.target_ion_mobility = 1.1;
    peak_map.rt_start = 95.0;
    peak_map.rt_end = 105.0;
    peak_map.mz = mz_values;
    peak_map.rt = std::vector<double>(mz_values.size(), 100.0);
    peak_map.ion_mobility = std::vector<double>(mz_values.size(), 1.05);
    peak_map.intensity = std::vector<double>(mz_values.size(), 1000.0);
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
}

START_TEST(XIPMParquetConsumer, "$Id$")

START_SECTION(XIPMParquetConsumer_basic_roundtrip)
{
  const auto light_exp = makeExperiment_();

  String tmp;
  NEW_TMP_FILE(tmp);
  const String out = tmp + ".xipm";
  {
    XIPMParquetConsumer consumer(out, light_exp);
    consumer.consumePeakMap(makeTransitionPeakMap_(), 7, "run1.mzML", 2);
    consumer.consumePeakMap(makePrecursorPeakMap_(), 7, "run1.mzML", 1);
  }

  TEST_EQUAL(File::exists(out), true)
  TEST_EQUAL(File::fileSize(out) > 0, true)

  XIPMParquetFile xipm(out);
  std::vector<XIPMParquetFile::XIPMPeakMap> peak_maps;
  xipm.load(peak_maps);
  TEST_EQUAL(peak_maps.size(), 2)
  TEST_EQUAL(peak_maps[0].has_target_rt, true)
  TEST_REAL_SIMILAR(peak_maps[0].target_rt, 100.0)
  TEST_EQUAL(peak_maps[1].has_target_rt, true)
  TEST_REAL_SIMILAR(peak_maps[1].target_rt, 100.0)
}
END_SECTION

START_SECTION(XIPMParquetConsumer_empty_file)
{
  const auto light_exp = makeExperiment_();

  String tmp;
  NEW_TMP_FILE(tmp);
  const String out = tmp + ".xipm";
  {
    XIPMParquetConsumer consumer(out, light_exp);
  }

  TEST_EQUAL(File::exists(out), true)
}
END_SECTION

START_SECTION(XIPMParquetConsumer_mixed_target_and_decoy_identifying_transitions_keep_target_precursor_decoy)
{
  const auto light_exp = makeMixedIPFExperiment_();

  String tmp;
  NEW_TMP_FILE(tmp);
  const String out = tmp + ".xipm";
  {
    XIPMParquetConsumer consumer(out, light_exp);
    consumer.consumePeakMap(makeTransitionPeakMap_("tr1", 500.2, {500.19, 500.20}), 7, "run1.mzML", 2);
    consumer.consumePeakMap(makeTransitionPeakMap_("tr1_ident_decoy", 501.2, {501.19, 501.20}), 7, "run1.mzML", 2);
  }

  XIPMParquetFile xipm(out);
  std::vector<XIPMParquetFile::XIPMPeakMap> peak_maps;
  xipm.load(peak_maps);
  TEST_EQUAL(peak_maps.size(), 2)

  Size precursor_decoy_zero = 0;
  Size product_decoy_zero = 0;
  Size product_decoy_one = 0;
  for (const auto& peak_map : peak_maps)
  {
    TEST_EQUAL(peak_map.has_precursor_decoy, true)
    TEST_EQUAL(peak_map.precursor_decoy, 0)
    if (peak_map.precursor_decoy == 0) ++precursor_decoy_zero;
    TEST_EQUAL(peak_map.has_product_decoy, true)
    if (peak_map.product_decoy == 0) ++product_decoy_zero;
    if (peak_map.product_decoy == 1) ++product_decoy_one;
  }
  TEST_EQUAL(precursor_decoy_zero, 2)
  TEST_EQUAL(product_decoy_zero, 1)
  TEST_EQUAL(product_decoy_one, 1)
}
END_SECTION

END_TEST
