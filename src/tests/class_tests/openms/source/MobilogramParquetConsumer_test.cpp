// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Justin Sing $
// $Authors: Justin Sing $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/test_config.h>

#include <OpenMS/FORMAT/DATAACCESS/MobilogramParquetConsumer.h>
#include <OpenMS/FORMAT/XIMParquetFile.h>
#include <OpenMS/KERNEL/Mobilogram.h>
#include <OpenMS/KERNEL/MobilityPeak1D.h>
#include <OpenMS/SYSTEM/File.h>

using namespace OpenMS;
using namespace std;

START_TEST(MobilogramParquetConsumer, "$Id$")

START_SECTION(MobilogramParquetConsumer_basic)
{
  // create a small mobilogram
  Mobilogram m;
  m.setRT(123.45);
  m.push_back(MobilityPeak1D(1.0, 10.0));
  m.push_back(MobilityPeak1D(2.0, 20.0));

  OpenSwath::LightTargetedExperiment light_exp; // empty

  String tmp;
  NEW_TMP_FILE(tmp);
  String out = tmp + ".xim";
  {
    MobilogramParquetConsumer consumer(out, 1, "test_source", light_exp);
    consumer.consumeMobilogram(m);
  }

  TEST_EQUAL(File::exists(out), true)
}
END_SECTION

START_SECTION(MobilogramParquetConsumer_const_mobilogram_not_modified)
{
  // consumeMobilogram takes const Mobilogram& — verify the input is unchanged afterwards
  Mobilogram m;
  m.setRT(1.0);
  m.push_back(MobilityPeak1D(0.5, 100.0));
  m.push_back(MobilityPeak1D(1.5, 200.0));
  const Size initial_size = m.size();

  OpenSwath::LightTargetedExperiment light_exp;
  String tmp;
  NEW_TMP_FILE(tmp);
  String out = tmp + ".xim";
  {
    MobilogramParquetConsumer consumer(out, 1, "src", light_exp);
    consumer.consumeMobilogram(m);
  }
  // Original mobilogram must be untouched
  TEST_EQUAL(m.size(), initial_size)
  TEST_REAL_SIMILAR(m[0].getMobility(), 0.5)
  TEST_REAL_SIMILAR(m[1].getMobility(), 1.5)
}
END_SECTION

START_SECTION(MobilogramParquetConsumer_round_trip_file_written)
{
  // Write several mobilograms and verify the output file is non-empty.
  // Decoded value verification is covered by XIMParquetFile_test.cpp.
  Mobilogram m;
  m.setRT(50.0);
  m.push_back(MobilityPeak1D(0.80, 1000.0));
  m.push_back(MobilityPeak1D(0.90, 2000.0));
  m.push_back(MobilityPeak1D(1.00, 3000.0));

  OpenSwath::LightTargetedExperiment light_exp;
  String tmp;
  NEW_TMP_FILE(tmp);
  String out = tmp + ".xim";
  {
    MobilogramParquetConsumer consumer(out, 42, "run_src", light_exp);
    consumer.consumeMobilogram(m, "ms2", 2, -1, "tr1", 50.0, 7);
    consumer.consumeMobilogram(m, "ms2", 2, -1, "tr2", 50.0, 7);
  }

  TEST_EQUAL(File::exists(out), true)
  // File must contain actual parquet data (not just an empty shell)
  TEST_EQUAL(File::fileSize(out) > 0, true)
}
END_SECTION

START_SECTION(MobilogramParquetConsumer_empty_mobilograms)
{
  OpenSwath::LightTargetedExperiment light_exp;

  String tmp;
  NEW_TMP_FILE(tmp);
  String out = tmp + ".xim";
  {
    MobilogramParquetConsumer consumer(out, 1, "test_source", light_exp);
    // Do not consume any mobilograms
  }

  TEST_EQUAL(File::exists(out), true)
}
END_SECTION

START_SECTION(MobilogramParquetConsumer_destructor_flushes)
{
  // Destructor should flush pending data even without explicit finalize()
  Mobilogram m;
  m.push_back(MobilityPeak1D(1.0, 100.0));

  OpenSwath::LightTargetedExperiment light_exp;
  String tmp;
  NEW_TMP_FILE(tmp);
  String out = tmp + ".xim";
  {
    MobilogramParquetConsumer consumer(out, 1, "test_source", light_exp);
    consumer.consumeMobilogram(m);
    // No explicit finalize() — destructor must flush
  }

  // File was flushed by destructor
  TEST_EQUAL(File::exists(out), true)
  TEST_EQUAL(File::fileSize(out) > 0, true)
}
END_SECTION

START_SECTION(MobilogramParquetConsumer_mixed_target_and_decoy_identifying_transitions_keep_single_target_analyte)
{
  OpenSwath::LightTargetedExperiment light_exp;

  OpenSwath::LightCompound compound;
  compound.id = "2407";
  compound.sequence = "TSHVENDYIDNPSLALT(UniMod:21)TGPK";
  compound.charge = 3;
  light_exp.compounds.push_back(compound);

  OpenSwath::LightTransition target_transition;
  target_transition.transition_name = "tr_target";
  target_transition.peptide_ref = "2407";
  target_transition.precursor_mz = 784.6967;
  target_transition.product_mz = 326.1459;
  target_transition.fragment_charge = 1;
  target_transition.fragment_nr = 3;
  target_transition.setFragmentType("b");
  target_transition.setDetectingTransition(true);
  target_transition.setIdentifyingTransition(false);
  target_transition.setDecoy(false);
  light_exp.transitions.push_back(target_transition);

  OpenSwath::LightTransition identifying_decoy = target_transition;
  identifying_decoy.transition_name = "tr_ident_decoy";
  identifying_decoy.product_mz = 302.1710;
  identifying_decoy.setDetectingTransition(false);
  identifying_decoy.setIdentifyingTransition(true);
  identifying_decoy.setDecoy(true);
  light_exp.transitions.push_back(identifying_decoy);

  Mobilogram mobilogram;
  mobilogram.setRT(100.0);
  mobilogram.push_back(MobilityPeak1D(1.0, 1000.0));

  String tmp;
  NEW_TMP_FILE(tmp);
  String out = tmp + ".xim";
  {
    MobilogramParquetConsumer consumer(out, 1, "test_source", light_exp);
    consumer.consumeMobilogram(mobilogram, "transition", 2, -1, "tr_target", 100.0, 7);
    consumer.consumeMobilogram(mobilogram, "transition", 2, -1, "tr_ident_decoy", 100.0, 7);
  }

  XIMParquetFile xim(out);
  std::vector<XIMMobilogram> mobilograms;
  xim.load(mobilograms);
  TEST_EQUAL(mobilograms.size(), 2)
  for (const auto& m : mobilograms)
  {
    TEST_EQUAL(m.has_precursor_decoy, true)
    TEST_EQUAL(m.precursor_decoy, 0)
  }

  std::vector<XIMAnalyte> analytes;
  xim.getAnalytes(analytes);
  TEST_EQUAL(analytes.size(), 1)
}
END_SECTION

END_TEST
