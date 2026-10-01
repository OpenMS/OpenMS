// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>

#include <OpenMS/ANALYSIS/OPENSWATH/PeptDeepLibraryPredictor.h>
#include <OpenMS/CHEMISTRY/Residue.h>
#include <OpenMS/IONMOBILITY/IMTypes.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepCCSInference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepMS2Inference.h>
#include <OpenMS/ML/PEPTDEEP/PeptDeepRTInference.h>

#include <cstdint>
#include <string>
#include <vector>

using namespace OpenMS;
using namespace std;

START_TEST(PeptDeepLibraryPredictor, "$Id$")

const string rt_model = "data/peptdeep_rt_dynamic.onnx";
const string ccs_model = "data/peptdeep_ccs_dynamic.onnx";
const string ms2_model = "data/peptdeep_ms2_dynamic.onnx";

START_SECTION((OpenSwath::LightTargetedExperiment predict(const std::vector<PeptDeepLibraryPrecursor>& precursors)))
{
  const AASequence peptide = AASequence::fromString("PEPTIDEK");
  const AASequence modified_peptide = AASequence::fromString("M(Oxidation)PEPTIDE");
  const vector<string> peptide_strings{peptide.toString(), modified_peptide.toString()};
  const vector<float> charges{2.0f, 3.0f};
  const vector<float> nces{30.0f, 27.0f};
  const vector<int64_t> instruments{0, 2};

  float expected_rt = 0.0f;
  float expected_ccs = 0.0f;
  vector<float> expected_ms2;
  {
    PeptDeepRTInference rt_predictor(rt_model, 1, 2);
    expected_rt = rt_predictor.predictRT(peptide_strings).front();
  }
  {
    PeptDeepCCSInference ccs_predictor(ccs_model, 1, 2);
    expected_ccs = ccs_predictor.predictCCS(peptide_strings, charges).front();
  }
  {
    PeptDeepMS2Inference ms2_predictor(ms2_model, 1, 2);
    expected_ms2 = ms2_predictor.predictMS2(peptide_strings, charges, nces, instruments).front();
  }

  PeptDeepLibraryPredictor::Config config;
  config.rt_model_path = rt_model;
  config.ccs_model_path = ccs_model;
  config.ms2_model_path = ms2_model;
  config.intra_op_threads = 1;
  config.batch_size = 2;

  PeptDeepLibraryPredictor predictor(config);
  PeptDeepLibraryPrecursor precursor;
  precursor.peptide = peptide;
  precursor.id = "PEPTIDEK/2";
  precursor.charge = 2;
  precursor.nce = 30.0f;
  precursor.instrument_index = 0;
  precursor.protein_refs = {"P01234"};

  PeptDeepLibraryPrecursor modified_precursor;
  modified_precursor.peptide = modified_peptide;
  modified_precursor.id = "M(UniMod:35)PEPTIDE/3";
  modified_precursor.charge = 3;
  modified_precursor.nce = 27.0f;
  modified_precursor.instrument_index = 2;
  modified_precursor.protein_refs = {"P01234", "P56789"};

  const auto experiment = predictor.predict({precursor, modified_precursor});

  TEST_EQUAL(experiment.compounds.size(), 2)
  TEST_EQUAL(experiment.transitions.size(), 56)
  TEST_EQUAL(experiment.proteins.size(), 2)

  const auto& compound = experiment.compounds.front();
  const double precursor_mz = peptide.getMZ(2);
  const double expected_im = IMTypes::ccsToOneOverK0(expected_ccs, precursor_mz, 2);

  TEST_STRING_EQUAL(compound.id, "PEPTIDEK/2")
  TEST_STRING_EQUAL(compound.sequence, peptide.toUniModString())
  TEST_EQUAL(compound.charge, 2)
  TEST_REAL_SIMILAR(compound.rt, expected_rt)
  TEST_REAL_SIMILAR(compound.drift_time, expected_im)
  TEST_EQUAL(compound.protein_refs.size(), 1)
  TEST_STRING_EQUAL(compound.protein_refs.front(), "P01234")
  TEST_STRING_EQUAL(experiment.proteins[0].id, "P01234")
  TEST_STRING_EQUAL(experiment.proteins[1].id, "P56789")

  const auto& modified_compound = experiment.compounds[1];
  TEST_STRING_EQUAL(modified_compound.sequence, modified_peptide.toUniModString())
  TEST_EQUAL(modified_compound.modifications.size(), 1)
  TEST_EQUAL(modified_compound.modifications.front().location, 0)
  TEST_EQUAL(modified_compound.modifications.front().unimod_id, 35)

  const auto& b1_z1 = experiment.transitions[0];
  const auto& b1_z2 = experiment.transitions[1];
  const auto& y7_z1 = experiment.transitions[2];
  const auto& y7_z2 = experiment.transitions[3];

  TEST_STRING_EQUAL(b1_z1.transition_name, "PEPTIDEK/2_b1^1")
  TEST_STRING_EQUAL(b1_z1.peptide_ref, "PEPTIDEK/2")
  TEST_STRING_EQUAL(b1_z1.getFragmentType(), "b")
  TEST_EQUAL(b1_z1.fragment_nr, 1)
  TEST_EQUAL(b1_z1.fragment_charge, 1)
  TEST_REAL_SIMILAR(b1_z1.precursor_mz, precursor_mz)
  TEST_REAL_SIMILAR(b1_z1.precursor_im, expected_im)
  TEST_REAL_SIMILAR(b1_z1.product_mz, peptide.getPrefix(1).getMZ(1, Residue::BIon))
  TEST_REAL_SIMILAR(b1_z1.library_intensity, expected_ms2[0])

  TEST_STRING_EQUAL(b1_z2.getFragmentType(), "b")
  TEST_EQUAL(b1_z2.fragment_nr, 1)
  TEST_EQUAL(b1_z2.fragment_charge, 2)
  TEST_REAL_SIMILAR(b1_z2.library_intensity, expected_ms2[1])

  TEST_STRING_EQUAL(y7_z1.getFragmentType(), "y")
  TEST_EQUAL(y7_z1.fragment_nr, 7)
  TEST_EQUAL(y7_z1.fragment_charge, 1)
  TEST_REAL_SIMILAR(y7_z1.product_mz, peptide.getSuffix(7).getMZ(1, Residue::YIon))
  TEST_REAL_SIMILAR(y7_z1.library_intensity, expected_ms2[2])

  TEST_STRING_EQUAL(y7_z2.getFragmentType(), "y")
  TEST_EQUAL(y7_z2.fragment_nr, 7)
  TEST_EQUAL(y7_z2.fragment_charge, 2)
  TEST_REAL_SIMILAR(y7_z2.library_intensity, expected_ms2[3])
}
END_SECTION

START_SECTION((input validation and optional CCS prediction))
{
  PeptDeepLibraryPredictor::Config config;
  config.rt_model_path = rt_model;
  config.ccs_model_path = "/this/model/should/not/be/loaded.onnx";
  config.ms2_model_path = ms2_model;
  config.intra_op_threads = 1;
  config.batch_size = 2;
  config.predict_ccs = false;

  PeptDeepLibraryPredictor predictor(config);

  PeptDeepLibraryPrecursor precursor;
  precursor.peptide = AASequence::fromString("PEPTIDEK");
  precursor.id = "PEPTIDEK/2";
  precursor.charge = 2;

  const auto experiment = predictor.predict({precursor});
  TEST_REAL_SIMILAR(experiment.compounds.front().drift_time, -1.0)
  TEST_REAL_SIMILAR(experiment.transitions.front().precursor_im, -1.0)

  TEST_EXCEPTION(Exception::IllegalArgument, predictor.predict({}))

  PeptDeepLibraryPrecursor invalid = precursor;
  invalid.id = "";
  TEST_EXCEPTION(Exception::IllegalArgument, predictor.predict({invalid}))

  invalid = precursor;
  invalid.charge = 0;
  TEST_EXCEPTION(Exception::IllegalArgument, predictor.predict({invalid}))

  invalid = precursor;
  invalid.peptide = AASequence::fromString("A");
  TEST_EXCEPTION(Exception::IllegalArgument, predictor.predict({invalid}))

  TEST_EXCEPTION(Exception::IllegalArgument, predictor.predict({precursor, precursor}))
}
END_SECTION

END_TEST
