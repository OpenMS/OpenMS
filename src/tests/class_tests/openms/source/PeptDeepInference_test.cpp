// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ML/PeptDeepRTInference.h>
#include <OpenMS/ML/PeptDeepMS2Inference.h>
#include <OpenMS/ML/AminoAcidVocabulary.h>

using namespace OpenMS;
using namespace std;

START_TEST(PeptDeepInference, "$Id$")

// Note: These paths will be resolved by the OpenMS CMake testing environment
const string rt_model = "data/peptdeep_rt_dynamic.onnx";
const string ms2_model = "data/peptdeep_ms2_dynamic.onnx";

START_SECTION(AminoAcidVocabulary and Utilities)
    STATUS("Testing vocabulary tokenization for 'AGHCEWQMKYR'...");
    string pep = "AGHCEWQMKYR";

    // Expected index vector matching Justin's code review
    vector<Int64> expected = {0, 1, 7, 8, 3, 5, 23, 17, 13, 11, 25, 18, 0};

    vector<Int64> actual;
    actual.push_back(0); // N-term padding
    for (char aa : pep) {
        actual.push_back(ML::getAAIndex(aa));
    }
    actual.push_back(0); // C-term padding

    TEST_EQUAL(actual.size(), expected.size());
    for (size_t i = 0; i < expected.size(); ++i) {
        TEST_EQUAL(actual[i], expected[i]);
    }
END_SECTION

START_SECTION(PeptDeepRTInference)
    STATUS("Verifying RT prediction for iRT peptide LGGNEQVTR...");
    PeptDeepRTInference rt_engine(rt_model);

    vector<string> peptides = {"LGGNEQVTR"};
    auto rt_preds = rt_engine.predictRT(peptides);

    TEST_EQUAL(rt_preds.size(), 1);

    // Check against expected AlphaPeptDeep Python Output
    // Reference: https://github.com/MannLabs/alphapeptdeep/blob/main/nbs_tests/pretrained_models.ipynb
    TEST_REAL_SIMILAR(rt_preds[0], -0.035467f);
END_SECTION

START_SECTION(PeptDeepMS2Inference)
    STATUS("Verifying MS2 fragment intensity array lengths and math...");
    PeptDeepMS2Inference ms2_engine(ms2_model);

    vector<string> ms2_peptides = {"PEPTIDEK"};
    vector<float> charges = {2.0f};
    vector<float> nces = {30.0f};
    vector<int64_t> instruments = {0}; // 0 = Lumos

    auto ms2_preds = ms2_engine.predictMS2(ms2_peptides, charges, nces, instruments);

    TEST_EQUAL(ms2_preds.size(), 1);

    // PEPTIDEK length = 8. Valid fragments = 7.
    // AlphaPeptDeep outputs 8 ion types per fragment (b1, b2, y1, y2, etc.)
    // Exact expected length = (8 - 1) * 8 = 56
    size_t expected_length = (ms2_peptides[0].length() - 1) * 8;
    TEST_EQUAL(ms2_preds[0].size(), expected_length);

    // Verify Base Peak Normalization
    // Since we applied Base Peak Normalization, the max value in the spectrum MUST be exactly 1.0
    float max_intensity = 0.0f;
    for (float val : ms2_preds[0]) {
        if (val > max_intensity) {
            max_intensity = val;
        }
    }
    TEST_REAL_SIMILAR(max_intensity, 1.0f);

END_SECTION

END_TEST
