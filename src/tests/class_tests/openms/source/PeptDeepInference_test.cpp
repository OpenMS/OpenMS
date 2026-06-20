// Copyright (c) 2002-present, OpenMS Team -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Timo Sachsenberg $
// $Authors: Satyam Yadav $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ML/PeptDeepRTInference.h>
#include <OpenMS/ML/PeptDeepMS2Inference.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>
#include <string>
#include <iostream>
#include <OpenMS/test_config.h>

using namespace OpenMS;
using namespace std;

START_TEST(PeptDeepInference, "$Id$")

const string rt_model = "data/peptdeep_rt_dynamic.onnx";
const string ms2_model = "data/peptdeep_ms2_dynamic.onnx";

// TEST CASE 1: Initialization & Failures
TOLERANCE_ABSOLUTE(1e-5);

START_SECTION((PeptDeepRTInference(const string& model_path)))
    STATUS("Checking invalid model path handling...");
    TEST_EXCEPTION(Exception::FileNotFound, PeptDeepRTInference("non_existent_path.onnx"));
END_SECTION

// TEST CASE 2: Tokenization & Constraints
START_SECTION((std::vector<float> predictMS2(const std::string&, float, float, int64_t)))
    PeptDeepMS2Inference ms2_engine(ms2_model);

    STATUS("Testing empty input sequence violation...");
    TEST_EXCEPTION(Exception::IllegalArgument, ms2_engine.predictMS2("", 2.0f, 30.0f, 0));

    STATUS("Testing unsupported residue exception...");
    TEST_EXCEPTION(Exception::InvalidValue, ms2_engine.predictMS2("PEPTIDEX", 2.0f, 30.0f, 0));
END_SECTION

// TEST CASE 3: Execution Sanity & Numeric Regression
START_SECTION((Mathematical regression checks))
    STATUS("Verifying RT prediction returns a valid float output...");
    PeptDeepRTInference rt_engine(rt_model);
    vector<string> peptides = {"PEPTIDEK"};
    vector<float> rt_predictions = rt_engine.predictRT(peptides);

    TEST_EQUAL(rt_predictions.size(), 1);

    TEST_REAL_SIMILAR(rt_predictions[0], -0.0202387f);

    STATUS("Verifying MS2 engine returns normalized fragment slice...");
    PeptDeepMS2Inference ms2_engine_exec(ms2_model);

    // NCE expects raw percentage values (0-100)
    auto ms2_preds = ms2_engine_exec.predictMS2("PEPTIDEK", 2.0f, 30.0f, 0);

    TEST_EQUAL(ms2_preds.size() > 0, true);

    // Verify Base Peak Normalization matches Python output
    TEST_REAL_SIMILAR(ms2_preds[0], 0.0f);
    TEST_REAL_SIMILAR(ms2_preds[1], 0.0f);
    TEST_REAL_SIMILAR(ms2_preds[2], 0.0644772f);

    STATUS("All PeptDeepInference tests passed safely within this section!");
END_SECTION

END_TEST
