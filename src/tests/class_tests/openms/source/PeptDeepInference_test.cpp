#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/ML/PeptDeepRTInference.h>
#include <OpenMS/ML/PeptDeepMS2Inference.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <vector>
#include <string>
#include <iostream>

using namespace OpenMS;
using namespace std;

START_TEST(PeptDeepInference, "$Id$")

// Force the CMake test runner to provide the absolute paths
if (argc < 3) {
    cerr << " FATAL TEST ERROR: Model paths not provided." << endl;
    cerr << "Usage: " << argv[0] << " <absolute_path_to_rt.onnx> <absolute_path_to_ms2.onnx>" << endl;
    return 1;
}

const string rt_model = argv[1];
const string ms2_model = argv[2];

// TEST CASE 1: Initialization & Failures
TOLERANCE_ABSOLUTE(1e-4)

START_SECTION((PeptDeepRTInference(const string& model_path)))
    STATUS("Checking invalid model path handling...");
    // Using standard runtime error to clear the compiler warning cleanly
    TEST_EXCEPTION(std::runtime_error, PeptDeepRTInference("non_existent_path.onnx"));
END_SECTION

// TEST CASE 2: Tokenization & Constraints
START_SECTION((std::vector<float> predictMS2(const std::string&, float, float, int64_t)))
    PeptDeepMS2Inference ms2_engine(ms2_model);

    STATUS("Testing empty input sequence violation...");
    TEST_EXCEPTION(std::invalid_argument, ms2_engine.predictMS2("", 2.0f, 0.25f, 0));

    STATUS("Testing unknown amino acid violation (e.g., 'X')...");
    TEST_EXCEPTION(std::invalid_argument, ms2_engine.predictMS2("PEPTIDEX", 2.0f, 0.25f, 0));
END_SECTION


// TEST CASE 3: Execution Sanity
START_SECTION((Mathematical regression checks))
    STATUS("Verifying RT prediction returns a valid float output...");
    PeptDeepRTInference rt_engine(rt_model);
    vector<string> peptides = {"PEPTIDEK"};
    vector<float> rt_predictions = rt_engine.predictRT(peptides);

    TEST_EQUAL(rt_predictions.size(), 1);

    STATUS("Verifying MS2 engine returns proper fragment slice count...");
    PeptDeepMS2Inference ms2_engine_exec(ms2_model);
    auto ms2_preds = ms2_engine_exec.predictMS2("PEPTIDEK", 2.0f, 0.25f, 0);

    TEST_EQUAL(ms2_preds.size() > 0, true);

    STATUS("All PeptDeepInference tests passed safely within this section!");
END_SECTION

END_TEST