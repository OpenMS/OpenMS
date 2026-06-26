// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ML/SVM/SimpleSVM.h>
///////////////////////////

#include <fstream>
#include <sstream>


using namespace OpenMS;
using namespace std;

START_TEST(SimpleSVM, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

SimpleSVM* ptr = nullptr;
SimpleSVM* null_ptr = nullptr;

START_SECTION((SimpleSVM()))
{
  ptr = new SimpleSVM();
  TEST_NOT_EQUAL(ptr, null_ptr);
}
END_SECTION

START_SECTION((~SimpleSVM()))
{
  delete ptr;
}
END_SECTION

SimpleSVM svm, untrained_svm;
SimpleSVM::PredictorMap predictors;
map<Size, double> labels;

// read test data:
ifstream pred_file(OPENMS_GET_TEST_DATA_PATH("SimpleSVM_test_predictors.txt"));
if (pred_file.is_open())
{
  string line;
  while (getline(pred_file, line))
  {
    stringstream ss(line);
    std::string name;
    ss >> name;
    while (ss.good())
    {
      double value;
      ss >> value;
      predictors[name].push_back(value);
    }
  }
  pred_file.close();
}

ifstream label_file(OPENMS_GET_TEST_DATA_PATH("SimpleSVM_test_labels.txt"));
if (label_file.is_open())
{
  while (label_file.good())
  {
    Size index;
    Int label;
    label_file >> index >> label;
    labels[index] = label;
  }
  label_file.close();
}

START_SECTION((~SimpleSVM))
{
  SimpleSVM* svm_new = new SimpleSVM();
  delete svm_new;
}
END_SECTION

START_SECTION((void setup(PredictorMap& predictors,
                          const map<Size, double>& labels)))
{
  ABORT_IF(predictors.empty());
  ABORT_IF(labels.empty());

  SimpleSVM::PredictorMap empty_pred;
  TEST_EXCEPTION(Exception::IllegalArgument, svm.setup(empty_pred, labels));
  map<Size, double> bad_labels;
  bad_labels[0] = 1;
  SimpleSVM::PredictorMap tmp(predictors); // copy predictors to prevent rescaling to 0..1
  TEST_EXCEPTION(Exception::MissingInformation,
                 svm.setup(tmp, bad_labels));
  bad_labels[100] = 0;
  tmp = predictors; // copy predictors to prevent rescaling to 0..1
  TEST_EXCEPTION(Exception::InvalidValue, svm.setup(tmp, bad_labels));

  svm.setup(predictors, labels);

  // check that data has been scaled:
  for (const auto& it : predictors)
  {
    vector<double>::const_iterator pos = min_element(it.second.begin(),
                                                     it.second.end());
    TEST_REAL_SIMILAR(*pos, 0);
    pos = max_element(it.second.begin(), it.second.end());
    TEST_REAL_SIMILAR(*pos, 1);
  }
}
END_SECTION

START_SECTION((void predict(vector<Prediction>& predictions,
                            vector<Size> indexes) const))
{
  vector<SimpleSVM::Prediction> predictions;
  TEST_EXCEPTION(Exception::Precondition, untrained_svm.predict(predictions));

  svm.predict(predictions);

  TEST_EQUAL(predictions.size(), predictors.begin()->second.size());
  for (vector<SimpleSVM::Prediction>::iterator it = predictions.begin();
       it != predictions.end(); ++it)
  {
    if (it->outcome == 0)
    {
      TEST_EQUAL((it->probabilities[0] > 0.5) && (it->probabilities[0] < 1.0) &&
                 (it->probabilities[1] < 0.5) && (it->probabilities[1] > 0.0),
                 true);
    }
    else if (it->outcome == 1)
    {
      TEST_EQUAL((it->probabilities[1] > 0.5) && (it->probabilities[1] < 1.0) &&
                 (it->probabilities[0] < 0.5) && (it->probabilities[0] > 0.0),
                 true);
    }
    else TEST_EQUAL((it->outcome == 0) || (it->outcome == 1), true);
  }

  vector<Size> indexes(1, 0);
  svm.predict(predictions, indexes);
  TEST_EQUAL(predictions.size(), 1);

  indexes.push_back(100);
  TEST_EXCEPTION(Exception::InvalidValue, svm.predict(predictions, indexes));
}
END_SECTION

START_SECTION((void getFeatureWeights(map<std::string, double> feature_weights)
               const))
{
  map<std::string, double> feat_weights;
  TEST_EXCEPTION(Exception::Precondition,
                 untrained_svm.getFeatureWeights(feat_weights));

  svm.getFeatureWeights(feat_weights);
  TEST_EQUAL(feat_weights.size(), predictors.size());
}
END_SECTION

START_SECTION((map<std::string, pair<double, double>> void getScaling()
               const))
{
  auto scaling = svm.getScaling();
   
  TEST_REAL_SIMILAR(scaling["main_var_xx_swath_prelim_score"].first, -8.88447);
  TEST_REAL_SIMILAR(scaling["main_var_xx_swath_prelim_score"].second, 4.96923);
  TEST_REAL_SIMILAR(scaling["peak_apices_sum"].first, 1333.54);
  TEST_REAL_SIMILAR(scaling["peak_apices_sum"].second, 16131200);
  TEST_REAL_SIMILAR(scaling["rt_delta"].first, 0.0);
  TEST_REAL_SIMILAR(scaling["rt_delta"].second, 308.866);
  TEST_REAL_SIMILAR(scaling["sn_ratio"].first, 0.460926);
  TEST_REAL_SIMILAR(scaling["sn_ratio"].second, 176.142);
  TEST_REAL_SIMILAR(scaling["var_elution_model_fit_score"].first, -0.0801376);
  TEST_REAL_SIMILAR(scaling["var_elution_model_fit_score"].second, 0.998876);
  TEST_REAL_SIMILAR(scaling["var_intensity_score"].first,  0.00847651);
  TEST_REAL_SIMILAR(scaling["var_intensity_score"].second, 0.994559);
  TEST_REAL_SIMILAR(scaling["var_isotope_correlation_score"].first, -0.407676);
  TEST_REAL_SIMILAR(scaling["var_isotope_correlation_score"].second, 0.999652);
  TEST_REAL_SIMILAR(scaling["var_isotope_overlap_score"].first, 0.0);
  TEST_REAL_SIMILAR(scaling["var_isotope_overlap_score"].second, 1.0);
  TEST_REAL_SIMILAR(scaling["var_library_sangle"].first, 0.00104554)
  TEST_REAL_SIMILAR(scaling["var_library_sangle"].second, 1.17082);
  TEST_REAL_SIMILAR(scaling["var_log_sn_score"].first, 0.0);
  TEST_REAL_SIMILAR(scaling["var_log_sn_score"].second, 5.17129);
  TEST_REAL_SIMILAR(scaling["var_massdev_score"].first, 0.0)
  TEST_REAL_SIMILAR(scaling["var_massdev_score"].second, 34.4887);
  TEST_REAL_SIMILAR(scaling["var_xcorr_coelution"].first, 0.0)
  TEST_REAL_SIMILAR(scaling["var_xcorr_coelution"].second, 73.397);
  TEST_REAL_SIMILAR(scaling["var_xcorr_shape"].first, 0.333333)
  TEST_REAL_SIMILAR(scaling["var_xcorr_shape"].second, 0.999965);
  TEST_REAL_SIMILAR(scaling["xx_lda_prelim_score"].first, -4.61057)
  TEST_REAL_SIMILAR(scaling["xx_lda_prelim_score"].second, 7.95003);

  TEST_EQUAL(scaling.size(), predictors.size()); // one min/max entry for every predictor
}
END_SECTION


START_SECTION((void writeXvalResults(const std::string& path) const))
{
  string xval_file;
  NEW_TMP_FILE(xval_file);
  svm.writeXvalResults(xval_file);
  // cross-validation results are somewhat random, so don't be too strict:
  TOLERANCE_ABSOLUTE(0.2);
  TOLERANCE_RELATIVE(1.2);
  TEST_FILE_SIMILAR(xval_file,
                    OPENMS_GET_TEST_DATA_PATH("SimpleSVM_test_xval.txt"));
}
END_SECTION


START_SECTION(regression_train_and_predict_on_all)
{
  // create some noisy sinus data (data with clear outlier)
  // python code:
  // X = np.sort(5 * np.random.rand(40, 1), axis=0)
  // y = np.sin(X).ravel()
  // # add noise to targets
  // y[::5] += 3 * (0.5 - np.random.rand(8))

  SimpleSVM::PredictorMap x{ {"x", {
    0.0604460,     0.0827345,    0.0860808,     0.1587253, 
    0.2340380,     0.2376890,    0.2420991,     0.3825909, 
    0.4185480,     0.4909027,    0.7375445,     0.9507514, 
    0.9566240,     1.0680931,    1.1448901,     1.2531016, 
    1.3127954,     1.6548730,    1.8522109,     1.9021369, 
    2.0112177,     2.2142165,    2.2331448,     2.3437173, 
    2.7221007,     2.8684072,    2.9972565,     3.0146882, 
    3.0934168,     3.2048384,    3.2171502,     3.2959385, 
    3.6047516,     4.0592506,    4.0725104,     4.1267237, 
    4.7097008,     4.7699635,    4.8069954,     4.8126888 }}};
  
  map<Size, double> y{ {
    {0, -0.0604637},  {1, 0.0826401},   {2, 0.0859745},   {3, 0.1580597},   {4, 0.2319073},   {5, 1.1990320},
    {6, 0.2397410},   {7, 0.3733253},   {8, 0.4064343},   {9, 0.4714222},  {10, 2.1056848},  {11, 0.8138523},
   {12, 0.8172507},  {13, 0.8762834},  {14, 0.9106647},  {15, 1.4817945},  {16, 0.9669019},  {17, 0.9964676},
   {18, 0.9606635},  {19, 0.9456070},  {20, 1.3416947},  {21, 0.8000485},  {22, 0.7885501},  {23, 0.7158742},
   {24, 0.4072964},  {25, -0.1396564}, {26, 0.1438354},  {27, 0.1265640},  {28, 0.0481571},  {29, -0.0632036},
   {30, -0.2082548}, {31, -0.1537337}, {32, -0.4467765}, {33, -0.7941806}, {34, -0.8021682}, {35, -1.9997169},
   {36, -0.9999963}, {37, -0.9983430}, {38, -0.9955281}, {39, -0.9949741 } } };

  auto param = svm.getParameters();
  param.setValue("kernel", "RBF");
  param.setValue("log2_C", ListUtils::create<double>("9,11,13"));
  param.setValue("log2_gamma", ListUtils::create<double>("-1,1,3"));
  param.setValue("log2_p", ListUtils::create<double>("-6,-3.32192809489,0,3.32192809489"));

  svm.setParameters(param);

  svm.setup(x, y, false); // set up regression

  vector<SimpleSVM::Prediction> predictions;
  svm.predict(predictions);

  // test a few inlier
  TEST_EQUAL(std::abs(predictions[0].outcome - y[0]) < 0.2, true);
  TEST_EQUAL(std::abs(predictions[23].outcome - y[23]) < 0.2, true);
  TEST_EQUAL(std::abs(predictions[36].outcome - y[36]) < 0.2, true);

  // test a few outlier
  TEST_EQUAL(std::abs(predictions[10].outcome - y[10]) > 0.2, true);
  TEST_EQUAL(std::abs(predictions[15].outcome - y[15]) > 0.2, true);
  TEST_EQUAL(std::abs(predictions[35].outcome - y[35]) > 0.2, true);

  /* debug code to produce prediction error
  size_t index{};
  for (const auto& p : predictions)
  {
    std::cout << "index: " << index << " y: " << y[index] << " predicted y: " << p.outcome << " abs. error: " << std::abs(y[index] - p.outcome) << std::endl;
    ++index;
  }
  */

}
END_SECTION

START_SECTION(regression_train_and_predict_on_separate)
{
  // Same data as above but with split out test data (inlier + outlier)

  SimpleSVM::PredictorMap x{ {"x", {
    0.0827345,    0.0860808,     0.1587253, 
    0.2340380,     0.2376890,    0.2420991,     0.3825909, 
    0.4185480,     0.4909027,    0.9507514, 
    0.9566240,     1.0680931,    1.1448901,   
    1.3127954,     1.6548730,    1.8522109,     1.9021369, 
    2.0112177,     2.2142165,    2.2331448,   
    2.7221007,     2.8684072,    2.9972565,     3.0146882, 
    3.0934168,     3.2048384,    3.2171502,     3.2959385, 
    3.6047516,     4.0592506,    4.0725104,   
    4.7699635,     4.8069954,    4.8126888 }}};
  
  map<Size, double> y{ {
    {0, 0.0826401},   {1, 0.0859745},   {2, 0.1580597},   {3, 0.2319073},   {4, 1.1990320},
    {5, 0.2397410},   {6, 0.3733253},   {7, 0.4064343},   {8, 0.4714222},  {9, 0.8138523},
   {10, 0.8172507},  {11, 0.8762834},  {12, 0.9106647},  {13, 0.9669019},  {14, 0.9964676},
   {15, 0.9606635},  {16, 0.9456070},  {17, 1.3416947},  {18, 0.8000485},  {19, 0.7885501}, 
   {20, 0.4072964},  {21, -0.1396564}, {22, 0.1438354},  {23, 0.1265640},  {24, 0.0481571},  {25, -0.0632036},
   {26, -0.2082548}, {27, -0.1537337}, {28, -0.4467765}, {29, -0.7941806}, {30, -0.8021682},
   {31, -0.9983430}, {32, -0.9955281}, {33, -0.9949741 } } };

  SimpleSVM::PredictorMap x_test{ {"x", {
    0.0604460,
    0.7375445,
    1.2531016,
    2.3437173,
    4.1267237,
    4.7097008 }}};
  
  map<Size, double> y_test{ {
   {0, -0.0604637},
   {1, 2.1056848},
   {2, 1.4817945},
   {3, 0.7158742},   
   {4, -1.9997169},
   {5, -0.9999963} } };

  auto param = svm.getParameters();
  param.setValue("kernel", "RBF");
  param.setValue("log2_C", ListUtils::create<double>("1,5"));
  param.setValue("log2_gamma", ListUtils::create<double>("-5,5"));
  param.setValue("log2_p", ListUtils::create<double>("-15,-3.32192809489"));
  svm.setParameters(param);
  
  svm.setParameters(param);

  svm.setup(x, y, false); // set up regression

  vector<SimpleSVM::Prediction> predictions;
  svm.predict(x_test, predictions);

  // test a few inlier
  TEST_EQUAL(std::abs(predictions[0].outcome - y_test[0]) < 0.2, true);
  TEST_EQUAL(std::abs(predictions[3].outcome - y_test[3]) < 0.2, true);
  TEST_EQUAL(std::abs(predictions[5].outcome - y_test[5]) < 0.2, true);

  // test a few outlier
  TEST_EQUAL(std::abs(predictions[1].outcome - y_test[1]) > 0.2, true);
  TEST_EQUAL(std::abs(predictions[2].outcome - y_test[2]) > 0.2, true);
  TEST_EQUAL(std::abs(predictions[4].outcome - y_test[4]) > 0.2, true);
}
END_SECTION

START_SECTION([EXTRA] predict() stays index-aligned when a predictor is constant during training but present at prediction (issue #9661))
{
  // A predictor that is constant during training is dropped by setup() (it gets
  // no LIBSVM feature index). If predict() indexed features by prediction-map
  // position, such a predictor -- still present at prediction -- would shift
  // every later feature's index and silently corrupt the result. Here 'a'
  // (which sorts before 'b','c') is constant during training, so the trained
  // model must behave exactly as if it had only ever seen 'b','c'.

  // Separable two-class problem in informative features a,c. The constant
  // predictor 'b' sorts BETWEEN them, so it is dropped during training yet
  // shifts 'c' if predict() indexed by position. (It must NOT sort first: the
  // alphabetically-first predictor seeds convertData_'s observation count, and
  // a constant one there would be a different, setup-side issue.)
  auto make_train = [](bool with_b) {
    SimpleSVM::PredictorMap m;
    std::vector<double> a, c;
    for (int k = 0; k < 10; ++k) { a.push_back(0.02 * k);        c.push_back(0.02 * k); }        // class 0: 0.00..0.18
    for (int k = 0; k < 10; ++k) { a.push_back(0.82 + 0.02 * k); c.push_back(0.82 + 0.02 * k); } // class 1: 0.82..1.00
    m["a"] = a; m["c"] = c;
    if (with_b) { m["b"] = std::vector<double>(20, 7.0); } // constant -> dropped by setup()
    return m;
  };

  std::map<Size, double> y;
  for (Size i = 0; i < 10; ++i) { y[i] = 0.0; }
  for (Size i = 10; i < 20; ++i) { y[i] = 1.0; }

  Param param;
  param.setValue("kernel", "linear");
  param.setValue("log2_C", ListUtils::create<double>("0,3,6"));

  SimpleSVM svm_with_b, svm_without_b;
  svm_with_b.setParameters(param);
  svm_without_b.setParameters(param);
  SimpleSVM::PredictorMap train_with_b = make_train(true);
  SimpleSVM::PredictorMap train_without_b = make_train(false);
  svm_with_b.setup(train_with_b, y);       // 'b' constant -> dropped; trains on a,c
  svm_without_b.setup(train_without_b, y); // trains on a,c

  // The constant predictor really was dropped from the trained model.
  TEST_TRUE(svm_with_b.getScaling().find("b") == svm_with_b.getScaling().end())

  // Prediction data: identical a,c; for the 'b' model, 'b' is PRESENT and varies
  // (the exact trigger for the index misalignment).
  auto make_test = [](bool with_b) {
    SimpleSVM::PredictorMap t;
    t["a"] = std::vector<double>{0.05, 0.15, 0.85, 0.95};
    t["c"] = std::vector<double>{0.05, 0.15, 0.85, 0.95};
    if (with_b) { t["b"] = std::vector<double>{1.0, 9.0, 2.0, 8.0}; }
    return t;
  };
  SimpleSVM::PredictorMap test_with_b = make_test(true);
  SimpleSVM::PredictorMap test_without_b = make_test(false);

  std::vector<SimpleSVM::Prediction> pred_with_b, pred_without_b;
  svm_with_b.predict(test_with_b, pred_with_b);
  svm_without_b.predict(test_without_b, pred_without_b);

  TEST_EQUAL(pred_with_b.size(), 4)
  TEST_EQUAL(pred_without_b.size(), 4)

  // Invariant (#9661): a constant-in-training predictor present at prediction
  // must not change anything -> identical to the model that never saw 'b'.
  for (Size i = 0; i < pred_with_b.size(); ++i)
  {
    TEST_EQUAL(pred_with_b[i].outcome, pred_without_b[i].outcome)
    TEST_REAL_SIMILAR(pred_with_b[i].probabilities[0.0], pred_without_b[i].probabilities[0.0])
    TEST_REAL_SIMILAR(pred_with_b[i].probabilities[1.0], pred_without_b[i].probabilities[1.0])
  }

  // Sanity: the problem is separable, so the classes are recovered as expected
  // (this is exactly what silently breaks without the fix).
  TEST_EQUAL(pred_with_b[0].outcome, 0.0)
  TEST_EQUAL(pred_with_b[1].outcome, 0.0)
  TEST_EQUAL(pred_with_b[2].outcome, 1.0)
  TEST_EQUAL(pred_with_b[3].outcome, 1.0)
}
END_SECTION

START_SECTION([EXTRA] setup() handles a predictor that is constant during training and sorts first)
{
  // scaleData_ empties the value vector of every constant predictor, so
  // convertData_ must not take the observation count from a (now-empty) first
  // predictor. A constant predictor 'a' that sorts BEFORE the informative
  // features 'x','y' previously yielded n_obs == 0, sized the node arrays to
  // zero, and crashed setup(). It must now train and predict identically to a
  // model that never saw it.
  auto make_train = [](bool with_a) {
    SimpleSVM::PredictorMap m;
    std::vector<double> x, y;
    for (int k = 0; k < 10; ++k) { x.push_back(0.02 * k);        y.push_back(0.02 * k); }
    for (int k = 0; k < 10; ++k) { x.push_back(0.82 + 0.02 * k); y.push_back(0.82 + 0.02 * k); }
    m["x"] = x; m["y"] = y;
    if (with_a) { m["a"] = std::vector<double>(20, 3.0); } // constant, sorts first
    return m;
  };

  std::map<Size, double> lab;
  for (Size i = 0; i < 10; ++i) { lab[i] = 0.0; }
  for (Size i = 10; i < 20; ++i) { lab[i] = 1.0; }

  Param param;
  param.setValue("kernel", "linear");
  param.setValue("log2_C", ListUtils::create<double>("0,3,6"));

  SimpleSVM svm_a, svm_plain;
  svm_a.setParameters(param);
  svm_plain.setParameters(param);
  SimpleSVM::PredictorMap train_a = make_train(true);
  SimpleSVM::PredictorMap train_plain = make_train(false);
  svm_a.setup(train_a, lab);        // must NOT crash (previously n_obs == 0)
  svm_plain.setup(train_plain, lab);

  TEST_TRUE(svm_a.getScaling().find("a") == svm_a.getScaling().end())

  auto make_test = [](bool with_a) {
    SimpleSVM::PredictorMap t;
    t["x"] = std::vector<double>{0.05, 0.95};
    t["y"] = std::vector<double>{0.05, 0.95};
    if (with_a) { t["a"] = std::vector<double>{3.0, 3.0}; }
    return t;
  };
  SimpleSVM::PredictorMap test_a = make_test(true);
  SimpleSVM::PredictorMap test_plain = make_test(false);

  std::vector<SimpleSVM::Prediction> pred_a, pred_plain;
  svm_a.predict(test_a, pred_a);
  svm_plain.predict(test_plain, pred_plain);

  TEST_EQUAL(pred_a.size(), 2)
  TEST_EQUAL(pred_a[0].outcome, pred_plain[0].outcome)
  TEST_EQUAL(pred_a[1].outcome, pred_plain[1].outcome)
  TEST_EQUAL(pred_a[0].outcome, 0.0)
  TEST_EQUAL(pred_a[1].outcome, 1.0)

  // All predictors constant -> a clean error, not a crash or empty model.
  SimpleSVM svm_bad;
  svm_bad.setParameters(param);
  SimpleSVM::PredictorMap all_const;
  all_const["a"] = std::vector<double>(20, 1.0);
  all_const["b"] = std::vector<double>(20, 2.0);
  TEST_EXCEPTION(Exception::MissingInformation, svm_bad.setup(all_const, lab))
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
