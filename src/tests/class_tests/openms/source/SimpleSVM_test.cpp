// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////
#include <OpenMS/ANALYSIS/SVM/SimpleSVM.h>
///////////////////////////

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
map<Size, Int> labels;

// read test data:
ifstream pred_file(OPENMS_GET_TEST_DATA_PATH("SimpleSVM_test_predictors.txt"));
if (pred_file.is_open())
{
  string line;
  while (getline(pred_file, line))
  {
    stringstream ss(line);
    String name;
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

START_SECTION((void setup(PredictorMap& predictors,
                          const map<Size, Int>& labels)))
{
  ABORT_IF(predictors.empty());
  ABORT_IF(labels.empty());

  SimpleSVM::PredictorMap empty_pred;
  TEST_EXCEPTION(Exception::IllegalArgument, svm.setup(empty_pred, labels));
  map<Size, Int> bad_labels;
  bad_labels[0] = 1;
  TEST_EXCEPTION(Exception::MissingInformation,
                 svm.setup(predictors, bad_labels));
  bad_labels[100] = 0;
  TEST_EXCEPTION(Exception::InvalidValue, svm.setup(predictors, bad_labels));

  svm.setup(predictors, labels);

  // check that data has been scaled:
  for (SimpleSVM::PredictorMap::const_iterator it = predictors.begin();
       it != predictors.end(); ++it)
  {
    vector<double>::const_iterator pos = min_element(it->second.begin(),
                                                     it->second.end());
    TEST_REAL_SIMILAR(*pos, 0);
    pos = max_element(it->second.begin(), it->second.end());
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
    if (it->label == 0)
    {
      TEST_EQUAL((it->probabilities[0] > 0.5) && (it->probabilities[0] < 1.0) &&
                 (it->probabilities[1] < 0.5) && (it->probabilities[1] > 0.0),
                 true);
    }
    else if (it->label == 1)
    {
      TEST_EQUAL((it->probabilities[1] > 0.5) && (it->probabilities[1] < 1.0) &&
                 (it->probabilities[0] < 0.5) && (it->probabilities[0] > 0.0),
                 true);
    }
    else TEST_EQUAL((it->label == 0) || (it->label == 1), true);
  }

  vector<Size> indexes(1, 0);
  svm.predict(predictions, indexes);
  TEST_EQUAL(predictions.size(), 1);

  indexes.push_back(100);
  TEST_EXCEPTION(Exception::InvalidValue, svm.predict(predictions, indexes));
}
END_SECTION

START_SECTION((void getFeatureWeights(map<String, double> feature_weights)
               const))
{
  map<String, double> feat_weights;
  TEST_EXCEPTION(Exception::Precondition,
                 untrained_svm.getFeatureWeights(feat_weights));

  svm.getFeatureWeights(feat_weights);
  TEST_EQUAL(feat_weights.size(), predictors.size());
}
END_SECTION

START_SECTION((void writeXvalResults(const String& path) const))
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

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
