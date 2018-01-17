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

#include <OpenMS/ANALYSIS/SVM/SimpleSVM.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/SVOutStream.h>

#include <fstream>

using namespace OpenMS;
using namespace std;


SimpleSVM::SimpleSVM():
  DefaultParamHandler("SimpleSVM"), data_(), model_(nullptr)
{
  defaults_.setValue("kernel", "RBF", "SVM kernel");
  defaults_.setValidStrings("kernel", ListUtils::create<String>("RBF,linear"));

  defaults_.setValue("xval", 5, "Number of partitions for cross-validation (parameter optimization)");
  defaults_.setMinInt("xval", 1);

  String values = "-5,-3,-1,1,3,5,7,9,11,13,15";
  defaults_.setValue("log2_C", ListUtils::create<double>(values), "Values to try for the SVM parameter 'C' during parameter optimization. A value 'x' is used as 'C = 2^x'.");

  values = "-15,-13,-11,-9,-7,-5,-3,-1,1,3";
  defaults_.setValue("log2_gamma", ListUtils::create<double>(values), "Values to try for the SVM parameter 'gamma' during parameter optimization (RBF kernel only). A value 'x' is used as 'gamma = 2^x'.");

  vector<String> advanced(1, "advanced");
  defaults_.setValue("epsilon", 0.001, "Stopping criterion", advanced);
  defaults_.setMinFloat("epsilon", 0.0);

  defaults_.setValue("cache_size", 100.0, "Size of the kernel cache (in MB)", 
                     advanced);
  defaults_.setMinFloat("cache_size", 1.0);

  defaults_.setValue("no_shrinking", "false",
                     "Disable the shrinking heuristics", advanced);
  defaults_.setValidStrings("no_shrinking",
                            ListUtils::create<String>("true,false"));

  defaultsToParam_();

  svm_set_print_string_function(&printNull_); // suppress output of LIBSVM
}


SimpleSVM::~SimpleSVM()
{
  if (model_ != nullptr) svm_free_model_content(model_);
  delete[] data_.x;
  delete[] data_.y;
}


void SimpleSVM::setup(PredictorMap& predictors, const map<Size, Int>& labels)
{
  if (predictors.empty() || predictors.begin()->second.empty())
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "Predictors for SVM must not be empty.");
  }
  Size n_obs = predictors.begin()->second.size();
  n_parts_ = param_.getValue("xval");

  scaleData_(predictors);
  convertData_(predictors);

  data_.l = labels.size();
  data_.x = new svm_node*[data_.l];
  data_.y = new double[data_.l];
  map<Int, Size> label_table;
  Size index = 0;
  for (map<Size, Int>::const_iterator it = labels.begin(); it != labels.end();
       ++it, ++index)
  {
    if (it->first >= n_obs)
    {
      String msg = "Invalid training index; there are only " + String(n_obs) +
        " observations.";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    msg, String(it->first));
    }
    data_.x[index] = &(nodes_[it->first][0]);
    data_.y[index] = double(it->second);
    label_table[it->second]++;
  }
  if (label_table.size() < 2)
  {
    throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                        "Need at least two classes (distinct "
                                        "labels) for SVM classification.");
  }
  String msg = "Training SVM on " + String(data_.l) + " observations. Classes:";
  for (map<Int, Size>::iterator it = label_table.begin(); 
       it != label_table.end(); ++it)
  {
    if (it->second < n_parts_)
    {
      msg = "Not enough observations of class " + String(it->first) + " for " +
        String(n_parts_) + "-fold cross-validation.";
      throw Exception::MissingInformation(__FILE__, __LINE__, 
                                          OPENMS_PRETTY_FUNCTION, msg);
    }
    msg += "\n- '" + String(it->first) + "': " + String(it->second) +
      " observations";
  }
  LOG_INFO << msg << endl;

  svm_params_.svm_type = C_SVC;
  String kernel = param_.getValue("kernel");
  svm_params_.kernel_type = (kernel == "RBF") ? RBF : LINEAR;
  svm_params_.eps = param_.getValue("epsilon");
  svm_params_.cache_size = param_.getValue("cache_size");
  svm_params_.shrinking = !param_.getValue("no_shrinking").toBool();
  svm_params_.nr_weight = 0; // weighting not supported for now
  svm_params_.probability = 0; // no prob. estimation during cross-validation

  optimizeParameters_();
  svm_params_.probability = 1;
  // in case "setup" was called before:
  if (model_ != nullptr) svm_free_model_content(model_);
  model_ = svm_train(&data_, &svm_params_);
  LOG_INFO << "Number of support vectors in the final model: " << model_->l
           << endl;
}


void SimpleSVM::predict(vector<Prediction>& predictions, vector<Size> indexes)
  const
{
  if (model_ == nullptr)
  {
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "SVM model has not been trained (use the "
                                  "'setup' method)");
  }

  Size n_obs = nodes_.size();
  if (indexes.empty())
  {
    indexes.reserve(n_obs);
    for (Size i = 0; i < n_obs; indexes.push_back(i++));
  }
  Size n_classes = svm_get_nr_class(model_);
  vector<Int> labels(n_classes);
  svm_get_labels(model_, &(labels[0]));
  vector<double> probabilities(n_classes);
  predictions.clear();
  predictions.reserve(indexes.size());
  for (vector<Size>::iterator it = indexes.begin(); it != indexes.end(); ++it)
  {
    if (*it >= n_obs)
    {
      String msg = "Invalid index for prediction; there are only " + 
        String(n_obs) + " observations.";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    msg, String(*it));
    }
    Prediction pred;
    pred.label = Int(svm_predict_probability(model_, &(nodes_[*it][0]), 
                                             &(probabilities[0])));
    for (Size i = 0; i < n_classes; ++i)
    {
      pred.probabilities[labels[i]] = probabilities[i];
    }
    predictions.push_back(pred);
  }
}


void SimpleSVM::getFeatureWeights(map<String, double>& feature_weights) const
{
  if (model_ == nullptr)
  {
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "SVM model has not been trained (use the "
                                  "'setup' method)");
  }
  Size k = model_->nr_class;
  if (k > 2)
  {
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "Output of feature weights is currently only "
                                  "supported for two-class classification");
  }

  feature_weights.clear();
  Size n_sv = model_->l; // number of support vectors
  for (Size l = 0; l < n_sv; ++l)
  {
    double sv_coef = model_->sv_coef[0][l];
    // LIBSVM uses a sparse representation for data (incl. support vectors):
    for (Size n = 0; ; ++n)
    {
      const struct svm_node& node = model_->SV[l][n];
      if (node.index == -1) break;
      const String& predictor_name = predictor_names_[node.index - 1];
      feature_weights[predictor_name] += sv_coef * node.value;
    }
  }
}


void SimpleSVM::scaleData_(PredictorMap& predictors) const
{
  for (PredictorMap::iterator pred_it = predictors.begin();
       pred_it != predictors.end(); ++pred_it)
  {
    // if (pred_it->second.empty()) continue;
    vector<double>::iterator val_begin = pred_it->second.begin();
    vector<double>::iterator val_end = pred_it->second.end();
    double vmin = *min_element(val_begin, val_end);
    double vmax = *max_element(val_begin, val_end);
    if (vmin == vmax)
    {
      LOG_INFO << "Predictor '" + pred_it->first + "' is uninformative." 
               << endl;
      pred_it->second.clear();
      continue;
    }
    double range = vmax - vmin;
    for (; val_begin != val_end; ++val_begin)
    {
      *val_begin = (*val_begin - vmin) / range;
    }
  }
}


void SimpleSVM::convertData_(const PredictorMap& predictors)
{
  Size n_obs = predictors.begin()->second.size();
  nodes_.clear();
  nodes_.resize(n_obs);
  predictor_names_.clear();
  int pred_index = 0; // "int" for use by LIBSVM
  for (PredictorMap::const_iterator pred_it = predictors.begin();
       pred_it != predictors.end(); ++pred_it)
  {
    if (pred_it->second.empty()) continue; // uninformative predictor
    pred_index++; // LIBSVM counts observations from 1
    predictor_names_.push_back(pred_it->first);
    for (Size obs_index = 0; obs_index < n_obs; ++obs_index)
    {
      double value = pred_it->second[obs_index];
      if (value > 0.0)
      {
        svm_node node = {pred_index, value};
        nodes_[obs_index].push_back(node);
      }
    }
  }
  LOG_DEBUG << "Number of predictors for SVM: " << pred_index << endl;
  svm_node final = {-1, 0.0};
  for (vector<vector<struct svm_node> >::iterator node_it = nodes_.begin();
       node_it != nodes_.end(); ++node_it)
  {
    node_it->push_back(final);
  }
}


void SimpleSVM::writeXvalResults(const String& path) const
{
  SVOutStream output(path);
  output.modifyStrings(false);
  output << "log2_C" << "log2_gamma" << "performance" << nl;
  for (Size g_index = 0; g_index < log2_gamma_.size(); ++g_index)
  {
    for (Size c_index = 0; c_index < log2_C_.size(); ++c_index)
    {
      output << log2_C_[c_index] << log2_gamma_[g_index] 
             << performance_[g_index][c_index] << nl;
    }
  }
}


pair<double, double> SimpleSVM::chooseBestParameters_() const
{
  // which parameter set(s) achieved best cross-validation performance?
  double best_value = 0.0;
  vector<pair<Size, Size> > best_indexes;
  for (Size g_index = 0; g_index < log2_gamma_.size(); ++g_index)
  {
    for (Size c_index = 0; c_index < log2_C_.size(); ++c_index)
    {
      double value = performance_[g_index][c_index];
      if (value == best_value)
      {
        best_indexes.push_back(make_pair(g_index, c_index));
      }
      else if (value > best_value)
      {
        best_value = value;
        best_indexes.clear();
        best_indexes.push_back(make_pair(g_index, c_index));
      }
    }
  }
  LOG_INFO << "Best cross-validation performance: " 
           << float(best_value * 100.0) << "% correct" << endl;
  if (best_indexes.size() == 1)
  {
    return make_pair(log2_C_[best_indexes[0].second],
                     log2_gamma_[best_indexes[0].first]);
  }
  // break ties between parameter sets - look at "neighboring" parameters:
  multimap<pair<double, Size>, Size> tiebreaker;
  for (Size i = 0; i < best_indexes.size(); ++i)
  {
    const pair<Size, Size>& indexes = best_indexes[i];
    Size n_neighbors = 0;
    double neighbor_value = 0.0;
    if (indexes.first > 0)
    {
      neighbor_value += performance_[indexes.first - 1][indexes.second];
      ++n_neighbors;
    }
    if (indexes.first + 1 < log2_gamma_.size())
    {
      neighbor_value += performance_[indexes.first + 1][indexes.second];
      ++n_neighbors;
    }
    if (indexes.second > 0)
    {
      neighbor_value += performance_[indexes.first][indexes.second - 1];
      ++n_neighbors;
    }
    if (indexes.second + 1 < log2_C_.size())
    {
      neighbor_value += performance_[indexes.first][indexes.second + 1];
      ++n_neighbors;
    }
    neighbor_value /= n_neighbors; // avg. performance of neighbors
    tiebreaker.insert(make_pair(make_pair(neighbor_value, n_neighbors), i));
  }
  const pair<Size, Size>& indexes = best_indexes[tiebreaker.rbegin()->second];
  return make_pair(log2_C_[indexes.second], log2_gamma_[indexes.first]);
}


void SimpleSVM::optimizeParameters_()
{
  log2_C_ = param_.getValue("log2_C");
  if (svm_params_.kernel_type == RBF)
  {
    log2_gamma_ = param_.getValue("log2_gamma");
  }
  else
  {
    log2_gamma_ = vector<double>(1, 0.0);
  }

  LOG_INFO << "Running cross-validation to find optimal SVM parameters..." 
           << endl;
  Size prog_counter = 0;
  ProgressLogger prog_log;
  prog_log.startProgress(1, log2_gamma_.size() * log2_C_.size(),
                         "testing SVM parameters");
  // classification performance for different parameter pairs:
  performance_.resize(log2_gamma_.size());
  // vary "C"s in inner loop to keep results for all "C"s in one vector:
  for (Size g_index = 0; g_index < log2_gamma_.size(); ++g_index)
  {
    svm_params_.gamma = pow(2.0, log2_gamma_[g_index]);
    performance_[g_index].resize(log2_C_.size());
    for (Size c_index = 0; c_index < log2_C_.size(); ++c_index)
    {
      svm_params_.C = pow(2.0, log2_C_[c_index]);
      vector<double> targets(data_.l);
      svm_cross_validation(&data_, &svm_params_, n_parts_, &(targets[0]));
      Size n_correct = 0;
      for (Size i = 0; i < Size(data_.l); ++i)
      {
        if (targets[i] == data_.y[i]) n_correct++;
      }
      double ratio = n_correct / double(data_.l);
      performance_[g_index][c_index] = ratio;
      prog_log.setProgress(++prog_counter);
      LOG_DEBUG << "Performance (log2_C = " << log2_C_[c_index] 
                << ", log2_gamma = " << log2_gamma_[g_index] << "): " 
                << n_correct << " correct (" << float(ratio * 100.0) << "%)"
                << endl;
    }
  }
  prog_log.endProgress();

  pair<double, double> best_params = chooseBestParameters_();
  LOG_INFO << "Best SVM parameters: log2_C = " << best_params.first
           << ", log2_gamma = " << best_params.second << endl;

  svm_params_.C = pow(2.0, best_params.first);
  svm_params_.gamma = pow(2.0, best_params.second);
}

