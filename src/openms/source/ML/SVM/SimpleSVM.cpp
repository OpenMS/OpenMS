// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/ML/SVM/SimpleSVM.h>
#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/SVOutStream.h>
#include <OpenMS/DATASTRUCTURES/ListUtils.h>

#include <OpenMS/MATH/StatisticFunctions.h>
#include <OpenMS/ML/GRIDSEARCH/GridSearch.h>

#include <cstdlib>

using namespace OpenMS;
using namespace std;


SimpleSVM::SimpleSVM():
  DefaultParamHandler("SimpleSVM"), data_(), model_(nullptr)
{
  defaults_.setValue("kernel", "RBF", "SVM kernel");
  defaults_.setValidStrings("kernel", {"RBF","linear"});

  defaults_.setValue("xval", 5, "Number of partitions for cross-validation (parameter optimization)");
  defaults_.setMinInt("xval", 1);

  String values = "-5,-3,-1,1,3,5,7,9,11,13,15";
  defaults_.setValue("log2_C", ListUtils::create<double>(values), "Values to try for the SVM parameter 'C' during parameter optimization. A value 'x' is used as 'C = 2^x'.");

  values = "-15,-13,-11,-9,-7,-5,-3,-1,1,3";
  defaults_.setValue("log2_gamma", ListUtils::create<double>(values), "Values to try for the SVM parameter 'gamma' during parameter optimization (RBF kernel only). A value 'x' is used as 'gamma = 2^x'.");

  values = "-15,-12,-9,-6,-3.32192809489,0,3.32192809489,6,9,12,15";
  defaults_.setValue("log2_p", ListUtils::create<double>(values), "Values to try for the SVM parameter 'epsilon' during parameter optimization (epsilon-SVR only). A value 'x' is used as 'epsilon = 2^x'.");

  vector<std::string> advanced(1, "advanced");
  defaults_.setValue("epsilon", 0.001, "Stopping criterion", advanced);
  defaults_.setMinFloat("epsilon", 0.0);

  defaults_.setValue("cache_size", 100.0, "Size of the kernel cache (in MB)", 
                     advanced);
  defaults_.setMinFloat("cache_size", 1.0);

  defaults_.setValue("no_shrinking", "false",
                     "Disable the shrinking heuristics", advanced);
  defaults_.setValidStrings("no_shrinking",
                            {"true","false"});

  defaultsToParam_();

  svm_set_print_string_function(&printNull_); // suppress output of LIBSVM
}


SimpleSVM::~SimpleSVM()
{
  clear_();
}

void SimpleSVM::clear_()
{
  if (model_ != nullptr) svm_free_and_destroy_model(&model_); // frees model *and* sets ptr to zero
  delete[] data_.x;
  delete[] data_.y;
}

void SimpleSVM::setup(PredictorMap& predictors, const map<Size, double>& outcomes, bool classification)
{
  if (predictors.empty() || predictors.begin()->second.empty())
  {
    throw Exception::IllegalArgument(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                     "Predictors for SVM must not be empty.");
  }

  // count elements for first feature dimension to determine number of observations
  Size n_obs = predictors.begin()->second.size();
  n_parts_ = param_.getValue("xval");

  // clear old models
  clear_();

  scaleData_(predictors);
  convertData_(predictors);

  data_.l = outcomes.size();
  data_.x = new svm_node*[data_.l];
  data_.y = new double[data_.l];
  map<double, Size> label_table;
  Size index = 0;
  for (auto it = outcomes.cbegin(); it != outcomes.cend();
       ++it, ++index)
  {
    const Size& training_index = it->first;
    const double& outcome = it->second;
    if (it->first >= n_obs)
    {
      String msg = "Invalid training index; there are only " + String(n_obs) +
        " observations.";
      throw Exception::InvalidValue(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                    msg, String(it->first));
    }
    data_.x[index] = &(nodes_[training_index][0]);
    data_.y[index] = outcome;
    label_table[outcome]++;
  }

  if (classification)
  {
    // check for 2 or more classes
    if (label_table.size() < 2)
    {
      throw Exception::MissingInformation(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                          "Need at least two classes (distinct "
                                          "labels) for SVM classification.");
    }

    String msg = "Training SVM on " + String(data_.l) + " observations. Classes:";
    for (map<double, Size>::iterator it = label_table.begin(); 
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
    OPENMS_LOG_INFO << msg << endl;

    svm_params_.svm_type = C_SVC;
  }
  else
  { // regression
    if ((unsigned int)data_.l < n_parts_) // TODO: check minimum amount of points needed for training a regression model. Assume 1 is enough for now.
    {
        String msg = "Not enough observations for " + String(n_parts_) + "-fold cross-validation.";
        throw Exception::MissingInformation(__FILE__, __LINE__, 
                                            OPENMS_PRETTY_FUNCTION, msg);
    }

    OPENMS_LOG_INFO << "Training SVR on " + String(data_.l) + " observations." << endl;

    svm_params_.svm_type = EPSILON_SVR;
    svm_params_.p = 0.1; // epsilon parameter of epsilon-SVR
  }

  std::string kernel = param_.getValue("kernel");
  svm_params_.kernel_type = (kernel == "RBF") ? RBF : LINEAR;
  svm_params_.eps = param_.getValue("epsilon");

  svm_params_.cache_size = param_.getValue("cache_size");
  svm_params_.shrinking = !param_.getValue("no_shrinking").toBool();
  svm_params_.nr_weight = 0; // weighting not supported for now
  svm_params_.probability = 0; // no prob. estimation during cross-validation

  optimizeParameters_(classification);
  svm_params_.probability = 1;

  model_ = svm_train(&data_, &svm_params_);
  OPENMS_LOG_INFO << "Number of support vectors in the final model: " << model_->l
           << endl;
}

// predict on (subset) of training data
void SimpleSVM::predict(vector<Prediction>& predictions, vector<Size> indexes) const
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
    for (Size i = 0; i < n_obs; indexes.push_back(i++)){};
  }
  Size n_classes = svm_get_nr_class(model_);
  vector<int> outcomes(n_classes);
  svm_get_labels(model_, &(outcomes[0]));
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
    pred.outcome = svm_predict_probability(model_, &(nodes_[*it][0]), 
                                             &(probabilities[0]));
    for (Size i = 0; i < n_classes; ++i)
    {
      pred.probabilities[outcomes[i]] = probabilities[i];
    }
    predictions.push_back(pred);
  }
}

void scaleDataUsingTrainingRanges(SimpleSVM::PredictorMap& predictors, const map<String, pair<double, double>>& scaling)
{
  // scale each feature dimension to the min-max-range
  for (auto pred_it = predictors.begin();
       pred_it != predictors.end(); ++pred_it)
  {
    if (pred_it->second.empty()) continue; // uninformative predictor
    auto val_begin = pred_it->second.begin();
    auto val_end = pred_it->second.end();
    for (; val_begin != val_end; ++val_begin)
    {
      if (scaling.count(pred_it->first) == 0)
      {
        //std::cout << "Predictor: '" << pred_it->first << "' not found in scale map because it was uninformative during training." << std::endl;
        continue;
      }      
      auto [min, max] = scaling.at(pred_it->first);
      double range = max - min;
       *val_begin = (*val_begin - min) / range;
    }
  }
  
}

// predict on novel e.g., test data
void SimpleSVM::predict(PredictorMap& predictors, vector<Prediction>& predictions) const
{
  if (model_ == nullptr)
  {
    throw Exception::Precondition(__FILE__, __LINE__, OPENMS_PRETTY_FUNCTION,
                                  "SVM/SVR model has not been trained (use the "
                                  "'setup' method)");
  }

  Size n_obs = predictors.begin()->second.size(); // length of the first feature ...
  Size feature_dim = predictors.size(); 

  scaleDataUsingTrainingRanges(predictors, scaling_);

  //std::cout << "Predicting on novel data with obs./feature dimensionality: " << n_obs << "/" << feature_dim << std::endl;

  Size n_classes = svm_get_nr_class(model_);
  vector<int> outcomes(n_classes);
  svm_get_labels(model_, &(outcomes[0]));
  vector<double> probabilities(n_classes);
  predictions.clear();
  predictions.reserve(n_obs);

  svm_node *x = new svm_node[feature_dim + 1];
  for (Size i = 0; i != n_obs; ++i)
  {
    size_t feature_index{0};
    for (auto p : predictors) 
    {
      x[feature_index].index = feature_index + 1;
      x[feature_index].value = p.second[i]; // feature value for observation i
      ++feature_index;
    }
    x[feature_dim].index = -1;
    x[feature_dim].value = 0;

    Prediction pred;
    pred.outcome = svm_predict_probability(model_, x, &(probabilities[0]));
    for (Size c = 0; c < n_classes; ++c)
    {
      pred.probabilities[outcomes[c]] = probabilities[c];
    }
    predictions.push_back(pred);
  }
  delete[] x;  
}

// only works in classification mode
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

void SimpleSVM::scaleData_(PredictorMap& predictors)
{
  scaling_.clear();
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
      OPENMS_LOG_INFO << "Predictor '" + pred_it->first + "' is uninformative. Ignoring." 
               << endl;
      pred_it->second.clear();
      continue;
    }
    double range = vmax - vmin;
    for (; val_begin != val_end; ++val_begin)
    {
      *val_begin = (*val_begin - vmin) / range;
    }
    scaling_[pred_it->first] = make_pair(vmin, vmax); // store old range
  }
}

const SimpleSVM::ScaleMap& SimpleSVM::getScaling() const
{
  return scaling_;
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
//      if (value > 0.0) // TODO: why > 0.0?
//      {
        svm_node node = {pred_index, value};
        nodes_[obs_index].push_back(node);
//      }
    }
  }
  OPENMS_LOG_DEBUG << "Number of predictors for SVM: " << pred_index << endl;
  svm_node sentinel = {-1, 0.0};
  for (auto node_it = nodes_.begin(); node_it != nodes_.end(); ++node_it)
  {
    node_it->push_back(sentinel);
  }
}

void SimpleSVM::writeXvalResults(const String& path) const
{
  SVOutStream output(path);
  output.modifyStrings(false);
  output << "log2_C" << "log2_gamma" << "log2_p" << "performance" << nl;
  for (Size g_index = 0; g_index < log2_gamma_.size(); ++g_index)
  {
    for (Size c_index = 0; c_index < log2_C_.size(); ++c_index)
    {
      for (Size p_index = 0; p_index < log2_p_.size(); ++p_index)
      {
        output << log2_C_[c_index] << log2_gamma_[g_index] << log2_p_[p_index]
               << performance_[g_index][c_index][p_index] << nl;
      }
    }
  }
}

tuple<double, double, double> SimpleSVM::chooseBestParameters_(bool higher_better) const
{
  auto is_better = [&higher_better](double l, double r)->bool { return higher_better ? (l > r) : (l < r); };

  // which parameter set(s) achieved best cross-validation performance?
  double best_value = higher_better ? std::numeric_limits<double>::lowest() : std::numeric_limits<double>::max();
  vector<tuple<Size, Size, Size>> best_indexes;
  for (Size g_index = 0; g_index < log2_gamma_.size(); ++g_index)
  {
    for (Size c_index = 0; c_index < log2_C_.size(); ++c_index)
    {
      for (Size p_index = 0; p_index < log2_p_.size(); ++p_index)
      {
        double value = performance_[g_index][c_index][p_index];
        // cout << "value " << value << " best_value " << best_value << endl;
        if (value == best_value)
        {
          best_indexes.emplace_back(g_index, c_index, p_index); // tie
        }
        else if (is_better(value, best_value))
        {
          best_value = value;
          best_indexes.clear();
          best_indexes.emplace_back(g_index, c_index, p_index);
        }
      }
    }
  }
  OPENMS_LOG_INFO << "Best cross-validation performance: " 
           << best_value << " (ties: " << best_indexes.size() << ")" << endl;

  if (best_indexes.size() == 1)
  {
    return make_tuple(
      log2_C_[std::get<1>(best_indexes[0])], // TODO: check why order changed
      log2_gamma_[std::get<0>(best_indexes[0])], 
      log2_p_[std::get<2>(best_indexes[0])]);
  }

  // break ties between parameter sets - look at "neighboring" parameters:
  multimap<pair<double, Size>, Size> tiebreaker;
  for (Size i = 0; i < best_indexes.size(); ++i)
  {
    const auto& indexes = best_indexes[i];
    Size n_neighbors = 0;
    double neighbor_value = 0.0;
    if (std::get<0>(indexes) > 0)
    {
      neighbor_value += performance_[std::get<0>(indexes) - 1][std::get<1>(indexes)][std::get<2>(indexes)];
      ++n_neighbors;
    }
    if (std::get<0>(indexes) + 1 < log2_gamma_.size())
    {
      neighbor_value += performance_[std::get<0>(indexes) + 1][std::get<1>(indexes)][std::get<2>(indexes)];
      ++n_neighbors;
    }
    if (std::get<1>(indexes) > 0)
    {
      neighbor_value += performance_[std::get<0>(indexes)][std::get<1>(indexes) - 1][std::get<2>(indexes)];
      ++n_neighbors;
    }
    if (std::get<1>(indexes) + 1 < log2_C_.size())
    {
      neighbor_value += performance_[std::get<0>(indexes)][std::get<1>(indexes) + 1][std::get<2>(indexes)];
      ++n_neighbors;
    }
    if (std::get<2>(indexes) > 0)
    {
      neighbor_value += performance_[std::get<0>(indexes)][std::get<1>(indexes)][std::get<2>(indexes) - 1];
      ++n_neighbors;
    }
    if (std::get<2>(indexes) + 1 < log2_p_.size())
    {
      neighbor_value += performance_[std::get<0>(indexes)][std::get<1>(indexes)][std::get<2>(indexes) + 1];
      ++n_neighbors;
    }
    neighbor_value /= n_neighbors; // avg. performance of neighbors
    tiebreaker.insert(make_pair(make_pair(neighbor_value, n_neighbors), i));
  }
  const auto& indexes = best_indexes[tiebreaker.rbegin()->second]; // TODO: use begin if higher_better == false

  return make_tuple(log2_C_[std::get<1>(indexes)],
      log2_gamma_[std::get<0>(indexes)], 
      log2_p_[std::get<2>(indexes)]);
}

void SimpleSVM::optimizeParameters_(bool classification)
{
  OPENMS_LOG_INFO << "Optimizing parameters." << endl;
  auto perFoldClassificationAccuracy = [&](const auto& d, const auto& targets)->double {
    Size n_correct = 0;
    for (Size i = 0; i < Size(d.l); ++i)
    {
      if (targets[i] == d.y[i]) n_correct++;
    }
    const double ratio = n_correct / double(d.l);
    return ratio;            
  };

  [[maybe_unused]]auto perFoldRegressionRSquared = [&](const auto& d, const auto& targets)->double {

    double targets_mean = Math::mean(std::begin(targets), std::end(targets)); // mean of truth y-values

    double u{}, v{u};
    for (Size i = 0; i < Size(d.l); ++i)
    {
      u += std::pow(targets[i] - d.y[i], 2.0);
      v += std::pow(targets[i] - targets_mean, 2.0);      
    }
    const double Rsquared = (v != 0.0) ? (1.0 - u/v) : -1.0;
    return Rsquared;
  };

  auto perFoldRMSE = [&](const auto& d, const auto& targets)->double {
    double err{};
    for (Size i = 0; i < Size(d.l); ++i)
    {
      err += std::pow(targets[i] - d.y[i], 2.0);
    }
    err /= (double)Size(d.l);
    err = std::sqrt(err);
    return err;
  };


  log2_C_ = param_.getValue("log2_C");
  if (svm_params_.kernel_type == RBF)
  {
    log2_gamma_ = param_.getValue("log2_gamma");
  }
  else
  {
    log2_gamma_ = vector<double>(1, 0.0);
  }
  log2_p_ = param_.getValue("log2_p");
  if (classification)
  {
    log2_p_ = vector<double >(1, log2(0.1));
  }

  OPENMS_LOG_INFO << "Running cross-validation to find optimal parameters..." 
          << endl;
  Size prog_counter = 0;
  ProgressLogger prog_log;
  prog_log.startProgress(1, log2_gamma_.size() * log2_C_.size() * log2_p_.size(), 
                        "testing parameters");

  const String& performance_type = classification ? "accuracy: " : "error: ";

  // classification performance for different parameter pairs:
  // vary "C"s in inner loop to keep results for all "C"s in one vector:

  performance_.resize(log2_gamma_.size());
  for (int g_index = 0; g_index < (int)log2_gamma_.size(); ++g_index)
  {
    svm_params_.gamma = pow(2.0, log2_gamma_[g_index]);
    performance_[g_index].resize(log2_C_.size());
    for (int c_index = 0; c_index < (int)log2_C_.size(); ++c_index)
    {
      svm_params_.C = pow(2.0, log2_C_[c_index]);
      performance_[g_index][c_index].resize(log2_p_.size());
      for (int p_index = 0; p_index < (int)log2_p_.size(); ++p_index)
      {
        svm_params_.p = pow(2.0, log2_p_[p_index]);

        double* targets = (double *)malloc(sizeof(double) * data_.l);
        svm_cross_validation(&data_, &svm_params_, n_parts_, &(targets[0]));

        double acc = classification ? perFoldClassificationAccuracy(data_, targets) : perFoldRMSE(data_, targets);

        performance_[g_index][c_index][p_index] = acc;
        prog_log.setProgress(++prog_counter);

        OPENMS_LOG_DEBUG << "Performance (log2_C = " << log2_C_[c_index] 
            << ", log2_gamma = " << log2_gamma_[g_index] << ") " 
            << ", log2_p = " << log2_p_[p_index] << ") "
            << performance_type << acc << endl;
        free(targets);
      }
    }
  }
  prog_log.endProgress();

  auto best_params = classification ? chooseBestParameters_(true) : chooseBestParameters_(false);
  //auto best_params = classification ? chooseBestParameters_(true) : chooseBestParameters_(true); // for Rsquared
  OPENMS_LOG_INFO << "Best SVM parameters: log2_C = " << std::get<0>(best_params)
          << ", log2_gamma = " << std::get<1>(best_params)
          << ", log2_p = " << std::get<2>(best_params)
          << endl;

  svm_params_.C = pow(2.0, std::get<0>(best_params));
  svm_params_.gamma = pow(2.0, std::get<1>(best_params));
  svm_params_.p = pow(2.0, std::get<2>(best_params));
  OPENMS_LOG_INFO << "... done." << endl;
}
