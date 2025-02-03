// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Douglas McCloskey, Pasquale Domenico Colaianni, Svetlana Kutuzova $
// $Authors: Douglas McCloskey, Pasquale Domenico Colaianni, Svetlana Kutuzova $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/OPENSWATH/MRMFeatureSelector.h>

#include <OpenMS/CONCEPT/LogStream.h>
#include <OpenMS/DATASTRUCTURES/LPWrapper.h>
#include <OpenMS/KERNEL/FeatureMap.h>
#include <algorithm>
#include <set>

namespace OpenMS
{
  Int MRMFeatureSelector::addVariable_(
    LPWrapper& problem,
    const String& name,
    const bool bounded,
    const double obj,
    const VariableType variableType
  ) const
  {
    const Int index = problem.addColumn();

    if (bounded)
    {
      problem.setColumnBounds(index, 0, 1, LPWrapper::DOUBLE_BOUNDED);
    }
    else
    {
      problem.setColumnBounds(index, 0, 1, LPWrapper::UNBOUNDED);
    }

    problem.setColumnName(index, name);

    if (variableType == VariableType::INTEGER)
    {
      problem.setColumnType(index, LPWrapper::INTEGER);
    }
    else if (variableType == VariableType::CONTINUOUS)
    {
      problem.setColumnType(index, LPWrapper::CONTINUOUS);
    }
    else
    {
      throw std::runtime_error("Variable type not supported\n");
    }

    problem.setObjective(index, obj);
    return index;
  }

  void MRMFeatureSelector::addConstraint_(
    LPWrapper& problem,
    const std::vector<Int>& indices,
    const std::vector<double>& values,
    const String& name,
    const double lb,
    const double ub,
    const LPWrapper::Type param
  ) const
  {
    problem.addRow(indices, values, name, lb, ub, param);
  }

  void MRMFeatureSelectorScore::optimize(
    const std::vector<std::pair<double, String>>& time_to_name,
    const std::map<String, std::vector<Feature>>& feature_name_map,
    std::vector<String>& result,
    const SelectorParameters& parameters
  ) const
  {
    result.clear();
    std::set<String> variables;
    LPWrapper problem;
    problem.setObjectiveSense(LPWrapper::MIN);
    for (const std::pair<double, String>& elem : time_to_name)
    {
      std::vector<Int> constraints;
      for (const Feature& feature : feature_name_map.at(elem.second))
      {
        const String name1 = elem.second + "_" + String(feature.getUniqueId());
        if (variables.count(name1) == 0)
        {
          const double score = computeScore_(feature, parameters.score_weights);
          const Int col_idx = addVariable_(problem, name1, true, score, parameters.variable_type);
          constraints.push_back(col_idx);
          variables.insert(name1);
        }
      }
      std::vector<double> constraints_values(constraints.size(), 1.0);
      addConstraint_(problem, constraints, constraints_values, elem.second + "_constraint", 1.0, 1.0, LPWrapper::DOUBLE_BOUNDED);
    }
    LPWrapper::SolverParam param;
    problem.solve(param);
    for (Int c = 0; c < problem.getNumberOfColumns(); ++c)
    {
      if (problem.getColumnValue(c) >= parameters.optimal_threshold)
      {
        result.push_back(problem.getColumnName(c));
      }
    }
  }

  String MRMFeatureSelector::removeSpaces_(String str) const
  {
    String::iterator end_pos = std::remove(str.begin(), str.end(), ' ');
    str.erase(end_pos, str.end());
    return str;
  }

  void MRMFeatureSelectorQMIP::optimize(
    const std::vector<std::pair<double, String>>& time_to_name,
    const std::map<String, std::vector<Feature>>& feature_name_map,
    std::vector<String>& result,
    const SelectorParameters& parameters
  ) const
  {
    result.clear();
    std::set<String> variables;
    LPWrapper problem;
    problem.setObjectiveSense(LPWrapper::MIN);
    for (Int cnt1 = 0; static_cast<Size>(cnt1) < time_to_name.size(); ++cnt1)
    {
      const Size start_iter = std::max(cnt1 - parameters.nn_threshold, 0);
      const Size stop_iter = std::min(static_cast<Size>(cnt1 + parameters.nn_threshold + 1), time_to_name.size()); // assuming nn_threshold >= -1
      std::vector<Int> constraints;
      const std::vector<Feature>& feature_row1 = feature_name_map.at(time_to_name[cnt1].second);

      for (Size i = 0; i < feature_row1.size(); ++i)
      {
        const String name1 = time_to_name[cnt1].second + "_" + String(feature_row1[i].getUniqueId());

        if (variables.count(name1) == 0)
        {
          constraints.push_back(addVariable_(problem, name1, true, 0, parameters.variable_type));
          variables.insert(name1);
        }
        else
        {
          constraints.push_back(problem.getColumnIndex(name1));
        }

        double score_1 = computeScore_(feature_row1[i], parameters.score_weights);
        const Size n_score_weights = parameters.score_weights.size();

        if (n_score_weights > 1)
        {
          score_1 = std::pow(score_1, 1.0 / n_score_weights);
        }

        const Int index1 = problem.getColumnIndex(name1);

        for (Size cnt2 = start_iter; cnt2 < stop_iter; ++cnt2)
        {
          if (static_cast<Size>(cnt1) == cnt2)
          {
            continue;
          }

          const std::vector<Feature>& feature_row2 = feature_name_map.at(time_to_name[cnt2].second);
          const double locality_weight = parameters.locality_weight
            ? 1.0 / (parameters.nn_threshold - std::abs(static_cast<Int>(start_iter + cnt2) - cnt1) + 1)
            : 1.0;
          const double tr_delta_expected = time_to_name[cnt1].first - time_to_name[cnt2].first;

          for (Size j = 0; j < feature_row2.size(); ++j)
          {
            const String name2 = time_to_name[cnt2].second + "_" + String(feature_row2[j].getUniqueId());
            if (variables.count(name2) == 0)
            {
              addVariable_(problem, name2, true, 0, parameters.variable_type);
              variables.insert(name2);
            }

            const String var_qp_name = time_to_name[cnt1].second + "_" + String(i) + "-" + time_to_name[cnt2].second + "_" + String(j);

            const Int index_var_qp = addVariable_(problem, var_qp_name, true, 0, VariableType::CONTINUOUS);
            const Int index_var_abs = addVariable_(problem, var_qp_name + "-ABS", false, 1, VariableType::CONTINUOUS);

            const Int index2 = problem.getColumnIndex(name2);

            double score_2 = computeScore_(feature_row2[j], parameters.score_weights);
            if (n_score_weights > 1)
            {
              score_2 = std::pow(score_2, 1.0 / n_score_weights);
            }

            const double tr_delta = feature_row1[i].getRT() - feature_row2[j].getRT();
            const double score = locality_weight * score_1 * score_2 * (tr_delta - tr_delta_expected);

            addConstraint_(problem, {index1, index_var_qp}, {1.0, -1.0}, var_qp_name + "-QP1", 0.0, 1.0, LPWrapper::LOWER_BOUND_ONLY);
            addConstraint_(problem, {index2, index_var_qp}, {1.0, -1.0}, var_qp_name + "-QP2", 0.0, 1.0, LPWrapper::LOWER_BOUND_ONLY);
            addConstraint_(problem, {index1, index2, index_var_qp}, {1.0, 1.0, -1.0}, var_qp_name + "-QP3", 0.0, 1.0, LPWrapper::UPPER_BOUND_ONLY);
            std::vector<Int> indices_abs = {index_var_abs, index_var_qp};
            addConstraint_(problem, indices_abs, {-1.0, score}, var_qp_name + "-obj+", -1.0, 0.0, LPWrapper::UPPER_BOUND_ONLY);
            addConstraint_(problem, indices_abs, {-1.0, -score}, var_qp_name + "-obj-", -1.0, 0.0, LPWrapper::UPPER_BOUND_ONLY);
          }
        }
      }
      std::vector<double> constraints_values(constraints.size(), 1.0);
      addConstraint_(problem, constraints, constraints_values, time_to_name[cnt1].second + "_constraint", 1.0, 1.0, LPWrapper::DOUBLE_BOUNDED);
      // addConstraint_(problem, constraints, constraints_values, time_to_name[cnt1].second + "_constraint", 1.0, 1.0, LPWrapper::FIXED); // glpk
    }
    LPWrapper::SolverParam param;
    problem.solve(param);
    for (Int c = 0; c < problem.getNumberOfColumns(); ++c)
    {
      const String name = problem.getColumnName(c);
      if (problem.getColumnValue(c) > parameters.optimal_threshold && variables.count(name))
      {
        result.push_back(name);
      }
    }
  }

  void MRMFeatureSelector::constructTargTransList_(
    const FeatureMap& features,
    std::vector<std::pair<double, String>>& time_to_name,
    std::map<String, std::vector<Feature>>& feature_name_map,
    const bool select_transition_group
  ) const
  {
    time_to_name.clear();
    feature_name_map.clear();
    std::set<String> names;
    for (const Feature& feature : features)
    {
      const String component_group_name = removeSpaces_(feature.getMetaValue("PeptideRef").toString());
      const double assay_retention_time = feature.getMetaValue("assay_rt");
      if (names.count(component_group_name) == 0)
      {
        time_to_name.emplace_back(assay_retention_time, component_group_name);
        names.insert(component_group_name);
      }
      if (feature_name_map.count(component_group_name) == 0)
      {
        feature_name_map[component_group_name] = std::vector<Feature>();
      }
      feature_name_map[component_group_name].push_back(feature);
      if (select_transition_group)
      {
        continue;
      }
      for (const Feature& subordinate : feature.getSubordinates())
      {
        const String component_name = removeSpaces_(subordinate.getMetaValue("native_id").toString());
        if (names.count(component_name))
        {
          time_to_name.emplace_back(assay_retention_time, component_name);
          names.insert(component_name);
        }
        if (feature_name_map.count(component_name) == 0)
        {
          feature_name_map[component_name] = std::vector<Feature>();
        }
        feature_name_map[component_name].push_back(subordinate);
      }
    }
  }

  void MRMFeatureSelector::selectMRMFeature(
    const FeatureMap& features,
    FeatureMap& selected_filtered,
    const SelectorParameters& parameters
  ) const
  {
    selected_filtered.clear();

    if (features.empty())
    {
      return;
    }

    std::vector<std::pair<double, String>> time_to_name;
    std::map<String, std::vector<Feature>> feature_name_map;
    constructTargTransList_(features, time_to_name, feature_name_map, parameters.select_transition_group);

    sort(time_to_name.begin(), time_to_name.end());
    Int window_length = parameters.segment_window_length;
    Int step_length = parameters.segment_step_length;
    if (window_length == -1 && step_length == -1)
    {
      window_length = step_length = time_to_name.size();
    }
    Size n_segments = time_to_name.size() / step_length;
    if (time_to_name.size() % step_length)
    {
      ++n_segments;
    }
    std::vector<String> result_names;

    for (Size i = 0; i < n_segments; ++i)
    {
      const Size start = step_length * i;
      const Size end = std::min(start + window_length, time_to_name.size());
      const std::vector<std::pair<double, String>> time_slice(time_to_name.begin() + start, time_to_name.begin() + end);
      std::vector<String> result;
      optimize(time_slice, feature_name_map, result, parameters);
      result_names.insert(result_names.end(), result.begin(), result.end());
    }
    const std::set<String> result_names_set(result_names.begin(), result_names.end());
    for (const Feature& feature : features)
    {
      std::vector<Feature> subordinates_filtered;
      for (const Feature& subordinate : feature.getSubordinates())
      {
        const String feature_name = parameters.select_transition_group
          ? removeSpaces_(feature.getMetaValue("PeptideRef").toString()) + "_" + String(feature.getUniqueId())
          : removeSpaces_(subordinate.getMetaValue("native_id").toString()) + "_" + String(feature.getUniqueId());

        if (result_names_set.count(feature_name))
        {
          subordinates_filtered.push_back(subordinate);
        }
      }
      if (!subordinates_filtered.empty())
      {
        Feature feature_filtered(feature);
        feature_filtered.setSubordinates(subordinates_filtered);
        selected_filtered.push_back(feature_filtered);
      }
    }
  }

  double MRMFeatureSelector::computeScore_(const Feature& feature, const std::map<String, MRMFeatureSelector::LambdaScore>& score_weights) const
  {
    double score_1 = 1.0;
    for (const std::pair<const String, LambdaScore>& score_weight : score_weights)
    {
      const String& metavalue_name = score_weight.first;
      const LambdaScore lambda_score = score_weight.second;
      if (!feature.metaValueExists(metavalue_name))
      {
        OPENMS_LOG_WARN << "computeScore_(): Metavalue \"" << metavalue_name << "\" not found.\n";
        continue;
      }
      const double value = weightScore_(feature.getMetaValue(metavalue_name), lambda_score);
      if (value > 0.0 && !std::isnan(value) && !std::isinf(value))
      {
        score_1 *= value;
      }
    }
    return score_1;
  }

  double MRMFeatureSelector::weightScore_(const double score, const LambdaScore lambda_score) const
  {
    if (lambda_score == LambdaScore::LINEAR)
    {
      return score;
    }
    else if (lambda_score == LambdaScore::INVERSE)
    {
      return 1.0 / score;
    }
    else if (lambda_score == LambdaScore::LOG)
    {
      return std::log(score);
    }
    else if (lambda_score == LambdaScore::INVERSE_LOG)
    {
      return 1.0 / std::log(score);
    }
    else if (lambda_score == LambdaScore::INVERSE_LOG10)
    {
      return 1.0 / std::log10(score);
    }
    else
    {
      throw Exception::IllegalArgument(__FILE__, __LINE__, __FUNCTION__,
        "`lambda_score`'s value is not handled by any current condition.");
    }
  }
}
