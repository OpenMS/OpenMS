#ifndef _BRUTEFORCEINFERENCEENGINE_HPP
#define _BRUTEFORCEINFERENCEENGINE_HPP

#include "TableDependency.hpp"
template <typename VARIABLE_KEY>
class BruteForceInferenceEngine;
#include "../Utility/inference_utilities.hpp"
#include "../Engine/InferenceEngine.hpp"
#include "AdditiveDependency.hpp"
#include "VariableBounds.hpp"
#include <map>

// forward declaration:
template <typename VARIABLE_KEY>
class BruteForceInferenceEngine : public InferenceEngine<VARIABLE_KEY> {
protected:
  LabeledPMF<VARIABLE_KEY> _joint;
  const double _p;

private:
  void multiply_in_additives(const std::vector<AdditiveDependency<VARIABLE_KEY> > & all_additive) {
    std::map<VARIABLE_KEY, std::vector<AdditiveDependency<VARIABLE_KEY>> > var_to_additives_using_it;
    build_variable_to_additives(all_additive, var_to_additives_using_it);
    std::map<VARIABLE_KEY, std::pair<long, long> > var_to_bounds = find_bounds_from_joint(_joint);
    add_additive_bounds(var_to_bounds, var_to_additives_using_it);
    for (const AdditiveDependency<VARIABLE_KEY> & additive : all_additive) {
      LabeledPMF<VARIABLE_KEY> lpmf = to_lpmf(additive, var_to_bounds); 

      _joint = _joint * lpmf;
    }
  }

  void add_to_map(const VARIABLE_KEY & var, const AdditiveDependency<VARIABLE_KEY> & additive, std::map<VARIABLE_KEY, std::vector<AdditiveDependency<VARIABLE_KEY>> > & var_to_additives_using_it) const {
    auto additive_iter = var_to_additives_using_it.find(var);
    if (additive_iter != var_to_additives_using_it.end())
      additive_iter->second.push_back(additive);         
    else {
      // first time this variable has been seen in an additive dependency:            
      var_to_additives_using_it.insert(std::pair<VARIABLE_KEY, std::vector<AdditiveDependency<VARIABLE_KEY>> >(var, {additive})); 
    }  
  }

  void build_variable_to_additives(const std::vector<AdditiveDependency<VARIABLE_KEY> > & all_additive, std::map<VARIABLE_KEY, std::vector<AdditiveDependency<VARIABLE_KEY>> > & var_to_additives_using_it) const {
    for (const AdditiveDependency<VARIABLE_KEY> & additive : all_additive) {
      for (const std::vector<VARIABLE_KEY> & vect : additive.get_inputs())
        for (const VARIABLE_KEY & var : vect)
	        add_to_map(var, additive, var_to_additives_using_it);
      for (const VARIABLE_KEY & var : additive.get_output())
      	add_to_map(var, additive, var_to_additives_using_it);
    }
  }

  LabeledPMF<VARIABLE_KEY> to_lpmf(const AdditiveDependency<VARIABLE_KEY> & additive, const std::map<VARIABLE_KEY, std::pair<long, long> > & var_to_bounds) const {
    std::vector<std::vector<VARIABLE_KEY> > inputs = additive.get_inputs();
    inputs.push_back(additive.get_output());
    std::vector<VARIABLE_KEY> flattened_inputs = additive.get_all_variables_used();
    int num_dimensions = inputs[0].size();
    LabeledPMF<VARIABLE_KEY> joint;
    std::vector<long> first_support;
    std::vector<unsigned long> result_table_dims;
    for(const VARIABLE_KEY & var : flattened_inputs) {
      std::pair<long, long> bounds;
      if (var_to_bounds.find(var) != var_to_bounds.end())
	      bounds = var_to_bounds.find(var)->second;
      long var_small = bounds.first;
      long var_large = bounds.second;
      first_support.push_back(var_small);
      result_table_dims.push_back(var_large - var_small + 1);
    }
    // Initialized with zeros:
    Tensor<double> result_table(result_table_dims);
    enumerate_apply_tensors([&first_support, num_dimensions](const_tup_t index, const int dim, double & res_val){
      std::vector<long> sum_val(num_dimensions, 0);
      for (int current_dim = 0; current_dim < num_dimensions; ++current_dim)
        for(int i = current_dim; i < dim-num_dimensions; i += num_dimensions)
          sum_val[current_dim] += index[i] + first_support [i];
      bool is_additive = true;        
      for (int i = 0; i < num_dimensions; ++i)
        if (sum_val[i] != long(index[dim-num_dimensions+i] + first_support[dim-num_dimensions+i]))
          is_additive = false;
      if (is_additive)
        res_val = 1.0;
      },
    result_table.data_shape(),
    result_table);
    LabeledPMF<VARIABLE_KEY> table_dependency(flattened_inputs, PMF(first_support, result_table));
    table_dependency.reset_log_normalization_constant();
    return table_dependency;
  }

public:
  BruteForceInferenceEngine(const std::vector<TableDependency<VARIABLE_KEY> > & all_tables, const std::vector<AdditiveDependency<VARIABLE_KEY> > & all_additive, double p):
    _p(p)
  {
    // Note: instead of !=, could check fabs(table.p - _p) > epsilon
    for (const TableDependency<VARIABLE_KEY> & table : all_tables)
      if (table.p != _p) {
      	std::cerr << "Cannot do brute force on non-homogeneous p norms" << std::endl;
	      assert(false);
      }
    for (const AdditiveDependency<VARIABLE_KEY> & additive : all_additive)
      if (additive.p != _p) {
	      std::cerr << "Cannot do brute force on non-homogeneous p norms" << std::endl;
	      assert(false);
      }

    // Note: This could be performed more efficiently by first
    // allocating the result table (using the union of all variables
    // and the intersection of their supports) and then performing a
    // single product over all distributions simultaneously.
    for (const TableDependency<VARIABLE_KEY> & table : all_tables)
      _joint = _joint * table.labeled_pmf();

    multiply_in_additives(all_additive);
  }

  std::vector<LabeledPMF<VARIABLE_KEY> > estimate_posteriors(const std::vector<std::vector<VARIABLE_KEY> > & joint_distributions_to_retrieve) {
    std::vector<LabeledPMF<VARIABLE_KEY> > results;
    for (const std::vector<VARIABLE_KEY> & ordered_vars : joint_distributions_to_retrieve)
      results.push_back(_joint.marginal(ordered_vars, _p));
    return results;
  }

  // Note: uses p=1:
  double log_normalization_constant() {
    return _joint.log_normalization_constant();
  }
};

#endif
