#include <map>
#ifndef _VARIABLEBOUNDS_HPP
#define _VARIABLEBOUNDS_HPP

template <typename VARIABLE_KEY> std::map<VARIABLE_KEY, std::pair<long, long> > find_bounds_from_joint(const LabeledPMF<VARIABLE_KEY> & _joint) {
  std::map<VARIABLE_KEY, std::pair<long, long> > var_to_bounds;
  for (const VARIABLE_KEY & var : _joint.ordered_variables()) {
    LabeledPMF<VARIABLE_KEY> var_marg(_joint.marginal({var}, 1));
    long var_small = var_marg.pmf().first_support()[0];
    long var_large = var_marg.pmf().last_support()[0]; 
    var_to_bounds[var] = std::make_pair(var_small, var_large);       
  }  
  return var_to_bounds;
}

template <typename VARIABLE_KEY> bool all_bounds_calculable(const AdditiveDependency<VARIABLE_KEY> additive, const std::map<VARIABLE_KEY, std::pair<long, long> > & var_to_bounds) {
  int dimensions = additive.get_output().size();
  int max_missing_bounds_in_dimension = 0;
  for (int dimension = 0; dimension < dimensions; ++dimension) {
    int missing_bounds_in_dimension = 0;
    for (const std::vector<VARIABLE_KEY> & input : additive.get_inputs()) {
      if (var_to_bounds.count(input[dimension]) == 0)
	missing_bounds_in_dimension++;
    }
    if (var_to_bounds.count(additive.get_output()[dimension]) == 0)
      missing_bounds_in_dimension++;
    if (missing_bounds_in_dimension > max_missing_bounds_in_dimension)
      max_missing_bounds_in_dimension = missing_bounds_in_dimension;
  }
  //not <= 1 b/c function used to test not only if bounds are calculable, but if they need to be calculated as well.
  return max_missing_bounds_in_dimension == 1;
}
 
template <typename VARIABLE_KEY> bool no_vars_missing_bounds(const std::map<VARIABLE_KEY, std::vector<AdditiveDependency<VARIABLE_KEY>> > & var_to_additives_using_it, const std::map<VARIABLE_KEY, std::pair<long, long> > & var_to_bounds) {
  for(auto var_iter = var_to_additives_using_it.begin(); var_iter != var_to_additives_using_it.end(); ++var_iter) {
    VARIABLE_KEY var =  var_iter->first;
    if (var_to_bounds.count(var) == 0)
      return false;    
  } 
  return true;
}

inline std::pair<long, long> compact_two_bounds (const std::pair<long, long> & new_bound, const std::pair<long, long> & old_bound) {
    
  std::pair <long, long> bound1(old_bound.first, old_bound.second);
  std::pair <long, long> bound2(new_bound.first, new_bound.second);
  std::pair <long, long> bound3(new_bound.first, old_bound.second);
  std::pair <long, long> bound4(old_bound.first, new_bound.second);

  std::pair <long, long> compacted_bound = bound1;
  if (compacted_bound.second - compacted_bound.first > bound2.second - bound2.first && bound2.second - bound2.first > 0)
    compacted_bound = bound2;
  if (compacted_bound.second - compacted_bound.first > bound3.second - bound3.first && bound3.second - bound3.first > 0)
    compacted_bound = bound3;  
  if (compacted_bound.second - compacted_bound.first > bound4.second - bound4.first && bound4.second - bound4.first > 0)
    compacted_bound = bound4;  
  return compacted_bound;
} 

template <typename VARIABLE_KEY> std::pair<long, long> find_bounds_from_additive(const VARIABLE_KEY & var, const std::vector<std::vector<VARIABLE_KEY> > & inputs, const std::vector<VARIABLE_KEY> & output, const std::map<VARIABLE_KEY, std::pair<long, long> > & var_to_bounds) {
  long low_bound = 0;
  long high_bound = 0;
  bool var_in_output = find(output.begin(), output.end(), var) != output.end();
  int var_dimension = 0;
  if (var_in_output)
    var_dimension = std::distance(output.begin(), find(output.begin(), output.end(), var));
  else {
    for (const std::vector<VARIABLE_KEY> & input : inputs) {
      auto input_it = find(input.begin(), input.end(), var);
      if (input_it != input.end())
        var_dimension = std::distance(input.begin(), input_it);
    }
  }
  for(const std::vector<VARIABLE_KEY> & input : inputs) {
    VARIABLE_KEY additive_var = input[var_dimension];
    if (additive_var != var && var_to_bounds.find(additive_var) != var_to_bounds.end()) {
      long additive_low_bound = var_to_bounds.find(additive_var)->second.first;
      long additive_high_bound = var_to_bounds.find(additive_var)->second.second;
      if (var_in_output) {
        low_bound +=  additive_low_bound;
        high_bound += additive_high_bound;
      } else {
        low_bound -= additive_high_bound;
        high_bound -= additive_low_bound;
      }
    }
  }
  if (!var_in_output && var_to_bounds.find(output[var_dimension]) != var_to_bounds.end()) {
    low_bound += var_to_bounds.find(output[var_dimension])->second.first;
    high_bound += var_to_bounds.find(output[var_dimension])->second.second;
  }
  return std::make_pair(low_bound, high_bound);
} 

template <typename VARIABLE_KEY> void add_additive_bounds(std::map<VARIABLE_KEY, std::pair<long, long> > & var_to_bounds, const std::map<VARIABLE_KEY, std::vector<AdditiveDependency<VARIABLE_KEY>> > & var_to_additives_using_it) {
  bool bounds_changed = true;
  while(bounds_changed) {
    bounds_changed = false;
    for(auto var_iter = var_to_additives_using_it.begin(); var_iter != var_to_additives_using_it.end(); ++var_iter) {
      VARIABLE_KEY var =  var_iter->first;
      std::vector<AdditiveDependency<VARIABLE_KEY> > additives = var_iter->second;
      for (const AdditiveDependency<VARIABLE_KEY> & additive : additives) {
	if (all_bounds_calculable(additive, var_to_bounds)) {
	  std::pair<long, long> new_var_bounds = find_bounds_from_additive(var, additive.get_inputs(), additive.get_output(), var_to_bounds);
	  if (var_to_bounds.count(var) == 0)
	    var_to_bounds[var] = new_var_bounds;
	  else
	    var_to_bounds[var] = compact_two_bounds(new_var_bounds, var_to_bounds[var]);     
	  bounds_changed = true;
	}
      }
    }
  }
  assert(no_vars_missing_bounds(var_to_additives_using_it, var_to_bounds));
}

#endif
