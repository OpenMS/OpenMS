#ifndef _LABELEDPMF_H
#define _LABELEDPMF_H

#include <unordered_map>
#include <vector>

#include "PMF.hpp"
#include "semi_outer_product_and_quotient.hpp"
#include "divergence.hpp"
#include "dampen.hpp"

template <typename VARIABLE_KEY>
class LabeledPMF {
protected:
  std::vector<VARIABLE_KEY> _ordered_variables;
  std::unordered_map<VARIABLE_KEY, unsigned char> _variable_to_index;
  PMF _pmf;

  void construct_var_to_index() {
    for (unsigned char i=0; i<_ordered_variables.size(); ++i) {
      #ifdef SHAPE_CHECK
      auto iter = _variable_to_index.find(_ordered_variables[i]);
      // The ordered variables must be unique:
      assert(iter == _variable_to_index.end() );
      #endif
      _variable_to_index[_ordered_variables[i]] = i;
    }
  }
    
public:
  LabeledPMF() { }

  LabeledPMF(const std::vector<VARIABLE_KEY> & ordered_variables, const PMF & pmf_param):
    _ordered_variables(ordered_variables),
    _pmf(pmf_param)
  {
    #ifdef SHAPE_CHECK
    assert(_ordered_variables.size() == _pmf.dimension());
    #endif
    
    construct_var_to_index();
  }

  LabeledPMF(const std::vector<VARIABLE_KEY> & ordered_variables, PMF && pmf_param):
    _ordered_variables(ordered_variables),
    _pmf(pmf_param)
  {
    #ifdef SHAPE_CHECK
    assert(_ordered_variables.size() == _pmf.dimension());
    #endif
    
    construct_var_to_index();
  }

  LabeledPMF(std::vector<VARIABLE_KEY> && ordered_variables, PMF && pmf_param):
    _ordered_variables(ordered_variables),
    _pmf(pmf_param)
  {
    #ifdef SHAPE_CHECK
    assert(_ordered_variables.size() == _pmf.dimension());
    #endif
    
    construct_var_to_index();
  }

  const PMF & pmf() const {
    return _pmf;
  }

  const unsigned char dimension() const {
    return _pmf.dimension();
  }

  const double log_normalization_constant() const {
    return _pmf.log_normalization_constant();
  }

  void add_to_log_normalization_constant(const double log_c) {
    _pmf.add_to_log_normalization_constant(log_c);
  }

  void reset_log_normalization_constant() {
      _pmf.reset_norm_constant();
  }

  const std::vector<VARIABLE_KEY> & ordered_variables() const {
    return _ordered_variables;
  }
  
  LabeledPMF marginal(const std::vector<VARIABLE_KEY> & ordered_vars_to_keep, double p) const {
    Vector<unsigned char> indices = lookup_indices(ordered_vars_to_keep);
    #ifdef SHAPE_CHECK
    verify_subpermutation(indices, dimension());
    #endif

    // When the vars are just a permutation (none are eliminated),
    // then simply transpose (this can be more efficient):
    if (ordered_vars_to_keep.size() == dimension())
      return transposed(ordered_vars_to_keep);
    
    return LabeledPMF(ordered_vars_to_keep, _pmf.marginal(indices,p));
  }

  LabeledPMF transposed(const Vector<unsigned char> & new_axis_order) const {
    #ifdef BOUNDSCHECK
    assert(new_variable_order.size() == dimension());
    verify_permutation(new_axis_order);
    #endif
    std::vector<VARIABLE_KEY> new_variable_order(dimension());
    for (unsigned char i=0; i<dimension(); ++i)
      new_variable_order[i] = _ordered_variables[ new_axis_order[i] ];
    return LabeledPMF(new_variable_order, _pmf.transposed(new_axis_order));
  }
  LabeledPMF transposed(const std::vector<VARIABLE_KEY> & new_variable_order) const {
    Vector<unsigned char> new_axis_order = lookup_indices(new_variable_order);
    // Note: there is code shared with the function above, but write
    // it from scratch because new_variable_order does not need to be
    // constructed for this version.
    #ifdef BOUNDSCHECK
    assert(new_variable_order.size() == dimension());
    verify_permutation(new_axis_order);
    #endif
    return LabeledPMF(new_variable_order, _pmf.transposed(new_axis_order));
  }
  void transpose(const std::vector<VARIABLE_KEY> & new_variable_order) {
    if (new_variable_order == _ordered_variables)
      return;
    
    Vector<unsigned char> new_axis_order = lookup_indices(new_variable_order);
    // Note: there is code shared with the function above, but write
    // it from scratch because new_variable_order does not need to be
    // constructed for this version.
    #ifdef BOUNDSCHECK
    assert(new_variable_order.size() == dimension());
    verify_permutation(new_axis_order);
    #endif

    _ordered_variables = new_variable_order;
    _pmf.transpose(new_axis_order);
  }

  // To avoid building sets (since the _variable_to_index table is
  // already built):
  int variable_index(const VARIABLE_KEY & var) const {
    auto iter = _variable_to_index.find(var);
    if (iter != _variable_to_index.end())
      return iter->second;
    return -1;
  }

  Vector<unsigned char> lookup_indices(const std::vector<VARIABLE_KEY> & vars) const {
    Vector<unsigned char> res(vars.size());
    for (unsigned char i=0; i<vars.size(); ++i) {
      auto iter = _variable_to_index.find(vars[i]);
      #ifdef SHAPE_CHECK
      assert(iter != _variable_to_index.end() && "Variable not found in LabeledPMF");
      #endif
      res[i] = iter->second;
    }

    #ifdef SHAPE_CHECK
    verify_subpermutation(res, dimension());
    #endif

    return res;
  }

  bool contains_variable(const VARIABLE_KEY & var) const {
    return variable_index(var) != -1;
  }

  // Helper function for efficiently computing products, quotients,
  // and other tasks that require aligned tables to be intersected:
  std::pair<TensorView<double>, Vector<long> > view_of_intersection_with(const LabeledPMF & rhs) const {
    // Uses existing dictionaries to avoid creating a set on the fly:
    unsigned char intersection_size = 0;
    for (unsigned char i=0; i<dimension(); ++i) {
      const VARIABLE_KEY & var = ordered_variables()[i];
      if (rhs.contains_variable(var))
	++intersection_size;
    }

    Vector<long> first_sup = _pmf.first_support();
    Vector<long> view_shape(dimension());
    for (unsigned char i=0; i<dimension(); ++i) {
      const VARIABLE_KEY & var = ordered_variables()[i];
      int index_rhs = rhs.variable_index(var);

      // Compute the intersection over the minumum supports:
      if (index_rhs != -1)
	first_sup[i] = std::max(first_sup[i], rhs._pmf.first_support()[index_rhs]);

      // Compute the intersection over the maximum supports (as a
      // shape):
      const long max_sup_plus_one = _pmf.first_support()[i] + _pmf.table().data_shape()[i];

      view_shape[i] = max_sup_plus_one;
      if (index_rhs != -1) {
	const long rhs_max_sup_plus_one = rhs._pmf.first_support()[index_rhs] + rhs._pmf.table().data_shape()[index_rhs];
	view_shape[i] = std::min(view_shape[i], rhs_max_sup_plus_one);
      }

      #ifdef SHAPE_CHECK
      if (view_shape[i] < first_sup[i]) {
	std::cerr << "Error: narrowing LabeledPMF would produce empty LabeledPMF" << std::endl;
	assert(false);
      }
      #endif
      
      view_shape[i] -= first_sup[i];
    }

    return std::make_pair(_pmf.table().start_at_const(first_sup - _pmf.first_support(), view_shape), first_sup);
  }

  bool has_same_variables(const LabeledPMF & rhs) const {
    // Checks that the variable sets are identical:
    for (unsigned char i=0; i<dimension(); ++i) {
      const VARIABLE_KEY & var = ordered_variables()[i];
      if ( ! rhs.contains_variable(var) )
	return false;
    }
    for (unsigned char i=0; i<rhs.dimension(); ++i) {
      const VARIABLE_KEY & var = rhs.ordered_variables()[i];
      if ( ! contains_variable(var) )
	return false;
    }
    return true;
  }
};

template <typename VARIABLE_KEY, bool MULT> // true is means multiply, false means divide
LabeledPMF<VARIABLE_KEY> mult_or_div(const LabeledPMF<VARIABLE_KEY> & lhs, const LabeledPMF<VARIABLE_KEY> & rhs) {
  #ifdef SHAPE_CHECK
  // Check that bounds intersect:
  for (unsigned int lhs_index=0; lhs_index<lhs.ordered_variables().size(); ++lhs_index) {
    const VARIABLE_KEY & var = lhs.ordered_variables()[lhs_index];
    int rhs_index = rhs.variable_index(var);
    if (rhs_index != -1) {
      // Variable is in both lhs and rhs:
      long min_lhs_outcome = lhs.pmf().first_support()[lhs_index];
      long max_lhs_outcome = min_lhs_outcome + lhs.pmf().table().view_shape()[lhs_index] - 1;
      long min_rhs_outcome = rhs.pmf().first_support()[rhs_index];
      long max_rhs_outcome = min_rhs_outcome + rhs.pmf().table().view_shape()[rhs_index] - 1;

      assert( ((min_rhs_outcome <= max_lhs_outcome && max_rhs_outcome >= min_lhs_outcome) || (min_lhs_outcome <= max_rhs_outcome && max_lhs_outcome >= min_rhs_outcome)) && "Error: multiplying LabeledPMFs would produce empty product");
    }
  }
  #endif

  std::pair<TensorView<double>, Vector<long> > lhs_view_and_first_sup = lhs.view_of_intersection_with(rhs);
  std::pair<TensorView<double>, Vector<long> > rhs_view_and_first_sup = rhs.view_of_intersection_with(lhs);

  unsigned char intersection_size=0;
    
  // Check if the shared variables are already in the same order and
  // in the inner-most indices; in that case, transposition is
  // unnecessary. Simultaneously compute the number of shared
  // variables.
  int last_shared_index = -1;
  bool already_in_order = true;
  int rhs_index = -1;

  for (unsigned char i=0; i<lhs.dimension(); ++i) {
    const VARIABLE_KEY & var = lhs.ordered_variables()[i];

    rhs_index = rhs.variable_index(var);

    if (rhs_index != -1) {
      // Unrelated to other code here; these loops are fused to
      // compute intersection_size without performing variable
      // lookup again.
      ++intersection_size;
      
      if (last_shared_index != -1 && last_shared_index != rhs_index-1) {
	// A block of consecutive in-order intersecting variables was
	// broken by a non-consecutive shared variable.
	already_in_order = false;
      }

      last_shared_index = rhs_index;
    }
    else {
      if (last_shared_index != -1) {
	// A block of consecutive in-order intersecting variables was
	// broken by a non-intersecting shared variable. This is valid
	// because shared variables must be inner-most, so any time an
	// shared variable is found, the rest of the variables must be
	// shared).
	already_in_order = false;
      }
    }
  }

  // To be in order, a block must be consecutive and the final
  // variable in that block (i.e., the final variable in lhs) must be
  // the final variable in rhs:
  already_in_order = already_in_order && rhs_index+1 == rhs.dimension();

  const unsigned char unique_lhs_dims=lhs.dimension()-intersection_size;
  const unsigned char unique_rhs_dims=rhs.dimension()-intersection_size;

  std::vector<VARIABLE_KEY> new_variable_order;
  // First, insert variables unique to lhs:
  for (unsigned char i=0; i<lhs.dimension(); ++i) {
    const VARIABLE_KEY & var = lhs.ordered_variables()[i];
    if ( ! rhs.contains_variable(var))
      new_variable_order.push_back(var);
  }
  // Second, insert variables unique to rhs:
  for (unsigned char i=0; i<rhs.dimension(); ++i) {
    const VARIABLE_KEY & var = rhs.ordered_variables()[i];
    if ( ! lhs.contains_variable(var))
      new_variable_order.push_back(var);
  }
  // Lastly, insert variables shared (use order of lhs):
  for (unsigned char i=0; i<lhs.dimension(); ++i) {
    const VARIABLE_KEY & var = lhs.ordered_variables()[i];
    if ( rhs.contains_variable(var))
      new_variable_order.push_back(var);
  }
    
  Vector<long> new_first_support(lhs.dimension() + rhs.dimension() - intersection_size);
  
  if (already_in_order) {
    // Work directly with the tensor views:

    // Compute the first support of the result distribution:
    for (unsigned char i=0; i<unique_lhs_dims; ++i)
      new_first_support[i] = lhs_view_and_first_sup.second[i];
    for (unsigned char i=0; i<unique_rhs_dims; ++i)
      new_first_support[unique_lhs_dims + i] = rhs_view_and_first_sup.second[i];
    for (unsigned char i=0; i<intersection_size; ++i)
      new_first_support[unique_lhs_dims + unique_rhs_dims + i] = lhs_view_and_first_sup.second[unique_lhs_dims + i];

    if (MULT) {
      PMF res( new_first_support, semi_outer_product(lhs_view_and_first_sup.first, rhs_view_and_first_sup.first, intersection_size) );
      res.add_to_log_normalization_constant( lhs.log_normalization_constant() + rhs.log_normalization_constant() );
      return LabeledPMF<VARIABLE_KEY>( new_variable_order, res );
    }
    else {
      PMF res(new_first_support, semi_outer_quotient(lhs_view_and_first_sup.first, rhs_view_and_first_sup.first, intersection_size));
      res.add_to_log_normalization_constant( lhs.log_normalization_constant() - rhs.log_normalization_constant() );
      return LabeledPMF<VARIABLE_KEY>( new_variable_order, res);
    }
  }
  else {
    // Transpose into the proper order:
    Tensor<double> lhs_part(lhs_view_and_first_sup.first);
    Tensor<double> rhs_part(rhs_view_and_first_sup.first);

    Vector<unsigned char> new_lhs_order(lhs.dimension());
    for (unsigned char i=0; i<unique_lhs_dims; ++i)
      new_lhs_order[i] = lhs.variable_index(new_variable_order[i]);
    for (unsigned char i=0; i<intersection_size; ++i)
      new_lhs_order[unique_lhs_dims+i] = lhs.variable_index(new_variable_order[unique_lhs_dims+unique_rhs_dims+i]);

    Vector<unsigned char> new_rhs_order(rhs.dimension());
    for (unsigned char i=0; i<unique_rhs_dims; ++i)
      new_rhs_order[i] = rhs.variable_index(new_variable_order[unique_lhs_dims + i]);
    for (unsigned char i=0; i<intersection_size; ++i)
      new_rhs_order[unique_rhs_dims+i] = rhs.variable_index(new_variable_order[unique_lhs_dims+unique_rhs_dims+i]);

    evergreen::transpose(lhs_part, new_lhs_order);
    evergreen::transpose(rhs_part, new_rhs_order);

    // Compute the first support of the result distribution:
    for (unsigned char i=0; i<unique_lhs_dims; ++i) {
      unsigned char new_index = new_lhs_order[i];
      new_first_support[i] = lhs_view_and_first_sup.second[new_index];
    }
    for (unsigned char i=0; i<unique_rhs_dims; ++i) {
      unsigned char new_index = new_rhs_order[i];
      new_first_support[unique_lhs_dims+i] = rhs_view_and_first_sup.second[new_index];
    }
    for (unsigned char i=0; i<intersection_size; ++i) {
      unsigned char new_index = new_lhs_order[unique_lhs_dims+i];
      new_first_support[unique_lhs_dims+unique_rhs_dims+i] = lhs_view_and_first_sup.second[new_index];
    }

    if (MULT) {
      PMF res( new_first_support, semi_outer_product(lhs_part, rhs_part, intersection_size));
      res.add_to_log_normalization_constant(lhs.log_normalization_constant() + rhs.log_normalization_constant());
      return LabeledPMF<VARIABLE_KEY>(new_variable_order, res);
    }
    else {
      PMF res(new_first_support, semi_outer_quotient(lhs_part, rhs_part, intersection_size));
      res.add_to_log_normalization_constant(lhs.log_normalization_constant() - rhs.log_normalization_constant());
      return LabeledPMF<VARIABLE_KEY>(new_variable_order, res);
    }
  }
}

template <typename VARIABLE_KEY>
LabeledPMF<VARIABLE_KEY> operator *(const LabeledPMF<VARIABLE_KEY> & lhs, const LabeledPMF<VARIABLE_KEY> & rhs) {
  if (rhs.dimension() == 0)
    return lhs;
  if (lhs.dimension() == 0)
    return rhs;

  return mult_or_div<VARIABLE_KEY, true>(lhs, rhs);
}

template <typename VARIABLE_KEY>
LabeledPMF<VARIABLE_KEY> operator /(const LabeledPMF<VARIABLE_KEY> & lhs, const LabeledPMF<VARIABLE_KEY> & rhs) {
  #ifdef SHAPE_CHECK
  // Dividing an empty LabeledPMF by a non-empty LabeledPMF makes no
  // sense (unless you were to make every element 1.0/old_value; but
  // it is unclear when that would ever be useful):
  if (rhs.dimension() > 0)
    assert(lhs.dimension() > 0);
  #endif

  if (rhs.dimension() == 0)
    return lhs;

  return mult_or_div<VARIABLE_KEY, false>(lhs, rhs);
}

template <typename VARIABLE_KEY>
std::ostream & operator << (std::ostream & os, const LabeledPMF<VARIABLE_KEY> & rhs) {
  for (unsigned char i=0; i<rhs.dimension(); ++i)
    os << rhs.ordered_variables()[i] << " ";
  os << rhs.pmf();
  return os;
}

#endif
