#ifndef _DIVERGENCE_HPP
#define _DIVERGENCE_HPP

#include "squared.hpp"

template <template <typename> class TENSOR_LHS, template <typename> class TENSOR_RHS>
double se(const TensorLike<double, TENSOR_LHS> & lhs, const TensorLike<double, TENSOR_RHS> & rhs) {
  #ifdef SHAPE_CHECK
  assert( lhs.view_shape() == rhs.view_shape() );
  #endif
  
  double tot = 0.0;

  for_each_tensors([&tot](double lhs_val, double rhs_val){
      tot += squared(lhs_val - rhs_val);
    },
    lhs.view_shape(),
    lhs, rhs);

  return tot;
}

template <typename VARIABLE_KEY>
class LabeledPMF;

template <typename VARIABLE_KEY>
double mse_divergence(const LabeledPMF<VARIABLE_KEY> & lhs, const LabeledPMF<VARIABLE_KEY> & rhs) {
  #ifdef SHAPE_CHECK
  assert(lhs.has_same_variables(rhs));
  #endif
  
  std::pair<TensorView<double>, Vector<long> > lhs_view_and_first_sup = lhs.view_of_intersection_with(rhs);
  std::pair<TensorView<double>, Vector<long> > rhs_view_and_first_sup = rhs.view_of_intersection_with(lhs);

  const TensorView<double> & lhs_view = lhs_view_and_first_sup.first;
  const TensorView<double> & rhs_view = rhs_view_and_first_sup.first;
  
  double lhs_view_mass = 0.0;
  for_each_tensors([&lhs_view_mass](double val){
      lhs_view_mass += val;
    },
    lhs_view.view_shape(),
    lhs_view
    );
  double rhs_view_mass = 0.0;
  for_each_tensors([&rhs_view_mass](double val){
      rhs_view_mass += val;
    },
    rhs_view.view_shape(),
    rhs_view
    );
  
  double nonintersecting_se = squared(1.0 - lhs_view_mass) + squared(1.0 - rhs_view_mass);

  double intersecting_se;
  if (lhs.ordered_variables() == rhs.ordered_variables()) {
    // variables are in the same order; no need to transpose:
    intersecting_se = se(lhs_view, rhs_view);
  }
  else {
    // transpose rhs to get variables in the same order:
    Tensor<double> rhs_part(rhs_view);
    Vector<unsigned int> new_rhs_order = rhs.lookup_indices(lhs.ordered_variables());
    transpose(rhs_part, new_rhs_order);
    
    intersecting_se = se(lhs_view_and_first_sup.first, rhs_part);
  }

  // Note: lhs_view.flat_size() == rhs_view.flat_size()
  return (nonintersecting_se + intersecting_se) / (lhs.pmf().table().flat_size() + rhs.pmf().table().flat_size() - lhs_view.flat_size());
}

#endif
