#ifndef _DAMPEN_HPP
#define _DAMPEN_HPP

// For computing a convex combination of two LabeledPMFs (using only the
// intersecting support):
template <typename VARIABLE_KEY>
LabeledPMF<VARIABLE_KEY> dampen(const LabeledPMF<VARIABLE_KEY> & lhs, const LabeledPMF<VARIABLE_KEY> & rhs, double lambda) {
  #ifdef SHAPE_CHECK
  assert(lhs.has_same_variables(rhs));
  #endif
  #ifdef NUMERIC_CHECK
  assert(lambda >= 0 && lambda <= 1);
  #endif

  // It is important to call this in the consistent order (lhs,
  // rhs) so that lambda and 1-lambda are multiplied with the
  // appropriate respective values:
  auto convex_combination = [lambda](double a, double b) {
    return lambda*a + (1-lambda)*b;
  };

  std::pair<TensorView<double>, Vector<long> > lhs_view_and_first_sup = lhs.view_of_intersection_with(rhs);
  std::pair<TensorView<double>, Vector<long> > rhs_view_and_first_sup = rhs.view_of_intersection_with(lhs);

  const TensorView<double> & lhs_view = lhs_view_and_first_sup.first;
  const TensorView<double> & rhs_view = rhs_view_and_first_sup.first;
  Vector<long> & first_support = lhs_view_and_first_sup.second;
  
  if (lhs.ordered_variables() == rhs.ordered_variables()) {
    // variables are in the same order; no need to transpose:
    Tensor<double> res_table(lhs_view);

    apply_tensors([&convex_combination](double & res_val, double rhs_val){
	res_val = convex_combination(res_val, rhs_val);
      },
      res_table.data_shape(),
      res_table, rhs_view);

    PMF pmf(first_support, std::move(res_table));
    return LabeledPMF<VARIABLE_KEY>(lhs.ordered_variables(), std::move(pmf));
  }
  else {
    // transpose rhs to get variables in the same order:
    Tensor<double> res_table(lhs_view);
    Vector<unsigned int> new_rhs_order = rhs.lookup_indices(lhs.ordered_variables());
    transpose(res_table, new_rhs_order);

    apply_tensors([&convex_combination](double & res_val, double rhs_val){
	res_val = convex_combination(res_val, rhs_val);
      },
      res_table.data_shape(),
      res_table, rhs_view);

    PMF pmf(first_support, std::move(res_table));
    return LabeledPMF<VARIABLE_KEY>(lhs.ordered_variables(), std::move(pmf));
  }
}

#endif
