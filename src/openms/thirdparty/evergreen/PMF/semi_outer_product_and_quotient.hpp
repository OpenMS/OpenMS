#ifndef _SEMI_OUTER_PRODUCT_AND_QUOTIENT_HPP
#define _SEMI_OUTER_PRODUCT_AND_QUOTIENT_HPP

// For performing semi_outer_... functions:
template <typename FUNCTION, template <typename> class TENSOR>
Tensor<double> semi_outer_apply(const TensorLike<double, TENSOR> & lhs, const TensorLike<double, TENSOR> & rhs, const unsigned char overlapping_inner_dims, FUNCTION semi_outer_function) {
  #ifdef SHAPE_CHECK
  assert(lhs.dimension() > 0 && rhs.dimension() > 0);
  #endif

  // semi_outer_function is either or semi_outer_product or semi_outer_quotient

  const unsigned char unique_lhs_dims = lhs.dimension() - overlapping_inner_dims;
  const unsigned char unique_rhs_dims = rhs.dimension() - overlapping_inner_dims;
  
  Vector<unsigned long> outer_shape_lhs = lhs.view_shape().start_at_const(0, unique_lhs_dims);
  Vector<unsigned long> outer_shape_rhs = rhs.view_shape().start_at_const(0, unique_rhs_dims);
  Vector<unsigned long> inner_shape_lhs = lhs.view_shape().start_at_const(unique_lhs_dims, overlapping_inner_dims);
  Vector<unsigned long> inner_shape_rhs = rhs.view_shape().start_at_const(unique_rhs_dims, overlapping_inner_dims);

  Vector<unsigned long> result_shape = concatenate(concatenate(outer_shape_lhs, outer_shape_rhs), inner_shape_lhs);

  #ifdef SHAPE_CHECK
  assert( lhs.dimension() >= overlapping_inner_dims );
  assert( rhs.dimension() >= overlapping_inner_dims );

  // Inner shapes must match:
  assert(inner_shape_lhs == inner_shape_rhs);
  #endif

  Tensor<double> result( result_shape );

  if (unique_lhs_dims > 0 || unique_rhs_dims > 0) {
    Vector<unsigned long> counter_lhs(lhs.dimension());
    Vector<unsigned long> counter_rhs(rhs.dimension());
    enumerate_apply_tensors([&counter_lhs, &counter_rhs, &lhs, &rhs, unique_lhs_dims, unique_rhs_dims, overlapping_inner_dims, &semi_outer_function](const_tup_t counter_result, const unsigned char result_dims, double & res_val) {
	// Note: This could be optimized to not use the counter Vectors:
	for (unsigned char i=0; i<unique_lhs_dims; ++i)
	  counter_lhs[i] = counter_result[i];
	for (unsigned char i=0; i<overlapping_inner_dims; ++i)
	  counter_lhs[unique_lhs_dims+i] = counter_result[unique_lhs_dims+unique_rhs_dims+i];
	
	for (unsigned char i=0; i<unique_rhs_dims; ++i)
	  counter_rhs[i] = counter_result[unique_lhs_dims+i];
	for (unsigned char i=0; i<overlapping_inner_dims; ++i)
	  counter_rhs[unique_rhs_dims+i] = counter_result[unique_lhs_dims+unique_rhs_dims+i];

	res_val = semi_outer_function(lhs[counter_lhs], rhs[counter_rhs]);
      },
      result.data_shape(),
      result);
  }
  else // unique_lhs_dims == 0 && unique_rhs_dims == 0 (compute element-wise product):
    apply_tensors([&semi_outer_function](double & res_val, double lhs_val, double rhs_val) {
	res_val = semi_outer_function(lhs_val, rhs_val);
      },
      result.data_shape(),
      result, lhs, rhs);

  return result;
}

template <template <typename> class TENSOR>
Tensor<double> semi_outer_product(const TensorLike<double, TENSOR> & lhs, const TensorLike<double, TENSOR> & rhs, const unsigned char overlapping_inner_dims) {
  return semi_outer_apply(lhs, rhs, overlapping_inner_dims, [](double x, double y) {
      return x * y;
    });
}

template <template <typename> class TENSOR>
Tensor<double> semi_outer_quotient(const TensorLike<double, TENSOR> & lhs, const TensorLike<double, TENSOR> & rhs, const unsigned char overlapping_inner_dims) {
  return semi_outer_apply(lhs, rhs, overlapping_inner_dims, [](double x, double y) {
      // Note: fabs not necessary for probabilistic problems (since
      // PMFs are >=0), but it's better to be tidy:
      if ( fabs(y) > tau_denom )
	return x / y;
      return 0.0;
    });
}

#endif
