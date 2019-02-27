#ifndef _NAIVE_CONVOLVE_HPP
#define _NAIVE_CONVOLVE_HPP

#include "../Tensor/Tensor.hpp"
#include "custom_pow.hpp"

// O(n^2); use for testing and on small problems:
template <typename T>
Tensor<T> naive_convolve(const Tensor<T> & lhs, const Tensor<T> & rhs) {
  #ifdef SHAPE_CHECK
  assert(lhs.dimension() == rhs.dimension());
  assert(lhs.data_shape() + rhs.data_shape() >= 1ul);
  #endif
  if (lhs.dimension() == 0)
    return Tensor<T>();

  Tensor<T> result(lhs.data_shape() + rhs.data_shape() - 1ul);
  Vector<unsigned long> counter_result(result.dimension());

  enumerate_for_each_tensors([&counter_result, &result, &rhs](const_tup_t counter_lhs, const unsigned char dim_lhs, T lhs_val) {
      enumerate_for_each_tensors([&counter_result, &result, &rhs, &counter_lhs, &lhs_val](const_tup_t counter_rhs, const unsigned char dim_rhs, T rhs_val) {
	  for (unsigned char i=0; i<dim_rhs; ++i)
	    counter_result[i] = counter_lhs[i] + counter_rhs[i];
	  unsigned long result_flat = tuple_to_index(counter_result, result.data_shape(), dim_rhs);
	  result[result_flat] += lhs_val * rhs_val;
	},
	rhs.data_shape(),
	rhs);
    },
    lhs.data_shape(),
    lhs);

  return result;
}

template <typename T>
Tensor<T> naive_max_convolve(const Tensor<T> & lhs, const Tensor<T> & rhs) {
  #ifdef SHAPE_CHECK
  assert(lhs.dimension() == rhs.dimension());
  assert(lhs.data_shape() + rhs.data_shape() >= 1ul);
  #endif
  if (lhs.dimension() == 0)
    return Tensor<T>();

  Tensor<T> result(lhs.data_shape() + rhs.data_shape() - 1ul);
  Vector<unsigned long> counter_result(result.dimension());

  enumerate_for_each_tensors([&counter_result, &result, &rhs](const_tup_t counter_lhs, const unsigned char dim_lhs, T lhs_val) {
      enumerate_for_each_tensors([&counter_result, &result, &rhs, &counter_lhs, &lhs_val](const_tup_t counter_rhs, const unsigned char dim_rhs, T rhs_val) {
	  for (unsigned char i=0; i<dim_rhs; ++i)
	    counter_result[i] = counter_lhs[i] + counter_rhs[i];
	  unsigned long result_flat = tuple_to_index(counter_result, result.data_shape(), dim_rhs);
	  result[result_flat] = std::max(result[result_flat], lhs_val * rhs_val);
	},
	rhs.data_shape(),
	rhs);
    },
    lhs.data_shape(),
    lhs);

  return result;
}

template <typename T>
Tensor<T> naive_p_convolve(const Tensor<T> & lhs, const Tensor<T> & rhs, double p_goal) {
  #ifdef SHAPE_CHECK
  assert(lhs.dimension() == rhs.dimension());
  assert(lhs.data_shape() + rhs.data_shape() >= 1ul);
  #endif
  if (lhs.dimension() == 0)
    return Tensor<T>();

  Tensor<T> max_result(lhs.data_shape() + rhs.data_shape() - 1ul);
  Vector<unsigned long> counter_result(max_result.dimension());

  // For numeric stability, perform three passes. On the first pass,
  // perform max max-convolution (this gets the largest element of
  // each u vector). On the second pass, apply p-norm. On third pass,
  // scale by max value.
  enumerate_for_each_tensors([&counter_result, &max_result, &rhs](const_tup_t counter_lhs, const unsigned char dim_lhs, T lhs_val) {
      enumerate_for_each_tensors([&counter_result, &max_result, &rhs, &counter_lhs, &lhs_val](const_tup_t counter_rhs, const unsigned char dim_rhs, T rhs_val) {
	  for (unsigned char i=0; i<dim_rhs; ++i)
	    counter_result[i] = counter_lhs[i] + counter_rhs[i];
	  unsigned long result_flat = tuple_to_index(counter_result, max_result.data_shape(), dim_rhs);
	  max_result[result_flat] = std::max( max_result[result_flat], lhs_val*rhs_val );
	},
	rhs.data_shape(),
	rhs);
    },
    lhs.data_shape(),
    lhs);

  Tensor<T> result(max_result.data_shape());

  // Apply p-norms:
  enumerate_for_each_tensors([&counter_result, &result, &rhs, &max_result, &p_goal](const_tup_t counter_lhs, const unsigned char dim_lhs, T lhs_val) {
      enumerate_for_each_tensors([&counter_result, &result, &rhs, &counter_lhs, &lhs_val, &max_result, &p_goal](const_tup_t counter_rhs, const unsigned char dim_rhs, T rhs_val) {
	  for (unsigned char i=0; i<dim_rhs; ++i)
	    counter_result[i] = counter_lhs[i] + counter_rhs[i];
	  unsigned long result_flat = tuple_to_index(counter_result, result.data_shape(), dim_rhs);
	  // Note: using tau_denom is too conservative, but it may be
	  // good to use some numeric epsilon here.
	  if (max_result[result_flat] > 0.0)
	    result[result_flat] += custom_pow(lhs_val * rhs_val / max_result[result_flat], p_goal);
	},
	rhs.data_shape(),
	rhs);
    },
    lhs.data_shape(),
    lhs);

  for (unsigned long k=0; k<result.flat_size(); ++k)
    result[k] = custom_pow(result[k], 1.0/p_goal);

  result.flat() *= max_result.flat();

  return result;
}

#endif
