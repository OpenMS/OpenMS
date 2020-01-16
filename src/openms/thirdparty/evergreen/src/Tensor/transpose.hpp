#ifndef _TRANSPOSE_HPP
#define _TRANSPOSE_HPP

#include "MatrixTranspose.hpp"

template <typename T>
class Tensor;

// Empirically chosen:
const unsigned long SIZE_WHERE_NAIVE_TRANSPOSE_BECOMES_SLOWER = 8;

template <typename T>
inline Tensor<T> naive_transposed(const Tensor<T> & ten, const Vector<unsigned char> & new_axis_order) {
#ifdef SHAPE_CHECK
  assert(ten.dimension() == new_axis_order.size());
  verify_permutation(new_axis_order);
#endif

  Vector<unsigned long> new_shape(ten.dimension());
  for (unsigned char i=0; i<ten.dimension(); ++i)
    new_shape[i] = ten.data_shape()[ new_axis_order[i] ];
  Tensor<T> result(new_shape);

  Vector<unsigned long> reordered_tup(ten.dimension());
  enumerate_for_each_tensors([&result, &reordered_tup, &new_axis_order](const_tup_t tup, const unsigned char dim, const T & val){
      for (unsigned char i=0; i<dim; ++i)
	reordered_tup[i] = tup[ new_axis_order[i] ];
      result[ tuple_to_index(reordered_tup, result.data_shape(), dim) ] = val;
    },
    ten.data_shape(),
    ten);

  return result;
}

template <typename T>
void naive_transpose(Tensor<T> & ten, const Vector<unsigned char> & new_axis_order) {
  ten = std::move(naive_transposed(ten, new_axis_order));
}

// Cache friendly version:
template <typename T>
void cache_friendly_transpose(Tensor<T> & ten, const Vector<unsigned char> & new_axis_order) {
#ifdef SHAPE_CHECK
  assert(ten.dimension() == new_axis_order.size());
  verify_permutation(new_axis_order);
#endif

  // For performance: when the prefix of the new order is already
  // partially in order, do not visit those indices.
  unsigned char already_ordered_prefix;
  for (already_ordered_prefix=0; already_ordered_prefix<new_axis_order.size(); ++already_ordered_prefix)
    if ( new_axis_order[already_ordered_prefix] != already_ordered_prefix)
      break;

  if (already_ordered_prefix < ten.dimension()) {
    T* __restrict buffer_from = &ten.flat()[0];
    Tensor<T> buffer(ten.data_shape());
    T* __restrict buffer_to = &buffer.flat()[0];

    /*
      This function performs tensor transposition in O(d) matrix
      transpositions. This allows it to use the cache-oblivious matrix
      transposition, and therefore be more performant..

      (a,b,c,d,e,f,g) @
      (3,1,5,0,6,4,2) -->
      (d,b,c,f,e,g,a)
    
      2D contiguous transposes only swap adjacent inner
      indices. Therefore, you can send a given index to the far right
      and shift the others left:
    
      (0,1,2,3,4,5,6) -->
      (0,1,2,4,5,6,3) -->
      (0,2,4,5,6,3,1) -->
      (0,2,4,6,3,1,5) -->
      (2,4,6,3,1,5,0) -->
      (2,4,3,1,5,0,6) -->
      (2,3,1,5,0,6,4) -->
      (3,1,5,0,6,4,2)
    */

    // Note: this method is in O( N d + d^2 ). It could likely be done
    // in O( N d + d log(d) ) or better, but O( N d + d^2 ) = O( d ( N +
    // d ) ), which will equal O( N d ) when d is in O(N). When d is not
    // in O(N), some axes must have a length of 1. Therefore, these axes
    // could also alternatively be suppressed, since they will not
    // effect the result. If a o(d^2) solution is necessary in general,
    // this is likely preferred over using a tree / map to somehow
    // reduce the d^2 --> d log(d).
    Vector<unsigned char> current_axis_order = seq<unsigned char>(ten.dimension());
    for (unsigned char i=already_ordered_prefix; i<ten.dimension(); ++i) {
      unsigned char next_axis = new_axis_order[i];

      unsigned char next_axis_index = 0;
      for (next_axis_index=0; next_axis_index<ten.dimension(); ++next_axis_index)
	if (current_axis_order[next_axis_index] == next_axis)
	  break;
      
      // Compute number of 2D transposes and the R,C:
      unsigned long number_of_2d_transposes = 1;
      for (unsigned char j=0; j<next_axis_index; ++j)
	number_of_2d_transposes *= ten.data_shape()[ current_axis_order[j] ];
      unsigned long R = ten.data_shape()[ current_axis_order[next_axis_index] ];
      unsigned long C = 1;
      for (unsigned char j=next_axis_index+1; j<ten.dimension(); ++j)
	C *= ten.data_shape()[ current_axis_order[j] ];

      // Note that this could be sped up by swapping in the largest
      // blocks possible (e.g., a contig of the first indices may
      // already be in order).

      // Perform 2D transposes if non-trivial:
      if (R > 1 && C > 1) {
	for (unsigned long j=0; j<number_of_2d_transposes; ++j)
	  MatrixTranspose<T>::apply_buffered(buffer_to + j*R*C, buffer_from + j*R*C, R, C);
	// The following does not change the vector pointer inside source:
	std::swap(buffer_from, buffer_to);
      }

      // Shift: axes to reflect the updated order:
      for (unsigned char j=next_axis_index; j<ten.dimension()-1; ++j)
	current_axis_order[j] = current_axis_order[j+1];
      current_axis_order[ten.dimension()-1] = next_axis;
    }

    // Data was last transposed into buffer_to, which was swapped with
    // buffer_from. If buffer_from is not the true source array, then
    // move it into ten. 
    if (buffer_from != &ten[0ul])
      ten = std::move(buffer);

    // Change the shape to correspond:
    Vector<unsigned long> old_shape = ten.data_shape();

    Vector<unsigned long> new_shape(ten.dimension());
    for (unsigned char i=0; i<ten.dimension(); ++i)
      new_shape[i] = old_shape[ new_axis_order[i] ];

    ten.reshape(new_shape);
  }
}

template <typename T>
Tensor<T> cache_friendly_transposed(const Tensor<T> & ten, const Vector<unsigned char> & new_axis_order) {
  Tensor<T> res = ten;
  transpose(res, new_axis_order);
  return res;
}

template <typename T>
void transpose(Tensor<T> & ten, const Vector<unsigned char> & new_axis_order) {
  if (ten.flat_size() < SIZE_WHERE_NAIVE_TRANSPOSE_BECOMES_SLOWER)
    naive_transpose(ten, new_axis_order);
  else
    cache_friendly_transpose(ten, new_axis_order);
}

template <typename T>
Tensor<T> transposed(const Tensor<T> & ten, const Vector<unsigned char> & new_axis_order) {
  if (ten.flat_size() < SIZE_WHERE_NAIVE_TRANSPOSE_BECOMES_SLOWER)
    return naive_transposed(ten, new_axis_order);
  return cache_friendly_transposed(ten, new_axis_order);
}

#endif
