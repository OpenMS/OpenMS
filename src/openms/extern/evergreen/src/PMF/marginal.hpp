#ifndef _MARGINAL_HPP
#define _MARGINAL_HPP

// Note: it may be possible to optimize these marginal routines for
// rvalue references, directly writing the results to the first
// element in each collapsed group, then collapsing values down, and
// then shrinking the tensor:


// Empirically chosen:
const unsigned long SIZE_WHERE_NAIVE_MARGINAL_BECOMES_SLOWER = 32;

// Naive marginal: more optimized than the obvious version; saves
// time by iterating over new tuple and then inside that iterating
// over remaining tuple. This allows the pow computation to be
// performed in a numerically stable manner (by dividing out the
// maximum) without computing a separate tensor full of the maximum
// marginals first.
inline Tensor<double> naive_marginal(const Tensor<double> & table, Vector<unsigned char> axes_to_keep, double p) {
  #ifdef SHAPE_CHECK
  verify_subpermutation(axes_to_keep, table.dimension());
  #endif

  unsigned char k;
  
  Vector<unsigned long> new_shape(axes_to_keep.size());
  for (k=0; k<axes_to_keep.size(); ++k)
    new_shape[k] = table.data_shape()[ axes_to_keep[k] ];
  
  std::vector<bool> axes_eliminated(table.dimension(), true);
  for (unsigned char i=0; i<axes_to_keep.size(); ++i)
    axes_eliminated[ axes_to_keep[i] ] = false;

  Vector<unsigned char> axes_to_remove( table.dimension() - axes_to_keep.size() );
  for (unsigned char i=0, j=0; i<axes_eliminated.size(); ++i)
    if (axes_eliminated[i]) {
      axes_to_remove[j] = i;
      ++j;
    }
  
  Vector<unsigned long> shape_removed( axes_to_remove.size() );
  for (unsigned char i=0; i<shape_removed.size(); ++i)
    shape_removed[i] = table.data_shape()[ axes_to_remove[i] ];
  
  Tensor<double> new_table(new_shape);
  Vector<unsigned long> full_counter(table.dimension());
  enumerate_apply_tensors([&axes_to_keep, &axes_to_remove, &full_counter, &table, p, &shape_removed](const_tup_t counter_kept, const unsigned char dim_kept, double & new_val){
      for (unsigned char i=0; i<dim_kept; ++i)
	full_counter[ axes_to_keep[i] ] = counter_kept[i];
      
      double max_val = 0.0;
      enumerate_for_each_tensors([&axes_to_remove, &full_counter, &table, p, &max_val, dim_kept](const_tup_t counter_removed, const unsigned char dim_removed){
	  for (unsigned char i=0; i<dim_removed; ++i)
	    full_counter[ axes_to_remove[i] ] = counter_removed[i];
	  
	  unsigned long full_index = tuple_to_index(full_counter, table.data_shape(), dim_kept + dim_removed);
	  
	  max_val = std::max(max_val, table[full_index]);
	},
	shape_removed);

      if ( max_val > tau_denom ) {
	enumerate_for_each_tensors([&axes_to_remove, &full_counter, &table, p, max_val, dim_kept, &new_val](const_tup_t counter_removed, const unsigned char dim_removed){
	    for (unsigned char i=0; i<dim_removed; ++i)
	      full_counter[ axes_to_remove[i] ] = counter_removed[i];
	  
	    unsigned long full_index = tuple_to_index(full_counter, table.data_shape(), dim_kept + dim_removed);
	    
	    new_val += custom_pow(table[full_index] / max_val, p);
	  },
	  shape_removed);
      }
      // Otherwise, let result = 0.0. 

      // Note: The numeric stability of this could possibly be
      // improved; e.g., when max is close to zero, the 1-norm may not
      // necessarily be 0.
      
      new_val = custom_pow(new_val, 1.0/p) * max_val;
    },
    new_table.data_shape(),
    new_table);

  return new_table;
}

// First transpose so that the innermost indices are the ones lost:
// (this makes marginalization more cache friendly):
inline Tensor<double> transposed_marginal(const Tensor<double> & table, Vector<unsigned char> axes_to_keep, double p) {
#ifdef SHAPE_CHECK
  verify_subpermutation(axes_to_keep, table.dimension());
#endif

  // Initialize variables:
  unsigned long k;

  Vector<unsigned long> new_shape(axes_to_keep.size());
  for (k=0; k<axes_to_keep.size(); ++k)
    new_shape[k] = table.data_shape()[ axes_to_keep[k] ];

  // Transpose so that the axes kept are first:
  Vector<unsigned char> new_axis_order(table.dimension());
  copy( new_axis_order, axes_to_keep );

  std::vector<bool> axes_eliminated(table.dimension(), true);
  for (unsigned char i=0; i<axes_to_keep.size(); ++i)
    axes_eliminated[ axes_to_keep[i] ] = false;
  for (unsigned char i=0, j=0; i<axes_eliminated.size(); ++i)
    if (axes_eliminated[i]) {
      new_axis_order[j+axes_to_keep.size()] = i;
      ++j;
    }

  Tensor<double> table_copy = table;
  transpose(table_copy, new_axis_order);

  if (axes_to_keep.size() == table.dimension())
    // all axes are kept:
    return table_copy;

  Tensor<double> new_table(new_shape);
  // Compute marginal:

  // Note: this strategy can be more efficient (unrolling final axes
  // into single, longer axis), but is not valid for TensorView since
  // the memory would not necessarily be contiguous. If this code were
  // generalized for TensorView, the following would therefore need to
  // be changed.
  unsigned long removed_axes_flat_length = flat_length( table_copy.data_shape().start_at_const(axes_to_keep.size() ) );
  enumerate_apply_tensors([&table_copy, &removed_axes_flat_length, p](const_tup_t counter, const unsigned char dim, double & new_val){
      unsigned long bias = tuple_to_index(counter, table_copy.data_shape(), dim) * removed_axes_flat_length;
      double max_val = 0.0;
      for (unsigned long k=0; k<removed_axes_flat_length; ++k)
	max_val = std::max(max_val, table_copy[bias + k]);
      if ( max_val > tau_denom ) {
	for (unsigned long k=0; k<removed_axes_flat_length; ++k)
	  new_val += custom_pow(table_copy[bias + k]/max_val, p);
	new_val = custom_pow(new_val, 1.0/p) * max_val;
      }
      // otherwise the result will be 0.0:
    },
    new_table.data_shape(),
    new_table);

  return new_table;
}

inline Tensor<double> marginal(const Tensor<double> & table, Vector<unsigned char> axes_to_keep, double p) {
  if (table.flat_size() < SIZE_WHERE_NAIVE_MARGINAL_BECOMES_SLOWER)
    return naive_marginal(table, axes_to_keep, p);
  return transposed_marginal(table, axes_to_keep, p);
}

#endif
