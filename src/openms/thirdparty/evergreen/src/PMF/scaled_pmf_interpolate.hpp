#ifndef _SCALED_PMF_INTERPOLATE_HPP
#define _SCALED_PMF_INTERPOLATE_HPP

#include "squared.hpp"

// Note: could be sped up to not make local Vector objects and to not
// use tuple indexing (using integer offset):
inline void add_scaled_outcome_interpolate(Tensor<double> & ten, const Vector<long> & new_first_support, const Vector<double> & scaled_tup, const Vector<double> & next_scaled_tup, double mass, const Vector<double> & factor) {
  // For performance, don't bother if the mass is 0.
  if (mass > 0.0) {
    Vector<unsigned long> start_index(ten.dimension());
    for (unsigned char i=0; i<ten.dimension(); ++i)
      start_index[i] = floor( std::min(scaled_tup[i], next_scaled_tup[i]) ) - new_first_support[i];
    
    Vector<unsigned long> scaled_bounding_box(ten.dimension());
    for (unsigned char i=0; i<ten.dimension(); ++i)
      scaled_bounding_box[i] = ceil( std::max(scaled_tup[i], next_scaled_tup[i]) - start_index[i] ) - new_first_support[i];

    // Split the mass over the partitioned boxes:
    for (unsigned char i=0; i<ten.dimension(); ++i)
      mass /= scaled_bounding_box[i];

    enumerate_apply_tensors([mass](const_tup_t tup, const unsigned char dim, double & val){
	val += mass;
      },
      scaled_bounding_box,
      ten.start_at(start_index));
  }
}

inline PMF scaled_pmf_interpolate(const PMF & pmf, const Vector<double> & factor) {
  Vector<double> extreme_a = pmf.first_support();
  extreme_a *= factor;
  Vector<double> extreme_b = pmf.last_support();
  extreme_b *= factor;

  Vector<long> new_first_support(pmf.dimension());
  Vector<unsigned long> new_shape(pmf.dimension());

  for (unsigned char i=0; i<pmf.dimension(); ++i) {
    new_first_support[i] = floor( std::min(extreme_a[i], extreme_b[i]) );
    new_shape[i] = long(ceil( std::max(extreme_a[i], extreme_b[i]) )) - new_first_support[i] + long(ceil(fabs(factor[i]))) ;
  }

  Tensor<double> result_table(new_shape);

  Vector<double> scaled_tup(pmf.dimension());
  Vector<double> next_scaled_tup(pmf.dimension());

  enumerate_for_each_tensors([&pmf, &result_table, &new_first_support, &scaled_tup, &next_scaled_tup, &factor](const_tup_t tup, const unsigned char dim, double mass){
      for (unsigned char i=0; i<dim; ++i) {
	scaled_tup[i] = (pmf.first_support()[i] + long(tup[i])) * factor[i];
	next_scaled_tup[i] = scaled_tup[i] + factor[i];
	
	// This hack is necessary in order to allow scaling by S and
	// then by 1/S come out the same as scaling by -S and -1/S. It
	// occurs because the continuous interpretation of bin 1 is
	// actually [1,2); however, this becomes inverted with
	// negative support: -1 indicates [-1, 0), which is not
	// symmetric. By shifting negatives left by 1 (meaning that
	// their scaled interpretations will shift left by factor[i]),
	// -1 will indicate (-2,-1].
	if (factor[i] < 0) {
	  scaled_tup[i] -= factor[i];
	  next_scaled_tup[i] -= factor[i];
	}
      }
      add_scaled_outcome_interpolate(result_table, new_first_support, scaled_tup, next_scaled_tup, mass, factor);
    },
    pmf.table().data_shape(),
    pmf.table());

  return PMF(new_first_support, std::move(result_table));
}

#endif
