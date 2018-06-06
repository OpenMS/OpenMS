#ifndef _SCALED_PMF_HPP
#define _SCALED_PMF_HPP

#include "squared.hpp"

// Note: could be sped up to not make local Vector objects and to not
// use tuple indexing (using integer offset):
inline void add_scaled_outcome(Tensor<double> & ten, const Vector<long> & new_first_support, const Vector<double> & scaled_tup, double mass) {
  // For performance, don't bother if the mass is 0.
  if (mass > 0.0) {
    Vector<unsigned long> start_index(ten.dimension());
    for (unsigned char i=0; i<ten.dimension(); ++i)
      start_index[i] = floor(scaled_tup[i]) - new_first_support[i];

    Vector<unsigned long> scaled_bounding_box(ten.dimension());
    for (unsigned char i=0; i<ten.dimension(); ++i)
      scaled_bounding_box[i] = (ceil(scaled_tup[i]) - new_first_support[i] + 1) - start_index[i];

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

inline PMF scaled_pmf(const PMF & pmf, const Vector<double> & factor) {
  Vector<double> extreme_a = pmf.first_support();
  extreme_a *= factor;
  Vector<double> extreme_b = pmf.last_support();
  extreme_b *= factor;

  Vector<long> new_first_support(pmf.dimension());
  Vector<long> new_last_support(pmf.dimension());

  for (unsigned char i=0; i<pmf.dimension(); ++i) {
    new_first_support[i] = floor( std::min(extreme_a[i], extreme_b[i]) );
    new_last_support[i] = ceil( std::max(extreme_a[i], extreme_b[i]) );
  }

  Tensor<double> result_table(new_last_support - new_first_support + 1L);
  Vector<double> scaled_tup(pmf.dimension());

  enumerate_for_each_tensors([&pmf, &result_table, &new_first_support, &scaled_tup, &factor](const_tup_t tup, const unsigned char dim, double mass){
      for (unsigned char i=0; i<dim; ++i)
	scaled_tup[i] = (pmf.first_support()[i] + long(tup[i])) * factor[i];
      add_scaled_outcome(result_table, new_first_support, scaled_tup, mass);
    },
    pmf.table().data_shape(),
    pmf.table());

  return PMF(new_first_support, std::move(result_table));
}

#endif
