#ifndef _SCALED_PMF_DITHER_HPP
#define _SCALED_PMF_DITHER_HPP

#include "squared.hpp"

// Note: For performance, it may be beneficial to manually add the
// offset and use the Tensor[unsigned long] operator instead;
// constructing a tensor view for each function call will construct
// a Vector, which may not get optimized out (unclear).
inline static void add_scaled_outcome_dither(Tensor<double> & ten, const Vector<double> & weighting_partition, const Vector<unsigned long> & scaled_counter_lower, const Vector<unsigned long> & scaled_bounding_box, double mass) {
  if (mass > 0.0) {
    enumerate_apply_tensors([mass, &weighting_partition](const_tup_t tup, const unsigned char dim, double & ten_value){
	double mass_partition = 1.0;
	for (unsigned char i=0; i<dim; ++i)
	  // When tup[i] == 0, use weighting_partition[i]
	  // When tup[i] == 1, use 1-weighting_partition[i]
	  mass_partition *= tup[i]*(1.0-weighting_partition[i]) + (1-tup[i])*weighting_partition[i];

	ten_value += mass_partition*mass;
      },
      scaled_bounding_box,
      ten.start_at(scaled_counter_lower));
  }
}

inline PMF scaled_pmf_dither(const PMF & pmf, const Vector<double> & factor, double sigma_squared) {
  // Largest index is shape - 1:
  Vector<double> abs_factor = factor;
  for (unsigned char i=0; i<factor.size(); ++i)
    abs_factor[i] = fabs(abs_factor[i]);

  Vector<long> res_shape = pmf.table().view_shape();
  res_shape -= 1L;
  for (unsigned char i=0; i<res_shape.size(); ++i)
    res_shape[i] = ceil( res_shape[i]*abs_factor[i] );

  // For the result shape, add +2 (the first support could round
  // down while last support could round up):
  Tensor<double> res_table(std::move(res_shape+2L));

  const Vector<long> & first_sup = pmf.first_support();
  Vector<long> last_sup = pmf.last_support();
  Vector<double> new_first_sup_double(pmf.dimension());
  for (unsigned char i=0; i<pmf.dimension(); ++i)
    new_first_sup_double[i] = std::min(first_sup[i]*factor[i], last_sup[i]*factor[i]);

  Vector<long> new_first_sup(pmf.dimension());
  for (unsigned char i=0; i<pmf.dimension(); ++i)
    new_first_sup[i] = floor(new_first_sup_double[i]);

  Vector<double> scaled_outcome(pmf.dimension());
  Vector<unsigned long> scaled_counter_lower(pmf.dimension());
  Vector<unsigned long> scaled_bounding_box(pmf.dimension());
  enumerate_for_each_tensors([&res_table, &scaled_counter_lower, &scaled_bounding_box, &factor, &first_sup, &new_first_sup, &scaled_outcome, sigma_squared](const_tup_t index, const unsigned char dim, double mass){
      for (unsigned char i=0; i<dim; ++i)
	scaled_outcome[i] = (long(index[i]) + first_sup[i])*factor[i];

      for (unsigned char i=0; i<dim; ++i)
	scaled_counter_lower[i] = floor(scaled_outcome[i]) - new_first_sup[i];

      for (unsigned char i=0; i<dim; ++i)
	scaled_bounding_box[i] = ceil(scaled_outcome[i]) - floor(scaled_outcome[i]) + 1;

      for (unsigned char i=0; i<dim; ++i) {
	if (scaled_bounding_box[i] == 1)
	  scaled_outcome[i] = 1.0;
	else {
	  // scaled_bounding_box[i] == 2:
	  scaled_outcome[i] -= floor(scaled_outcome[i]);
	  // Outcome is either +0 or +1. These are then smoothed to
	  // weight how much of the mass is partitioned into the +0
	  // outcome.
	  double smoothed_0 = exp(-squared(scaled_outcome[i])/sigma_squared);
	  double smoothed_1 = exp(-squared(scaled_outcome[i]-1.0)/sigma_squared);
	  scaled_outcome[i] = smoothed_0 / (smoothed_0 + smoothed_1);
	}
	  
	// scaled_outcome[i] is now the partition in the +0 category.
      }

      add_scaled_outcome_dither(res_table, scaled_outcome, scaled_counter_lower, scaled_bounding_box, mass);
    },
    pmf.table().view_shape(),
    pmf.table());

  auto result = PMF(new_first_sup, std::move(res_table));

  return result;
}

#endif
