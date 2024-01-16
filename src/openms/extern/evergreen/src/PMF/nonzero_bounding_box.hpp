#ifndef _NONZERO_BOUNDING_BOX_HPP
#define _NONZERO_BOUNDING_BOX_HPP

inline std::array<Vector<unsigned long>, 2> nonzero_bounding_box(const Tensor<double> & rhs, const double relative_mass_threshold_for_bounding_box) {
  // Initialize min with value greater than maximum possible, and
  // initialize max with minimum possible:
  Vector<unsigned long> min_tup  = rhs.data_shape(), max_tup(rhs.dimension());

  double max_mass = max(rhs.flat());
  const double epsilon = max_mass*relative_mass_threshold_for_bounding_box;

  bool exist_any_nonzero = false;
  enumerate_for_each_tensors([&min_tup, &max_tup, &exist_any_nonzero, epsilon](const_tup_t counter, const unsigned char dim, double val){
      if (val > epsilon) {
	exist_any_nonzero = true;
	for (unsigned char i=0; i<dim; ++i) {
	  min_tup[i] = std::min(min_tup[i], counter[i]);
	  max_tup[i] = std::max(max_tup[i], counter[i]);
	}
      }
    },
    rhs.data_shape(),
    rhs);

  assert(exist_any_nonzero && "PMF must be constructed from a tensor with at least one nonzero entry; this model has a contradiction in it (or is numerically very close to a contradiction).");
  return {{min_tup, max_tup}};
}

#endif
