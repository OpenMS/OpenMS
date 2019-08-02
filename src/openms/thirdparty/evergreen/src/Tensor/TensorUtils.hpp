#ifndef _TENSOR_UTILS_HPP
#define _TENSOR_UTILS_HPP

#include <vector>

#include "product.hpp"
#include "sum.hpp"
#include "Vector.hpp"

#ifndef MAX_TENSOR_DIMENSION
#define MAX_TENSOR_DIMENSION 24
#endif

// Tuple types:
typedef unsigned long* __restrict const tup_t;
typedef const unsigned long* __restrict const const_tup_t;

inline void print_tuple(const_tup_t tup, unsigned char dim) {
  for (unsigned char i=0; i<dim; ++i)
    std::cout << tup[i] << " ";
  std::cout << std::endl;
}

inline unsigned long flat_length(const_tup_t shape, unsigned char dimension) {
  if (dimension > 0)
    return product(shape, dimension);
  return 0;
}

inline unsigned long flat_length(const Vector<unsigned long> & shape) {
  return flat_length(static_cast<const unsigned long*const>(shape), shape.size());
}

inline void advance_tuple(unsigned long* __restrict tup, const_tup_t shape, unsigned char dimension) {
  ++tup[dimension-1];
    
  for (unsigned char k=dimension-1; k>=1; --k)
    if ( tup[k] >= shape[k] ) {
      ++tup[k-1];
      tup[k] = 0;
    }
    else
      // No more carry operations:
      return;
}

inline unsigned long tuple_to_index(const_tup_t tup, const_tup_t shape, unsigned char dimension) {
  unsigned long res = 0;
  unsigned char k;

  for (k=1; k<dimension; ++k) {
    res += tup[k-1];
    res *= shape[k];
  }
  res += tup[k-1];
  return res;
}

template <unsigned int DIMENSION>
inline unsigned long tuple_to_index_fixed_dimension(const_tup_t tup, const_tup_t shape) {
    unsigned long res = 0;
    unsigned int k;
    for (k=0; k<DIMENSION-1; ++k) {
        res += tup[k];
        res *= shape[k+1]; }
    res += tup[k];
    return res;
}

// Note: This is not very efficient, but is useful for debugging.
inline unsigned long* index_to_tuple(unsigned long index, const unsigned long* __restrict const shape, unsigned int dimension) {
  unsigned long* __restrict result = aligned_calloc<unsigned long>(dimension);
  for (int i=dimension-1; index>0 && i>=0; --i) {
    unsigned long next_axis = shape[i];
    
    // Note: There may be a speedup lurking in here where shared work
    // between index / next_axis and index % next_axis can be reused;
    // however, this code will only be used for bounds checking, so
    // speed is not very important.
    unsigned long next_value = index % next_axis;
    result[i] = next_value;
    
    index /= next_axis;
  }
  return result;
}


// No assertions when there are no duplicate indices and all are in
// range. Could also be implemented a set in O(n log(n)), but this
// version is in O(n).
inline void verify_subpermutation(const Vector<unsigned char> & permutation, unsigned char dim) {
  std::vector<bool> indices(dim, false);
  for (unsigned char i=0; i<permutation.size(); ++i) {
    // All values must be in 0, 1, ... n-1, where n is the number of
    // dimensions allowed:
    assert(permutation[i] < dim);

    indices[ permutation[i] ] = true;
  }
  unsigned char cardinality = 0;
  for (unsigned char i=0; i<permutation.size(); ++i)
    cardinality += indices[ permutation[i] ];

  // All indices must be included exactly once (therefore, there must
  // be no duplicates). Given all indices must also be in range (by
  // the assertion above), this means it is a valid subpermutation.
  assert(cardinality == permutation.size());
}

inline void verify_permutation(const Vector<unsigned char> & permutation) {
  verify_subpermutation(permutation, permutation.size());
}

#endif
