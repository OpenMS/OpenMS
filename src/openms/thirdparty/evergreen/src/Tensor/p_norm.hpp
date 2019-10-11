#ifndef _P_NORM_HPP
#define _P_NORM_HPP

template <typename T, template <typename> class VECTOR>
T p_norm(const VectorLike<T, VECTOR> & rhs, T p) {
  #ifdef SHAPE_CHECK
  assert(rhs.size() > 0);
  #endif

  T max_val = max(rhs);

  T res = pow((rhs[0]/max_val), p);

  for (unsigned long k=1; k<rhs.size(); ++k)
    res += pow(rhs[k]/max_val, p);
  return max_val*pow(res, T(1.0)/p);
}

#endif
