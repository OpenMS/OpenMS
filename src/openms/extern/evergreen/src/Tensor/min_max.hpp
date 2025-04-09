#ifndef _MIN_MAX_HPP
#define _MIN_MAX_HPP

template <typename T, template <typename> class VECTOR>
T min(const VectorLike<T, VECTOR> & rhs) {
  #ifdef SHAPE_CHECK
  assert(rhs.size() > 0);
  #endif
  
  T res = rhs[0];
  for (unsigned long k=1; k<rhs.size(); ++k)
    if (res > rhs[k])
      res = rhs[k];
  return res;
}

template <typename T, template <typename> class VECTOR>
T max(const VectorLike<T, VECTOR> & rhs) {
  #ifdef SHAPE_CHECK
  assert(rhs.size() > 0);
  #endif
  
  T res = rhs[0];
  for (unsigned long k=1; k<rhs.size(); ++k)
    if (res < rhs[k])
      res = rhs[k];
  return res;
}

#endif
