#ifndef _SUM_HPP
#define _SUM_HPP

template <typename T>
inline T sum(const T* __restrict const v, unsigned long length) {
  T res = 0;
  for (unsigned long k=0; k<length; ++k)
    res += v[k];
  return res;
}

template <typename T, template <typename> class VECTOR>
inline T sum(const VectorLike<T, VECTOR> & v) {
  return sum(static_cast<const T*const>(v), v.size());
}

#endif

