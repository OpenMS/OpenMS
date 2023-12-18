#ifndef _PRODUCT_HPP
#define _PRODUCT_HPP

template <typename T>
inline unsigned long product(const T* __restrict const v, unsigned long length) {
  T res = 1ul;
  for (unsigned long k=0; k<length; ++k)
    res *= v[k];
  return res;
}

template <typename T, template <typename> class VECTOR>
inline T product(const VectorLike<T, VECTOR> & v) {
  return product(static_cast<const T*const>(v), v.size());
}

#endif

