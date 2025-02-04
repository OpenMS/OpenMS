#ifndef _ANY_ALL_HPP
#define _ANY_ALL_HPP

template <template <typename> class VECTOR>
bool any(const VectorLike<bool, VECTOR> & rhs) {
  for (unsigned long k=0; k<rhs.size(); ++k)
    if (rhs[k])
      return true;
  return false;
}

template <template <typename> class VECTOR>
bool all(const VectorLike<bool, VECTOR> & rhs) {
  for (unsigned long k=0; k<rhs.size(); ++k)
    if ( ! rhs[k])
      return false;
  return true;
}

#endif
