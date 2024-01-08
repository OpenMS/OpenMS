#ifndef _VECTORCOMPARISON_HPP
#define _VECTORCOMPARISON_HPP

// Note: Vector comparison is currently performed overall rather than
// element-wise. This could be changed by having every function return
// a Vector<bool>, which could be used with any(...) and all(...)
// functions, which would aggregate. The downside of that approach is
// that [1,2] == [1,4,5] would no longer return false; instead, it
// would assert(false).

// However, a downside to the current implementation is that x <= y
// can be false and x > y can be false (typically, exactly one of
// these must be true).

template <typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
bool operator ==(const VectorLike<T, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
  if (lhs.size() != rhs.size())
    return false;
  for (unsigned long k=0; k<lhs.size(); ++k)
    if (lhs[k] != rhs[k])
      return false;
  return true;
}

template <typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
bool operator !=(const VectorLike<T, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
  return !( lhs == rhs );
}

template <typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
bool operator <(const VectorLike<T, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
  if (lhs.size() != rhs.size())
    return false;
  for (unsigned long k=0; k<lhs.size(); ++k)
    if (lhs[k] >= rhs[k])
      return false;
  return true;
}

template <typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
bool operator >(const VectorLike<T, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
  if (lhs.size() != rhs.size())
    return false;
  for (unsigned long k=0; k<lhs.size(); ++k)
    if (lhs[k] <= rhs[k])
      return false;
  return true;
}

template <typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
bool operator <=(const VectorLike<T, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
  if (lhs.size() != rhs.size())
    return false;
  for (unsigned long k=0; k<lhs.size(); ++k)
    if (lhs[k] > rhs[k])
      return false;
  return true;
}

template <typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
bool operator >=(const VectorLike<T, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
  if (lhs.size() != rhs.size())
    return false;
  for (unsigned long k=0; k<lhs.size(); ++k)
    if (lhs[k] < rhs[k])
      return false;
  return true;
}

template <typename T, template <typename> class VECTOR_A>
bool operator ==(const VectorLike<T, VECTOR_A> & lhs, T rhs) {
  for (unsigned long k=0; k<lhs.size(); ++k)
    if (lhs[k] != rhs)
      return false;
  return true;
}

template <typename T, template <typename> class VECTOR_A>
bool operator !=(const VectorLike<T, VECTOR_A> & lhs, T rhs) {
  return ! (lhs == rhs);
}

template <typename T, template <typename> class VECTOR_A>
bool operator <(const VectorLike<T, VECTOR_A> & lhs, T rhs) {
  for (unsigned long k=0; k<lhs.size(); ++k)
    if (lhs[k] >= rhs)
      return false;
  return true;
}

template <typename T, template <typename> class VECTOR_A>
bool operator >(const VectorLike<T, VECTOR_A> & lhs, T rhs) {
  for (unsigned long k=0; k<lhs.size(); ++k)
    if (lhs[k] <= rhs)
      return false;
  return true;
}

template <typename T, template <typename> class VECTOR_A>
bool operator <=(const VectorLike<T, VECTOR_A> & lhs, T rhs) {
  for (unsigned long k=0; k<lhs.size(); ++k)
    if (lhs[k] > rhs)
      return false;
  return true;
}

template <typename T, template <typename> class VECTOR_A>
bool operator >=(const VectorLike<T, VECTOR_A> & lhs, T rhs) {
  for (unsigned long k=0; k<lhs.size(); ++k)
    if (lhs[k] < rhs)
      return false;
  return true;
}

#endif
