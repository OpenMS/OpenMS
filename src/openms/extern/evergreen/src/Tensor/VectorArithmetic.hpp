#ifndef _VECTORARITHMETIC_HPP
#define _VECTORARITHMETIC_HPP

#include "VectorTRIOT.hpp"

template <typename T>
class Vector;

template <typename T>
class VectorView;

// +=, -=, ... with & lhs:
template <typename S, typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
const WritableVectorLike<S, VECTOR_A> & operator +=(WritableVectorLike<S, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
#ifdef SHAPE_CHECK
  assert(lhs.size() == rhs.size());
#endif

  apply_vectors([](S & vL, T vR){ vL += vR; }, lhs.size(), lhs, rhs);
  return lhs;
}

template <typename S, typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
const WritableVectorLike<S, VECTOR_A> & operator -=(WritableVectorLike<S, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
#ifdef SHAPE_CHECK
  assert(lhs.size() == rhs.size());
#endif

  apply_vectors([](S & vL, T vR){ vL -= vR; }, lhs.size(), lhs, rhs);
  return lhs;
}

template <typename S, typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
const WritableVectorLike<S, VECTOR_A> & operator *=(WritableVectorLike<S, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
#ifdef SHAPE_CHECK
  assert(lhs.size() == rhs.size());
#endif

  apply_vectors([](S & vL, T vR){ vL *= vR; }, lhs.size(), lhs, rhs);
  return lhs;
}

template <typename S, typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
const WritableVectorLike<S, VECTOR_A> & operator /=(WritableVectorLike<S, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
#ifdef SHAPE_CHECK
  assert(lhs.size() == rhs.size());
#endif

  apply_vectors([](S & vL, T vR){ vL /= vR; }, lhs.size(), lhs, rhs);
  return lhs;
}

// +=, -=, ... with & lhs, value rhs:
template <typename T, template <typename> class VECTOR_A>
const WritableVectorLike<T, VECTOR_A> & operator +=(WritableVectorLike<T, VECTOR_A> & lhs, T rhs) {
  apply_vectors([&rhs](T & vL){ vL += rhs; }, lhs.size(), lhs);
  return lhs;
}

template <typename T, template <typename> class VECTOR_A>
const WritableVectorLike<T, VECTOR_A> & operator -=(WritableVectorLike<T, VECTOR_A> & lhs, T rhs) {
  apply_vectors([&rhs](T & vL){ vL -= rhs; }, lhs.size(), lhs);
  return lhs;
}

template <typename T, template <typename> class VECTOR_A>
const WritableVectorLike<T, VECTOR_A> & operator *=(WritableVectorLike<T, VECTOR_A> & lhs, T rhs) {
  apply_vectors([&rhs](T & vL){ vL *= rhs; }, lhs.size(), lhs);
  return lhs;
}

template <typename T, template <typename> class VECTOR_A>
const WritableVectorLike<T, VECTOR_A> & operator /=(WritableVectorLike<T, VECTOR_A> & lhs, T rhs) {
  apply_vectors([&rhs](T & vL){ vL /= rhs; }, lhs.size(), lhs);
  return lhs;
}

// +=, -=, ... with && lhs:
template <typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
const WritableVectorLike<T, VECTOR_A> & operator +=(WritableVectorLike<T, VECTOR_A> && lhs, const VectorLike<T, VECTOR_B> & rhs) {
  return lhs += rhs;
}

template <typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
const WritableVectorLike<T, VECTOR_A> & operator -=(WritableVectorLike<T, VECTOR_A> && lhs, const VectorLike<T, VECTOR_B> & rhs) {
  return lhs -= rhs;
}

template <typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
const WritableVectorLike<T, VECTOR_A> & operator *=(WritableVectorLike<T, VECTOR_A> && lhs, const VectorLike<T, VECTOR_B> & rhs) {
  return lhs *= rhs;
}

template <typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
const WritableVectorLike<T, VECTOR_A> & operator /=(WritableVectorLike<T, VECTOR_A> && lhs, const VectorLike<T, VECTOR_B> & rhs) {
  return lhs /= rhs;
}

// +=, -=, ... with & lhs, value rhs:
template <typename T, template <typename> class VECTOR_A>
const WritableVectorLike<T, VECTOR_A> & operator +=(WritableVectorLike<T, VECTOR_A> && lhs, T rhs) {
  return lhs += rhs;
}

template <typename T, template <typename> class VECTOR_A>
const WritableVectorLike<T, VECTOR_A> & operator -=(WritableVectorLike<T, VECTOR_A> && lhs, T rhs) {
  return lhs -= rhs;
}

template <typename T, template <typename> class VECTOR_A>
const WritableVectorLike<T, VECTOR_A> & operator *=(WritableVectorLike<T, VECTOR_A> && lhs, T rhs) {
  return lhs *= rhs;
}

template <typename T, template <typename> class VECTOR_A>
const WritableVectorLike<T, VECTOR_A> & operator /=(WritableVectorLike<T, VECTOR_A> && lhs, T rhs) {
  return lhs /= rhs;
}

// +, -, ... 
template <typename S, typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
Vector<S> operator +(const VectorLike<S, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
  Vector<S> result = lhs;
  return result += rhs;
}

template <typename S, typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
Vector<S> operator -(const VectorLike<S, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
  Vector<S> result = lhs;
  return result -= rhs;
}

template <typename S, typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
Vector<S> operator *(const VectorLike<S, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
  Vector<S> result = lhs;
  return result *= rhs;
}

template <typename S, typename T, template <typename> class VECTOR_A, template <typename> class VECTOR_B>
Vector<S> operator /(const VectorLike<S, VECTOR_A> & lhs, const VectorLike<T, VECTOR_B> & rhs) {
  Vector<S> result = lhs;
  return result /= rhs;
}

// +, -, ... Vector with value:
template <typename T, template <typename> class VECTOR_A>
Vector<T> operator +(const VectorLike<T, VECTOR_A> & lhs, T rhs) {
  Vector<T> result = lhs;
  return result += rhs;
}

template <typename T, template <typename> class VECTOR_A>
Vector<T> operator -(const VectorLike<T, VECTOR_A> & lhs, T rhs) {
  Vector<T> result = lhs;
  return result -= rhs;
}

template <typename T, template <typename> class VECTOR_A>
Vector<T> operator *(const VectorLike<T, VECTOR_A> & lhs, T rhs) {
  Vector<T> result = lhs;
  return result *= rhs;
}

template <typename T, template <typename> class VECTOR_A>
Vector<T> operator /(const VectorLike<T, VECTOR_A> & lhs, T rhs) {
  Vector<T> result = lhs;
  return result /= rhs;
}

template <typename T, template <typename> class VECTOR_A>
Vector<T> operator /(T lhs, const VectorLike<T, VECTOR_A> & rhs) {
  Vector<T> result(rhs.size(), lhs);
  return result /= rhs;
}

template <typename T>
Vector<T> seq(long length) {
  Vector<T> result(length);
  for (unsigned long k=0; k<result.size(); ++k)
    result[k] = T(k);
  return result;
}

#endif
