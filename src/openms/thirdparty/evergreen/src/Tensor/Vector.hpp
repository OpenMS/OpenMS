#ifndef _VECTOR_HPP
#define _VECTOR_HPP

#include <iostream>
#include <vector>

#include "alloc.hpp"
#include "VectorLike.hpp"
#include "VectorView.hpp"
#include "VectorArithmetic.hpp"
#include "VectorComparison.hpp"
#include "min_max.hpp"
#include "p_norm.hpp"
#include "any_all.hpp"


// Note: Vector<T> is for simple numeric T types; uses aligned_malloc
// rather than new[], so no constructor is called:
template <typename T>
class Vector : public WritableVectorLike<T, Vector> {
protected:
  unsigned long _length;
  T* __restrict _data;
public:
  template <typename S>
  friend class Vector;
  
  Vector():
    _length(0),
    _data(NULL)
  { }
  explicit Vector(unsigned long length_param):
    _length(length_param)
  {
    _data = aligned_calloc<T>(_length);
  }
  explicit Vector(unsigned long length_param, T fill_value):
    _length(length_param)
  {
    _data = aligned_malloc<T>(_length);
    this->fill(fill_value);
  }
  explicit Vector(unsigned long length_param, const T*const fill_vec):
    _length(length_param)
  {
    _data = aligned_malloc<T>(_length);
    for (unsigned long k=0; k<_length; ++k)
      _data[k] = fill_vec[k];
  }

  // Note: the following two constructors would be duplicates; this is
  // necessary because copy ctor is deleted when defining move
  // ctor. Delegate the constructor to the most general type to avoid
  // duplicate code:

  // Casts element S to element T:
  template <typename S, template <typename> class VECTOR>
  Vector(const VectorLike<S, VECTOR> & rhs):
    _length(rhs.size())
  {
    _data = aligned_malloc<T>(_length);
    for (unsigned long k=0; k<_length; ++k)
      _data[k] = (T)rhs[k];
  }
  Vector(const Vector<T> & rhs): 
    Vector( static_cast<const VectorLike<T, evergreen::Vector> &>(rhs) )
  { }

  Vector(Vector<T> && rhs):
    _length(rhs._length)
  {
    _data = rhs._data;
    rhs._data = NULL;
    rhs._length = 0;
  }
  Vector(const std::initializer_list<T> & rhs):
    _length(rhs.size())
  {
    _data = aligned_malloc<T>(_length);
    unsigned long k=0;
    for (auto iter = rhs.begin(); iter != rhs.end(); ++k, ++iter)
      _data[k] = *iter;
  }

  Vector(const std::vector<T> & rhs):
      _length(rhs.size())
  {
    _data = aligned_malloc<T>(_length);
    for (unsigned long k=0; k<_length; ++k)
      _data[k] = rhs[k];
  }


  // Casts element S to element T:
  template <typename S, template <typename> class VECTOR>
  const Vector & operator =(const VectorLike<S, VECTOR> & rhs) {
    // Ensure the vectors do not refer to the same memory (otherwise
    // freeing _data would also free rhs._data):
    bool no_overlap = (_data + _length <= &rhs[0ul]) || (&rhs[0ul] + rhs.size() <= _data);
    assert( no_overlap );
    clear();

    _length = rhs.size();
    _data = aligned_malloc<T>(_length);
    for (unsigned long k=0; k<_length; ++k)
      _data[k] = (T)rhs[k];

    return *this;
  }
  // Necessary because = (const&) operator is deleted when = (&&)
  // operator is defined:
  const Vector & operator =(const Vector & rhs) {
    (*this) = static_cast<const VectorLike<T, evergreen::Vector> &>(rhs);
    return *this;
  }

  const Vector & operator =(Vector && rhs) {
    // Ensure no overlap (as above):
    bool no_overlap = (_data + _length <= rhs._data) || (rhs._data + rhs._length <= _data);
    assert( no_overlap );
    clear();

    std::swap(_length, rhs._length);
    std::swap(_data, rhs._data);

    return *this;
  }
  ~Vector() {
    clear();
  }
  const T & operator [](unsigned long i) const {
    #ifdef BOUNDS_CHECK
    assert(i < size());
    #endif

    return _data[i];
  }
  T & operator [](unsigned long i) {
    #ifdef BOUNDS_CHECK
    assert(i < size());
    #endif

    return _data[i];
  }
  unsigned long size() const {
    return _length;
  }
  void clear() {
    _length = 0;
    if (_data != NULL) {
      free(_data);
      _data = NULL;
    }
  }

  operator const T*const() const {
    return _data;
  }
  operator T*const() const {
    return _data;
  }

  WritableVectorView<T> start_at(unsigned long start) {
    #ifdef SHAPE_CHECK
    assert(start < _length);
    #endif

    return WritableVectorView<T>(*this, start);
  }
  WritableVectorView<T> start_at(unsigned long start, unsigned long length) {
    #ifdef SHAPE_CHECK
    assert(start + length <= _length);
    #endif

    return WritableVectorView<T>(*this, start, length);
  }

  // TODO: should the const versions have different names? The
  // compiler should be able to choose the const version if necessary,
  // so the only reason to say so explicitly is when a Vector is
  // passed by & (non const) and you explicitly want a non-writable
  // view...
  VectorView<T> start_at_const(unsigned long start) const {
    #ifdef SHAPE_CHECK
    assert(start < _length);
    #endif

    return VectorView<T>(*this, start);
  }
  VectorView<T> start_at_const(unsigned long start, unsigned long length) const {
    #ifdef SHAPE_CHECK
    assert(start + length <= _length);
    #endif

    return VectorView<T>(*this, start, length);
  }

  void shrink(unsigned long new_length) {
    #ifdef SHAPE_CHECK
    assert(new_length <= _length);
    #endif

    // Note: This may not be cache aligned; it may be worth
    // investigating if there is a method to perform an
    // aligned_realloc.
    _data = (T*) realloc(_data, new_length*sizeof(T));
    _length = new_length;
  }

  // Only bother creating a && version of this function; if rhs cannot
  // be destroyed by directly moving the pointer, then it would be
  // trivial to simply copy element by element in O(n) time. This O(1)
  // solution is necessarily destructive.

  // Don't use this unless you know what you're doing.
  template <typename S>
  static Vector<T> create_reinterpreted(Vector<S> && rhs) {
    #ifdef SHAPE_CHECK
    assert(rhs._length * sizeof(S) % sizeof(T) == 0);
    #endif

    Vector<T> res;
    res._data = (T*)rhs._data;
    rhs._data = NULL;
    res._length = (rhs._length * sizeof(S)) /sizeof(T);
    rhs._length = 0ul;
    return res;
  }
};

template <typename T>
Vector<T> reversed(const Vector<T> & rhs) {
  Vector<T> result(rhs.size());
  for (unsigned long k=0; k<rhs.size(); ++k)
    result[rhs.size()-1-k] = rhs[k];
  return result;
}

template <typename T>
Vector<T> concatenate(const Vector<T> & lhs, const Vector<T> & rhs) {
  Vector<T> result(lhs.size() + rhs.size());
  for (unsigned long k=0; k<lhs.size(); ++k)
    result[k] = lhs[k];
  for (unsigned long k=0; k<rhs.size(); ++k)
    result[k+lhs.size()] = rhs[k];
  return result;
}

#endif
