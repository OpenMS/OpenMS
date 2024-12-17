#ifndef _VECTORVIEW_HPP
#define _VECTORVIEW_HPP

template <typename T>
class Vector;


template <typename T>
class VectorView : public VectorLike<T, VectorView> {
protected:
  const Vector<T> & _vec_ref;
  const unsigned long _start;
  const unsigned long _length;
public:
  explicit VectorView(const Vector<T> & vec, unsigned long start):
    _vec_ref(vec),
    _start(start),
    _length(vec.size() - _start)
  {
    #ifdef SHAPE_CHECK
    // Allows views of size 0:
    assert( start <= vec.size() );
    #endif
  }
  explicit VectorView(const Vector<T> & vec, unsigned long start, unsigned long length):
    _vec_ref(vec),
    _start(start),
    _length(length)
  {
    #ifdef SHAPE_CHECK
    // Allows views of size 0:
    assert( _start + _length <= vec.size() );
    #endif
  }
  const T & operator [] (unsigned long i) const {
    #ifdef BOUNDS_CHECK
    assert(i < size());
    #endif

    return _vec_ref[_start + i];
  }
  operator const T*const() const {
    return (const T*const)(_vec_ref) + _start;
  }
  unsigned long size() const {
    return _length;
  }
  VectorView<T> start_at_const(unsigned long start) const {
    return VectorView(_vec_ref, start+_start);
  }
  VectorView<T> start_at_const(unsigned long start, unsigned long length) const {
    return VectorView(_vec_ref, start+_start, length);
  }
};

template <typename T>
class WritableVectorView : public WritableVectorLike<T, WritableVectorView> {
protected:
  Vector<T> & _vec_ref;
  const unsigned long _start;
  const unsigned long _length;
public:
  explicit WritableVectorView(Vector<T> & vec, unsigned long start):
    _vec_ref(vec),
    _start(start),
    _length(vec.size() - _start)
  {
    #ifdef SHAPE_CHECK
    // Allows views of size 0:
    assert( start <= vec.size() );
    #endif
  }
  explicit WritableVectorView(Vector<T> & vec, unsigned long start, unsigned long length):
    _vec_ref(vec),
    _start(start),
    _length(length)
  {
    #ifdef SHAPE_CHECK
    // Allows views of size 0:
    assert( _start + _length <= vec.size() );
    #endif
  }
  template <typename S, template <typename> class VECTOR>
  const WritableVectorView & operator =(const VectorLike<S, VECTOR> & rhs) {
    copy(*this, rhs);
    return *this;
  }
  const T & operator [] (unsigned long i) const {
    #ifdef BOUNDS_CHECK
    assert(i < size());
    #endif

    return _vec_ref[_start + i];
  }
  T & operator [] (unsigned long i) {
    #ifdef BOUNDS_CHECK
    assert(i < size());
    #endif

    return _vec_ref[_start + i];
  }
  operator const T*const() const {
    return (const T*const)(_vec_ref) + _start;
  }
  operator T*const() const {
    return (T*const)(_vec_ref) + _start;
  }
  unsigned long size() const {
    return _length;
  }
  void fill(T value) {
    for (unsigned long k=0; k<_length; ++k)
      (*this)[k] = value;
  }

  WritableVectorView start_at(unsigned long start) {
    return WritableVectorView(_vec_ref, start+_start);
  }
  WritableVectorView start_at(unsigned long start, unsigned long length) {
    return WritableVectorView(_vec_ref, start+_start, length);
  }
  VectorView<T> start_at_const(unsigned long start) const {
    return VectorView<T>(_vec_ref, start+_start);
  }
  VectorView<T> start_at_const(unsigned long start, unsigned long length) const {
    return VectorView<T>(_vec_ref, start+_start, length);
  }
};

#endif
