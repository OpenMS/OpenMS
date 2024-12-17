#ifndef _TENSOR_HPP
#define _TENSOR_HPP

// #include this file to import Vector, Tensor, TRIOT, and everything
// else from this subdirectory.

#include <cmath>
#include "Vector.hpp"
#include "TensorUtils.hpp"
#include "TensorLike.hpp"
#include "TensorView.hpp"
#include "TRIOT.hpp"
#include "transpose.hpp"
#include "embed.hpp"
#include "ArrayShape.hpp"

// Note: Tensor<T> is for simple numeric T types; underneath, it uses
// Vector<T>, which uses aligned_malloc rather than new[], so no
// constructor is called:
template <typename T>
class Tensor : public WritableTensorLike<T, Tensor> {
protected:
  Vector<unsigned long> _data_shape;
  Vector<T> _flat_vector;
public:
  template <typename S>
  friend class Tensor;

  // Inits to dimension 0:
  Tensor()
  { }

  template <typename ARRAY>
  static Tensor<T> from_array(const ARRAY & arr) {
    // ARRAY should be T[A][B]...[Z]:

    // Cast will be unsafe when ARRAY is not T[A][B]...[Z], but if so,
    // the ArrayShape call below should not compiler. // FIXME: compiler -> compile?
    const T*arr_head = (const T*)arr;
    constexpr unsigned long flat_length = sizeof(arr) / sizeof(T);
    
    return Tensor<T>(std::move(ArrayShape<const T &,decltype(arr)>::eval(arr)), std::move(Vector<T>(flat_length, arr_head)));
  }
  explicit Tensor(Vector<unsigned long> && shape):
    _data_shape(std::move(shape)),
    _flat_vector(flat_length(_data_shape, _data_shape.size()))
  {
    #ifdef SHAPE_CHECK
    assert(dimension() <= MAX_TENSOR_DIMENSION && "Tensor dimension is too large; adjust MAX_TENSOR_DIMENSION value");
    #endif
  }
  template <template <typename> class VECTOR>
  explicit Tensor(const VectorLike<unsigned long, VECTOR> & shape):
    _data_shape(shape),
    _flat_vector(flat_length(_data_shape, _data_shape.size()))
  {
    #ifdef SHAPE_CHECK
    assert(dimension() <= MAX_TENSOR_DIMENSION && "Tensor dimension is too large; adjust MAX_TENSOR_DIMENSION value");
    #endif
  }
  template <template <typename> class VECTOR_A, template <typename> class VECTOR_B>
  explicit Tensor(const VectorLike<unsigned long, VECTOR_A> & shape, const VectorLike<T, VECTOR_B> & data):
    _data_shape(shape),
    _flat_vector(data)
  {
    #ifdef SHAPE_CHECK
    assert( flat_size() == flat_length(_data_shape, _data_shape.size()) );
    assert(dimension() <= MAX_TENSOR_DIMENSION && "Tensor dimension is too large; adjust MAX_TENSOR_DIMENSION value");
    #endif
  }
  template <template <typename> class VECTOR>
  explicit Tensor(const VectorLike<unsigned long, VECTOR> & shape, Vector<T> && data):
    _data_shape(shape),
    _flat_vector(std::move(data))
  {
    #ifdef SHAPE_CHECK
    assert( flat_size() == flat_length(_data_shape, _data_shape.size()) );
    assert(dimension() <= MAX_TENSOR_DIMENSION && "Tensor dimension is too large; adjust MAX_TENSOR_DIMENSION value");
    #endif
  }
  template <template <typename> class VECTOR>
  explicit Tensor(const VectorLike<unsigned long, VECTOR> & shape, const Vector<T> & data):
    _data_shape(shape),
    _flat_vector(data)
  {
    #ifdef SHAPE_CHECK
    assert( flat_size() == flat_length(_data_shape, _data_shape.size()) );
    assert(dimension() <= MAX_TENSOR_DIMENSION && "Tensor dimension is too large; adjust MAX_TENSOR_DIMENSION value");
    #endif
  }
  // Note: the following constructor is used to help the compiler
  // detect when you are creating using an initializer list Tensor<T>({1,2}):
  explicit Tensor(Vector<unsigned long> && shape, const Vector<T> & data):
    _data_shape(std::move(shape)),
    _flat_vector(data)
  {
    #ifdef SHAPE_CHECK
    assert( flat_size() == flat_length(_data_shape, _data_shape.size()) );
    assert(dimension() <= MAX_TENSOR_DIMENSION && "Tensor dimension is too large; adjust MAX_TENSOR_DIMENSION value");
    #endif
  }
  explicit Tensor(Vector<unsigned long> && shape, Vector<T> && data):
    _data_shape(std::move(shape)),
    _flat_vector(std::move(data))
  {
    #ifdef SHAPE_CHECK
    assert( flat_size() == flat_length(_data_shape, _data_shape.size()) );
    assert(dimension() <= MAX_TENSOR_DIMENSION && "Tensor dimension is too large; adjust MAX_TENSOR_DIMENSION value");
    #endif
  }
  template <template <typename> class TENSOR>
  Tensor(const TensorLike<T, TENSOR> & rhs):
    _data_shape(rhs.view_shape()),
    _flat_vector(rhs.flat_size()) // FIXME: figure out why allign_malloc results in indirect loss of memory
  {
    embed(*this, rhs);
  }
  Tensor(const Tensor & rhs):
    Tensor( static_cast<const TensorLike<T, evergreen::Tensor> &>(rhs) )
  { }

  Tensor(Tensor<T> && rhs):
    _data_shape( std::move(rhs._data_shape) ),
    _flat_vector( std::move(rhs._flat_vector) )
  { }
  template <template <typename> class TENSOR>
  const Tensor & operator =(const TensorLike<T, TENSOR> & rhs) {
    _data_shape = rhs.data_shape();
    _flat_vector = Vector<T>(flat_length(_data_shape, _data_shape.size()));
    embed(*this, rhs);
    return *this;
  }
  const Tensor & operator =(Tensor<T> && rhs) {
    _data_shape = std::move(rhs._data_shape);
    _flat_vector = std::move(rhs._flat_vector);
    return *this;
  }
  const Tensor & operator =(const Tensor<T> & rhs) {
    *this = static_cast<const TensorLike<T, evergreen::Tensor> &>(rhs);
    return *this;
  }

  // Providing access as a view essentially gives access as a raw flat
  // vector, but prevents assigning tensor.flat() = new_vector, which
  // could allow the flat length and the shape to become inconsistent.
  WritableVectorView<T> flat() {
    return _flat_vector.start_at(0);
  }
  VectorView<T> flat_const() const {
    return _flat_vector.start_at_const(0);
  }
  VectorView<T> flat() const {
    return flat_const();
  }

  T & operator[] (unsigned long index) {
    return _flat_vector[index];
  }
  const T & operator[] (unsigned long index) const {
    return _flat_vector[index];
  }

  // Rewiring other [] operators to TensorLike (unfortunately,
  // the compiler cannot detect the appropriate one on its own):
  const T & operator[](const_tup_t tuple) const {
    return static_cast<const TensorLike<T, evergreen::Tensor> &>(*this)[tuple];
  }
  T & operator[](const_tup_t tuple) {
    return static_cast<WritableTensorLike<T, evergreen::Tensor> &>(*this)[tuple];
  }
  template <template <typename> class VECTOR>
  const T & operator[](const VectorLike<unsigned long, VECTOR> & tuple) const {
    return static_cast<const TensorLike<T, evergreen::Tensor> &>(*this)[tuple];
  }
  template <template <typename> class VECTOR>
  T & operator[](const VectorLike<unsigned long, VECTOR> & tuple) {
    return static_cast<WritableTensorLike<T, evergreen::Tensor> &>(*this)[tuple];
  }

  unsigned char dimension() const {
    return _data_shape.size();
  }

  WritableTensorView<T> start_at(const Vector<unsigned long> & start) {
    #ifdef SHAPE_CHECK
    assert(start.size() == dimension());
    // Shape bounds will be checked by TensorView constructor.
    #endif
    return WritableTensorView<T>(*this, start);
  }
  TensorView<T> start_at_const(const Vector<unsigned long> & start) const {
    #ifdef SHAPE_CHECK
    assert(start.size() == dimension());
    // Shape bounds will be checked by TensorView constructor.
    #endif
    return TensorView<T>(*this, start);
  }
  WritableTensorView<T> start_at(const Vector<unsigned long> & start, const Vector<unsigned long> & new_view_shape) {
    #ifdef SHAPE_CHECK
    assert(start.size() == dimension());
    // Shape bounds will be checked by TensorView constructor.
    #endif
    return WritableTensorView<T>(*this, start, new_view_shape);
  }
  TensorView<T> start_at_const(const Vector<unsigned long> & start, const Vector<unsigned long> & new_view_shape) const {
    #ifdef SHAPE_CHECK
    assert(start.size() == dimension());
    // Shape bounds will be checked by TensorView constructor.
    #endif
    return TensorView<T>(*this, start, new_view_shape);
  }

  const Vector<unsigned long> & data_shape() const {
    return _data_shape;
  }
  const Vector<unsigned long> & view_shape() const {
    return data_shape();
  }
  unsigned long flat_size() const {
    return _flat_vector.size();
  }

  void shrink(const Vector<unsigned long> & new_shape) {
    #ifdef SHAPE_CHECK
    assert(new_shape <= data_shape());
    #endif

    // Moved elements from larger (or equal) indices into smaller
    // (guaranteed because new_shape <= data_shape(), so each dest
    // flattened index must be <= source flattened index; therefore,
    // should not overwrite any data until it's already been copied to
    // the new location.
    enumerate_for_each_tensors([this, &new_shape](const_tup_t counter, const unsigned long dim) {
	unsigned long old_index = tuple_to_index(counter, _data_shape, dim);
	unsigned long new_index = tuple_to_index(counter, new_shape, dim);
	_flat_vector[new_index] = _flat_vector[old_index];
      },
      new_shape
      );
    _data_shape = new_shape;
    _flat_vector.shrink( flat_length(_data_shape, _data_shape.size()) );
  }
  void shrink(const Vector<unsigned long> & start, const Vector<unsigned long> & new_shape) {
    #ifdef SHAPE_CHECK
    assert(new_shape <= data_shape());
    #endif

    // As above but with a start index; the tensor will be seen twice
    // from two different perspectives, but this should be compatible
    // with the restrict view in the underlying vector since since the
    // TensorView stores a reference to the flat vector rather than
    // storing an alternate pointer:
    TensorView<T> view = start_at_const(start);
    enumerate_for_each_tensors([this, &view, &new_shape](const_tup_t counter, const unsigned long dim) {
	unsigned long old_index = tuple_to_index(counter, _data_shape, dim);
	unsigned long new_index = tuple_to_index(counter, new_shape, dim);
	_flat_vector[new_index] = view[old_index];
      },
      new_shape
      );
    _data_shape = new_shape;
    _flat_vector.shrink( flat_length(_data_shape, _data_shape.size()) );
  }
  void reshape(const Vector<unsigned long> & new_shape) {
    #ifdef SHAPE_CHECK
    assert( flat_length(new_shape, new_shape.size()) == flat_size() );
    #endif
    
    _data_shape = new_shape;
  }
  void clear() {
    _flat_vector.clear();
    _data_shape.fill(0ul);
  }
  
  // See Vector<T>::create_reinterpreted for description:
  template <typename S>
  static Tensor<T> create_reinterpreted(Tensor<S> && rhs) {
    #ifdef SHAPE_CHECK
    assert(rhs.flat_size() * sizeof(S) % sizeof(T) == 0);
    #endif

    Tensor<T> res;
    res._flat_vector = Vector<T>::create_reinterpreted(std::move(rhs._flat_vector)); 
    res._data_shape = std::move(rhs._data_shape);
    res._data_shape[res._data_shape.size()-1] *= sizeof(S);
    res._data_shape[res._data_shape.size()-1] /= sizeof(T);
    return res;
  }
};

// TODO: these operators should be written for TensorLike (as in
// VectorComparison.hpp):
template <typename T>
bool operator ==(const Tensor<T> & lhs, const Tensor<T> & rhs) {
  if (lhs.data_shape() != rhs.data_shape())
    return false;
  return lhs.flat() == rhs.flat();
}

#endif
