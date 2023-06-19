#ifndef _TENSORVIEW_HPP
#define _TENSORVIEW_HPP

template <typename T>
class Tensor;

// Note: TensorView types are for local, temporary use only (e.g.,
// to allow TRIOT expressions that do not start at the same tuple
// indices). But they should not be stored long term; the tensor
// pointer may change, so long-term use is unsafe. 

template <typename T>
class TensorView : public TensorLike<T, TensorView> {
protected:
  const Tensor<T> & _tensor_ref;
  const unsigned long _flat_start;
  const Vector<unsigned long> _view_shape;
  const unsigned long _flat_size;

  // For constructing a TensorView from another TensorView:
  explicit TensorView(const TensorView<T> & ten_con, const Vector<unsigned long> & start):
    _tensor_ref(ten_con._tensor_ref),
    _flat_start( ten_con._flat_start + tuple_to_index(start, data_shape(), ten_con.dimension()) ),
    _view_shape(ten_con.data_shape() - start),
    _flat_size(flat_length(_view_shape))
  {
    #ifdef SHAPE_CHECK
    // Allows views of size 0:
    assert( start <= ten_con.data_shape() );
    #endif
  }

public:
  // View full Tensor:
  explicit TensorView(const Tensor<T> & ten, const Vector<unsigned long> & start):
    _tensor_ref(ten),
    _flat_start(tuple_to_index(start, data_shape(), ten.dimension())),
    _view_shape(ten.data_shape() - start),
    _flat_size(flat_length(_view_shape))
  {
    #ifdef SHAPE_CHECK
    // Allows views of size 0:
    assert( start <= ten.data_shape() );
    #endif
  }

  // View partial Tensor:
  template <template <typename> class VECTOR_A, template <typename> class VECTOR_B>
  explicit TensorView(const Tensor<T> & ten, const VectorLike<unsigned long, VECTOR_A> & start, const VectorLike<unsigned long, VECTOR_B> & new_view_shape):
    _tensor_ref(ten),
    _flat_start(tuple_to_index(start, data_shape(), ten.dimension())),
    _view_shape(new_view_shape),
    _flat_size(flat_length(_view_shape))
  {
    #ifdef SHAPE_CHECK
    // Allows views of size 0:
    assert( start + new_view_shape <= ten.data_shape() );
    #endif
  }

  // Note: This is used by TRIOT expressions of views, but is
  // deceptive for external use (the view is not guaranteed to be
  // contiguous). One option to address this would be to make the []
  // operator private and making the TRIOT helper classes friend
  // classes; however, access to the [] operator is still desirable to
  // the informed user, so it has not yet been changed.
  const T & operator[] (unsigned long index) const {
    #ifdef BOUNDS_CHECK
    // Note: This is expensive; it may be preferred to simply not
    // allow [] operation outside of TRIOT expressions (as described
    // above).

    // Compute starting tuple and starting index + current index tuple
    // in reference tensor:
    unsigned long* __restrict start_tuple = index_to_tuple(_flat_start, data_shape(), dimension());
    unsigned long* __restrict current_tuple = index_to_tuple(_flat_start+index, data_shape(), dimension());

    // Subtract starting tuple to get tuple corresponding to view:
    for (unsigned char k=0; k<dimension(); ++k)
      current_tuple[k] -= start_tuple[k];

    // Check that every axis is in range for this view:
    for (unsigned char k=0; k<dimension(); ++k)
      assert(current_tuple[k] < view_shape()[k]);
    
    free(start_tuple);
    free(current_tuple);
    #endif

    return _tensor_ref[_flat_start + index];
  }

  // Rewiring other [] operators to TensorLike (unfortunately,
  // the compiler cannot detect the appropriate one on its own):
  const T & operator[](const_tup_t tuple) const {
    return static_cast<const TensorLike<T, evergreen::TensorView> &>(*this)[tuple];
  }
  template <template <typename> class VECTOR>
  const T & operator[](const VectorLike<unsigned long, VECTOR> & tuple) const {
    return static_cast<const TensorLike<T, evergreen::TensorView> &>(*this)[tuple];
  }

  unsigned char dimension() const {
    return _tensor_ref.dimension();
  }
  
  TensorView<T> start_at_const(const Vector<unsigned long> & start) const {
    #ifdef SHAPE_CHECK
    assert(start.size() == dimension());
    assert(start < view_shape());
    #endif

    return TensorView<T>(*this, start);
  }

  const Vector<unsigned long> & data_shape() const {
    return _tensor_ref.data_shape();
  }
  const Vector<unsigned long> & view_shape() const {
    return _view_shape;
  }
  unsigned long flat_size() const {
    return _flat_size;
  }
};

template <typename T>
class WritableTensorView : public WritableTensorLike<T, WritableTensorView> {
protected:
  Tensor<T> & _tensor_ref;
  const unsigned long _flat_start;
  const Vector<unsigned long> _view_shape;
  const unsigned long _flat_size;

  // For constructing a TensorView from another TensorView:
  explicit WritableTensorView(WritableTensorView<T> & ten_con, const Vector<unsigned long> & start):
    _tensor_ref(ten_con._tensor_ref),
    _flat_start( ten_con._flat_start + tuple_to_index(start, data_shape(), ten_con.dimension()) ),
    _view_shape(ten_con.data_shape() - start),
    _flat_size(flat_length(_view_shape))
  {
    #ifdef SHAPE_CHECK
    // Allows views of size 0:
    assert( start <= ten_con.data_shape() );
    #endif
  }

public:
  // View full Tensor:
  explicit WritableTensorView(Tensor<T> & ten, const Vector<unsigned long> & start):
    _tensor_ref(ten),
    _flat_start(tuple_to_index(start, data_shape(), ten.dimension())),
    _view_shape(ten.data_shape() - start),
    _flat_size(flat_length(_view_shape))
  {
    #ifdef SHAPE_CHECK
    // Allows views of size 0:
    assert( start <= ten.data_shape() );
    #endif
  }

  // View partial Tensor:
  template <template <typename> class VECTOR_A, template <typename> class VECTOR_B>
  explicit WritableTensorView(Tensor<T> & ten, const VectorLike<unsigned long, VECTOR_A> & start, const VectorLike<unsigned long, VECTOR_B> & new_view_shape):
    _tensor_ref(ten),
    _flat_start(tuple_to_index(start, data_shape(), ten.dimension())),
    _view_shape(new_view_shape),
    _flat_size(flat_length(_view_shape))
  {
    #ifdef SHAPE_CHECK
    // Allows views of size 0:
    assert( start + new_view_shape <= ten.data_shape() );
    #endif
  }

  // Note: These are used by TRIOT expressions of views, but is
  // deceptive for external use (the view is not guaranteed to be
  // contiguous). See above.
  T & operator[] (unsigned long index) {
    return _tensor_ref[_flat_start + index];
  }
  const T & operator[] (unsigned long index) const {
    return _tensor_ref[_flat_start + index];
  }

  // Rewiring other [] operators to TensorLike (unfortunately,
  // the compiler cannot detect the appropriate one on its own):
  const T & operator[](const_tup_t tuple) const {
    return static_cast<const TensorLike<T, evergreen::WritableTensorView> &>(*this)[tuple];
  }
  T & operator[](const_tup_t tuple) {
    return static_cast<WritableTensorLike<T, evergreen::WritableTensorView> &>(*this)[tuple];
  }
  template <template <typename> class VECTOR>
  const T & operator[](const VectorLike<unsigned long, VECTOR> & tuple) const {
    return static_cast<const TensorLike<T, evergreen::WritableTensorView> &>(*this)[tuple];
  }
  template <template <typename> class VECTOR>
  T & operator[](const VectorLike<unsigned long, VECTOR> & tuple) {
    return static_cast<WritableTensorLike<T, evergreen::WritableTensorView> &>(*this)[tuple];
  }

  unsigned char dimension() const {
    return _tensor_ref.dimension();
  }
  
  WritableTensorView<T> start_at(const Vector<unsigned long> & start) {
    #ifdef SHAPE_CHECK
    assert(start.size() == dimension());
    assert(start < view_shape());
    #endif

    return WritableTensorView<T>(*this, start);
  }
  TensorView<T> start_at_const(const Vector<unsigned long> & start) const {
    #ifdef SHAPE_CHECK
    assert(start.size() == dimension());
    assert(start < view_shape());
    #endif

    return TensorView<T>(*this, start);
  }

  const Vector<unsigned long> & data_shape() const {
    return _tensor_ref.data_shape();
  }
  const Vector<unsigned long> & view_shape() const {
    return _view_shape;
  }
  unsigned long flat_size() const {
    return _flat_size;
  }
};

#endif
