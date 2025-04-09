#ifndef _TENSORLIKE_HPP
#define _TENSORLIKE_HPP

template <typename T>
class TensorView;

template <typename T>
class WritableTensorView;

// Never instantiate these; always pass by reference & or &&:

template <typename T, template <typename> class TENSOR >
class TensorLike {
public:
  unsigned char dimension() const {
    return static_cast<const TENSOR<T> &>(*this).dimension();
  }
  unsigned long flat_size() const {
    return static_cast<const TENSOR<T> &>(*this).flat_size();
  }
  const T & operator[](unsigned long i) const {
    return static_cast<const TENSOR<T> &>(*this)[i];
  }
  const T & operator[](const_tup_t tuple) const {
    #ifdef BOUNDS_CHECK
    for (unsigned char k=0; k<dimension(); ++k)
      assert( tuple[k] < view_shape()[k] );
    #endif
    return (*this)[ tuple_to_index(tuple, this->data_shape(), this->dimension()) ];
  }
  template <template <typename> class VECTOR>
  const T & operator[](const VectorLike<unsigned long, VECTOR> & tuple) const {
    return (*this)[ static_cast<const_tup_t>(tuple) ];
  }
  const Vector<unsigned long> & data_shape() const {
    return static_cast<const TENSOR<T> &>(*this).data_shape();
  }
  const Vector<unsigned long> & view_shape() const {
    return static_cast<const TENSOR<T> &>(*this).view_shape();
  }
  template <template <typename> class VECTOR>
  TensorView<T> start_at_const(const VectorLike<unsigned long, VECTOR> & start) const {
    return static_cast<const TENSOR<T> &>(*this).start_at_const(start);
  }
  template <template <typename> class VECTOR>
  TensorView<T> start_at_const(const VectorLike<unsigned long, VECTOR> & start, const VectorLike<unsigned long, VECTOR> & new_view_shape) const {
    return static_cast<const TENSOR<T> &>(*this).start_at_const(start, new_view_shape);
  }

  static void print_helper(std::ostream & os, const T*const rhs, const_tup_t data_shape, const_tup_t view_shape, unsigned char dimension) {
    os << "[";
    if (dimension > 1) {
      unsigned long flat_size_without_first = flat_length(data_shape+1, dimension-1);
      for (unsigned long i=0; i<view_shape[0]; ++i) {
	print_helper(os, rhs + i*flat_size_without_first, data_shape+1, view_shape+1, dimension-1);
	if (i != (view_shape[0]-1))
	  os << ", ";
      }
    }
    else {
      for (unsigned long i=0; i<view_shape[0]; ++i) {
	os << rhs[i];
	if (i != (view_shape[0]-1))
	  os << ", ";
      }
    }
    os << "]";
  }
};

template <typename T, template <typename> class TENSOR >
class WritableTensorLike : public TensorLike<T, TENSOR> {
public:
  T & operator[](unsigned long i) {
    return static_cast<TENSOR<T> &>(*this)[i];
  }
  T & operator[](const_tup_t tuple) {
    #ifdef BOUNDS_CHECK
    for (unsigned char k=0; k<this->dimension(); ++k)
      assert( tuple[k] < this->view_shape()[k] );
    #endif
    return (*this)[ tuple_to_index(tuple, this->data_shape(), this->dimension()) ];
  }
  template <template <typename> class VECTOR>
  T & operator[](const VectorLike<unsigned long, VECTOR> & tuple) {
    return (*this)[ static_cast<const_tup_t>(tuple) ];
  }

  template <template <typename> class VECTOR>
  WritableTensorView<T> start_at(const VectorLike<unsigned long, VECTOR> & start) {
    return static_cast<const TENSOR<T> &>(*this).start_at(start);
  }
  template <template <typename> class VECTOR>
  WritableTensorView<T> start_at(const VectorLike<unsigned long, VECTOR> & start, const VectorLike<unsigned long, VECTOR> & new_view_shape) {
    return static_cast<const TENSOR<T> &>(*this).start_at(start, new_view_shape);
  }
};

template <typename T, template <typename> class TENSOR>
std::ostream & operator<<(std::ostream & os, const TensorLike<T, TENSOR> & rhs) {
  // To distinguish 1D Tensor from Vector:
  os << "t:";
  if ( rhs.flat_size() == 0 ) {
    for (unsigned char k=0; k<rhs.dimension(); ++k)
      os << "[";
    for (unsigned char k=0; k<rhs.dimension(); ++k)
      os << "]";
  }
  else
    rhs.print_helper(os, &rhs[0ul], rhs.data_shape(), rhs.view_shape(), rhs.dimension());
  return os;
}

#endif
