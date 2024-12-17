#ifndef _EMBED_HPP
#define _EMBED_HPP

template <typename T>
class Tensor;

template <typename S, typename T, template <typename> class TENSOR_A, template <typename> class TENSOR_B>
void embed(WritableTensorLike<S, TENSOR_A> & dest, const TensorLike<T, TENSOR_B> & source) {
  #ifdef SHAPE_CHECK
  assert( dest.view_shape() >= source.view_shape() );
  #endif

  apply_tensors([](S & lhs, const T & rhs) {
      lhs = (S)rhs;
    },
    source.view_shape(),
    dest, source);
}

template <typename S, typename T, template <typename> class TENSOR_A, template <typename> class TENSOR_B>
void embed(WritableTensorLike<S, TENSOR_A> && dest, const TensorLike<T,  TENSOR_B> & source) {
  embed(dest, source);
}

#endif
