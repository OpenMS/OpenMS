#ifndef _FFT_CONVOLVE_HPP
#define _FFT_CONVOLVE_HPP

#include "naive_convolve.hpp"
#include "../FFT/FFT.hpp"

inline unsigned char log2_ceiling(unsigned long len){
  return (unsigned char)ceil(log2(len));
};

inline unsigned long power_of_2_ceiling(unsigned long len){
  return 1ul<<log2_ceiling(len);
};

inline Tensor<cpx> fft_convolve(const Tensor<cpx> & lhs, const Tensor<cpx> & rhs) {
  #ifdef SHAPE_CHECK
  assert(lhs.dimension() == rhs.dimension());
  assert(lhs.data_shape() + rhs.data_shape() >= 1ul);
  #endif
  if (lhs.dimension() == 0)
    return Tensor<cpx>();

  unsigned long k;

  Vector<unsigned long> conv_shape(lhs.dimension());
  for (k=0; k<lhs.dimension(); ++k) {
    unsigned long larger = std::max(lhs.data_shape()[k], rhs.data_shape()[k]);
    conv_shape[k] = power_of_2_ceiling(larger) * 2;
  }

  Tensor<cpx> lhs_padded(conv_shape);
  embed(lhs_padded, lhs);
  Tensor<cpx> rhs_padded(conv_shape);
  embed(rhs_padded, rhs);

  apply_fft<DIF, false, false, true>(lhs_padded);
  apply_fft<DIF, false, false, true>(rhs_padded);

  lhs_padded.flat() *= rhs_padded.flat();
  
  // Allow rhs_padded to deallocate:
  rhs_padded.clear();

  // Perform in-place inverse FFT on packed values:
  apply_ifft<DIT, false, false>(lhs_padded);

  lhs_padded.shrink(lhs.data_shape() + rhs.data_shape() - 1ul);
  return lhs_padded;
}

inline Vector<unsigned long> padded_convolution_shape(const Tensor<double> & lhs, const Tensor<double> & rhs) {

  unsigned long k;

  Vector<unsigned long> conv_shape_doubles(lhs.dimension());
  #ifdef SHAPE_CHECK
  assert(lhs.dimension() > 0);
  #endif
  for (k=0; k<lhs.dimension()-1u; ++k) {
    unsigned long larger = std::max(lhs.data_shape()[k], rhs.data_shape()[k]);
    conv_shape_doubles[k] = power_of_2_ceiling(larger) * 2;
  }
  // Final axis is n/2+1 cpx values, after *2 it becomes n+1 cpx
  // values, which will be 2*(n+1) double values:
  conv_shape_doubles[k] = 2*(power_of_2_ceiling(std::max(lhs.data_shape()[k], rhs.data_shape()[k])) + 1);
  return conv_shape_doubles;
}

inline Tensor<double> fft_convolve_already_padded_rvalue(Tensor<double> && lhs_padded_doubles, Tensor<double> && rhs_padded_doubles, Vector<unsigned long> result_shape) {
  #ifdef SHAPE_CHECK
  assert(lhs_padded_doubles.dimension() == rhs_padded_doubles.dimension());
  assert(lhs_padded_doubles.data_shape() + rhs_padded_doubles.data_shape() >= 1ul);
  #endif
  if (lhs_padded_doubles.dimension() == 0)
    return Tensor<double>();

  Tensor<cpx> lhs_padded = Tensor<cpx>::create_reinterpreted(std::move(lhs_padded_doubles));
  Tensor<cpx> rhs_padded = Tensor<cpx>::create_reinterpreted(std::move(rhs_padded_doubles));

  apply_real_fft_packed<DIF, false, false, true>(lhs_padded);
  apply_real_fft_packed<DIF, false, false, true>(rhs_padded);

  lhs_padded.flat() *= rhs_padded.flat();
  
  // Allow rhs_padded to deallocate:
  rhs_padded.clear();

  // Perform in-place inverse FFT on packed values:
  apply_real_ifft_packed<DIT, false, false>(lhs_padded);

  // Unpack:
  Tensor<double> result = Tensor<double>::create_reinterpreted(std::move(lhs_padded));

  result.shrink(result_shape);
  return result;
}

inline Tensor<double> fft_convolve(const Tensor<double> & lhs, const Tensor<double> & rhs) {
  #ifdef SHAPE_CHECK
  assert(lhs.dimension() == rhs.dimension());
  assert(lhs.data_shape() + rhs.data_shape() >= 1ul);
  #endif
  if (lhs.dimension() == 0)
    return Tensor<double>();

  Vector<unsigned long> conv_shape_doubles = padded_convolution_shape(lhs, rhs);

  Tensor<double> lhs_padded_doubles(conv_shape_doubles);
  embed(lhs_padded_doubles, lhs);
  Tensor<double> rhs_padded_doubles(conv_shape_doubles);
  embed(rhs_padded_doubles, rhs);

  return fft_convolve_already_padded_rvalue(std::move(lhs_padded_doubles), std::move(rhs_padded_doubles), lhs.data_shape() + rhs.data_shape() - 1ul);
}

#endif
