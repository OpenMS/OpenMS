#ifndef _FFT_HPP
#define _FFT_HPP

// #include this file to import all FFT utilities from this
// subdirectory.

const unsigned char FFT1D_MAX_LOG_N = 31u;

#include "../BitReversedShuffle/RecursiveShuffle.hpp"

#include "DIF.hpp"
#include "DIT.hpp"
#include "../Tensor/Tensor.hpp"
#include "shape_to_log_shape.hpp"

#include <assert.h>

// Note: template-recursive method may no longer be faster because you
// could use the OpenMP SIMD pragma in an iterative implementation

// Note that these are chosen to use Tensor rather than a TensorView
// or WritableTensorView; even though views do not support striding,
// they do allow the memory to be non-contiguous (in the case of
// multidimensional FFT). For FFT performance, the data should be in a
// contiguous block of memory.

// FFT1D should be either DIF or DIT:
template <template <unsigned char, bool> class FFT1D, bool SHUFFLE, bool UNDO_TRANSPOSE>
class NDFFTEnvironment {
public:
  
  // Wrapping in NDFFTEnvironment allows the following classes to only
  // depend on a single template parameter, LOG_NEXT_FROM_RIGHT;
  // therefore, it can be searched with LinearTemplateSearch:

  template<unsigned char LOG_N>
  class SingleFFT1D {
  public:
    inline static void apply(cpx* __restrict data) {
      FFT1D<LOG_N, SHUFFLE>::fft1d(data);
    }
  };

  template<unsigned char LOG_N>
  class SingleIFFT1D {
  public:
    inline static void apply(cpx* __restrict data) {
      for (unsigned long k=0; k<(1ul<<LOG_N); ++k)
	data[k] = data[k].conj();
      FFT1D<LOG_N, SHUFFLE>::fft1d(data);
      const double scale = 1.0 / (1ul<<LOG_N);
      for (unsigned long k=0; k<(1ul<<LOG_N); ++k) {
	data[k] = data[k].conj();
	data[k] *= scale;
      }
    }
  };

  template<unsigned char LOG_N>
  class RowFFTs {
  public:
    inline static void apply(cpx* __restrict data, const unsigned long flat, const bool freshly_zero_padded) {
      // Row FFTs:
      unsigned long k;
      for (k=0; k<(flat>>1); k+=(1ul<<LOG_N))
      	FFT1D<LOG_N, SHUFFLE>::fft1d(data+k);
      if ( ! freshly_zero_padded )
	for (; k<flat; k+=(1ul<<LOG_N))
	  FFT1D<LOG_N, SHUFFLE>::fft1d(data+k);
    }
  };

  // Note: not needed for ND IFFT (which conjugates, perform ND FFTs,
  // conjugates, and scales), but used for consistency so that 1D FFTs
  // have an easy access point for performing 1D FFTs and IFFTs (with
  // only one template parameter, LOG_N, allowing use of
  // LinearTemplateSearch); this can reduce overhead when a 1D FFT or
  // iFFT is needed.
  template<unsigned char LOG_N>
  class RowIFFTs {
  public:
    inline static void apply(cpx* __restrict data, const unsigned long flat) {
      for (unsigned long k=0; k<flat; ++k)
	data[k] = data[k].conj();

      RowFFTs<LOG_N>::apply(data, flat, false);

      const double scale = 1.0 / flat;
      for (unsigned long k=0; k<flat; ++k) {
	data[k] = data[k].conj();
	data[k] *= scale;
      }
    }
  };

  // Note: it may be possible to exploit freshly_zero_padded for a
  // speedup in this transposition code as well:
  template<unsigned char LOG_NEXT_FROM_RIGHT>
  static void transpose_so_next_dimension_becomes_row(cpx* __restrict & data, cpx* __restrict & buffer, const unsigned long flat, const unsigned long prod_shape_from_right) {
    if ((1ul<<LOG_NEXT_FROM_RIGHT) > 1 && prod_shape_from_right > 1) {
      for (unsigned long k=0; k<flat; k+=(1ul<<LOG_NEXT_FROM_RIGHT)*prod_shape_from_right)
	MatrixTranspose<cpx>::apply_buffered(buffer + k, data + k, 1ul<<LOG_NEXT_FROM_RIGHT, prod_shape_from_right);
      std::swap(data, buffer);
    }
  }

  template<unsigned char LOG_NEXT_FROM_RIGHT>
  static void undo_transpositions(cpx* __restrict & data, cpx* __restrict & buffer, const unsigned long flat, const unsigned long prod_shape_from_right) {
    if ((1ul<<LOG_NEXT_FROM_RIGHT) > 1 && prod_shape_from_right > 1) {
      for (unsigned long k=0; k<flat; k+=(1ul<<LOG_NEXT_FROM_RIGHT)*prod_shape_from_right)
	MatrixTranspose<cpx>::apply_buffered(buffer + k, data + k, prod_shape_from_right, 1ul<<LOG_NEXT_FROM_RIGHT);
      std::swap(data, buffer);
    }
  }

  template<unsigned char LOG_NEXT_FROM_RIGHT>
  class RowFFTsAndTransposes {
  public:
    inline static void apply(cpx* __restrict & data, cpx* __restrict & buffer, const unsigned long flat, const unsigned long prod_shape_from_right) {

      // Transpose:
      transpose_so_next_dimension_becomes_row<LOG_NEXT_FROM_RIGHT>(data, buffer, flat, prod_shape_from_right);

      // Row FFTs:
      RowFFTs<LOG_NEXT_FROM_RIGHT>::apply(data, flat, false);

      // Undo transpose:
      if (UNDO_TRANSPOSE)
	undo_transpositions<LOG_NEXT_FROM_RIGHT>(data, buffer, flat, prod_shape_from_right);
    }
  };

  template<unsigned char LOG_N>
  class SingleRealFFT1D {
  public:
    inline static void apply(cpx* __restrict data) {
      FFT1D<LOG_N, SHUFFLE>::real_fft1d_packed(data);
    }
  };

  template<unsigned char LOG_N>
  class SingleRealIFFT1D {
  public:
    inline static void apply(cpx* __restrict data) {
      FFT1D<LOG_N, SHUFFLE>::real_ifft1d_packed(data);
    }
  };

  template<unsigned char LOG_N>
  class RealRowFFTs {
  public:
    inline static void apply(cpx* __restrict data, const unsigned long flat, const bool freshly_zero_padded) {
      // This will only be performed on final axis, so do not bother
      // transposing.

      // Row FFTs:
      const unsigned long SINGLE_FFT_LEN = real_length_to_packed_length(1ul<<LOG_N);

      unsigned long k;
      for (k=0; k<flat>>1; k+=SINGLE_FFT_LEN)
	FFT1D<LOG_N, SHUFFLE>::real_fft1d_packed(data+k);
      if ( ! freshly_zero_padded )
	for (; k<flat; k+=SINGLE_FFT_LEN)
	  FFT1D<LOG_N, SHUFFLE>::real_fft1d_packed(data+k);
    }
  };

  template<unsigned char LOG_N>
  class RealRowIFFTs {
  public:
    inline static void apply(cpx* __restrict data, const unsigned long flat) {
      // This will only be performed on final axis, so do not bother
      // transposing.

      // Row IFFTs:
      const unsigned long SINGLE_FFT_LEN = real_length_to_packed_length(1ul<<LOG_N);
      for (unsigned long k=0; k<flat; k+=SINGLE_FFT_LEN)
	FFT1D<LOG_N, SHUFFLE>::real_ifft1d_packed(data+k);
    }
  };

  // This class is for performing a full multidimensional FFT by
  // moving right to left, performing row FFTs and transposing as
  // necessary:
  class NDFFT {
  public:
    inline static void fft(cpx* __restrict & data, cpx* __restrict & buffer, const unsigned char* __restrict const log_shape, int dimension, const bool freshly_zero_padded) {
      unsigned long flat = 1ul<<sum(log_shape, dimension);
      unsigned long prod_shape_from_right=1;
      if (dimension > 0) {
	// Apply last axis using freshly_zero_padded:
	LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, RowFFTs >::apply(log_shape[dimension-1], data, flat, freshly_zero_padded);
	prod_shape_from_right *= (1ul<<log_shape[dimension-1]);

	// Other axes no longer guarantee freshly_zero_padded:
	--dimension;
	for (; dimension>0; --dimension) {
	  LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, RowFFTsAndTransposes >::apply(log_shape[dimension-1], data, buffer, flat, prod_shape_from_right);
	  prod_shape_from_right *= (1ul<<log_shape[dimension-1]);
	}
      }
    }
    inline static void ifft(cpx* __restrict & data, cpx* __restrict & buffer, const unsigned char* __restrict const log_shape, int dimension) {
      unsigned long flat = 1ul<<sum(log_shape, dimension);
      
      // Conjugate:
      for (unsigned long k=0; k<flat; ++k)
	data[k] = data[k].conj();
      
      // FFT:
      fft(data, buffer, log_shape, dimension, false);
      
      // Conjugate and scale:
      const double scale = 1.0 / flat;
      for (unsigned long k=0; k<flat; ++k) {
	data[k] = data[k].conj();
	data[k] *= scale;
      }
    }

    inline static void real_fft_packed(cpx* __restrict & data, cpx* __restrict & buffer, const unsigned char* __restrict const log_shape, int dimension, const bool freshly_zero_padded) {
      unsigned long prod_shape_from_right = real_length_to_packed_length(1ul<<log_shape[dimension-1]);
      unsigned long flat = (1ul<<sum(log_shape, dimension-1))*prod_shape_from_right;

      // axis dimension-1 performed on reals:
      if (dimension > 0) {
	--dimension;
	// Final axis can use freshly_zero_padded:

	// Force SHUFFLE=true:
	LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, NDFFTEnvironment<FFT1D, true, UNDO_TRANSPOSE>::template RealRowFFTs>::apply(log_shape[dimension], data, flat, freshly_zero_padded);

	// Remaining axes use freshly_zero_padded = false:
	for (; dimension>0; --dimension) {
	  LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, RowFFTsAndTransposes >::apply(log_shape[dimension-1], data, buffer, flat, prod_shape_from_right);
	  prod_shape_from_right *= (1ul<<log_shape[dimension-1]);
	}
      }
    }

    inline static void real_ifft_packed(cpx* __restrict & data, cpx* __restrict & buffer, const unsigned char* __restrict const log_shape, int dimension) {
      unsigned long prod_shape_from_right, real_axis, flat;
      if (UNDO_TRANSPOSE) {
	real_axis = real_length_to_packed_length(1ul<<log_shape[dimension-1]);
	prod_shape_from_right = real_axis;
	flat = (1ul<<sum(log_shape, dimension-1))*prod_shape_from_right;
      }
      else {
	real_axis = real_length_to_packed_length(1ul<<log_shape[0]);
	prod_shape_from_right = 1ul;
	flat = real_axis * (1ul<<sum(log_shape+1, dimension-1));
      }

      // Real axis will been scaled; therefore, scale with length
      // flat/real_axis to fix remaining axes:
      const double scale = real_axis / double(flat); // 1.0 / length

      for (unsigned long k=0; k<flat; ++k)
	data[k] = data[k].conj();

      if (UNDO_TRANSPOSE) {
	// Since inverse on real axis is not performed first,
	// conjugate:
	for (unsigned char dim=dimension-1; dim>=1; --dim) {
	  LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, RowFFTsAndTransposes >::apply(log_shape[dim-1], data, buffer, flat, prod_shape_from_right);
	  prod_shape_from_right *= (1ul<<log_shape[dim-1]);
	}

	for (unsigned long k=0; k<flat; ++k) {
	  data[k] = data[k].conj();
	  data[k] *= scale;
	}

	// Final axis was originally reals; however, it should be
	// performed last:

	// Force SHUFFLE=true:
	LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, NDFFTEnvironment<FFT1D, true, UNDO_TRANSPOSE>::template RealRowIFFTs>::apply(log_shape[dimension-1], data, flat);
      }
      else {
	// When UNDO_TRANSPOSE is false, the first axis will contain
	// packed reals:
	for (unsigned char dim=dimension-1; dim>=1; --dim) {
	  LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, RowFFTsAndTransposes >::apply(log_shape[dim], data, buffer, flat, prod_shape_from_right);
	  prod_shape_from_right *= (1ul<<log_shape[dim]);
	}
	for (unsigned long k=0; k<flat; ++k) {
	  data[k] = data[k].conj();
	  data[k] *= scale;
	}

	// The order of all other axes has been reversed; swap the
	// first axis (the packed reals) with the block of all other
	// axes:
	if (real_axis > 1 && prod_shape_from_right > 1) {
	  MatrixTranspose<cpx>::apply_buffered(buffer, data, real_axis, prod_shape_from_right);
	  std::swap(data, buffer);
	}

	// Force SHUFFLE=true:
	LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, NDFFTEnvironment<FFT1D, true, UNDO_TRANSPOSE>::template RealRowIFFTs>::apply(log_shape[0], data, flat);
      }
    }
  };
};

template <template <unsigned char, bool> class FFT1D, bool SHUFFLE, bool UNDO_TRANSPOSE, bool FORWARD_FFT, bool FRESHLY_ZERO_PADDED>
inline void execute_fft(Tensor<cpx> & ten) {
  Vector<unsigned char> log_shape = shape_to_log_shape(ten.data_shape());
  cpx* __restrict buffer_a = &ten[0ul];

  // A buffer is necessary to perform transpositions:
  Tensor<cpx> buffer(ten.data_shape());
    
  cpx* __restrict buffer_b = &buffer[0ul];

  if (FORWARD_FFT)
    NDFFTEnvironment<FFT1D, SHUFFLE, UNDO_TRANSPOSE>::NDFFT::fft(buffer_a, buffer_b, &log_shape[0], ten.dimension(), FRESHLY_ZERO_PADDED);
  else
    NDFFTEnvironment<FFT1D, SHUFFLE, UNDO_TRANSPOSE>::NDFFT::ifft(buffer_a, buffer_b, &log_shape[0], ten.dimension());

  // buffer_a is the destination; if, after swapping, buffer_a is
  // not the same as the pointer used in ten, then swap:
  if (buffer_a != &ten[0ul])
    ten = std::move(buffer);

  if ( ! UNDO_TRANSPOSE )
    // Axes have been reversed:
    ten.reshape( reversed(ten.data_shape()) );
}

template <template <unsigned char, bool> class FFT1D, bool SHUFFLE, bool UNDO_TRANSPOSE, bool FORWARD_FFT, bool FRESHLY_ZERO_PADDED>
void execute_real_fft_packed(Tensor<cpx> & ten) {
  Vector<unsigned char> log_shape;
  if (UNDO_TRANSPOSE || FORWARD_FFT)
    // packed reals are on final axis:
    log_shape = packed_shape_to_log_shape(ten.data_shape());
  else
    // packed reals are on first axis:
    log_shape = reversed_packed_shape_to_log_shape(ten.data_shape());
    
  cpx* __restrict buffer_a = &ten[0ul];
  // A buffer is necessary to perform transpositions:
  Tensor<cpx> buffer(ten.data_shape());
    
  cpx* __restrict buffer_b = &buffer[0ul];

  if (FORWARD_FFT)
    NDFFTEnvironment<FFT1D, SHUFFLE, UNDO_TRANSPOSE>::NDFFT::real_fft_packed(buffer_a, buffer_b, &log_shape[0], ten.dimension(), FRESHLY_ZERO_PADDED);
  else {
    NDFFTEnvironment<FFT1D, SHUFFLE, UNDO_TRANSPOSE>::NDFFT::real_ifft_packed(buffer_a, buffer_b, &log_shape[0], ten.dimension());
  }

  if (buffer_a != &ten[0ul])
    ten = std::move(buffer);

  if ( ! UNDO_TRANSPOSE )
    // Axes have been reversed:
    ten.reshape( reversed(ten.data_shape()) );
}

template <template <unsigned char, bool> class FFT1D, bool SHUFFLE, bool UNDO_TRANSPOSE, bool FRESHLY_ZERO_PADDED=false>
inline void apply_fft(Tensor<cpx> & ten) {
  if (ten.dimension() == 0 || ten.flat_size() == 0) {
  }
  else if (ten.dimension() == 1)
    LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, NDFFTEnvironment<FFT1D, SHUFFLE, false>::template SingleFFT1D >::apply(log2(ten.flat_size()), &ten[0ul]);
  else
    execute_fft<FFT1D, SHUFFLE, UNDO_TRANSPOSE, true, FRESHLY_ZERO_PADDED>(ten);
}

template <template <unsigned char, bool> class FFT1D, bool SHUFFLE, bool UNDO_TRANSPOSE>
void apply_ifft(Tensor<cpx> & ten) {
  if (ten.dimension() == 0 || ten.flat_size() == 0) {
  }
  else if (ten.dimension() == 1)
    LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, NDFFTEnvironment<FFT1D, SHUFFLE, false>::template SingleIFFT1D >::apply(log2(ten.flat_size()), &ten[0ul]);
  else
    execute_fft<FFT1D, SHUFFLE, UNDO_TRANSPOSE, false, false>(ten);
}

template <template <unsigned char, bool> class FFT1D, bool SHUFFLE, bool UNDO_TRANSPOSE, bool FRESHLY_ZERO_PADDED=false>
Tensor<cpx> fft(Tensor<cpx> ten) {
  Tensor<cpx> ten_prime = ten;
  apply_fft<FFT1D, SHUFFLE, UNDO_TRANSPOSE, FRESHLY_ZERO_PADDED>(ten);
  return ten;
}

template <template <unsigned char, bool> class FFT1D, bool SHUFFLE, bool UNDO_TRANSPOSE>
Tensor<cpx> ifft(Tensor<cpx> ten) {
  apply_ifft<FFT1D, SHUFFLE, UNDO_TRANSPOSE>(ten);
  return ten;
}

template <template <unsigned char, bool> class FFT1D, bool SHUFFLE, bool UNDO_TRANSPOSE, bool FRESHLY_ZERO_PADDED=false>
void apply_real_fft_packed(Tensor<cpx> & ten) {
  if (ten.dimension() == 0 || ten.flat_size() == 0) {
  }
  else if (ten.dimension() == 1)
    LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, NDFFTEnvironment<FFT1D, true, UNDO_TRANSPOSE>::template SingleRealFFT1D >::apply(integer_log2( packed_length_to_real_length(ten.flat_size()) ), &ten[0ul]);
  else
    execute_real_fft_packed<FFT1D, SHUFFLE, UNDO_TRANSPOSE, true, FRESHLY_ZERO_PADDED>(ten);
}

template <template <unsigned char, bool> class FFT1D, bool SHUFFLE, bool UNDO_TRANSPOSE>
void apply_real_ifft_packed(Tensor<cpx> & ten) {
  if (ten.dimension() == 0 || ten.flat_size() == 0) {
  }
  else if (ten.dimension() == 1)
    LinearTemplateSearch<0u, FFT1D_MAX_LOG_N, NDFFTEnvironment<FFT1D, true, UNDO_TRANSPOSE>::template SingleRealIFFT1D >::apply(integer_log2( packed_length_to_real_length(ten.flat_size()) ), &ten[0ul]);
  else
    execute_real_fft_packed<FFT1D, SHUFFLE, UNDO_TRANSPOSE, false, false>(ten);
}

template <template <unsigned char, bool> class FFT1D, bool SHUFFLE, bool UNDO_TRANSPOSE, bool FRESHLY_ZERO_PADDED=false>
Tensor<cpx> real_fft(const Tensor<double> & ten) {
  if (ten.dimension() == 0)
    return Tensor<cpx>();

  Vector<unsigned long> shape = ten.data_shape();
  shape[shape.size()-1] = real_length_to_packed_length(shape[shape.size()-1])*2;
  Tensor<double> larger(std::move(shape));

  embed(larger, ten);
  Tensor<cpx> packed = Tensor<cpx>::create_reinterpreted(std::move(larger));

  apply_real_fft_packed<FFT1D, SHUFFLE, UNDO_TRANSPOSE, FRESHLY_ZERO_PADDED>(packed);
  return packed;
}

template <template <unsigned char, bool> class FFT1D, bool SHUFFLE, bool UNDO_TRANSPOSE>
Tensor<double> real_ifft(const Tensor<cpx> & ten) {
  if (ten.dimension() == 0)
    return Tensor<double>();

  Tensor<cpx> larger = ten;
  apply_real_ifft_packed<FFT1D, SHUFFLE, UNDO_TRANSPOSE>(larger);
  Tensor<double> larger_reals = Tensor<double>::create_reinterpreted(std::move(larger));

  Vector<unsigned long> shape;
  if (UNDO_TRANSPOSE)
    shape = ten.data_shape();
  else
    shape = reversed(ten.data_shape());
  shape[shape.size()-1] = packed_length_to_real_length(shape[shape.size()-1]);
  Tensor<double> smaller(std::move(shape));

  // Cut slice of larger into smaller (this could also be done by
  // creating a TensorView of larger from {0,0,...} to
  // smaller.data_shape() and then calling embed):
  apply_tensors([](double & small_val, double large_val){
      small_val = large_val;
    },
    smaller.data_shape(),
    smaller, larger_reals);

  return smaller;
}

// Note: there is a small unexploited speedup that would detect when
// any axes have length 1, and would perform the equivalent FFT with
// those axes removed. The result will be the same, but it will be
// faster.

#endif
