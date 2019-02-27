#ifndef _DIF_HPP
#define _DIF_HPP

#include <iostream>
#include <cmath>

#include "RealFFTPostprocessor.hpp"
#include "DIFButterfly.hpp"

template<unsigned char LOG_N, bool SHUFFLE>
class DIF {
public:
  inline static void fft1d(cpx* __restrict const data) {
    DIFButterfly< 1ul << LOG_N >::apply(data);
    if (SHUFFLE)
      RecursiveShuffle<cpx, LOG_N>::apply(data);
  }
  // Note: Here N refers to the full length of the equivalent complex
  // FFT (i.e., the FFT of N cpx values where all imaginary components
  // are 0). The actual allocation should be N/2+1 packed cpx values.
  inline static void real_fft1d_packed(cpx* __restrict const data) {
    const unsigned long int LOG_N_PACKED = LOG_N-1;

    DIF<LOG_N_PACKED, SHUFFLE>::fft1d(data);

    static_assert(SHUFFLE, "DIF on reals must be used with reordered data (the performance cost of inlining the shuffle operations is worse than simply shuffling)");
    RealFFTPostprocessor<LOG_N>::apply(data);
  }
  inline static void real_ifft1d_packed(cpx* __restrict const data) {
    const unsigned long int LOG_N_PACKED = LOG_N-1;
    const unsigned long int N_PACKED = 1ul<<LOG_N_PACKED;

    RealFFTPostprocessor<LOG_N>::apply_inverse(data);
    // Conj.:
    for (unsigned long int k=0; k<=N_PACKED; ++k)
      data[k] = data[k].conj();
    // FFT:
    DIF<LOG_N_PACKED, SHUFFLE>::fft1d(data);
    // Conj.:
    for (unsigned long int k=0; k<=N_PACKED; ++k)
      data[k] = data[k].conj();
    // Scale:
    for (unsigned long int k=0; k<=N_PACKED; ++k)
      data[k] *= 1.0/N_PACKED;
  }
};

template<bool SHUFFLE>
class DIF<0u, SHUFFLE> {
public:
  inline static void fft1d(cpx* __restrict const data) {
  }
  inline static void real_fft1d_packed(cpx* __restrict const data) {
  }
  inline static void real_ifft1d_packed(cpx* __restrict const data) {
  }
};

#endif
