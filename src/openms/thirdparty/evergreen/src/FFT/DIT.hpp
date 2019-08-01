#ifndef _DIT_HPP
#define _DIT_HPP

#include <iostream>
#include <cmath>

#include "RealFFTPostprocessor.hpp"
#include "DITButterfly.hpp"

template<unsigned char LOG_N, bool SHUFFLE>
class DIT {
public:
  inline static void fft1d(cpx* __restrict const data) {
    if (SHUFFLE)
      RecursiveShuffle<cpx, LOG_N>::apply(data);
    DITButterfly< 1ul << LOG_N >::apply(data);
  }

  // Note: Here N refers to the full length of the equivalent complex
  // FFT (i.e., the FFT of N cpx values where all imaginary components
  // are 0). The actual allocation should be N/2+1/2 packed cpx values.
  inline static void real_fft1d_packed(cpx* __restrict const data) {
    const unsigned long LOG_N_PACKED = LOG_N-1;
    DIT<LOG_N_PACKED, SHUFFLE>::fft1d(data);

    // Data must be in order, so real FFT cleanup is safe:
    RealFFTPostprocessor<LOG_N>::apply(data);
  }
  inline static void real_ifft1d_packed(cpx* __restrict const data) {
    const unsigned long LOG_N_PACKED = LOG_N-1;
    const unsigned long N_PACKED = 1ul<<LOG_N_PACKED;

    static_assert(SHUFFLE, "Inverse DIT on reals must be used with reordered data (the performance cost of inlining the shuffle operations is worse than simply shuffling)");

    RealFFTPostprocessor<LOG_N>::apply_inverse(data);
    // Conj.:
    for (unsigned long k=0; k<=N_PACKED; ++k)
      data[k] = data[k].conj();
    // FFT:
    DIT<LOG_N_PACKED, SHUFFLE>::fft1d(data);
    // Conj.:
    for (unsigned long k=0; k<=N_PACKED; ++k)
      data[k] = data[k].conj();
    // Scale:
    for (unsigned long k=0; k<=N_PACKED; ++k)
      data[k] *= 1.0/N_PACKED;
  }
};

template<bool SHUFFLE>
class DIT<0u, SHUFFLE> {
public:
  inline static void fft1d(cpx* __restrict const data) {
  }
  inline static void real_fft1d_packed(cpx* __restrict const data) {
  }
  inline static void real_ifft1d_packed(cpx* __restrict const data) {
  }
};

#endif
