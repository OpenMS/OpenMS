#ifndef _REALFFTPOSTPROCESSOR_HPP
#define _REALFFTPOSTPROCESSOR_HPP

#include "cpx.hpp"
#include "Twiddles.hpp"

template<unsigned int LOG_N>
class RealFFTPostprocessor {
public:
  // Note: Here N refers to the full length of the equivalent complex
  // FFT (as above):
  inline static void apply(cpx* const data) {
    const unsigned long N = 1ul<<LOG_N;

    // Postprocessing to convert FFT to equivalent real FFT:
    cpx bias = data[0];
    data[0].r = bias.r + bias.i;
    data[0].i = 0.0;

    data[N/2].r = bias.r - bias.i;
    data[N/2].i = 0.0;

    cpx current_twiddle{1.0,0.0};
    Twiddles<N/2>::advance(current_twiddle);
    unsigned long k;
    for (k=1; k<=N/4; ++k) {
      cpx x1 = 0.5*(data[k] + data[N/2-k].conj());
      cpx x2 = 0.5*(data[k] - data[N/2-k].conj());

      cpx temp = x2 * cpx{ current_twiddle.i, -current_twiddle.r };

      data[k] = x1 + temp;
      data[N/2-k] = (x1 - temp).conj();

      Twiddles<N/2>::advance(current_twiddle);
      //current_twiddle = current_twiddle * first_twiddle;
    }
  }
  
  inline static void apply_inverse(cpx* const data) {
    const unsigned long N = 1ul<<LOG_N;

    // Use 1D FFT result on packed to compute full FFT:
    cpx bias = data[0];
    cpx last = data[N/2];

    // Note: final element is not used, and so does not need to be
    // created; only create the first element.
    data[0].r = (bias.r + last.r)/2.0;
    data[0].i = (bias.r - last.r)/2.0;

    // Unnecessary, but tidy:
    data[N/2] = cpx{0.0,0.0};

    cpx current_twiddle{1.0,0.0};
    Twiddles<N/2>::advance(current_twiddle);
    unsigned long k;
    for (k=1; k<=N/4; ++k) {
      cpx from_back = data[N/2-k].conj();
      cpx x1 = 0.5*(data[k] + from_back);
      cpx temp = 0.5*(data[k] - from_back);

      // Perform cpx x2 = temp / currentTwiddle (where twiddle is
      // first swapped to (i,-r) to match the forward version
      // apply(...), implemented above):
      cpx x2 = temp * cpx{current_twiddle.i, current_twiddle.r};

      // Important: store data[k] after data[N/2-k] so that when k=N/4
      // and both indices are the same, the data[k] version is used
      // (it will give the correct imaginary sign):
      data[N/2-k] = (x1 - x2).conj();
      data[k] = x1 + x2;

      Twiddles<N/2>::advance(current_twiddle);
      //      current_twiddle = current_twiddle * first_twiddle;
    }
  }
};

template<>
class RealFFTPostprocessor<0u> {
public:
  // Note: Here N refers to the full length of the equivalent complex
  // FFT (as above):
  inline static void apply(cpx* const data) {
  }

  inline static void apply_inverse(cpx* const data) {
  }
};

// Note: could specialize real postprocessing for small lengths. The
// compiler may be able to do this automatically by simply unrolling
// the loop, since the cpx operations are forced to be inline;
// however, with FFT butterflies, explicitly specializing yielded a
// small performance boost, so it may be similar here.

#endif
