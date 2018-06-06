#ifndef TWIDDLES_H
#define TWIDDLES_H

#define _USE_MATH_DEFINES
#include <cmath>

#include "cpx.hpp"

template<unsigned long int N>
class Twiddles {
public:
  inline static void advance(cpx & current) {
    current += current * delta();
  }
  inline static constexpr cpx delta() {
    // Used in the recurrence
    // current_twiddle += current_twiddle * delta,
    // which will produce the next twiddle on the unit circle.

    // Note that it would be mathematically equivalent to perform
    // current_twiddle *= simple_delta
    // where simple_delta = cpx{cos(), -sin()}.

    // However, when the FFT size is large, cos(2*pi / N) \approx 1,
    // and so signal is lost by using the
    // cosine. -2.0*squared(Twiddles<N*2>::sin()) equals 1-cos(2*pi/N)
    // (see below). However, to get the numeric benefit, the 1 must
    // not be added (leaving -cos(2*pi/N)). Thus advancing
    // current_twiddle = current_twiddle * naive_delta is equivalent
    // to assigning current_twiddle = current_twiddle * (delta-1) +
    // current_twiddle which is the same as current_twiddle +=
    // current_twiddle * delta.

    // This final form leaves the real and imag parts of the delta
    // close to zero but applies a recurrence equivalent to the naive.
    return cpx{-2.0*squared(Twiddles<N*2>::sin()), -sin()};

    // Derivation of why cos(theta)-1 = -2*sin(theta/2)**2:
    // Using naive recurrence above:
    // ( cos(theta/2) + sin(theta/2)j ) **2 -->
    // ( cos(theta) + sin(theta)j ).
    // Thus cos(theta) = cos(theta/2)**2 - sin(theta/2)**2.
    // Also, 1 = cos(theta/2)**2 + sin(theta/2)**2.
    // subtracting the 1 =... equation from the cos(theta)
    // =... equation yields
    // cos(theta)-1 = -2*sin(theta/2)**2.
  }
  inline static constexpr double sin() {
    return ::sin(M_PI/N);
  }
  inline static constexpr double cos() {
    //    return ::cos(M_PI/N);
    // Reduces to sin() calls:
    return 1-2.0*squared(Twiddles<N*2>::sin());
  }
protected:
  inline static constexpr double squared(double x) {
    return x*x;
  }
};

#endif
