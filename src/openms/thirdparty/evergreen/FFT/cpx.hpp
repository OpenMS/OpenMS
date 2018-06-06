#ifndef _CPX_HPP
#define _CPX_HPP

#include <ostream>
#include <math.h>

struct cpx {
  static constexpr double PrintEpsilon = 1e-12;
  double r;
  double i;

  constexpr cpx():
    r(0.0),
    i(0.0)
  { }

  constexpr cpx(double rVal):
    r(rVal),
    i(0.0)
  { }

  constexpr cpx(double rVal, double iVal):
    r(rVal),
    i(iVal)
  { }
  
  __attribute__((always_inline))
  inline cpx operator+=(cpx rhs){
    r += rhs.r;
    i += rhs.i;
    return *this;
  }
  __attribute__((always_inline))
  inline cpx operator-=(cpx rhs){
    r -= rhs.r;
    i -= rhs.i;
    return *this;
  }
  // Slightly faster than * operator (needs only one temporary double):
  __attribute__((always_inline))
  inline const cpx & operator *=(cpx rhs) {
    double temp = r;
    r *= rhs.r;
    r -= i*rhs.i;
    i = temp*rhs.i+i*rhs.r;
    return *this;
  }
  __attribute__((always_inline))
  inline const cpx & operator *=(double scale) {
    r *= scale;
    i *= scale;
    return *this;
  }
  __attribute__((always_inline))
  inline const cpx & operator /=(double denom) {
    denom = 1.0/denom;
    r *= denom;
    i *= denom;
    return *this;
  }
  __attribute__((always_inline))
  inline cpx conj() const {
    return cpx{r, -i};
  }
};

__attribute__((always_inline))
inline cpx operator *(cpx lhs, cpx rhs) {
  // Gauss' method for multiplying complex numbers turns out not to be
  // faster after compiler optimizations; here is the naive method:
  return cpx{lhs.r*rhs.r-lhs.i*rhs.i, lhs.r*rhs.i+lhs.i*rhs.r};
}

__attribute__((always_inline))
inline cpx operator *(double lhs, cpx rhs) {
  rhs.r *= lhs;
  rhs.i *= lhs;
  return rhs;
}

__attribute__((always_inline))
inline cpx operator -(cpx lhs, cpx rhs){
  return cpx{lhs.r-rhs.r, lhs.i-rhs.i};
}

__attribute__((always_inline))
inline cpx operator +(cpx lhs, cpx rhs){
  return cpx{lhs.r+rhs.r, lhs.i+rhs.i};
}

__attribute__((always_inline))
inline cpx operator /(cpx lhs, double rhs) {
  lhs.r /= rhs;
  lhs.i /= rhs;
  return lhs;
}

__attribute__((always_inline))
inline bool operator ==(cpx lhs, cpx rhs){
  return (lhs.r == rhs.r) && (lhs.i == rhs.i);
}

inline std::ostream & operator << (std::ostream & os, cpx cmplx) {
  if ( fabs(cmplx.r) >= cpx::PrintEpsilon && fabs(cmplx.i) >= cpx::PrintEpsilon ) {
    os << cmplx.r;
    if (cmplx.i > 0)
      os << '+';
    os << cmplx.i << 'j' ;
    return os;
  }
  if ( fabs(cmplx.r) >= cpx::PrintEpsilon )
    return ( os << cmplx.r );
  if ( fabs(cmplx.i) >= cpx::PrintEpsilon )
    return (os << cmplx.i << 'j' );
  return (os << 0.0);
}

#endif
