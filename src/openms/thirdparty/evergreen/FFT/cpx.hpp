#ifndef _CPX_HPP
#define _CPX_HPP

#include <ostream>
// in MSVC M_PI constant is non-standard
#define _USE_MATH_DEFINES 
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
  
  EVERGREEN_ALWAYS_INLINE cpx operator+=(cpx rhs){
    r += rhs.r;
    i += rhs.i;
    return *this;
  }
  EVERGREEN_ALWAYS_INLINE cpx operator-=(cpx rhs){
    r -= rhs.r;
    i -= rhs.i;
    return *this;
  }
  // Slightly faster than * operator (needs only one temporary double):
  EVERGREEN_ALWAYS_INLINE const cpx & operator *=(cpx rhs) {
    double temp = r;
    r *= rhs.r;
    r -= i*rhs.i;
    i = temp*rhs.i+i*rhs.r;
    return *this;
  }
  EVERGREEN_ALWAYS_INLINE const cpx & operator *=(double scale) {
    r *= scale;
    i *= scale;
    return *this;
  }
  EVERGREEN_ALWAYS_INLINE const cpx & operator /=(double denom) {
    denom = 1.0/denom;
    r *= denom;
    i *= denom;
    return *this;
  }
  EVERGREEN_ALWAYS_INLINE cpx conj() const {
    return cpx{r, -i};
  }
};

EVERGREEN_ALWAYS_INLINE cpx operator *(cpx lhs, cpx rhs) {
  // Gauss' method for multiplying complex numbers turns out not to be
  // faster after compiler optimizations; here is the naive method:
  return cpx{lhs.r*rhs.r-lhs.i*rhs.i, lhs.r*rhs.i+lhs.i*rhs.r};
}

EVERGREEN_ALWAYS_INLINE cpx operator *(double lhs, cpx rhs) {
  rhs.r *= lhs;
  rhs.i *= lhs;
  return rhs;
}

EVERGREEN_ALWAYS_INLINE cpx operator -(cpx lhs, cpx rhs){
  return cpx{lhs.r-rhs.r, lhs.i-rhs.i};
}

EVERGREEN_ALWAYS_INLINE cpx operator +(cpx lhs, cpx rhs){
  return cpx{lhs.r+rhs.r, lhs.i+rhs.i};
}

EVERGREEN_ALWAYS_INLINE cpx operator /(cpx lhs, double rhs) {
  lhs.r /= rhs;
  lhs.i /= rhs;
  return lhs;
}

EVERGREEN_ALWAYS_INLINE bool operator ==(cpx lhs, cpx rhs){
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
