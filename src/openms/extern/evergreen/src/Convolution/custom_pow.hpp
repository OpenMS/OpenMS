#ifndef _CUSTOM_POW_HPP
#define _CUSTOM_POW_HPP

// greater numeric stability than the one listed below:
inline double fast_pow(double a, const double b) {
  // calculate approximation with just the fraction of the exponent
  int exp = (int) b;
  union {
    double d;
    int x[2];
  } u = { a };
  u.x[1] = (int)((b - exp) * (u.x[1] - 1072632447) + 1072632447);
  u.x[0] = 0;
 
  // exponentiation by squaring with the exponent's integer part
  // double r = u.d makes everything much slower, not sure why
  double r = 1.0;
  while (exp) {
    if (exp & 1) {
      r *= a;
    }
    a *= a;
    exp >>= 1;
  }
 
  return r * u.d;
}


inline double faster_pow(const double a, const double b) {
  union {
    double d;
    struct {
      int a;
      int b;
    } s;
  } u = { a };
  u.s.b = (int)(b * (u.s.b - 1072632447) + 1072632447);
  u.s.a = 0;
  return u.d;
}

inline double custom_pow(double a, const double b) {
  #ifdef FASTER_POW
  #pragma message( "using faster pow" )
  return faster_pow(a,b);
  #else
  #ifdef FAST_POW
  #pragma message( "using fast pow" )
  return fast_pow(a,b);
  #else
  //  #pragma message( "using pow" )
  return pow(a,b);
  #endif
  #endif
}

#endif
