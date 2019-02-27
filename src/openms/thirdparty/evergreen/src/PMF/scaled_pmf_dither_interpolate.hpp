#ifndef _SCALED_PMF_DITHER_INTERPOLATE_HPP
#define _SCALED_PMF_DITHER_INTERPOLATE_HPP

#include "scaled_pmf_dither.hpp"

inline PMF scaled_pmf_dither_interpolate(const PMF & pmf, const Vector<double> & factor, double sigma_squared) {
  // TODO: implement more general form that simultaneously dithers and
  // interpolates. If fabs of all scaling factors are <= 1, then
  // interpolation is unnecessary:
  if ( factor <= 1.0 && factor >= -1.0 )
    return scaled_pmf_dither(pmf, factor, sigma_squared);
  return scaled_pmf_interpolate(pmf, factor);
}

#endif
