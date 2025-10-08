// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

/*
 * This code below is C code, obtained from the common lisp stat project under the BSD licence:
 * https://raw.githubusercontent.com/blindglobe/common-lisp-stat/3bdd28c4ae3de28dce32d8b9158c1f8d1b2e3924/lib/lowess.c
 *
 * Like much lowess code, it is derived from the initial FORTRAN code by W. S.
 * Cleveland published at NETLIB. The original FORTRAN code can be found at
 * http://www.netlib.org/go/lowess.f
 *
 * Other implementations and ports of the same code can be found at the R project: http://svn.r-project.org/R/trunk/src/library/stats/src/lowess.c
 * while a Cython version is available at: https://github.com/statsmodels/statsmodels/blob/master/statsmodels/nonparametric/_smoothers_lowess.pyx
 *
*/

#include <OpenMS/PROCESSING/SMOOTHING/FastLowessSmoothing.h>

#include <cmath>
#include <algorithm>    // std::min, std::max
#include <cstdlib>
#include <vector>
#include <functional>

namespace c_lowess
{

  /*
    The lowess code below is Translated from RATFOR lowess code of W. S.
    Cleveland as obtained from NETLIB.

    It is based on two functions written in ratfor (see below), namely lowest
    and lowess. The code has since been refactored and commented further.
  */

  /* ratfor code for lowest:
  *
  *  subroutine lowest(x,y,n,xs,ys,nleft,nright,w,userw,rw,ok)
  *  real x(n),y(n),w(n),rw(n)
  *  logical userw,ok
  *  range = x(n)-x(1)
  *  h = amax1(xs-x(nleft),x(nright)-xs)
  *  h9 = .999*h
  *  h1 = .001*h
  *  a = 0.0        # sum of weights
  *  for(j=nleft; j<=n; j=j+1){     # compute weights (pick up all ties on right)
  *         w(j)=0.
  *         r = abs(x(j)-xs)
  *         if (r<=h9) {    # small enough for non-zero weight
  *                 if (r>h1) w(j) = (1.0-(r/h)**3)**3
  *                 else      w(j) = 1.
  *                 if (userw) w(j) = rw(j)*w(j)
  *                 a = a+w(j)
  *                 }
  *         else if(x(j)>xs)break   # get out at first zero wt on right
  *         }
  *  nrt=j-1        # rightmost pt (may be greater than nright because of ties)
  *  if (a<=0.0) ok = FALSE
  *  else { # weighted least squares
  *         ok = TRUE
  *         do j = nleft,nrt
  *                 w(j) = w(j)/a   # make sum of w(j) == 1
  *         if (h>0.) {     # use linear fit
  *                 a = 0.0
  *                 do j = nleft,nrt
  *                         a = a+w(j)*x(j) # weighted center of x values
  *                 b = xs-a
  *                 c = 0.0
  *                 do j = nleft,nrt
  *                         c = c+w(j)*(x(j)-a)**2
  *                 if(sqrt(c)>.001*range) {
  *  # points are spread out enough to compute slope
  *                         b = b/c
  *                         do j = nleft,nrt
  *                                 w(j) = w(j)*(1.0+b*(x(j)-a))
  *                         }
  *                 }
  *         ys = 0.0
  *         do j = nleft,nrt
  *                 ys = ys+w(j)*y(j)
  *         }
  *  return
  *  end

  */

  /* ratfor code for lowess:
  *
  *  subroutine lowess(x,y,n,f,nsteps,delta,ys,rw,res)
  *  real x(n),y(n),ys(n),rw(n),res(n)
  *  logical ok
  *  if (n<2){ ys(1) = y(1); return }
  *  ns = max0(min0(ifix(f*float(n)),n),2)  # at least two, at most n points
  *  for(iter=1; iter<=nsteps+1; iter=iter+1){      # robustness iterations
  *         nleft = 1; nright = ns
  *         last = 0        # index of prev estimated point
  *         i = 1   # index of current point
  *         repeat{
  *                 while(nright<n){
  *  # move nleft, nright to right if radius decreases
  *                         d1 = x(i)-x(nleft)
  *                         d2 = x(nright+1)-x(i)
  *  # if d1<=d2 with x(nright+1)==x(nright), lowest fixes
  *                         if (d1<=d2) break
  *  # radius will not decrease by move right
  *                         nleft = nleft+1
  *                         nright = nright+1
  *                         }
  *                 call lowest(x,y,n,x(i),ys(i),nleft,nright,res,iter>1,rw,ok)
  *  # fitted value at x(i)
  *                 if (!ok) ys(i) = y(i)
  *  # all weights zero - copy over value (all rw==0)
  *                 if (last<i-1) { # skipped points -- interpolate
  *                         denom = x(i)-x(last)    # non-zero - proof?
  *                         for(j=last+1; j<i; j=j+1){
  *                                 alpha = (x(j)-x(last))/denom
  *                                 ys(j) = alpha*ys(i)+(1.0-alpha)*ys(last)
  *                                 }
  *                         }
  *                 last = i        # last point actually estimated
  *                 cut = x(last)+delta     # x coord of close points
  *                 for(i=last+1; i<=n; i=i+1){     # find close points
  *                         if (x(i)>cut) break     # i one beyond last pt within cut
  *                         if(x(i)==x(last)){      # exact match in x
  *                                 ys(i) = ys(last)
  *                                 last = i
  *                                 }
  *                         }
  *                 i=max0(last+1,i-1)
  *  # back 1 point so interpolation within delta, but always go forward
  *                 } until(last>=n)
  *         do i = 1,n      # residuals
  *                 res(i) = y(i)-ys(i)
  *         if (iter>nsteps) break  # compute robustness weights except last time
  *         do i = 1,n
  *                 rw(i) = abs(res(i))
  *         call sort(rw,n)
  *         m1 = 1+n/2; m2 = n-m1+1
  *         cmad = 3.0*(rw(m1)+rw(m2))      # 6 median abs resid
  *         c9 = .999*cmad; c1 = .001*cmad
  *         do i = 1,n {
  *                 r = abs(res(i))
  *                 if(r<=c1) rw(i)=1.      # near 0, avoid underflow
  *                 else if(r>c9) rw(i)=0.  # near 1, avoid underflow
  *                 else rw(i) = (1.0-(r/cmad)**2)**2
  *                 }
  *         }
  *  return
  *  end

  */

  /// Templated lowess class, call with template container (can be anything
  /// that supports random access)
  template <typename ContainerType, typename ValueType>
  class TemplatedLowess
  {

    inline ValueType pow2(ValueType x) { return x * x;  }
    inline ValueType pow3(ValueType x) { return x * x * x;  }

    /// Calculate weights for weighted regression.
    bool calculate_weights(const ContainerType& x,
                           const size_t n,
                           const ValueType current_x,
                           const bool use_resid_weights,
                           const size_t nleft,
                           const ContainerType& resid_weights,
                           ContainerType& weights,
                           size_t& nrt,
                           const ValueType h)
    {
      ValueType r;
      size_t j;

      ValueType h9 = .999 * h;
      ValueType h1 = .001 * h;
      ValueType a = 0.0; // sum of weights

      // compute weights (pick up all ties on right)
      for (j = nleft; j < n; j++)
      {

        // Compute the distance measure, then apply the tricube
        // function on the distance to get the weight.
        // use_resid_weights will be False on the first iteration, then True
        // on the subsequent ones, after some residuals have been calculated.
        weights[j] = 0.0;
        r = std::abs(x[j] - current_x);
        if (r <= h9)
        {
          if (r > h1)
          {
            // small enough for non-zero weight
            // compute tricube function: ( 1 - (r/h)^3 )^3
            weights[j] = pow3(1.0 - pow3(r / h));
          }
          else
          {
            weights[j] = 1.0;
          }

          if (use_resid_weights)
          {
            weights[j] = resid_weights[j] * weights[j];
          }

          a += weights[j];
        }
        else if (x[j] > current_x)
        {
          // get out at first zero wt on right
          break;
        }
      }

      // rightmost pt (may be greater than nright because of ties)
      nrt = j - 1;
      if (a <= 0.0)
      {
        return false;
      }
      else
      {

        // normalize weights (make sum of w[j] == 1)
        for (j = nleft; j <= nrt; j++)
        {
          weights[j] = weights[j] / a;
        }

        return true;

      }
    }

    /// Calculate smoothed/fitted y-value by weighted regression.
    void calculate_y_fit(const ContainerType& x,
                         const ContainerType& y,
                         const ValueType current_x,
                         const size_t n,
                         const size_t nleft,
                         const size_t nrt,
                         const ValueType h,
                         ValueType& ys,
                         ContainerType& weights)
    {
      ValueType range = x[n - 1] - x[0];

      if (h > 0.0)
      {
        // use linear fit

        // No regression function (e.g. lstsq) is called. Instead a "projection
        // vector" p_i_j is calculated, and y_fit[i] = sum(p_i_j * y[j]) = y_fit[i]
        // for j s.t. x[j] is in the neighborhood of x[i]. p_i_j is a function of
        // the weights, x[i], and its neighbors.
        // To save space, p_i_j is computed in place using the weight vector.

        // find weighted center of x values
        ValueType sum_weighted_x = 0.0; // originally variable a
        for (size_t j = nleft; j <= nrt; j++)
        {
          sum_weighted_x += weights[j] * x[j];
        }

        ValueType b = current_x - sum_weighted_x; // originally variable b
        ValueType weighted_sqdev = 0.0; // originally variable c
        for (size_t j = nleft; j <= nrt; j++)
        {
          weighted_sqdev += weights[j] *
                            (x[j] - sum_weighted_x) * (x[j] - sum_weighted_x);
        }

        if (sqrt(weighted_sqdev) > .001 * range)
        {
          // points are spread out enough to compute slope
          b = b / weighted_sqdev;
          for (size_t j = nleft; j <= nrt; j++)
          {
            // Compute p_i_j in place
            weights[j] = weights[j] * (1.0 + b * (x[j] - sum_weighted_x));
          }
        }
      }

      ys = 0.0;
      for (size_t j = nleft; j <= nrt; j++)
      {
        ys += weights[j] * y[j];
      }
    }

    bool lowest(const ContainerType& x,
                const ContainerType& y,
                size_t n,
                ValueType current_x, //xs
                ValueType& ys,
                size_t nleft,
                size_t nright,
                ContainerType& weights, // vector w
                bool use_resid_weights,  // userw
                const ContainerType& resid_weights)
    {
      ValueType h;
      size_t nrt; // rightmost pt (may be greater than nright because of ties)

      h = std::max(current_x - x[nleft], x[nright] - current_x);

      // Calculate the weights for the regression in this neighborhood.
      // Determine if at least some weights are positive, so a regression
      // is ok.
      bool fit_ok = calculate_weights(x, n, current_x, use_resid_weights,
                                      nleft, resid_weights,
                                      weights, nrt, h);
      if (!fit_ok)
      {
        return fit_ok;
      }

      // If it is ok to fit, run the weighted least squares regression
      calculate_y_fit(x, y, current_x, n, nleft, nrt, h, ys, weights);

      return fit_ok;
    }

    /// Find the indices bounding the k-nearest-neighbors of the current point.
    void update_neighborhood(const ContainerType& x,
                             const size_t n,
                             const size_t i,
                             size_t& nleft,
                             size_t& nright)
    {
      ValueType d1, d2;
      // A subtle loop. Start from the current neighborhood range:
      // [nleft, nright). Shift both ends rightwards by one
      // (so that the neighborhood still contains ns points), until
      // the current point is in the center (or just to the left of
      // the center) of the neighborhood. This neighborhood will
      // contain the ns-nearest neighbors of x[i].
      //
      // Once the right end hits the end of the data, hold the
      // neighborhood the same for the remaining x[i]s.
      while (nright < n - 1)
      {
        // move nleft, nright to right if radius decreases
        d1 = x[i] - x[nleft];
        d2 = x[nright + 1] - x[i];
        // if d1 <= d2 with x[nright+1] == x[nright], lowest fixes
        if (d1 <= d2) break;
        // radius will not decrease by move right
        nleft++;
        nright++;
      }
    }

    /// Update the counters of the local regression.
    void update_indices(const ContainerType& x,
                        const size_t n,
                        const ValueType delta,
                        size_t& i,
                        size_t& last,
                        ContainerType& ys)
    {
      // For most points within delta of the current point, we skip the
      // weighted linear regression (which save much computation of
      // weights and fitted points). Instead, we'll jump to the last
      // point within delta, fit the weighted regression at that point,
      // and linearly interpolate in between.

      // the last point actually estimated
      last = i;

      // This loop increments until we fall just outside of delta distance,
      // copying the results for any repeated x's along the way.
      ValueType cut = x[last] + delta;
      for (i = last + 1; i < n; i++)
      {
        // find close points
        if (x[i] > cut) break;

        // i one beyond last pt within cut
        if (x[i] == x[last])
        {
          // exact match in x
          // if tied with previous x-value, just use the already
          // fitted y, and update the last-fit counter.
          ys[i] = ys[last];
          last = i;
        }
      }


      // the next point to fit the regression at is either one prior to i (since
      // i should be the first point outside of delta) or it is "last + 1" in the
      // case that i never got incremented. This insures we always step forward.
      // -> back 1 point so interpolation within delta, but always go forward
      i = std::max(last + 1, i - 1);
    }

    /// Calculate smoothed/fitted y by linear interpolation between the current
    /// and previous y fitted by weighted regression.
    void interpolate_skipped_fits(const ContainerType& x,
                                  const size_t i,
                                  const size_t last,
                                  ContainerType& ys)
    {
      // skipped points -- interpolate
      ValueType alpha;
      ValueType denom = x[i] - x[last]; // non-zero - proof?
      for (size_t j = last + 1; j < i; j = j + 1)
      {
        alpha = (x[j] - x[last]) / denom;
        ys[j] = alpha * ys[i] + (1.0 - alpha) * ys[last];
      }
    }

    /// Calculate residual weights for the next `robustifying` iteration.
    void calculate_residual_weights(const size_t n,
                                    const ContainerType& weights,
                                    ContainerType& resid_weights)
    {
      ValueType r;

      for (size_t i = 0; i < n; i++)
      {
        resid_weights[i] = std::abs(weights[i]);
      }

      // ***********************************
      // Compute pseudo-median (take average even if we have an odd number of
      // elements), following the original implementation. We could also use a
      // true median calculation here:
      // ValueType cmad = 6.0 * median(resid_weights.begin(), resid_weights.end());
      // ***********************************

      size_t m1 = n / 2; // FORTRAN starts with one, CPP with zero
      // size_t m1 = 1 + n / 2; // original FORTRAN code
      // size_t m2 = n - m1 + 1; // see below, we don't explicitly sort but use max_element

      // Use nth element to find element m1, which produces a partially sorted
      // vector. This means we can get element m2 by looking for the maximum in the
      // remainder.
      typename ContainerType::iterator it_m1 = resid_weights.begin() + m1;
      std::nth_element(resid_weights.begin(), it_m1, resid_weights.end());
      typename ContainerType::iterator it_m2 = std::max_element(
        resid_weights.begin(), it_m1);
      ValueType cmad = 3.0 * (*it_m1 + *it_m2);
      ValueType c9 = .999 * cmad;
      ValueType c1 = .001 * cmad;

      for (size_t i = 0; i < n; i++)
      {
        r = std::abs(weights[i]);
        if (r <= c1)
        {
          // near 0, avoid underflow
          resid_weights[i] = 1.0;
        }
        else if (r > c9)
        {
          // near 1, avoid underflow
          resid_weights[i] = 0.0;
        }
        else
        {
          resid_weights[i] = pow2(1.0 - pow2(r / cmad));
        }
      }
    }

public:

    int lowess(const ContainerType& x,
               const ContainerType& y,
               double frac,    // parameter f
               int nsteps,
               ValueType delta,
               ContainerType& ys,
               ContainerType& resid_weights,   // vector rw
               ContainerType& weights   // vector res
               )
    {
      bool fit_ok;

      size_t ns, n(x.size());
      if (n < 2)
      {
        ys[0] = y[0];
        return 1;
      }

      // how many points around estimation point should be used for regression:
      // at least two, at most n points
      size_t tmp = (size_t)(frac * (double)n);
      ns = std::max(std::min(tmp, n), (size_t)2);

      // robustness iterations
      for (int iter = 1; iter <= nsteps + 1; iter++)
      {
        // start of array in C++ at 0 / in FORTRAN at 1
        // last: index of prev estimated point
        // i: index of current point
        size_t i(0), last(-1), nleft(0), nright(ns -1);

        // Fit all data points y[i] until the end of the array
        do
        {
          // Identify the neighborhood around the current x[i]
          // -> get the nearest ns points
          update_neighborhood(x, n, i, nleft, nright);

          // Calculate weights and apply fit (original lowest function)
          fit_ok = lowest(x, y, n, x[i], ys[i], nleft, nright,
                          weights, (iter > 1), resid_weights);

          // if something went wrong during the fit, use y[i] as the
          // fitted value at x[i]
          if (!fit_ok) ys[i] = y[i];

          // If we skipped some points (because of how delta was set), go back
          // and fit them by linear interpolation.
          if (last < i - 1)
          {
            interpolate_skipped_fits(x, i, last, ys);
          }

          // Update the last fit counter to indicate we've now fit this point.
          // Find the next i for which we'll run a regression.
          update_indices(x, n, delta, i, last, ys);

        } while (last < n - 1);

        // compute current residuals
        for (i = 0; i < n; i++)
        {
          weights[i] = y[i] - ys[i];
        }

        // compute robustness weights except last time
        if (iter > nsteps) break;

        calculate_residual_weights(n, weights, resid_weights);
      }
      return 0;
    }

  };
}

namespace OpenMS::FastLowessSmoothing
  {

    int lowess(const std::vector<double>& x, const std::vector<double>& y,
               double f, int nsteps,
               double delta, std::vector<double>& result)
    {
      OPENMS_PRECONDITION(delta >= 0.0, "lowess: parameter delta must be zero or larger")
      OPENMS_PRECONDITION(f > 0.0, "lowess: parameter f must be larger than 0")
      OPENMS_PRECONDITION(f <= 1.0, "lowess: parameter f must be smaller or equal to 1")
      OPENMS_PRECONDITION(nsteps >= 0, "lowess: parameter nstesp must be zero or larger")
      OPENMS_PRECONDITION(x.size() == y.size(), "Vectors x and y must have the same length")
      OPENMS_PRECONDITION(x.size() >= 2, "Need at least two points for smoothing")
      OPENMS_PRECONDITION(std::adjacent_find(x.begin(), x.end(), std::greater<double>()) == x.end(),
                          "The vector x needs to be sorted")

      size_t n = x.size();

      // result as well as working vectors need to have the correct size
      result.clear();
      result.resize(n);
      std::vector<double> resid_weights(n);
      std::vector<double> weights(n);

      c_lowess::TemplatedLowess<std::vector<double>, double> clowess;

      int retval = clowess.lowess(x, y, f, nsteps, delta, result,
                                  resid_weights, weights);

      return retval;
    }

    int lowess(const std::vector<double>& x, const std::vector<double>& y,
               std::vector<double>& result)
    {
      OPENMS_PRECONDITION(x.size() == y.size(), "Vectors x and y must have the same length")
      OPENMS_PRECONDITION(x.size() >= 2, "Need at least two points for smoothing")
      OPENMS_PRECONDITION(std::adjacent_find(x.begin(), x.end(), std::greater<double>()) == x.end(),
          "The vector x needs to be sorted")

      double delta = 0.01 * (x[ x.size()-1 ] - x[0]); // x is sorted
      return lowess(x, y, 2.0/3, 3, delta, result);
    }

  }
