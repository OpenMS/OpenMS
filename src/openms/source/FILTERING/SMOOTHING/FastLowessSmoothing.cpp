// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2015.
//
// This software is released under a three-clause BSD license:
//  * Redistributions of source code must retain the above copyright
//    notice, this list of conditions and the following disclaimer.
//  * Redistributions in binary form must reproduce the above copyright
//    notice, this list of conditions and the following disclaimer in the
//    documentation and/or other materials provided with the distribution.
//  * Neither the name of any author or any participating institution
//    may be used to endorse or promote products derived from this software
//    without specific prior written permission.
// For a full list of authors, refer to the file AUTHORS.
// --------------------------------------------------------------------------
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
// AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
// IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
// ARE DISCLAIMED. IN NO EVENT SHALL ANY OF THE AUTHORS OR THE CONTRIBUTING
// INSTITUTIONS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
// OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY,
// WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR
// OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF
// ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// --------------------------------------------------------------------------
// $Maintainer: Hannes Roest $
// $Authors: Hannes Roest $
// --------------------------------------------------------------------------

// This code below is C code, obtained from the common lisp stat project under the BSD licence: 
// https://raw.githubusercontent.com/blindglobe/common-lisp-stat/3bdd28c4ae3de28dce32d8b9158c1f8d1b2e3924/lib/lowess.c

namespace c_lowess
{

/*
  Translated from RATFOR lowess code of W. S. Cleveland as obtained from NETLIB
*/

#include "math.h"
#include <stdlib.h>

#define FALSE 0
#define TRUE 1

extern long max(long, long);
extern long min(long, long);

static double pow2(double x) { return(x * x); }
static double pow3(double x) { return(x * x * x); }
static double local_fmax(double x, double y) { return (x > y ? x : y); }

int 
static compar(const void *aa, const void *bb)
{
  const double* a = (double*)aa;
  const double* b = (double*)bb;

  if (*a < *b) return(-1);
  else if (*a > *b) return(1);
  else return(0);
}

static void
lowest(double *x, double *y, size_t n, double xs, double *ys, long nleft, long nright,
       double *w, int userw, double *rw, int *ok)
{
  double range, h, h1, h9, a, b, c, r;
  long j, nrt;

  range = x[n - 1] - x[0];
  h = local_fmax(xs - x[nleft], x[nright] - xs);
  h9 = .999 * h;
  h1 = .001 * h;

  a = 0.0; // sum of weights

  // compute weights (pick up all ties on right)
  for (j = nleft; j < n; j++) 
  {

    w[j]=0.0;
    r = fabs(x[j] - xs);
    if (r <= h9)
    {  
      if (r > h1) 
      {
        // small enough for non-zero weight
        w[j] = pow3(1.0-pow3(r/h));
      }
      else w[j] = 1.0;

      if (userw) w[j] = rw[j] * w[j];
      a += w[j];
    }
    else if (x[j] > xs) 
    {
      // get out at first zero wt on right
      break;
    }
  }

  // rightmost pt (may be greater than nright because of ties)
  nrt = j - 1;
  if (a <= 0.0) 
  {
    *ok = FALSE;
  }
  else 
  { 
    *ok = TRUE;

    // weighted least squares
    
    // make sum of w[j] == 1
    for (j = nleft; j <= nrt; j++) w[j] = w[j] / a;

    if (h > 0.0) 
    { 
      // use linear fit

      // find weighted center of x values
      for (j = nleft, a = 0.0; j <= nrt; j++) a += w[j] * x[j];

      b = xs - a;
      for (j = nleft, c = 0.0; j <= nrt; j++) 
      {
        c += w[j] * (x[j] - a) * (x[j] - a);
      }

      if(sqrt(c) > .001 * range) 
      {
        // points are spread out enough to compute slope
        b = b/c;
        for (j = nleft; j <= nrt; j++) 
        {
          w[j] = w[j] * (1.0 + b*(x[j] - a));
        }
      }
    }

    for (j = nleft, *ys = 0.0; j <= nrt; j++)
    {
      *ys += w[j] * y[j];
    }

  }
}

static void
sort(double *x, size_t n)
{
  extern void qsort();

  qsort(x, n, sizeof(double), compar);
}



int
lowess(double *x, double *y, size_t n,
       double f, size_t nsteps,
       double delta, double *ys, double *rw, double *res)
{
  int iter, ok;
  long i, j, last, m1, m2, nleft, nright, ns;
  double d1, d2, denom, alpha, cut, cmad, c9, c1, r;
  
  if (n < 2) 
  { 
    ys[0] = y[0]; 
    return(1);
  }

  // at least two, at most n points
  ns = max(min( (long) (f * n), n), 2); 

  // robustness iterations
  for (iter = 1; iter <= nsteps + 1; iter++)
  {

    // start of array in C++ at 0 / in FORTRAN at 1
    nleft = 0;
    nright = ns - 1;
    last = -1;        // index of prev estimated point
    i = 0;            // index of current point

    do 
    {

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

      lowest(x, y,
             n, x[i],
             &ys[i],
             nleft, nright,
             res, (iter > 1), rw, &ok);

      // fitted value at x[i]
      if (! ok) ys[i] = y[i];

      // all weights zero - copy over value (all rw==0)
      if (last < i - 1) 
      {
        // skipped points -- interpolate
        // non-zero - proof?
        denom = x[i] - x[last];
        for (j = last + 1; j < i; j = j + 1)
        {
          alpha = (x[j] - x[last]) / denom;
          ys[j] = alpha * ys[i] + (1.0 - alpha) * ys[last];
        }
      }
        
      // last point actually estimated
      last = i;

      // x coord of close points
      cut = x[last] + delta;
      for(i=last + 1; i < n; i++) 
      { 
        // find close points 
        if (x[i] > cut) break;

        // i one beyond last pt within cut
        if (x[i] == x[last]) 
        {
          // exact match in x
          ys[i] = ys[last];
          last = i;
        }
      }

      // back 1 point so interpolation within delta, but always go forward
      i = max(last + 1,i - 1);

    } while (last < n - 1);

    // residuals
    for (i = 0; i < n; i++)
    {
      res[i] = y[i] - ys[i];
    }

    // compute robustness weights except last time
    if (iter > nsteps) break;

    for (i = 0; i < n; i++) 
    {
      rw[i] = fabs(res[i]);
    }

    sort(rw,n);
    m1 = 1 + n / 2;
    m2 = n - m1 + 1;
    // cmad = 6 * median abs resid
    cmad = 3.0 * (rw[m1] + rw[m2]);
    c9 = .999 * cmad; 
    c1 = .001 * cmad;

    for (i = 0; i < n; i++) 
    {
      r = fabs(res[i]);
      if (r <= c1) 
      {
        // near 0, avoid underflow
        rw[i] = 1.0;
      }
      else if (r > c9) 
      {
        // near 1, avoid underflow
        rw[i] = 0.0;
      }
      else 
      {
        rw[i] = pow2(1.0 - pow2(r / cmad));
      }
    }

  }
  return(0);
}


}
 

