// -*- mode: c++; c-basic-offset: 4; -*-
/************************************************************************
 * Copyright 2009 University Corporation for Atmospheric Research.
 * All rights reserved.
 *
 * Use of this code is subject to the standard BSD license:
 *
 *  http://www.opensource.org/licenses/bsd-license.html
 *
 * See the COPYRIGHT file in the source distribution for the license text,
 * or see this web page:
 *
 *  http://www.eol.ucar.edu/homes/granger/bspline/doc/
 *
 *************************************************************************/

/**
 * @file
 *
 * This file defines the implementation for the BSpline and BSplineBase
 * templates.
 **/
#include "BSpline.h"
#include "BandedMatrix.h"

#include <vector>
#include <algorithm>
#include <iterator>
#include <iostream>
#include <iomanip>
#include <map>
#include <assert.h>


//////////////////////////////////////////////////////////////////////
//////////////////////////////////////////////////////////////////////
// BSpline Class
//////////////////////////////////////////////////////////////////////

template<class T> struct BSplineP {
        std::vector<T> spline;
        std::vector<T> A;
};

//////////////////////////////////////////////////////////////////////
// Construction/Destruction
//////////////////////////////////////////////////////////////////////

/*
 * This BSpline constructor constructs and sets up a new base and 
 * solves for the spline curve coeffiecients all at once.
 */
template<class T> BSpline<T>::BSpline(const T *x,
                                      int nx,
                                      const T *y,
                                      double wl,
                                      int bc_type,
                                      int num_nodes) :
    BSplineBase<T>(x, nx, wl, bc_type, num_nodes), s(new BSplineP<T>) {
    solve(y);
}
//////////////////////////////////////////////////////////////////////
/*
 * Create a new spline given a BSplineBase.
 */
template<class T> BSpline<T>::BSpline(BSplineBase<T> &bb,
                                      const T *y) :
    BSplineBase<T>(bb), s(new BSplineP<T>) {
    solve(y);
}
//////////////////////////////////////////////////////////////////////
/*
 * (Re)calculate the spline for the given set of y values.
 */
template<class T> bool BSpline<T>::solve(const T *y) {
    if (!OK)
        return false;

    // Any previously calculated curve is now invalid.
    s->spline.clear();
    OK = false;

    // Given an array of data points over x and its precalculated
    // P+Q matrix, calculate the b vector and solve for the coefficients.
    std::vector<T> &B = s->A;
    std::vector<T> &A = s->A;
    A.clear();
    A.resize(M+1);

    if (Debug())
        std::cerr << "Solving for B..." << std::endl;

    // Find the mean of these data
    mean = 0.0;
    int i;
    for (i = 0; i < NX; ++i) {
        mean += y[i];
    }
    mean = mean / (double)NX;
    if (Debug())
        std::cerr << "Mean for y: " << mean << std::endl;

    int mx, m, j;
    for (j = 0; j < NX; ++j) {
        // Which node does this put us in?
        T &xj = base->X[j];
        T yj = y[j] - mean;
        mx = (int)((xj - xmin) / DX);

        for (m = my::max(0, mx-1); m <= my::min(mx+2, M); ++m) {
            B[m] += yj * this->Basis(m, xj);
        }
    }

    if (Debug() && M < 30) {
        std::cerr << "Solution a for (P+Q)a = b" << std::endl;
        std::cerr << " b: " << B << std::endl;
    }

    // Now solve for the A vector in place.
    if (LU_solve_banded(base->Q, A, 3) != 0) {
        if (Debug())
            std::cerr << "LU_solve_banded() failed." << std::endl;
    } else {
        OK = true;
        if (Debug())
            std::cerr << "Done." << std::endl;
        if (Debug() && M < 30) {
            std::cerr << " a: " << A << std::endl;
            std::cerr << "LU factor of (P+Q) = " << std::endl << base->Q
                    << std::endl;
        }
    }
    return (OK);
}
//////////////////////////////////////////////////////////////////////
template<class T> BSpline<T>::~BSpline() {
    delete s;
}
//////////////////////////////////////////////////////////////////////
template<class T> T BSpline<T>::coefficient(int n) {
    if (OK)
        if (0 <= n && n <= M)
            return s->A[n];
    return 0;
}
//////////////////////////////////////////////////////////////////////
template<class T> T BSpline<T>::evaluate(T x) {
    T y = 0;
    if (OK) {
        int n = (int)((x - xmin)/DX);
        for (int i = my::max(0, n-1); i <= my::min(M, n+2); ++i) {
            y += s->A[i] * this->Basis(i, x);
        }
        y += mean;
    }
    return y;
}
//////////////////////////////////////////////////////////////////////
template<class T> T BSpline<T>::slope(T x) {
    T dy = 0;
    if (OK) {
        int n = (int)((x - xmin)/DX);
        for (int i = my::max(0, n-1); i <= my::min(M, n+2); ++i) {
            dy += s->A[i] * this->DBasis(i, x);
        }
    }
    return dy;
}
