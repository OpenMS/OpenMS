// -*- mode: C++; tab-width: 2; -*-
// vi: set ts=2:
//
// --------------------------------------------------------------------------
//                   OpenMS Mass Spectrometry Framework 
// --------------------------------------------------------------------------
//  Copyright (C) 2003-2011 -- Oliver Kohlbacher, Knut Reinert
//
//  This library is free software; you can redistribute it and/or
//  modify it under the terms of the GNU Lesser General Public
//  License as published by the Free Software Foundation; either
//  version 2.1 of the License, or (at your option) any later version.
//
//  This library is distributed in the hope that it will be useful,
//  but WITHOUT ANY WARRANTY; without even the implied warranty of
//  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
//  Lesser General Public License for more details.
//
//  You should have received a copy of the GNU Lesser General Public
//  License along with this library; if not, write to the Free Software
//  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#ifndef OPENMS_MATH_MISC_NNLS_NNLS_H
#define OPENMS_MATH_MISC_NNLS_NNLS_H

#include <OpenMS/config.h>

namespace OpenMS 
{

  namespace NNLS
  {
    typedef int integer;
    typedef double doublereal;

    /*     SUBROUTINE NNLS  (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE) */

    /*  Algorithm NNLS: NONNEGATIVE LEAST SQUARES */

    /*  The original version of this code was developed by */
    /*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory */
    /*  1973 JUN 15, and published in the book */
    /*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-HalL, 1974. */
    /*  Revised FEB 1995 to accompany reprinting of the book by SIAM. */

    /*     GIVEN AN M BY N MATRIX, A, AND AN M-VECTOR, B,  COMPUTE AN */
    /*     N-VECTOR, X, THAT SOLVES THE LEAST SQUARES PROBLEM */

    /*                      A * X = B  SUBJECT TO X .GE. 0 */
    /*     ------------------------------------------------------------------ */
    /*                     Subroutine Arguments */

    /*     A(),MDA,M,N     MDA IS THE FIRST DIMENSIONING PARAMETER FOR THE */
    /*                     ARRAY, A().   ON ENTRY A() CONTAINS THE M BY N */
    /*                     MATRIX, A.           ON EXIT A() CONTAINS */
    /*                     THE PRODUCT MATRIX, Q*A , WHERE Q IS AN */
    /*                     M BY M ORTHOGONAL MATRIX GENERATED IMPLICITLY BY */
    /*                     THIS SUBROUTINE. */
    /*     B()     ON ENTRY B() CONTAINS THE M-VECTOR, B.   ON EXIT B() CON- */
    /*             TAINS Q*B. */
    /*     X()     ON ENTRY X() NEED NOT BE INITIALIZED.  ON EXIT X() WILL */
    /*             CONTAIN THE SOLUTION VECTOR. */
    /*     RNORM   ON EXIT RNORM CONTAINS THE EUCLIDEAN NORM OF THE */
    /*             RESIDUAL VECTOR. */
    /*     W()     AN N-ARRAY OF WORKING SPACE.  ON EXIT W() WILL CONTAIN */
    /*             THE DUAL SOLUTION VECTOR.   W WILL SATISFY W(I) = 0. */
    /*             FOR ALL I IN SET P  AND W(I) .LE. 0. FOR ALL I IN SET Z */
    /*     ZZ()     AN M-ARRAY OF WORKING SPACE. */
    /*     INDEX()     AN INTEGER WORKING ARRAY OF LENGTH AT LEAST N. */
    /*                 ON EXIT THE CONTENTS OF THIS ARRAY DEFINE THE SETS */
    /*                 P AND Z AS FOLLOWS.. */

    /*                 INDEX(1)   THRU INDEX(NSETP) = SET P. */
    /*                 INDEX(IZ1) THRU INDEX(IZ2)   = SET Z. */
    /*                 IZ1 = NSETP + 1 = NPP1 */
    /*                 IZ2 = N */
    /*     MODE    THIS IS A SUCCESS-FAILURE FLAG WITH THE FOLLOWING */
    /*             MEANINGS. */
    /*             1     THE SOLUTION HAS BEEN COMPUTED SUCCESSFULLY. */
    /*             2     THE DIMENSIONS OF THE PROBLEM ARE BAD. */
    /*                   EITHER M .LE. 0 OR N .LE. 0. */
    /*             3    ITERATION COUNT EXCEEDED.  MORE THAN 3*N ITERATIONS. */
    int OPENMS_DLLAPI nnls_(doublereal *a, integer *mda, integer *m, integer *
      n, doublereal *b, doublereal *x, doublereal *rnorm, doublereal *w, 
      doublereal *zz, integer *index, integer *mode);

    /* Subroutine */
    int OPENMS_DLLAPI g1_(doublereal *, doublereal *, doublereal *, doublereal *, doublereal *);

    /* Subroutine */
    int OPENMS_DLLAPI h12_(integer *, integer *, integer *, integer *, doublereal *, integer *, doublereal *, doublereal *, integer *, integer *, integer *);

    /* Subroutine */
    doublereal OPENMS_DLLAPI diff_(doublereal *, doublereal *);
    
    /* Subroutine */
    double OPENMS_DLLAPI d_sign_(double& a, double& b);
  }

}

#endif //OPENMS_MATH_MISC_NNLS_NNLS_H
