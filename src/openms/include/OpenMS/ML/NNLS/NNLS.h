// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Chris Bielow $
// $Authors: Chris Bielow $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/config.h>

namespace OpenMS
{

  namespace NNLS
  {
    typedef int integer;

    /*     SUBROUTINE NNLS  (A,MDA,M,N,B,X,RNORM,W,ZZ,INDEX,MODE) */

    /*  Algorithm NNLS: NONNEGATIVE LEAST SQUARES */

    /*  The original version of this code was developed by */
    /*  Charles L. Lawson and Richard J. Hanson at Jet Propulsion Laboratory */
    /*  1973 JUN 15, and published in the book */
    /*  "SOLVING LEAST SQUARES PROBLEMS", Prentice-Hall, 1974. */
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
    int OPENMS_DLLAPI nnls_(double * a, integer * mda, integer * m, integer *
                            n, double * b, double * x, double * rnorm, double * w,
                            double * zz, integer * index, integer * mode);

    /* Subroutine */
    int OPENMS_DLLAPI g1_(double *, double *, double *, double *, double *);

    /* Subroutine */
    int OPENMS_DLLAPI h12_(integer *, integer *, integer *, integer *, double *, integer *, double *, double *, integer *, integer *, integer *);

    /* Subroutine */
    double OPENMS_DLLAPI diff_(double *, double *);

    /* Subroutine */
    double OPENMS_DLLAPI d_sign_(double & a, double & b);
  }

}

