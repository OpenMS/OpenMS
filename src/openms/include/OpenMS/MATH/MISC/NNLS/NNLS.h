// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2020.
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

