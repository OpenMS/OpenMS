// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/MATH/STATISTICS/CumulativeBinomial.h>
#include <numeric>

using namespace OpenMS;
using namespace OpenMS::Math;

/////////////////////////////////////////////////////////////

START_TEST( CumulativeBinomial, "$Id$" );

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

START_SECTION(static double compute(Size n, Size k, double p))
{
  // Extreme cases
  TEST_EQUAL(CumulativeBinomial::compute(10, 5, 0.0), 0.0)
  TEST_EQUAL(CumulativeBinomial::compute(10, 5, 1.0), 1.0)
  TEST_EQUAL(CumulativeBinomial::compute(10, 15, 0.5), 1.0)

  TEST_EQUAL(CumulativeBinomial::compute(10, 5, 1e-10), nexttoward(1.0, 0.0))
   TEST_REAL_SIMILAR(CumulativeBinomial::compute(500, 5, 1e-10), 1.0)
  TEST_EQUAL(CumulativeBinomial::compute(500, 5, 0.9999999), 0.0)

  // Normal cases
  TEST_REAL_SIMILAR(CumulativeBinomial::compute(10, 5, 0.5), 0.37695)
  TEST_REAL_SIMILAR(CumulativeBinomial::compute(500, 5, 0.1), 5.58143e-18)
  TEST_REAL_SIMILAR(CumulativeBinomial::compute(500, 5, 0.9), 1e-05)
  TEST_REAL_SIMILAR(CumulativeBinomial::compute(500, 5, 0.5), 7.92409e-142)
  TEST_REAL_SIMILAR(CumulativeBinomial::compute(500, 50, 0.1), 0.478198)
  TEST_REAL_SIMILAR(CumulativeBinomial::compute(500, 50, 0.9), 1e-05)
  TEST_REAL_SIMILAR(CumulativeBinomial::compute(500, 50, 0.5), 8.78876e-83)
  TEST_REAL_SIMILAR(CumulativeBinomial::compute(500, 499, 1e-10), 0.99999)
  TEST_REAL_SIMILAR(CumulativeBinomial::compute(500, 499, 0.1), 1.0)
  TEST_REAL_SIMILAR(CumulativeBinomial::compute(500, 499, 0.9), 1.0)
}
END_SECTION

END_TEST
