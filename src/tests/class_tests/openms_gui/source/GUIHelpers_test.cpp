// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

#include <OpenMS/CONCEPT/ClassTest.h>

///////////////////////////

#include <OpenMS/VISUAL/MISC/GUIHelpers.h>

///////////////////////////

using namespace OpenMS;
using namespace std;

START_TEST(GUIHelpers, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////


START_SECTION(([EXTRA] size_t OverlapDetector::placeItem(double x_start, double x_end)))
	GUIHelpers::OverlapDetector od(3);
	TEST_EQUAL(od.placeItem(1, 3), 0);
  TEST_EQUAL(od.placeItem(1, 3), 1);
  TEST_EQUAL(od.placeItem(1, 3), 2);
  TEST_EQUAL(od.placeItem(1, 3), 0);

  TEST_EQUAL(od.placeItem(4, 8), 0);
  TEST_EQUAL(od.placeItem(5, 11), 1);
  TEST_EQUAL(od.placeItem(9, 10), 0);

  TEST_EQUAL(od.placeItem(12, 20), 0);
  TEST_EQUAL(od.placeItem(12, 18), 1);
  TEST_EQUAL(od.placeItem(12, 19), 2);

  TEST_EQUAL(od.placeItem(16, 25), 1);
  TEST_EQUAL(od.placeItem(16, 25), 2);
  TEST_EQUAL(od.placeItem(16, 25), 0);
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST



