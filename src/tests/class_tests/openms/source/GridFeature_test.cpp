// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2014.
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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/DATASTRUCTURES/GridFeature.h>
#include <OpenMS/KERNEL/BaseFeature.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/CHEMISTRY/AASequence.h>

///////////////////////////

START_TEST(GridFeature, "$Id$")

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////

using namespace OpenMS;
using namespace std;

GridFeature* gf_ptr = 0;
GridFeature* gf_nullPointer = 0;

START_SECTION((GridFeature(const BaseFeature& feature, Size map_index, Size feature_index)))
{
  BaseFeature bf;
  gf_ptr = new GridFeature(bf, 0, 0);
  TEST_NOT_EQUAL(gf_ptr, gf_nullPointer);
}
END_SECTION

START_SECTION((~GridFeature()))
{
  delete gf_ptr;
}
END_SECTION

START_SECTION((const BaseFeature& getFeature() const))
{
  BaseFeature bf;
  bf.setRT(1.1);
  bf.setMZ(2.2);
  bf.setCharge(3);
  const BaseFeature bf_const(bf);
  GridFeature gf(bf_const, 0, 0);
  TEST_EQUAL(gf.getFeature() == bf_const, true);
}
END_SECTION

START_SECTION((Size getMapIndex() const))
{
  BaseFeature bf;
  GridFeature gf(bf, 123, 0);
  TEST_EQUAL(gf.getMapIndex(), 123);
}
END_SECTION

START_SECTION((Size getFeatureIndex() const))
{
  BaseFeature bf;
  GridFeature gf(bf, 0, 123);
  TEST_EQUAL(gf.getFeatureIndex(), 123);
}
END_SECTION

START_SECTION((Int getID() const))
{
  BaseFeature bf;
  GridFeature gf(bf, 0, 123);
  TEST_EQUAL(gf.getID(), 123);
}
END_SECTION

START_SECTION((const std::set<AASequence>& getAnnotations() const))
{
  BaseFeature bf;
  GridFeature gf(bf, 0, 0);
  TEST_EQUAL(gf.getAnnotations().size(), 0);
  bf.getPeptideIdentifications().resize(2);
  PeptideHit hit;
  hit.setSequence(AASequence::fromString("AAA"));
  bf.getPeptideIdentifications()[0].insertHit(hit);
  hit.setSequence(AASequence::fromString("CCC"));
  bf.getPeptideIdentifications()[1].insertHit(hit);
  GridFeature gf2(bf, 0, 0);
  TEST_EQUAL(gf2.getAnnotations().size(), 2);
  TEST_EQUAL(*(gf2.getAnnotations().begin()), AASequence::fromString("AAA"));
  TEST_EQUAL(*(gf2.getAnnotations().rbegin()), AASequence::fromString("CCC"));
}
END_SECTION

START_SECTION((double getRT() const))
{
  BaseFeature bf;
  bf.setRT(4.56);
  GridFeature gf(bf, 0, 123);
  TEST_REAL_SIMILAR(gf.getRT(), 4.56);
}
END_SECTION

START_SECTION((double getMZ() const))
{
  BaseFeature bf;
  bf.setMZ(4.56);
  GridFeature gf(bf, 0, 123);
  TEST_REAL_SIMILAR(gf.getMZ(), 4.56);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
