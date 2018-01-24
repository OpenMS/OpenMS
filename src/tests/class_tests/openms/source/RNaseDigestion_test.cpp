// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2017.
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

#include <OpenMS/CHEMISTRY/RNaseDigestion.h>

using namespace OpenMS;
using namespace std;

///////////////////////////

START_TEST(RNaseDigestion, "$Id$")

/////////////////////////////////////////////////////////////
RNaseDigestion* rd_ptr = 0;
RNaseDigestion* rd_null = 0;

START_SECTION(([EXTRA] RNaseDigestion()))
{
  rd_ptr = new RNaseDigestion;
  TEST_NOT_EQUAL(rd_ptr, rd_null);
}
END_SECTION

START_SECTION([EXTRA] ~RNaseDigestion())
{
  delete rd_ptr;
}
END_SECTION

START_SECTION((void setEnzyme(const String& enzyme_name)))
{
  RNaseDigestion rd;
  rd.setEnzyme("RNase_T1");
  TEST_EQUAL(rd.getEnzymeName(), "RNase_T1");
  rd.setEnzyme("cusativin");
  TEST_EQUAL(rd.getEnzymeName(), "cusativin");
}
END_SECTION

START_SECTION((void digest(const NASequence& rna, set<NASequence>& output, Size min_length, Size max_length) const))
{
  RNaseDigestion rd;
  rd.setEnzyme("RNase_T1"); // cuts after G and leaves a 3'-phosphate
  set<NASequence> out;

  rd.digest(NASequence::fromString("AUC"), out);
  TEST_EQUAL(out.size(), 1);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("AUC"))==out.end(),true);
  out.clear();

  rd.digest(NASequence::fromString("AGUC"), out);
  TEST_EQUAL(out.size(), 2);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("AGp"))==out.end(),true);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("UC"))==out.end(),true);
  out.clear();

  rd.digest(NASequence::fromString("pAUGUCGCAG"), out);
  TEST_EQUAL(out.size(), 3);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("pAUGp"))==out.end(),true);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("UCGp"))==out.end(),true);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("CAG"))==out.end(),true);
  out.clear();

  rd.setMissedCleavages(2);
  rd.digest(NASequence::fromString("pAUGUCGCAG"), out);
  TEST_EQUAL(out.size(), 6);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("pAUGp"))==out.end(),true);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("pAUGUCGp"))==out.end(),true);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("pAUGUCGCAG"))==out.end(),true);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("UCGp"))==out.end(),true);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("UCGCAG"))==out.end(),true);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("CAG"))==out.end(),true);
  out.clear();

  rd.setEnzyme("cusativin");
  rd.setMissedCleavages(0);
  rd.digest(NASequence::fromString("CCCAUCCG"), out);
  TEST_EQUAL(out.size(), 3);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("CCCp"))==out.end(),true);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("AUCCp"))==out.end(),true);
  TEST_NOT_EQUAL(out.find(NASequence::fromString("G"))==out.end(),true);
}
END_SECTION


/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
