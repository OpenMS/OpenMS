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

START_SECTION((void digest(const NASequence& rna, vector<NASequence>& output, Size min_length, Size max_length) const))
{
  RNaseDigestion rd;
  rd.setEnzyme("RNase_T1"); // cuts after G and leaves a 3'-phosphate
  vector<NASequence> out;

  rd.digest(NASequence::fromString("AUC"), out);
  TEST_EQUAL(out.size(), 1);
  TEST_STRING_EQUAL(out[0].toString(), "AUC");
  out.clear();

  rd.digest(NASequence::fromString("AGUC"), out);
  TEST_EQUAL(out.size(), 2);
  TEST_STRING_EQUAL(out[0].toString(), "AGp");
  TEST_STRING_EQUAL(out[1].toString(), "UC");
  out.clear();

  rd.digest(NASequence::fromString("pAUGUCGCAG"), out);
  TEST_EQUAL(out.size(), 3);
  TEST_STRING_EQUAL(out[0].toString(), "pAUGp");
  TEST_STRING_EQUAL(out[1].toString(), "UCGp");
  TEST_STRING_EQUAL(out[2].toString(), "CAG");
  out.clear();

  // RNase T1 should cut after G and m1G, but not after Gm:
  rd.digest(NASequence::fromString("G[m1G][Gm]A"), out);
  TEST_EQUAL(out.size(), 3);
  TEST_STRING_EQUAL(out[0].toString(), "Gp");
  TEST_STRING_EQUAL(out[1].toString(), "[m1G]p");
  TEST_STRING_EQUAL(out[2].toString(), "[Gm]A");
  out.clear();

  rd.setMissedCleavages(2);
  rd.digest(NASequence::fromString("pAUGUCGCAG"), out);
  TEST_EQUAL(out.size(), 6);
  TEST_STRING_EQUAL(out[0].toString(), "pAUGp");
  TEST_STRING_EQUAL(out[1].toString(), "pAUGUCGp");
  TEST_STRING_EQUAL(out[2].toString(), "pAUGUCGCAG");
  TEST_STRING_EQUAL(out[3].toString(), "UCGp");
  TEST_STRING_EQUAL(out[4].toString(), "UCGCAG");
  TEST_STRING_EQUAL(out[5].toString(), "CAG");
  out.clear();

  rd.setEnzyme("cusativin");
  rd.setMissedCleavages(0);
  rd.digest(NASequence::fromString("CCCAUCCG"), out);
  TEST_EQUAL(out.size(), 3);
  TEST_STRING_EQUAL(out[0].toString(), "CCCp");
  TEST_STRING_EQUAL(out[1].toString(), "AUCCp");
  TEST_STRING_EQUAL(out[2].toString(), "G");

  rd.setEnzyme("no cleavage");
  rd.setMissedCleavages(3);
  rd.digest(NASequence::fromString("CCCAUCCG"), out);
  TEST_EQUAL(out.size(), 1);
  TEST_STRING_EQUAL(out[0].toString(), "CCCAUCCG");

  rd.setEnzyme("unspecific cleavage");
  rd.setMissedCleavages(0); // shouldn't matter for the result
  rd.digest(NASequence::fromString("ACGU"), out);
  TEST_EQUAL(out.size(), 10);
  TEST_STRING_EQUAL(out[0].toString(), "A");
  TEST_STRING_EQUAL(out[1].toString(), "AC");
  TEST_STRING_EQUAL(out[2].toString(), "ACG");
  TEST_STRING_EQUAL(out[3].toString(), "ACGU");
  TEST_STRING_EQUAL(out[4].toString(), "C");
  TEST_STRING_EQUAL(out[5].toString(), "CG");
  TEST_STRING_EQUAL(out[6].toString(), "CGU");
  TEST_STRING_EQUAL(out[7].toString(), "G");
  TEST_STRING_EQUAL(out[8].toString(), "GU");
  TEST_STRING_EQUAL(out[9].toString(), "U");
}
END_SECTION

START_SECTION((void digest(IdentificationData& id_data, Size min_length = 0,
                Size max_length = 0) const))
{
  IdentificationData id_data;
  IdentificationData::ParentMolecule rna("test", IdentificationData::MoleculeType::RNA, "pAUGUCGCAG");
  id_data.registerParentMolecule(rna);

  RNaseDigestion rd;
  rd.setEnzyme("RNase_T1"); // cuts after G and leaves a 3'-phosphate
  rd.digest(id_data);

  TEST_EQUAL(id_data.getIdentifiedOligos().size(), 3);

  /// multiple occurrences of the same oligo:
  IdentificationData id_data2;
  rna.sequence = "ACUGACUGG";
  id_data2.registerParentMolecule(rna);

  rd.digest(id_data2, 2);

  TEST_EQUAL(id_data2.getIdentifiedOligos().size(), 1);
  ABORT_IF(id_data2.getIdentifiedOligos().empty());
  IdentificationData::IdentifiedOligoRef ref = id_data2.getIdentifiedOligos().begin();
  TEST_EQUAL(ref->parent_matches.size(), 1);
  ABORT_IF(ref->parent_matches.empty());
  // oligo sequence matches in two locations:
  const set<IdentificationData::MoleculeParentMatch>& matches =
    ref->parent_matches.begin()->second;
  TEST_EQUAL(matches.size(), 2);
  ABORT_IF(matches.size() < 2);
  auto match_it = matches.begin();
  TEST_EQUAL(match_it->start_pos, 0);
  ++match_it;
  TEST_EQUAL(match_it->start_pos, 4);
}
END_SECTION

/////////////////////////////////////////////////////////////
/////////////////////////////////////////////////////////////
END_TEST
