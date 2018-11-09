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

#include <OpenMS/FORMAT/XQuestResultXMLFile.h>

#include <QStringList>

using namespace OpenMS;

START_TEST(XQuestResultXMLFile, "$Id$")


START_SECTION(void store(const String& filename, const std::vector<ProteinIdentification>& poid, const std::vector<PeptideIdentification>& peid) const)

  std::vector<ProteinIdentification> protein_ids;
  std::vector<PeptideIdentification> peptide_ids;

  String xquest_input_file= OPENMS_GET_TEST_DATA_PATH("XQuestResultXMLFile_test_data.xquest.xml");
  XQuestResultXMLFile().load(xquest_input_file, peptide_ids, protein_ids);

  String out_file;
  NEW_TMP_FILE(out_file)

  XQuestResultXMLFile().store(out_file, protein_ids, peptide_ids);

  std::vector< PeptideIdentification > peptide_id_vector;
  std::vector< ProteinIdentification > protein_id_vector;
  XQuestResultXMLFile().load(out_file, peptide_id_vector, protein_id_vector);

  for (Size i = 0; i < peptide_id_vector.size(); ++i)
  {
    std::vector<PeptideHit> hits = peptide_id_vector[i].getHits();
    for (Size k = 0; k < hits.size(); ++k)
    {
      TEST_REAL_SIMILAR(hits[k].getScore(), peptide_ids[i].getHits()[k].getScore())
      TEST_EQUAL(hits[k].getCharge(), peptide_ids[i].getHits()[k].getCharge())

      // only in alpha peptide hit
      if (k == 0)
      {
        TEST_EQUAL(hits[k].getMetaValue("xl_rank"), peptide_ids[i].getHits()[k].getMetaValue("xl_rank"))
        TEST_EQUAL(hits[k].getMetaValue("xl_type"), peptide_ids[i].getHits()[k].getMetaValue("xl_type"))
      }
      TEST_EQUAL(hits[k].getMetaValue("xl_pos"), peptide_ids[i].getHits()[k].getMetaValue("xl_pos"))
      TEST_EQUAL(hits[k].getMetaValue("xl_chain"), peptide_ids[i].getHits()[k].getMetaValue("xl_chain"))
      TEST_REAL_SIMILAR(hits[k].getMetaValue("OpenXQuest:match-odds"), peptide_ids[i].getHits()[k].getMetaValue("OpenXQuest:match-odds"))
      TEST_REAL_SIMILAR(hits[k].getMetaValue("OpenXQuest:intsum"), peptide_ids[i].getHits()[k].getMetaValue("OpenXQuest:intsum"))
    }
  }

END_SECTION

END_TEST
