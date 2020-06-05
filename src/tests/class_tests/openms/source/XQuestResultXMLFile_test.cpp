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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/CONCEPT/Constants.h>

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

  for (Size i = 0; i < peptide_id_vector.size(); i+=20)
  {
    std::vector<PeptideHit> hits = peptide_id_vector[i].getHits();
    for (Size k = 0; k < hits.size(); ++k)
    {
      TEST_REAL_SIMILAR(hits[k].getScore(), peptide_ids[i].getHits()[k].getScore())
      TEST_EQUAL(hits[k].getCharge(), peptide_ids[i].getHits()[k].getCharge())

      TEST_EQUAL(hits[k].getMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK), peptide_ids[i].getHits()[k].getMetaValue(Constants::UserParam::OPENPEPXL_XL_RANK))
      TEST_EQUAL(hits[k].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), peptide_ids[i].getHits()[k].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE))
      TEST_EQUAL(hits[k].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), peptide_ids[i].getHits()[k].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1))
      TEST_EQUAL(hits[k].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), peptide_ids[i].getHits()[k].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2))
      TEST_EQUAL(hits[k].getSequence().toString(), peptide_ids[i].getHits()[k].getSequence().toString())
      TEST_EQUAL(hits[k].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), peptide_ids[i].getHits()[k].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE))
      TEST_REAL_SIMILAR(hits[k].getMetaValue("OpenPepXL:match-odds"), peptide_ids[i].getHits()[k].getMetaValue("OpenPepXL:match-odds"))
      TEST_REAL_SIMILAR(hits[k].getMetaValue("OpenPepXL:intsum"), peptide_ids[i].getHits()[k].getMetaValue("OpenPepXL:intsum"))
    }
  }

  TEST_EQUAL(peptide_id_vector.size(), 296)
  TEST_EQUAL(peptide_id_vector[0].getHits().size(), 1)
  TEST_EQUAL(peptide_id_vector[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "cross-link")
  TEST_EQUAL(peptide_id_vector[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 14)
  TEST_EQUAL(peptide_id_vector[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), 5)
  TEST_EQUAL(peptide_id_vector[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "ANYWHERE")
  TEST_EQUAL(peptide_id_vector[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA), "ANYWHERE")
  TEST_EQUAL(peptide_id_vector[0].getHits()[0].getSequence().toString(), "LTEIISHDPNIELHKK")
  TEST_EQUAL(peptide_id_vector[0].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "VEGCPKHPK")

  TEST_EQUAL(peptide_id_vector[17].getHits().size(), 1)
  TEST_EQUAL(peptide_id_vector[17].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "cross-link")
  TEST_EQUAL(peptide_id_vector[17].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 15)
  TEST_EQUAL(peptide_id_vector[17].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), 11)
  TEST_EQUAL(peptide_id_vector[17].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "C_TERM")
  TEST_EQUAL(peptide_id_vector[17].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA), "ANYWHERE")
  TEST_EQUAL(peptide_id_vector[17].getHits()[0].getSequence().toString(), "VILHLKEDQTEYLEER")
  TEST_EQUAL(peptide_id_vector[17].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "EYGCAPWPMVEKLIK")

  TEST_EQUAL(peptide_id_vector[289].getHits().size(), 1)
  TEST_EQUAL(peptide_id_vector[289].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "cross-link")
  TEST_EQUAL(peptide_id_vector[289].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 15)
  TEST_EQUAL(peptide_id_vector[289].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), 0)
  TEST_EQUAL(peptide_id_vector[289].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "ANYWHERE")
  TEST_EQUAL(peptide_id_vector[289].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA), "ANYWHERE")
  TEST_EQUAL(peptide_id_vector[289].getHits()[0].getSequence().toString(), "DYHFVNATEESDALAKLR")
  TEST_EQUAL(peptide_id_vector[289].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "KETFDDLPK")

  TEST_EQUAL(peptide_id_vector[279].getHits().size(), 2)
  TEST_EQUAL(peptide_id_vector[279].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TYPE), "cross-link")
  TEST_EQUAL(peptide_id_vector[279].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS1), 0)
  TEST_EQUAL(peptide_id_vector[279].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_POS2), 6)
  TEST_EQUAL(peptide_id_vector[279].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_ALPHA), "N_TERM")
  TEST_EQUAL(peptide_id_vector[279].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_XL_TERM_SPEC_BETA), "ANYWHERE")
  TEST_EQUAL(peptide_id_vector[279].getHits()[0].getSequence().toString(), "MASGSCQGCEEDEETLKK")
  TEST_EQUAL(peptide_id_vector[279].getHits()[0].getMetaValue(Constants::UserParam::OPENPEPXL_BETA_SEQUENCE), "NTEGTQKQK")

END_SECTION

END_TEST
