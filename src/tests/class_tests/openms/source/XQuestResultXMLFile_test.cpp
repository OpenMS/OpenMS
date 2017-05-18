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
// $Maintainer: Eugen Netz $
// $Authors: Eugen Netz $
// --------------------------------------------------------------------------

#include <OpenMS/CONCEPT/ClassTest.h>
#include <OpenMS/test_config.h>

///////////////////////////

#include <OpenMS/FORMAT/XQuestResultXMLFile.h>
#include <OpenMS/FORMAT/MzIdentMLFile.h>
#include <OpenMS/METADATA/ProteinIdentification.h>
#include <OpenMS/ANALYSIS/XLMS/OPXLDataStructs.h>

#include <QStringList>

using namespace OpenMS;

START_TEST(XQuestResultXMLFile, "$Id$")


START_SECTION(static void writeXQuestXML(String out_file, String base_name, const std::vector< PeptideIdentification >& peptide_ids, const std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > >& all_top_csms, const PeakMap& spectra,
                                                  String precursor_mass_tolerance_unit, String fragment_mass_tolerance_unit, double precursor_mass_tolerance, double fragment_mass_tolerance, double fragment_mass_tolerance_xlinks, String cross_link_name,
                                                  double cross_link_mass_light, DoubleList cross_link_mass_mono_link, String in_fasta, String in_decoy_fasta, StringList cross_link_residue1, StringList cross_link_residue2, double cross_link_mass_iso_shift, String enzyme_name, Size missed_cleavages))

  // TODO
  TEST_EQUAL(1, 1)

  std::vector<ProteinIdentification> protein_ids;
  std::vector<PeptideIdentification> peptide_ids;

  String mzid_input_file= OPENMS_GET_TEST_DATA_PATH("MzIdentML_XLMS_labelled.mzid");
  MzIdentMLFile().load(mzid_input_file, protein_ids, peptide_ids);

  std::vector< std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > > all_top_csms;
  std::vector< OPXLDataStructs::CrossLinkSpectrumMatch > spectrum_csms;
  PeakMap spectra;

  for (Size i = 0; i < peptide_ids.size(); ++i)
  {
    std::vector<PeptideHit> hits = peptide_ids[i].getHits();
    std::cout << peptide_ids[i].getMetaValue("spectrum_reference") << std::endl;

    OPXLDataStructs::CrossLinkSpectrumMatch csm;
    OPXLDataStructs::ProteinProteinCrossLink cross_link;
    cross_link.cross_linker_name = "MyLinker";
    std::pair <SignedSize, SignedSize> positions;

    PeakSpectrum spectrum;
    spectrum.setRT(peptide_ids[i].getRT());
    std::vector<Precursor> precs;
    Precursor prec;

    csm.scan_index_light = i;
    csm.scan_index_heavy = i;
    csm.matched_common_alpha = 0;
    csm.matched_common_beta = 0;
    csm.matched_xlink_alpha = 0;
    csm.matched_xlink_beta = 0;
    csm.HyperCommon = 0;
    csm.HyperXlink = 0;
    csm.HyperAlpha = 0;
    csm.HyperBeta = 0;
    csm.HyperBoth = 0;
    csm.PScoreCommon = 0;
    csm.PScoreXlink = 0;
    csm.PScoreAlpha = 0;
    csm.PScoreBeta = 0;
    csm.PScoreBoth = 0;

    for (Size j = 0; j < hits.size(); ++j)
    {
      if (hits[j].metaValueExists("xl_chain"))
      {
        std::cout << hits[j].getMetaValue("xl_chain") << std::endl;

        if (hits[j].getMetaValue("xl_chain") == "MS:1002509")
        {
          cross_link.alpha = hits[j].getSequence();
          positions.first = hits[j].getMetaValue("xl_pos").toString().toInt();

          if (hits[j].getMetaValue("xl_type") == "mono-link")
          {
            positions.second = -1;
          }
          else if (hits[j].getMetaValue("xl_type") == "loop-link")
          {
            positions.second = hits[j].getMetaValue("xl_pos2").toString().toInt();
          }

          std::cout << "POSITIONS true: " << hits[j].getMetaValue("xl_pos").toString().toInt() << " | written: " << positions.first << std::endl;

          prec.setCharge(hits[j].getCharge());
          prec.setMZ(peptide_ids[i].getMZ());
          precs.push_back(prec);

          csm.score = hits[j].getScore();
          csm.rank  = hits[j].getMetaValue("xl_rank");
          csm.pre_score = i*0.1;
          csm.percTIC = i*0.2;
          csm.wTIC = hits[j].getMetaValue("OpenXQuest:wTIC");
          csm.int_sum = hits[j].getMetaValue("OpenXQuest:intsum");
          csm.match_odds = hits[j].getMetaValue("OpenXQuest:match-odds");
          csm.xcorrx_max = hits[j].getMetaValue("OpenXQuest:xcorr xlink");
          csm.xcorrc_max = hits[j].getMetaValue("OpenXQuest:xcorr common");
          csm.peptide_id_index = i;

          if (hits[j].getMetaValue("xl_type") == "cross-link" || hits[j].getMetaValue("xl_type") == "loop-link")
          {
            cross_link.cross_linker_mass = 150;
          }
          else if (hits[j].getMetaValue("xl_type") == "mono-link")
          {
            cross_link.cross_linker_mass = 50;
          }
        }
        else
        {
          cross_link.beta = hits[j].getSequence();
          positions.second = hits[j].getMetaValue("xl_pos").toString().toInt();
          std::cout << "POSITIONS true: " << hits[j].getMetaValue("xl_pos").toString().toInt() << " | written: " << positions.second << std::endl;
        }
      }
      else
      {
        std::cout << "WTF?" << std::endl;
       }
    }
    spectrum.setPrecursors(precs);
    spectra.addSpectrum(spectrum);

    cross_link.cross_link_position = positions;
    csm.cross_link = cross_link;
    spectrum_csms.push_back(csm);
  }
  all_top_csms.push_back(spectrum_csms);

  std::cout << "after loops" << std::endl;

  QStringList q_str_list1;
  QStringList q_str_list2;
  q_str_list1 << "K" << "E";
  q_str_list2 << "D" << "E";
  StringList cross_link_residue1 = StringListUtils::fromQStringList(q_str_list1);
  StringList cross_link_residue2 = StringListUtils::fromQStringList(q_str_list2);

  std::vector<String> mono_masses;
  mono_masses.push_back("50.0");
  DoubleList cross_link_mass_mono_link = ListUtils::create<double>(mono_masses);

  double cross_link_mass_light = 150.0;
  double cross_link_mass_iso_shift = 15.0;
  String base_name = "TEST_SPECTRUM";
  String in_fasta = "database.fasta";
  String in_decoy_fasta = "decoy_database.fasta";
  String enzyme_name = "Trypsin";
  Size missed_cleavages = 2;
  String cross_link_name = "MyLinker";

  String precursor_mass_tolerance_unit = "ppm";
  String fragment_mass_tolerance_unit = "Da";
  double precursor_mass_tolerance = 10;
  double fragment_mass_tolerance = 0.2;
  double fragment_mass_tolerance_xlinks = 0.3;

  String out_file;
  NEW_TMP_FILE(out_file)

  XQuestResultXMLFile::writeXQuestXML(out_file, base_name, peptide_ids, all_top_csms, spectra,
                                                  precursor_mass_tolerance_unit, fragment_mass_tolerance_unit, precursor_mass_tolerance, fragment_mass_tolerance, fragment_mass_tolerance_xlinks, cross_link_name,
                                                  cross_link_mass_light, cross_link_mass_mono_link, in_fasta, in_decoy_fasta, cross_link_residue1, cross_link_residue2, cross_link_mass_iso_shift, enzyme_name, missed_cleavages);

  std::vector< std::vector< PeptideIdentification > > peptide_id_vector_vector;
  std::vector< ProteinIdentification > protein_id_vector;
  XQuestResultXMLFile().load(out_file, peptide_id_vector_vector, protein_id_vector, 0, true);

  for (Size j = 0; j < peptide_id_vector_vector.size(); ++j)
  {
    std::vector< PeptideIdentification > peptide_id_vector = peptide_id_vector_vector[j];
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
  }


END_SECTION

END_TEST


