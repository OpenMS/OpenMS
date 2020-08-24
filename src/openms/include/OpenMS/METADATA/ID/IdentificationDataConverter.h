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
// $Maintainer: Hendrik Weisser $
// $Authors: Hendrik Weisser $
// --------------------------------------------------------------------------

#pragma once

#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/FORMAT/FASTAFile.h>
#include <OpenMS/FORMAT/MzTab.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <OpenMS/METADATA/ProteinIdentification.h>

namespace OpenMS
{
  class OPENMS_DLLAPI IdentificationDataConverter
  {
  public:

    /// Import from legacy peptide/protein identifications
    static void importIDs(IdentificationData& id_data,
                          const std::vector<ProteinIdentification>& proteins,
                          const std::vector<PeptideIdentification>& peptides);

    /// Export to legacy peptide/protein identifications
    static void exportIDs(const IdentificationData& id_data,
                          std::vector<ProteinIdentification>& proteins,
                          std::vector<PeptideIdentification>& peptides,
                          bool export_oligonucleotides = false);

    /// Export to mzTab format
    static MzTab exportMzTab(const IdentificationData& id_data);

    /// Import FASTA sequences as parent molecules
    static void importSequences(IdentificationData& id_data,
                                const std::vector<FASTAFile::FASTAEntry>& fasta,
                                IdentificationData::MoleculeType type =
                                IdentificationData::MoleculeType::PROTEIN,
                                const String& decoy_pattern = "");

  protected:

    /// Export a parent molecule (protein or nucleic acid) to mzTab
    template <typename MzTabSectionRow>
    static void exportParentMoleculeToMzTab_(
      const IdentificationData::ParentMolecule& parent,
      std::vector<MzTabSectionRow>& output,
      std::map<IdentificationData::ScoreTypeRef, Size>& score_map)
    {
      MzTabSectionRow row;
      row.accession.set(parent.accession);
      exportStepsAndScoresToMzTab_(parent.steps_and_scores, row.search_engine,
                                   row.best_search_engine_score, score_map);
      row.description.set(parent.description);
      row.coverage.set(parent.coverage);
      if (!parent.sequence.empty())
      {
        MzTabOptionalColumnEntry opt_seq;
        opt_seq.first = "opt_sequence";
        opt_seq.second.set(parent.sequence);
        row.opt_.push_back(opt_seq);
      }
      output.push_back(row);
    }

    /// Export an identified sequence (peptide or oligonucleotide, but not small molecule/compound) to mzTab
    template <typename MzTabSectionRow, typename IdentSeq>
    static void exportPeptideOrOligoToMzTab_(
      const IdentSeq& identified, std::vector<MzTabSectionRow>& output,
      std::map<IdentificationData::ScoreTypeRef, Size>& score_map)
    {
      MzTabSectionRow row;
      // @TODO: handle modifications properly
      row.sequence.set(identified.sequence.toString());
      exportStepsAndScoresToMzTab_(identified.steps_and_scores,
                                   row.search_engine,
                                   row.best_search_engine_score, score_map);
      if (identified.parent_matches.empty()) // no parent information given
      {
        // row.unique.set(false); // leave this unset?
        output.push_back(row);
      }
      else // generate entries (with duplicated data) for every accession
      {
        bool unique = (identified.parent_matches.size() == 1);
        for (const auto& match_pair : identified.parent_matches)
        {
          const String& accession = match_pair.first->accession;
          MzTabSectionRow copy = row;
          copy.accession.set(accession);
          copy.unique.set(unique);
          addMzTabMoleculeParentContext_(match_pair.second, copy);
          output.push_back(copy);
        }
      }
    }

    /// Export a molecule-query match (peptide- or oligonucleotide-spectrum match) to mzTab
    template <typename MzTabSectionRow>
    static void exportQueryMatchToMzTab_(
      const String& sequence,
      const IdentificationData::MoleculeQueryMatch& match, double calc_mass,
      std::vector<MzTabSectionRow>& output,
      std::map<IdentificationData::ScoreTypeRef, Size>& score_map,
      std::map<IdentificationData::InputFileRef, Size>& file_map)
    {
      MzTabSectionRow xsm; // PSM or OSM
      // @TODO: handle modifications properly
      xsm.sequence.set(sequence);
      exportStepsAndScoresToMzTab_(match.steps_and_scores, xsm.search_engine,
                                   xsm.search_engine_score, score_map);
      const IdentificationData::DataQuery& query = *match.data_query_ref;
      std::vector<MzTabDouble> rts(1);
      rts[0].set(query.rt);
      xsm.retention_time.set(rts);
      xsm.charge.set(match.charge);
      xsm.exp_mass_to_charge.set(query.mz);
      xsm.calc_mass_to_charge.set(calc_mass / abs(match.charge));
      if (query.input_file_opt)
      {
        xsm.spectra_ref.setMSFile(file_map[*query.input_file_opt]);
      }
      xsm.spectra_ref.setSpecRef(query.data_id);
      // @TODO: find a way of passing in the names of relevant meta values
      // (e.g. from NucleicAcidSearchEngine), instead of hard-coding them here
      static const std::vector<String> meta_out({"adduct", "isotope_offset"});
      for (const String& meta : meta_out)
      {
        if (match.metaValueExists(meta))
        {
          MzTabOptionalColumnEntry opt_meta;
          opt_meta.first = "opt_" + meta;
          opt_meta.second.set(match.getMetaValue(meta));
          xsm.opt_.push_back(opt_meta);
        }
      }
      // don't repeat data from the peptide section (e.g. accessions)
      // why are "pre"/"post"/"start"/"end" not in the peptide section?!
      output.push_back(xsm);
    }

    /// Helper function to add processing steps (search engines) and their scores to MzTab
    static void exportStepsAndScoresToMzTab_(
      const IdentificationData::AppliedProcessingSteps& steps_and_scores,
      MzTabParameterList& steps_out, std::map<Size, MzTabDouble>& scores_out,
      std::map<IdentificationData::ScoreTypeRef, Size>& score_map);

    /// Helper function to add search engine score entries to MzTab's meta data section
    static void addMzTabSEScores_(
      const std::map<IdentificationData::ScoreTypeRef, Size>& scores,
      std::map<Size, MzTabParameter>& output);

    /// Helper function for @ref exportPeptideOrOligoToMzTab_() - oligonucleotide variant
    static void addMzTabMoleculeParentContext_(
      const std::set<IdentificationData::MoleculeParentMatch>& matches,
      MzTabOligonucleotideSectionRow& row);

    /// Helper function for @ref exportPeptideOrOligoToMzTab_() - peptide variant
    static void addMzTabMoleculeParentContext_(
      const std::set<IdentificationData::MoleculeParentMatch>& matches,
      MzTabPeptideSectionRow& row);

    /// Helper function to import DB search parameters from legacy format
    static IdentificationData::SearchParamRef importDBSearchParameters_(
      const ProteinIdentification::SearchParameters& pisp,
      IdentificationData& id_data);

    /// Helper function to export DB search parameters to legacy format
    static ProteinIdentification::SearchParameters exportDBSearchParameters_(
      IdentificationData::SearchParamRef ref);
  };
}
