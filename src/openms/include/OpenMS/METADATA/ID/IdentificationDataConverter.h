// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
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
  class FeatureMap;

  class OPENMS_DLLAPI IdentificationDataConverter
  {
  public:

    /// Import from legacy peptide/protein identifications
    static void importIDs(IdentificationData& id_data,
                          const std::vector<ProteinIdentification>& proteins,
                          const std::vector<PeptideIdentification>& peptides);

    /*!
      @brief Export to legacy peptide/protein identifications

      Results are added to existing data (if any) in @p proteins and @p peptides.
    */
    static void exportIDs(const IdentificationData& id_data,
                          std::vector<ProteinIdentification>& proteins,
                          std::vector<PeptideIdentification>& peptides,
                          bool export_ids_wo_scores = false);

    /// Export to mzTab format
    static MzTab exportMzTab(const IdentificationData& id_data);

    /// Import FASTA sequences as parent sequences
    static void importSequences(IdentificationData& id_data,
                                const std::vector<FASTAFile::FASTAEntry>& fasta,
                                IdentificationData::MoleculeType type =
                                IdentificationData::MoleculeType::PROTEIN,
                                const String& decoy_pattern = "");

    /// Convert parent matches to peptide evidences
    static void exportParentMatches(
      const IdentificationData::ParentMatches& parent_matches, PeptideHit& hit);

    /*!
      @brief Convert IDs from legacy peptide/protein identifications in a feature map

      @param features Feature map containing IDs in legacy format
      @param clear_original Clear original IDs after conversion?
    */
    static void importFeatureIDs(FeatureMap& features, bool clear_original = true);

    /*!
      @brief Convert IDs in a feature map to legacy peptide/protein identifications

      @param features Feature map containing IDs in new format
      @param clear_original Clear original IDs after conversion?
    */
    static void exportFeatureIDs(FeatureMap& features, bool clear_original = true);

    /*!
      @brief Convert IDs from legacy peptide/protein identifications in a consensus map

      @param consensus Consensus map containing IDs in legacy format
      @param clear_original Clear original IDs after conversion?
    */
    static void importConsensusIDs(ConsensusMap& consensus, bool clear_original = true);

    /*!
      @brief Convert IDs in a consensus map to legacy peptide/protein identifications

      @param consensus Consensus map containing IDs in new format
      @param clear_original Clear original IDs after conversion?
    */
    static void exportConsensusIDs(ConsensusMap& consensus, bool clear_original = true);

  protected:

    using StepOpt = std::optional<IdentificationData::ProcessingStepRef>;

    /// Functor for ordering @p StepOpt (by date of the steps, if available):
    struct StepOptCompare
    {
      bool operator()(const StepOpt& left, const StepOpt& right) const
      {
        // @TODO: should runs without associated step go first or last?
        if (!left) return bool(right);
        if (!right) return false;
        return **left < **right;
      }
    };

    /// Functor for ordering peptide IDs by RT and m/z (if available)
    struct PepIDCompare
    {
      bool operator()(const PeptideIdentification& left,
                      const PeptideIdentification& right) const
      {
        // @TODO: should IDs without RT go first or last?
        if (left.hasRT())
        {
          if (right.hasRT())
          {
            if (right.getRT() != left.getRT())
            {
              return left.getRT() < right.getRT();
            } // else: compare by m/z (below)
          }
          else
          {
            return false;
          }
        }
        else if (right.hasRT())
        {
          return true;
        }
        // no RTs or same RTs -> try to compare by m/z:
        if (left.hasMZ())
        {
          if (right.hasMZ())
          {
            return left.getMZ() < right.getMZ();
          }
          else
          {
            return false;
          }
        }
        // if both PI's have nothing, return false (to ensure 'x < x' is false for strict weak ordering)
        return right.hasMZ();
      }
    };

    /// Export a parent sequence (protein or nucleic acid) to mzTab
    template <typename MzTabSectionRow>
    static void exportParentSequenceToMzTab_(
      const IdentificationData::ParentSequence& parent,
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
        // in mzTab, "unique" means "peptide is unique for this protein"
        row.unique.set(identified.parent_matches.size() == 1);
        for (const auto& match_pair : identified.parent_matches)
        {
          row.accession.set(match_pair.first->accession);
          for (const IdentificationData::ParentMatch& match :
                 match_pair.second)
          {
            MzTabSectionRow copy = row;
            addMzTabMoleculeParentContext_(match, copy);
            output.push_back(copy);
          }
        }
      }
    }

    /// Export an input match (peptide- or oligonucleotide-spectrum match) to mzTab
    template <typename MzTabSectionRow>
    static void exportObservationMatchToMzTab_(
      const String& sequence,
      const IdentificationData::ObservationMatch& match, double calc_mass,
      std::vector<MzTabSectionRow>& output,
      std::map<IdentificationData::ScoreTypeRef, Size>& score_map,
      std::map<IdentificationData::InputFileRef, Size>& file_map)
    {
      MzTabSectionRow xsm; // PSM or OSM
      // @TODO: handle modifications properly
      xsm.sequence.set(sequence);
      exportStepsAndScoresToMzTab_(match.steps_and_scores, xsm.search_engine,
                                   xsm.search_engine_score, score_map);
      const IdentificationData::Observation& query = *match.observation_ref;
      std::vector<MzTabDouble> rts(1);
      rts[0].set(query.rt);
      xsm.retention_time.set(rts);
      xsm.charge.set(match.charge);
      xsm.exp_mass_to_charge.set(query.mz);
      xsm.calc_mass_to_charge.set(calc_mass / abs(match.charge));
      xsm.spectra_ref.setMSFile(file_map[query.input_file]);
      xsm.spectra_ref.setSpecRef(query.data_id);
      // optional column for adduct:
      if (match.adduct_opt)
      {
        MzTabOptionalColumnEntry opt_adduct;
        opt_adduct.first = "opt_adduct";
        opt_adduct.second.set((*match.adduct_opt)->getName());
        xsm.opt_.push_back(opt_adduct);
      }
      // optional columns for isotope offset:
      // @TODO: find a way of passing in the names of relevant meta values
      // (e.g. from NucleicAcidSearchEngine), instead of hard-coding them here
      if (match.metaValueExists("isotope_offset"))
      {
        MzTabOptionalColumnEntry opt_meta;
        opt_meta.first = "opt_isotope_offset";
        opt_meta.second.set(match.getMetaValue("isotope_offset"));
        xsm.opt_.push_back(opt_meta);
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
      const IdentificationData::ParentMatch& match,
      MzTabOligonucleotideSectionRow& row);

    /// Helper function for @ref exportPeptideOrOligoToMzTab_() - peptide variant
    static void addMzTabMoleculeParentContext_(
      const IdentificationData::ParentMatch& match,
      MzTabPeptideSectionRow& row);

    /// Helper function to import DB search parameters from legacy format
    static IdentificationData::SearchParamRef importDBSearchParameters_(
      const ProteinIdentification::SearchParameters& pisp,
      IdentificationData& id_data);

    /// Helper function to export DB search parameters to legacy format
    static ProteinIdentification::SearchParameters exportDBSearchParameters_(
      IdentificationData::SearchParamRef ref);

    /// Helper function to export (primary) MS run information to legacy format
    static void exportMSRunInformation_(
      IdentificationData::ProcessingStepRef step_ref,
      ProteinIdentification& protein);

    static void handleFeatureImport_(Feature& feature, const IntList& indexes,
                                     std::vector<PeptideIdentification>& peptides,
                                     Size& id_counter, bool clear_original);

    static void handleFeatureExport_(Feature& feature, const IntList& indexes,
                                     IdentificationData& id_data, Size& id_counter);
  };
}
