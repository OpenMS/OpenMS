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

#include <OpenMS/METADATA/ID/IdentificationDataConverter.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <OpenMS/CONCEPT/Constants.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/FORMAT/FileHandler.h>

using namespace std;

using ID = OpenMS::IdentificationData;

namespace OpenMS
{
  void IdentificationDataConverter::importIDs(
    IdentificationData& id_data, const vector<ProteinIdentification>& proteins,
    const vector<PeptideIdentification>& peptides)
  {
    map<String, ID::ProcessingStepRef> id_to_step;
    ProgressLogger progresslogger;
    progresslogger.setLogType(ProgressLogger::CMD);

    // ProteinIdentification:
    progresslogger.startProgress(0, proteins.size(),
                                 "converting protein identification runs");
    Size proteins_counter = 0;
    for (const ProteinIdentification& prot : proteins)
    {
      proteins_counter++;
      progresslogger.setProgress(proteins_counter);
      ID::ScoreType score_type(prot.getScoreType(), prot.isHigherScoreBetter());
      ID::ScoreTypeRef prot_score_ref = id_data.registerScoreType(score_type);

      ID::DataProcessingSoftware software(prot.getSearchEngine(),
                                          prot.getSearchEngineVersion());
      software.assigned_scores.push_back(prot_score_ref);
      ID::ProcessingSoftwareRef software_ref =
        id_data.registerDataProcessingSoftware(software);

      ID::SearchParamRef search_ref =
        importDBSearchParameters_(prot.getSearchParameters(), id_data);

      ID::DataProcessingStep step(software_ref);
      // ideally, this should give us the raw files:
      prot.getPrimaryMSRunPath(step.primary_files, true);
      // ... and this should give us mzML files:
      vector<String> spectrum_files;
      prot.getPrimaryMSRunPath(spectrum_files);
      for (const String& path : spectrum_files)
      {
        ID::InputFileRef file_ref = id_data.registerInputFile(path);
        step.input_file_refs.push_back(file_ref);
      }
      step.date_time = prot.getDateTime();
      ID::ProcessingStepRef step_ref =
        id_data.registerDataProcessingStep(step, search_ref);
      id_to_step[prot.getIdentifier()] = step_ref;
      id_data.setCurrentProcessingStep(step_ref);

      ProgressLogger sublogger;
      sublogger.setLogType(ProgressLogger::CMD);
      String run_label = "(run " + String(proteins_counter) + "/" +
        String(proteins.size()) + ")";

      // ProteinHit:
      sublogger.startProgress(0, prot.getHits().size(),
                              "converting protein hits " + run_label);
      Size hits_counter = 0;
      for (const ProteinHit& hit : prot.getHits())
      {
        ++hits_counter;
        sublogger.setProgress(hits_counter);
        ID::ParentMolecule parent(hit.getAccession());
        parent.sequence = hit.getSequence();
        parent.description = hit.getDescription();
        // coverage comes in percents, -1 for missing; we want 0 to 1:
        parent.coverage = max(hit.getCoverage(), 0.0) / 100.0;
        static_cast<MetaInfoInterface&>(parent) = hit;
        ID::AppliedProcessingStep applied(step_ref);
        applied.scores[prot_score_ref] = hit.getScore();
        parent.steps_and_scores.push_back(applied);
        id_data.registerParentMolecule(parent);
      }
      sublogger.endProgress();

      // indistinguishable protein groups:
      if (!prot.getIndistinguishableProteins().empty())
      {
        sublogger.startProgress(0, prot.getIndistinguishableProteins().size(),
                                "converting indistinguishable proteins " +
                                run_label);
        Size groups_counter = 0;

        ID::ScoreType score("probability", true);
        ID::ScoreTypeRef score_ref = id_data.registerScoreType(score);

        ID::ParentMoleculeGrouping grouping;
        grouping.label = "indistinguishable proteins";

        for (const auto& group : prot.getIndistinguishableProteins())
        {
          ++groups_counter;
          sublogger.setProgress(groups_counter);
          ID::ParentMoleculeGroup new_group;
          new_group.scores[score_ref] = group.probability;
          for (const String& acc : group.accessions)
          {
            ID::ParentMolecule parent(acc);
            ID::ParentMoleculeRef ref = id_data.registerParentMolecule(parent);
            new_group.parent_molecule_refs.insert(ref);
          }
          grouping.groups.insert(new_group);
        }
        id_data.registerParentMoleculeGrouping(grouping);
        sublogger.endProgress();
      }
      // other protein groups:
      if (!prot.getProteinGroups().empty())
      {
        sublogger.startProgress(0, prot.getProteinGroups().size(),
                                "converting protein groups " + run_label);
        Size groups_counter = 0;

        ID::ScoreType score("probability", true);
        ID::ScoreTypeRef score_ref = id_data.registerScoreType(score);

        ID::ParentMoleculeGrouping grouping;
        grouping.label = "protein groups";

        for (const auto& group : prot.getProteinGroups())
        {
          ++groups_counter;
          sublogger.setProgress(groups_counter);
          ID::ParentMoleculeGroup new_group;
          new_group.scores[score_ref] = group.probability;
          for (const String& acc : group.accessions)
          {
            ID::ParentMolecule parent(acc);
            ID::ParentMoleculeRef ref = id_data.registerParentMolecule(parent);
            new_group.parent_molecule_refs.insert(ref);
          }
          grouping.groups.insert(new_group);
        }
        id_data.registerParentMoleculeGrouping(grouping);
        sublogger.endProgress();
      }

      id_data.clearCurrentProcessingStep();
    }
    progresslogger.endProgress();

    // PeptideIdentification:
    Size unknown_query_counter = 1;
    progresslogger.startProgress(0, peptides.size(),
                                 "converting peptide identifications");
    Size peptides_counter = 0;
    for (const PeptideIdentification& pep : peptides)
    {
      peptides_counter++;
      progresslogger.setProgress(peptides_counter);
      const String& id = pep.getIdentifier();
      ID::ProcessingStepRef step_ref = id_to_step[id];
      ID::DataQuery query(""); // fill in "data_id" later
      if (!step_ref->input_file_refs.empty())
      {
        // @TODO: what if there's more than one input file?
        query.input_file_opt = step_ref->input_file_refs[0];
      }
      else
      {
        String file = "UNKNOWN_INPUT_FILE_" + id;
        ID::InputFileRef file_ref = id_data.registerInputFile(file);
        query.input_file_opt = file_ref;
      }
      query.rt = pep.getRT();
      query.mz = pep.getMZ();
      static_cast<MetaInfoInterface&>(query) = pep;
      if (pep.metaValueExists("spectrum_reference"))
      {
        query.data_id = pep.getMetaValue("spectrum_reference");
        query.removeMetaValue("spectrum_reference");
      }
      else
      {
        if (pep.hasRT() && pep.hasMZ())
        {
          query.data_id = String("RT=") + String(float(query.rt)) + "_MZ=" +
            String(float(query.mz));
        }
        else
        {
          query.data_id = "UNKNOWN_QUERY_" + String(unknown_query_counter);
          ++unknown_query_counter;
        }
      }
      ID::DataQueryRef query_ref = id_data.registerDataQuery(query);

      ID::ScoreType score_type(pep.getScoreType(), pep.isHigherScoreBetter());
      ID::ScoreTypeRef score_ref = id_data.registerScoreType(score_type);

      // PeptideHit:
      for (const PeptideHit& hit : pep.getHits())
      {
        if (hit.getSequence().empty()) continue;
        ID::IdentifiedPeptide peptide(hit.getSequence());
        peptide.addProcessingStep(step_ref);
        for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
        {
          const String& accession = evidence.getProteinAccession();
          if (accession.empty()) continue;
          ID::ParentMolecule parent(accession);
          parent.addProcessingStep(step_ref);
          // this will merge information if the protein already exists:
          ID::ParentMoleculeRef parent_ref =
            id_data.registerParentMolecule(parent);
          ID::MoleculeParentMatch match(evidence.getStart(), evidence.getEnd(),
                                        evidence.getAABefore(),
                                        evidence.getAAAfter());
          peptide.parent_matches[parent_ref].insert(match);
        }
        ID::IdentifiedPeptideRef peptide_ref =
          id_data.registerIdentifiedPeptide(peptide);

        ID::MoleculeQueryMatch match(peptide_ref, query_ref);
        match.charge = hit.getCharge();
        static_cast<MetaInfoInterface&>(match) = hit;
        if (!hit.getPeakAnnotations().empty())
        {
          match.peak_annotations[step_ref] = hit.getPeakAnnotations();
        }
        ID::AppliedProcessingStep applied(step_ref);
        applied.scores[score_ref] = hit.getScore();

        // analysis results from pepXML:
        for (const PeptideHit::PepXMLAnalysisResult& ana_res :
               hit.getAnalysisResults())
        {
          ID::DataProcessingSoftware software;
          software.setName(ana_res.score_type); // e.g. "peptideprophet"
          ID::AppliedProcessingStep sub_applied;
          ID::ScoreType main_score;
          main_score.cv_term.setName(ana_res.score_type + "_probability");
          main_score.higher_better = ana_res.higher_is_better;
          ID::ScoreTypeRef main_score_ref =
            id_data.registerScoreType(main_score);
          software.assigned_scores.push_back(main_score_ref);
          sub_applied.scores[main_score_ref] = ana_res.main_score;
          for (const pair<const String, double>& sub_pair : ana_res.sub_scores)
          {
            ID::ScoreType sub_score;
            sub_score.cv_term.setName(sub_pair.first);
            ID::ScoreTypeRef sub_score_ref =
              id_data.registerScoreType(sub_score);
            software.assigned_scores.push_back(sub_score_ref);
            sub_applied.scores[sub_score_ref] = sub_pair.second;
          }
          ID::ProcessingSoftwareRef software_ref =
            id_data.registerDataProcessingSoftware(software);
          ID::DataProcessingStep sub_step(software_ref);
          if (query.input_file_opt)
          {
            sub_step.input_file_refs.push_back(*query.input_file_opt);
          }
          ID::ProcessingStepRef sub_step_ref =
            id_data.registerDataProcessingStep(sub_step);
          sub_applied.processing_step_opt = sub_step_ref;
          match.addProcessingStep(sub_applied);
        }

        // most recent step (with primary score) goes last:
        match.addProcessingStep(applied);
        id_data.registerMoleculeQueryMatch(match);
      }
    }
    progresslogger.endProgress();
  }


  void IdentificationDataConverter::exportIDs(
    const IdentificationData& id_data, vector<ProteinIdentification>& proteins,
    vector<PeptideIdentification>& peptides, bool export_oligonucleotides)
  {
    proteins.clear();
    peptides.clear();

    // "DataQuery" roughly corresponds to "PeptideIdentification",
    // "DataProcessingStep" roughly corresponds to "ProteinIdentification";
    // score type is stored in "PeptideIdent.", not "PeptideHit":
    map<pair<ID::DataQueryRef, boost::optional<ID::ProcessingStepRef>>,
        pair<vector<PeptideHit>, ID::ScoreTypeRef>> psm_data;
    // we only export peptides and proteins (or oligos and RNAs), so start by
    // getting the PSMs (or OSMs):
    const String& ppm_error_name = Constants::UserParam::PRECURSOR_ERROR_PPM_USERPARAM;

    for (const ID::MoleculeQueryMatch& query_match :
           id_data.getMoleculeQueryMatches())
    {
      PeptideHit hit;
      const ID::ParentMatches* parent_matches_ptr;
      if (!export_oligonucleotides) // export peptides
      {
        if (query_match.getMoleculeType() != ID::MoleculeType::PROTEIN) continue;
        static_cast<MetaInfoInterface&>(hit) = query_match;
        ID::IdentifiedPeptideRef peptide_ref =
          query_match.getIdentifiedPeptideRef();
        hit.setSequence(peptide_ref->sequence);
        parent_matches_ptr = &(peptide_ref->parent_matches);
      }
      else
      {
        if (query_match.getMoleculeType() != ID::MoleculeType::RNA) continue;
        static_cast<MetaInfoInterface&>(hit) = query_match;
        ID::IdentifiedOligoRef oligo_ref = query_match.getIdentifiedOligoRef();
        hit.setMetaValue("label", oligo_ref->sequence.toString());
        parent_matches_ptr = &(oligo_ref->parent_matches);
      }
      hit.setCharge(query_match.charge);
      if (query_match.metaValueExists(ppm_error_name))
      {
        hit.setMetaValue(ppm_error_name,
                         query_match.getMetaValue(ppm_error_name));
      }
      for (const auto& pair : *parent_matches_ptr)
      {
        ID::ParentMoleculeRef parent_ref = pair.first;
        for (const ID::MoleculeParentMatch& parent_match : pair.second)
        {
          PeptideEvidence evidence;
          evidence.setProteinAccession(parent_ref->accession);
          evidence.setStart(parent_match.start_pos);
          evidence.setEnd(parent_match.end_pos);
          if (!parent_match.left_neighbor.empty())
          {
            evidence.setAABefore(parent_match.left_neighbor[0]);
          }
          if (!parent_match.right_neighbor.empty())
          {
            evidence.setAAAfter(parent_match.right_neighbor[0]);
          }
          hit.addPeptideEvidence(evidence);
        }
      }
      // sort the evidences:
      vector<PeptideEvidence> evidences = hit.getPeptideEvidences();
      sort(evidences.begin(), evidences.end());
      hit.setPeptideEvidences(evidences);

      // generate hits in different ID runs for different processing steps:
      for (ID::AppliedProcessingStep applied : query_match.steps_and_scores)
      {
        if (applied.scores.empty()) continue;
        PeptideHit hit_copy = hit;
        vector<pair<ID::ScoreTypeRef, double>> scores =
          applied.getScoresInOrder();
        // "primary" score comes first:
        hit_copy.setScore(scores[0].second);
        // add meta values for "secondary" scores:
        for (auto it = ++scores.begin(); it != scores.end(); ++it)
        {
          hit_copy.setMetaValue(it->first->cv_term.getName(), it->second);
        }
        auto pos =
          query_match.peak_annotations.find(applied.processing_step_opt);
        if (pos != query_match.peak_annotations.end())
        {
          hit_copy.setPeakAnnotations(pos->second);
        }
        auto key = make_pair(query_match.data_query_ref,
                             applied.processing_step_opt);
        psm_data[key].first.push_back(hit_copy);
        psm_data[key].second = scores[0].first; // primary score type
      }
    }

    set<boost::optional<ID::ProcessingStepRef>> steps;
    for (const auto& psm : psm_data)
    {
      const ID::DataQuery& query = *psm.first.first;
      PeptideIdentification peptide;
      static_cast<MetaInfoInterface&>(peptide) = query;
      peptide.setRT(query.rt);
      peptide.setMZ(query.mz);
      peptide.setMetaValue("spectrum_reference", query.data_id);
      peptide.setHits(psm.second.first);
      const ID::ScoreType& score_type = *psm.second.second;
      peptide.setScoreType(score_type.cv_term.getName());
      peptide.setHigherScoreBetter(score_type.higher_better);
      if (psm.first.second) // processing step given
      {
        peptide.setIdentifier(String(Size(&(**psm.first.second))));
      }
      else
      {
        peptide.setIdentifier("dummy");
      }
      peptides.push_back(peptide);
      steps.insert(psm.first.second);
    }

    map<boost::optional<ID::ProcessingStepRef>,
        pair<vector<ProteinHit>, ID::ScoreTypeRef>> prot_data;
    for (const auto& parent : id_data.getParentMolecules())
    {
      bool right_type =
        parent.molecule_type == (export_oligonucleotides ?
                                 ID::MoleculeType::RNA :
                                 ID::MoleculeType::PROTEIN);
      if (!right_type) continue;
      ProteinHit hit;
      hit.setAccession(parent.accession);
      hit.setSequence(parent.sequence);
      hit.setDescription(parent.description);
      if (parent.coverage > 0.0)
      {
        hit.setCoverage(parent.coverage * 100.0); // convert to percents
      }
      else // zero coverage means coverage is unknown
      {
        hit.setCoverage(ProteinHit::COVERAGE_UNKNOWN);
      }
      static_cast<MetaInfoInterface&>(hit) = parent;
      if (!parent.metaValueExists("target_decoy"))
      {
        hit.setMetaValue("target_decoy", parent.is_decoy ? "decoy" : "target");
      }

      // generate hits in different ID runs for different processing steps:
      for (ID::AppliedProcessingStep applied : parent.steps_and_scores)
      {
        if (applied.scores.empty() && !steps.count(applied.processing_step_opt))
        {
          continue; // no scores and no associated peptides -> skip
        }
        ProteinHit hit_copy = hit;
        if (!applied.scores.empty())
        {
          vector<pair<ID::ScoreTypeRef, double>> scores =
            applied.getScoresInOrder();
          // "primary" score comes first:
          hit_copy.setScore(scores[0].second);
          // add meta values for "secondary" scores:
          for (auto it = ++scores.begin(); it != scores.end(); ++it)
          {
            hit_copy.setMetaValue(it->first->cv_term.getName(), it->second);
          }
          prot_data[applied.processing_step_opt].first.push_back(hit_copy);
          prot_data[applied.processing_step_opt].second = scores[0].first;
          continue;
        }
        // always include steps that have generated peptides:
        auto pos = prot_data.find(applied.processing_step_opt);
        if (pos != prot_data.end())
        {
          pos->second.first.push_back(hit);
          // existing entry, don't overwrite score type
        }
        else
        {
          prot_data[applied.processing_step_opt].first.push_back(hit);
          prot_data[applied.processing_step_opt].second =
            id_data.getScoreTypes().end(); // no score given
        }
      }
    }

    for (const auto& step_ref_opt : steps)
    {
      ProteinIdentification protein;
      if (!step_ref_opt) // no processing step given
      {
        protein.setIdentifier("dummy");
      }
      else
      {
        ID::ProcessingStepRef step_ref = *step_ref_opt;
        protein.setIdentifier(String(Size(&(*step_ref))));
        protein.setDateTime(step_ref->date_time);
        exportMSRunInformation_(step_ref, protein);
        const Software& software = *step_ref->software_ref;
        protein.setSearchEngine(software.getName());
        protein.setSearchEngineVersion(software.getVersion());
        ID::DBSearchSteps::const_iterator ss_pos =
          id_data.getDBSearchSteps().find(step_ref);
        if (ss_pos != id_data.getDBSearchSteps().end())
        {
          protein.setSearchParameters(exportDBSearchParameters_(ss_pos->
                                                                second));
        }
      }
      auto pd_pos = prot_data.find(step_ref_opt);
      if (pd_pos != prot_data.end())
      {
        protein.setHits(pd_pos->second.first);
        if (pd_pos->second.second != id_data.getScoreTypes().end())
        {
          const ID::ScoreType& score_type = *pd_pos->second.second;
          protein.setScoreType(score_type.cv_term.getName());
          protein.setHigherScoreBetter(score_type.higher_better);
        }
      }

      // protein groups:
      for (const auto& grouping : id_data.getParentMoleculeGroupings())
      {
        // do these protein groups belong to the current search run?
        if (grouping.getStepsAndScoresByStep().find(step_ref_opt) !=
            grouping.getStepsAndScoresByStep().end())
        {
          for (const auto& group : grouping.groups)
          {
            ProteinIdentification::ProteinGroup new_group;
            if (!group.scores.empty())
            {
              // @TODO: what if there are several scores?
              new_group.probability = group.scores.begin()->second;
            }
            for (const auto& parent_ref : group.parent_molecule_refs)
            {
              new_group.accessions.push_back(parent_ref->accession);
            }
            sort(new_group.accessions.begin(), new_group.accessions.end());
            if (grouping.label == "indistinguishable proteins")
            {
              protein.insertIndistinguishableProteins(new_group);
            }
            else
            {
              protein.insertProteinGroup(new_group);
            }
          }
        }
      }
      sort(protein.getIndistinguishableProteins().begin(),
           protein.getIndistinguishableProteins().end());
      sort(protein.getProteinGroups().begin(),
           protein.getProteinGroups().end());
      proteins.push_back(protein);
    }
  }


  MzTab IdentificationDataConverter::exportMzTab(const IdentificationData&
                                                 id_data)
  {
    MzTabMetaData meta;
    Size counter = 1;
    for (const auto& software : id_data.getDataProcessingSoftwares())
    {
      MzTabSoftwareMetaData sw_meta;
      sw_meta.software.setName(software.getName());
      sw_meta.software.setValue(software.getVersion());
      meta.software[counter] = sw_meta;
      ++counter;
    }
    counter = 1;
    map<ID::InputFileRef, Size> file_map;
    for (auto it = id_data.getInputFiles().begin();
         it != id_data.getInputFiles().end(); ++it)
    {
      MzTabMSRunMetaData run_meta;
      run_meta.location.set(*it);
      meta.ms_run[counter] = run_meta;
      file_map[it] = counter;
      ++counter;
    }
    set<String> fixed_mods, variable_mods;
    for (const auto& search_param : id_data.getDBSearchParams())
    {
      fixed_mods.insert(search_param.fixed_mods.begin(),
                        search_param.fixed_mods.end());
      variable_mods.insert(search_param.variable_mods.begin(),
                           search_param.variable_mods.end());
    }
    counter = 1;
    for (const String& mod : fixed_mods)
    {
      MzTabModificationMetaData mod_meta;
      mod_meta.modification.setName(mod);
      meta.fixed_mod[counter] = mod_meta;
      ++counter;
    }
    counter = 1;
    for (const String& mod : variable_mods)
    {
      MzTabModificationMetaData mod_meta;
      mod_meta.modification.setName(mod);
      meta.variable_mod[counter] = mod_meta;
      ++counter;
    }

    map<ID::ScoreTypeRef, Size> protein_scores, peptide_scores,
      psm_scores, nucleic_acid_scores, oligonucleotide_scores, osm_scores;
    // compound_scores;

    MzTabProteinSectionRows proteins;
    MzTabNucleicAcidSectionRows nucleic_acids;
    for (const auto& parent : id_data.getParentMolecules())
    {
      if (parent.molecule_type == ID::MoleculeType::PROTEIN)
      {
        exportParentMoleculeToMzTab_(parent, proteins, protein_scores);
      }
      else if (parent.molecule_type == ID::MoleculeType::RNA)
      {
        exportParentMoleculeToMzTab_(parent, nucleic_acids,
                                     nucleic_acid_scores);
      }
    }

    MzTabPeptideSectionRows peptides;
    for (const auto& peptide : id_data.getIdentifiedPeptides())
    {
      exportPeptideOrOligoToMzTab_(peptide, peptides, peptide_scores);
    }

    MzTabOligonucleotideSectionRows oligos;
    for (const auto& oligo : id_data.getIdentifiedOligos())
    {
      exportPeptideOrOligoToMzTab_(oligo, oligos, oligonucleotide_scores);
    }

    MzTabPSMSectionRows psms;
    MzTabOSMSectionRows osms;
    for (const auto& query_match : id_data.getMoleculeQueryMatches())
    {
      ID::IdentifiedMoleculeRef molecule_ref =
        query_match.identified_molecule_ref;
      // @TODO: what about small molecules?
      ID::MoleculeType molecule_type = query_match.getMoleculeType();
      if (molecule_type == ID::MoleculeType::PROTEIN)
      {
        const AASequence& seq = query_match.getIdentifiedPeptideRef()->sequence;
        double calc_mass = seq.getMonoWeight(Residue::Full, query_match.charge);
        exportQueryMatchToMzTab_(seq.toString(), query_match, calc_mass, psms,
                                 psm_scores, file_map);
        // "PSM_ID" field is set at the end, after sorting
      }
      else if (molecule_type == ID::MoleculeType::RNA)
      {
        const NASequence& seq = query_match.getIdentifiedOligoRef()->sequence;
        double calc_mass = seq.getMonoWeight(NASequence::Full,
                                             query_match.charge);
        exportQueryMatchToMzTab_(seq.toString(), query_match, calc_mass, osms,
                                 osm_scores, file_map);
      }
    }

    addMzTabSEScores_(protein_scores, meta.protein_search_engine_score);
    addMzTabSEScores_(peptide_scores, meta.peptide_search_engine_score);
    addMzTabSEScores_(psm_scores, meta.psm_search_engine_score);
    addMzTabSEScores_(nucleic_acid_scores,
                      meta.nucleic_acid_search_engine_score);
    addMzTabSEScores_(oligonucleotide_scores,
                      meta.oligonucleotide_search_engine_score);
    addMzTabSEScores_(osm_scores, meta.osm_search_engine_score);

    // sort rows:
    sort(proteins.begin(), proteins.end(),
         MzTabProteinSectionRow::RowCompare());
    sort(peptides.begin(), peptides.end(),
         MzTabPeptideSectionRow::RowCompare());
    sort(psms.begin(), psms.end(), MzTabPSMSectionRow::RowCompare());
    // set "PSM_ID" fields - the IDs are repeated for peptides with multiple
    // protein accessions; however, we don't write out the accessions on PSM
    // level (only on peptide level), so just number consecutively:
    for (Size i = 0; i < psms.size(); ++i)
    {
      psms[i].PSM_ID.set(i + 1);
    }
    sort(nucleic_acids.begin(), nucleic_acids.end(),
         MzTabNucleicAcidSectionRow::RowCompare());
    sort(oligos.begin(), oligos.end(),
         MzTabOligonucleotideSectionRow::RowCompare());
    sort(osms.begin(), osms.end(), MzTabOSMSectionRow::RowCompare());

    MzTab output;
    output.setMetaData(meta);
    output.setProteinSectionRows(proteins);
    output.setPeptideSectionRows(peptides);
    output.setPSMSectionRows(psms);
    output.setNucleicAcidSectionRows(nucleic_acids);
    output.setOligonucleotideSectionRows(oligos);
    output.setOSMSectionRows(osms);

    return output;
  }


  void IdentificationDataConverter::importSequences(
    IdentificationData& id_data, const vector<FASTAFile::FASTAEntry>& fasta,
    ID::MoleculeType type, const String& decoy_pattern)
  {
    for (const FASTAFile::FASTAEntry& entry : fasta)
    {
      ID::ParentMolecule parent(entry.identifier, type, entry.sequence,
                                entry.description);
      if (!decoy_pattern.empty() &&
          entry.identifier.hasSubstring(decoy_pattern))
      {
        parent.is_decoy = true;
      }
      id_data.registerParentMolecule(parent);
    }
  }


  void IdentificationDataConverter::exportStepsAndScoresToMzTab_(
    const ID::AppliedProcessingSteps& steps_and_scores,
    MzTabParameterList& steps_out, map<Size, MzTabDouble>& scores_out,
    map<ID::ScoreTypeRef, Size>& score_map)
  {
    vector<MzTabParameter> search_engines;
    for (const ID::AppliedProcessingStep& applied : steps_and_scores)
    {
      MzTabParameter param;
      if (applied.processing_step_opt)
      {
        ID::ProcessingSoftwareRef sw_ref =
          (*applied.processing_step_opt)->software_ref;
        param.setName(sw_ref->getName());
        param.setValue(sw_ref->getVersion());
      }
      else
      {
        param.setName("unknown");
      }
      search_engines.push_back(param);

      for (const pair<ID::ScoreTypeRef, double>& score_pair :
             applied.getScoresInOrder())
      {
        if (!score_map.count(score_pair.first)) // new score type
        {
          score_map.insert(make_pair(score_pair.first, score_map.size() + 1));
        }
        Size index = score_map[score_pair.first];
        scores_out[index].set(score_pair.second);
      }
    }
    steps_out.set(search_engines);
  }


  void IdentificationDataConverter::addMzTabSEScores_(
    const map<ID::ScoreTypeRef, Size>& scores,
    map<Size, MzTabParameter>& output)
  {
    for (const pair<const ID::ScoreTypeRef, Size>& score_pair : scores)
    {
      const ID::ScoreType& score_type = *score_pair.first;
      MzTabParameter param;
      param.setName(score_type.cv_term.getName());
      param.setAccession(score_type.cv_term.getAccession());
      param.setCVLabel(score_type.cv_term.getCVIdentifierRef());
      output[score_pair.second] = param;
    }
  }


  void IdentificationDataConverter::addMzTabMoleculeParentContext_(
    const ID::MoleculeParentMatch& match,
    MzTabOligonucleotideSectionRow& row)
  {
    if (match.left_neighbor == String(ID::MoleculeParentMatch::LEFT_TERMINUS))
    {
      row.pre.set("-");
    }
    else if (match.left_neighbor !=
             String(ID::MoleculeParentMatch::UNKNOWN_NEIGHBOR))
    {
      row.pre.set(match.left_neighbor);
    }
    if (match.right_neighbor == String(ID::MoleculeParentMatch::RIGHT_TERMINUS))
    {
      row.post.set("-");
    }
    else if (match.right_neighbor !=
             String(ID::MoleculeParentMatch::UNKNOWN_NEIGHBOR))
    {
      row.post.set(match.right_neighbor);
    }
    if (match.start_pos != ID::MoleculeParentMatch::UNKNOWN_POSITION)
    {
      row.start.set(match.start_pos + 1); // mzTab counts from 1
    }
    if (match.end_pos != ID::MoleculeParentMatch::UNKNOWN_POSITION)
    {
      row.end.set(match.end_pos + 1); // mzTab counts from 1
    }
  }


  void IdentificationDataConverter::addMzTabMoleculeParentContext_(
    const ID::MoleculeParentMatch& /* match */,
    MzTabPeptideSectionRow& /* row */)
  {
    // nothing to do here
  }


  ID::SearchParamRef
  IdentificationDataConverter::importDBSearchParameters_(
    const ProteinIdentification::SearchParameters& pisp,
    IdentificationData& id_data)
  {
    ID::DBSearchParam dbsp;
    dbsp.molecule_type = ID::MoleculeType::PROTEIN;
    dbsp.mass_type = ID::MassType(pisp.mass_type);
    dbsp.database = pisp.db;
    dbsp.database_version = pisp.db_version;
    dbsp.taxonomy = pisp.taxonomy;
    vector<Int> charges;
    try
    {
      charges = ListUtils::create<Int>(pisp.charges);
    }
    catch (Exception::ConversionError& /*e*/) { // X! Tandem notation, e.g. "+1-+4"?
      charges = ListUtils::create<Int>(pisp.charges, '-');
      if ((charges.size() == 2) && (charges[0] < charges[1]))
      {
        for (Int z = charges[0] + 1; z < charges[1]; ++z)
        {
          charges.push_back(z);
        }
        sort(charges.begin(), charges.end());
      }
    }
    dbsp.charges.insert(charges.begin(), charges.end());
    dbsp.fixed_mods.insert(pisp.fixed_modifications.begin(),
                           pisp.fixed_modifications.end());
    dbsp.variable_mods.insert(pisp.variable_modifications.begin(),
                              pisp.variable_modifications.end());
    dbsp.precursor_mass_tolerance = pisp.precursor_mass_tolerance;
    dbsp.fragment_mass_tolerance = pisp.fragment_mass_tolerance;
    dbsp.precursor_tolerance_ppm = pisp.precursor_mass_tolerance_ppm;
    dbsp.fragment_tolerance_ppm = pisp.fragment_mass_tolerance_ppm;
    const String& enzyme_name = pisp.digestion_enzyme.getName();
    if (ProteaseDB::getInstance()->hasEnzyme(enzyme_name))
    {
      dbsp.digestion_enzyme = ProteaseDB::getInstance()->getEnzyme(enzyme_name);
    }
    dbsp.missed_cleavages = pisp.missed_cleavages;
    dbsp.enzyme_term_specificity = pisp.enzyme_term_specificity;
    static_cast<MetaInfoInterface&>(dbsp) = pisp;

    return id_data.registerDBSearchParam(dbsp);
  }


  ProteinIdentification::SearchParameters
  IdentificationDataConverter::exportDBSearchParameters_(
    ID::SearchParamRef ref)
  {
    const ID::DBSearchParam& dbsp = *ref;
    ProteinIdentification::SearchParameters pisp;
    pisp.mass_type = ProteinIdentification::PeakMassType(dbsp.mass_type);
    pisp.db = dbsp.database;
    pisp.db_version = dbsp.database_version;
    pisp.taxonomy = dbsp.taxonomy;
    pisp.charges = ListUtils::concatenate(dbsp.charges, ", ");
    pisp.fixed_modifications.insert(pisp.fixed_modifications.end(),
                                    dbsp.fixed_mods.begin(),
                                    dbsp.fixed_mods.end());
    pisp.variable_modifications.insert(pisp.variable_modifications.end(),
                                       dbsp.variable_mods.begin(),
                                       dbsp.variable_mods.end());
    pisp.precursor_mass_tolerance = dbsp.precursor_mass_tolerance;
    pisp.fragment_mass_tolerance = dbsp.fragment_mass_tolerance;
    pisp.precursor_mass_tolerance_ppm = dbsp.precursor_tolerance_ppm;
    pisp.fragment_mass_tolerance_ppm = dbsp.fragment_tolerance_ppm;
    if (dbsp.digestion_enzyme &&
        (dbsp.molecule_type == ID::MoleculeType::PROTEIN))
    {
      pisp.digestion_enzyme =
        *(static_cast<const DigestionEnzymeProtein*>(dbsp.digestion_enzyme));
    }
    else
    {
      pisp.digestion_enzyme = DigestionEnzymeProtein("unknown_enzyme", "");
    }
    pisp.missed_cleavages = dbsp.missed_cleavages;
    static_cast<MetaInfoInterface&>(pisp) = dbsp;

    return pisp;
  }


  void IdentificationDataConverter::exportMSRunInformation_(
    ID::ProcessingStepRef step_ref, ProteinIdentification& protein)
  {
    // are input files mzMLs?
    // @TODO: what if there's a mix of mzMLs and other files?
    bool mzml_inputs = false;
    vector<String> mzml_files;
    for (ID::InputFileRef input_ref : step_ref->input_file_refs)
    {
      FileTypes::Type type = FileHandler::getTypeByFileName(*input_ref);
      if (type == FileTypes::MZML)
      {
        mzml_inputs = true;
        mzml_files.push_back(*input_ref);
      }
      else
      {
        mzml_inputs = false;
        break;
      }
    }
    if (mzml_inputs)
    {
      protein.setPrimaryMSRunPath(mzml_files);
      // also store raw files (or equivalent):
      protein.setPrimaryMSRunPath(step_ref->primary_files, true);
      return;
    }
    // alternatively, are the primary files mzMLs?
    bool mzml_primaries = false;
    for (const String& file : step_ref->primary_files)
    {
      FileTypes::Type type = FileHandler::getTypeByFileName(file);
      if (type == FileTypes::MZML)
      {
        mzml_primaries = true;
      }
      else
      {
        mzml_primaries = false;
        break;
      }
    }
    // store as (raw) primary files depending on file type:
    protein.setPrimaryMSRunPath(step_ref->primary_files, !mzml_primaries);
  }

} // end namespace OpenMS
