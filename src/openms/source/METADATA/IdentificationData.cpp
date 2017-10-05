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

#include <OpenMS/METADATA/IdentificationData.h>

using namespace std;

namespace OpenMS
{
  void IdentificationData::importIDs(
    const vector<ProteinIdentification>& proteins,
    const vector<PeptideIdentification>& peptides)
  {
    map<String, ProcessingStepKey> id_to_step;
    for (vector<ProteinIdentification>::const_iterator prot_it =
           proteins.begin(); prot_it != proteins.end(); ++prot_it)
    {
      DataProcessingParameters params;
      params.tool.setName(prot_it->getSearchEngine());
      params.tool.setVersion(prot_it->getSearchEngineVersion());
      ProcessingParamsKey params_key =
        insertIntoBimap_(params, processing_params).first;

      ScoreType score_type;
      score_type.name = prot_it->getScoreType();
      score_type.params_key = params_key;
      score_type.higher_better = prot_it->isHigherScoreBetter();
      ScoreTypeKey score_key =
        insertIntoBimap_(score_type, score_types).first;

      DataProcessingStep step;
      prot_it->getPrimaryMSRunPath(step.ms_data_path);
      for (vector<String>::const_iterator path_it = step.ms_data_path.begin();
           path_it != step.ms_data_path.end(); ++ path_it)
      {
        InputFileKey file_key = insertIntoBimap_(*path_it, input_files).first;
        step.input_files.push_back(file_key);
      }
      step.date_time = prot_it->getDateTime();
      step.params_key = params_key;
      // no uniqueness check necessary here, as we don't expect duplicates:
      ProcessingStepKey step_key =
        ProcessingStepKey(UniqueIdGenerator::getUniqueId());
      processing_steps.insert(make_pair(step_key, step));
      id_to_step[prot_it->getIdentifier()] = step_key;

      for (vector<ProteinHit>::const_iterator hit_it =
             prot_it->getHits().begin(); hit_it != prot_it->getHits().end();
           ++hit_it)
      {
        pair<ParentMoleculeKey, bool> result =
          insertIntoBimap_(hit_it->getAccession(), parent_molecules);
        if (result.second) // new protein
        {
          ParentMetaData meta;
          meta.molecule_type = MT_PROTEIN;
          meta.sequence = hit_it->getSequence();
          meta.description = hit_it->getDescription();
          meta.coverage = hit_it->getCoverage();
          meta.scores.insert(make_pair(score_key, hit_it->getScore()));
          meta.processing_steps.push_back(step_key);
          static_cast<MetaInfoInterface&>(meta) = *hit_it;
          parent_meta_data.insert(make_pair(result.first, meta));
        }
        else // previously seen protein
        {
          ParentMetaData& meta = parent_meta_data.at(result.first);
          // this won't overwrite:
          meta.scores.insert(make_pair(score_key, hit_it->getScore()));
          meta.processing_steps.push_back(step_key);
        }
      }
    }

    for (vector<PeptideIdentification>::const_iterator pep_it =
           peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      const String& id = pep_it->getIdentifier();
      const DataProcessingStep& step = processing_steps.at(id_to_step[id]);
      DataQuery query;
      if (!step.input_files.empty())
      {
        query.input_file_key = step.input_files[0];
      }
      else
      {
        String file = "UNKNOWN_INPUT_FILE_" + id;
        InputFileKey file_key = insertIntoBimap_(file, input_files).first;
        query.input_file_key = file_key;
      }
      query.rt = pep_it->getRT();
      query.mz = pep_it->getMZ();
      static_cast<MetaInfoInterface&>(query) = *pep_it;
      if (pep_it->metaValueExists("spectrum_reference"))
      {
        query.data_id = pep_it->getMetaValue("spectrum_reference");
        query.removeMetaValue("spectrum_reference");
      }
      else
      {
        query.data_id = String("RT=") + String(float(query.rt)) + "_MZ=" +
          String(float(query.mz));
      }
      DataQueryKey query_key = insertIntoBimap_(query, data_queries).first;

      ScoreType score_type;
      score_type.name = pep_it->getScoreType();
      score_type.params_key = step.params_key;
      score_type.higher_better = pep_it->isHigherScoreBetter();
      ScoreTypeKey score_key =
        insertIntoBimap_(score_type, score_types).first;

      for (vector<PeptideHit>::const_iterator hit_it =
             pep_it->getHits().begin(); hit_it != pep_it->getHits().end();
           ++hit_it)
      {
        const AASequence& seq = hit_it->getSequence();
        pair<IdentifiedMoleculeKey, bool> result =
          insertIntoBimap_(seq, identified_peptides);
        if (result.second && !hit_it->getPeptideEvidences().empty())
        {
          identified_meta_data[result.first] = hit_it->getPeptideEvidences();
        }
        // @TODO: check/merge evidences if peptide already exists?

        MatchMetaData match;
        match.molecule_type = MT_PROTEIN;
        match.molecule_key = result.first;
        match.scores.insert(make_pair(score_key, hit_it->getScore()));
        match.rank = hit_it->getRank();
        match.charge = hit_it->getCharge();
        match.processing_steps.push_back(id_to_step[id]);
        static_cast<MetaInfoInterface&>(match) = *hit_it;

        matches.insert(make_pair(query_key, match));
      }
    }
  }


  void IdentificationData::exportIDs(vector<ProteinIdentification>& proteins,
                                     vector<PeptideIdentification>& peptides)
    const
  {
    proteins.clear();
    peptides.clear();

    map<pair<DataQueryKey, ProcessingStepKey>,
        pair<vector<PeptideHit>, ScoreTypeKey>> psm_data;
    // we only export peptides and proteins, so start by getting the PSMs:
    for (MatchMultimap::const_iterator match_it = matches.begin();
         match_it != matches.end(); ++match_it)
    {
      const MatchMetaData& match = match_it->second;
      if (match.molecule_type != MT_PROTEIN) continue;
      // "DataQuery" roughly corresponds to "PeptideIdentification":
      DataQueryKey query_key = match_it->first;
      IdentifiedMoleculeKey molecule_key = match.molecule_key;
      PeptideHit hit;
      hit.setSequence(identified_peptides.left.at(molecule_key));
      hit.setCharge(match.charge);
      hit.setRank(match.rank);
      // "DataProcessingStep" roughly corresponds to "ProteinIdentification";
      // find the last step that assigned a score:
      ProcessingStepKey step_key = match.processing_steps.back(); // most recent
      bool done = false;
      for (vector<ProcessingStepKey>::const_reverse_iterator step_it =
             match.processing_steps.rbegin();
           !done && (step_it != match.processing_steps.rend()); ++step_it)
      {
        const DataProcessingStep& step = processing_steps.at(*step_it);
        for (unordered_map<ScoreTypeKey, double>::const_iterator score_it =
               match.scores.begin(); score_it != match.scores.end();
             ++score_it)
        {
          const ScoreType& score = score_types.left.at(score_it->first);
          if (score.params_key == step.params_key)
          {
            hit.setScore(score_it->second);
            step_key = *step_it;
            psm_data[make_pair(query_key, step_key)].second = score_it->first;
            done = true;
            break;
          }
        }
      }
      EvidenceMap::const_iterator pos =
        identified_meta_data.find(molecule_key);
      if (pos != identified_meta_data.end())
      {
        hit.setPeptideEvidences(pos->second);
      }
      static_cast<MetaInfoInterface&>(hit) = match;

      psm_data[make_pair(query_key, step_key)].first.push_back(hit);
    }

    set<ProcessingStepKey> steps;
    for (auto psm_it = psm_data.begin(); psm_it != psm_data.end(); ++psm_it)
    {
      const DataQuery& query = data_queries.left.at(psm_it->first.first);
      PeptideIdentification peptide;
      static_cast<MetaInfoInterface&>(peptide) = query;
      peptide.setRT(query.rt);
      peptide.setMZ(query.mz);
      peptide.setMetaValue("spectrum_reference", query.data_id);
      peptide.setHits(psm_it->second.first);
      const ScoreType& score_type = score_types.left.at(psm_it->second.second);
      peptide.setScoreType(score_type.name);
      peptide.sortByRank();
      peptide.setIdentifier(String(psm_it->first.second));
      peptides.push_back(peptide);
      steps.insert(psm_it->first.second);
    }

    map<ProcessingStepKey, pair<vector<ProteinHit>, ScoreTypeKey>> prot_data;
    for (unordered_map<ParentMoleculeKey, ParentMetaData>::const_iterator
           par_it = parent_meta_data.begin(); par_it != parent_meta_data.end();
         ++par_it)
    {
      if (par_it->second.molecule_type != MT_PROTEIN) continue;
      const ParentMetaData& meta = par_it->second;
      ProteinHit hit;
      hit.setSequence(meta.sequence);
      hit.setDescription(meta.description);
      hit.setCoverage(meta.coverage);
      hit.setAccession(parent_molecules.left.at(par_it->first));
      static_cast<MetaInfoInterface&>(hit) = meta;
      // find all steps that assigned a score:
      for (vector<ProcessingStepKey>::const_iterator step_it =
             meta.processing_steps.begin(); step_it !=
             meta.processing_steps.end(); ++step_it)
      {
        const DataProcessingStep& step = processing_steps.at(*step_it);
        for (unordered_map<ScoreTypeKey, double>::const_iterator score_it =
               meta.scores.begin(); score_it != meta.scores.end();
             ++score_it)
        {
          const ScoreType& score = score_types.left.at(score_it->first);
          if (score.params_key == step.params_key)
          {
            hit.setScore(score_it->second);
            prot_data[*step_it].first.push_back(hit);
            prot_data[*step_it].second = score_it->first;
            steps.insert(*step_it);
            break;
          }
        }
      }
    }

    for (set<ProcessingStepKey>::const_iterator step_it = steps.begin();
         step_it != steps.end(); ++step_it)
    {
      ProteinIdentification protein;
      protein.setIdentifier(String(*step_it));
      const DataProcessingStep& step = processing_steps.at(*step_it);
      protein.setDateTime(step.date_time);
      protein.setPrimaryMSRunPath(step.ms_data_path);
      const DataProcessingParameters& params =
        processing_params.left.at(step.params_key);
      protein.setSearchEngine(params.tool.getName());
      protein.setSearchEngineVersion(params.tool.getVersion());
      map<ProcessingStepKey, pair<vector<ProteinHit>, ScoreTypeKey>>::
        const_iterator pos = prot_data.find(*step_it);
      if (pos != prot_data.end())
      {
        protein.setHits(pos->second.first);
        const ScoreType& score_type = score_types.left.at(pos->second.second);
        protein.setScoreType(score_type.name);
      }
      proteins.push_back(protein);
    }
  }
}
