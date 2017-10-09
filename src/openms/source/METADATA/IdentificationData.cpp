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
#include <OpenMS/CHEMISTRY/ProteaseDB.h>

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

      SearchParamsKey search_key =
        importDBSearchParameters_(prot_it->getSearchParameters());
      db_search_steps.insert(make_pair(step_key, search_key));

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
          meta.scores.push_back(make_pair(score_key, hit_it->getScore()));
          meta.processing_steps.push_back(step_key);
          static_cast<MetaInfoInterface&>(meta) = *hit_it;
          parent_meta_data.insert(make_pair(result.first, meta));
        }
        else // previously seen protein
        {
          ParentMetaData& meta = parent_meta_data.at(result.first);
          // this won't overwrite:
          meta.scores.push_back(make_pair(score_key, hit_it->getScore()));
          meta.processing_steps.push_back(step_key);
        }
      }
    }

    Size unknown_query_counter = 1;
    for (vector<PeptideIdentification>::const_iterator pep_it =
           peptides.begin(); pep_it != peptides.end(); ++pep_it)
    {
      const String& id = pep_it->getIdentifier();
      ProcessingStepKey step_key = id_to_step[id];
      const DataProcessingStep& step = processing_steps.at(step_key);
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
        if (pep_it->hasRT() && pep_it->hasMZ())
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
        if (result.second)
        {
          if (!hit_it->getPeptideEvidences().empty())
          {
            parent_evidence[result.first] = hit_it->getPeptideEvidences();
          }
          IdentifiedMetaData meta;
          meta.molecule_type = MT_PROTEIN;
          identified_meta_data.insert(make_pair(result.first, meta));
        }
        // @TODO: check/merge evidences if peptide already exists?
        identified_meta_data[result.first].processing_steps.push_back(step_key);

        pair<DataQueryKey, IdentifiedMoleculeKey> psm_key =
          make_pair(query_key, result.first);
        MatchMap::iterator pos = matches.find(psm_key);
        if (pos == matches.end()) // new PSM
        {
          MatchMetaData match;
          match.rank = hit_it->getRank();
          match.charge = hit_it->getCharge();
          match.peak_annotations = hit_it->getPeakAnnotations();
          static_cast<MetaInfoInterface&>(match) = *hit_it;
          pos = matches.insert(make_pair(psm_key, match)).first;
        }
        pos->second.scores.push_back(make_pair(score_key, hit_it->getScore()));
        pos->second.processing_steps.push_back(step_key);
      }
    }
  }


  void IdentificationData::exportIDs(vector<ProteinIdentification>& proteins,
                                     vector<PeptideIdentification>& peptides)
    const
  {
    proteins.clear();
    peptides.clear();

    // "DataQuery" roughly corresponds to "PeptideIdentification",
    // "DataProcessingStep" roughly corresponds to "ProteinIdentification":
    map<pair<DataQueryKey, ProcessingStepKey>,
        pair<vector<PeptideHit>, ScoreTypeKey>> psm_data;
    // we only export peptides and proteins, so start by getting the PSMs:
    for (MatchMap::const_iterator match_it = matches.begin();
         match_it != matches.end(); ++match_it)
    {
      DataQueryKey query_key = match_it->first.first;
      IdentifiedMoleculeKey molecule_key = match_it->first.second;
      const MatchMetaData& match = match_it->second;
      const IdentifiedMetaData& meta =
        identified_meta_data.at(molecule_key);
      if (meta.molecule_type != MT_PROTEIN) continue;
      PeptideHit hit;
      hit.setSequence(identified_peptides.left.at(molecule_key));
      hit.setCharge(match.charge);
      hit.setRank(match.rank);
      hit.setPeakAnnotations(match.peak_annotations);
      EvidenceMap::const_iterator pos =
        parent_evidence.find(molecule_key);
      if (pos != parent_evidence.end())
      {
        hit.setPeptideEvidences(pos->second);
      }
      static_cast<MetaInfoInterface&>(hit) = match;
      // find all steps that assigned a score:
      for (vector<ProcessingStepKey>::const_iterator step_it =
             match.processing_steps.begin(); step_it !=
             match.processing_steps.end(); ++step_it)
      {
        const DataProcessingStep& step = processing_steps.at(*step_it);
        // give priority to "later" scores:
        for (ScoreList::const_reverse_iterator score_it = match.scores.rbegin();
             score_it != match.scores.rend(); ++score_it)
        {
          const ScoreType& score = score_types.left.at(score_it->first);
          if (score.params_key == step.params_key)
          {
            hit.setScore(score_it->second);
            pair<DataQueryKey, ProcessingStepKey> key =
              make_pair(query_key, *step_it);
            psm_data[key].first.push_back(hit);
            psm_data[key].second = score_it->first;
            break;
          }
        }
      }
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
        // give priority to "later" scores:
        for (ScoreList::const_reverse_iterator score_it = meta.scores.rbegin();
             score_it != meta.scores.rend(); ++score_it)
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
        const_iterator pd_pos = prot_data.find(*step_it);
      if (pd_pos != prot_data.end())
      {
        protein.setHits(pd_pos->second.first);
        const ScoreType& score_type = score_types.left.at(pd_pos->
                                                          second.second);
        protein.setScoreType(score_type.name);
      }
      unordered_map<ProcessingStepKey, SearchParamsKey>::const_iterator ss_pos =
        db_search_steps.find(*step_it);
      if (ss_pos != db_search_steps.end())
      {
        protein.setSearchParameters(exportDBSearchParameters_(ss_pos->second));
      }

      proteins.push_back(protein);
    }
  }


  double IdentificationData::findScore_(ScoreTypeKey key,
                                        const ScoreList& scores)
  {
    // give priority to "later" scores in the list:
    for (ScoreList::const_reverse_iterator it = scores.rbegin();
         it != scores.rend(); ++it)
    {
      if (it->first == key) return it->second;
    }
    return numeric_limits<double>::quiet_NaN(); // or throw an exception?
  }


  IdentificationData::SearchParamsKey IdentificationData::importDBSearchParameters_(
      const ProteinIdentification::SearchParameters& pisp)
  {
    DBSearchParameters dbsp;
    dbsp.molecule_type = MT_PROTEIN;
    dbsp.peak_mass_type = pisp.mass_type;
    dbsp.database = pisp.db;
    dbsp.database_version = pisp.db_version;
    dbsp.taxonomy = pisp.taxonomy;
    vector<Int> charges = ListUtils::create<Int>(pisp.charges);
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
    static_cast<MetaInfoInterface&>(dbsp) = pisp;

    return insertIntoBimap_(dbsp, db_search_params).first;
  }


  ProteinIdentification::SearchParameters
  IdentificationData::exportDBSearchParameters_(SearchParamsKey key) const
  {
    const DBSearchParameters& dbsp = db_search_params.left.at(key);
    if (dbsp.molecule_type != MT_PROTEIN)
    {
      String msg = "only proteomics search parameters can be exported";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    ProteinIdentification::SearchParameters pisp;
    pisp.mass_type = dbsp.peak_mass_type;
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
    if (dbsp.digestion_enzyme)
    {
      pisp.digestion_enzyme = *(static_cast<const DigestionEnzymeProtein*>(
                                  dbsp.digestion_enzyme));
    }
    else
    {
      pisp.digestion_enzyme = DigestionEnzymeProtein("unknown_enzyme", "");
    }
    pisp.missed_cleavages = dbsp.missed_cleavages;
    static_cast<MetaInfoInterface&>(pisp) = dbsp;

    return pisp;
  }

} // end namespace OpenMS
