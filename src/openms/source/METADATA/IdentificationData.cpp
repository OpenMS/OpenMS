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
  const Size IdentificationData::MoleculeParentMatch::UNKNOWN_POSITION =
    Size(-1);
  const char IdentificationData::MoleculeParentMatch::UNKNOWN_NEIGHBOR = 'X';
  const char IdentificationData::MoleculeParentMatch::LEFT_TERMINUS = '[';
  const char IdentificationData::MoleculeParentMatch::RIGHT_TERMINUS = ']';


  void IdentificationData::importIDs(
    const vector<ProteinIdentification>& proteins,
    const vector<PeptideIdentification>& peptides)
  {
    map<String, ProcessingStepKey> id_to_step;

    // ProteinIdentification:
    for (const ProteinIdentification& prot : proteins)
    {
      DataProcessingSoftware software(prot.getSearchEngine(),
                                      prot.getSearchEngineVersion());
      ProcessingSoftwareKey software_key =
        registerDataProcessingSoftware(software).first;

      ScoreType score_type(prot.getScoreType(),
                           prot.isHigherScoreBetter(), software_key);
      ScoreTypeKey score_key = registerScoreType(score_type).first;

      SearchParamsKey search_key =
        importDBSearchParameters_(prot.getSearchParameters());

      DataProcessingStep step;
      step.software_key = software_key;
      prot.getPrimaryMSRunPath(step.primary_files);
      for (const String& path : step.primary_files)
      {
        InputFileKey file_key = registerInputFile(path).first;
        step.input_files.push_back(file_key);
      }
      step.date_time = prot.getDateTime();
      ProcessingStepKey step_key =
        registerDataProcessingStep(step, search_key).first;
      id_to_step[prot.getIdentifier()] = step_key;
      setCurrentProcessingStep(step_key);

      // ProteinHit:
      for (const ProteinHit& hit : prot.getHits())
      {
        pair<ParentMoleculeKey, bool> result =
          registerParentMolecule(hit.getAccession());
        ParentMetaData& meta = parent_meta_data.at(result.first);
        if (result.second) // new protein
        {
          meta.sequence = hit.getSequence();
          meta.description = hit.getDescription();
          meta.coverage = hit.getCoverage();
          static_cast<MetaInfoInterface&>(meta) = hit;
        }
        // add this information also for previously seen proteins:
        meta.scores.push_back(make_pair(score_key, hit.getScore()));
      }
      clearCurrentProcessingStep();
    }

    // PeptideIdentification:
    Size unknown_query_counter = 1;
    for (const PeptideIdentification& pep : peptides)
    {
      const String& id = pep.getIdentifier();
      ProcessingStepKey step_key = id_to_step[id];
      const DataProcessingStep& step = processing_steps.left.at(step_key);
      DataQuery query;
      if (!step.input_files.empty())
      {
        // @TODO: what if there's more than one input file?
        query.input_file_key = step.input_files[0];
      }
      else
      {
        String file = "UNKNOWN_INPUT_FILE_" + id;
        InputFileKey file_key = insertIntoBimap_(file, input_files).first;
        query.input_file_key = file_key;
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
      DataQueryKey query_key = registerDataQuery(query).first;

      ScoreType score_type(pep.getScoreType(), pep.isHigherScoreBetter(),
                           step.software_key);
      ScoreTypeKey score_key = registerScoreType(score_type).first;

      // PeptideHit:
      for (const PeptideHit& hit : pep.getHits())
      {
        if (hit.getSequence().empty()) continue;
        pair<IdentifiedMoleculeKey, bool> result =
          registerPeptide(hit.getSequence());
        for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
        {
          const String& accession = evidence.getProteinAccession();
          if (accession.empty()) continue;
          // this won't overwrite if the protein already exists:
          ParentMoleculeKey parent_key =
            registerParentMolecule(accession).first;
          vector<ProcessingStepKey>& steps =
            parent_meta_data.at(parent_key).processing_steps;
          if (find(steps.begin(), steps.end(), step_key) == steps.end())
          {
            steps.push_back(step_key);
          }
          MoleculeParentMatch match(evidence.getStart(), evidence.getEnd(),
                                    evidence.getAABefore(),
                                    evidence.getAAAfter());
          parent_matches[result.first][parent_key].insert(match);
        }
        IdentifiedMetaData& meta = identified_meta_data.at(result.first);
        meta.processing_steps.push_back(step_key);

        QueryMatchKey psm_key = make_pair(result.first, query_key);
        QueryMatchMap::iterator pos = query_matches.find(psm_key);
        if (pos == query_matches.end()) // new PSM
        {
          MoleculeQueryMatch match;
          match.charge = hit.getCharge();
          match.peak_annotations = hit.getPeakAnnotations();
          static_cast<MetaInfoInterface&>(match) = hit;
          pos = query_matches.insert(make_pair(psm_key, match)).first;
        }
        MoleculeQueryMatch& match = pos->second;
        match.processing_steps.push_back(step_key);

        // analysis results from pepXML:
        for (const PeptideHit::PepXMLAnalysisResult& ana_res :
               hit.getAnalysisResults())
        {
          DataProcessingSoftware software;
          software.tool.setName(ana_res.score_type); // e.g. "peptideprophet"
          ProcessingSoftwareKey software_key =
            registerDataProcessingSoftware(software).first;
          DataProcessingStep sub_step;
          sub_step.software_key = software_key;
          if (query.input_file_key != 0)
          {
            sub_step.input_files.push_back(query.input_file_key);
          }
          ProcessingStepKey sub_step_key =
            insertIntoBimap_(sub_step, processing_steps).first;
          match.processing_steps.push_back(sub_step_key);
          for (const pair<String, double>& sub_pair : ana_res.sub_scores)
          {
            ScoreType sub_score;
            sub_score.name = sub_pair.first;
            sub_score.software_key = software_key;
            ScoreTypeKey sub_score_key =
              insertIntoBimap_(sub_score, score_types).first;
            match.scores.push_back(make_pair(sub_score_key, sub_pair.second));
          }
          ScoreType main_score;
          main_score.name = ana_res.score_type + "_probability";
          main_score.higher_better = ana_res.higher_is_better;
          main_score.software_key = software_key;
          ScoreTypeKey main_score_key =
            insertIntoBimap_(main_score, score_types).first;
          match.scores.push_back(make_pair(main_score_key, ana_res.main_score));
        }
        // primary score goes last:
        match.scores.push_back(make_pair(score_key, hit.getScore()));
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
    for (QueryMatchMap::const_iterator match_it = query_matches.begin();
         match_it != query_matches.end(); ++match_it)
    {
      IdentifiedMoleculeKey molecule_key = match_it->first.first;
      DataQueryKey query_key = match_it->first.second;
      const MoleculeQueryMatch& match = match_it->second;
      const IdentifiedMetaData& meta =
        identified_meta_data.at(molecule_key);
      if (meta.molecule_type != MT_PROTEIN) continue;
      PeptideHit hit;
      hit.setSequence(identified_peptides.left.at(molecule_key));
      hit.setCharge(match.charge);
      hit.setPeakAnnotations(match.peak_annotations);
      ParentMatchMap::const_iterator pos = parent_matches.find(molecule_key);
      if (pos != parent_matches.end())
      {
        for (const ParentSubMap::value_type& sub : pos->second)
        {
          ParentMoleculeKey parent_key = sub.first;
          for (const MoleculeParentMatch& match : sub.second)
          {
            PeptideEvidence evidence;
            evidence.setProteinAccession(parent_molecules.left.at(parent_key));
            evidence.setStart(match.start_pos);
            evidence.setEnd(match.end_pos);
            evidence.setAABefore(match.left_neighbor);
            evidence.setAAAfter(match.right_neighbor);
            hit.addPeptideEvidence(evidence);
          }
        }
      }
      static_cast<MetaInfoInterface&>(hit) = match;
      // find all steps that assigned a score:
      for (vector<ProcessingStepKey>::const_iterator step_it =
             match.processing_steps.begin(); step_it !=
             match.processing_steps.end(); ++step_it)
      {
        const DataProcessingStep& step = processing_steps.left.at(*step_it);
        // give priority to "later" scores:
        for (ScoreList::const_reverse_iterator score_it = match.scores.rbegin();
             score_it != match.scores.rend(); ++score_it)
        {
          const ScoreType& score = score_types.left.at(score_it->first);
          if (score.software_key == step.software_key)
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
        const DataProcessingStep& step = processing_steps.left.at(*step_it);
        // give priority to "later" scores:
        for (ScoreList::const_reverse_iterator score_it = meta.scores.rbegin();
             score_it != meta.scores.rend(); ++score_it)
        {
          const ScoreType& score = score_types.left.at(score_it->first);
          if (score.software_key == step.software_key)
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
      const DataProcessingStep& step = processing_steps.left.at(*step_it);
      protein.setDateTime(step.date_time);
      protein.setPrimaryMSRunPath(step.primary_files);
      const DataProcessingSoftware& software =
        processing_software.left.at(step.software_key);
      protein.setSearchEngine(software.tool.getName());
      protein.setSearchEngineVersion(software.tool.getVersion());
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


  IdentificationData::SearchParamsKey
  IdentificationData::importDBSearchParameters_(
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

    return registerDBSearchParameters(dbsp).first;
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


  void IdentificationData::checkScoreTypes_(const ScoreList& scores)
  {
    for (const pair<ScoreTypeKey, double>& score_pair : scores)
    {
      if (!UniqueIdInterface::isValid(score_pair.first) ||
          (score_types.left.count(score_pair.first) == 0))
      {
        String msg = "invalid reference to a score type - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
  }


  void IdentificationData::checkProcessingSteps_(
    const std::vector<ProcessingStepKey>& steps)
  {
    for (ProcessingStepKey step : steps)
    {
      if (!UniqueIdInterface::isValid(step) ||
          (processing_steps.left.count(step) == 0))
      {
        String msg = "invalid reference to a data processing step - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
  }


  bool IdentificationData::addCurrentProcessingStep_(
    vector<ProcessingStepKey>& processing_steps)
  {
    if (UniqueIdInterface::isValid(current_step_key_) &&
        (processing_steps.empty() ||
         (processing_steps.back() != current_step_key_)))
    {
      processing_steps.push_back(current_step_key_);
      return true;
    }
    return false;
  }


  pair<IdentificationData::InputFileKey, bool>
  IdentificationData::registerInputFile(const String& file)
  {
    return insertIntoBimap_(file, input_files);
  }


  pair<IdentificationData::ProcessingSoftwareKey, bool>
  IdentificationData::registerDataProcessingSoftware(
    const DataProcessingSoftware& software)
  {
    return insertIntoBimap_(software, processing_software);
  }


  pair<IdentificationData::SearchParamsKey, bool>
  IdentificationData::registerDBSearchParameters(
    const DBSearchParameters& params)
  {
    // @TODO: any required information that should be checked?
    return insertIntoBimap_(params, db_search_params);
  }


  pair<IdentificationData::ProcessingStepKey, bool>
  IdentificationData::registerDataProcessingStep(const DataProcessingStep& step,
                                                 SearchParamsKey search_key)
  {
    // valid reference to software is required:
    if (!UniqueIdInterface::isValid(step.software_key) ||
        (processing_software.left.count(step.software_key) == 0))
    {
      String msg = "invalid reference to data processing parameters - register those first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    // if given, references to input files must be valid:
    for (InputFileKey key : step.input_files)
    {
      if (!UniqueIdInterface::isValid(key) ||
          (input_files.left.count(key) == 0))
      {
      String msg = "invalid reference to input file - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
      }
    }

    pair<ProcessingStepKey, bool> result =
      insertIntoBimap_(step, processing_steps);
    if (result.second && UniqueIdInterface::isValid(search_key))
    {
      // if given, reference to DB search params must be valid:
      if (db_search_params.left.count(search_key) == 0)
      {
        String msg = "invalid reference to DB search parameters - register those first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      db_search_steps.insert(make_pair(result.first, search_key));
    }
    return result;
  }


  pair<IdentificationData::ScoreTypeKey, bool>
  IdentificationData::registerScoreType(const ScoreType& score)
  {
    // reference to software may be missing, but otherwise must be valid:
    if (UniqueIdInterface::isValid(score.software_key))
    {
      if (processing_software.left.count(score.software_key) == 0)
      {
        String msg = "invalid reference to data processing software - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
    else if (UniqueIdInterface::isValid(current_step_key_))
    {
      // transfer the software key from the current data processing step:
      const DataProcessingStep& step =
        processing_steps.left.at(current_step_key_);
      ScoreType copy(score); // need a copy so we can modify it
      copy.software_key = step.software_key;
      return insertIntoBimap_(copy, score_types);
    }
    return insertIntoBimap_(score, score_types);
  }


  pair<IdentificationData::DataQueryKey, bool>
  IdentificationData::registerDataQuery(const DataQuery& query)
  {
    // reference to spectrum or feature is required:
    if (query.data_id.empty())
    {
      String msg = "missing identifier in data query";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    // reference to input file may be missing, but otherwise must be valid:
    if (UniqueIdInterface::isValid(query.input_file_key) &&
        (input_files.left.count(query.input_file_key) == 0))
    {
      String msg = "invalid reference to input file - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    return insertIntoBimap_(query, data_queries);
  }


  pair<IdentificationData::IdentifiedMoleculeKey, bool>
  IdentificationData::registerPeptide(const AASequence& seq,
                                      IdentifiedMetaData meta_data)
  {
    if (seq.empty())
    {
      String msg = "missing sequence for peptide";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    if (meta_data.molecule_type != MT_PROTEIN)
    {
      String msg = "invalid molecule type setting for peptide";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    pair<IdentifiedMoleculeKey, bool> result =
      insertIntoBimap_(seq, identified_peptides);
    if (result.second)
    {
      checkScoreTypes_(meta_data.scores);
      checkProcessingSteps_(meta_data.processing_steps);
      // @TODO: this prevents passing "meta_data" by const ref. - any
      // workaround?
      addCurrentProcessingStep_(meta_data.processing_steps);
      identified_meta_data.insert(make_pair(result.first, meta_data));
    }
    return result;
  }


  pair<IdentificationData::IdentifiedMoleculeKey, bool>
  IdentificationData::registerCompound(const String& id,
                                       const CompoundMetaData& compound_meta,
                                       IdentifiedMetaData id_meta)
  {
    if (id.empty())
    {
      String msg = "missing identifier for compound";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    if (id_meta.molecule_type != MT_COMPOUND)
    {
      String msg = "invalid molecule type setting for compound";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    pair<IdentifiedMoleculeKey, bool> result =
      insertIntoBimap_(id, identified_compounds);
    if (result.second)
    {
      checkScoreTypes_(id_meta.scores);
      checkProcessingSteps_(id_meta.processing_steps);
      // @TODO: this prevents passing "id_meta" by const ref. - any workaround?
      addCurrentProcessingStep_(id_meta.processing_steps);
      identified_meta_data.insert(make_pair(result.first, id_meta));
      // @TODO: insert "compound_meta" even if it's empty?
      compound_meta_data.insert(make_pair(result.first, compound_meta));
    }
    return result;
  }


  pair<IdentificationData::IdentifiedMoleculeKey, bool>
  IdentificationData::registerOligo(const NASequence& seq,
                                    IdentifiedMetaData meta_data)
  {
    if (seq.empty())
    {
      String msg = "missing sequence for oligonucleotide";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    if (meta_data.molecule_type != MT_RNA)
    {
      String msg = "invalid molecule type setting for oligonucleotide";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    pair<IdentifiedMoleculeKey, bool> result =
      insertIntoBimap_(seq, identified_oligos);
    if (result.second)
    {
      checkScoreTypes_(meta_data.scores);
      checkProcessingSteps_(meta_data.processing_steps);
      // @TODO: this prevents passing "meta_data" by const ref. - any
      // workaround?
      addCurrentProcessingStep_(meta_data.processing_steps);
      identified_meta_data.insert(make_pair(result.first, meta_data));
    }
    return result;
  }


  pair<IdentificationData::ParentMoleculeKey, bool>
  IdentificationData::registerParentMolecule(const String& accession,
                                             ParentMetaData meta_data)
  {
    if (accession.empty())
    {
      String msg = "missing accession for parent molecule";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    pair<ParentMoleculeKey, bool> result =
      insertIntoBimap_(accession, parent_molecules);
    if (result.second)
    {
      checkScoreTypes_(meta_data.scores);
      checkProcessingSteps_(meta_data.processing_steps);
      addCurrentProcessingStep_(meta_data.processing_steps);
      parent_meta_data.insert(make_pair(result.first, meta_data));
    }
    return result;
  }


  bool IdentificationData::addMoleculeParentMatch(
    IdentifiedMoleculeKey molecule_key, ParentMoleculeKey parent_key,
    const MoleculeParentMatch& meta_data)
  {
    if (!UniqueIdInterface::isValid(molecule_key) ||
        // don't know which "identified_[type]" to check, so check meta data:
        (identified_meta_data.count(molecule_key) == 0))
    {
      String msg = "invalid reference to an identified molecule - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    if (!UniqueIdInterface::isValid(parent_key) ||
        (parent_molecules.left.count(parent_key) == 0))
    {
      String msg = "invalid reference to a parent molecule - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    const ParentMetaData& parent_meta = parent_meta_data.at(parent_key);
    const IdentifiedMetaData& id_meta = identified_meta_data.at(molecule_key);
    if (parent_meta.molecule_type != id_meta.molecule_type)
    {
      String msg = "molecule types of parent and child molecule don't match";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    return parent_matches[molecule_key][parent_key].insert(meta_data).second;
  }


  bool IdentificationData::addMoleculeQueryMatch(
    IdentifiedMoleculeKey molecule_key, DataQueryKey query_key,
    MoleculeQueryMatch meta_data)
  {
    if (!UniqueIdInterface::isValid(molecule_key) ||
        // don't know which "identified_[type]" to check, so check meta data:
        (identified_meta_data.count(molecule_key) == 0))
    {
      String msg = "invalid reference to an identified molecule - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    if (!UniqueIdInterface::isValid(query_key) ||
        (data_queries.left.count(query_key) == 0))
    {
      String msg = "invalid reference to a data query - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    checkScoreTypes_(meta_data.scores);
    checkProcessingSteps_(meta_data.processing_steps);
    // @TODO: this prevents passing "meta_data" by const ref. - any workaround?
    addCurrentProcessingStep_(meta_data.processing_steps);
    // @TODO: disallow charge zero in "meta_data"?

    QueryMatchKey match_key = make_pair(molecule_key, query_key);
    return query_matches.insert(make_pair(match_key, meta_data)).second;
  }


  void IdentificationData::setCurrentProcessingStep(ProcessingStepKey step_key)
  {
    if (!UniqueIdInterface::isValid(step_key) ||
        (processing_steps.left.count(step_key) == 0))
    {
      String msg = "invalid reference to a processing step - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    current_step_key_ = step_key;
  }


  void IdentificationData::clearCurrentProcessingStep()
  {
    current_step_key_ = 0;
  }

} // end namespace OpenMS
