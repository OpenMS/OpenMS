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
    map<String, ProcessingStepRef> id_to_step;

    // ProteinIdentification:
    for (const ProteinIdentification& prot : proteins)
    {
      Software software(prot.getSearchEngine(), prot.getSearchEngineVersion());
      ProcessingSoftwareRef software_ref =
        registerDataProcessingSoftware(software);

      ScoreType score_type(prot.getScoreType(),
                           prot.isHigherScoreBetter(), software_ref);
      ScoreTypeRef score_ref = registerScoreType(score_type);

      SearchParamRef search_ref =
        importDBSearchParameters_(prot.getSearchParameters());

      DataProcessingStep step(software_ref);
      prot.getPrimaryMSRunPath(step.primary_files);
      for (const String& path : step.primary_files)
      {
        InputFileRef file_ref = registerInputFile(path);
        step.input_file_refs.push_back(file_ref);
      }
      step.date_time = prot.getDateTime();
      ProcessingStepRef step_ref = registerDataProcessingStep(step, search_ref);
      id_to_step[prot.getIdentifier()] = step_ref;
      setCurrentProcessingStep(step_ref);

      // ProteinHit:
      for (const ProteinHit& hit : prot.getHits())
      {
        ParentMolecule parent(hit.getAccession());
        parent.sequence = hit.getSequence();
        parent.description = hit.getDescription();
        parent.coverage = hit.getCoverage();
        static_cast<MetaInfoInterface&>(parent) = hit;
        parent.scores.push_back(make_pair(score_ref, hit.getScore()));
        registerParentMolecule(parent);
      }
      clearCurrentProcessingStep();
    }

    // PeptideIdentification:
    Size unknown_query_counter = 1;
    for (const PeptideIdentification& pep : peptides)
    {
      const String& id = pep.getIdentifier();
      ProcessingStepRef step_ref = id_to_step[id];
      DataQuery query(""); // fill in "data_id" later
      if (!step_ref->input_file_refs.empty())
      {
        // @TODO: what if there's more than one input file?
        query.input_file_ref = step_ref->input_file_refs[0];
      }
      else
      {
        String file = "UNKNOWN_INPUT_FILE_" + id;
        InputFileRef file_ref = registerInputFile(file);
        query.input_file_ref = file_ref;
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
      DataQueryRef query_ref = registerDataQuery(query);

      ScoreType score_type(pep.getScoreType(), pep.isHigherScoreBetter(),
                           step_ref->software_ref);
      ScoreTypeRef score_ref = registerScoreType(score_type);

      // PeptideHit:
      for (const PeptideHit& hit : pep.getHits())
      {
        if (hit.getSequence().empty()) continue;
        IdentifiedPeptide peptide(hit.getSequence());
        peptide.processing_step_refs.push_back(step_ref);
        for (const PeptideEvidence& evidence : hit.getPeptideEvidences())
        {
          const String& accession = evidence.getProteinAccession();
          if (accession.empty()) continue;
          ParentMolecule parent(accession);
          parent.processing_step_refs.push_back(step_ref);
          // this will merge information if the protein already exists:
          ParentMoleculeRef parent_ref = registerParentMolecule(parent);
          MoleculeParentMatch match(evidence.getStart(), evidence.getEnd(),
                                    evidence.getAABefore(),
                                    evidence.getAAAfter());
          peptide.parent_matches[parent_ref].insert(match);
        }
        IdentifiedPeptideRef peptide_ref = registerPeptide(peptide);

        MoleculeQueryMatch match(peptide_ref, query_ref);
        match.charge = hit.getCharge();
        static_cast<MetaInfoInterface&>(match) = hit;
        if (!hit.getPeakAnnotations().empty())
        {
          match.peak_annotations[step_ref] = hit.getPeakAnnotations();
        }
        match.processing_step_refs.push_back(step_ref);

        // analysis results from pepXML:
        for (const PeptideHit::PepXMLAnalysisResult& ana_res :
               hit.getAnalysisResults())
        {
          Software software;
          software.setName(ana_res.score_type); // e.g. "peptideprophet"
          ProcessingSoftwareRef software_ref =
            registerDataProcessingSoftware(software);
          DataProcessingStep sub_step(software_ref);
          if (query.input_file_ref != nullptr)
          {
            sub_step.input_file_refs.push_back(query.input_file_ref);
          }
          ProcessingStepRef sub_step_ref = registerDataProcessingStep(sub_step);
          match.processing_step_refs.push_back(sub_step_ref);
          for (const pair<String, double>& sub_pair : ana_res.sub_scores)
          {
            ScoreType sub_score;
            sub_score.name = sub_pair.first;
            sub_score.software_ref = software_ref;
            ScoreTypeRef sub_score_ref = registerScoreType(sub_score);
            match.scores.push_back(make_pair(sub_score_ref, sub_pair.second));
          }
          ScoreType main_score;
          main_score.name = ana_res.score_type + "_probability";
          main_score.higher_better = ana_res.higher_is_better;
          main_score.software_ref = software_ref;
          ScoreTypeRef main_score_ref = registerScoreType(main_score);
          match.scores.push_back(make_pair(main_score_ref, ana_res.main_score));
        }

        // primary score goes last:
        match.scores.push_back(make_pair(score_ref, hit.getScore()));
        registerMoleculeQueryMatch(match);
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
    map<pair<DataQueryRef, ProcessingStepRef>,
        pair<vector<PeptideHit>, ScoreTypeRef>> psm_data;
    // we only export peptides and proteins, so start by getting the PSMs:
    for (const MoleculeQueryMatch& query_match : query_matches)
    {
      if (query_match.getMoleculeType() != MT_PROTEIN) continue;
      IdentifiedPeptideRef peptide_ref = query_match.getIdentifiedPeptideRef();
      DataQueryRef query_ref = query_match.data_query_ref;
      PeptideHit hit;
      hit.setSequence(peptide_ref->sequence);
      hit.setCharge(query_match.charge);
      for (const auto& pair : peptide_ref->parent_matches)
      {
        ParentMoleculeRef parent_ref = pair.first;
        for (const MoleculeParentMatch& parent_match : pair.second)
        {
          PeptideEvidence evidence;
          evidence.setProteinAccession(parent_ref->accession);
          evidence.setStart(parent_match.start_pos);
          evidence.setEnd(parent_match.end_pos);
          evidence.setAABefore(parent_match.left_neighbor);
          evidence.setAAAfter(parent_match.right_neighbor);
          hit.addPeptideEvidence(evidence);
        }
      }
      static_cast<MetaInfoInterface&>(hit) = query_match;
      // find all steps that assigned a score:
      for (ProcessingStepRef step_ref : query_match.processing_step_refs)
      {
        auto pos = query_match.peak_annotations.find(step_ref);
        if (pos != query_match.peak_annotations.end())
        {
          hit.setPeakAnnotations(pos->second);
        }
        // give priority to "later" scores:
        for (ScoreList::const_reverse_iterator score_it =
               query_match.scores.rbegin(); score_it !=
               query_match.scores.rend(); ++score_it)
        {
          ScoreTypeRef score_ref = score_it->first;
          if (score_ref->software_ref == step_ref->software_ref)
          {
            hit.setScore(score_it->second);
            pair<DataQueryRef, ProcessingStepRef> key = make_pair(query_ref,
                                                                  step_ref);
            psm_data[key].first.push_back(hit);
            psm_data[key].second = score_ref;
            break;
          }
        }
      }
    }

    set<ProcessingStepRef> steps;
    for (auto psm_it = psm_data.begin(); psm_it != psm_data.end(); ++psm_it)
    {
      const DataQuery& query = *psm_it->first.first;
      PeptideIdentification peptide;
      static_cast<MetaInfoInterface&>(peptide) = query;
      peptide.setRT(query.rt);
      peptide.setMZ(query.mz);
      peptide.setMetaValue("spectrum_reference", query.data_id);
      peptide.setHits(psm_it->second.first);
      const ScoreType& score_type = *psm_it->second.second;
      peptide.setScoreType(score_type.name);
      peptide.setIdentifier(String(Size(&(*psm_it->first.second))));
      peptides.push_back(peptide);
      steps.insert(psm_it->first.second);
    }

    map<ProcessingStepRef, pair<vector<ProteinHit>, ScoreTypeRef>> prot_data;
    for (const auto& parent : parent_molecules)
    {
      if (parent.molecule_type != MT_PROTEIN) continue;
      ProteinHit hit;
      hit.setAccession(parent.accession);
      hit.setSequence(parent.sequence);
      hit.setDescription(parent.description);
      hit.setCoverage(parent.coverage);
      static_cast<MetaInfoInterface&>(hit) = parent;
      // find all steps that assigned a score:
      for (ProcessingStepRef step_ref : parent.processing_step_refs)
      {
        // give priority to "later" scores:
        for (ScoreList::const_reverse_iterator score_it =
               parent.scores.rbegin(); score_it != parent.scores.rend();
             ++score_it)
        {
          ScoreTypeRef score_ref = score_it->first;
          if (score_ref->software_ref == step_ref->software_ref)
          {
            hit.setScore(score_it->second);
            prot_data[step_ref].first.push_back(hit);
            prot_data[step_ref].second = score_ref;
            steps.insert(step_ref);
            break;
          }
        }
      }
    }

    for (ProcessingStepRef step_ref : steps)
    {
      ProteinIdentification protein;
      protein.setIdentifier(String(Size(step_ref)));
      protein.setDateTime(step_ref->date_time);
      protein.setPrimaryMSRunPath(step_ref->primary_files);
      const Software& software = *step_ref->software_ref;
      protein.setSearchEngine(software.getName());
      protein.setSearchEngineVersion(software.getVersion());
      map<ProcessingStepRef, pair<vector<ProteinHit>, ScoreTypeRef>>::
        const_iterator pd_pos = prot_data.find(step_ref);
      if (pd_pos != prot_data.end())
      {
        protein.setHits(pd_pos->second.first);
        const ScoreType& score_type = *pd_pos->second.second;
        protein.setScoreType(score_type.name);
      }
      map<ProcessingStepRef, SearchParamRef>::const_iterator ss_pos =
        db_search_steps.find(step_ref);
      if (ss_pos != db_search_steps.end())
      {
        protein.setSearchParameters(exportDBSearchParameters_(ss_pos->second));
      }

      proteins.push_back(protein);
    }
  }


  MzTab IdentificationData::exportMzTab() const
  {
    MzTabMetaData meta;
    Size counter = 1;
    for (const auto& software : processing_software)
    {
      MzTabSoftwareMetaData sw_meta;
      sw_meta.software.setName(software.getName());
      sw_meta.software.setValue(software.getVersion());
      meta.software[counter] = sw_meta;
      ++counter;
    }
    counter = 1;
    map<InputFileRef, Size> file_map;
    for (const String& input_file : input_files)
    {
      MzTabMSRunMetaData run_meta;
      run_meta.location.set(input_file);
      meta.ms_run[counter] = run_meta;
      file_map[&input_file] = counter;
      ++counter;
    }
    set<String> fixed_mods, variable_mods;
    for (const auto& search_param : db_search_params)
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

    map<ScoreTypeRef, Size> protein_scores, peptide_scores, psm_scores,
      nucleic_acid_scores, oligonucleotide_scores, osm_scores;
    // compound_scores;

    MzTabProteinSectionRows proteins;
    MzTabNucleicAcidSectionRows nucleic_acids;
    for (const auto& parent : parent_molecules)
    {
      if (parent.molecule_type == MT_PROTEIN)
      {
        exportParentMoleculeToMzTab_(parent, proteins, protein_scores);
      }
      else if (parent.molecule_type == MT_RNA)
      {
        exportParentMoleculeToMzTab_(parent, nucleic_acids,
                                     nucleic_acid_scores);
      }
    }

    MzTabPeptideSectionRows peptides;
    for (const auto& peptide : identified_peptides)
    {
      exportPeptideOrOligoToMzTab_(peptide, peptides, peptide_scores);
    }

    MzTabOligonucleotideSectionRows oligos;
    for (const auto& oligo : identified_oligos)
    {
      exportPeptideOrOligoToMzTab_(oligo, oligos, oligonucleotide_scores);
    }

    MzTabPSMSectionRows psms;
    MzTabOSMSectionRows osms;
    counter = 1;
    for (const auto& query_match : query_matches)
    {
      IdentifiedMoleculeRef molecule_ref = query_match.identified_molecule_ref;
      // @TODO: what about small molecules?
      enum MoleculeType molecule_type = query_match.getMoleculeType();
      if (molecule_type == MT_PROTEIN)
      {
        const AASequence& seq = query_match.getIdentifiedPeptideRef()->sequence;
        double calc_mass = seq.getMonoWeight(Residue::Full, query_match.charge);
        exportQueryMatchToMzTab_(seq.toString(), query_match, calc_mass, psms,
                                 psm_scores, file_map);
        psms.back().PSM_ID.set(counter);
        ++counter;
      }
      else if (molecule_type == MT_RNA)
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


  void IdentificationData::exportScoresToMzTab_(
    const ScoreList& scores, map<Size, MzTabDouble>& output,
    map<ScoreTypeRef, Size>& score_map) const
  {
    for (const pair<ScoreTypeRef, double>& score_pair : scores)
    {
      if (!score_map.count(score_pair.first)) // new score type
      {
        score_map.insert(make_pair(score_pair.first, score_map.size() + 1));
      }
      Size index = score_map[score_pair.first];
      output[index].set(score_pair.second);
    }
  }


  void IdentificationData::exportProcessingStepsToMzTab_(
    const vector<ProcessingStepRef>& steps, MzTabParameterList& output) const
  {
    vector<MzTabParameter> search_engines;
    for (const ProcessingStepRef& step_ref : steps)
    {
      const Software& sw = *step_ref->software_ref;
      MzTabParameter param;
      param.setName(sw.getName());
      param.setValue(sw.getVersion());
      search_engines.push_back(param);
    }
    output.set(search_engines);
  }


  void IdentificationData::addMzTabSEScores_(
    const map<ScoreTypeRef, Size>& scores, map<Size, MzTabParameter>& output)
    const
  {
    for (const pair<ScoreTypeRef, Size>& score_pair : scores)
    {
      const ScoreType& score_type = *score_pair.first;
      MzTabParameter param;
      param.setName(score_type.name); // or "score_type.cv_term.getName()"?
      param.setAccession(score_type.cv_term.getAccession());
      param.setCVLabel(score_type.cv_term.getCVIdentifierRef());
      output[score_pair.second] = param;
    }
  }


  void IdentificationData::addMzTabMoleculeParentContext_(
    const set<MoleculeParentMatch>& matches,
    const MzTabOligonucleotideSectionRow& row,
    vector<MzTabOligonucleotideSectionRow>& output) const
  {
    for (const MoleculeParentMatch& match : matches)
    {
      MzTabOligonucleotideSectionRow copy = row;
      if (match.left_neighbor == MoleculeParentMatch::LEFT_TERMINUS)
      {
        copy.pre.set("-");
      }
      else if (match.left_neighbor != MoleculeParentMatch::UNKNOWN_NEIGHBOR)
      {
        copy.pre.set(String(match.left_neighbor));
      }
      if (match.right_neighbor == MoleculeParentMatch::RIGHT_TERMINUS)
      {
        copy.post.set("-");
      }
      else if (match.right_neighbor != MoleculeParentMatch::UNKNOWN_NEIGHBOR)
      {
        copy.post.set(String(match.right_neighbor));
      }
      if (match.start_pos != MoleculeParentMatch::UNKNOWN_POSITION)
      {
        copy.start.set(String(match.start_pos + 1));
      }
      if (match.end_pos != MoleculeParentMatch::UNKNOWN_POSITION)
      {
        copy.end.set(String(match.end_pos + 1));
      }
      output.push_back(copy);
    }
  }


  void IdentificationData::addMzTabMoleculeParentContext_(
    const set<MoleculeParentMatch>& /* matches */,
    const MzTabPeptideSectionRow& /* row */,
    vector<MzTabPeptideSectionRow>& /* output */) const
  {
    // do nothing here
  }


  IdentificationData::SearchParamRef
  IdentificationData::importDBSearchParameters_(
    const ProteinIdentification::SearchParameters& pisp)
  {
    DBSearchParam dbsp;
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

    return registerDBSearchParam(dbsp);
  }


  ProteinIdentification::SearchParameters
  IdentificationData::exportDBSearchParameters_(SearchParamRef ref) const
  {
    const DBSearchParam& dbsp = *ref;
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
    for (const pair<ScoreTypeRef, double>& score_pair : scores)
    {
      if (!isValidReference_(score_pair.first, score_types))
      {
        String msg = "invalid reference to a score type - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
  }


  void IdentificationData::checkProcessingSteps_(
    const std::vector<ProcessingStepRef>& step_refs)
  {
    for (ProcessingStepRef step_ref : step_refs)
    {
      if (!isValidReference_(step_ref, processing_steps))
      {
        String msg = "invalid reference to a data processing step - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
  }


  void IdentificationData::checkParentMatches_(const ParentMatches& matches,
                                               enum MoleculeType expected_type)
  {
    for (const auto& pair : matches)
    {
      if (!isValidReference_(pair.first, parent_molecules))
      {
        String msg = "invalid reference to a parent molecule - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      if (pair.first->molecule_type != expected_type)
      {
        String msg = "unexpected molecule type for parent molecule";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
  }


  IdentificationData::InputFileRef
  IdentificationData::registerInputFile(const String& file)
  {
    return &(*input_files.insert(file).first);
  }


  IdentificationData::ProcessingSoftwareRef
  IdentificationData::registerDataProcessingSoftware(const Software& software)
  {
    return &(*processing_software.insert(software).first);
  }


  IdentificationData::SearchParamRef
  IdentificationData::registerDBSearchParam(const DBSearchParam& param)
  {
    // @TODO: any required information that should be checked?
    return &(*db_search_params.insert(param).first);
  }


  IdentificationData::ProcessingStepRef
  IdentificationData::registerDataProcessingStep(
    const DataProcessingStep& step, SearchParamRef search_ref)
  {
    // valid reference to software is required:
    if (!isValidReference_(step.software_ref, processing_software))
    {
      String msg = "invalid reference to data processing software - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    // if given, references to input files must be valid:
    for (InputFileRef ref : step.input_file_refs)
    {
      if (!isValidReference_(ref, input_files))
      {
        String msg = "invalid reference to input file - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }

    ProcessingStepRef step_ref = &(*processing_steps.insert(step).first);
    // if given, reference to DB search param. must be valid:
    if (search_ref != nullptr)
    {
      if (!isValidReference_(search_ref, db_search_params))
      {
        String msg = "invalid reference to database search parameters - register those first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
      }
      db_search_steps.insert(make_pair(step_ref, search_ref));
    }
    return step_ref;
  }


  IdentificationData::ScoreTypeRef
  IdentificationData::registerScoreType(const ScoreType& score)
  {
    pair<ScoreTypes::iterator, bool> result;
    if ((score.software_ref == nullptr) && (current_step_ref_ != nullptr))
    {
      // transfer the software ref. from the current data processing step:
      const DataProcessingStep& step = *current_step_ref_;
      ScoreType copy(score); // need a copy so we can modify it
      copy.software_ref = step.software_ref;
      result = score_types.insert(copy);
    }
    else
    {
      // ref. to software may be missing, but must otherwise be valid:
      if ((score.software_ref != nullptr) &&
          !isValidReference_(score.software_ref, processing_software))
      {
        String msg = "invalid reference to data processing software - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      result = score_types.insert(score);
    }
    if (!result.second && (score.higher_better != result.first->higher_better))
    {
      String msg = "score type already exists with opposite orientation";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    return &(*result.first);
  }


  IdentificationData::DataQueryRef
  IdentificationData::registerDataQuery(const DataQuery& query)
  {
    // reference to spectrum or feature is required:
    if (query.data_id.empty())
    {
      String msg = "missing identifier in data query";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    // ref. to input file may be missing, but must otherwise be valid:
    if ((query.input_file_ref != nullptr) &&
        !isValidReference_(query.input_file_ref, input_files))
    {
      String msg = "invalid reference to an input file - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    return &(*data_queries.insert(query).first);
  }


  IdentificationData::IdentifiedPeptideRef
  IdentificationData::registerPeptide(const IdentifiedPeptide& peptide)
  {
    if (peptide.sequence.empty())
    {
      String msg = "missing sequence for peptide";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    checkParentMatches_(peptide.parent_matches, MT_PROTEIN);

    return insertIntoMultiIndex_(identified_peptides, peptide);
  }


  IdentificationData::IdentifiedCompoundRef
  IdentificationData::registerCompound(const IdentifiedCompound& compound)
  {
    if (compound.identifier.empty())
    {
      String msg = "missing identifier for compound";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    return insertIntoMultiIndex_(identified_compounds, compound);
  }


  IdentificationData::IdentifiedOligoRef
  IdentificationData::registerOligo(const IdentifiedOligo& oligo)
  {
    if (oligo.sequence.empty())
    {
      String msg = "missing sequence for oligonucleotide";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    checkParentMatches_(oligo.parent_matches, MT_RNA);

    return insertIntoMultiIndex_(identified_oligos, oligo);
  }


  IdentificationData::ParentMoleculeRef
  IdentificationData::registerParentMolecule(const ParentMolecule& parent)
  {
    if (parent.accession.empty())
    {
      String msg = "missing accession for parent molecule";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    return insertIntoMultiIndex_(parent_molecules, parent);
  }


  IdentificationData::QueryMatchRef
  IdentificationData::registerMoleculeQueryMatch(const MoleculeQueryMatch&
                                                 match)
  {
    if (const IdentifiedPeptideRef* ref_ptr =
        boost::get<IdentifiedPeptideRef>(&match.identified_molecule_ref))
    {
      if (!isValidReference_(*ref_ptr, identified_peptides))
      {
        String msg = "invalid reference to an identified peptide - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
    else if (const IdentifiedCompoundRef* ref_ptr =
             boost::get<IdentifiedCompoundRef>(&match.identified_molecule_ref))
    {
      if (!isValidReference_(*ref_ptr, identified_compounds))
      {
        String msg = "invalid reference to an identified compound - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
    else if (const IdentifiedOligoRef* ref_ptr =
             boost::get<IdentifiedOligoRef>(&match.identified_molecule_ref))
    {
      if (!isValidReference_(*ref_ptr, identified_oligos))
      {
        String msg = "invalid reference to an identified oligonucleotide - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
    if (!isValidReference_(match.data_query_ref, data_queries))
    {
      String msg = "invalid reference to a data query - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    return insertIntoMultiIndex_(query_matches, match);
  }


  void IdentificationData::setCurrentProcessingStep(ProcessingStepRef step_ref)
  {
    if (!isValidReference_(step_ref, processing_steps))
    {
      String msg = "invalid reference to a processing step - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    current_step_ref_ = step_ref;
  }


  IdentificationData::ProcessingStepRef
  IdentificationData::getCurrentProcessingStep()
  {
    return current_step_ref_;
  }


  void IdentificationData::clearCurrentProcessingStep()
  {
    current_step_ref_ = nullptr;
  }


  IdentificationData::ScoreTypeRef IdentificationData::findScoreType(
    const String& score_name, ProcessingSoftwareRef software_ref) const
  {
    for (ScoreTypes::iterator it = score_types.begin(); it != score_types.end();
         ++it)
    {
      if ((it->name == score_name) &&
          ((software_ref == nullptr) || (it->software_ref == software_ref)))
      {
        return &(*it);
      }
    }
    return nullptr;
  }


  vector<IdentificationData::QueryMatchRef>
  IdentificationData::getBestMatchPerQuery(ScoreTypeRef score_ref) const
  {
    vector<QueryMatchRef> results;
    bool higher_better = score_ref->higher_better;
    pair<double, bool> best_score = make_pair(0.0, false);
    QueryMatchRef best_ref = nullptr;
    for (const MoleculeQueryMatch& match : query_matches)
    {
      pair<double, bool> current_score = match.getScore(score_ref);
      if (match.data_query_ref != best_ref->data_query_ref)
      {
        // finalize previous query:
        if (best_score.second) results.push_back(best_ref);
        best_score = current_score;
        best_ref = &match;
      }
      else if (current_score.second &&
               (!best_score.second ||
                isBetterScore(current_score.first, best_score.first,
                              higher_better)))
      {
        // new best score for the current query:
        best_score = current_score;
        best_ref = &match;
      }
    }
    // finalize last query:
    if (best_score.second) results.push_back(best_ref);

    return results;
  }

} // end namespace OpenMS
