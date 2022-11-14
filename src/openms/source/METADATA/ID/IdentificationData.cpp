// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2022.
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

#include <OpenMS/METADATA/ID/IdentificationData.h>
#include <OpenMS/CHEMISTRY/ProteaseDB.h>
#include <numeric>

using namespace std;

namespace OpenMS
{
  void IdentificationData::checkScoreTypes_(const map<ScoreTypeRef, double>&
                                            scores) const
  {
    for (const auto& pair : scores)
    {
      if (!isValidReference_(pair.first, score_types_))
      {
        String msg = "invalid reference to a score type - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
  }

  void IdentificationData::checkAppliedProcessingSteps_(
    const AppliedProcessingSteps& steps_and_scores) const
  {
    for (const auto& step : steps_and_scores)
    {
      if ((step.processing_step_opt != std::nullopt) &&
          (!isValidReference_(*step.processing_step_opt, processing_steps_)))
      {
        String msg = "invalid reference to a data processing step - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      checkScoreTypes_(step.scores);
    }
  }


  void IdentificationData::checkParentMatches_(const ParentMatches& matches,
                                               MoleculeType expected_type) const
  {
    for (const auto& pair : matches)
    {
      if (!isValidHashedReference_(pair.first, parent_lookup_))
      {
        String msg = "invalid reference to a parent sequence - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      if (pair.first->molecule_type != expected_type)
      {
        String msg = "unexpected molecule type for parent sequence";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
  }


  IdentificationData::InputFileRef
  IdentificationData::registerInputFile(const InputFile& file)
  {
    if (!no_checks_ && file.name.empty()) // key may not be empty
    {
      String msg = "input file must have a name";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    auto result = input_files_.insert(file);
    if (!result.second) // existing element - merge in new information
    {
      input_files_.modify(result.first, [&file](InputFile& existing)
                          {
                            existing.merge(file);
                          });
    }

    return result.first;
  }


  IdentificationData::ProcessingSoftwareRef
  IdentificationData::registerProcessingSoftware(
    const ProcessingSoftware& software)
  {
    if (!no_checks_)
    {
      for (ScoreTypeRef score_ref : software.assigned_scores)
      {
        if (!isValidReference_(score_ref, score_types_))
        {
          String msg = "invalid reference to a score type - register that first";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                           OPENMS_PRETTY_FUNCTION, msg);
        }
      }
    }
    return processing_softwares_.insert(software).first;
  }


  IdentificationData::SearchParamRef
  IdentificationData::registerDBSearchParam(const DBSearchParam& param)
  {
    // @TODO: any required information that should be checked?
    return db_search_params_.insert(param).first;
  }


  IdentificationData::ProcessingStepRef
  IdentificationData::registerProcessingStep(
    const ProcessingStep& step)
  {
    return registerProcessingStep(step, db_search_params_.end());
  }


  IdentificationData::ProcessingStepRef
  IdentificationData::registerProcessingStep(
    const ProcessingStep& step, SearchParamRef search_ref)
  {
    if (!no_checks_)
    {
      // valid reference to software is required:
      if (!isValidReference_(step.software_ref, processing_softwares_))
      {
        String msg = "invalid reference to data processing software - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      // if given, references to input files must be valid:
      for (InputFileRef ref : step.input_file_refs)
      {
        if (!isValidReference_(ref, input_files_))
        {
          String msg = "invalid reference to input file - register that first";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                           OPENMS_PRETTY_FUNCTION, msg);
        }
      }
    }

    ProcessingStepRef step_ref = processing_steps_.insert(step).first;
    // if given, reference to DB search param. must be valid:
    if (search_ref != db_search_params_.end())
    {
      if (!no_checks_ && !isValidReference_(search_ref, db_search_params_))
      {
        String msg = "invalid reference to database search parameters - register those first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      db_search_steps_.insert(make_pair(step_ref, search_ref));
    }
    return step_ref;
  }


  IdentificationData::ScoreTypeRef
  IdentificationData::registerScoreType(const ScoreType& score)
  {
    // @TODO: allow just an accession? (all look-ups are currently by name)
    if (!no_checks_ && score.cv_term.getName().empty())
    {
      String msg = "score type must have a name (as part of its CV term)";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    pair<ScoreTypes::iterator, bool> result;
    result = score_types_.insert(score);
    if (!result.second && (score.higher_better != result.first->higher_better))
    {
      String msg = "score type already exists with opposite orientation";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    return result.first;
  }

  IdentificationData::ObservationRef
  IdentificationData::registerObservation(const Observation& obs)
  {
    if (!no_checks_)
    {
      // reference to spectrum or feature is required:
      if (obs.data_id.empty())
      {
        String msg = "missing identifier in observation";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      // ref. to input file must be valid:
      if (!isValidReference_(obs.input_file, input_files_))
      {
        String msg = "invalid reference to an input file - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }

    // can't use "insertIntoMultiIndex_" because Observation doesn't have the
    // "steps_and_scores" member (from ScoredProcessingResult)
    auto result = observations_.insert(obs);
    if (!result.second) // existing element - merge in new information
    {
      observations_.modify(result.first, [&obs](Observation& existing)
                           {
                             existing.merge(obs);
                           });
    }
    // add address of new element to look-up table (for existence checks):
    observation_lookup_.insert(uintptr_t(&(*result.first)));

    // @TODO: add processing step? (currently not supported by Observation)
    return result.first;
  }


  IdentificationData::IdentifiedPeptideRef
  IdentificationData::registerIdentifiedPeptide(const IdentifiedPeptide&
                                                peptide)
  {
    if (!no_checks_)
    {
      if (peptide.sequence.empty())
      {
        String msg = "missing sequence for peptide";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      checkParentMatches_(peptide.parent_matches, MoleculeType::PROTEIN);
    }

    return insertIntoMultiIndex_(identified_peptides_, peptide,
                                 identified_peptide_lookup_);
  }


  IdentificationData::IdentifiedCompoundRef
  IdentificationData::registerIdentifiedCompound(const IdentifiedCompound&
                                                 compound)
  {
    if (!no_checks_ && compound.identifier.empty())
    {
      String msg = "missing identifier for compound";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    return insertIntoMultiIndex_(identified_compounds_, compound,
                                 identified_compound_lookup_);
  }


  IdentificationData::IdentifiedOligoRef
  IdentificationData::registerIdentifiedOligo(const IdentifiedOligo& oligo)
  {
    if (!no_checks_)
    {
      if (oligo.sequence.empty())
      {
        String msg = "missing sequence for oligonucleotide";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      checkParentMatches_(oligo.parent_matches, MoleculeType::RNA);
    }

    return insertIntoMultiIndex_(identified_oligos_, oligo,
                                 identified_oligo_lookup_);
  }


  IdentificationData::ParentSequenceRef
  IdentificationData::registerParentSequence(const ParentSequence& parent)
  {
    if (!no_checks_)
    {
      if (parent.accession.empty())
      {
        String msg = "missing accession for parent sequence";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      if ((parent.coverage < 0.0) || (parent.coverage > 1.0))
      {
        String msg = "parent sequence coverage must be between 0 and 1";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }

    return insertIntoMultiIndex_(parents_, parent,
                                 parent_lookup_);
  }


  void IdentificationData::registerParentGroupSet(const ParentGroupSet& groups)
  {
    if (!no_checks_)
    {
      checkAppliedProcessingSteps_(groups.steps_and_scores);

      for (const auto& group : groups.groups)
      {
        checkScoreTypes_(group.scores); // are the score types registered?

        for (const auto& ref : group.parent_refs)
        {
          if (!isValidHashedReference_(ref, parent_lookup_))
          {
            String msg = "invalid reference to a parent sequence - register that first";
            throw Exception::IllegalArgument(__FILE__, __LINE__,
                                             OPENMS_PRETTY_FUNCTION, msg);
          }
        }
      }
    }

    parent_groups_.push_back(groups);

    // add the current processing step?
    if ((current_step_ref_ != processing_steps_.end()) &&
        (groups.steps_and_scores.get<1>().find(current_step_ref_) ==
         groups.steps_and_scores.get<1>().end()))
    {
      parent_groups_.back().steps_and_scores.push_back(
        IdentificationDataInternal::AppliedProcessingStep(current_step_ref_));
    }
  }


  IdentificationData::AdductRef
  IdentificationData::registerAdduct(const AdductInfo& adduct)
  {
    // @TODO: require non-empty name? (auto-generate from formula?)
    auto result = adducts_.insert(adduct);
    if (!result.second && (result.first->getName() != adduct.getName()))
    {
      OPENMS_LOG_WARN << "Warning: adduct '" << adduct.getName()
                      << "' is already known under the name '"
                      << result.first->getName() << "'";
    }
    return result.first;
  }


  IdentificationData::ObservationMatchRef
  IdentificationData::registerObservationMatch(const ObservationMatch& match)
  {
    if (!no_checks_)
    {
      if (const IdentifiedPeptideRef* ref_ptr =
          std::get_if<IdentifiedPeptideRef>(&match.identified_molecule_var))
      {
        if (!isValidHashedReference_(*ref_ptr, identified_peptide_lookup_))
        {
          String msg = "invalid reference to an identified peptide - register that first";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                           OPENMS_PRETTY_FUNCTION, msg);
        }
      }
      else if (const IdentifiedCompoundRef* ref_ptr =
               std::get_if<IdentifiedCompoundRef>(&match.identified_molecule_var))
      {
        if (!isValidHashedReference_(*ref_ptr, identified_compound_lookup_))
        {
          String msg = "invalid reference to an identified compound - register that first";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                           OPENMS_PRETTY_FUNCTION, msg);
        }
      }
      else if (const IdentifiedOligoRef* ref_ptr =
               std::get_if<IdentifiedOligoRef>(&match.identified_molecule_var))
      {
        if (!isValidHashedReference_(*ref_ptr, identified_oligo_lookup_))
        {
          String msg = "invalid reference to an identified oligonucleotide - register that first";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                           OPENMS_PRETTY_FUNCTION, msg);
        }
      }

      if (!isValidHashedReference_(match.observation_ref, observation_lookup_))
      {
        String msg = "invalid reference to an observation - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }

      if (match.adduct_opt && !isValidReference_(*match.adduct_opt, adducts_))
      {
        String msg = "invalid reference to an adduct - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }

    return insertIntoMultiIndex_(observation_matches_, match,
                                 observation_match_lookup_);
  }


  IdentificationData::MatchGroupRef
  IdentificationData::registerObservationMatchGroup(const ObservationMatchGroup& group)
  {
    if (!no_checks_)
    {
      for (const auto& ref : group.observation_match_refs)
      {
        if (!isValidHashedReference_(ref, observation_match_lookup_))
        {
          String msg = "invalid reference to an input match - register that first";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                           OPENMS_PRETTY_FUNCTION, msg);
        }
      }
    }

    return insertIntoMultiIndex_(observation_match_groups_, group);
  }


  void IdentificationData::addScore(ObservationMatchRef match_ref,
                                    ScoreTypeRef score_ref, double value)
  {
    if (!no_checks_ && !isValidReference_(score_ref, score_types_))
    {
      String msg = "invalid reference to a score type - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    ModifyMultiIndexAddScore<ObservationMatch> modifier(score_ref, value);
    observation_matches_.modify(match_ref, modifier);
  }


  void IdentificationData::setCurrentProcessingStep(ProcessingStepRef step_ref)
  {
    if (!no_checks_ && !isValidReference_(step_ref, processing_steps_))
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
    current_step_ref_ = processing_steps_.end();
  }


  IdentificationData::ScoreTypeRef
  IdentificationData::findScoreType(const String& score_name) const
  {
    for (ScoreTypeRef it = score_types_.begin(); it != score_types_.end(); ++it)
    {
      if (it->cv_term.getName() == score_name)
      {
        return it;
      }
    }
    return score_types_.end();
  }


  vector<IdentificationData::ObservationMatchRef>
  IdentificationData::getBestMatchPerObservation(ScoreTypeRef score_ref,
                                                 bool require_score) const
  {
    vector<ObservationMatchRef> results;
    pair<double, bool> best_score = make_pair(0.0, false);
    ObservationMatchRef best_ref = observation_matches_.end();
    Size n_matches = 1; // number of matches for current observation
    // matches for same observation appear consecutively, so just iterate:
    for (ObservationMatchRef ref = observation_matches_.begin();
         ref != observation_matches_.end(); ++ref, ++n_matches)
    {
      pair<double, bool> current_score = ref->getScore(score_ref);
      if (current_score.second && (!best_score.second ||
                                   score_ref->isBetterScore(current_score.first,
                                                           best_score.first)))
      {
        // new best score for the current observation:
        best_score = current_score;
        best_ref = ref;
      }
      // peek ahead:
      ObservationMatchRef next = ref;
      ++next;
      if ((next == observation_matches_.end()) ||
          (next->observation_ref != ref->observation_ref))
      {
        // last match for this observation - finalize:
        if (best_score.second)
        {
          results.push_back(best_ref);
        }
        else if (!require_score && (n_matches == 1))
        {
          results.push_back(ref); // only match for this observation
        }
        best_score.second = false;
        n_matches = 0; // will be incremented by for-loop
      }
    }

    return results;
  }


  pair<IdentificationData::ObservationMatchRef, IdentificationData::ObservationMatchRef>
  IdentificationData::getMatchesForObservation(ObservationRef obs_ref) const
  {
    return observation_matches_.equal_range(obs_ref);
  }


  void IdentificationData::calculateCoverages(bool check_molecule_length)
  {
    // aggregate parent matches by parent:
    struct ParentData
    {
      Size length = 0;
      double coverage = 0.0;
      vector<pair<Size, Size>> fragments;
    };
    map<ParentSequenceRef, ParentData> parent_info;

    // go through all peptides:
    for (const auto& molecule : identified_peptides_)
    {
      Size molecule_length = check_molecule_length ?
        molecule.sequence.size() : 0;
      for (const auto& pair : molecule.parent_matches)
      {
        auto pos = parent_info.find(pair.first);
        if (pos == parent_info.end()) // new parent sequence
        {
          ParentData pd;
          pd.length = AASequence::fromString(pair.first->sequence).size();
          if (pd.length == 0)
          {
            break; // sequence not available
          }
          pos = parent_info.insert(make_pair(pair.first, pd)).first;
        }
        Size parent_length = pos->second.length; // always check this
        for (const auto& match : pair.second)
        {
          if (match.hasValidPositions(molecule_length, parent_length))
          {
            pos->second.fragments.emplace_back(match.start_pos,
                                                      match.end_pos);
          }
        }
      }
    }
    // go through all oligonucleotides:
    for (const auto& molecule : identified_oligos_)
    {
      Size molecule_length = check_molecule_length ?
        molecule.sequence.size() : 0;
      for (const auto& pair : molecule.parent_matches)
      {
        auto pos = parent_info.find(pair.first);
        if (pos == parent_info.end()) // new parent sequence
        {
          ParentData pd;
          pd.length = NASequence::fromString(pair.first->sequence).size();
          if (pd.length == 0)
          {
            break; // sequence not available
          }
          pos = parent_info.insert(make_pair(pair.first, pd)).first;
        }
        Size parent_length = pos->second.length; // always check this
        for (const auto& match : pair.second)
        {
          if (match.hasValidPositions(molecule_length, parent_length))
          {
            pos->second.fragments.emplace_back(match.start_pos,
                                                      match.end_pos);
          }
        }
      }
    }

    // calculate coverage for each parent:
    for (auto& pair : parent_info)
    {
      vector<bool> covered(pair.second.length, false);
      for (const auto& fragment : pair.second.fragments)
      {
        fill(covered.begin() + fragment.first,
             covered.begin() + fragment.second + 1, true);
      }
      pair.second.coverage = (accumulate(covered.begin(), covered.end(), 0) /
                              double(pair.second.length));
    }
    // set coverage:
    for (ParentSequenceRef ref = parents_.begin();
         ref != parents_.end(); ++ref)
    {
      auto pos = parent_info.find(ref);
      double coverage = (pos == parent_info.end()) ? 0.0 : pos->second.coverage;
      parents_.modify(ref, [coverage](ParentSequence& parent)
                               {
                                 parent.coverage = coverage;
                               });
    }
  }


  void IdentificationData::cleanup(bool require_observation_match,
                                   bool require_identified_sequence,
                                   bool require_parent_match,
                                   bool require_parent_group,
                                   bool require_match_group)
  {
    // we expect that only "primary results" (stored in classes derived from
    // "ScoredProcessingResult") will be directly removed (by filters) - not
    // meta data (incl. score types, processing steps etc.)

    // remove parent sequences based on parent groups:
    if (require_parent_group)
    {
      parent_lookup_.clear(); // will become invalid anyway
      for (const auto& groups: parent_groups_)
      {
        for (const auto& group : groups.groups)
        {
          for (const auto& ref : group.parent_refs)
          {
            parent_lookup_.insert(ref);
          }
        }
      }
      removeFromSetIfNotHashed_(parents_, parent_lookup_);
    }
    // update look-up table of parent sequence addresses (in case parent
    // molecules were removed):
    updateAddressLookup_(parents_, parent_lookup_);

    // remove parent matches based on parent sequences:
    ModifyMultiIndexRemoveParentMatches<IdentifiedPeptide>
      pep_modifier(parent_lookup_);
    for (auto it = identified_peptides_.begin();
         it != identified_peptides_.end(); ++it)
    {
      identified_peptides_.modify(it, pep_modifier);
    }
    ModifyMultiIndexRemoveParentMatches<IdentifiedOligo>
      oli_modifier(parent_lookup_);
    for (auto it = identified_oligos_.begin();
         it != identified_oligos_.end(); ++it)
    {
      identified_oligos_.modify(it, oli_modifier);
    }

    // remove identified molecules based on parent matches:
    if (require_parent_match)
    {
      removeFromSetIf_(identified_peptides_, [](IdentifiedPeptides::iterator it)
                       {
                         return it->parent_matches.empty();
                       });
      removeFromSetIf_(identified_oligos_, [](IdentifiedOligos::iterator it)
                       {
                         return it->parent_matches.empty();
                       });
    }

    // remove observation matches based on identified molecules:
    set<IdentifiedMolecule> id_vars;
    for (IdentifiedPeptideRef it = identified_peptides_.begin();
         it != identified_peptides_.end(); ++it)
    {
      id_vars.insert(it);
    }
    for (IdentifiedCompoundRef it = identified_compounds_.begin();
         it != identified_compounds_.end(); ++it)
    {
      id_vars.insert(it);
    }
    for (IdentifiedOligoRef it = identified_oligos_.begin();
         it != identified_oligos_.end(); ++it)
    {
      id_vars.insert(it);
    }
    removeFromSetIf_(observation_matches_, [&](ObservationMatches::iterator it)
                     {
                       return !id_vars.count(it->identified_molecule_var);
                     });

    // remove observation matches based on observation match groups:
    if (require_match_group)
    {
      observation_match_lookup_.clear(); // will become invalid anyway
      for (const auto& group : observation_match_groups_)
      {
        for (const auto& ref : group.observation_match_refs)
        {
          observation_match_lookup_.insert(ref);
        }
      }
      removeFromSetIfNotHashed_(observation_matches_, observation_match_lookup_);
    }
    // update look-up table of input match addresses:
    updateAddressLookup_(observation_matches_, observation_match_lookup_);

    // remove id'd molecules, observations and adducts based on observation matches:
    if (require_observation_match)
    {
      observation_lookup_.clear();
      identified_peptide_lookup_.clear();
      identified_compound_lookup_.clear();
      identified_oligo_lookup_.clear();
      set<AdductRef> adduct_refs;
      for (const auto& match : observation_matches_)
      {
        observation_lookup_.insert(match.observation_ref);
        const IdentifiedMolecule& molecule_var = match.identified_molecule_var;
        switch (molecule_var.getMoleculeType())
        {
          case IdentificationData::MoleculeType::PROTEIN:
            identified_peptide_lookup_.insert(molecule_var.getIdentifiedPeptideRef());
            break;
          case IdentificationData::MoleculeType::COMPOUND:
            identified_compound_lookup_.insert(molecule_var.getIdentifiedCompoundRef());
            break;
          case IdentificationData::MoleculeType::RNA:
            identified_oligo_lookup_.insert(molecule_var.getIdentifiedOligoRef());
        }
        if (match.adduct_opt) adduct_refs.insert(*match.adduct_opt);
      }
      removeFromSetIfNotHashed_(observations_, observation_lookup_);
      removeFromSetIfNotHashed_(identified_peptides_,
                                identified_peptide_lookup_);
      removeFromSetIfNotHashed_(identified_compounds_,
                                identified_compound_lookup_);
      removeFromSetIfNotHashed_(identified_oligos_, identified_oligo_lookup_);
      removeFromSetIf_(adducts_, [&](Adducts::iterator it)
      {
        return !adduct_refs.count(it);
      });
    }
    // update look-up tables of addresses:
    updateAddressLookup_(observations_, observation_lookup_);
    updateAddressLookup_(identified_peptides_, identified_peptide_lookup_);
    updateAddressLookup_(identified_compounds_, identified_compound_lookup_);
    updateAddressLookup_(identified_oligos_, identified_oligo_lookup_);

    // remove parent sequences based on identified molecules:
    if (require_identified_sequence)
    {
      parent_lookup_.clear(); // will become invalid anyway
      for (const auto& peptide : identified_peptides_)
      {
        for (const auto& parent_pair : peptide.parent_matches)
        {
          parent_lookup_.insert(parent_pair.first);
        }
      }
      for (const auto& oligo : identified_oligos_)
      {
        for (const auto& parent_pair : oligo.parent_matches)
        {
          parent_lookup_.insert(parent_pair.first);
        }
      }
      removeFromSetIfNotHashed_(parents_, parent_lookup_);
      // update look-up table of parent sequence addresses (again):
      updateAddressLookup_(parents_, parent_lookup_);
    }

    // remove entries from parent sequence groups based on parent sequences
    // (if a parent sequence doesn't exist anymore, remove it from any groups):
    bool warn = false;
    for (auto& group_set : parent_groups_)
    {
      for (auto group_it = group_set.groups.begin();
           group_it != group_set.groups.end(); )
      {
        Size old_size = group_it->parent_refs.size();
        group_set.groups.modify(group_it, [&](ParentGroup& group)
        {
          removeFromSetIfNotHashed_(group.parent_refs, parent_lookup_);
        });
        if (group_it->parent_refs.empty())
        {
          group_it = group_set.groups.erase(group_it);
        }
        else
        {
          if (group_it->parent_refs.size() != old_size)
          {
            warn = true;
          }
          ++group_it;
        }
      }
      // @TODO: if no group is left, remove the whole grouping?
    }
    if (warn)
    {
      OPENMS_LOG_WARN << "Warning: filtering removed elements from parent sequence groups - associated scores may not be valid any more" << endl;
    }

    // remove entries from input match groups based on input matches:
    warn = false;
    for (auto group_it = observation_match_groups_.begin();
         group_it != observation_match_groups_.end(); )
    {
      Size old_size = group_it->observation_match_refs.size();
      observation_match_groups_.modify(group_it, [&](ObservationMatchGroup& group)
      {
        removeFromSetIfNotHashed_(group.observation_match_refs, observation_match_lookup_);
      });
      if (group_it->observation_match_refs.empty())
      {
        group_it = observation_match_groups_.erase(group_it);
      }
      else
      {
        if (group_it->observation_match_refs.size() != old_size)
        {
          warn = true;
        }
        ++group_it;
      }
    }
    if (warn)
    {
      OPENMS_LOG_WARN << "Warning: filtering removed elements from observation match groups - associated scores may not be valid any more" << endl;
    }
  }


  bool IdentificationData::empty() const
  {
    return (input_files_.empty() && processing_softwares_.empty() &&
            processing_steps_.empty() && db_search_params_.empty() &&
            db_search_steps_.empty() && score_types_.empty() &&
            observations_.empty() && parents_.empty() &&
            parent_groups_.empty() &&
            identified_peptides_.empty() && identified_compounds_.empty() &&
            identified_oligos_.empty() && adducts_.empty() &&
            observation_matches_.empty() && observation_match_groups_.empty());
  }


  void IdentificationData::mergeScoredProcessingResults_(
    IdentificationData::ScoredProcessingResult& result,
    const IdentificationData::ScoredProcessingResult& other,
    const RefTranslator& trans)
  {
    result.MetaInfoInterface::operator=(other);
    for (const AppliedProcessingStep& applied : other.steps_and_scores)
    {
      AppliedProcessingStep copy;
      if (applied.processing_step_opt)
      {
        // need to reference a processing step in 'result', not the original one
        // from 'other', so find the corresponding one:
        copy.processing_step_opt = trans.processing_step_refs.at(*applied.processing_step_opt);
      }
      for (const auto& pair : applied.scores)
      {
        // need to reference a score type in 'result', not the original one from
        // 'other', so find the corresponding one:
        ScoreTypeRef score_ref = trans.score_type_refs.at(pair.first);
        copy.scores[score_ref] = pair.second;
      }
      result.addProcessingStep(copy);
    }
  }


  IdentificationData::RefTranslator
  IdentificationData::merge(const IdentificationData& other)
  {
    RefTranslator trans;
    // incoming data (stored in IdentificationData) is guaranteed to be consistent,
    // so no need to check for consistency again:
    no_checks_ = true;
    // input files:
    for (InputFileRef other_ref = other.getInputFiles().begin();
         other_ref != other.getInputFiles().end(); ++other_ref)
    {
      trans.input_file_refs[other_ref] = registerInputFile(*other_ref);
    }
    // score types:
    for (ScoreTypeRef other_ref = other.getScoreTypes().begin();
         other_ref != other.getScoreTypes().end(); ++other_ref)
    {
      trans.score_type_refs[other_ref] = registerScoreType(*other_ref);
    }
    // processing software:
    for (ProcessingSoftwareRef other_ref = other.getProcessingSoftwares().begin();
         other_ref != other.getProcessingSoftwares().end(); ++other_ref)
    {
      // update internal references:
      ProcessingSoftware copy = *other_ref;
      for (ScoreTypeRef& score_ref : copy.assigned_scores)
      {
        score_ref = trans.score_type_refs[score_ref];
      }
      trans.processing_software_refs[other_ref] = registerProcessingSoftware(copy);
    }
    // search params:
    for (SearchParamRef other_ref = other.getDBSearchParams().begin();
         other_ref != other.getDBSearchParams().end(); ++other_ref)
    {
      trans.search_param_refs[other_ref] = registerDBSearchParam(*other_ref);
    }
    // processing steps:
    for (ProcessingStepRef other_ref = other.getProcessingSteps().begin();
         other_ref != other.getProcessingSteps().end(); ++other_ref)
    {
      // update internal references:
      ProcessingStep copy = *other_ref;
      copy.software_ref = trans.processing_software_refs[copy.software_ref];
      for (InputFileRef& file_ref : copy.input_file_refs)
      {
        file_ref = trans.input_file_refs[file_ref];
      }
      trans.processing_step_refs[other_ref] = registerProcessingStep(copy);
    }
    // search steps:
    for (const auto& pair : other.getDBSearchSteps())
    {
      ProcessingStepRef step_ref = trans.processing_step_refs[pair.first];
      SearchParamRef param_ref = trans.search_param_refs[pair.second];
      db_search_steps_[step_ref] = param_ref;
    }
    // observations:
    for (ObservationRef other_ref = other.getObservations().begin();
         other_ref != other.getObservations().end(); ++other_ref)
    {
      // update internal references:
      Observation copy = *other_ref;
      copy.input_file = trans.input_file_refs[copy.input_file];
      trans.observation_refs[other_ref] = registerObservation(copy);
    }
    // parent sequences:
    for (ParentSequenceRef other_ref = other.getParentSequences().begin();
         other_ref != other.getParentSequences().end(); ++other_ref)
    {
      // don't copy processing steps and scores yet:
      ParentSequence copy(other_ref->accession, other_ref->molecule_type,
                          other_ref->sequence, other_ref->description,
                          other_ref->coverage, other_ref->is_decoy);
      // now copy precessing steps and scores while updating references:
      mergeScoredProcessingResults_(copy, *other_ref, trans);
      trans.parent_sequence_refs[other_ref] = registerParentSequence(copy);
    }
    // identified peptides:
    for (IdentifiedPeptideRef other_ref = other.getIdentifiedPeptides().begin();
         other_ref != other.getIdentifiedPeptides().end(); ++other_ref)
    {
      // don't copy parent matches, steps/scores yet:
      IdentifiedPeptide copy(other_ref->sequence, ParentMatches());
      // now copy steps/scores and parent matches while updating references:
      mergeScoredProcessingResults_(copy, *other_ref, trans);
      for (const auto& pair : other_ref->parent_matches)
      {
        ParentSequenceRef parent_ref = trans.parent_sequence_refs[pair.first];
        copy.parent_matches[parent_ref] = pair.second;
      }
      trans.identified_peptide_refs[other_ref] = registerIdentifiedPeptide(copy);
    }
    // identified oligonucleotides:
    for (IdentifiedOligoRef other_ref = other.getIdentifiedOligos().begin();
         other_ref != other.getIdentifiedOligos().end(); ++other_ref)
    {
      // don't copy parent matches, steps/scores yet:
      IdentifiedOligo copy(other_ref->sequence, ParentMatches());
      // now copy steps/scores and parent matches while updating references:
      mergeScoredProcessingResults_(copy, *other_ref, trans);
      for (const auto& pair : other_ref->parent_matches)
      {
        ParentSequenceRef parent_ref = trans.parent_sequence_refs[pair.first];
        copy.parent_matches[parent_ref] = pair.second;
      }
      trans.identified_oligo_refs[other_ref] = registerIdentifiedOligo(copy);
    }
    // identified compounds:
    for (IdentifiedCompoundRef other_ref = other.getIdentifiedCompounds().begin();
         other_ref != other.getIdentifiedCompounds().end(); ++other_ref)
    {
      IdentifiedCompound copy(other_ref->identifier, other_ref->formula,
                              other_ref->name, other_ref->smile, other_ref->inchi);
      mergeScoredProcessingResults_(copy, *other_ref, trans);
      trans.identified_compound_refs[other_ref] = registerIdentifiedCompound(copy);
    }
    // adducts:
    for (AdductRef other_ref = other.getAdducts().begin();
         other_ref != other.getAdducts().end(); ++other_ref)
    {
      trans.adduct_refs[other_ref] = registerAdduct(*other_ref);
    }
    // observation matches:
    for (ObservationMatchRef other_ref = other.getObservationMatches().begin();
         other_ref != other.getObservationMatches().end(); ++other_ref)
    {
      IdentifiedMolecule molecule_var =
        trans.translate(other_ref->identified_molecule_var);
      ObservationRef obs_ref = trans.observation_refs[other_ref->observation_ref];
      ObservationMatch copy(molecule_var, obs_ref, other_ref->charge);
      if (other_ref->adduct_opt)
      {
        copy.adduct_opt = trans.adduct_refs[*other_ref->adduct_opt];
      }
      for (const auto& pair : other_ref->peak_annotations)
      {
        std::optional<ProcessingStepRef> opt_ref;
        if (pair.first)
        {
          opt_ref = trans.processing_step_refs[*pair.first];
        }
        copy.peak_annotations[opt_ref] = pair.second;
      }
      mergeScoredProcessingResults_(copy, *other_ref, trans);
      trans.observation_match_refs[other_ref] = registerObservationMatch(copy);
    }
    // parent sequence groups:
    // @TODO: does this need to be more sophisticated?
    for (const ParentGroupSet& groups : other.parent_groups_)
    {
      ParentGroupSet copy(groups.label);
      mergeScoredProcessingResults_(copy, groups, trans);
      for (const ParentGroup& group : groups.groups)
      {
        ParentGroup group_copy;
        for (const auto& pair : group.scores)
        {
          ScoreTypeRef score_ref = trans.score_type_refs[pair.first];
          group_copy.scores[score_ref] = pair.second;
        }
        for (ParentSequenceRef parent_ref : group.parent_refs)
        {
          group_copy.parent_refs.insert(trans.parent_sequence_refs[parent_ref]);
        }
        copy.groups.insert(group_copy);
      }
      registerParentGroupSet(copy);
    }
    no_checks_ = false;

    return trans;
  }


  IdentificationData::IdentificationData(const IdentificationData& other):
    MetaInfoInterface(other)
  {
    // don't add a processing step during merging:
    current_step_ref_ = processing_steps_.end();
    RefTranslator trans = merge(other);
    if (other.current_step_ref_ != other.processing_steps_.end())
    {
      current_step_ref_ = trans.processing_step_refs[other.current_step_ref_];
    }
    no_checks_ = other.no_checks_;
  }


  void IdentificationData::swap(IdentificationData& other)
  {
    MetaInfoInterface::swap(other);
    input_files_.swap(other.input_files_);
    processing_softwares_.swap(other.processing_softwares_);
    processing_steps_.swap(other.processing_steps_);
    db_search_params_.swap(other.db_search_params_);
    db_search_steps_.swap(other.db_search_steps_);
    score_types_.swap(other.score_types_);
    observations_.swap(other.observations_);
    parents_.swap(other.parents_);
    parent_groups_.swap(other.parent_groups_);
    identified_peptides_.swap(other.identified_peptides_);
    identified_compounds_.swap(other.identified_compounds_);
    identified_oligos_.swap(other.identified_oligos_);
    adducts_.swap(other.adducts_);
    observation_matches_.swap(other.observation_matches_);
    observation_match_groups_.swap(other.observation_match_groups_);
    std::swap(current_step_ref_, other.current_step_ref_);
    std::swap(no_checks_, other.no_checks_);
    // look-up tables:
    observation_lookup_.swap(other.observation_lookup_);
    parent_lookup_.swap(other.parent_lookup_);
    identified_peptide_lookup_.swap(other.identified_peptide_lookup_);
    identified_compound_lookup_.swap(other.identified_compound_lookup_);
    identified_oligo_lookup_.swap(other.identified_oligo_lookup_);
    observation_match_lookup_.swap(other.observation_match_lookup_);
  }


  void IdentificationData::clear()
  {
    IdentificationData tmp;
    swap(tmp);
  }


  void IdentificationData::setMetaValue(const ObservationMatchRef ref, const String& key,
                                        const DataValue& value)
  {
    setMetaValue_(ref, key, value, observation_matches_, observation_match_lookup_);
  }


  void IdentificationData::setMetaValue(const ObservationRef ref, const String& key,
                                        const DataValue& value)
  {
    setMetaValue_(ref, key, value, observations_, observation_lookup_);
  }


  void IdentificationData::setMetaValue(const IdentifiedMolecule& var, const String& key,
                                        const DataValue& value)
  {
    switch (var.getMoleculeType())
    {
      case MoleculeType::PROTEIN:
        setMetaValue_(var.getIdentifiedPeptideRef(), key, value,
                      identified_peptides_, identified_peptide_lookup_);
        break;
      case MoleculeType::COMPOUND:
        setMetaValue_(var.getIdentifiedCompoundRef(), key, value,
                      identified_compounds_, identified_compound_lookup_);
        break;
      case MoleculeType::RNA:
        setMetaValue_(var.getIdentifiedOligoRef(), key, value,
                      identified_oligos_, identified_oligo_lookup_);
    }
  }

} // end namespace OpenMS
