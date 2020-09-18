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
      if ((step.processing_step_opt != boost::none) &&
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
  IdentificationData::registerInputFile(const InputFile& file)
  {
    if (!no_checks_ && file.name.empty())
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
                            existing += file;
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


  IdentificationData::InputItemRef
  IdentificationData::registerInputItem(const InputItem& query)
  {
    // reference to spectrum or feature is required:
    if (!no_checks_ && query.data_id.empty())
    {
      String msg = "missing identifier in data query";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    // ref. to input file may be missing, but must otherwise be valid:
    if (!no_checks_ && query.input_file_opt &&
        !isValidReference_(*query.input_file_opt, input_files_))
    {
      String msg = "invalid reference to an input file - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    // can't use "insertIntoMultiIndex_" because InputItem doesn't have the
    // "steps_and_scores" member (from ScoredProcessingResult)
    auto result = input_items_.insert(query);
    if (!result.second) // existing element - merge in new information
    {
      input_items_.modify(result.first, [&query](InputItem& existing)
                           {
                             existing += query;
                           });
    }

    input_item_lookup_.insert(uintptr_t(&(*result.first)));

    // @TODO: add processing step (currently not supported by InputItem)
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
        String msg = "missing accession for parent molecule";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
      if ((parent.coverage < 0.0) || (parent.coverage > 1.0))
      {
        String msg = "parent molecule coverage must be between 0 and 1";
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
        checkScoreTypes_(group.scores);

        for (const auto& ref : group.parent_refs)
        {
          if (!isValidHashedReference_(ref, parent_lookup_))
          {
            String msg = "invalid reference to a parent molecule - register that first";
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


  IdentificationData::InputMatchRef
  IdentificationData::registerInputMatch(const InputMatch&
                                                 match)
  {
    if (!no_checks_)
    {
      if (const IdentifiedPeptideRef* ref_ptr =
          boost::get<IdentifiedPeptideRef>(&match.identified_molecule_var))
      {
        if (!isValidHashedReference_(*ref_ptr, identified_peptide_lookup_))
        {
          String msg = "invalid reference to an identified peptide - register that first";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                           OPENMS_PRETTY_FUNCTION, msg);
        }
      }
      else if (const IdentifiedCompoundRef* ref_ptr =
               boost::get<IdentifiedCompoundRef>(&match.identified_molecule_var))
      {
        if (!isValidHashedReference_(*ref_ptr, identified_compound_lookup_))
        {
          String msg = "invalid reference to an identified compound - register that first";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                           OPENMS_PRETTY_FUNCTION, msg);
        }
      }
      else if (const IdentifiedOligoRef* ref_ptr =
               boost::get<IdentifiedOligoRef>(&match.identified_molecule_var))
      {
        if (!isValidHashedReference_(*ref_ptr, identified_oligo_lookup_))
        {
          String msg = "invalid reference to an identified oligonucleotide - register that first";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                           OPENMS_PRETTY_FUNCTION, msg);
        }
      }

      if (!isValidHashedReference_(match.input_item_ref, input_item_lookup_))
      {
        String msg = "invalid reference to a data query - register that first";
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

    return insertIntoMultiIndex_(input_matches_, match, input_match_lookup_);
  }


  IdentificationData::MatchGroupRef
  IdentificationData::registerInputMatchGroup(const InputMatchGroup& group)
  {
    if (!no_checks_)
    {
      for (const auto& ref : group.input_match_refs)
      {
        if (!isValidHashedReference_(ref, input_match_lookup_))
        {
          String msg = "invalid reference to a molecule-query match - register that first";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                           OPENMS_PRETTY_FUNCTION, msg);
        }
      }
    }

    return insertIntoMultiIndex_(input_match_groups_, group);
  }


  void IdentificationData::addScore(InputMatchRef match_ref,
                                    ScoreTypeRef score_ref, double value)
  {
    if (!no_checks_ && !isValidReference_(score_ref, score_types_))
    {
      String msg = "invalid reference to a score type - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    ModifyMultiIndexAddScore<InputMatch> modifier(score_ref, value);
    input_matches_.modify(match_ref, modifier);
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


  vector<IdentificationData::InputMatchRef>
  IdentificationData::getBestMatchPerQuery(ScoreTypeRef score_ref) const
  {
    vector<InputMatchRef> results;
    pair<double, bool> best_score = make_pair(0.0, false);
    InputMatchRef best_ref = input_matches_.end();
    for (InputMatchRef ref = input_matches_.begin();
         ref != input_matches_.end(); ++ref)
    {
      pair<double, bool> current_score = ref->getScore(score_ref);
      if ((best_ref != input_matches_.end()) &&
          (ref->input_item_ref != best_ref->input_item_ref))
      {
        // finalize previous query:
        if (best_score.second) results.push_back(best_ref);
        best_score = current_score;
        best_ref = ref;
      }
      else if (current_score.second &&
               (!best_score.second ||
                score_ref->isBetterScore(current_score.first,
                                         best_score.first)))
      {
        // new best score for the current query:
        best_score = current_score;
        best_ref = ref;
      }
    }
    // finalize last query:
    if (best_score.second) results.push_back(best_ref);

    return results;
  }


  void IdentificationData::calculateCoverages(bool check_molecule_length)
  {
    // aggregate molecule-parent matches by parent:
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
        if (pos == parent_info.end()) // new parent molecule
        {
          ParentData pd;
          pd.length = AASequence::fromString(pair.first->sequence).size();
          if (pd.length == 0) break; // sequence not available
          pos = parent_info.insert(make_pair(pair.first, pd)).first;
        }
        Size parent_length = pos->second.length; // always check this
        for (const auto& match : pair.second)
        {
          if (match.hasValidPositions(molecule_length, parent_length))
          {
            pos->second.fragments.push_back(make_pair(match.start_pos,
                                                      match.end_pos));
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
        if (pos == parent_info.end()) // new parent molecule
        {
          ParentData pd;
          pd.length = NASequence::fromString(pair.first->sequence).size();
          if (pd.length == 0) break; // sequence not available
          pos = parent_info.insert(make_pair(pair.first, pd)).first;
        }
        Size parent_length = pos->second.length; // always check this
        for (const auto& match : pair.second)
        {
          if (match.hasValidPositions(molecule_length, parent_length))
          {
            pos->second.fragments.push_back(make_pair(match.start_pos,
                                                      match.end_pos));
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


  void IdentificationData::cleanup(bool require_input_match,
                                   bool require_identified_sequence,
                                   bool require_parent_match,
                                   bool require_parent_group,
                                   bool require_match_group)
  {
    // we expect that only "primary results" (stored in classes derived from
    // "ScoredProcessingResult") will be directly removed (by filters) - not
    // meta data (incl. data queries, score types, processing steps etc.)

    // remove parent molecules based on parent groups:
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
    // update look-up table of parent molecule addresses (in case parent
    // molecules were removed):
    updateAddressLookup_(parents_, parent_lookup_);

    // remove parent matches based on parent molecules:
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

    // remove molecule-query matches based on identified molecules:
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
    removeFromSetIf_(input_matches_, [&](InputMatches::iterator it)
                     {
                       return !id_vars.count(it->identified_molecule_var);
                     });

    // remove molecule-query matches based on query match groups:
    if (require_match_group)
    {
      input_match_lookup_.clear(); // will become invalid anyway
      for (const auto& group : input_match_groups_)
      {
        for (const auto& ref : group.input_match_refs)
        {
          input_match_lookup_.insert(ref);
        }
      }
      removeFromSetIfNotHashed_(input_matches_, input_match_lookup_);
    }
    // update look-up table of query match addresses:
    updateAddressLookup_(input_matches_, input_match_lookup_);

    // remove id'd molecules and data queries based on molecule-query matches:
    if (require_input_match)
    {
      input_item_lookup_.clear();
      identified_peptide_lookup_.clear();
      identified_compound_lookup_.clear();
      identified_oligo_lookup_.clear();
      for (const auto& match : input_matches_)
      {
        input_item_lookup_.insert(match.input_item_ref);
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
            break;
          default: // shouldn't happen
            break;
        }
      }
      removeFromSetIfNotHashed_(input_items_, input_item_lookup_);
      removeFromSetIfNotHashed_(identified_peptides_,
                                identified_peptide_lookup_);
      removeFromSetIfNotHashed_(identified_compounds_,
                                identified_compound_lookup_);
      removeFromSetIfNotHashed_(identified_oligos_, identified_oligo_lookup_);
    }
    // update look-up tables of addresses:
    updateAddressLookup_(input_items_, input_item_lookup_);
    updateAddressLookup_(identified_peptides_, identified_peptide_lookup_);
    updateAddressLookup_(identified_compounds_, identified_compound_lookup_);
    updateAddressLookup_(identified_oligos_, identified_oligo_lookup_);

    // remove parent molecules based on identified molecules:
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
      // update look-up table of parent molecule addresses (again):
      updateAddressLookup_(parents_, parent_lookup_);
    }

    // remove entries from parent molecule groups based on parent molecules:
    bool warn = false;
    for (auto& groups : parent_groups_)
    {
      for (auto group_it = groups.groups.begin();
           group_it != groups.groups.end(); )
      {
        Size old_size = group_it->parent_refs.size();
        groups.groups.modify(group_it, [&](ParentGroup& group)
        {
          removeFromSetIfNotHashed_(group.parent_refs, parent_lookup_);
        });
        if (group_it->parent_refs.empty())
        {
          group_it = groups.groups.erase(group_it);
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
      OPENMS_LOG_WARN << "Warning: filtering removed elements from parent molecule groups - associated scores may not be valid any more" << endl;
    }

    // remove entries from query match groups based on molecule-query matches:
    warn = false;
    for (auto group_it = input_match_groups_.begin();
         group_it != input_match_groups_.end(); )
    {
      Size old_size = group_it->input_match_refs.size();
      input_match_groups_.modify(group_it, [&](InputMatchGroup& group)
      {
        removeFromSetIfNotHashed_(group.input_match_refs, input_match_lookup_);
      });
      if (group_it->input_match_refs.empty())
      {
        group_it = input_match_groups_.erase(group_it);
      }
      else
      {
        if (group_it->input_match_refs.size() != old_size)
        {
          warn = true;
        }
        ++group_it;
      }
    }
    if (warn)
    {
      OPENMS_LOG_WARN << "Warning: filtering removed elements from query match groups - associated scores may not be valid any more" << endl;
    }
  }


  bool IdentificationData::empty() const
  {
    return (input_files_.empty() && processing_softwares_.empty() &&
            processing_steps_.empty() && db_search_params_.empty() &&
            db_search_steps_.empty() && score_types_.empty() &&
            input_items_.empty() && parents_.empty() &&
            parent_groups_.empty() &&
            identified_peptides_.empty() && identified_compounds_.empty() &&
            identified_oligos_.empty() && adducts_.empty() &&
            input_matches_.empty() && input_match_groups_.empty());
  }


  void IdentificationData::mergeScoredProcessingResults_(
    IdentificationData::ScoredProcessingResult& result,
    const IdentificationData::ScoredProcessingResult& other,
    const map<ProcessingStepRef, ProcessingStepRef>& step_refs,
    const map<ScoreTypeRef, ScoreTypeRef>& score_refs)
  {
    result.MetaInfoInterface::operator=(other);
    for (const AppliedProcessingStep& applied : other.steps_and_scores)
    {
      AppliedProcessingStep copy;
      if (applied.processing_step_opt)
      {
        copy.processing_step_opt = step_refs.at(*applied.processing_step_opt);
      }
      for (const auto& pair : applied.scores)
      {
        ScoreTypeRef score_ref = score_refs.at(pair.first);
        copy.scores[score_ref] = pair.second;
      }
      result.addProcessingStep(copy);
    }
  }


  IdentificationData::ProcessingStepRef
  IdentificationData::merge(const IdentificationData& other)
  {
    no_checks_ = true;
    // input files:
    map<InputFileRef, InputFileRef> file_refs;
    for (InputFileRef other_ref = other.getInputFiles().begin();
         other_ref != other.getInputFiles().end(); ++other_ref)
    {
      file_refs[other_ref] = registerInputFile(*other_ref);
    }
    // score types:
    map<ScoreTypeRef, ScoreTypeRef> score_refs;
    for (ScoreTypeRef other_ref = other.getScoreTypes().begin();
         other_ref != other.getScoreTypes().end(); ++other_ref)
    {
      score_refs[other_ref] = registerScoreType(*other_ref);
    }
    // processing software:
    map<ProcessingSoftwareRef, ProcessingSoftwareRef> sw_refs;
    for (ProcessingSoftwareRef other_ref = other.getProcessingSoftwares().begin();
         other_ref != other.getProcessingSoftwares().end(); ++other_ref)
    {
      // update internal references:
      ProcessingSoftware copy = *other_ref;
      for (ScoreTypeRef& score_ref : copy.assigned_scores)
      {
        score_ref = score_refs[score_ref];
      }
      sw_refs[other_ref] = registerProcessingSoftware(copy);
    }
    // search params:
    map<SearchParamRef, SearchParamRef> param_refs;
    for (SearchParamRef other_ref = other.getDBSearchParams().begin();
         other_ref != other.getDBSearchParams().end(); ++other_ref)
    {
      param_refs[other_ref] = registerDBSearchParam(*other_ref);
    }
    // processing steps:
    map<ProcessingStepRef, ProcessingStepRef> step_refs;
    for (ProcessingStepRef other_ref = other.getProcessingSteps().begin();
         other_ref != other.getProcessingSteps().end(); ++other_ref)
    {
      // update internal references:
      ProcessingStep copy = *other_ref;
      copy.software_ref = sw_refs[copy.software_ref];
      for (InputFileRef& file_ref : copy.input_file_refs)
      {
        file_ref = file_refs[file_ref];
      }
      step_refs[other_ref] = registerProcessingStep(copy);
    }
    // search steps:
    for (const auto& pair : other.getDBSearchSteps())
    {
      ProcessingStepRef step_ref = step_refs[pair.first];
      SearchParamRef param_ref = param_refs[pair.second];
      db_search_steps_[step_ref] = param_ref;
    }
    // data queries:
    map<InputItemRef, InputItemRef> query_refs;
    for (InputItemRef other_ref = other.getInputItems().begin();
         other_ref != other.getInputItems().end(); ++other_ref)
    {
      // update internal references:
      InputItem copy = *other_ref;
      if (copy.input_file_opt)
      {
        copy.input_file_opt = file_refs[*copy.input_file_opt];
      }
      query_refs[other_ref] = registerInputItem(copy);
    }
    // parent molecules:
    map<ParentSequenceRef, ParentSequenceRef> parent_refs;
    for (ParentSequenceRef other_ref = other.getParentSequences().begin();
         other_ref != other.getParentSequences().end(); ++other_ref)
    {
      // don't copy processing steps and scores yet:
      ParentSequence copy(other_ref->accession, other_ref->molecule_type,
                          other_ref->sequence, other_ref->description,
                          other_ref->coverage, other_ref->is_decoy);
      // now copy precessing steps and scores while updating references:
      mergeScoredProcessingResults_(copy, *other_ref, step_refs, score_refs);
      parent_refs[other_ref] = registerParentSequence(copy);
    }
    // identified peptides:
    map<IdentifiedPeptideRef, IdentifiedPeptideRef> peptide_refs;
    for (IdentifiedPeptideRef other_ref = other.getIdentifiedPeptides().begin();
         other_ref != other.getIdentifiedPeptides().end(); ++other_ref)
    {
      // don't copy parent matches, steps/scores yet:
      IdentifiedPeptide copy(other_ref->sequence, ParentMatches());
      // now copy steps/scores and parent matches while updating references:
      mergeScoredProcessingResults_(copy, *other_ref, step_refs, score_refs);
      for (const auto& pair : other_ref->parent_matches)
      {
        ParentSequenceRef parent_ref = parent_refs[pair.first];
        copy.parent_matches[parent_ref] = pair.second;
      }
      // @TODO: with C++17, 'map::extract' offers a better way to update keys
      peptide_refs[other_ref] = registerIdentifiedPeptide(copy);
    }
    // identified oligonucleotides:
    map<IdentifiedOligoRef, IdentifiedOligoRef> oligo_refs;
    for (IdentifiedOligoRef other_ref = other.getIdentifiedOligos().begin();
         other_ref != other.getIdentifiedOligos().end(); ++other_ref)
    {
      // don't copy parent matches, steps/scores yet:
      IdentifiedOligo copy(other_ref->sequence, ParentMatches());
      // now copy steps/scores and parent matches while updating references:
      mergeScoredProcessingResults_(copy, *other_ref, step_refs, score_refs);
      for (const auto& pair : other_ref->parent_matches)
      {
        ParentSequenceRef parent_ref = parent_refs[pair.first];
        copy.parent_matches[parent_ref] = pair.second;
      }
     // @TODO: with C++17, 'map::extract' offers a better way to update keys
      oligo_refs[other_ref] = registerIdentifiedOligo(copy);
    }
    // identified compounds:
    map<IdentifiedCompoundRef, IdentifiedCompoundRef> compound_refs;
    for (IdentifiedCompoundRef other_ref = other.getIdentifiedCompounds().begin();
         other_ref != other.getIdentifiedCompounds().end(); ++other_ref)
    {
      IdentifiedCompound copy(other_ref->identifier, other_ref->formula,
                              other_ref->name, other_ref->smile, other_ref->inchi);
      mergeScoredProcessingResults_(copy, *other_ref, step_refs, score_refs);
      compound_refs[other_ref] = registerIdentifiedCompound(copy);
    }
    // adducts:
    map<AdductRef, AdductRef> adduct_refs;
    for (AdductRef other_ref = other.getAdducts().begin();
         other_ref != other.getAdducts().end(); ++other_ref)
    {
      adduct_refs[other_ref] = registerAdduct(*other_ref);
    }
    // molecule-query matches:
    map<InputMatchRef, InputMatchRef> match_refs;
    for (InputMatchRef other_ref = other.getInputMatches().begin();
         other_ref != other.getInputMatches().end(); ++other_ref)
    {
      IdentifiedMolecule molecule_var;
      const IdentifiedMolecule& other_var = other_ref->identified_molecule_var;
      switch (other_var.getMoleculeType())
      {
        case MoleculeType::PROTEIN:
          molecule_var = peptide_refs[other_var.getIdentifiedPeptideRef()];
          break;
        case MoleculeType::COMPOUND:
          molecule_var = compound_refs[other_var.getIdentifiedCompoundRef()];
          break;
        case MoleculeType::RNA:
          molecule_var = oligo_refs[other_var.getIdentifiedOligoRef()];
          break;
        default: // avoid compiler warning
          break;
      }
      InputItemRef query_ref = query_refs[other_ref->input_item_ref];
      InputMatch copy(molecule_var, query_ref, other_ref->charge);
      if (other_ref->adduct_opt)
      {
        copy.adduct_opt = adduct_refs[*other_ref->adduct_opt];
      }
      for (const auto& pair : other_ref->peak_annotations)
      {
        boost::optional<ProcessingStepRef> opt_ref;
        if (pair.first)
        {
          opt_ref = step_refs[*pair.first];
        }
        copy.peak_annotations[opt_ref] = pair.second;
      }
      mergeScoredProcessingResults_(copy, *other_ref, step_refs, score_refs);
      match_refs[other_ref] = registerInputMatch(copy);
    }
    // parent molecule groups:
    // @TODO: does this need to be more sophisticated?
    for (const ParentGroupSet& groups : other.parent_groups_)
    {
      ParentGroupSet copy(groups.label);
      mergeScoredProcessingResults_(copy, groups, step_refs, score_refs);
      for (const ParentGroup& group : groups.groups)
      {
        ParentGroup group_copy;
        for (const auto& pair : group.scores)
        {
          ScoreTypeRef score_ref = score_refs[pair.first];
          group_copy.scores[score_ref] = pair.second;
        }
        for (ParentSequenceRef parent_ref : group.parent_refs)
        {
          group_copy.parent_refs.insert(parent_refs[parent_ref]);
        }
        copy.groups.insert(group_copy);
      }
      registerParentGroupSet(copy);
    }
    no_checks_ = false;

    // this is only needed for reusing "merge" in the copy c'tor:
    if (other.current_step_ref_ == other.processing_steps_.end())
    {
      return processing_steps_.end();
    }
    return step_refs[other.current_step_ref_];
  }


  IdentificationData::IdentificationData(const IdentificationData& other):
    MetaInfoInterface(other)
  {
    // don't add a processing step during merging:
    current_step_ref_ = processing_steps_.end();
    current_step_ref_ = merge(other);
    no_checks_ = other.no_checks_;
  }
} // end namespace OpenMS
