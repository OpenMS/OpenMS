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
      if (!isValidHashedReference_(pair.first, parent_molecule_lookup_))
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
    if (file.name.empty())
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
  IdentificationData::registerDataProcessingSoftware(
    const DataProcessingSoftware& software)
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
    return processing_softwares_.insert(software).first;
  }


  IdentificationData::SearchParamRef
  IdentificationData::registerDBSearchParam(const DBSearchParam& param)
  {
    // @TODO: any required information that should be checked?
    return db_search_params_.insert(param).first;
  }


  IdentificationData::ProcessingStepRef
  IdentificationData::registerDataProcessingStep(
    const DataProcessingStep& step)
  {
    return registerDataProcessingStep(step, db_search_params_.end());
  }


  IdentificationData::ProcessingStepRef
  IdentificationData::registerDataProcessingStep(
    const DataProcessingStep& step, SearchParamRef search_ref)
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

    ProcessingStepRef step_ref = processing_steps_.insert(step).first;
    // if given, reference to DB search param. must be valid:
    if (search_ref != db_search_params_.end())
    {
      if (!isValidReference_(search_ref, db_search_params_))
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
    if (score.cv_term.getName().empty())
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
    if (query.input_file_opt && !isValidReference_(*query.input_file_opt,
                                                   input_files_))
    {
      String msg = "invalid reference to an input file - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    // can't use "insertIntoMultiIndex_" because DataQuery doesn't have the
    // "steps_and_scores" member (from ScoredProcessingResult)
    auto result = data_queries_.insert(query);
    if (!result.second) // existing element - merge in new information
    {
      data_queries_.modify(result.first, [&query](DataQuery& existing)
                           {
                             existing += query;
                           });
    }

    data_query_lookup_.insert(uintptr_t(&(*result.first)));

    // @TODO: add processing step (currently not supported by DataQuery)
    return result.first;
  }


  IdentificationData::IdentifiedPeptideRef
  IdentificationData::registerIdentifiedPeptide(const IdentifiedPeptide&
                                                peptide)
  {
    if (peptide.sequence.empty())
    {
      String msg = "missing sequence for peptide";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    checkParentMatches_(peptide.parent_matches, MoleculeType::PROTEIN);

    return insertIntoMultiIndex_(identified_peptides_, peptide,
                                 identified_peptide_lookup_);
  }


  IdentificationData::IdentifiedCompoundRef
  IdentificationData::registerIdentifiedCompound(const IdentifiedCompound&
                                                 compound)
  {
    if (compound.identifier.empty())
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
    if (oligo.sequence.empty())
    {
      String msg = "missing sequence for oligonucleotide";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }
    checkParentMatches_(oligo.parent_matches, MoleculeType::RNA);

    return insertIntoMultiIndex_(identified_oligos_, oligo,
                                 identified_oligo_lookup_);
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
    if ((parent.coverage < 0.0) || (parent.coverage > 1.0))
    {
      String msg = "parent molecule coverage must be between 0 and 1";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    return insertIntoMultiIndex_(parent_molecules_, parent,
                                 parent_molecule_lookup_);
  }


  void IdentificationData::registerParentMoleculeGrouping(
    const ParentMoleculeGrouping& grouping)
  {
    checkAppliedProcessingSteps_(grouping.steps_and_scores);

    for (const auto& group : grouping.groups)
    {
      checkScoreTypes_(group.scores);

      for (const auto& ref : group.parent_molecule_refs)
      {
        if (!isValidHashedReference_(ref, parent_molecule_lookup_))
        {
          String msg = "invalid reference to a parent molecule - register that first";
          throw Exception::IllegalArgument(__FILE__, __LINE__,
                                           OPENMS_PRETTY_FUNCTION, msg);
        }
      }
    }

    parent_molecule_groupings_.push_back(grouping);

    // add the current processing step?
    if ((current_step_ref_ != processing_steps_.end()) &&
        (grouping.steps_and_scores.get<1>().find(current_step_ref_) ==
         grouping.steps_and_scores.get<1>().end()))
    {
      parent_molecule_groupings_.back().steps_and_scores.push_back(
        IdentificationDataInternal::AppliedProcessingStep(current_step_ref_));
    }
  }


  IdentificationData::QueryMatchRef
  IdentificationData::registerMoleculeQueryMatch(const MoleculeQueryMatch&
                                                 match)
  {
    if (const IdentifiedPeptideRef* ref_ptr =
        boost::get<IdentifiedPeptideRef>(&match.identified_molecule_ref))
    {
      if (!isValidHashedReference_(*ref_ptr, identified_peptide_lookup_))
      {
        String msg = "invalid reference to an identified peptide - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
    else if (const IdentifiedCompoundRef* ref_ptr =
             boost::get<IdentifiedCompoundRef>(&match.identified_molecule_ref))
    {
      if (!isValidHashedReference_(*ref_ptr, identified_compound_lookup_))
      {
        String msg = "invalid reference to an identified compound - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }
    else if (const IdentifiedOligoRef* ref_ptr =
             boost::get<IdentifiedOligoRef>(&match.identified_molecule_ref))
    {
      if (!isValidHashedReference_(*ref_ptr, identified_oligo_lookup_))
      {
        String msg = "invalid reference to an identified oligonucleotide - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }

    if (!isValidHashedReference_(match.data_query_ref, data_query_lookup_))
    {
      String msg = "invalid reference to a data query - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    return insertIntoMultiIndex_(query_matches_, match, query_match_lookup_);
  }


  IdentificationData::MatchGroupRef
  IdentificationData::registerQueryMatchGroup(const QueryMatchGroup& group)
  {
    for (const auto& ref : group.query_match_refs)
    {
      if (!isValidHashedReference_(ref, query_match_lookup_))
      {
        String msg = "invalid reference to a molecule-query match - register that first";
        throw Exception::IllegalArgument(__FILE__, __LINE__,
                                         OPENMS_PRETTY_FUNCTION, msg);
      }
    }

    return insertIntoMultiIndex_(query_match_groups_, group);
  }


  void IdentificationData::addScore(QueryMatchRef match_ref,
                                    ScoreTypeRef score_ref, double value)
  {
    if (!isValidReference_(score_ref, score_types_))
    {
      String msg = "invalid reference to a score type - register that first";
      throw Exception::IllegalArgument(__FILE__, __LINE__,
                                       OPENMS_PRETTY_FUNCTION, msg);
    }

    ModifyMultiIndexAddScore<MoleculeQueryMatch> modifier(score_ref, value);
    query_matches_.modify(match_ref, modifier);
  }


  void IdentificationData::setCurrentProcessingStep(ProcessingStepRef step_ref)
  {
    if (!isValidReference_(step_ref, processing_steps_))
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


  pair<IdentificationData::ScoreTypeRef, bool>
  IdentificationData::findScoreType(const String& score_name) const
  {
    for (ScoreTypeRef it = score_types_.begin(); it != score_types_.end(); ++it)
    {
      if (it->cv_term.getName() == score_name)
      {
        return make_pair(it, true);
      }
    }
    return make_pair(score_types_.end(), false);
  }


  vector<IdentificationData::QueryMatchRef>
  IdentificationData::getBestMatchPerQuery(ScoreTypeRef score_ref) const
  {
    vector<QueryMatchRef> results;
    pair<double, bool> best_score = make_pair(0.0, false);
    QueryMatchRef best_ref = query_matches_.end();
    for (QueryMatchRef ref = query_matches_.begin();
         ref != query_matches_.end(); ++ref)
    {
      pair<double, bool> current_score = ref->getScore(score_ref);
      if ((best_ref != query_matches_.end()) &&
          (ref->data_query_ref != best_ref->data_query_ref))
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
    map<ParentMoleculeRef, ParentData> parent_info;

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
    for (ParentMoleculeRef ref = parent_molecules_.begin();
         ref != parent_molecules_.end(); ++ref)
    {
      auto pos = parent_info.find(ref);
      double coverage = (pos == parent_info.end()) ? 0.0 : pos->second.coverage;
      parent_molecules_.modify(ref, [coverage](ParentMolecule& parent)
                               {
                                 parent.coverage = coverage;
                               });
    }
  }


  void IdentificationData::cleanup(bool require_query_match,
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
      parent_molecule_lookup_.clear(); // will become invalid anyway
      for (const auto& grouping: parent_molecule_groupings_)
      {
        for (const auto& group : grouping.groups)
        {
          for (const auto& ref : group.parent_molecule_refs)
          {
            parent_molecule_lookup_.insert(ref);
          }
        }
      }
      removeFromSetIfNotHashed_(parent_molecules_, parent_molecule_lookup_);
    }
    // update look-up table of parent molecule addresses (in case parent
    // molecules were removed):
    updateAddressLookup_(parent_molecules_, parent_molecule_lookup_);

    // remove parent matches based on parent molecules:
    ModifyMultiIndexRemoveParentMatches<IdentifiedPeptide>
      pep_modifier(parent_molecule_lookup_);
    for (auto it = identified_peptides_.begin();
         it != identified_peptides_.end(); ++it)
    {
      identified_peptides_.modify(it, pep_modifier);
    }
    ModifyMultiIndexRemoveParentMatches<IdentifiedOligo>
      oli_modifier(parent_molecule_lookup_);
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
    set<IdentifiedMoleculeRef> id_refs;
    for (auto it = identified_peptides_.begin();
         it != identified_peptides_.end(); ++it)
    {
      id_refs.insert(it);
    }
    for (auto it = identified_compounds_.begin();
         it != identified_compounds_.end(); ++it)
    {
      id_refs.insert(it);
    }
    for (auto it = identified_oligos_.begin();
         it != identified_oligos_.end(); ++it)
    {
      id_refs.insert(it);
    }
    removeFromSetIf_(query_matches_, [&](MoleculeQueryMatches::iterator it)
                     {
                       return !id_refs.count(it->identified_molecule_ref);
                     });

    // remove molecule-query matches based on query match groups:
    if (require_match_group)
    {
      query_match_lookup_.clear(); // will become invalid anyway
      for (const auto& group : query_match_groups_)
      {
        for (const auto& ref : group.query_match_refs)
        {
          query_match_lookup_.insert(ref);
        }
      }
      removeFromSetIfNotHashed_(query_matches_, query_match_lookup_);
    }
    // update look-up table of query match addresses:
    updateAddressLookup_(query_matches_, query_match_lookup_);

    // remove id'd molecules and data queries based on molecule-query matches:
    if (require_query_match)
    {
      data_query_lookup_.clear();
      identified_peptide_lookup_.clear();
      identified_compound_lookup_.clear();
      identified_oligo_lookup_.clear();
      for (const auto& match : query_matches_)
      {
        data_query_lookup_.insert(match.data_query_ref);
        IdentificationData::MoleculeType molecule_type =
          match.getMoleculeType();
        if (molecule_type == IdentificationData::MoleculeType::PROTEIN)
        {
          identified_peptide_lookup_.insert(match.getIdentifiedPeptideRef());
        }
        else if (molecule_type == IdentificationData::MoleculeType::COMPOUND)
        {
          identified_compound_lookup_.insert(match.getIdentifiedCompoundRef());
        }
        else if (molecule_type == IdentificationData::MoleculeType::RNA)
        {
          identified_oligo_lookup_.insert(match.getIdentifiedOligoRef());
        }
      }
      removeFromSetIfNotHashed_(data_queries_, data_query_lookup_);
      removeFromSetIfNotHashed_(identified_peptides_,
                                identified_peptide_lookup_);
      removeFromSetIfNotHashed_(identified_compounds_,
                                identified_compound_lookup_);
      removeFromSetIfNotHashed_(identified_oligos_, identified_oligo_lookup_);
    }
    // update look-up tables of addresses:
    updateAddressLookup_(data_queries_, data_query_lookup_);
    updateAddressLookup_(identified_peptides_, identified_peptide_lookup_);
    updateAddressLookup_(identified_compounds_, identified_compound_lookup_);
    updateAddressLookup_(identified_oligos_, identified_oligo_lookup_);

    // remove parent molecules based on identified molecules:
    if (require_identified_sequence)
    {
      parent_molecule_lookup_.clear(); // will become invalid anyway
      for (const auto& peptide : identified_peptides_)
      {
        for (const auto& parent_pair : peptide.parent_matches)
        {
          parent_molecule_lookup_.insert(parent_pair.first);
        }
      }
      for (const auto& oligo : identified_oligos_)
      {
        for (const auto& parent_pair : oligo.parent_matches)
        {
          parent_molecule_lookup_.insert(parent_pair.first);
        }
      }
      removeFromSetIfNotHashed_(parent_molecules_, parent_molecule_lookup_);
      // update look-up table of parent molecule addresses (again):
      updateAddressLookup_(parent_molecules_, parent_molecule_lookup_);
    }

    // remove entries from parent molecule groups based on parent molecules:
    bool warn = false;
    for (auto& grouping : parent_molecule_groupings_)
    {
      for (auto group_it = grouping.groups.begin();
           group_it != grouping.groups.end(); )
      {
        Size old_size = group_it->parent_molecule_refs.size();
        grouping.groups.modify(
          group_it, [&](ParentMoleculeGroup& group)
          {
            removeFromSetIfNotHashed_(group.parent_molecule_refs,
                                      parent_molecule_lookup_);
          });
        if (group_it->parent_molecule_refs.empty())
        {
          group_it = grouping.groups.erase(group_it);
        }
        else
        {
          if (group_it->parent_molecule_refs.size() != old_size)
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
    for (auto group_it = query_match_groups_.begin();
         group_it != query_match_groups_.end(); )
    {
      Size old_size = group_it->query_match_refs.size();
      query_match_groups_.modify(
        group_it, [&](QueryMatchGroup& group)
        {
          removeFromSetIfNotHashed_(group.query_match_refs,
                                    query_match_lookup_);
        });
      if (group_it->query_match_refs.empty())
      {
        group_it = query_match_groups_.erase(group_it);
      }
      else
      {
        if (group_it->query_match_refs.size() != old_size)
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
            data_queries_.empty() && parent_molecules_.empty() &&
            parent_molecule_groupings_.empty() &&
            identified_peptides_.empty() && identified_compounds_.empty() &&
            identified_oligos_.empty() && query_matches_.empty() &&
            query_match_groups_.empty());
  }

} // end namespace OpenMS
