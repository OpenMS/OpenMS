// Copyright (c) 2002-present, OpenMS Inc. -- EKU Tuebingen, ETH Zurich, and FU Berlin
// SPDX-License-Identifier: BSD-3-Clause
//
// --------------------------------------------------------------------------
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <unordered_map>
#include <array>

using namespace std;
namespace OpenMS
{


  //TODO parameterize so it only adds/keeps best per peptide, peptide charge, modified peptide
  // How? Maybe keep a map here about the best scores and lookup before adding and update and insert only if better
  // proteins of this peptide could be skipped (if we assume same database as we do currently, it has to be there already)
  IDMergerAlgorithm::IDMergerAlgorithm(const String& runIdentifier) :
      IDMergerAlgorithm::DefaultParamHandler("IDMergerAlgorithm"),
      prot_result_(),
      pep_result_(),
      collected_protein_hits_(0, accessionHash_, accessionEqual_),
      id_(runIdentifier)
  {
    defaults_.setValue("annotate_origin",
                       "true",
                       "If true, adds a map_index MetaValue to the PeptideIDs to annotate the IDRun they came from.");
    defaults_.setValidStrings("annotate_origin", {"true","false"});
    defaults_.setValue("allow_disagreeing_settings",
                       "false",
                       "Force merging of disagreeing runs. Use at your own risk.");
    defaults_.setValidStrings("allow_disagreeing_settings", {"true","false"});
    defaultsToParam_();
    prot_result_.setIdentifier(getNewIdentifier_());
  }

  //TODO overload to accept a set of specific runIDs only
  void IDMergerAlgorithm::insertRuns(
      std::vector<ProteinIdentification>&& prots,
      std::vector<PeptideIdentification>&& peps
      )
  {
    if (prots.empty() || peps.empty()) return; //error?

    //TODO instead of only checking consistency, merge if possible (especially for SILAC mods)
    if (!filled_)
    {
      if (prots.size() > 1)
      {
        //Without any exp. design we assume label-free for checking mods
        checkOldRunConsistency_(prots, "label-free");
      }
      copySearchParams_(prots[0], prot_result_);
      filled_ = true;
    }
    else
    {
      //Without any exp. design we assume label-free for checking mods
      checkOldRunConsistency_(prots, this->prot_result_, "label-free");
    }
    // move proteins and move peps
    movePepIDsAndRefProteinsToResultFaster_(std::move(peps), std::move(prots));
  }

  void IDMergerAlgorithm::insertRuns(
      const std::vector<ProteinIdentification>& prots,
      const std::vector<PeptideIdentification>& peps
  )
  {
    //copy
    std::vector<ProteinIdentification> pr = prots;
    std::vector<PeptideIdentification> pep = peps;
    if (prots.empty() || peps.empty()) return; //error?

    //TODO instead of only checking consistency, merge if possible (especially for SILAC mods)
    if (!filled_)
    {
      if (prots.size() > 1)
      {
        //Without any exp. design we assume label-free for checking mods
        checkOldRunConsistency_(prots, "label-free");
      }
      copySearchParams_(prots[0], prot_result_);
      filled_ = true;
    }
    else
    {
      //Without any exp. design we assume label-free for checking mods
      checkOldRunConsistency_(prots, this->prot_result_, "label-free");
    }
    movePepIDsAndRefProteinsToResultFaster_(std::move(pep), std::move(pr));
  }

  void IDMergerAlgorithm::returnResultsAndClear(
      ProteinIdentification& prots,
      vector<PeptideIdentification>& peps)
  {
    // convert the map from file origin to idx into
    // a vector
    StringList newOrigins(file_origin_to_idx_.size());
    for (auto& entry : file_origin_to_idx_)
    {
      newOrigins[entry.second] = entry.first;
    }
    // currently setPrimaryMSRunPath does not support move (const ref)
    prot_result_.setPrimaryMSRunPath(newOrigins);
    std::swap(prots, prot_result_);
    std::swap(peps, pep_result_);
    //reset so the new this class is reuseable
    prot_result_ = ProteinIdentification{};
    prot_result_.setIdentifier(getNewIdentifier_());
    //clear, if user gave non-empty vector
    pep_result_.clear();
    //reset internals
    file_origin_to_idx_.clear();

    for (auto& p : collected_protein_hits_)
      prots.getHits().push_back(std::move(const_cast<ProteinHit&>(p)));
    // above invalidates set but we clear right after

    collected_protein_hits_.clear();
  }

  String IDMergerAlgorithm::getNewIdentifier_() const
  {
    std::array<char, 64> buffer;
    buffer.fill(0);
    time_t rawtime;
    time(&rawtime);
    const auto timeinfo = localtime(&rawtime);
    strftime(buffer.data(), sizeof(buffer), "%d-%m-%Y %H-%M-%S", timeinfo);
    return id_ + String(buffer.data());
  }


  void IDMergerAlgorithm::insertProteinIDs_(
      vector<ProteinIdentification>&& old_protRuns
  )
  {
    typedef std::vector<ProteinHit>::iterator iter_t;
    for (auto& protRun : old_protRuns) //TODO check run ID when option is added
    {
      auto& hits = protRun.getHits();
      collected_protein_hits_.insert(
          std::move_iterator<iter_t>(hits.begin()),
          std::move_iterator<iter_t>(hits.end())
      );
      hits.clear();
    }
  }

  void IDMergerAlgorithm::updateAndMovePepIDs_(
      vector<PeptideIdentification>&& pepIDs,
      const map<String, Size>& runID_to_runIdx,
      const vector<StringList>& originFiles,
      bool annotate_origin)
  {
    //TODO if we allow run IDs, we should do a remove_if,
    // then use the iterator to update and move
    // the IDs, then erase them so we don't encounter them in
    // subsequent calls of this function
    for (auto &pid : pepIDs)
    {
      const String &runID = pid.getIdentifier();

      const auto& runIdxIt = runID_to_runIdx.find(runID);

      if (runIdxIt == runID_to_runIdx.end())
      {
        //This is an easy way to just merge peptides from a certain run
        continue;
        /*
        throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Old IdentificationRun not found for PeptideIdentification "
            "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ").");
        */
      }

      bool annotated = pid.metaValueExists(Constants::UserParam::ID_MERGE_INDEX);
      if (annotate_origin || annotated)
      {
        Size oldFileIdx(0);
        const StringList& origins = originFiles[runIdxIt->second];
        if (annotated)
        {
          oldFileIdx = pid.getMetaValue(Constants::UserParam::ID_MERGE_INDEX);
        }
        else if (origins.size() > 1)
        {
          // If there is more than one possible file it might be from
          // and it is not annotated -> fail
          throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Trying to annotate new id_merge_index for PeptideIdentification "
              "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ") but"
              "no old id_merge_index present");
        }

        if (oldFileIdx >= origins.size())
        {
          throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Trying to annotate new id_merge_index for PeptideIdentification "
              "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ") but"
              " the index exceeds the number of files in the run.");
        }
        pid.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, file_origin_to_idx_[origins[oldFileIdx]]);
      }
      pid.setIdentifier(prot_result_.getIdentifier());
      //move peptides into right vector
      pep_result_.emplace_back(std::move(pid));
    }
  }


  // this merges without checking the existence of a parent protein for the PeptideHits
  // therefore it can merge peptides and proteins separately and a bit faster.
  void IDMergerAlgorithm::movePepIDsAndRefProteinsToResultFaster_(
      vector<PeptideIdentification>&& pepIDs,
      vector<ProteinIdentification>&& old_protRuns
  )
  {
    bool annotate_origin(param_.getValue("annotate_origin").toBool());
    vector<StringList> originFiles{};
    //TODO here check run ID if we allow this option
    for (const auto& protRun : old_protRuns)
    {
      StringList toFill{};
      protRun.getPrimaryMSRunPath(toFill);
      if (toFill.empty() && annotate_origin)
      {
        throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Annotation of origin requested during merge, but no origin present in run "
            + protRun.getIdentifier() + ".");
      }
      //TODO this will make multiple runs from the same file appear multiple times.
      // should be ok but check all possibilities at some point
      originFiles.push_back(toFill);
      for (String& f : toFill)
      {
        file_origin_to_idx_.emplace(std::move(f), file_origin_to_idx_.size());
      }
      toFill.clear();
    }

    std::map<String, Size> runIDToRunIdx;
    for (Size oldProtRunIdx = 0; oldProtRunIdx < old_protRuns.size(); ++oldProtRunIdx)
    {
      ProteinIdentification &protIDRun = old_protRuns[oldProtRunIdx];
      runIDToRunIdx[protIDRun.getIdentifier()] = oldProtRunIdx;
    }

    updateAndMovePepIDs_(std::move(pepIDs), runIDToRunIdx, originFiles, annotate_origin);
    insertProteinIDs_(std::move(old_protRuns));
    pepIDs.clear();
    old_protRuns.clear();
  }

  /* Old version. Quite slower but only copies actually referenced proteins
  void IDMergerAlgorithm::movePepIDsAndRefProteinsToResult_(
      vector<PeptideIdentification>&& pepIDs,
      vector<ProteinIdentification>&& oldProtRuns
  )
  {
    bool annotate_origin(param_.getValue("annotate_origin").toBool());
    vector<StringList> originFiles{};
    //TODO here check run ID if we allow this option
    for (const auto& protRun : oldProtRuns)
    {
      StringList toFill{};
      protRun.getPrimaryMSRunPath(toFill);
      if (toFill.empty() && annotate_origin)
      {
        throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Annotation of origin requested during merge, but no origin present in run "
            + protRun.getIdentifier() + ".");
      }
      //TODO this will make multiple runs from the same file appear multiple times.
      // should be ok but check all possibilities at some point
      originFiles.push_back(toFill);
      for (String& f : toFill)
      {
        file_origin_to_idx_.emplace(std::move(f), file_origin_to_idx_.size());
      }
      toFill.clear();
    }

    for (auto &pid : pepIDs)
    {
      const String &runID = pid.getIdentifier();

      //TODO maybe create lookup table in the beginning runIDToRunRef
      Size oldProtRunIdx = 0;
      for (; oldProtRunIdx < oldProtRuns.size(); ++oldProtRunIdx)
      {
        ProteinIdentification &protIDRun = oldProtRuns[oldProtRunIdx];
        if (protIDRun.getIdentifier() == runID)
        {
          break;
        }
      }
      if (oldProtRunIdx == oldProtRuns.size())
      {
        throw Exception::MissingInformation(
            __FILE__,
            __LINE__,
            OPENMS_PRETTY_FUNCTION,
            "Old IdentificationRun not found for PeptideIdentification "
            "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ").");

      }

      bool annotated = pid.metaValueExists(Constants::UserParam::ID_MERGE_INDEX);
      if (annotate_origin || annotated)
      {
        Size oldFileIdx(0);
        if (annotated)
        {
          oldFileIdx = pid.getMetaValue(Constants::UserParam::ID_MERGE_INDEX);
        }
          // If there is more than one possible file it might be from
          // and it is not annotated -> fail
        else if (originFiles[oldProtRunIdx].size() > 1)
        {
          throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Trying to annotate new id_merge_index for PeptideIdentification "
              "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ") but"
              "no old id_merge_index present");
        }
        pid.setMetaValue(Constants::UserParam::ID_MERGE_INDEX, file_origin_to_idx_[originFiles[oldProtRunIdx].at(oldFileIdx)]);
      }
      pid.setIdentifier(prot_result_.getIdentifier());
      for (auto &phit : pid.getHits())
      {
        for (auto &acc : phit.extractProteinAccessionsSet())
        {
          const auto &it = proteinsCollected.emplace(acc);
          if (it.second) // was newly inserted
          {
            //TODO this linear findHit is not so nice: Maybe we can just insert all proteins into a
            // unordered_set member
            prot_result_.getHits().emplace_back(std::move(*oldProtRuns[oldProtRunIdx].findHit(acc)));
          } // else it was there already
        }
      }
      //move peptides into right vector
      pep_result_.emplace_back(std::move(pid));
    }
    pepIDs.clear();
    oldProtRuns.clear();
  }*/

  void IDMergerAlgorithm::copySearchParams_(const ProteinIdentification& from, ProteinIdentification& to)
  {
      to.setSearchEngine(from.getSearchEngine());
      to.setSearchEngineVersion(from.getSearchEngineVersion());
      to.setSearchParameters(from.getSearchParameters());
  }

  bool IDMergerAlgorithm::checkOldRunConsistency_(const vector<ProteinIdentification>& protRuns, const String& experiment_type) const
  {
    return checkOldRunConsistency_(protRuns, protRuns[0], experiment_type);
  }

  bool IDMergerAlgorithm::checkOldRunConsistency_(const vector<ProteinIdentification>& protRuns, const ProteinIdentification& ref, const String& experiment_type) const
  {
    bool ok = true;
    for (const auto& idRun : protRuns)
    {
      // collect warnings and throw at the end if at least one failed
      ok = ok && ref.peptideIDsMergeable(idRun, experiment_type);
    }
    if (!ok && !param_.getValue("allow_disagreeing_settings").toBool())
    {
      throw Exception::MissingInformation(__FILE__,
                                          __LINE__,
                                          OPENMS_PRETTY_FUNCTION,
                                          "Search settings are not matching across IdentificationRuns. "
                                          "See warnings. Aborting..");
    }
    return ok;
  }
} // namespace OpenMS
