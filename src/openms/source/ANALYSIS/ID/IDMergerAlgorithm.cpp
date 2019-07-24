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
// $Maintainer: Julianus Pfeuffer $
// $Authors: Julianus Pfeuffer $
// --------------------------------------------------------------------------

#include <OpenMS/ANALYSIS/ID/IDMergerAlgorithm.h>
#include <OpenMS/METADATA/PeptideIdentification.h>
#include <unordered_map>

using namespace std;
namespace OpenMS
{


  //TODO parameterize so it only adds/keeps best per peptide, peptide charge, modpeptide
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
    defaults_.setValidStrings("annotate_origin", ListUtils::create<String>("true,false"));
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
    // the IDs, then erase them so we dont encounter them in
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

      bool annotated = pid.metaValueExists("map_index");
      if (annotate_origin || annotated)
      {
        Size oldFileIdx(0);
        const StringList& origins = originFiles[runIdxIt->second];
        if (annotated)
        {
          oldFileIdx = pid.getMetaValue("map_index");
        }
        else if (origins.size() > 1)
        {
          // If there is more than one possible file it might be from
          // and it is not annotated -> fail
          throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Trying to annotate new map_index for PeptideIdentification "
              "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ") but"
              "no old map_index present");
        }

        if (oldFileIdx >= origins.size())
        {
          throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Trying to annotate new map_index for PeptideIdentification "
              "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ") but"
              " the index exceeds the number of files in the run.");
        }
        pid.setMetaValue("map_index", file_origin_to_idx_[origins[oldFileIdx]]);
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

      bool annotated = pid.metaValueExists("map_index");
      if (annotate_origin || annotated)
      {
        Size oldFileIdx(0);
        if (annotated)
        {
          oldFileIdx = pid.getMetaValue("map_index");
        }
          // If there is more than one possible file it might be from
          // and it is not annotated -> fail
        else if (originFiles[oldProtRunIdx].size() > 1)
        {
          throw Exception::MissingInformation(
              __FILE__,
              __LINE__,
              OPENMS_PRETTY_FUNCTION,
              "Trying to annotate new map_index for PeptideIdentification "
              "(" + String(pid.getMZ()) + ", " + String(pid.getRT()) + ") but"
              "no old map_index present");
        }
        pid.setMetaValue("map_index", file_origin_to_idx_[originFiles[oldProtRunIdx].at(oldFileIdx)]);
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
    const String& engine = ref.getSearchEngine();
    const String& version = ref.getSearchEngineVersion();
    ProteinIdentification::SearchParameters params = ref.getSearchParameters();
    set<String> fixed_mods(params.fixed_modifications.begin(), params.fixed_modifications.end());
    set<String> var_mods(params.variable_modifications.begin(), params.variable_modifications.end());
    bool ok = false;
    unsigned runID = 0;
    //TODO if one equals the ref, continue
    for (const auto& idRun : protRuns)
    {
      ok = true;
      if (idRun.getSearchEngine() != engine || idRun.getSearchEngineVersion() != version)
      {
        ok = false;
       OPENMS_LOG_WARN << "Search engine " + idRun.getSearchEngine() + " from IDRun " + String(runID) + " does not match "
        "with the others. You probably do not want to merge the results with this tool." << std::endl;
        break;
      }
      const ProteinIdentification::SearchParameters& sp = idRun.getSearchParameters();
      if (params.precursor_mass_tolerance != sp.precursor_mass_tolerance ||
          params.precursor_mass_tolerance_ppm != sp.precursor_mass_tolerance_ppm ||
          params.db != sp.db ||
          params.db_version != sp.db_version ||
          params.fragment_mass_tolerance != sp.fragment_mass_tolerance ||
          params.fragment_mass_tolerance_ppm != sp.fragment_mass_tolerance_ppm ||
          params.charges != sp.charges ||
          params.digestion_enzyme != sp.digestion_enzyme ||
          params.taxonomy != sp.taxonomy)
      {
        ok = false;
       OPENMS_LOG_WARN << "Searchengine settings from IDRun " + String(runID) + " does not match with the others."
        " You probably do not want to merge the results with this tool if they differ significantly." << std::endl;
        break;
      }

      set<String> curr_fixed_mods(sp.fixed_modifications.begin(), sp.fixed_modifications.end());
      set<String> curr_var_mods(sp.variable_modifications.begin(), sp.variable_modifications.end());
      if (fixed_mods != curr_fixed_mods ||
          var_mods != curr_var_mods)
      {
        if (experiment_type != "labeled_MS1")
        {
          ok = false;
         OPENMS_LOG_WARN << "Used modification settings from IDRun " + String(runID) + " does not match with the others."
          " Since the experiment is not annotated as MS1-labeled "
          "you probably do not want to merge the results with this tool." << std::endl;
          break;
        }
        else
        {
          //TODO actually introduce a flag for labelling modifications in the Mod datastructures?
          //OR put a unique ID for the used mod as a UserParam to the mapList entries (consensusHeaders)
          //TODO actually you would probably need an experimental design here, because
          //settings have to agree exactly in a FractionGroup but can slightly differ across runs.
         OPENMS_LOG_WARN << "Used modification settings from IDRun " + String(runID) + " does not match with the others."
          " Although it seems to be an MS1-labeled experiment,"
          " check carefully that only non-labelling mods differ." << std::endl;
        }
      }
    }
    if (!ok /*&& TODO and no force flag*/)
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
