// --------------------------------------------------------------------------
//                   OpenMS -- Open-Source Mass Spectrometry
// --------------------------------------------------------------------------
// Copyright The OpenMS Team -- Eberhard Karls University Tuebingen,
// ETH Zurich, and Freie Universitaet Berlin 2002-2018.
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

#pragma once

#include <OpenMS/DATASTRUCTURES/DefaultParamHandler.h>
#include <OpenMS/CONCEPT/Exception.h>
#include <OpenMS/CONCEPT/ProgressLogger.h>
#include <OpenMS/KERNEL/ConsensusMap.h>
#include <OpenMS/METADATA/ExperimentalDesign.h>

#include <unordered_set>

namespace OpenMS
{
  /**
   * @brief Merges identification data in ConsensusMaps
   * Has some stuff in common with IDMergerAlgorithm and therefore could be merged
   * but you can save some overhead by only going through the CMap once. Therefore
   * the extra class.
   */
  class OPENMS_DLLAPI ConsensusMapMergerAlgorithm:
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:
    ConsensusMapMergerAlgorithm ();
    /// Takes a cmap with one IDRun per column and merges them to one proteinIDRun per Condition
    /// while reassociating the PeptideEvidences
    /// Constructs the mapping based on the exp. design and then uses mergeProteinIDRuns
    /// @throws MissingInformationException for e.g. missing map_indices in PeptideIDs
    void mergeProteinsAcrossFractionsAndReplicates(ConsensusMap& cmap, const ExperimentalDesign& exp_design) const;

    /// Similar to above, merges every ID Run into one big run. Proteins get only inserted once but Peptides stay unfiltered
    /// i.e. might occur in several PeptideIdentifications afterwards
    /// @throws MissingInformationException for e.g. missing map_indices in PeptideIDs
    void mergeAllIDRuns(ConsensusMap& cmap) const;

    /// Takes a ConsensusMap and a mapping between ConsensusMap column index (map index) and
    /// the new ProteinIdentificationRun index and merges them. If you know the number of resulting
    /// ProteinIDRuns, consider specifying new_size.
    /// @throws MissingInformationException for e.g. missing map_indices in PeptideIDs
    void mergeProteinIDRuns(ConsensusMap &cmap,
                            const std::map<unsigned, unsigned>& mapIdx_to_new_protIDRun) const;

    /// Takes a vector of old protein ID runs and old peptide ID runs, which will be moved or overwritten
    /// and a map from old run to new run, as well as a to-be-filled vector of peptide IDs
    /// It merges the proteins from runs that map to the same new run (by moving the first occurence to it)
    /// It concatenates and moves the peptides of those runs into the vector at the according index while updating their
    /// run references.
    void mergeIDRunsAndSplitPeptides(
        std::vector<ProteinIdentification>& oldProtRuns,
        std::vector<PeptideIdentification>& pepIDs,
        const std::map<Size, Size>& oldrunToNewrun,
        std::vector<std::vector<PeptideIdentification>>& splitPepIDs) const;

    void mergeIDRunsAndSplitPeptides(
        std::vector<ProteinIdentification>& oldProtRuns,
        std::vector<std::vector<PeptideIdentification>>& pepIDs,
        const std::map<Size, Size>& oldrunToNewrun,
        std::vector<std::vector<PeptideIdentification>>& splitPepIDs) const;

  private:

    /// Moves multiple ID vectors into a long one
    template<class Identification>
    static void concatenateIdentifications_(
        std::vector<std::vector<Identification>>&& oldIDs,
        std::vector<Identification>& newIDs)
    {
      for (auto& IDs : oldIDs)
      {
        newIDs.reserve(newIDs.size() + IDs.size());
        std::copy(std::make_move_iterator(begin(IDs)),
                  std::make_move_iterator(end(IDs)),
                  std::back_inserter(newIDs));
      }
    }

    /// Checks consistency of search engines and settings across runs before merging.
    /// @return all same? TODO: return a merged RunDescription about what to put in the new runs (e.g. for SILAC)
    /// @throws BaseException for disagreeing settings
    bool checkOldRunConsistency_(const std::vector<ProteinIdentification>& protRuns, const String& experiment_type) const;
    bool checkOldRunConsistency_(const std::vector<ProteinIdentification>& protRuns, const ProteinIdentification& ref, const String& experiment_type) const;
    bool checkRunSettings_(const ProteinIdentification& idRun, const ProteinIdentification& ref, const String& experiment_type) const;

    void initNewRunsAndFileMappings_(
        const std::vector<ProteinIdentification>& oldProtRuns,
        const std::map<Size,Size>& oldrunToNewrun,
        std::vector<std::map<Size, Size>>& oldToNewFileIdx,
        std::vector<ProteinIdentification>& newProtIDRuns) const;

    /// In (will be moved and cleared): pepIDs, oldProtRuns
    /// Out: newProtIDRuns, splitPepIDs
    /// Helpers: oldrunToNewrun (mapping of the run itself), oldToNewFileIdx (map of fileidx in old Run to fileidx in
    /// new run), proteinsCollected (keeps track of already inserted proteins)
    void movePepIDsAndRefProteinsToResult_(
        std::vector<PeptideIdentification>& pepIDs,
        std::vector<ProteinIdentification>& oldProtRuns,
        std::vector<ProteinIdentification>& newProtIDRuns,
        std::vector<std::vector<PeptideIdentification>>& splitPepIDs,
        const std::map<Size,Size>& oldrunToNewrun,
        const std::vector<std::map<Size,Size>>& oldToNewFileIdx,
        std::vector<std::unordered_set<std::string>> proteinsCollected
    ) const;


    static size_t accessionHash(const ProteinHit& p)
    {
      return std::hash<String>()(p.getAccession());
    }
    static bool accessionEqual(const ProteinHit& p1, const ProteinHit& p2)
    {
      return p1.getAccession() == p2.getAccession();
    }
    using hash_type = std::size_t (*)(const ProteinHit&);
    using equal_type = bool (*)(const ProteinHit&, const ProteinHit&);

  };
} // namespace OpenMS
