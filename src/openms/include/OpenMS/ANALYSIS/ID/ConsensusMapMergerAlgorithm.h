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
   * @todo This could be merged in the future with the general IDMergerAlgorithm since it shares a lot.
   *     IDMergerAlgorithm needs additional methods to have multiple runs as output. It also needs to store
   *     an extended mapping internally to distribute the PeptideIDs to the right output run according to origin and
   *     label.
   *     And should have non-copying/moving
   *     overloads for inserting PeptideIDs since we probably do not want to distribute the PeptideIDs to the features
   *     again. In general detaching IDs from features would be of great help here.
   * @todo Untested for TMT/iTraq data where you usually have one Identification run per File but in one File you
   *     might have multiple conditions multiplexed, that you might want to split for inference. Problem:
   *     There is only one PeptideIdentification object per Feature that is representative for all "submaps" (in this
   *     case the labels/reporter ions). -> A lookup is necessary if the reporter ion had non-zero intensity and if
   *     so, the peptide ID needs to be duplicated for every new (condition-based) IdentificationRun it is supposed
   *     to be used in, according to the mapping.
   */
  class OPENMS_DLLAPI ConsensusMapMergerAlgorithm:
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:
    ConsensusMapMergerAlgorithm ();

    /// Takes a ConsensusMap (with usually one IdentificationRun per column [= sub map]
    /// and merges them to one IdentificationRun per Condition (unique entries in Sample section when removing
    /// replicate columns) while reassociating the PeptideHits accordingly.
    /// It just does not make sense to have every protein duplicated.
    /// And IdentificationRuns are used to guide inference methods on what identifications to perform inference on
    /// @note Constructs the mapping based on the exp. design and then uses mergeProteinIDRuns
    /// @note Groups are not carried over during merging.
    /// @throws MissingInformationException for e.g. missing map_indices in PeptideIDs
    void mergeProteinsAcrossFractionsAndReplicates(ConsensusMap& cmap, const ExperimentalDesign& exp_design) const;

    /// Similar to above, merges every ID Run into one big run. Proteins get only inserted once but Peptides stay unfiltered
    /// i.e. might occur in several PeptideIdentifications afterwards
    /// @note Groups are not carried over during merging.
    /// @throws MissingInformationException for e.g. missing map_indices in PeptideIDs
    void mergeAllIDRuns(ConsensusMap& cmap) const;

    /// Takes a ConsensusMap and a mapping between ConsensusMap column index (sub map index in the map_index meta value)
    /// and the new ProteinIdentification run index and merges them based on this.
    /// @throws MissingInformationException for e.g. missing map_indices in PeptideIDs
    /// @todo Do we need to consider the old IDRun identifier in addition to the sub map index
    void mergeProteinIDRuns(ConsensusMap &cmap,
                            const std::map<unsigned, unsigned>& mapIdx_to_new_protIDRun) const;

  private:

    /// Checks consistency of search engines and settings across runs before merging.
    /// Uses the first run as reference and compares all to it.
    /// @return all same? TODO: return a merged RunDescription about what to put in the new runs (e.g. for SILAC)
    /// @throws BaseException for disagreeing settings
    bool checkOldRunConsistency_(const std::vector<ProteinIdentification>& protRuns, const String& experiment_type) const;
    /// Same as above but with specific reference run
    bool checkOldRunConsistency_(const std::vector<ProteinIdentification>& protRuns, const ProteinIdentification& ref, const String& experiment_type) const;
    /// Compares exactly two runs @todo refactor with above
    bool checkRunSettings_(const ProteinIdentification& idRun, const ProteinIdentification& ref, const String& experiment_type) const;


    static size_t accessionHash_(const ProteinHit& p)
    {
      return std::hash<String>()(p.getAccession());
    }
    static bool accessionEqual_(const ProteinHit& p1, const ProteinHit& p2)
    {
      return p1.getAccession() == p2.getAccession();
    }
    using hash_type = std::size_t (*)(const ProteinHit&);
    using equal_type = bool (*)(const ProteinHit&, const ProteinHit&);

  };
} // namespace OpenMS
