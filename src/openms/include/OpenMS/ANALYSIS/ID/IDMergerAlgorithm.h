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

namespace OpenMS
{
  class OPENMS_DLLAPI IDMergerAlgorithm:
    public DefaultParamHandler,
    public ProgressLogger
  {
  public:
    IDMergerAlgorithm ();
    /// Takes a cmap with one IDRun per column and merges them to one proteinIDRun per Condition
    /// while reassociating the PeptideEvidences
    /// Constructs the mapping based on the exp. design and then uses mergeProteinIDRuns
    /// @throws MissingInformationException for e.g. missing map_indices in PeptideIDs
    void mergeProteinsAcrossFractionsAndReplicates(ConsensusMap& cmap, const ExperimentalDesign& exp_design) const;

    /// Similar to above, merges every ID Run into one big run. Proteins get only inserted once but Peptides stay unfiltered
    /// i.e. might occur in several PeptideIdentifications afterwards
    /// @throws MissingInformationException for e.g. missing map_indices in PeptideIDs
    void mergeAllIDRuns(ConsensusMap& cmap) const;
    void mergeAllIDRuns(std::vector<ProteinIdentification>& protRuns, std::vector<PeptideIdentification>& pepIDs) const;

    /// Takes a ConsensusMap and a mapping between ConsensusMap column index (map index) and
    /// the new ProteinIdentificationRun index and merges them. If you know the number of resulting
    /// ProteinIDRuns, consider specifying new_size.
    /// @throws MissingInformationException for e.g. missing map_indices in PeptideIDs
    void mergeProteinIDRuns(ConsensusMap &cmap,
                            std::map<unsigned, unsigned> const &mapIdx_to_new_protIDRun) const;

    //TODO Add methods for to merge vectors of PepIDs based on experimental design

  private:
    struct RunDescription
    {
      String engine;
      String version;
      ProteinIdentification::SearchParameters params;
    };

    /// Checks consistency of search engines and settings across runs before merging.
    /// @return a merged RunDescription about what to put in the new runs
    /// @throws BaseException for disagreeing settings
    bool checkOldRunConsistency_(const std::vector<ProteinIdentification> protRuns, String experiment_type) const;


  };
} // namespace OpenMS
